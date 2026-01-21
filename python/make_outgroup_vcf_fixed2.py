
# 13 Making outgroup vcf usinf peregrin falcon assembly

"""
make_outgroup_vcf_fixed2.py

Robust SAM -> single-sample outgroup VCF generator.

Inputs:
 - SAM file (from minimap2 mapping of ref windows -> outgroup assembly)
 - sites_lookup.tsv : CHROM<TAB>POS<TAB>REF<TAB>ALT  (one line per site)
 - original focal VCF (used to copy header/meta lines)
Outputs:
 - out_vcf (uncompressed). Use bgzip+tabix afterwards.

Example:
python3 make_outgroup_vcf_fixed2.py --sam queries_on_outgroup.sam \
  --outgroup_fa /path/to/peregrine_falcon.fasta \
  --sites_lookup sites_lookup.tsv \
  --orig_vcf /path/to/focal.vcf.gz \
  --out_vcf outgroup_calls.vcf --min_dp 3 --min_mapq 0 --half_window 100
"""
import argparse, re, subprocess, sys
from collections import namedtuple, defaultdict

Aln = namedtuple("Aln", "qname rname pos mapq cigar flag qlen seq")

def parse_sites_lookup(fn):
    sites = []
    qname_map = {}  # multiple possible qname forms -> (chrom,pos,ref,alt)
    with open(fn) as fh:
        for line in fh:
            line=line.strip()
            if not line: continue
            parts=line.split('\t')
            if len(parts) < 4: continue
            chrom = parts[0]; pos=int(parts[1]); ref=parts[2]; alt=parts[3]
            sites.append((chrom,pos,ref,alt))
        
            qname_map[f"{chrom}_{pos}"] = (chrom,pos,ref,alt)
            qname_map[f"{chrom}:{pos}"] = (chrom,pos,ref,alt)
            qname_map[f"{chrom}|{pos}"] = (chrom,pos,ref,alt)
            
            qname_map[f"{chrom}_{pos}_center"] = (chrom,pos,ref,alt)
            qname_map[f"{chrom}_center_{pos}"] = (chrom,pos,ref,alt)
    return sites, qname_map

def parse_sam_best(samfile):
    best = {}
    with open(samfile) as fh:
        for L in fh:
            if L.startswith('@'): continue
            cols = L.rstrip("\n").split("\t")
            if len(cols) < 11: continue
            qname = cols[0]
            flag = int(cols[1])
            rname = cols[2]
            pos = int(cols[3]) if cols[3].isdigit() else 0
            mapq = int(cols[4]) if cols[4].isdigit() else 0
            cigar = cols[5]
            seq = cols[9]
            if rname == '*' or cigar == '*':
                continue
            # compute aligned length
            ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
            aligned_len = sum(int(n) for n,o in ops if o in ('M','=','X'))
            cur = best.get(qname)
            aln = Aln(qname,rname,pos,mapq,cigar,flag,len(seq),seq)
            if cur is None:
                best[qname] = (aln, mapq, aligned_len)
            else:
                if mapq > cur[1] or (mapq == cur[1] and aligned_len > cur[2]):
                    best[qname] = (aln, mapq, aligned_len)
    # return only the Aln objects
    return {k:v[0] for k,v in best.items()}

def parse_cigar_ops(cigar):
    return [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]

def qidx_to_refpos(cigar, ref_start1, qlen, qidx1):
    # build mapping qpos -> refpos (1-based).
    ops = parse_cigar_ops(cigar)
    qpos = 1
    rpos = ref_start1
    q_to_r = [None]*(qlen+1)

    for op_len, op in ops:
        if op in ('M','=','X'):
            for i in range(op_len):
                if qpos <= qlen:
                    q_to_r[qpos] = rpos
                qpos += 1
                rpos += 1
        elif op == 'I':
            for i in range(op_len):
                if qpos <= qlen:
                    q_to_r[qpos] = None
                qpos += 1
        elif op in ('D','N'):
            rpos += op_len
        elif op == 'S':
            for i in range(op_len):
                if qpos <= qlen:
                    q_to_r[qpos] = None
                qpos += 1
        else:
            # H,P etc ignore
            pass

    # if exact qidx maps:
    if 1 <= qidx1 <= qlen and q_to_r[qidx1] is not None:
        return q_to_r[qidx1]

    # else search outwards
    radius = max(20, int(0.05 * qlen))
    for d in range(1, radius+1):
        for j in (qidx1 - d, qidx1 + d):
            if 1 <= j <= qlen and q_to_r[j] is not None:
                return q_to_r[j]
    return None

def fetch_base(fa, contig, pos1):
    try:
        p = subprocess.run(['samtools','faidx',fa,f"{contig}:{pos1}-{pos1}"],
                          capture_output=True, text=True, check=True)
        lines = p.stdout.splitlines()
        if len(lines) >= 2:
            return lines[1].strip().upper()
    except subprocess.CalledProcessError:
        return None
    return None

def get_vcf_header(orig_vcf):
    try:
        # Try bcftools first
        out = subprocess.run(['bcftools','view','-h', orig_vcf],
                            capture_output=True, text=True, check=True)
        return out.stdout
    except subprocess.CalledProcessError:
        # If bcftools fails, try zgrep for compressed VCF
        try:
            import gzip
            if orig_vcf.endswith('.gz'):
                with gzip.open(orig_vcf, 'rt') as f:
                    header_lines = []
                    for line in f:
                        if line.startswith('#'):
                            header_lines.append(line)
                        else:
                            break
                    return ''.join(header_lines)
            else:
                with open(orig_vcf) as f:
                    header_lines = []
                    for line in f:
                        if line.startswith('#'):
                            header_lines.append(line)
                        else:
                            break
                    return ''.join(header_lines)
        except:
            # fallback minimal header
            return "##fileformat=VCFv4.2\n"

def write_out_vcf(out_vcf, header_meta, sites, qname_map, best_map, out_fa, min_dp, min_mapq, half_window):
    called = 0
    total = 0
    found_map = 0

    # regex to extract chrom and pos from query names
    qname_extract = re.compile(r'^([^_:\|]+)[_:\|](\d+)')

    with open(out_vcf,'w') as out:
        meta_lines = header_meta.rstrip('\n').split('\n')
        # clean any existing #CHROM line to avoid duplication
        meta_lines = [l for l in meta_lines if not l.startswith("#CHROM")]
        # add OG_BASE INFO if missing
        if not any(l.startswith('##INFO=<ID=OG_BASE') for l in meta_lines):
            meta_lines.append('##INFO=<ID=OG_BASE,Number=1,Type=String,Description="Outgroup base at mapped position (NA if not available)">')
        if not any(l.startswith('##FORMAT=<ID=DP') for l in meta_lines):
            meta_lines.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth for outgroup heuristic">')
        for l in meta_lines:
            out.write(l + '\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tOUTGROUP\n")

        # go through sites in order and try to determine OG base
        for chrom, pos, ref, alt in sites:
            total += 1
            og_base = "NA"
            info = "OG_BASE=NA"
            gt = "./.:0"

            candidates = [f"{chrom}_{pos}", f"{chrom}:{pos}", f"{chrom}|{pos}"]

            # see if any of those were mapped
            best_aln = None
            for q in candidates:
                if q in best_map:
                    best_aln = best_map[q]
                    break

            if best_aln is None:
                # try to find SAM qname that encodes chrom and pos in other formats
                # search best_map keys for those matching pattern
                for q in best_map.keys():
                    m = qname_extract.match(q)
                    if not m: continue
                    qc = m.group(1)
                    qp = m.group(2)
                    try:
                        if qc == chrom and int(qp) == pos:
                            best_aln = best_map[q]
                            break
                    except:
                        continue

            if best_aln:
                found_map += 1
                aln = best_aln
                if aln.mapq >= min_mapq:
                    center_qpos = (aln.qlen // 2) + 1
                    refpos = qidx_to_refpos(aln.cigar, aln.pos, aln.qlen, center_qpos)
                    if refpos:
                        base = fetch_base(out_fa, aln.rname, refpos)
                        if base:
                            og_base = base
                            info = f"OG_BASE={og_base}"
                            dp = min_dp if min_dp else 10
                            if og_base == ref:
                                gt = f"0/0:{dp}"
                            elif alt != "." and og_base == alt:
                                gt = f"1/1:{dp}"
                            else:
                                gt = "./.:0"
                            called += 1

            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\tGT:DP\t{gt}\n")

    sys.stderr.write(f"[info] outgroup total_sites={total} mapped_candidates={found_map} called={called}\n")

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sam', required=True)
    p.add_argument('--outgroup_fa', required=True)
    p.add_argument('--sites_lookup', required=True)
    p.add_argument('--orig_vcf', required=True)
    p.add_argument('--out_vcf', required=True)
    p.add_argument('--min_dp', type=int, default=3)
    p.add_argument('--min_mapq', type=int, default=0)
    p.add_argument('--half_window', type=int, default=100)
    args = p.parse_args()

    sites, qname_map = parse_sites_lookup(args.sites_lookup)
    best_map = parse_sam_best(args.sam)
    hdr = get_vcf_header(args.orig_vcf)
    write_out_vcf(args.out_vcf, hdr, sites, qname_map, best_map, args.outgroup_fa,
                  args.min_dp, args.min_mapq, args.half_window)
    print("WROTE OUTGROUP VCF:", args.out_vcf, file=sys.stderr)

if __name__ == '__main__':
    main()
