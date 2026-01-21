
# 14 Polarizing focal vcf using outgroup vcf
"""
polarize_and_mutload_fixed.py

Polarize focal VCF using outgroup VCF (OG_BASE in INFO) and compute per-sample mutation load.
Outputs:
 - per_sample_mutload_summary.csv
 - per_sample_mutload_bootstrap_summary.csv
 - per_location_aggregates.csv
"""
import argparse, os
import numpy as np, pandas as pd
from cyvcf2 import VCF
from tqdm import tqdm

IMPACTS = ['HIGH','MODERATE','LOW','MODIFIER']

def read_callable_map(callable_dir):
    mapping={}
    summary = os.path.join(callable_dir, "summary.txt")
    if os.path.exists(summary):
        for line in open(summary):
            parts=line.strip().split()
            if len(parts)>=2:
                fname=os.path.basename(parts[0])
                sample=fname.split('.')[0]
                try:
                    mapping[sample]=int(parts[1])
                except:
                    pass
    else:
        for fn in os.listdir(callable_dir):
            if fn.endswith('.counts.txt'):
                sample=fn.split('.')[0]
                try:
                    mapping[sample]=int(open(os.path.join(callable_dir,fn)).read().strip())
                except:
                    pass
    return mapping

def load_annotation_lookup(annot_vcf):
    lookup={}
    if not annot_vcf or not os.path.exists(annot_vcf):
        return lookup
    v = VCF(annot_vcf)
    for rec in v:
        if len(rec.ALT)!=1: continue
        key=(rec.CHROM, rec.POS, rec.REF, rec.ALT[0])
        ann = rec.INFO.get('ANN')
        worst=None
        if ann:
            annstr = ann if isinstance(ann,str) else (ann[0] if isinstance(ann,(list,tuple)) else str(ann))
            for entry in str(annstr).split(','):
                fields = entry.split('|')
                if len(fields) > 2:
                    impact = fields[2].upper()
                    if impact in IMPACTS:
                        if worst is None or IMPACTS.index(impact) < IMPACTS.index(worst):
                            worst = impact
        if worst is None: worst='MODIFIER'
        lookup[key]=worst
    v.close()
    return lookup

def load_outgroup_lookup(outgroup_vcf):
    """Load outgroup base calls from VCF INFO field OG_BASE"""
    og_lookup = {}
    if not outgroup_vcf or not os.path.exists(outgroup_vcf):
        return og_lookup

    try:
        v = VCF(outgroup_vcf)
        total = 0
        called = 0

        for rec in v:
            total += 1
            if len(rec.ALT) != 1:
                key = (rec.CHROM, rec.POS, rec.REF, rec.ALT[0] if rec.ALT else None)
                og_lookup[key] = None
                continue

            # Get OG_BASE from INFO
            og = rec.INFO.get('OG_BASE')
            key = (rec.CHROM, rec.POS, rec.REF, rec.ALT[0])

            if og is None:
                og_lookup[key] = None
            else:
                # Handle different formats
                if isinstance(og, (list, tuple)):
                    ogval = og[0]
                else:
                    ogval = og

                if ogval is None or ogval == 'NA':
                    og_lookup[key] = None
                else:
                    og_lookup[key] = ogval.upper()
                    called += 1

        v.close()
        print(f"[info] outgroup entries total={total} called={called}")

    except Exception as e:
        print(f"[warn] could not open outgroup vcf {outgroup_vcf}: {e}")
        og_lookup = {}

    return og_lookup

def remap_location(sample):
    if sample.startswith('EG') or sample.startswith('EV') or sample.startswith('BKN'):
        return 'BKN'
    if sample.startswith('WR') or sample.startswith('HYD'):
        return 'HYD'
    return sample.split('_')[0] if '_' in sample else sample

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--vcf', required=True)
    p.add_argument('--outgroup_vcf', default="")
    p.add_argument('--annot_vcf', default="")
    p.add_argument('--callable_dir', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--bootstrap', type=int, default=2000)
    p.add_argument('--min_dp', type=int, default=3)
    p.add_argument('--min_mapq', type=int, default=20)
    p.add_argument('--exclude', default="")
    p.add_argument('--weights', default="1.0,0.5,0.1,0.0")
    args=p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    callable_map = read_callable_map(args.callable_dir)
    annot = load_annotation_lookup(args.annot_vcf) if args.annot_vcf else {}
    weights_list = [float(x) for x in args.weights.split(',')]
    weights = {'HIGH':weights_list[0], 'MODERATE':weights_list[1],
               'LOW':weights_list[2], 'MODIFIER':weights_list[3]}

    # Load outgroup data
    og_lookup = load_outgroup_lookup(args.outgroup_vcf)

    # Read focal VCF
    v = VCF(args.vcf)
    samples = v.samples[:]
    excludes = set([x.strip() for x in args.exclude.split(',') if x.strip()])
    included = [s for s in samples if s not in excludes]

    print(f"[info] total samples: {len(samples)}, included: {len(included)}, excluded: {excludes}")

    site_counts = []
    site_hom = []
    site_impacts = []
    n_total = n_bial = n_pol = 0

    for rec in tqdm(v, desc="iter vcf"):
        n_total += 1
        if len(rec.ALT) != 1:
            continue
        n_bial += 1

        key = (rec.CHROM, rec.POS, rec.REF, rec.ALT[0])
        og_base = og_lookup.get(key, None) if og_lookup else None

        # Skip if outgroup provided but no og_base
        if og_lookup and og_base is None:
            continue

        imp = annot.get(key, 'MODIFIER')
        gtlist = rec.genotypes
        counts = []
        homs = []

        for s in included:
            idx = samples.index(s)
            try:
                gt = gtlist[idx]
            except:
                gt = None

            if gt is None:
                counts.append(0)
                homs.append(0)
                continue

            a1, a2, _ = gt[:3]
            if a1 is None or a2 is None or int(a1) < 0 or int(a2) < 0:
                counts.append(0)
                homs.append(0)
                continue

            # Count derived alleles
            alleles = [rec.REF] + list(rec.ALT)
            derived = 0

            for ai in (int(a1), int(a2)):
                if ai < 0:
                    continue
                base = alleles[ai]

                if og_lookup:  # Outgroup available
                    if base != og_base:
                        derived += 1
                else:  # No outgroup, treat ALT as derived
                    if ai != 0:
                        derived += 1

            counts.append(derived)
            homs.append(1 if derived == 2 else 0)

        site_counts.append(counts)
        site_hom.append(homs)
        site_impacts.append(imp)
        n_pol += 1

    v.close()

    print(f"[info] totals: total_sites={n_total}, biallelic={n_bial}, polarized_and_callable={n_pol}")

    if n_pol == 0:
        print("[error] no polarized sites found; exiting")
        return

    # Convert to numpy arrays
    counts_mat = np.array(site_counts, dtype=int)
    hom_mat = np.array(site_hom, dtype=int)
    totals_nonref = counts_mat.sum(axis=0)
    totals_hom = hom_mat.sum(axis=0)
    n_samples = len(included)

    # Count homozygous sites by impact
    per_imp_hom = {imp: np.zeros(n_samples, dtype=int) for imp in IMPACTS}
    for i, imp in enumerate(site_impacts):
        if imp not in per_imp_hom:
            imp = 'MODIFIER'
        per_imp_hom[imp] += hom_mat[i, :]

    # Write per-sample summary
    rows = []
    for i, s in enumerate(included):
        cbp = callable_map.get(s, np.nan)
        perMb = totals_nonref[i] / (cbp / 1e6) if (not np.isnan(cbp) and cbp > 0) else np.nan

        r = {
            'sample': s,
            'location': remap_location(s),
            'callable_bp': int(cbp) if (not np.isnan(cbp) and cbp) else None,
            'total_nonref': int(totals_nonref[i]),
            'total_hom': int(totals_hom[i]),
            'total_nonref_perMb': perMb
        }

        for imp in IMPACTS:
            r[f'hom_{imp}'] = int(per_imp_hom[imp][i])

        rows.append(r)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.outdir, 'per_sample_mutload_summary.csv'), index=False)
    print("[info] wrote per_sample_mutload_summary.csv")

    # Bootstrap resampling
    if args.bootstrap > 0:
        rng = np.random.default_rng(123456)
        nboot = args.bootstrap
        n_sites = counts_mat.shape[0]

        boot_totals = np.zeros((nboot, n_samples), dtype=int)
        boot_hom_imp = {imp: np.zeros((nboot, n_samples), dtype=int) for imp in IMPACTS}

        for b in tqdm(range(nboot), desc="bootstrap"):
            idxs = rng.choice(n_sites, size=n_sites, replace=True)
            bc = counts_mat[idxs, :].sum(axis=0)
            bh = hom_mat[idxs, :].sum(axis=0)
            boot_totals[b, :] = bc

            for imp in IMPACTS:
                mask = np.array([1 if site_impacts[ii] == imp else 0 for ii in idxs], dtype=int)
                arr = hom_mat[idxs, :] * mask[:, None]
                boot_hom_imp[imp][b, :] = arr.sum(axis=0)

        # Write bootstrap summary
        boot_rows = []
        for i, s in enumerate(included):
            mean_tot = float(boot_totals[:, i].mean())
            ci = np.percentile(boot_totals[:, i], [2.5, 97.5]).astype(int)

            r = {
                'sample': s,
                'total_nonref_boot_mean': mean_tot,
                'total_nonref_ci2.5': int(ci[0]),
                'total_nonref_ci97.5': int(ci[1])
            }

            cbp = callable_map.get(s, np.nan)
            if not np.isnan(cbp) and cbp > 0:
                perMb_mean = float(np.mean(boot_totals[:, i] / (cbp / 1e6)))
                perMb_ci = np.percentile(boot_totals[:, i] / (cbp / 1e6), [2.5, 97.5])
                r['total_nonref_perMb_boot_mean'] = perMb_mean
                r['total_nonref_perMb_ci2.5'] = float(perMb_ci[0])
                r['total_nonref_perMb_ci97.5'] = float(perMb_ci[1])
            else:
                r['total_nonref_perMb_boot_mean'] = None
                r['total_nonref_perMb_ci2.5'] = None
                r['total_nonref_perMb_ci97.5'] = None

            for imp in IMPACTS:
                arr = boot_hom_imp[imp][:, i]
                r[f'{imp}_hom_ci2.5'] = int(np.percentile(arr, 2.5))
                r[f'{imp}_hom_ci97.5'] = int(np.percentile(arr, 97.5))

            boot_rows.append(r)

        bdf = pd.DataFrame(boot_rows)
        bdf.to_csv(os.path.join(args.outdir, 'per_sample_mutload_bootstrap_summary.csv'), index=False)
        print("[info] wrote per_sample_mutload_bootstrap_summary.csv")

    print("[info] DONE")

if __name__ == '__main__':
    main()
