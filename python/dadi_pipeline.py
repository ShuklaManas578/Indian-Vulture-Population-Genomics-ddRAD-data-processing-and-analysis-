
# 16 Resume capable dadi_pipeline for 4 demogrpahic models: AM, SI, SC and IM
"""
dadi_pipeline.py
- Resume logic: if --resume and model_{X}/popt.npy exists, the script loads popt and computes
  model SFS, theta and log-likelihood, then runs bootstrap+plots (no re-optimization).
"""
from __future__ import annotations
import os, sys, argparse, time, math, traceback
import numpy as np
import multiprocessing as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    import dadi
except Exception as e:
    sys.stderr.write("ERROR: dadi import failed. Install dadi in this python environment.\n")
    raise

# utilities 
def now(): return time.strftime("%Y-%m-%d %H:%M:%S")
def safe_print(*a, **k): print(*a, **k); sys.stdout.flush()

# robust file load/write
def robust_load_2d_matrix(path):
    if not os.path.exists(path): raise FileNotFoundError(path)
    try:
        mat = np.genfromtxt(path, delimiter=None, autostrip=True)
        mat = np.atleast_2d(mat)
        return mat
    except Exception:
        with open(path) as fh:
            txt = fh.read().strip()
        toks = txt.replace(",", " ").split()
        arr = np.array([float(x) for x in toks], dtype=float)
        return arr

def write_dadi_fs(int_mat, out_fs):
    rows, cols = int_mat.shape
    with open(out_fs, 'w') as fh:
        fh.write(f"{rows} {cols}\n")
        for r in range(rows):
            fh.write(" ".join(str(int(int_mat[r, c])) for c in range(cols)))
            fh.write("\n")

# ensure clean numeric array
def ensure_spectrum(obj):
    arr = np.asarray(obj)
    if isinstance(arr, np.ma.MaskedArray): arr = arr.filled(0)
    arr = np.array(arr, dtype=float)
    return dadi.Spectrum(arr)

# phi wrapper compatibility
def call_phi_1d_to_2d_safe(xx, phi):
    fn = dadi.PhiManip.phi_1D_to_2D
    tries = [lambda: fn(xx, phi), lambda: fn(phi, xx), lambda: fn(phi), lambda: fn(xx)]
    last_exc = None
    for t in tries:
        try:
            return t()
        except TypeError as e:
            last_exc = e
            continue
    raise RuntimeError(f"phi_1D_to_2D signature mismatch; last error: {last_exc}")

# plotting helpers
def plot_obs_exp_res(obs_mat, exp_mat, outpng, title, pop_labels=("BKN","HYD")):
    resid = obs_mat - exp_mat
    fig, axs = plt.subplots(2,2, figsize=(10,8))
    im0 = axs[0,0].imshow(obs_mat, origin='lower', aspect='auto')
    axs[0,0].set_title('Observed SFS'); axs[0,0].set_xlabel(pop_labels[1]); axs[0,0].set_ylabel(pop_labels[0])
    plt.colorbar(im0, ax=axs[0,0])
    im1 = axs[0,1].imshow(exp_mat, origin='lower', aspect='auto')
    axs[0,1].set_title('Model SFS (scaled)'); axs[0,1].set_xlabel(pop_labels[1])
    plt.colorbar(im1, ax=axs[0,1])
    im2 = axs[1,0].imshow(resid, origin='lower', aspect='auto')
    axs[1,0].set_title('Residuals (obs - exp)'); axs[1,0].set_xlabel(pop_labels[1]); axs[1,0].set_ylabel(pop_labels[0])
    plt.colorbar(im2, ax=axs[1,0])
    axs[1,1].hist(resid.flatten(), bins=80)
    axs[1,1].set_title('Residual histogram')
    fig.suptitle(title)
    plt.tight_layout(rect=[0,0,1,0.96])
    fig.savefig(outpng, dpi=1200); plt.close(fig)

def plot_param_violin(popt, boots, param_names, outpng):
    fig = plt.figure(figsize=(12,6))
    P = len(popt)
    ax1 = fig.add_subplot(121)
    ax1.bar(range(P), popt)
    ax1.set_xticks(range(P)); ax1.set_xticklabels(param_names, rotation=45)
    ax1.set_title('Fitted params (point estimates)')
    ax2 = fig.add_subplot(122)
    boots_clean = boots[~np.isnan(boots).any(axis=1)] if boots.size else np.empty((0,P))
    if boots_clean.shape[0] > 0:
        ax2.violinplot([boots_clean[:,i] for i in range(P)], showmeans=True)
        ax2.scatter(np.arange(1,P+1), popt, marker='x')
    ax2.set_xticks(range(1,P+1)); ax2.set_xticklabels(param_names, rotation=45)
    ax2.set_title('Bootstrap parameter distributions')
    plt.tight_layout(); fig.savefig(outpng, dpi=1200); plt.close(fig)

def plot_param_histograms(popt, boots, param_names, outpng):
    P = len(popt)
    boots_clean = boots[~np.isnan(boots).any(axis=1)] if boots.size else np.empty((0,P))
    ncols = min(4, P); nrows = int(math.ceil(P / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows))
    axs = np.atleast_2d(axs).reshape(-1)
    for i in range(P):
        ax = axs[i]; name = param_names[i]
        ax.axvline(popt[i], linestyle='-', label='popt')
        if boots_clean.shape[0] > 0:
            vals = boots_clean[:, i]; ax.hist(vals, bins=30, alpha=0.8)
            lo = np.nanpercentile(vals, 2.5); hi = np.nanpercentile(vals, 97.5)
            med = np.nanmedian(vals); mean = np.nanmean(vals)
            ax.axvline(lo, color='red', linestyle='--'); ax.axvline(hi, color='red', linestyle='--')
            ax.axvline(med, color='green', linestyle=':', label='median')
            ax.text(0.98, 0.95, f"mean={mean:.3g}\nmed={med:.3g}\n95%CI=[{lo:.3g},{hi:.3g}]", transform=ax.transAxes, va='top', ha='right', fontsize=8, bbox=dict(facecolor='white', alpha=0.6))
        else:
            ax.text(0.5,0.5,"no bootstrap", ha='center', va='center')
        ax.set_title(name)
    for j in range(P, len(axs)): axs[j].axis('off')
    plt.tight_layout(); fig.savefig(outpng, dpi=1200); plt.close(fig)

def plot_pairwise_scatter(boots, param_names, outpng):
    boots_clean = boots[~np.isnan(boots).any(axis=1)] if boots.size else np.empty((0,boots.shape[1] if boots.size else 0))
    P = len(param_names)
    if boots_clean.shape[0] == 0 or P < 2:
        fig = plt.figure(figsize=(4,3)); plt.text(0.5,0.5,"No bootstrap draws to show pairwise scatter", ha='center', va='center'); plt.savefig(outpng, dpi=1200); plt.close(fig); return
    fig, axes = plt.subplots(P,P, figsize=(2.5*P, 2.5*P))
    for i in range(P):
        for j in range(P):
            ax = axes[i,j]
            if i == j: ax.hist(boots_clean[:, i], bins=30)
            else: ax.scatter(boots_clean[:, j], boots_clean[:, i], s=6, alpha=0.6)
            if i == P-1: ax.set_xlabel(param_names[j])
            else: ax.set_xticks([])
            if j == 0: ax.set_ylabel(param_names[i])
            else: ax.set_yticks([])
    plt.tight_layout(); plt.savefig(outpng, dpi=1200); plt.close(fig)

def plot_bootstrap_corr_heatmap(boots, param_names, outpng):
    boots_clean = boots[~np.isnan(boots).any(axis=1)] if boots.size else np.empty((0,boots.shape[1] if boots.size else 0))
    fig = plt.figure(figsize=(6,5))
    if boots_clean.shape[0] == 0:
        plt.text(0.5,0.5,"No bootstrap draws for correlation heatmap", ha='center', va='center')
    else:
        corr = np.corrcoef(boots_clean.T)
        im = plt.imshow(corr, origin='lower', vmin=-1, vmax=1); plt.colorbar(im)
        plt.xticks(range(len(param_names)), param_names, rotation=45); plt.yticks(range(len(param_names)), param_names)
        plt.title('Bootstrap parameter correlation')
    plt.tight_layout(); plt.savefig(outpng, dpi=1200); plt.close(fig)

# fallback model implementations
def model_SI_fallback(params, ns, pts):
    _, nu1, nu2, T = params
    xx = dadi.Numerics.default_grid(pts)
    phi1 = dadi.PhiManip.phi_1D(xx)
    phi2 = call_phi_1d_to_2d_safe(xx, phi1)
    phi2 = dadi.Integration.two_pops(phi2, xx, T, nu1, nu2, m12=0.0, m21=0.0)
    return dadi.Spectrum.from_phi(phi2, ns, (xx, xx))

def model_IM_fallback(params, ns, pts):
    _, nu1, nu2, T, m12, m21 = params
    xx = dadi.Numerics.default_grid(pts)
    phi1 = dadi.PhiManip.phi_1D(xx)
    phi2 = call_phi_1d_to_2d_safe(xx, phi1)
    phi2 = dadi.Integration.two_pops(phi2, xx, T, nu1, nu2, m12=m12, m21=m21)
    return dadi.Spectrum.from_phi(phi2, ns, (xx, xx))

def model_SC_fallback(params, ns, pts):
    _, nu1, nu2, Tsc, Tsec, m12, m21 = params
    xx = dadi.Numerics.default_grid(pts)
    phi1 = dadi.PhiManip.phi_1D(xx)
    phi2 = call_phi_1d_to_2d_safe(xx, phi1)
    phi2 = dadi.Integration.two_pops(phi2, xx, Tsc, nu1, nu2, m12=0.0, m21=0.0)
    phi2 = dadi.Integration.two_pops(phi2, xx, Tsec, nu1, nu2, m12=m12, m21=m21)
    return dadi.Spectrum.from_phi(phi2, ns, (xx, xx))

def model_AM_fallback(params, ns, pts):
    _, nu1, nu2, Tanc, Tam, m12, m21 = params
    xx = dadi.Numerics.default_grid(pts)
    phi1 = dadi.PhiManip.phi_1D(xx)
    phi2 = call_phi_1d_to_2d_safe(xx, phi1)
    phi2 = dadi.Integration.two_pops(phi2, xx, Tanc, nu1, nu2, m12=m12, m21=m21)
    phi2 = dadi.Integration.two_pops(phi2, xx, Tam, nu1, nu2, m12=0.0, m21=0.0)
    return dadi.Spectrum.from_phi(phi2, ns, (xx, xx))

FALLBACK_MODELS = {'SI': model_SI_fallback, 'IM': model_IM_fallback, 'SC': model_SC_fallback, 'AM': model_AM_fallback}

# optimization helpers (coarse->refine)
def make_model_func(model_name, force_fallback=False):
    if not force_fallback:
        try:
            model_func = getattr(dadi.Demographics2D, model_name)
            safe_print(f"[{now()}] Using dadi.Demographics2D.{model_name}")
            return model_func, False
        except Exception:
            pass
    if model_name in FALLBACK_MODELS:
        safe_print(f"[{now()}] Using fallback implementation for model {model_name}")
        return FALLBACK_MODELS[model_name], True
    raise RuntimeError(f"Model {model_name} not available and no fallback implemented.")

def get_extrap_func(model_func):
    try: return dadi.Numerics.make_extrap_func(model_func)
    except Exception:
        try: return dadi.Numerics.make_extrap_log_func(model_func)
        except Exception as e: raise RuntimeError(f"Could not create extrapolation wrapper: {e}")

def _opt_attempt(start, data_spec, func_ex, pts_l, lower, upper, maxiter):
    popt = dadi.Inference.optimize_log(start, data_spec, func_ex, pts_l, lower_bound=lower, upper_bound=upper, verbose=0, maxiter=maxiter)
    model_raw = func_ex(popt, data_spec.sample_sizes, pts_l)
    model_arr = np.asarray(model_raw)
    if np.any(np.isnan(model_arr)) or np.any(model_arr < 0):
        raise RuntimeError('Invalid model fs produced (NaN or negative)')
    model_spec = dadi.Spectrum(model_arr)
    theta = dadi.Inference.optimal_sfs_scaling(model_spec, data_spec)
    ll = dadi.Inference.ll_multinom(model_spec * theta, data_spec)
    return popt, model_arr, theta, ll

def optimize_two_stage(data_fs, ns, func_ex, p0, lower, upper, pts_l, n_restarts=8, n_coarse=None, coarse_maxiter=30, refine_maxiter=200, verbose=False, debug=False):
    data_spec = ensure_spectrum(data_fs)
    n_coarse = n_coarse if n_coarse is not None else max(4, n_restarts-2)
    starts = [p0]
    for i in range(n_coarse-1):
        starts.append(p0 * (1.0 + 0.2*(np.random.rand(len(p0)) - 0.5)))
    safe_print(f"[{now()}]  Coarse stage: {len(starts)} starts (maxiter={coarse_maxiter})")
    coarse_results = []
    for i, st in enumerate(starts):
        try:
            t0 = time.time()
            popt_c, model_arr_c, theta_c, ll_c = _opt_attempt(st, data_spec, func_ex, pts_l, lower, upper, coarse_maxiter)
            elapsed = time.time() - t0
            safe_print(f"[{now()}]   coarse attempt {i+1} ok: ll={ll_c:.3f} (time={elapsed:.1f}s)")
            coarse_results.append((popt_c, model_arr_c, theta_c, ll_c))
        except Exception as e:
            safe_print(f"[{now()}]   coarse attempt {i+1} failed: {e}")
            if debug: traceback.print_exc()
            continue
    if not coarse_results:
        raise RuntimeError('All coarse attempts failed')
    coarse_results.sort(key=lambda x: x[3], reverse=True)
    best_coarse = coarse_results[0]
    to_refine = [cr[0] for cr in coarse_results[:max(1, min(n_restarts, len(coarse_results)) )]]
    safe_print(f"[{now()}]  Refine stage: {len(to_refine)} starts (maxiter={refine_maxiter})")
    best = None; best_ll = -np.inf
    for i, st in enumerate(to_refine):
        try:
            t0 = time.time()
            popt_r = dadi.Inference.optimize_log(st, data_spec, func_ex, pts_l, lower_bound=lower, upper_bound=upper, verbose=1 if verbose else 0, maxiter=refine_maxiter)
            model_raw = func_ex(popt_r, data_spec.sample_sizes, pts_l)
            model_arr = np.asarray(model_raw)
            if np.any(np.isnan(model_arr)) or np.any(model_arr < 0):
                safe_print(f"[{now()}]   refine attempt {i+1} produced invalid model (NaN or negative) — skipping")
                continue
            model_spec = dadi.Spectrum(model_arr)
            theta_r = dadi.Inference.optimal_sfs_scaling(model_spec, data_spec)
            ll_r = dadi.Inference.ll_multinom(model_spec * theta_r, data_spec)
            elapsed = time.time() - t0
            safe_print(f"[{now()}]   refine attempt {i+1} done: ll={ll_r:.3f} (time={elapsed:.1f}s)")
            if np.isfinite(ll_r) and ll_r > best_ll:
                best_ll = ll_r; best = (popt_r, model_spec, theta_r, ll_r)
        except Exception as e:
            safe_print(f"[{now()}]   refine attempt {i+1} failed: {e}")
            if debug: traceback.print_exc()
            continue
    if best is None:
        safe_print(f"[{now()}] No refine succeeded; using best coarse result")
        popt_c, model_arr_c, theta_c, ll_c = best_coarse
        return popt_c, dadi.Spectrum(model_arr_c), theta_c, ll_c
    return best

# bootstrap parallel
def _bootstrap_worker(args):
    """
    args: (
      i, total_sites, probs (1D np array), shape (tuple), p0_guess (1d array),
      model_name (str), force_fallback (bool), pts_l (list), lower (list/array), upper (list/array), maxiter (int)
    )
    """
    i, total_sites, probs, shape, p0_guess, model_name, force_fallback, pts_l, lower, upper, maxiter = args
    try:
        # reconstruct func_ex inside worker
        model_func, _ = make_model_func(model_name, force_fallback=force_fallback)
        func_ex = get_extrap_func(model_func)
        samp = np.random.multinomial(total_sites, probs)
        samp_mat = samp.reshape(shape)
        fs = dadi.Spectrum(samp_mat)
        popt_b = dadi.Inference.optimize_log(p0_guess, fs, func_ex, pts_l, lower_bound=lower, upper_bound=upper, verbose=0, maxiter=maxiter)
        return popt_b
    except Exception:
        return np.full(len(p0_guess), np.nan)

def parametric_bootstrap_parallel(model_fs_prob, total_sites, model_name, force_fallback, p0_guess, pts_l, lower, upper, nboot=50, threads=1, verbose=False, worker_maxiter=40):
    prob_arr = np.array(model_fs_prob, dtype=float)
    probs = prob_arr.flatten()
    if probs.sum() <= 0 or np.isnan(probs).any(): raise ValueError('Model SFS probabilities invalid for bootstrap.')
    probs = probs / probs.sum()
    shape = prob_arr.shape
    args_list = []
    for i in range(nboot):
        args_list.append((i, total_sites, probs, shape, p0_guess, model_name, force_fallback, pts_l, list(lower), list(upper), worker_maxiter))
    boots = []
    if threads > 1:
        try:
            safe_print(f"[{now()}]  starting parallel bootstrap with {threads} workers")
            with mp.Pool(processes=threads) as pool:
                for res in pool.imap_unordered(_bootstrap_worker, args_list):
                    boots.append(res)
            return np.array(boots)
        except Exception as e:
            safe_print(f"[{now()}] Parallel bootstrap failed: {e}; falling back to serial")
    # serial fallback
    for i, args_ in enumerate(args_list):
        if (i+1) % 10 == 0 or verbose:
            safe_print(f"[{now()}]  Bootstrap {i+1}/{nboot}")
        res = _bootstrap_worker(args_)
        boots.append(res)
    return np.array(boots)

# read existing model summaries in outdir 
def collect_existing_models(outdir):
    collected = []
    if not os.path.isdir(outdir):
        return collected
    for name in os.listdir(outdir):
        if not name.startswith("model_"): continue
        sub = os.path.join(outdir, name)
        if not os.path.isdir(sub): continue
        summary_file = os.path.join(sub, "fit_summary.txt")
        if not os.path.exists(summary_file): continue
        model_name = name.replace("model_","")
        ll = None; aic = None
        try:
            with open(summary_file) as fh:
                for line in fh:
                    l = line.strip()
                    if not l: continue
                    low = l.lower()
                    if 'loglikelihood' in low or low.startswith('ll') or low.startswith('loglik'):
                        if ':' in l:
                            try:
                                ll = float(l.split(':',1)[1].strip())
                            except Exception:
                                pass
                        else:
                            toks = l.split()
                            try:
                                ll = float(toks[-1])
                            except Exception:
                                pass
                    if 'aic' in low:
                        if ':' in l:
                            try:
                                aic = float(l.split(':',1)[1].strip())
                            except Exception:
                                pass
                        else:
                            toks = l.split()
                            try:
                                aic = float(toks[-1])
                            except Exception:
                                pass
            if ll is not None and aic is not None:
                collected.append({'model': model_name, 'll': ll, 'aic': aic, 'resdir': sub})
        except Exception:
            continue
    return collected

def write_combined_model_summary(outdir, results_local):
    existing = collect_existing_models(outdir)
    combined = []
    # add existing but not in results_local
    for e in existing:
        if e['model'] in results_local: continue
        combined.append((e['model'], e['ll'], e['aic'], e['resdir']))
    # add local results
    for m, info in results_local.items():
        combined.append((m, info['ll'], info['aic'], info['resdir']))
    if not combined:
        safe_print(f"[{now()}] No models available for combined summary.")
        return
    tsv = os.path.join(outdir, "model_summary.tsv")
    with open(tsv, 'w') as fh:
        fh.write("model\tll\taic\tdeltaAIC\tAkaikeWeight\tresdir\n")
        aics = np.array([x[2] for x in combined])
        best = aics.min()
        delta = aics - best
        weights = np.exp(-0.5*delta) / np.exp(-0.5*delta).sum()
        for i,(m,ll,aic,resdir) in enumerate(combined):
            fh.write(f"{m}\t{ll}\t{aic}\t{delta[i]:.6g}\t{weights[i]:.6g}\t{resdir}\n")
    models = [x[0] for x in combined]
    aics = np.array([x[2] for x in combined])
    fig, ax = plt.subplots(figsize=(max(6,len(models)*1.2),4))
    ax.bar(models, aics)
    ax.set_ylabel('AIC'); ax.set_title('Model comparison (AIC)')
    for i,m in enumerate(models):
        ax.text(i, aics[i], f" w={weights[i]:.2f}", va='bottom', ha='center')
    plt.tight_layout(); plt.savefig(os.path.join(outdir,"model_comparison_AIC.png"), dpi=1200); plt.close(fig)
    safe_print(f"[{now()}] Wrote combined model summary to {tsv} and model_comparison_AIC.png")

# CLI
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--matrix', required=True)
    p.add_argument('--nBKN', type=int, required=True)
    p.add_argument('--nHYD', type=int, required=True)
    p.add_argument('--outdir', default='dadi_results_v3')
    p.add_argument('--pts', type=int, default=None)
    p.add_argument('--ptslist', default=None)
    p.add_argument('--maxiter', type=int, default=200)
    p.add_argument('--coarse-maxiter', type=int, default=30)
    p.add_argument('--nrestarts', type=int, default=8)
    p.add_argument('--ncoarse', type=int, default=None)
    p.add_argument('--nboot', type=int, default=50)
    p.add_argument('--no-bootstrap', action='store_true')
    p.add_argument('--folded', action='store_true')
    p.add_argument('--poplabels', nargs=2, default=["BKN","HYD"])
    p.add_argument('--models', default='SI,IM,SC,AM')
    p.add_argument('--force-fallback', action='store_true')
    p.add_argument('--ignore-monomorphic', action='store_true')
    p.add_argument('--verbose', action='store_true')
    p.add_argument('--debug', action='store_true')
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--im-p0', type=str, default=None, help="Comma-separated 6 numbers to set initial IM p0: s,nu1,nu2,T,m12,m21")
    p.add_argument('--ncandidates', type=int, default=3, help="Max number of IM starting candidates to try (default 3)")
    p.add_argument('--resume', action='store_true', help="If set, resume from existing model_{X}/popt.npy (skip re-optimizing that model)")
    p.add_argument('--force-rerun', action='store_true', help="Force re-running bootstrap+plots even if outputs exist")
    return p.parse_args()

# main 
def main():
    args = parse_args()
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    safe_print(f"[{now()}] START; verbose={args.verbose}, debug={args.debug}, threads={args.threads}, ncandidates={args.ncandidates}, resume={args.resume}")

    raw = robust_load_2d_matrix(args.matrix)
    R = 2*args.nBKN + 1; C = 2*args.nHYD + 1
    if isinstance(raw, np.ndarray) and raw.ndim == 1:
        if raw.size != R*C: raise SystemExit(f"ERROR: flattened input length {raw.size} != expected {R*C}")
        mat = raw.reshape((R,C))
    else:
        mat = np.atleast_2d(raw)
        if mat.shape == (R,C): pass
        elif mat.shape == (C,R): mat = mat.T
        elif mat.size == R*C: mat = mat.reshape((R,C))
        else: raise SystemExit(f"ERROR reading matrix: got shape {mat.shape}, expected ({R},{C})")

    float_mat = np.array(mat, dtype=float)
    total = float_mat.sum()
    int_mat = np.rint(float_mat).astype(int)
    diff = int(round(total)) - int(int_mat.sum())
    if diff != 0:
        idx = np.unravel_index(np.argmax(float_mat), float_mat.shape); int_mat[idx] += diff

    fs_file = os.path.join(outdir, 'data.fs'); write_dadi_fs(int_mat, fs_file)
    safe_print(f"[{now()}] Wrote dadi fs: {fs_file}")

    try:
        data_fs = dadi.Spectrum.from_file(fs_file)
    except Exception:
        data_fs = ensure_spectrum(int_mat)

    if args.folded:
        try: data_fs = data_fs.fold()
        except Exception:
            try: data_fs = ensure_spectrum(data_fs).fold()
            except Exception as e: safe_print(f"[{now()}] Warning: could not fold data_fs: {e}")

    data_fs = ensure_spectrum(data_fs)
    if args.debug:
        safe_print(f"[{now()}] DEBUG data_fs shape={np.asarray(data_fs).shape}, folded={getattr(data_fs,'folded',None)}, sample_sizes={getattr(data_fs,'sample_sizes',None)}")

    ns = (2*args.nBKN, 2*args.nHYD); mnss = max(ns)
    if args.pts is not None: pts_l = [args.pts]
    elif args.ptslist: pts_l = [int(x) for x in args.ptslist.split(',')]
    else: pts_l = [mnss + 10]
    safe_print(f"[{now()}] pts list: {pts_l}, ns={ns}")

    # default p0s
    model_p0 = {
        'SI': np.array([0.5, 1.0, 1.0, 0.1]),
        'IM': np.array([0.5, 1.0, 1.0, 0.1, 0.5, 0.5]),
        'SC': np.array([0.5, 1.0, 1.0, 0.05, 0.05, 0.5, 0.5]),
        'AM': np.array([0.5, 1.0, 1.0, 0.05, 0.05, 0.5, 0.5]),
    }
    param_names = {
        'SI': ['s','nu1','nu2','T'],
        'IM': ['s','nu1','nu2','T','m12','m21'],
        'SC': ['s','nu1','nu2','Tsc','Tsec','m12','m21'],
        'AM': ['s','nu1','nu2','Tanc','Tam','m12','m21'],
    }

    im_user_p0 = None
    if args.im_p0:
        toks = args.im_p0.replace(" ", "").split(",")
        if len(toks) != 6:
            raise SystemExit("ERROR: --im-p0 requires 6 comma-separated numbers: s,nu1,nu2,T,m12,m21")
        im_user_p0 = np.array([float(x) for x in toks], dtype=float)
        safe_print(f"[{now()}] Using user-supplied IM initial p0: {im_user_p0}")

    models_want = [m.strip() for m in args.models.split(",") if m.strip()]
    results_local = {}

    for model in models_want:
        if model not in model_p0:
            safe_print(f"[{now()}] Skipping unknown model {model}")
            continue
        safe_print(f"\n[{now()}] START model {model}")
        start_t = time.time()
        p0 = model_p0[model].copy()

        # default bounds
        lower = []; upper = []
        for name in param_names[model]:
            if name == 's': lower.append(1e-4); upper.append(0.9999)
            elif name.startswith('nu') or name in ('nu1','nu2'): lower.append(1e-5); upper.append(1e3)
            elif name.startswith('T'): lower.append(1e-6); upper.append(50.0)
            elif name.startswith('m'): lower.append(1e-8); upper.append(500.0)
            else: lower.append(1e-6); upper.append(1e3)
        lower = np.array(lower); upper = np.array(upper)

        # IM-specific tighter bounds
        used_fallback = False
        if model == 'IM':
            safe_print(f"[{now()}] Applying IM-specific tighter bounds (nu<=100, T<=10, m<=50)")
            for i, name in enumerate(param_names['IM']):
                if name in ('nu1','nu2'):
                    upper[i] = min(upper[i], 100.0)
                if name.startswith('T'):
                    upper[i] = min(upper[i], 10.0)
                if name.startswith('m'):
                    upper[i] = min(upper[i], 50.0)

        # prepare model function (native if available)
        try:
            model_func, is_fallback = make_model_func(model, force_fallback=args.force_fallback)
            func_ex = get_extrap_func(model_func)
            used_fallback = is_fallback
        except Exception as e:
            safe_print(f"[{now()}] ERROR selecting model function for {model}: {e}")
            if args.debug: traceback.print_exc()
            continue

        # RESUME mode: if --resume and a popt.npy exists, skip optimization and use saved popt
        resdir = os.path.join(outdir, f"model_{model}")
        os.makedirs(resdir, exist_ok=True)
        popt_path = os.path.join(resdir, "popt.npy")
        fit_summary_path = os.path.join(resdir, "fit_summary.txt")

        if args.resume and os.path.exists(popt_path):
            safe_print(f"[{now()}] Resuming for model {model}: loading existing popt from {popt_path}")
            try:
                popt = np.load(popt_path)
                # recompute model SFS and theta and ll
                try:
                    model_arr = np.asarray(func_ex(popt, data_fs.sample_sizes, pts_l))
                    if np.any(np.isnan(model_arr)) or np.any(model_arr < 0):
                        raise ValueError("native func_ex produced invalid model")
                except Exception:
                    # try fallback
                    safe_print(f"[{now()}] native model function failed to evaluate saved popt; trying fallback implementation")
                    try:
                        model_func_fb, _ = make_model_func(model, force_fallback=True)
                        func_ex_fb = get_extrap_func(model_func_fb)
                        model_arr = np.asarray(func_ex_fb(popt, data_fs.sample_sizes, pts_l))
                        used_fallback = True
                    except Exception as e:
                        safe_print(f"[{now()}] Could not evaluate popt with fallback either: {e}")
                        continue
                model_spec = dadi.Spectrum(model_arr)
                theta = dadi.Inference.optimal_sfs_scaling(model_spec, data_fs)
                ll = dadi.Inference.ll_multinom(model_spec * theta, data_fs)
                safe_print(f"[{now()}] Resumed model {model}: ll={ll:.3f}, theta={theta}")
            except Exception as e:
                safe_print(f"[{now()}] Failed to load/resume popt for {model}: {e}")
                if args.debug: traceback.print_exc()
                continue

        else:
            # Normal optimization path
            if model == 'IM':
                safe_print(f"[{now()}] IM: robust multi-candidate optimization (limited to {args.ncandidates} candidates)")
                prior_good = np.array([0.5, 0.8, 0.8, 0.25, 0.5, 0.5])
                candidates = [
                    np.array([0.5, 0.5, 0.5, 0.05, 0.001, 0.001]),
                    np.array([0.5, 5.0, 0.05, 0.2, 0.001, 0.001]),
                    np.array([0.5, 1.0, 1.0, 0.2, 0.1, 0.1]),
                    np.array([0.5, 2.0, 0.5, 1.0, 0.01, 0.01]),
                    np.array([0.5, 10.0, 0.2, 0.5, 0.001, 0.2]),
                ]
                priority = []
                if im_user_p0 is not None: priority.append(im_user_p0)
                if not any(np.allclose(prior_good, c) for c in candidates): priority.append(prior_good)
                master = []
                for item in priority + candidates + [p0]:
                    if not any(np.allclose(item, x) for x in master):
                        master.append(item)
                nc = max(1, int(args.ncandidates))
                if len(master) > nc:
                    safe_print(f"[{now()}] Trimming candidate list from {len(master)} to {nc} (keeping highest-priority starts first).")
                candidates = master[:nc]

                best = None; best_ll = -np.inf
                for idx, cand in enumerate(candidates):
                    safe_print(f"[{now()}] IM candidate {idx+1}/{len(candidates)} start={np.round(cand,6)}")
                    try:
                        res = optimize_two_stage(
                            data_fs, ns, func_ex, cand, lower, upper, pts_l,
                            n_restarts=max(1, min(args.nrestarts, 6)),
                            n_coarse=args.ncoarse if args.ncoarse is not None else 3,
                            coarse_maxiter=min(args.coarse_maxiter, 25),
                            refine_maxiter=min(args.maxiter, 200),
                            verbose=args.verbose, debug=args.debug
                        )
                        popt_c, model_spec_c, theta_c, ll_c = res
                        if theta_c is None or theta_c <= 0 or not np.isfinite(ll_c): continue
                        if ll_c < -5e6: continue
                        if ll_c > best_ll:
                            best_ll = ll_c; best = (popt_c, model_spec_c, theta_c, ll_c)
                            used_fallback = False
                            safe_print(f"[{now()}]  candidate {idx+1} new best (ll={ll_c:.3f})")
                    except Exception as e:
                        safe_print(f"[{now()}]  candidate {idx+1} failed: {e}")
                        if args.debug: traceback.print_exc()
                        continue

                if best is None and not args.force_fallback:
                    safe_print(f"[{now()}] IM: native impl failed — trying fallback impl with same {len(candidates)} candidates")
                    try:
                        model_func_fb, _ = make_model_func(model, force_fallback=True)
                        func_ex_fb = get_extrap_func(model_func_fb)
                        for idx, cand in enumerate(candidates):
                            safe_print(f"[{now()}] IM(fallback) candidate {idx+1}/{len(candidates)} start={np.round(cand,6)}")
                            try:
                                res = optimize_two_stage(
                                    data_fs, ns, func_ex_fb, cand, lower, upper, pts_l,
                                    n_restarts=max(1, min(args.nrestarts, 6)),
                                    n_coarse=args.ncoarse if args.ncoarse is not None else 3,
                                    coarse_maxiter=min(args.coarse_maxiter, 25),
                                    refine_maxiter=min(args.maxiter, 200),
                                    verbose=args.verbose, debug=args.debug
                                )
                                popt_c, model_spec_c, theta_c, ll_c = res
                                if theta_c is None or theta_c <= 0 or not np.isfinite(ll_c): continue
                                if ll_c < -5e6: continue
                                if ll_c > best_ll:
                                    best_ll = ll_c; best = (popt_c, model_spec_c, theta_c, ll_c)
                                    used_fallback = True
                                    safe_print(f"[{now()}]  fallback candidate {idx+1} new best (ll={ll_c:.3f})")
                            except Exception as e:
                                safe_print(f"[{now()}]  fallback candidate {idx+1} failed: {e}")
                                if args.debug: traceback.print_exc()
                                continue
                    except Exception as e:
                        safe_print(f"[{now()}] Could not instantiate fallback IM: {e}")
                        if args.debug: traceback.print_exc()

                if best is None:
                    safe_print(f"[{now()}] ERROR: All IM attempts failed (native + fallback). Try different starts or smaller pts for testing.")
                    continue

                popt, model_spec, theta, ll = best

            else:
                # non-IM standard model fit
                try:
                    popt, model_spec, theta, ll = optimize_two_stage(
                        data_fs, ns, func_ex, p0, lower, upper, pts_l,
                        n_restarts=args.nrestarts, n_coarse=args.ncoarse,
                        coarse_maxiter=args.coarse_maxiter, refine_maxiter=args.maxiter,
                        verbose=args.verbose, debug=args.debug
                    )
                    used_fallback = is_fallback
                except Exception as e:
                    safe_print(f"[{now()}] Model {model} optimization failed: {e}")
                    if args.debug: traceback.print_exc()
                    continue

        # post-fit checks and save
        obs_mat = np.asarray(data_fs)
        total_sites = int(np.rint(obs_mat.sum())) if obs_mat.size else 1
        if theta is None or theta <= 0 or not np.isfinite(ll):
            safe_print(f"[{now()}] Post-fit sanity failed for {model}: theta={theta}, ll={ll} -> skipping model")
            continue
        per_site_ll = ll / total_sites if total_sites > 0 else float('nan')
        exp_mat = np.asarray(model_spec * theta)
       
        np.save(popt_path, popt)
        safe_print(f"[{now()}] Saved popt -> {popt_path}")

        # decide whether to run bootstrap
        boots = np.empty((0, len(popt)))
        boot_out_path = os.path.join(resdir, "bootstrap_params.npy")
        if args.no_bootstrap:
            safe_print(f"[{now()}] Skipping bootstrap for {model} (--no-bootstrap)")
        else:
            if os.path.exists(boot_out_path) and not args.force_rerun:
                try:
                    boots = np.load(boot_out_path)
                    safe_print(f"[{now()}] Found existing bootstrap file and --force-rerun not set: loaded {boot_out_path}")
                except Exception:
                    safe_print(f"[{now()}] Could not load existing bootstrap file; will run bootstrap afresh")
            else:
                try:
                    model_prob = np.asarray(model_spec, dtype=float)
                    if np.any(np.isnan(model_prob)): raise ValueError('Model fs contains NaNs')
                    model_prob = model_prob / model_prob.sum()
                    safe_print(f"[{now()}] Running bootstrap n={args.nboot} (total_sites={total_sites}) with threads={args.threads}")
                    boots = parametric_bootstrap_parallel(
                        model_prob, total_sites, model, used_fallback,
                        popt, pts_l, lower, upper,
                        nboot=args.nboot, threads=args.threads, verbose=args.verbose, worker_maxiter=40
                    )
                    safe_print(f"[{now()}] Bootstrap finished for {model}")
                except Exception as e:
                    safe_print(f"[{now()}] Bootstrap failed for {model}: {e}")
                    if args.debug: traceback.print_exc()
                    boots = np.empty((0, len(popt)))

        # save bootstrap and summary
        np.save(boot_out_path, boots)
        if boots.size:
            np.savetxt(os.path.join(resdir, "bootstrap_params.csv"), boots, delimiter=",", fmt="%.6e", header=",".join(param_names[model]), comments='')
        else:
            with open(os.path.join(resdir, "bootstrap_params.csv"), 'w') as fh: fh.write(",".join(param_names[model]) + "\n")

        k = len(popt); aic = 2*k - 2*ll
        with open(fit_summary_path, 'w') as fh:
            fh.write(f"Model: {model}\n")
            fh.write(f"param_names: {param_names[model]}\n")
            fh.write(f"popt: {popt.tolist()}\n")
            fh.write(f"theta: {theta}\n")
            fh.write(f"loglikelihood: {ll}\n")
            fh.write(f"per_site_ll: {per_site_ll}\n")
            fh.write(f"AIC: {aic}\n")
            fh.write(f"total_sites: {total_sites}\n")
            fh.write(f"used_fallback: {used_fallback}\n")
            if boots.size:
                boots_clean = boots[~np.isnan(boots).any(axis=1)]
                fh.write("Bootstrap summary (per parameter):\nparam\tmean\tmedian\tstd\t2.5%\t97.5%\n")
                for i, name in enumerate(param_names[model]):
                    if boots_clean.shape[0] > 0:
                        vals = boots_clean[:, i]
                        fh.write(f"{name}\t{np.mean(vals):.6e}\t{np.median(vals):.6e}\t{np.std(vals):.6e}\t{np.percentile(vals,2.5):.6e}\t{np.percentile(vals,97.5):.6e}\n")
                    else:
                        fh.write(f"{name}\tNA\tNA\tNA\tNA\tNA\n")
            else:
                fh.write("No bootstrap results available.\n")

        # plots
        try:
            plot_obs_exp_res(obs_mat, exp_mat, os.path.join(resdir, f"{model}_obs_exp_res.png"), title=f"{model} Observed vs Expected (LL={ll:.2f})", pop_labels=tuple(args.poplabels))
            plot_param_violin(popt, boots, param_names[model], os.path.join(resdir, f"{model}_params_violin.png"))
            plot_param_histograms(popt, boots, param_names[model], os.path.join(resdir, f"{model}_params_histograms.png"))
            plot_pairwise_scatter(boots, param_names[model], os.path.join(resdir, f"{model}_params_pairwise.png"))
            plot_bootstrap_corr_heatmap(boots, param_names[model], os.path.join(resdir, f"{model}_params_corr_heatmap.png"))
        except Exception as e:
            safe_print(f"[{now()}] Warning: plotting failed for {model}: {e}")
            if args.debug: traceback.print_exc()

        safe_print(f"[{now()}] Saved results to {resdir} (ll={ll:.3f}, per_site_ll={per_site_ll:.6e})")
        safe_print(f"[{now()}] Finished {model} in {(time.time()-start_t):.1f}s")
        results_local[model] = {'popt':popt, 'll':ll, 'aic':aic, 'resdir':resdir}

    # combine existing model_* folders + local results for summary/plot
    write_combined_model_summary(outdir, results_local)
    safe_print(f"[{now()}] DONE. Results in: {outdir}")

if __name__ == '__main__':
    main()
