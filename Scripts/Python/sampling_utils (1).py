import os
from copy import deepcopy
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import ks_2samp, gaussian_kde
from cobra.io import read_sbml_model
from cobra.sampling.optgp import OptGPSampler

try:
    from cobra.core.gene import GPR
    def _gpr_copy(self):
        
        """ Fixing GPR copy bug"""
        
        new = GPR()
        try:
            new._genes = self._genes.copy()
        except Exception:
            new._genes = self._genes
        new.body = deepcopy(self.body)
        return new
    GPR.__copy__ = _gpr_copy
except Exception:
    pass


def run_flux_sampling(model_wt, model_ros, candidate_reactions, n_samples=5000, thinning=100,
                       chains=1, output_dir=Path("sampling_results")):

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    """
    Takes as input two pre-loaded metabolic models and a set of candidate reactions that want to be inspected.
    Performs flux sampling onto the two models and creates a dataframe that appends all results with a label column (wt vs ros)
    Default parameters are 10,000 samples and thinning factor of k=1000, which can be lowered to reduce computational burden
    Creates a directory called ./sampling_results to store outputs.
    Returns two dataframes: all reactions results, and only the ones corresponding to the given candidate reactions.
    Also calculates Gelman-Rubin statistic which gives a measure of the convergence of the sampled chains.
    """

    all_samples = []
    candidate_samples = []

    models = [("wt", model_wt), ("ros", model_ros)]

    for label, model_orig in models:

        chain_samples = []  # store per-chain samples for R-hat

        for chain in range(chains):

            model = deepcopy(model_orig)

            try:
                sampler = OptGPSampler(model, thinning=thinning)
                samples = sampler.sample(n=n_samples)
                samples = samples.iloc[-n_samples:]  # discard initial samples (warmup-like behavior)
            except Exception as e:
                print(f"Sampling failed for {label}, chain {chain}: {e}")
                continue

            chain_samples.append(samples)

            # all reactions
            df_all = samples.reset_index(drop=True).copy()
            df_all["label"] = label
            df_all["chain"] = chain
            all_samples.append(df_all)

            # only candidate reactions
            df_cand = samples.loc[:, candidate_reactions].reset_index(drop=True).copy()
            df_cand["label"] = label
            df_cand["chain"] = chain
            candidate_samples.append(df_cand)

        valid_rxns = [
            rxn for rxn in candidate_reactions
            if all(rxn in df.columns for df in chain_samples)
        ]

        if len(chain_samples) < 2 or len(valid_rxns) == 0:
            print(f"Gelman-Rubin (R̂) for {label}: not computed")
        else:
            rhat_values = []
            for rxn in valid_rxns:
                chains_rxn = [df[rxn].values for df in chain_samples]
                n = min(len(c) for c in chains_rxn)
                chains_rxn = [c[:n] for c in chains_rxn]

                chain_means = np.array([c.mean() for c in chains_rxn])
                chain_vars = np.array([c.var(ddof=1) for c in chains_rxn])

                B = n * chain_means.var(ddof=1)
                W = chain_vars.mean()

                if W == 0:
                    rhat = np.nan
                else:
                    var_hat = (n - 1) / n * W + B / n
                    rhat = np.sqrt(var_hat / W)

                rhat_values.append(rhat)

            max_rhat = np.nanmax(rhat_values)
            print(f"Gelman-Rubin (R̂) for {label}: {max_rhat:.3f}  (max over {len(valid_rxns)} candidate reactions)") # prints in the console the mean value of Rhat for both models. 

    df_all = pd.concat(all_samples, ignore_index=True)
    df_candidates = pd.concat(candidate_samples, ignore_index=True)

    df_all.to_csv(output_dir / "all_reactions_flux_samples.csv", index=False)
    df_candidates.to_csv(output_dir / "candidate_reactions_flux_samples.csv", index=False)

    return df_all, df_candidates


def plot_pca_samples(df_all, labels=["wt", "ros"]):
    """
    Simple 2D PCA scatter of samples in df_all. Expects a 'label' column with values in `labels`.
    Saves 'pca_samples.png' in the working directory.
    """

    # numeric reaction columns (drop common metadata)
    drop = {"label", "chain", "model"}
    cols = [c for c in df_all.columns if c not in drop and np.issubdtype(df_all[c].dtype, np.number)]
    X = df_all[cols].values
    pca = PCA(n_components=2)
    Z = pca.fit_transform(X)
    df_plot = pd.DataFrame(Z, columns=("PC1", "PC2"))
    df_plot["label"] = df_all["label"].values

    plt.figure(figsize=(6,5))
    for lbl in labels:
        sel = df_plot[df_plot["label"] == lbl]
        plt.scatter(sel["PC1"], sel["PC2"], s=10, label=lbl, alpha=0.6)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.title("PCA of flux samples")
    plt.tight_layout()
    plt.savefig("pca_samples.png", dpi=150)
    plt.close()

def run_ks_test(df_candidates, labels=["wt", "ros"]):

    """
    Runs 2-sample KS test for every reaction in df_candidates between labels[0] and labels[1].
    Returns a dataframe (reaction, pvalue) and saves 'ks_test_results.xlsx'.
    Expects 'label' column in df_candidates.
    """
    drop = {"label", "chain", "model"}
    rxns = [c for c in df_candidates.columns if c not in drop]
    rows = []
    a_lbl, b_lbl = labels[0], labels[1]
    for rxn in rxns:
        a = df_candidates.loc[df_candidates["label"] == a_lbl, rxn].dropna().values
        b = df_candidates.loc[df_candidates["label"] == b_lbl, rxn].dropna().values
        if a.size == 0 or b.size == 0:
            p = np.nan
        else:
            p = float(ks_2samp(a, b).pvalue)
        rows.append({"reaction": rxn, "pvalue": p})
    res = pd.DataFrame(rows)
    res.to_excel("ks_test_results.xlsx", index=False)
    return res

def plot_distributions(df_candidates, labels=["wt", "ros"]):


    """
    For each reaction in df_candidates plots KDEs for the two labels and saves PNGs.
    
    """
    out_dir = Path("distribution_plots")
    out_dir.mkdir(exist_ok=True)
    drop = {"label", "chain", "model"}
    rxns = [c for c in df_candidates.columns if c not in drop]

    for rxn in rxns:
        a = df_candidates.loc[df_candidates["label"] == labels[0], rxn].dropna().values
        b = df_candidates.loc[df_candidates["label"] == labels[1], rxn].dropna().values
        if a.size < 2 and b.size < 2:
            continue

        plt.figure(figsize=(5,4))
        plotted = False
        if a.size >= 2:
            kde_a = gaussian_kde(a)
            xs = np.linspace(min(a.min(), b.min() if b.size>=1 else a.min()),
                             max(a.max(), b.max() if b.size>=1 else a.max()), 200)
            plt.plot(xs, kde_a(xs), label=labels[0], lw=1.2)
            plotted = True
        if b.size >= 2:
            kde_b = gaussian_kde(b)
            if not plotted:
                xs = np.linspace(b.min(), b.max(), 200)
            plt.plot(xs, kde_b(xs), label=labels[1], lw=1.2)
        plt.title(rxn)
        plt.xlabel("flux")
        plt.legend()
        plt.tight_layout()
        fname = f"dist_{rxn.replace(os.sep, '_').replace(' ', '_')}.png"
        plt.savefig(out_dir / fname, dpi=150)
        plt.close()


