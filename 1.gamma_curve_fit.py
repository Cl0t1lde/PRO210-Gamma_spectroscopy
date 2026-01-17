import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# === Function to load a spectrum from a single file ===
def load_spectrum(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    start_row = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Channel Data"):
            start_row = i + 2
            break

    if start_row is None:
        raise ValueError(f"Could not find 'Channel Data' in {filename}")

    df = pd.read_csv(
        filename,
        skiprows=start_row,
        usecols=[0, 2],
        names=["channel", "counts"]
    )

    df["channel"] = pd.to_numeric(df["channel"], errors="coerce")
    df["counts"] = pd.to_numeric(df["counts"], errors="coerce")
    df = df.dropna()

    return df


SAMPLE_LABELS = {
    "Tungsten": "Tungsten",
    "Banana": "Banana",
    "Potassium": "Potassium (K-40)",
    "Uranium": "Uranium",
    "Rocks": "Rock sample",
    "Background": "Background"
}


def makegraph(filename):
    df = load_spectrum(filename)

    plt.figure(figsize=(10,5))
    plt.plot(df["channel"], df["counts"], drawstyle="steps-mid")
    plt.xlabel("Channel")
    plt.ylabel("Counts")
    plt.title(f"Full Gamma Spectrum: {filename}")
    plt.grid(alpha=0.4)
    plt.show()

    return df


def peak_analyser(t_k, filename, df_bg):

    # >>> NEW / MODIFIED <<<
    is_background = "background" in filename.lower()

    sample_label = "Sample"
    for key, label in SAMPLE_LABELS.items():
        if key.lower() in filename.lower():
            sample_label = label
            break

    t_bg = 4838
    df_p = load_spectrum(filename)

    scale = t_k / t_bg
    print(f"\nFile: {filename}")
    print(f"Background scale factor: {scale:.3f}")

    # >>> NEW / MODIFIED <<<
    if not is_background:
        df_bg = df_bg.copy()
        df_bg["counts_scaled"] = df_bg["counts"] * scale
        df = df_p.merge(df_bg[["channel", "counts_scaled"]], on="channel", how="inner")
        df["counts_sub"] = df["counts"] - df["counts_scaled"]
    else:
        df = df_p.copy()
        df["counts_sub"] = df["counts"]
        df["counts_scaled"] = 0.0

    # --- Plot ---
    plt.figure(figsize=(10,5))
    plt.plot(df["channel"], df["counts"], drawstyle="steps-mid", label=sample_label)
    if not is_background:
        plt.plot(df["channel"], df["counts_sub"], drawstyle="steps-mid", label="Background subtracted")
    plt.xlabel("Channel")
    plt.ylabel("Counts")
    plt.legend()
    plt.grid(alpha=0.4)
    plt.show()

    x = df["channel"].to_numpy()
    y_sub = df["counts_sub"].to_numpy()

    peak_windows_dict = {
        "PHA PreAmp Tungsten 4196sec; 14-1-2026.csv": [(700, 740), (524, 563), (348, 411), (243, 285), (145, 175), (90, 132),(16, 44)],
        "PreAmp Background measurement 3671sec;14-1-2026.csv": [(480, 540), (240, 300)],
        "PHA PreAmp Banana 7422 sec; 8-1-2026.csv": [(490, 530)],
        "PHA PreAmp Potassium 4835 sec; 8-1-2026.csv": [(485, 540)],
        "PHA PreAmp Rocks 7200sec; 15-1-2026.csv": [(495, 530), (250, 286), (152, 179), (105, 130)],
        "PHA PreAmp Vase Uranium 3638sec; 14-1-2026.csv": [(81, 97), (37, 52), (25, 37), (499, 530), (62, 78.5)],
    }

    peak_windows = peak_windows_dict.get(filename, [])
    all_results = []

    for low, high in peak_windows:

        mask = (x > low) & (x < high)
        x_peak = x[mask]
        y_peak = y_sub[mask]

        if len(x_peak) == 0:
            continue

        print(f"\nAnalyzing window {low}-{high}")

        # --- Baseline estimation ---
        side_width = 5
        left_bg = (x >= low-side_width) & (x < low)
        right_bg = (x > high) & (x <= high+side_width)
        baseline_mask = left_bg | right_bg

        baseline = np.mean(y_sub[baseline_mask]) if baseline_mask.any() else 0.0

        # --- Baseline subtraction (Compton removal) ---
        y_corr = y_peak - baseline

        # --- Statistics ---
        N_net = y_corr.sum()
        N_gross = y_peak.sum()

        # >>> NEW / MODIFIED <<<
        if not is_background:
            y_bg_scaled = df["counts_scaled"][mask].to_numpy()
            N_bg = y_bg_scaled.sum()
            sigma_net = np.sqrt(abs(N_gross) + N_bg * scale**2)
        else:
            N_bg = 0.0
            sigma_net = np.sqrt(abs(N_gross))

        significance = N_net / sigma_net if sigma_net > 0 else 0.0

        print(f"N_net = {N_net:.1f}")
        print(f"Significance = {significance:.2f} σ")

        # --- Gaussian fit ---
        def gaussian(x, A, mu, sigma, C):
            return A * np.exp(-(x - mu)**2 / (2*sigma**2)) + C

        if y_corr.max() <= 0:
            print("No positive peak, skipping fit")
            continue

        p0 = [y_corr.max(), x_peak[np.argmax(y_corr)], (high-low)/4, 0]

        try:
            params, cov = curve_fit(gaussian, x_peak, y_corr, p0=p0, maxfev=10000)
            A, mu, sigma, C = params
            
            perr = np.sqrt(np.diag(cov))
            FWHM = 2.355 * sigma  # Full Width at Half Maximum
            residuals = y_corr - gaussian(x_peak, *params)
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((y_corr - np.mean(y_corr))**2)
            r_squared = 1 - (ss_res / ss_tot)
            chi2_red = np.sum((residuals**2)/np.maximum(y_corr,1)) / (len(y_corr)-len(params))

            
            x_fit = np.linspace(low, high, 500)
            y_fit = gaussian(x_fit, *params)

            plt.figure(figsize=(8,4))
            plt.scatter(x_peak, y_peak, s=10, label="Data")
            plt.plot(x_fit, y_fit + baseline, 'r', label="Gaussian fit")
            plt.hlines(baseline, low, high, colors="orange", linestyles="--", label="Baseline")
            plt.xlabel("Channel")
            plt.ylabel("Counts")
            plt.legend()
            plt.grid(alpha=0.4)
            plt.show()
            
            print(f"μ     = {mu:.2f} ± {perr[1]:.2f}")
            print(f"σ     = {sigma:.2f}")
            print(f"FWHM  = {FWHM:.2f}")
            print(f"R²    = {r_squared:.4f}")
            print(f"χ²red = {chi2_red:.2f}")

        except RuntimeError:
            print("Fit failed")

    return all_results


# ==========================================================
# Run analysis
# ==========================================================
graph_bg_1 = makegraph("PHA PreAmp Background 4838 sec; 7-1-2026.csv")

peak_analyser(3638, "PHA PreAmp Vase Uranium 3638sec; 14-1-2026.csv", graph_bg_1)
peak_analyser(7200, "PHA PreAmp Rocks 7200sec; 15-1-2026.csv", graph_bg_1)
peak_analyser(4196, "PHA PreAmp Tungsten 4196sec; 14-1-2026.csv", graph_bg_1)
peak_analyser(3671, "PreAmp Background measurement 3671sec;14-1-2026.csv", graph_bg_1)
peak_analyser(7422, "PHA PreAmp Banana 7422 sec; 8-1-2026.csv", graph_bg_1)
peak_analyser(4835, "PHA PreAmp Potassium 4835 sec; 8-1-2026.csv", graph_bg_1)
