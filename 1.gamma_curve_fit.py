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

    is_background = "background" in filename.lower()

    # Determine sample label
    sample_label = "Sample"
    for key, label in SAMPLE_LABELS.items():
        if key.lower() in filename.lower():
            sample_label = label
            break

    # Reference time for background scaling
    t_bg = 4838 if df_bg is graph_bg_1 else 3671

    # Load sample spectrum
    df_p = load_spectrum(filename)
    scale = t_k / t_bg
    print(f"\nFile: {filename}")
    print(f"Background scale factor: {scale:.3f}")

    # Merge with scaled background
    if not is_background:
        df_bg = df_bg.copy()
        df_bg["counts_scaled"] = df_bg["counts"] * scale
        df = df_p.merge(df_bg[["channel", "counts_scaled"]], on="channel", how="inner")
        df["counts_sub"] = df["counts"] - df["counts_scaled"]
    else:
        df = df_p.copy()
        df["counts_sub"] = df["counts"]
        df["counts_scaled"] = 0.0

    x = df["channel"].to_numpy()
    y_sub = df["counts_sub"].to_numpy()
    
    plt.figure(figsize=(10,5))
    plt.plot(df["channel"], df["counts"], drawstyle="steps-mid", label="Sample")
    
    if not is_background:
        plt.plot(df["channel"], df["counts_scaled"], drawstyle="steps-mid", label="Background (scaled)")
        plt.plot(df["channel"], df["counts_sub"], drawstyle="steps-mid", label="Subtracted")
    plt.xlabel("Channel")
    plt.ylabel("Counts")
    plt.legend()
    plt.grid(alpha=0.4)
    num_ticks = 20
    plt.xticks(np.linspace(df["channel"].min(), df["channel"].max(), num_ticks))
    plt.show()
    
    
    plt.figure(figsize=(10,5))
    plt.plot(df["channel"], df["counts"], linestyle="None", marker=".", markersize=4, label="Sample")
    if not is_background:
        plt.plot(df["channel"], df["counts_scaled"], linestyle="None", marker=".", markersize=4, label="Background (scaled)")
        plt.plot(df["channel"], df["counts_sub"], linestyle="None", marker=".", markersize=4, label="Subtracted")
    
    plt.yscale("log")                # logarithmic y-axis
    plt.xlim(0, 600)                 # stop at channel 600
    plt.xlabel("Channel")
    plt.ylabel("Counts (log scale)")
    plt.legend()
    plt.grid(alpha=0.4, which="both")

    # Optional: nicer x ticks
    num_ticks = 15
    plt.xticks(np.linspace(0, 600, num_ticks))

    plt.show()


    # Nominal peak windows
    peak_windows_dict = {
        "PHA PreAmp Tungsten 4196sec; 14-1-2026.csv": [(700, 740), (524, 563), (348, 411), (243, 285), (145, 175), (90, 132),(16, 44)],
        "PreAmp Background measurement 3671sec;14-1-2026.csv": [(480, 540), (240, 300)],
        "PHA PreAmp Banana 7422 sec; 8-1-2026.csv": [(480, 540)],
        "PHA PreAmp Potassium 4835 sec; 8-1-2026.csv": [(485, 540)],
        "PHA PreAmp Rocks 7200sec; 15-1-2026.csv": [(495, 530), (250, 286), (152, 179), (105, 130)],
        "PHA PreAmp Vase Uranium 3638sec; 14-1-2026.csv": [(81, 97), (37, 52), (25, 37), (499, 530), (62, 78.5)],
    }

    peak_windows = peak_windows_dict.get(filename, [])
    all_results = []

    for nominal_low, nominal_high in peak_windows:

        # Nominal width
        nominal_width = int(nominal_high - nominal_low)
        
        # Candidate shifts (around center)
        shifts = [-4, -2, -1, 0, 1, 2, 4]  # integer channels

        # Candidate width variations: ±10% of nominal width
        width_variants = [nominal_width, int(nominal_width * 1.1)]

        candidate_windows = []

        for w in width_variants:
            for s in shifts:
                center = (nominal_low + nominal_high) // 2 + s
                low = center - w // 2
                high = center + w // 2
                # Avoid negative or zero-width windows
                if low < 0 or high <= low:
                    continue
                candidate_windows.append((low, high))

        best_r2 = -np.inf
        best_window = None
        candidate_results = []

        for low, high in candidate_windows:
            mask = (x > low) & (x < high)
            x_peak = x[mask]
            y_peak = y_sub[mask]
            if len(x_peak) == 0:
                continue

            # Baseline estimation
            side_width = 4
            baseline_mask = ((x >= low-side_width) & (x < low)) | ((x > high) & (x <= high+side_width))
            baseline = np.min(y_sub[baseline_mask]) if baseline_mask.any() else 0.0
            baseline = max(baseline, 0.0)
            y_corr = y_peak - baseline
            y_corr = np.clip(y_corr, 0, None)

            if y_corr.max() <= 0:
                continue

            # Gaussian fit
            # Gaussian fit
            def gaussian(x, A, mu, sigma, C):
                return A * np.exp(-(x-mu)**2/(2*sigma**2)) + C

            p0 = [y_corr.max(), x_peak[np.argmax(y_corr)], (high-low)/4, 0]

            try:
                # Fit and get covariance matrix
                params, cov = curve_fit(gaussian, x_peak, y_corr, p0=p0, maxfev=10000)
                y_fit = gaussian(x_peak, *params)

                residuals = y_corr - y_fit
                ss_res = np.sum(residuals**2)
                ss_tot = np.sum((y_corr - np.mean(y_corr))**2)
                r2 = 1 - ss_res/ss_tot

                # Uncertainties from covariance matrix
                sigma_params = np.sqrt(np.diag(cov))  # [σ_A, σ_μ, σ_σ, σ_C]
                sigma_A, sigma_mu, sigma_sigma, sigma_C = sigma_params

                # FWHM and its uncertainty
                FWHM = 2.355 * params[2]
                sigma_FWHM = 2.355 * sigma_sigma

                # Save candidate result
                candidate_results.append({"low": low, "high": high, "r2": r2})

                # Best R² check
                if r2 > best_r2:
                    best_r2 = r2
                    best_window = (low, high)
                    best_params = params
                    best_cov = cov
                    best_sigma_params = sigma_params
                    best_FWHM = FWHM
                    best_sigma_FWHM = sigma_FWHM
                    best_baseline = baseline
                    best_x_peak = x_peak
                    best_y_peak = y_peak
                    best_y_corr = y_corr

            except RuntimeError:
                continue


        # Print all tried ROIs for this peak
        print(f"\nCandidate ROIs for peak near {int((nominal_low+nominal_high)/2)} channels:")
        for cr in candidate_results:
            print(f"  {cr['low']:.1f} – {cr['high']:.1f}, R² = {cr['r2']:.4f}")

        if best_window is None:
            continue

        # Compute statistics for best ROI
        low, high = best_window
        A, mu, sigma, C = best_params
        sigma_A, sigma_mu, sigma_sigma, sigma_C = best_sigma_params
        FWHM = best_FWHM
        sigma_FWHM = best_sigma_FWHM
        N_gross = best_y_peak.sum()
        N_net = best_y_corr.sum()
        N_bg = df["counts_scaled"][(x>low)&(x<high)].sum() if not is_background else 0.0
        s_b = N_net/N_bg if N_bg>0 else np.nan
        sigma_net = np.sqrt(N_gross + N_bg)
        sigma_bg = np.sqrt(N_bg) if N_bg>0 else 0.0
        if N_net > 0 and N_bg > 0:
            sigma_s_b = s_b * np.sqrt( (sigma_net / N_net)**2 + (sigma_bg / N_bg)**2 )
        else:
            sigma_s_b = np.nan
        significance = N_net / sigma_net if sigma_net>0 else 0.0

        print(f"\nBest ROI for peak near {int((nominal_low+nominal_high)/2)} channels: {low:.1f} – {high:.1f}")
        print(f"N_net = {N_net:.2f} ± {sigma_net:.2f}")
        print(f"N_bg  = {N_bg:.2f} ± {sigma_bg:.2f}")
        print(f"S/B ratio = {s_b:.2f} ± {sigma_s_b:.2f}")
        print(f"Significance = {significance:.2f} σ")
        print(f"μ = {mu:.2f} ± {sigma_mu:.2f}")
        print(f"σ = {sigma:.2f} ± {sigma_sigma:.2f}")
        print(f"FWHM = {FWHM:.2f} ± {sigma_FWHM:.2f}")
        print(f"R² = {best_r2:.4f}")
        
        # Plot everything
        # Plot everything
        plot_margin = 50
        x_min_plot = max(0, low - plot_margin)  # don't go below channel 0
        x_max_plot = min(x.max(), high + plot_margin)
        plot_mask = (x >= x_min_plot) & (x <= x_max_plot)

        # Extend Gaussian for plotting: mu ± 3 sigma
        x_fit_extended = np.linspace(mu - 3*sigma, mu + 3*sigma, 500)
        y_fit_extended = gaussian(x_fit_extended, *best_params) + best_baseline  # include baseline offset

        plt.figure(figsize=(10,5))
        plt.plot(x[plot_mask], df["counts"][plot_mask], drawstyle="steps-mid", label="Raw data")
        if not is_background:
            plt.plot(x[plot_mask], df["counts_scaled"][plot_mask], drawstyle="steps-mid", label="Scaled background")
            plt.plot(x[plot_mask], df["counts_sub"][plot_mask], drawstyle="steps-mid", label="Background removed")
        plt.plot(x_fit_extended, y_fit_extended, 'r', label="Gaussian fit")  # now full Gaussian
        plt.hlines(best_baseline, x[plot_mask].min(), x[plot_mask].max(), linestyles="--", color="gray", label="Baseline")
        plt.xlabel("Channel")
        plt.ylabel("Counts")
        plt.title(f"Peak {int(low)}-{int(high)}")
        plt.legend()
         # data stays positive, Gaussian may go below baseline visually
        plt.grid(alpha=0.4)
        plt.show()




        # Store results
        all_results.append({
            "window": best_window,
            "N_net": N_net,
            "N_bg": N_bg,
            "N_gross": N_gross,
            "s_b": s_b,
            "significance": significance,
            "mu": mu,
            "sigma": sigma,
            "FWHM": FWHM,
            "R2": best_r2
        })

    return all_results


# ==========================================================
# Run analysis
# ==========================================================
graph_bg_1 = makegraph("PHA PreAmp Background 4838 sec; 7-1-2026.csv")
graph_bg_2 = makegraph("PreAmp Background measurement 3671sec;14-1-2026.csv")


peak_analyser(3638, "PHA PreAmp Vase Uranium 3638sec; 14-1-2026.csv", graph_bg_2)
peak_analyser(7200, "PHA PreAmp Rocks 7200sec; 15-1-2026.csv", graph_bg_2)

peak_analyser(3671, "PreAmp Background measurement 3671sec;14-1-2026.csv", graph_bg_2)


peak_analyser(7422, "PHA PreAmp Banana 7422 sec; 8-1-2026.csv", graph_bg_1)
peak_analyser(4835, "PHA PreAmp Potassium 4835 sec; 8-1-2026.csv", graph_bg_1)
