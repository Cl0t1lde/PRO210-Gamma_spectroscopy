import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os

# ==========================================================
# Plot saving utilities
# ==========================================================
OUTPUT_DIR = "plots_2"  # base directory where all plots will be saved

def save_figure(folder, filename):
    """
    Save the current matplotlib figure to the specified folder with high resolution.
    Creates the folder if it doesn't exist.
    """
    os.makedirs(folder, exist_ok=True)
    path = os.path.join(folder, filename)
    plt.savefig(path, dpi=300, bbox_inches="tight")
    plt.close()  # close figure to avoid overlap

# ==========================================================
# Load spectrum
# ==========================================================
def load_spectrum(filename):
    """
    Load a gamma spectrum CSV file.
    - Automatically finds the line where channel data starts ("Channel Data")
    - Reads channel numbers and counts
    - Converts to numeric and drops invalid rows
    """
    with open(filename, "r") as f:
        lines = f.readlines()

    start_row = None
    for i, line in enumerate(lines):
        if line.strip().startswith("Channel Data"):
            start_row = i + 2  # data starts two lines after header
            break

    if start_row is None:
        raise ValueError(f"Could not find 'Channel Data' in {filename}")

    df = pd.read_csv(
        filename,
        skiprows=start_row,
        usecols=[0, 2],
        names=["channel", "counts"]
    )

    # Ensure data is numeric
    df["channel"] = pd.to_numeric(df["channel"], errors="coerce")
    df["counts"] = pd.to_numeric(df["counts"], errors="coerce")
    df = df.dropna()  # remove any invalid rows

    return df

# ==========================================================
# Sample labels
# ==========================================================
SAMPLE_LABELS = {
    "Tungsten": "Tungsten",
    "Banana": "Banana",
    "Potassium": "Potassium",
    "Uranium": "Uranium",
    "Rocks": "Rocks",
    "Background": "Background"
}

# ==========================================================
# Overview graph
# ==========================================================
def makegraph(filename):
    """
    Create an overview plot of the full spectrum for a given file.
    - Saves plot in the OUTPUT_DIR/Overview folder
    """
    df = load_spectrum(filename)

    plt.figure(figsize=(10,5))
    plt.plot(df["channel"], df["counts"], drawstyle="steps-mid")
    plt.xlabel("Channel")
    plt.ylabel("Counts")
    plt.title(f"Full Gamma Spectrum: {filename}")
    plt.grid(alpha=0.4)

    # Save figure
    save_figure(
        os.path.join(OUTPUT_DIR, "Overview"),
        f"full_spectrum_{filename.replace('.csv','')}.png"
    )

    return df

# ==========================================================
# Peak analyser
# ==========================================================
def peak_analyser(t_k, filename, df_bg):
    """
    Analyse peaks in a gamma spectrum:
    - Scales background according to acquisition time
    - Subtracts background
    - Generates linear and logarithmic plots
    - Fits Gaussian functions to peaks in defined ROIs
    - Saves all plots in structured folders
    """
    # Determine if the file is a background measurement
    is_background = "background" in filename.lower()

    # Determine sample label for plot naming
    sample_label = "Sample"
    for key, label in SAMPLE_LABELS.items():
        if key.lower() in filename.lower():
            sample_label = label
            break

    # Folder to save plots for this sample
    sample_folder = os.path.join(OUTPUT_DIR, sample_label)

    # Determine reference background acquisition time
    t_bg = 4838 if df_bg is graph_bg_1 else 3671

    # Load spectrum
    df_p = load_spectrum(filename)
    scale = t_k / t_bg  # scaling factor to match measurement time
    print(f"\nFile: {filename}")
    print(f"Background scale factor: {scale:.3f}")

    # Merge sample spectrum with scaled background
    if not is_background:
        df_bg = df_bg.copy()
        df_bg["counts_scaled"] = df_bg["counts"] * scale  # rescale background
        df = df_p.merge(df_bg[["channel", "counts_scaled"]], on="channel", how="inner")
        df["counts_sub"] = df["counts"] - df["counts_scaled"]  # background-subtracted spectrum
    else:
        # For pure background, no subtraction needed
        df = df_p.copy()
        df["counts_sub"] = df["counts"]
        df["counts_scaled"] = 0.0

    x = df["channel"].to_numpy()
    y_sub = df["counts_sub"].to_numpy()

    # ------------------------------------------------------
    # Linear spectrum plot
    # ------------------------------------------------------
    plt.figure(figsize=(10,5))
    plt.plot(df["channel"], df["counts"], drawstyle="steps-mid", label="Sample")
    plt.plot(df["channel"], df["counts_scaled"], drawstyle="steps-mid", label="Background (scaled)")
    plt.plot(df["channel"], df["counts_sub"], drawstyle="steps-mid", label="Subtracted")
    plt.xlabel("Channel")
    plt.ylabel("Counts")
    plt.legend()
    plt.grid(alpha=0.4)
    num_ticks = 20
    plt.xticks(np.linspace(df["channel"].min(), df["channel"].max(), num_ticks))
    
    save_figure(
        sample_folder,
        f"{filename.replace('.csv','')}_linear_spectrum.png"
    )

    # ------------------------------------------------------
    # Logarithmic spectrum plot
    # ------------------------------------------------------
    plt.figure(figsize=(10,5))
    plt.plot(df["channel"], df["counts"], ".", markersize=4, label="Sample")
    plt.plot(df["channel"], df["counts_scaled"], ".", markersize=4, label="Background")
    plt.plot(df["channel"], df["counts_sub"], ".", markersize=4, label="Subtracted")
    plt.yscale("log")
    plt.xlim(0, 600)
    plt.xlabel("Channel")
    plt.ylabel("Counts (log)")
    plt.legend()
    plt.grid(alpha=0.4, which="both")
    step = 30 
    plt.xticks(np.arange(0, 601, step))

    save_figure(
        sample_folder,
        f"{filename.replace('.csv','')}_log_spectrum.png"
    )

    # ======================================================
    # Peak windows (pre-defined ROIs in channel space)
    # ======================================================
    peak_windows_dict = {
        "PHA PreAmp Tungsten 4196sec; 14-1-2026.csv": [(700,740),(524,563),(348,411),(243,285),(145,175),(90,132),(16,44)],
        "PreAmp Background measurement 3671sec;14-1-2026.csv": [(480,540),(240,300)],
        "PHA PreAmp Banana 7422 sec; 8-1-2026.csv": [(480,540)],
        "PHA PreAmp Potassium 4835 sec; 8-1-2026.csv": [(485,540)],
        "PHA PreAmp Rocks 7200sec; 15-1-2026.csv": [(495,530),(250,286),(152,179),(105,130)],
        "PHA PreAmp Vase Uranium 3638sec; 14-1-2026.csv": [(81,97),(37,52),(25,37),(499,530),(62,78.5)]
    }

    peak_windows = peak_windows_dict.get(filename, [])
    all_results = []

    # ------------------------------------------------------
    # Iterate over each ROI
    # ------------------------------------------------------
    for nominal_low, nominal_high in peak_windows:
        nominal_width = int(nominal_high - nominal_low)
        shifts = [-4, -2, -1, 0, 1, 2, 4]  # candidate shifts to find best ROI
        width_variants = [nominal_width, int(nominal_width * 1.1)]  # candidate widths ±10%

        candidate_windows = []
        for w in width_variants:
            for s in shifts:
                center = (nominal_low + nominal_high) // 2 + s
                low = center - w // 2
                high = center + w // 2
                if low >= 0 and high > low:
                    candidate_windows.append((low, high))

        best_r2 = -np.inf

        # ------------------------------------------------------
        # Iterate candidate windows, estimate baseline, fit Gaussian
        # ------------------------------------------------------
        for low, high in candidate_windows:
            mask = (x > low) & (x < high)
            x_peak = x[mask]
            y_peak = y_sub[mask]
            if len(x_peak) == 0:
                continue

            # Baseline estimation using edges of the ROI (Compton continuum)
            side_width = 4
            baseline_mask = ((x >= low-side_width)&(x<low)) | ((x>high)&(x<=high+side_width))
            baseline = np.min(y_sub[baseline_mask]) if baseline_mask.any() else 0.0
            baseline = max(baseline, 0.0)

            # Subtract baseline, clip negative values
            y_corr = np.clip(y_peak - baseline, 0, None)
            if y_corr.max() <= 0:
                continue

            # Gaussian function: A*exp(-(x-mu)^2/(2*sigma^2)) + C
            def gaussian(x, A, mu, sigma, C):
                return A*np.exp(-(x-mu)**2/(2*sigma**2)) + C

            # Initial guess
            p0 = [y_corr.max(), x_peak[np.argmax(y_corr)], (high-low)/4, 0]

            # Fit Gaussian to peak
            try:
                params, _ = curve_fit(gaussian, x_peak, y_corr, p0=p0, maxfev=10000)
                y_fit = gaussian(x_peak, *params)
                # Goodness of fit (R²)
                r2 = 1 - np.sum((y_corr-y_fit)**2) / np.sum((y_corr-np.mean(y_corr))**2)

                if r2 > best_r2:
                    best_r2 = r2
                    best_params = params
                    best_baseline = baseline
                    best_window = (low, high)

            except RuntimeError:
                continue  # skip if fit fails

        # Skip if no valid fit found
        if best_r2 < 0:
            continue

        # ------------------------------------------------------
        # Plot final peak with Gaussian overlay
        # ------------------------------------------------------
        low, high = best_window
        A, mu, sigma, C = best_params
        plot_margin = 50
        mask_plot = (x >= max(0,low-plot_margin)) & (x <= min(x.max(),high+plot_margin))
        x_fit = np.linspace(mu-3*sigma, mu+3*sigma, 500)
        y_fit = gaussian(x_fit, *best_params) + best_baseline

        plt.figure(figsize=(10,5))
        plt.plot(x[mask_plot], df["counts"][mask_plot], drawstyle="steps-mid", label="Raw")
        plt.plot(x[mask_plot], df["counts_scaled"][mask_plot], drawstyle="steps-mid", label="Background")
        plt.plot(x[mask_plot], df["counts_sub"][mask_plot], drawstyle="steps-mid", label="Subtracted")
        plt.plot(x_fit, y_fit, "r", label="Gaussian fit")
        plt.hlines(best_baseline, x[mask_plot].min(), x[mask_plot].max(), linestyles="--", label="Baseline")
        plt.legend()
        plt.grid(alpha=0.4)
        plt.title(f"Peak {low}-{high}")

        # Save plot
        save_figure(
            sample_folder,
            f"{filename.replace('.csv','')}_peak_{int(mu)}_{low}_{high}.png"
        )

        # Store peak results
        all_results.append({"ROI":best_window,"mu":mu,"sigma":sigma,"R2":best_r2})

    return all_results

# ==========================================================
# Run analysis
# ==========================================================
graph_bg_1 = makegraph("PHA PreAmp Background 4838 sec; 7-1-2026.csv")
graph_bg_2 = makegraph("PreAmp Background measurement 3671sec;14-1-2026.csv")

peak_analyser(3638, "PHA PreAmp Vase Uranium 3638sec; 14-1-2026.csv", graph_bg_2)
peak_analyser(7200, "PHA PreAmp Rocks 7200sec; 15-1-2026.csv", graph_bg_2)
peak_analyser(7422, "PHA PreAmp Banana 7422 sec; 8-1-2026.csv", graph_bg_1)
peak_analyser(4835, "PHA PreAmp Potassium 4835 sec; 8-1-2026.csv", graph_bg_1)
