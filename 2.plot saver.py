import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Create folder to save plots
PLOT_DIR = "plots"
os.makedirs(PLOT_DIR, exist_ok=True)

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
    
    df = pd.read_csv(filename, skiprows=start_row, usecols=[0,2], names=["channel","counts"])
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

# === Plot full spectrum ===
def makegraph(filename):
    df = load_spectrum(filename)
    plt.figure(figsize=(10,5))
    plt.plot(df["channel"], df["counts"], drawstyle="steps-mid")
    plt.xlabel("Channel")
    plt.ylabel("Counts")
    plt.title(f"Full Gamma Spectrum: {filename}")
    plt.grid(alpha=0.4)
    plt.xticks(np.arange(df["channel"].min(), df["channel"].max()+1, 50))
    outname = filename.replace(".csv","_full_spectrum.png")
    plt.savefig(os.path.join(PLOT_DIR, outname), dpi=300, bbox_inches="tight")
    plt.close()
    return df

# === Background subtraction + full spectrum plot ===
def plot_background_removed(t_k, filename, df_bg):
    # Determine sample label
    sample_label = "Sample"
    for key, label in SAMPLE_LABELS.items():
        if key.lower() in filename.lower():
            sample_label = label
            break

    df_p = load_spectrum(filename)
    
    # Scale background
    t_bg = 3671
    scale = t_k / t_bg
    df_bg_scaled = df_bg.copy()
    df_bg_scaled["counts_scaled"] = df_bg_scaled["counts"] * scale
    
    # Merge and subtract
    df = df_p.merge(df_bg_scaled[["channel","counts_scaled"]], on="channel", how="inner")
    df["counts_sub"] = df["counts"] - df["counts_scaled"]

    # Plot full spectrum + background removed
    plt.figure(figsize=(12,6))
    plt.plot(df["channel"], df["counts"], drawstyle="steps-mid", label=sample_label)
    plt.plot(df["channel"], df["counts_sub"], drawstyle="steps-mid", label="Background removed")
    plt.xlabel("Channel")
    plt.ylabel("Counts")
    plt.title(f"{sample_label} - Full Spectrum with Background Removed")
    plt.legend()
    plt.xticks(np.arange(df["channel"].min(), df["channel"].max()+1, 50))
    plt.grid(alpha=0.4)
    
    # Save figure
    outname = f"{filename.replace('.csv','')}_full_background_removed.png"
    plt.savefig(os.path.join(PLOT_DIR, outname), dpi=300, bbox_inches="tight")
    plt.close()


# ==========================================================
# Run analysis
# ==========================================================
graph_bg_1 = makegraph("PHA PreAmp Background 4838 sec; 7-1-2026.csv")
graph_bg_2 = makegraph("PreAmp Background measurement 3671sec;14-1-2026.csv")
#graph_b  = makegraph("PHA PreAmp Banana 7422 sec; 8-1-2026.csv")
#graph_p  = makegraph("PHA PreAmp Potassium 4835 sec; 8-1-2026.csv")
#graph_t  = makegraph("PHA PreAmp Tungsten 4196sec; 14-1-2026.csv")
#graph_u  = makegraph("PHA PreAmp Vase Uranium 3638sec; 14-1-2026.csv")
#graph_r  = makegraph("PHA PreAmp Rocks 7200sec; 15-1-2026.csv")

print("\nResults for Uranium:")
#results_t = plot_background_removed(3638, "PHA PreAmp Vase Uranium 3638sec; 14-1-2026.csv", graph_bg_1)

print("\nResults for Rocks:")
#results_t = plot_background_removed(7200, "PHA PreAmp Rocks 7200sec; 15-1-2026.csv", graph_bg_1)

print("\nResults for Tungsten:")
#results_t = plot_background_removed(4196, "PHA PreAmp Tungsten 4196sec; 14-1-2026.csv", graph_bg_1)

print("\nResults for Background 2:")
results_t = plot_background_removed(4838, "PHA PreAmp Background 4838 sec; 7-1-2026.csv", graph_bg_2)

print("\nResults for Banana:")
#results_b = plot_background_removed(7422, "PHA PreAmp Banana 7422 sec; 8-1-2026.csv", graph_bg_1)

print("\nResults for Potassium:")
#results_k = plot_background_removed(4835, "PHA PreAmp Potassium 4835 sec; 8-1-2026.csv", graph_bg_1)