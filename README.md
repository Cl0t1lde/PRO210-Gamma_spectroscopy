# Gamma Spectroscopy Analysis

## Overview
This project analyzes gamma spectroscopy data. It can:  
- Plot full gamma spectra  
- Subtract background spectra  
- Identify and quantify peaks using Gaussian fits  

All generated plots are saved in the `plots_2` folder for easy review.

---

## Folder Structure

        /<folder>
        │
        ├─ plots_2/ # Saved plots organized by sample
        │ ├─ Overview/ # Full spectrum overview plots
        │ └─ <SampleName>/ # Linear, log, and peak plots for each sample
        │
        ├─ gamma_data.csv # Example gamma spectroscopy data
        ├─ script_peak_analysis_1.py # Original analysis script (interactive plots)
        ├─ script_peak_analysis_2.py # Improved script with saved plots
        ├─ .DS_Store # macOS system file (can be ignored)

---

## Files

- **`script_peak_analysis_1.py`** – Loads spectra, plots them, subtracts background, and fits Gaussian peaks. Plots are shown interactively.  
- **`script_peak_analysis_2.py`** – Same as above, but automatically saves all plots in `plots_2`.  
- **`gamma_data.csv`** – Gamma spectroscopy data (channels and counts).  

  ### File Naming Convention

  All data files follow this pattern:  

  **`PHA PreAmp <SampleName> <AcquisitionTime>sec; <Date>.csv`**  

  Where:
    - `<SampleName>` – Name of the sample (e.g., Uranium, Banana, Potassium)  
    - `<AcquisitionTime>` – Duration of the spectroscopy measurement in seconds
    - `<Date>` – Date of the data acquisition (format: day-month-year)

  **Example:**  
        PHA PreAmp Vase Uranium 3638sec; 14-1-2026.csv

- **`.DS_Store`** – macOS system file; safe to ignore.  
- **`plots_2/`** – Folder containing all saved plots from the analysis.

---

## Usage

1. Ensure the data file(s) are in the same folder as the scripts.  
2. Run the script with Python 3.x:  
   ```bash
   python script_peak_analysis_2.py
3. Generated plots will be saved in plots_2/ in subfolders:

## Notes

Background spectra are automatically scaled by acquisition time before subtraction.
Peak regions are predefined; Gaussian fits extract peak position, width, and significance.
.DS_Store files appear on macOS but do not affect scripts or data.
