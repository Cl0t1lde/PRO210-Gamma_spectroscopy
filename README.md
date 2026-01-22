Gamma Spectrum Analysis-README

This project analyses gamma-ray spectra using a NaI(Tl) scintillator connected to PASCO UCS-30 universal computer spectrometer. The analysis includes background subtraction, peak identification, Gaussian fitting, and basic statistical tests to characterise gamma-ray peaks.


Two scripts are provided:

Final_plot_saver.py runs the analysis and saves the plots in the folder ''plots_2'' without orinting the results.

gamma_curve_fit.py runs the analysis and prints numerical results while also showing plots interactively.




Input data (CSV files)



CSV files follow the pattern:
```
PHA PreAmp <Sample> <AcquisitionTime>sec; <Date>.csv
```
Example:

```
PHA PreAmp Rocks 7200sec; 15-1-2026.csv
```


Each CSV file contains:

(detector settings, acquisition time, voltage, gains, etc.)

A section starting with:

```
Channel Data:
Channel,Energy,Counts
```

Only Channel and Counts are used in the analysis.

The Energy column is ignored and also not present


What the analysis does:

For each spectrum, the scripts perform:

Spectrum loading
    Reads channel–count data from the CSV file.

Background scaling and subtraction
   Background spectra are scaled by measuring time and subtracted from sample spectra.

Spectrum visualisation

   Linear-scale spectrum (raw, background, and background-subtracted)
   Log-scale spectrum (channels 0–600)

Peak analysis

   The best regions of interest (ROIs) are selected for each sample and its peaks using the R^2 value
   A Gaussian function is fitted to each candidate ROI

Statistical tests (gamma_curve_fit.py only)
   For each accepted peak, the following are calculated and printed:

   Net counts and background counts
   Signal-to-background ratio (S/B)
   Peak significance (in σ)
   Gaussian centroid (μ), width (σ), and FWHM
   Fit quality (R²)
