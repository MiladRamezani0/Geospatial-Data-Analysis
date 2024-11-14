# Geospatial Data Analysis


This project is part of the "Geospatial Data Analysis" course at Politecnico di Milano, utilizing real-world geospatial data provided by ARPA Veneto. The analysis focuses on piezometer measurements from various boreholes within the Padua district, aiming to model and interpret groundwater levels using polynomial and spectral analysis techniques.

Usage
Execute the Main Code: Run the script LAB_Assignment.m to initiate the analysis.
Supporting Files: Additional .m files contain required functions.
Dataset: PZ.mat includes the time series data used in this analysis.
Comparison: The analysis is extended to additional boreholes LAB_515 and LAB_518.
Report: For a full description of the methodology and output, refer to CARLASSARA_LAB_Assignment_Report.pdf.
Introduction
The provided dataset includes borehole measurements spanning from 1999 to 2021, covering piezometer data, borehole coordinates, and elevations. Key features of the dataset include:

Timestamps converted to integer epochs.
Static water levels measured at approximately quarterly intervals.
Borehole locations and elevations within Padua's district.
Analysis Overview
Polynomial Modeling: Fit a low-degree polynomial using Smoothing Least Squares to model the deterministic component of the time series.
Spline Interpolation: Perform Exact Cubic Spline interpolation to handle missing values and create equidistant data points.
Spectral Analysis: Decompose the signal into sin and cosine components to examine its power spectrum.
Further Considerations
Due to variations in borehole data, results may vary across samples. Future improvements could include:

Clustering analysis for better sample categorization.
Enhanced data layers for pattern recognition.
Additional methods such as Linearized Least Squares and Kriging, with model evaluation using R-squared, AIC, and BIC metrics.
