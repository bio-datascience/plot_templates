# Correlation Plot Notebook

## Overview
This Jupyter notebook (`Correlation_plot.ipynb`) computes and visualizes the pairwise correlation matrix of bacterial abundance data. It loads two DataFrames, (1) abundance, and (2) taxonomy information. It calculates the correlation coefficients, and displays them in a heatmap ordered based on the taxonomy.

## Dependencies
- Python 3.7 or higher
- pandas
- seaborn
- matplotlib
- ete3
- Bio

Install the required packages via pip:
```bash
pip install pandas seaborn matplotlib ete3 Toytree