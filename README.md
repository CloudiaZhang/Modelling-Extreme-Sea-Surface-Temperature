# Sample of Work 1

This repository contains code and report materials for a project on modelling extreme sea surface temperatures (SSTs) in the Red Sea using Generalized Extreme Value (GEV) theory.

## Project overview

The project studies extreme SSTs in the Red Sea using GEV models. The analysis includes both a simulation study and an application to real SST data from 108 locations across the Red Sea.

The main goals of the project are:
- to examine the behaviour of GEV fitting under different block sizes and sample lengths
- to model extreme SSTs using geographical and environmental covariates
- to compare candidate models using AIC
- to estimate return levels and assess model fit through diagnostic plots

## Methods used

- simulation of SST data from a normal distribution
- extraction of block maxima
- GEV fitting using the `ismev` and `texmex` package
- non-stationary GEV modelling with spatial covariates
- model comparison using AIC
- diagnostic assessment using Q-Q plots, density plots, and return level plots

## Files

- `01simulation_study.R` — simulation study for block maxima, GEV fitting, and return-level comparisons under simulated SST data
- `02gev_model_fitting.R` — real-data analysis, covariate construction, candidate GEV models, AIC comparison, and final model fitting
- `03gev_diagnostics_return_levels.R` — diagnostic plots, return-level estimation, and exploratory dependence checks for selected locations
- `Allcodes.Rmd` — full R Markdown source containing the combined analysis workflow
- `Sample_of_Work_1_GEV_Modelling_of_Extreme_SST.pdf` — polished report

## Data

This project uses external data files for the empirical analysis, including:
- `Red_Sea_data.RData`
- `bathymetry.xlsx`

## Note

This repository contains supporting code for the accompanying sample of work. The code includes both the main analysis and some exploratory steps taken during model development.
