# Testing fkPT against MG-Copter

Analysis code for my MSc dissertation, *Testing Galaxy Clustering Predictions in Modified Gravity Theories* (MSc Theoretical Physics, University of Edinburgh, 2025). Supervisors: Dr Benjamin Bose and Dr Alkistis Pourtsidou.

## What this project tests
Upcoming galaxy surveys such as Euclid need fast, accurate models of how galaxies cluster. Modified-gravity theories change this clustering, but numerical predictions are expensive, so faster perturbative approximations are used instead. This project asks: **how well do two perturbative approximations reproduce the numerical MG-Copter predictions for the redshift-space galaxy power spectrum?**

## Approach
- **Approximations compared:** Einstein-de Sitter (EdS) and fk-perturbation theory (fkPT)
- **Benchmark:** numerical predictions from MG-Copter
- **Gravity models:** General Relativity / LCDM and Hu-Sawicki f(R) gravity (F4 and F5 configurations)
- **Setup:** CDM and halo bias prescriptions at redshift z = 1
- **Observables:** redshift-space power-spectrum multipoles (monopole P0, quadrupole P2, hexadecapole P4), compared through ratio diagnostics
- **Statistics:** covariance-weighted reduced chi-squared as the maximum fitted wavenumber (k_max) is varied
- **Extras:** refitting of effective-field-theory counterterms, different nuisance-parameter configurations, and Euclid-like statistical uncertainties

## Main result
EdS was generally the more consistent approximation across the scales tested; fkPT gave limited improvements at selected scales. These conclusions apply only to the cosmology, bias prescriptions, parameter setup and numerical configuration used here.

## Limitations
- A single redshift (z = 1) and a limited set of f(R) configurations
- Results depend on the chosen bias model and nuisance-parameter setup
- Statistical uncertainties are Euclid-like approximations, not a full survey analysis

## Tools
Python (NumPy, Matplotlib, Jupyter), C++ analysis scripts, MG-Copter, fkPT / ACTio-ReACTio.

## Acknowledgements
Thanks to Dr Benjamin Bose and Dr Alkistis Pourtsidou for supervision, and to the authors of MG-Copter and fkPT / ACTio-ReACTio.

## Contact
Siddhi S. Gandhi - [LinkedIn](https://www.linkedin.com/in/siddhi0802/)
