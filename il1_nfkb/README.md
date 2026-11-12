# IL-1–NF-κB signaling module — Python reconstruction and Streamlit app

**Author** — Nataša A. Kablar  
**Version** — 3.0  
**Release DOI** — [10.5281/zenodo.24148450](https://doi.org/10.5281/zenodo.24148450)

This package reconstructs the model reported by Nataša A. Kablar and Vladimir Kvrgić in [“Mathematical Model of IL -1- NFkB Biological Module”](https://www.ifnaworld.org/ojs/index.php/GJMS/article/download/16/7), Global Journal of Mathematical Sciences 1(1), 31–37 (2012), using the available model tables, reaction scheme and author-reviewed derivations. The app includes a 22-state ODE model, a 16-state reduced ODE model and a 16-differential + 6-algebraic DAE formulation.

## Run the app

From the unpacked folder:

```bash
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

To enable the IDA solver and comparison of all three formulations:

```bash
python -m pip install -r requirements-dae.txt
python -m streamlit run app.py
```

To regenerate the equilibrium results and run the numerical verification:

```bash
python run_equilibrium_analysis.py
python verify_equilibria.py
```

The release includes the 22-state transient, equilibrium concentrations, Jacobian spectra, and parameter-scan tables and plots in `results_equilibria/`. See `README_SR.md`, `FORMULACIJE.md` and `NALAZI_RAVNOTEZE.md` for the formulations, analysis and parameter ranges.

## Reconstruction decisions and scope

The reconstruction uses x1 for active IKK and x21 for inactive IKK. The k6 term is a basal source of IκB transcript (μM/s). The k19 term describes a first-order reaction of x22 and uses 1.2 × 10⁻⁴ s⁻¹. The reference equilibrium is locally asymptotically stable on the corresponding conservation manifold. Nine one-parameter scans identified no Hopf candidate on the tracked branch over the documented ranges. Kinetic parameter provenance and experimental calibration remain open. The models and outputs describe signaling dynamics under the specified assumptions.

## Cite this release

Kablar, N. A. (2026). *IL-1–NF-κB Signaling Module — Python Reconstruction, ODE and DAE Formulations, and Streamlit App* (Version 3.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.24148450

Please also cite the original article when referring to the model. See `CITATION.cff` for machine-readable metadata.

## Licenses

The code and Streamlit application are available under MIT. The original documentation and generated simulation results are available under CC BY 4.0. See [LICENSES.md](LICENSES.md).

## Acknowledgments

ChatGPT (OpenAI) reconstructed the model and its numerical responses in Python, identified omissions and unit inconsistencies in the available material, and developed the Streamlit interface in collaboration with Nataša A. Kablar. Nataša A. Kablar reviewed the derived equations and confirmed the corrections used in this reconstruction.
