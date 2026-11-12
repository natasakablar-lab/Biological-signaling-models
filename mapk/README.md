# MAPK signaling module — Python reconstruction and Streamlit app

**Author** — Nataša A. Kablar  
**Version** — 5.0  
**Release DOI** — [10.5281/zenodo.24148152](https://doi.org/10.5281/zenodo.24148152)

This package reconstructs the eight-state MAPK model in Nataša A. Kablar, “MAPK Module: Biological Basis, Structure, Mathematical Model and Dynamical Analyse,” Proceedings of the 19th International Symposium on Mathematical Theory of Networks and Systems (MTNS 2010), pp. 2363–2368. The reconstruction uses the published equations and figures as its reference. [Original conference paper](https://www.conferences.hu/mtns2010/proceedings/Papers/422_259.pdf).

## Run the app

From the unpacked folder, with Python 3.10 or later:

```bash
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

To regenerate the numerical results, run these commands in order. They overwrite matching files in the results folders.

```bash
python mapk_model.py
python hopf_analysis.py
python lyapunov_analysis.py
python sensitivity_analysis.py
```

## Results and scope

The baseline equilibrium at V1 = 0.25 nM/s agrees with the concentrations printed in the paper to their reported precision. The simulations reproduce the damped baseline response and sustained oscillations at V1 = 5 nM/s over the simulated interval. The package includes eight- and five-state formulations, equilibrium continuation in V1, Hopf analysis, first Lyapunov coefficients and sensitivity scenarios. The reference scenario sets k11 and k12 to zero. The nonzero k11 and k12 and the phosphatase totals used in the additional scenarios are documented modeling assumptions. See `README_SR.md`, `HOPF_SR.md`, `LYAPUNOV_SR.md`, `SENSITIVITY_SR.md` and `STREAMLIT_SR.md` for details.

## Cite this release

Kablar, N. A. (2026). *MAPK Signaling Module — Python Reconstruction, Dynamical Analysis and Streamlit App* (Version 5.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.24148152

Please also cite the original model publication when referring to the model. See `CITATION.cff` for machine-readable metadata.

## Licenses

The code and Streamlit application are available under MIT. The original documentation and generated simulation results are available under CC BY 4.0. See [LICENSES.md](LICENSES.md).

## Acknowledgments

ChatGPT (OpenAI) reconstructed the published model and its numerical responses in Python and developed the Streamlit interface with guidance from Nataša A. Kablar. Nataša A. Kablar reviewed the reconstructed responses against the original study.
