# Biological signaling models

Streamlit applications for the MAPK and IL-1–NF-κB models by Nataša A. Kablar.

| App | Streamlit entrypoint | Archived release |
| --- | --- | --- |
| MAPK | `mapk/app.py` | [Zenodo 5.0](https://doi.org/10.5281/zenodo.24148152) |
| IL-1–NF-κB | `il1_nfkb/app.py` | [Zenodo 3.0](https://doi.org/10.5281/zenodo.24148450) |

## Deploy on Streamlit Community Cloud

Push the contents of this folder to a public GitHub repository. In Streamlit Community Cloud, create one app per entrypoint using the same repository and branch. Select Python 3.12 under Advanced settings. Each entrypoint has its own `requirements.txt`. The NF-κB app installs scikit-sundae for the SUNDIALS IDA solver.

Each subdirectory contains the release citation, numerical results, and license notices. The Python code and applications are MIT-licensed. Original documentation and generated simulation results are CC BY 4.0.
