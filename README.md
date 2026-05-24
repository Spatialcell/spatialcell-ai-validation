SpatialCell AI Validation Framework
This repository contains the validation code and analysis for the manuscript:

SpatialCell AI achieves reference-free single-cell resolution from spot-based spatial transcriptomics through morphology-guided enhancement
Abdalla Elbialy, Scientific Reports (2026).
DOI: 10.1038/s41598-026-55246-w

Overview
Distribution-based validation framework comparing SpatialCell AI outputs against Xenium single-cell ground truth data. The framework evaluates per-cell expression reconstruction across multiple spatial transcriptomics input platforms (Visium 55 µm, Visium HD 16 µm, and Visium HD 8 µm) using a panel of distribution-based metrics.
Requirements

Python 3.8+
See requirements.txt for dependencies

Install dependencies with:
pip install -r requirements.txt
Usage
validation/Validation_script_v2.ipynb
Data Availability

Validation framework code: https://github.com/Spatialcell/spatialcell-ai-validation
Transformed single-cell resolution datasets: https://zenodo.org/records/15794331
Raw spatial transcriptomics data (Oliveira et al.): https://www.10xgenomics.com/platforms/visium/product-family/dataset-human-crc
SpatialCell AI web platform: https://spatialcell.tech

Citation
If you use this validation framework, please cite:

Elbialy, A. SpatialCell AI achieves reference-free single-cell resolution from spot-based spatial transcriptomics through morphology-guided enhancement. Scientific Reports (2026). https://doi.org/10.1038/s41598-026-55246-w

Scope of this repository
This repository contains only the validation and analysis code used to evaluate SpatialCell AI against Xenium ground-truth data for the manuscript above. It does not contain the SpatialCell AI software platform itself. The SpatialCell AI product is a separate, proprietary, commercially licensed service available at https://spatialcell.tech and is not included in or covered by this repository or its license.
License
The validation code in this repository is licensed under the MIT License — see the LICENSE file for details. This license applies to the validation code only, not to the SpatialCell AI platform/product.
