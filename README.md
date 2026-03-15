# Adaptive Respiratory Motion Correction Method for Hepatic 4D Flow MRI


This repository provides MATLAB implementations of Adaptive Respiratory Motion Correction(AMOCO) reconstruction algorithms presented in:

> **Ruiyu Cao, Shaohang Li, Aiqi Sun, et al.**  
> *Adaptive Respiratory Motion Correction Framework for Ungated Hepatic 4D Flow MRI.* 
> DOI: **XXX**

The primary contribution is **AMOCO**—a image-based motion-corrected reconstruction method designed to correct motion artifacts abdominal **4D flow** MRI .

---

## Overview

AMOCO is a framwork of optimized acquisition and motion corrected reconstruction to mitigate motion artifacts and improve quantification accuracy of abdominal 4D flow.

This repository includes source code of VPD pattern design and AMOCO reconstrction with intra-bin and inter-bin motino correction.

---

## Key Features

- **VPD pattern design**: 
  - implementations of variable-probability Cartesian sampling pattern with variable Poisson disk constraints

- **amoco_code**: 
  - **ADMM_LLR_TV_uncorrected**: conventional CS reconstruction algorithm with locally low rank and total variation regularization for uncorrected hepatic 4D flow imaging data
  - **ADMM_LLR_TV_amoco**: AMOCO reconstruction algorithm with intra-bin soft-gating weights and inter-bin motion fields for motion resolved hepatic 4D flow imaging data
  
---

## Getting Started

### Requirements
- **MATLAB R2019a or later**

### Installation
1. Clone this repository:
   - `git clone <repo_url>`
2. Open MATLAB and add the repository to your path (or run scripts from within the repository folder).

---

## Usage

### AMOCO
AMOCO reconstruction: run main_recon.main amoco_code directory

### VPD pattern design
在VPD_pattern_design中的generate_3D_mask.m中设置相关参数后运行

---

## Datasets

The example dataset for AMOCO reconstruction consist of:
- **Ungated Cartesian undersampled** acquisitions using VPD sampling pattern with respiratory bellows recording
- `.mat` files containing
   - **k-space data**
   - **sampling mask**
   - **soft-gating weight**
   - **coil sensitive map**
   - **motion fields**

---

## Citation

If you use this code, please cite:

> Ruiyu Cao, Shaohang Li, Aiqi Sun, et al. Adaptive Respiratory Motion Correction Framework for Ungated Hepatic 4D Flow MRI.

---

## Authors

- **Ruiyu Cao** (Primary Developer) — rycao22@m.fudan.edu.cn
- **Shaohang Li** — shli23@m.fudan.edu.cn
- **Aiqi Sun** —
- **Hao Li** -

---

## Notes

- Please ensure the dataset files are placed in the expected directories referenced by the study scripts.
- If you encounter path or dependency issues, confirm MATLAB has the repository folders added via `addpath(genpath(pwd))`.
