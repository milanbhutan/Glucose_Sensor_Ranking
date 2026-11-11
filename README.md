## Overview

This framework simulates photon transport within a three-layer tissue model to determine the optimal wavelength and source-detector spacing for glucose detection. Each wavelength and source-detector spacing is given a wavelength performance index (WPI) based on glucose-sensitivity, dermal absorption sensitivity, detected optical power, water interference.

## Software Layout
<img width="1066" height="311" alt="image" src="https://github.com/user-attachments/assets/353b5136-5a76-4b43-8a77-f6be7d7efb95" />



## GPU Requirements
Tested hardware: NVIDIA GeForce RTX 4060.

Generating Monte Carlo simulations using sim_gen.py requires:
- GPU: NVIDIA GPU with CUDA support.
- GPU driver: A compatible, up-to-date NVIDIA graphics driver.
- GPU memory: Sufficient VRAM for the selected tissue volume, time gates, and simulation settings.


## Setup
The recommended development tool for running this repo is VS Code. VS Code provides an intuitive way to navigate through directories and run Python modules.

Download VS Code Here:
https://code.visualstudio.com/download

1. In terminal, clone the repository
```powershell
git clone "https://github.com/milanbhutan/Glucose_Sensor_Ranking.git"
```
2. Open the repository in VS Code
```powershell
code Glucose_Sensor_Ranking
```
3. Open the VS Code terminal with, Ctrl + `, and create a virtual environment
```powershell
py -3.12 -m venv .MCVenv
```
4. Activate the environment
```powershell
.\.MCVenv\Scripts\Activate.ps1
```
5. Read the .toml file. This will install all the necessary libraries. Ignore the VS Code warning
```powershell
python -m pip install -e .
```

## Running .py files
To run a .py file in VS Code, click the play button in the top right.

### main.py
This will generate plots of the inputs for the WPI, and plot the WPI for a range of wavelengths and source detector spacings.

### sim_gen.py
This will run photon transport simulations and regenerate the inputs used for the WPI. These inputs depend on the selected tissue layer thicknesses, tissue optical properties, photon count, and batch size. Currently this module will take about 1.5 hours to run. The runtime can be reduced by decreasing photon count and batch size, but this will reduce the accuracy of the simulation.

### 600_1000nm.py
This will plot the constituent-based absorption and the reduced scattering spectra within 600 - 1000 nm. The absorption spectra is modeled from Jacques 2013 [1] and Lapitan 2024 [2]. The scattering spectra is provided by Jacques 2013 [1]

### 1000_1700nm.py
This will plot the absorption and reduced scattering spectra within 1001 - 1700 nm, both of which are provided by Iino 2003[3].

### 600_1700nm.py
This will plot the entire constituent-based absorption spectra and the Salomatina-based[4] reduced scattering spectra from 600 - 1700 nm. 

### exp_data.py
This will plot the tissue constituent absorption spectra.

## References
[1] S. L. Jacques, "Optical properties of biological tissues: A review," Physics in Medicine & Biology, vol. 58, no. 11, pp. R37–R61, 2013, doi: 10.1088/0031-9155/58/11/R37.

[2] D. G. Lapitan, A. P. Tarasov, M. E. Shtyflyuk, and D. A. Rogatkin, "On the origin of the photoplethysmography signal: Modeling of volumetric and aggregation effects," Photonics, vol. 11, no. 7, Art. no. 637, 2024, doi: 10.3390/photonics11070637.

[3] K. Iino, K. Maruo, H. Arimoto, K. Hyodo, T. Nakatani, and Y. Yamada, "Monte Carlo simulation of near infrared reflectance spectroscopy in the wavelength range from 1000 nm to 1900 nm," Optical Review, vol. 10, no. 6, pp. 600–606, 2003, doi: 10.1007/s10043-003-0600-z.

[4] E. Salomatina, B. Jiang, J. Novak, and A. N. Yaroslavsky, "Optical properties of normal and cancerous human skin in the visible and near-infrared spectral range," Journal of Biomedical Optics, vol. 11, no. 6, 064026, 2006, doi: 10.1117/1.2398928.

