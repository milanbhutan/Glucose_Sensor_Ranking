## Overview

This framework simulates photon transport within a three-layer tissue model to determine the optimal wavelength and source-detector spacing for glucose detection. Each wavelength and source-detector spacing is given a wavelength performance index (WPI) based on glucose-sensitivity, signal-to-noise ratio, water interference, and dermal absorption sensitivity.

## Software Layout
<img width="1061" height="302" alt="image" src="https://github.com/user-attachments/assets/54c81aad-afe4-4c1a-9e6f-35fbf6c3fb4e" />



## GPU Requirements
Tested hardware: NVIDIA GeForce RTX 4060.

Generating Monte Carlo simulations using sim_gen.py requires:
- GPU: NVIDIA GPU with CUDA support.
- GPU driver: A compatible, up-to-date NVIDIA graphics driver.
- GPU memory: Sufficient VRAM for the selected tissue volume, time gates, and simulation settings.


## Running the Repo
The recommended development tool for running this repo is VS Code. VS Code provides an intuitive way to navigate through directories and run python modules.




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
1. Open the VSCode terminal with, Ctrl + `, and create a virtual environment
```powershell
py -3.12 -m venv .MCVenv
```
2. Activate the environment
```powershell
.\.MCVenv\Scripts\Activate.ps1
```
3. Read the .toml file. This will install all the necessary libraries. Ignore the VS Code warning
```powershell
python -m pip install -e .
```

Now you have everything you need to run the .py files. To run a module in VS code, click the play button in the top right.

### main.py
This will generate a plots of the inputs for the WPI, aswell as a plot of the WPI for a range of wavelengths and source detector spacings.

### sim_gen.py
This will run photon transport simulations and regenerate the inputs used for the WPI. These inputs depend on the selected tissue layer thicknesses, tissue optical properties, photon count, and batch size. Currently this module will ran for apporximately one hour. The run time can be reduced by decreasing photon count and batch size but this will reduce the accuracy of the simulation.

### 600_1000nm.py
This will plot the absorption and reduced scattering spectra within 600 - 1000 nm. The absorption spectra is modeled from Jacques 2013[] and Laptian 2024[]. The scattering spectra is provided by Jacques 2013[]

### 1001_1700nm.py
This will plot the absorption and reduced scattering spectra within 1001 - 1700 nm, both of which are provided by Iino 2003.

### exp_data.py
This will plot the tissue constituent absorption spectra
