## Overview

This framework simulates photon transport within a three-layer tissue model to determine the optimal wavelength and source-detector spacing for glucose detection. Each wavelength and source-detector spacing is given a wavelength performance index (WPI) based on glucose-sensitivity, signal-to-noise ratio, water interference, and dermal absorption sensitivity.

## Software Layout
<img width="1061" height="302" alt="image" src="https://github.com/user-attachments/assets/54c81aad-afe4-4c1a-9e6f-35fbf6c3fb4e" />

## Installing the Framework
1. Create a virtual environment
```powershell
py -3.12 -m venv .MCVenv
```
2. Activate the environment
```powershell
.\.MCVenv\Scripts\Activate.ps1
```
3. Read the .toml file. This will install all the necessary libraries
```powershell
python -m pip install -e .
```
