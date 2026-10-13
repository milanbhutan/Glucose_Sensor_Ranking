import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# from pysr import PySRRegressor as SR
# import physo
# import torch
# from scipy.optimize import least_squares
# from pathlib import Path

import tissue_model.layer_calc as calc
import tissue_model.tissue_params as params
import analysis.error_analysis as err
import analysis.uncertainty_analysis as uncertain
import photon_sim.monte_carlo as mc


# uncomment to run monte carlo simulation


sim_wavelength_nm = 1314

mu_s_values = np.arange(5, 14, 1)
A = []
for mu_s in mu_s_values:
    val_A = mc.run_simulation(
            [
            params.epi_mu_a_sim[sim_wavelength_nm], 
            params.epi_mu_s_sim[sim_wavelength_nm], 
            params.epi_g_sim[sim_wavelength_nm], 
            params.epi_n_sim ],

            [
            params.derm_mu_a_sim[sim_wavelength_nm], 
            mu_s, 
            params.derm_g_sim[sim_wavelength_nm], 
            params.derm_n_sim],

            [
            params.sub_mu_a_sim[sim_wavelength_nm], 
            params.sub_mu_s_sim[sim_wavelength_nm], 
            params.sub_g_sim[sim_wavelength_nm], 
            params.sub_n_sim]
    )
    A.append(val_A)




plt.figure()
plt.scatter(mu_s_values, A)
plt.xlabel("$\mu_s$ $(mm^-1)$")
plt.ylabel("$A$ $(a.u.)$")

m, b = np.polyfit(mu_s_values, A, 1)
print(m)

plt.show()