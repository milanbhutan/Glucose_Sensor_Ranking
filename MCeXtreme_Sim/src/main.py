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


df_dmu_a_dCg = (
    pd.read_csv("Glucose_Absorption_Data/ge2024_dmu_a_dCg.csv")
    .set_index("wavelength_nm")
    .sort_index()
)
dmu_a_dCg = df_dmu_a_dCg.loc[: ,"dmu_a_dCg_mm_inv_per_mmol_L"]

dmu_s_dCg = -0.0029



plt.figure()
plt.plot(df_dmu_a_dCg)
plt.title("$d\mu_a/dC_g$ from Ge 2024", fontsize = 18)
plt.xlabel("Wavelength (nm)", fontsize = 18)
plt.xticks(fontsize = 18)
plt.ylabel("$d\mu_a/dC_g$", fontsize = 18)
plt.yticks(fontsize = 18)


wavelength_range = range(1000, 1701, 10)
sds = [1.7, 2.0, 2.3, 2.6, 2.9]

dA_dCg = pd.DataFrame(
    np.nan,
    index = wavelength_range,
    columns = sds
)


for sim_wavelength_nm in wavelength_range:
#Following Ge 2024 procedure for computing dA/dCg
    #baseline
    A_array_base = mc.run_simulation(
        [
        params.glass_mu_a_sim,
        params.glass_mu_s_sim,
        params.glass_g_sim,
        params.glass_n_sim
        ],

        [
        params.epi_mu_a_sim[sim_wavelength_nm], 
        params.epi_mu_s_sim[sim_wavelength_nm], 
        params.epi_g_sim[sim_wavelength_nm], 
        params.epi_n_sim],

        [
        params.derm_mu_a_sim[sim_wavelength_nm], 
        params.derm_mu_s_sim[sim_wavelength_nm], 
        params.derm_g_sim[sim_wavelength_nm], 
        params.derm_n_sim],

        [
        params.sub_mu_a_sim[sim_wavelength_nm], 
        params.sub_mu_s_sim[sim_wavelength_nm], 
        params.sub_g_sim[sim_wavelength_nm], 
        params.sub_n_sim]
    )

    #Increment mu_a by 1
    A_array_dmu_a = mc.run_simulation(
        [
        params.glass_mu_a_sim,
        params.glass_mu_s_sim,
        params.glass_g_sim,
        params.glass_n_sim
        ],

        [
        params.epi_mu_a_sim[sim_wavelength_nm], 
        params.epi_mu_s_sim[sim_wavelength_nm], 
        params.epi_g_sim[sim_wavelength_nm], 
        params.epi_n_sim ],

        [
        params.derm_mu_a_sim[sim_wavelength_nm] + 1, 
        params.derm_mu_s_sim[sim_wavelength_nm], 
        params.derm_g_sim[sim_wavelength_nm], 
        params.derm_n_sim],

        [
        params.sub_mu_a_sim[sim_wavelength_nm], 
        params.sub_mu_s_sim[sim_wavelength_nm], 
        params.sub_g_sim[sim_wavelength_nm], 
        params.sub_n_sim]
    )

    #Increment mu_s by 1
    A_array_dmu_s = mc.run_simulation(
        [
        params.glass_mu_a_sim,
        params.glass_mu_s_sim,
        params.glass_g_sim,
        params.glass_n_sim
        ],
        
        [
        params.epi_mu_a_sim[sim_wavelength_nm], 
        params.epi_mu_s_sim[sim_wavelength_nm], 
        params.epi_g_sim[sim_wavelength_nm], 
        params.epi_n_sim 
        ],

        [
        params.derm_mu_a_sim[sim_wavelength_nm], 
        params.derm_mu_s_sim[sim_wavelength_nm] + 1, 
        params.derm_g_sim[sim_wavelength_nm], 
        params.derm_n_sim
        ],

        [
        params.sub_mu_a_sim[sim_wavelength_nm], 
        params.sub_mu_s_sim[sim_wavelength_nm], 
        params.sub_g_sim[sim_wavelength_nm], 
        params.sub_n_sim
        ]
    )
    dA_dCg.loc[sim_wavelength_nm] = (
        (A_array_dmu_s - A_array_base) * dmu_s_dCg * params.derm_mu_s_sim[sim_wavelength_nm] 
        +
        (A_array_dmu_a - A_array_base) * dmu_a_dCg.loc[sim_wavelength_nm]
        )
  


dA_dCg.plot()
plt.xlabel("Wavelength (nm)", fontsize = 18)
plt.xticks(fontsize = 18)
plt.ylabel("$dA/dC_g$", fontsize = 18)
plt.yticks(fontsize = 18)
plt.legend(["SDS_1.7mm", "SDS_2.0mm", "SDS_2.3mm", "SDS_2.6mm", "SDS_2.9mm"])


plt.show()