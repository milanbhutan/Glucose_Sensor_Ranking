import pmcx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import layer
import error_analysis as err
from sklearn.metrics import mean_squared_error
from scipy.optimize import least_squares


lower_wavelength = 600 # original: 600
upper_wavelength = 1060

# Generate random tissue compositions here

rng = np.random.default_rng(45)

epi_mu_a_sim_set = []
derm_mu_a_sim_set = []
sub_mu_a_sim_set = []

derm_params = [
    "blood",
    "oxygen_saturation",
    "water",
    "gamma"
]

sub_params = [
    "blood",
    "oxygen_saturation",
    "water",
    "fat",
    "gamma",
]

derm_comp_upper = {
    "blood" : 0.015,
    "oxygen_saturation" : 0.65,    
    "water" : 0.75,
    "gamma" : 0.8

}

derm_comp_lower = {
    "blood" : 0.002,
    "oxygen_saturation" : 0.39,    
    "water" : 0.65,
    "gamma" : 0.2

}

sub_comp_upper = {
    "blood" : 0.010,
    "oxygen_saturation" : 0.80,    
    "water" : 0.30,
    "fat" : 0.80,
    "gamma" : 0.3,

}

sub_comp_lower = {
    "blood" : 0.002,
    "oxygen_saturation" : 0.50,    
    "water" : 0.10,
    "fat" : 0.55,
    "gamma" : 0.2,

}

for i in range(200):
    derm_blood = rng.uniform(derm_comp_lower["blood"], derm_comp_upper["blood"])
    derm_oxygen_saturation = rng.uniform(derm_comp_lower["oxygen_saturation"], derm_comp_upper["oxygen_saturation"])
    derm_water = rng.uniform(derm_comp_lower["water"], derm_comp_upper["water"])
    derm_gamma = rng.uniform(derm_comp_lower["gamma"], derm_comp_upper["gamma"])
    derm_mu_a_sim_set.append(
        layer.mu_a_derm(
            derm_blood, 
            derm_oxygen_saturation, 
            derm_water, 
            layer.derm_Jacques["fat"], 
            layer.derm_Jacques["melanosome"],
            derm_gamma,
    )
    )

for i in range(200):
    sub_blood = rng.uniform(sub_comp_lower["blood"], sub_comp_upper["blood"])
    sub_oxygen_saturation = rng.uniform(sub_comp_lower["oxygen_saturation"], sub_comp_upper["oxygen_saturation"])
    sub_water = rng.uniform(sub_comp_lower["water"], sub_comp_upper["water"])
    sub_fat = rng.uniform(sub_comp_lower["fat"], sub_comp_upper["fat"])
    sub_gamma = rng.uniform(sub_comp_lower["gamma"], sub_comp_upper["gamma"])
    sub_mu_a_sim_set.append(
        layer.mu_a_sub(
            sub_blood, 
            sub_oxygen_saturation, 
            sub_water, 
            sub_fat, 
            layer.sub_Jacques["melanosome"],
            sub_gamma
    )
    )

# End of generating random tissue compositions


norm_wavelength_nm = np.arange(0, 2001, 1)


epi_mu_a_sim = layer.mu_a_epi(layer.epi_Jacques["blood"], 
                                    layer.epi_Jacques["oxygen_saturation"], 
                                    layer.epi_Jacques["water"], 
                                    layer.epi_Jacques["fat"], 
                                    layer.epi_Jacques["melanosome"]
                                    )

epi_mu_a_sim_Jolivot_SPT_II = layer.mu_a_epi(layer.epi_Jolivot_SPT_II["blood"], 
                                    layer.epi_Jolivot_SPT_II["oxygen_saturation"], 
                                    layer.epi_Jolivot_SPT_II["water"], 
                                    layer.epi_Jolivot_SPT_II["fat"], 
                                    layer.epi_Jolivot_SPT_II["melanosome"]
                                    )

epi_mu_a_sim_Jolivot_SPT_III = layer.mu_a_epi(layer.epi_Jolivot_SPT_III["blood"], 
                                    layer.epi_Jolivot_SPT_III["oxygen_saturation"], 
                                    layer.epi_Jolivot_SPT_III["water"], 
                                    layer.epi_Jolivot_SPT_III["fat"], 
                                    layer.epi_Jolivot_SPT_III["melanosome"]
                                    )

epi_mu_a_sim_Jolivot_SPT_IV = layer.mu_a_epi(layer.epi_Jolivot_SPT_IV["blood"], 
                                    layer.epi_Jolivot_SPT_IV["oxygen_saturation"], 
                                    layer.epi_Jolivot_SPT_IV["water"], 
                                    layer.epi_Jolivot_SPT_IV["fat"], 
                                    layer.epi_Jolivot_SPT_IV["melanosome"]
                                    )

epi_mu_a_sim_Jolivot_SPT_V = layer.mu_a_epi(layer.epi_Jolivot_SPT_V["blood"], 
                                    layer.epi_Jolivot_SPT_V["oxygen_saturation"], 
                                    layer.epi_Jolivot_SPT_V["water"], 
                                    layer.epi_Jolivot_SPT_V["fat"], 
                                    layer.epi_Jolivot_SPT_V["melanosome"]
                                    )

epi_mu_a_sim_Jolivot_SPT_VI = layer.mu_a_epi(layer.epi_Jolivot_SPT_VI["blood"], 
                                    layer.epi_Jolivot_SPT_VI["oxygen_saturation"], 
                                    layer.epi_Jolivot_SPT_VI["water"], 
                                    layer.epi_Jolivot_SPT_VI["fat"], 
                                    layer.epi_Jolivot_SPT_VI["melanosome"]
                                    )

derm_mu_a_sim = layer.mu_a_derm(layer.derm_Jacques["blood"], 
                                  layer.derm_Jacques["oxygen_saturation"], 
                                  layer.derm_Jacques["water"], 
                                  layer.derm_Jacques["fat"], 
                                  layer.derm_Jacques["melanosome"]
                                  )

sub_mu_a_sim_Jacques = layer.mu_a_sub(layer.sub_Jacques["blood"], 
                                       layer.sub_Jacques["oxygen_saturation"], 
                                       layer.sub_Jacques["water"], 
                                       layer.sub_Jacques["fat"], 
                                       layer.sub_Jacques["melanosome"]
                                       )

sub_mu_a_sim = layer.mu_a_sub(layer.sub_Ganesan["blood"], 
                                               layer.sub_Ganesan["oxygen_saturation"],
                                               layer.sub_Ganesan["water"],
                                               layer.sub_Ganesan["fat"],
                                               layer.sub_Ganesan["melanosome"]
                                               )



# For scattering coefficent
epi_g = 0.9
derm_g = 0.9
sub_g = 0.9

epi_mu_s_sim = layer.mu_s_layer("epidermis", epi_g)
derm_mu_s_sim = layer.mu_s_layer("dermis", derm_g)
sub_mu_s_sim = layer.mu_s_layer("subcutaneous", sub_g)

# For reduced scattering coefficent
epi_g_red = 0
derm_g_red = 0
sub_g_red = 0

epi_mu_s_red_sim = layer.mu_s_layer("epidermis", epi_g_red)
derm_mu_s_red_sim = layer.mu_s_layer("dermis", derm_g_red)
sub_mu_s_red_sim = layer.mu_s_layer("subcutaneous", sub_g_red)


epi_n = 1.4
derm_n = 1.4
sub_n = 1.4


# Obtain data for experimental absorption and scattering spectras for three-layer tissue model from Shimojo 2020

# Absorption Data
df_Wan_epi_mu_a = (
    pd.read_csv("Exp_Absorption/wan_epidermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Shimojo_epi_mu_a = (
    pd.read_csv("Exp_Absorption/shimojo_epidermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)


df_Simpson_derm_mu_a = (
    pd.read_csv("Exp_Absorption/simpson_dermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
mu_a_derm_Simpson = df_Simpson_derm_mu_a.loc[lower_wavelength:upper_wavelength, "caucasian_dermis_mu_a_mm^-1"]

df_Shimojo_derm_mu_a = (
    pd.read_csv("Exp_Absorption/shimojo_dermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Ganesan_sub_mu_a = (
    pd.read_csv("Exp_Absorption/ganesan_subcutaneous_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
mu_a_sub_Ganesan = df_Ganesan_sub_mu_a.loc[lower_wavelength:upper_wavelength, "mu_a_mm^-1"]

df_Simpson_sub_mu_a = (
    pd.read_csv("Exp_Absorption/simpson_subcutaneous_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
mu_a_sub_Simpson = df_Simpson_sub_mu_a.loc[lower_wavelength:upper_wavelength, "caucasian_subdermis_mu_a_mm^-1"]

# Scattering Data
df_Shimojo_epi_red_mu_s = (
    pd.read_csv("Exp_Scattering/shimojo_epidermis_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Shimojo_derm_red_mu_s = (
    pd.read_csv("Exp_Scattering/shimojo_dermis_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Simpson_derm_red_mu_s = (
    pd.read_csv("Exp_Scattering/simpson_dermis_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Shimojo_sub_red_mu_s = (
    pd.read_csv("Exp_Scattering/shimojo_subcutaneous_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Ganesan_sub_red_mu_s = (
    pd.read_csv("Exp_Scattering/ganesan_subcutaneous_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Simpson_sub_red_mu_s = (
    pd.read_csv("Exp_Scattering/simpson_subcutaneous_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)






# ===============================================
# Tissue spectra from 600 - 1100nm
# ===============================================

# Note: This subcutaneous spectra is only valid within 429-1100 nm where reliable tissue constituent abosroption data is availible
fig_mu_a, axes_mu_a = plt.subplots(1, 3, layout = "constrained")
fig_mu_a.get_layout_engine().set(h_pad = 0.1, w_pad = 0.1)
fig_mu_a.suptitle("Absorption Spectras for Layers in Three-Layer Tissue Model", fontsize = 18)

for ax in axes_mu_a:
    ax.set_xlabel("Wavelength (nm)", fontsize = 18)
    ax.set_ylabel("Absorption (1/mm)", fontsize = 18)
    ax.tick_params(axis = "both", labelsize = 18)


# Simulated absorption spectra plots
axes_mu_a[0].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  epi_mu_a_sim[lower_wavelength:upper_wavelength], 
                  label = "Analytical Models", 
                  color = "red")
axes_mu_a[0].text(lower_wavelength, 
                  epi_mu_a_sim[lower_wavelength], 
                  "M=2.5%", 
                  color = "red")

axes_mu_a[0].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  epi_mu_a_sim_Jolivot_SPT_II[lower_wavelength:upper_wavelength], 
                  color = "red")


axes_mu_a[0].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  epi_mu_a_sim_Jolivot_SPT_III[lower_wavelength:upper_wavelength], 
                  color = "red")

axes_mu_a[0].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  epi_mu_a_sim_Jolivot_SPT_IV[lower_wavelength:upper_wavelength], 
                  color = "red")

axes_mu_a[0].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  epi_mu_a_sim_Jolivot_SPT_V[lower_wavelength:upper_wavelength], 
                  color = "red")

axes_mu_a[0].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  epi_mu_a_sim_Jolivot_SPT_VI[lower_wavelength:upper_wavelength], 
                  color = "red")
axes_mu_a[0].text(lower_wavelength - 22, 
                  epi_mu_a_sim_Jolivot_SPT_VI[lower_wavelength], 
                  "M=44.7%", 
                  color = "red")

for i in epi_mu_a_sim_set:
    axes_mu_a[0].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
            i[lower_wavelength:upper_wavelength], 
            color = "red",
            alpha = 0.08
            )

axes_mu_a[0].set_title(r"$\mu_a$ epidermis", fontsize = 18)

axes_mu_a[1].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  derm_mu_a_sim[lower_wavelength:upper_wavelength], 
                  color = "red")

mu_a_derm_660 = []
mu_a_derm_940 = []

for i in derm_mu_a_sim_set:
    axes_mu_a[1].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
            i[lower_wavelength:upper_wavelength], 
            color = "red",
            alpha = 0.02
            )
    mu_a_derm_660.append(i[660])
    mu_a_derm_940.append(i[940])
axes_mu_a[1].set_title(r"$\mu_a$ dermis",  fontsize = 18)



axes_mu_a[2].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  sub_mu_a_sim[lower_wavelength:upper_wavelength], 
                  color = "red")

mu_a_sub_650 = []
mu_a_sub_RMSEs = []
for i in sub_mu_a_sim_set:
    axes_mu_a[2].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
            i[lower_wavelength:upper_wavelength], 
            color = "red",
            alpha = 0.02
            )
    mu_a_sub_RMSE = err.RMSE(i[lower_wavelength:upper_wavelength], mu_a_sub_Ganesan)
    mu_a_sub_RMSEs.append(mu_a_sub_RMSE)
    mu_a_sub_650.append(i[650])

axes_mu_a[2].set_title(r"$\mu_a$ subcutaneous", fontsize = 18)


# Experimental absorption spectra plots
axes_mu_a[0].plot(df_Wan_epi_mu_a.loc[lower_wavelength:upper_wavelength], 
                  color = "orange", 
                  label = "Wan 1981: Dark Skin",
                  linestyle = "--")
axes_mu_a[0].plot(df_Shimojo_epi_mu_a.loc[lower_wavelength:upper_wavelength], 
                  color = "black", 
                  label = "Shimojo 2020",
                  linestyle = "--")

axes_mu_a[1].plot(df_Simpson_derm_mu_a.loc[lower_wavelength:upper_wavelength], 
                  color = "green",
                  linestyle = "--")
# axes_mu_a[1].plot(df_Shimojo_derm_mu_a.loc[lower_wavelength:upper_wavelength], 
#                   color = "black",
#                   linestyle = "--")

axes_mu_a[2].plot(df_Ganesan_sub_mu_a.loc[lower_wavelength:upper_wavelength], 
                  color = "black", 
                  label = "Ganesan 2016", 
                  linestyle = "--")
axes_mu_a[2].plot(df_Simpson_sub_mu_a.loc[lower_wavelength:upper_wavelength], 
                  color = "green", 
                  label = "Simpson 1998", 
                  linestyle = "--")
axes_mu_a[2].scatter(x = 650, 
                     y = 0.006,
                     color = "pink",
                     label = "Weber 2009",
                     s = 100)


fig_mu_a.legend(fontsize = 10)




fig_mu_s, axes_mu_s = plt.subplots(1, 3, layout = "constrained")
fig_mu_s.get_layout_engine().set(h_pad = 0.1, w_pad = 0.1)
fig_mu_s.suptitle("Reduced Scattering Spectras for Layers in Three-Layer Tissue Model", fontsize = 18)

for ax in axes_mu_s:
    ax.set_xlabel("Wavelength (nm)", fontsize = 18)
    ax.set_ylabel("Reduced Scattering (1/mm)", fontsize = 18)
    ax.tick_params(axis = "both", labelsize = 18)

# Simulated reduced scattering coefficients
axes_mu_s[0].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  epi_mu_s_red_sim[lower_wavelength:upper_wavelength], 
                  label = "Jacques 2013 Model", 
                  color = "red")
axes_mu_s[0].set_title(r"$\mu_s^'$ epidermis", fontsize = 18)



axes_mu_s[1].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  derm_mu_s_red_sim[lower_wavelength:upper_wavelength], 
                  color = "red")
axes_mu_s[1].set_title(r"$\mu^'_s$ dermis", fontsize = 18)

axes_mu_s[2].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  sub_mu_s_red_sim[lower_wavelength:upper_wavelength], 
                  color = "red")
axes_mu_s[2].set_title(r"$\mu^'_s$ subcutaneous tissue", fontsize = 18)

# Experimental reduced scattering coefficients
axes_mu_s[0].plot(df_Shimojo_epi_red_mu_s.loc[lower_wavelength:upper_wavelength], 
                  color = "cyan", 
                  label = "Shimojo 2020", 
                  linestyle = "--")

axes_mu_s[1].plot(df_Shimojo_derm_red_mu_s.loc[lower_wavelength:upper_wavelength], 
                  color = "cyan", 
                  linestyle = "--")
axes_mu_s[1].plot(df_Simpson_derm_red_mu_s.loc[lower_wavelength:upper_wavelength], 
                  color = "green", 
                  label = "Simpson 1998", 
                  linestyle = "--")

axes_mu_s[2].plot(df_Shimojo_sub_red_mu_s.loc[lower_wavelength:upper_wavelength], 
                  color = "cyan", 
                  linestyle = "--")
axes_mu_s[2].plot(df_Simpson_sub_red_mu_s.loc[lower_wavelength:upper_wavelength], 
                  color = "green",  
                  linestyle = "--")
axes_mu_s[2].plot(df_Ganesan_sub_red_mu_s.loc[lower_wavelength:upper_wavelength], 
                  color = "purple", 
                  label = "Ganesan 2016", 
                  linestyle = "--")


fig_mu_s.legend(fontsize = 10)

# print(f"Maximum mu_a_derm_660nm: {max(mu_a_derm_660)}")
# print(f"Minimum mu_a_derm_660nm: {min(mu_a_derm_660)}")
# print(f"Maximum mu_a_derm_940nm: {max(mu_a_derm_940)}")
# print(f"Minimum mu_a_derm_940nm: {min(mu_a_derm_940)}")
# print(f"Maximum mu_a_sub_650nm: {max(mu_a_sub_650)}")
# print(f"Minimum mu_a_sub_650nm: {min(mu_a_sub_650)}")
# print(f"Maximum mu_a_sub_RMSE: {max(mu_a_sub_RMSEs)}")
# print(f"Minimum mu_a_sub_RMSE: {min(mu_a_sub_RMSEs)}")

# Find tissue compostion which matches Ganesan 2016
result_Ganesan = least_squares(
    err.residual_sub, 
    [layer.sub_Ganesan["blood"],
     layer.sub_Ganesan["oxygen_saturation"],
     layer.sub_Ganesan["water"],
     layer.sub_Ganesan["fat"],
     0.25],
    args = (mu_a_sub_Ganesan, layer.mu_a_sub,),
    bounds =([sub_comp_lower[p] for p in sub_params],
             [sub_comp_upper[p] for p in sub_params])
    )

print(
    f"===Tissue Composition Prediction of Ganesan 2016===\n"
    f"Blood: {result_Ganesan.x[0]:.4f}\n"
    f"Oxygen Saturation: {result_Ganesan.x[1]:.4f}\n"
    f"Water: {result_Ganesan.x[2]:.4f}\n"
    f"Fat: {result_Ganesan.x[3]:.4f}\n"
    f"Melanin: 0.0\n"
    f"Gamma: {result_Ganesan.x[4]:.4f}\n"
    )

mu_a_sub_Ganesan_fitted = layer.mu_a_sub(
    B = result_Ganesan.x[0], 
    S = result_Ganesan.x[1],
    W = result_Ganesan.x[2],
    F = result_Ganesan.x[3],
    M = layer.sub_Ganesan["melanosome"],
    gamma = result_Ganesan.x[4]
    )

axes_mu_a[2].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  mu_a_sub_Ganesan_fitted[lower_wavelength:upper_wavelength],
                  label = "model guess",
                  color = "brown")

# Find tissue compostion which matches Simpson 1998, subcutaneous layer
# Simpson 1998 compostion not know so using Ganesan 2016 composition as an initial guess
result_Simpson_sub = least_squares(
    err.residual_sub, 
    [layer.sub_Ganesan["blood"],
     layer.sub_Ganesan["oxygen_saturation"],
     layer.sub_Ganesan["water"],
     layer.sub_Ganesan["fat"],
     0.25
     ],
    args = (mu_a_sub_Simpson, layer.mu_a_sub,),
    bounds =([sub_comp_lower[p] for p in sub_params],
             [sub_comp_upper[p] for p in sub_params])
    )

print(
    f"===Tissue Composition Prediction of Simpson 1998, Subcutaneous===\n"
    f"Blood: {result_Simpson_sub.x[0]:.4f}\n"
    f"Oxygen Saturation: {result_Simpson_sub.x[1]:.4f}\n"
    f"Water: {result_Simpson_sub.x[2]:.4f}\n"
    f"Fat: {result_Simpson_sub.x[3]:.4f}\n"
    f"Melanin: 0.0\n"
    f"Gamma: {result_Simpson_sub.x[4]:.4f}\n"
    )

mu_a_sub_Simpson_fitted = layer.mu_a_sub(
    B = result_Simpson_sub.x[0], 
    S = result_Simpson_sub.x[1],
    W = result_Simpson_sub.x[2],
    F = result_Simpson_sub.x[3],
    M = 0.0,
    gamma = result_Simpson_sub.x[4]
    )

axes_mu_a[2].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  mu_a_sub_Simpson_fitted[lower_wavelength:upper_wavelength],
                  color = "brown")


# Find tissue compostion which matches Simpson 1998, dermis layer
# Simpson 1998 compostion not know so using Jacques 2013 composition as an initial guess
result_Simpson_derm = least_squares(
    err.residual_derm, 
    [layer.derm_Jacques["blood"],
     layer.derm_Jacques["oxygen_saturation"],
     layer.derm_Jacques["water"],
     layer.derm_Jacques["gamma"],
     ],
    args = (mu_a_derm_Simpson, layer.mu_a_derm,),
    bounds =([derm_comp_lower[p] for p in derm_params],
             [derm_comp_upper[p] for p in derm_params])
    )

print(
    f"===Tissue Composition Prediction of Simpson 1998, Dermis===\n"
    f"Blood: {result_Simpson_derm.x[0]:.4f}\n"
    f"Oxygen Saturation: {result_Simpson_derm.x[1]:.4f}\n"
    f"Water: {result_Simpson_derm.x[2]:.4f}\n"
    f"Melanin: 0.0\n"
    f"Gamma: {result_Simpson_sub.x[3]:.4f}\n"
    )

mu_a_derm_Simpson_fitted = layer.mu_a_derm(
    B = result_Simpson_derm.x[0], 
    S = result_Simpson_derm.x[1],
    W = result_Simpson_derm.x[2],
    F = 0.0,
    M = 0.0,
    gamma = result_Simpson_derm.x[3]
    )

axes_mu_a[1].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  mu_a_derm_Simpson_fitted[lower_wavelength:upper_wavelength],
                  color = "brown")


# # Comment out to plot tissue spectras
# # plt.close("all")

# # =========================================
# # MCeXtreme Simulation
# # =========================================


# sim_wavelength_nm = 600

# # Format: [mu_a, mu_s, g, n]
# air_optic_params= [0.01, 0.01, 1, 1]
# epi_optic_params = [epi_mu_a_sim[sim_wavelength_nm], epi_mu_s_sim[sim_wavelength_nm], epi_g, epi_n]
# derm_optic_params = [derm_mu_a_sim[sim_wavelength_nm], derm_mu_s_sim[sim_wavelength_nm], derm_g, derm_n]
# sub_optic_params = [sub_mu_a_sim[sim_wavelength_nm], sub_mu_s_sim[sim_wavelength_nm], sub_g, sub_n]


# # Define simulation values
# prop = [
#     air_optic_params,
#     epi_optic_params,
#     derm_optic_params,
#     sub_optic_params
# ]


# Nx, Ny, Nz = 160, 160, 160

# vol = np.zeros((Nx, Ny, Nz), dtype = np.uint8)
# vol[:,:,0:3] = 1
# vol[:,:,3:13] = 2
# vol[:,:,13:160] = 3


# sim_config = {
#     "nphoton": 1e7, # 1e9: 57 seconds
#     "vol": vol,
#     "prop": prop,
#     "srcpos": [40, 30, 0],
#     "srcdir": [0, 0, 1],
#     "detdir": [0, 0, 1, 20],
#     "detpos": [[120, 30, 0, 20]],
#     "tstart": 0,
#     "tend": 10e-9,
#     "tstep": 1e-9,
#     "unitinmm": 0.1,
#     "autopilot": 1,
#     "gpuid": 1,
#     "issrcfrom0": 1,
#     "issavedet": 1,
#     "issaveseed": 1,
#     "savedetflag" : "dp"
# }

# # Run simulation


# sim_out = pmcx.run(sim_config)

# replay_config = sim_config.copy()
# replay_config["outputtype"] = "jacobian"
# replay_config["detphotons"] = sim_out["detp"]
# replay_config["seed"] = sim_out["seeds"]

# replay_out = pmcx.run(replay_config)

# jacobian_PMDF = replay_out["flux"]



# flux = sim_out["flux"]

# fig_flux, axes_flux = plt.subplots(1, 2, layout = "constrained")

# # select one time gate
# flux_1 = flux[:, :, :, 0]
# flux_1_yslice = flux_1[:, 30, :]
# axes_flux[0].set_title("Fluence Map")
# axes_flux[0].set_xlabel("x-axis (voxels)")
# axes_flux[0].set_ylabel("z-axis (voxels)")
# axes_flux[0].axhline(3, linestyle = "--", color = "black")
# axes_flux[0].axhline(13, linestyle = "--", color = "black")
# axes_flux[0].text(1, 1.5, "epidermis")
# axes_flux[0].text(1, 8, "dermis")
# axes_flux[0].text(1, 30, "subcutaneous")
# axes_flux[0].arrow(sim_config["srcpos"][0], sim_config["srcpos"][2] - 1, 0, 2, head_width = 1, head_length = 1, color = "red", label = "source")
# axes_flux[0].arrow(sim_config["detpos"][0][0], sim_config["detpos"][0][2] - 1, 0, 2, head_width = 1, head_length = 1, color = "blue", label = "detector" )

# axes_flux[0].legend()
# axes_flux[0].imshow(np.log10(np.abs(flux_1_yslice.T) + 1e-20))

# jacobian_1 = np.sum(jacobian_PMDF, axis = 3)
# jacobian_1_slice = jacobian_1[:, 30, :]
# axes_flux[1].set_title("Spatial Absorption Sensitivity")
# axes_flux[1].set_xlabel("x-axis (voxels)")
# axes_flux[1].set_ylabel("z-axis (voxels)")
# axes_flux[1].axhline(3, linestyle = "--", color = "black")
# axes_flux[1].axhline(13, linestyle = "--", color = "black")
# axes_flux[1].text(1, 1.5, "epidermis")
# axes_flux[1].text(1, 8, "dermis")
# axes_flux[1].text(1, 30, "subcutaneous")
# axes_flux[1].arrow(sim_config["srcpos"][0], sim_config["srcpos"][2] - 1, 0, 2, head_width = 1, head_length = 1, color = "red", label = "source")
# axes_flux[1].arrow(sim_config["detpos"][0][0], sim_config["detpos"][0][2] - 1, 0, 2, head_width = 1, head_length = 1, color = "blue", label = "detector" )

# axes_flux[1].legend()
# image = axes_flux[1].imshow(np.log10(np.abs(jacobian_1_slice.T) + 1e-20))

# # fig_flux.colorbar(
# #     image,
# #     ax = axes_flux[0],
# # )

# fig_flux.colorbar(
#     image,
#     ax = axes_flux[1],
# )

# PMDF = np.sum(jacobian_PMDF, axis = 3)
# PMDF_epi = np.sum(PMDF[vol == 1])
# PMDF_derm = np.sum(PMDF[vol == 2])
# PMDF_sub = np.sum(PMDF[vol == 3])

# PMDF_total = PMDF_epi + PMDF_derm + PMDF_sub

# PMDF_epi_per = 100 * (PMDF_epi / PMDF_total)
# PMDF_derm_per = 100 * (PMDF_derm / PMDF_total)
# PMDF_sub_per = 100 * (PMDF_sub / PMDF_total)

# print(f"ES(%): {PMDF_epi_per:.2f} %")
# print(f"DS(%): {PMDF_derm_per:.2f} %")
# print(f"SS(%): {PMDF_sub_per:.2f} %")




plt.show()