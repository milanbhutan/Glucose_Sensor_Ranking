import pmcx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import layer
from sklearn.metrics import mean_squared_error


norm_wavelength_nm = np.arange(0, 2001, 1)

epidermis_mu_a_sim = layer.mu_a_epi(layer.epidermis_Jacques["blood"], 
                                    layer.epidermis_Jacques["oxygen_saturation"], 
                                    layer.epidermis_Jacques["water"], 
                                    layer.epidermis_Jacques["fat"], 
                                    layer.epidermis_Jacques["melanosome"]
                                    )

epidermis_mu_a_sim_Jolivot_SPT_II = layer.mu_a_epi(layer.epidermis_Jolivot_SPT_II["blood"], 
                                    layer.epidermis_Jolivot_SPT_II["oxygen_saturation"], 
                                    layer.epidermis_Jolivot_SPT_II["water"], 
                                    layer.epidermis_Jolivot_SPT_II["fat"], 
                                    layer.epidermis_Jolivot_SPT_II["melanosome"]
                                    )

epidermis_mu_a_sim_Jolivot_SPT_III = layer.mu_a_epi(layer.epidermis_Jolivot_SPT_III["blood"], 
                                    layer.epidermis_Jolivot_SPT_III["oxygen_saturation"], 
                                    layer.epidermis_Jolivot_SPT_III["water"], 
                                    layer.epidermis_Jolivot_SPT_III["fat"], 
                                    layer.epidermis_Jolivot_SPT_III["melanosome"]
                                    )

epidermis_mu_a_sim_Jolivot_SPT_IV = layer.mu_a_epi(layer.epidermis_Jolivot_SPT_IV["blood"], 
                                    layer.epidermis_Jolivot_SPT_IV["oxygen_saturation"], 
                                    layer.epidermis_Jolivot_SPT_IV["water"], 
                                    layer.epidermis_Jolivot_SPT_IV["fat"], 
                                    layer.epidermis_Jolivot_SPT_IV["melanosome"]
                                    )

epidermis_mu_a_sim_Jolivot_SPT_V = layer.mu_a_epi(layer.epidermis_Jolivot_SPT_V["blood"], 
                                    layer.epidermis_Jolivot_SPT_V["oxygen_saturation"], 
                                    layer.epidermis_Jolivot_SPT_V["water"], 
                                    layer.epidermis_Jolivot_SPT_V["fat"], 
                                    layer.epidermis_Jolivot_SPT_V["melanosome"]
                                    )

epidermis_mu_a_sim_Jolivot_SPT_VI = layer.mu_a_epi(layer.epidermis_Jolivot_SPT_VI["blood"], 
                                    layer.epidermis_Jolivot_SPT_VI["oxygen_saturation"], 
                                    layer.epidermis_Jolivot_SPT_VI["water"], 
                                    layer.epidermis_Jolivot_SPT_VI["fat"], 
                                    layer.epidermis_Jolivot_SPT_VI["melanosome"]
                                    )

dermis_mu_a_sim = layer.mu_a_derm(layer.dermis_Jacques["blood"], 
                                  layer.dermis_Jacques["oxygen_saturation"], 
                                  layer.dermis_Jacques["water"], 
                                  layer.dermis_Jacques["fat"], 
                                  layer.dermis_Jacques["melanosome"]
                                  )

lower_dermis_mu_a_sim_Jonasson = layer.mu_a_derm(layer.lower_dermis_Jonasson["blood"], 
                                  layer.lower_dermis_Jonasson["oxygen_saturation"], 
                                  layer.lower_dermis_Jonasson["water"], 
                                  layer.lower_dermis_Jonasson["fat"], 
                                  layer.lower_dermis_Jonasson["melanosome"]
                                  )

subcutaneous_mu_a_sim = layer.mu_a_sub(layer.subcutaneous_Jacques["blood"], 
                                       layer.subcutaneous_Jacques["oxygen_saturation"], 
                                       layer.subcutaneous_Jacques["water"], 
                                       layer.subcutaneous_Jacques["fat"], 
                                       layer.subcutaneous_Jacques["melanosome"]
                                       )

subcutaneous_mu_a_sim_Ganesan = layer.mu_a_sub(layer.subcutaneous_Ganesan["blood"], 
                                               layer.subcutaneous_Ganesan["oxygen_saturation"],
                                               layer.subcutaneous_Ganesan["water"],
                                               layer.subcutaneous_Ganesan["fat"],
                                               layer.subcutaneous_Ganesan["melanosome"]
                                               )



# For scattering coefficent
epidermis_g = 0.9
dermis_g = 0.9
subcutaneous_g = 0.9

epidermis_mu_s_sim = layer.mu_s_layer("epidermis", epidermis_g)
dermis_mu_s_sim = layer.mu_s_layer("dermis", dermis_g)
subcutaneous_mu_s_sim = layer.mu_s_layer("subcutaneous", subcutaneous_g)

# For reduced scattering coefficent
epidermis_g_red = 0
dermis_g_red = 0
subcutaneous_g_red = 0

epidermis_mu_s_red_sim = layer.mu_s_layer("epidermis", epidermis_g_red)
dermis_mu_s_red_sim = layer.mu_s_layer("dermis", dermis_g_red)
subcutaneous_mu_s_red_sim = layer.mu_s_layer("subcutaneous", subcutaneous_g_red)


epidermis_n = 1.4
dermis_n = 1.4
subcutaneous_n = 1.4


# Obtain data for experimental absorption and scattering spectras for three-layer tissue model from Shimojo 2020

# Absorption Data
df_Shimojo_epi_mu_a = (
    pd.read_csv("Exp_Absorption/shimojo_epidermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
epidermis_mu_a_exp = df_Shimojo_epi_mu_a["coefficient_mm^-1"]

df_Salomatina_epi_mu_a = (
    pd.read_csv("Exp_Absorption/salomatina_epidermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Shimojo_derm_mu_a = (
    pd.read_csv("Exp_Absorption/shimojo_dermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
dermis_mu_a_exp = df_Shimojo_derm_mu_a["coefficient_mm^-1"]

df_Jonasson_derm_mu_a = (
    pd.read_csv("Exp_Absorption/jonasson_lower_dermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Salomatina_derm_mu_a = (
    pd.read_csv("Exp_Absorption/salomatina_dermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Shimojo_sub_mu_a = (
    pd.read_csv("Exp_Absorption/shimojo_subcutaneous_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
subcutaneous_mu_a_exp = df_Shimojo_sub_mu_a["coefficient_mm^-1"]

df_Ganesan_2016_mu_a = (
    pd.read_csv("Exp_Absorption/ganesan_subcutaneous_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
gan_subcutaneous_mu_a_exp = df_Ganesan_2016_mu_a["mu_a_mm^-1"]

df_Salomatina_sub_mu_a = (
    pd.read_csv("Exp_Absorption/salomatina_subcutaneous_fat_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

# Scattering Data
df_Shimojo_epi_mu_s = (
    pd.read_csv("Exp_Scattering/shimojo_epidermis_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
epidermis_mu_s_exp = df_Shimojo_epi_mu_s["coefficient_mm^-1"]

df_Salomatina_epi_mu_s = (
    pd.read_csv("Exp_Scattering/salomatina_2006_epidermis_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Shimojo_derm_mu_s = (
    pd.read_csv("Exp_Scattering/shimojo_dermis_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
dermis_mu_s_exp = df_Shimojo_derm_mu_s["coefficient_mm^-1"]

df_Salomatina_derm_mu_s = (
    pd.read_csv("Exp_Scattering/salomatina_2006_dermis_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Shimojo_sub_mu_s = (
    pd.read_csv("Exp_Scattering/shimojo_subcutaneous_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
subcutaneous_mu_s_exp = df_Shimojo_sub_mu_s["coefficient_mm^-1"]

df_Ganesan_2016_mu_s = (
    pd.read_csv("Exp_Scattering/ganesan_subcutaneous_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
gan_subcutaneous_mu_s_exp = df_Ganesan_2016_mu_s["mu_s_prime_mm^-1"]

df_Salomatina_sub_mu_s = (
    pd.read_csv("Exp_Scattering/salomatina_2006_subcutaneous_fat_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

lower_wavelength = 600
upper_wavelength = 1060

# Calculating RMSE
epidermis_mu_a_exp_align, epidermis_mu_a_sim_align = epidermis_mu_a_exp.align(epidermis_mu_a_sim[lower_wavelength:upper_wavelength], join = "inner")
dermis_mu_a_exp_align, dermis_mu_a_sim_align = dermis_mu_a_exp.align(dermis_mu_a_sim[lower_wavelength:upper_wavelength], join = "inner")
subcutaneous_mu_a_exp_align, subcutaneous_mu_a_sim_align = subcutaneous_mu_a_exp.align(subcutaneous_mu_a_sim[lower_wavelength:upper_wavelength], join = "inner")
epidermis_mu_s_exp_align, epidermis_mu_s_red_sim_align = epidermis_mu_s_exp.align(epidermis_mu_s_red_sim[lower_wavelength:upper_wavelength], join = "inner")
dermis_mu_s_exp_align, dermis_mu_s_red_sim_align = dermis_mu_s_exp.align(dermis_mu_s_red_sim[lower_wavelength:upper_wavelength], join = "inner")
subcutaneous_mu_s_exp_align, subcutaneous_mu_s_red_sim_align = subcutaneous_mu_s_exp.align(subcutaneous_mu_s_red_sim[lower_wavelength:upper_wavelength], join = "inner")

epi_mu_a_RMSE = np.sqrt(mean_squared_error(epidermis_mu_a_exp_align, epidermis_mu_a_sim_align))
derm_mu_a_RMSE = np.sqrt(mean_squared_error(dermis_mu_a_exp_align, dermis_mu_a_sim_align))
sub_mu_a_RMSE = np.sqrt(mean_squared_error(subcutaneous_mu_a_exp_align, subcutaneous_mu_a_sim_align))
epi_mu_s_RMSE = np.sqrt(mean_squared_error(epidermis_mu_s_exp_align, epidermis_mu_s_red_sim_align))
derm_mu_s_RMSE = np.sqrt(mean_squared_error(dermis_mu_s_exp_align, dermis_mu_s_red_sim_align))
sub_mu_s_RMSE = np.sqrt(mean_squared_error(subcutaneous_mu_s_exp_align, subcutaneous_mu_s_red_sim_align))

epi_mu_a_NRMSE = epi_mu_a_RMSE / (epidermis_mu_a_exp_align.max() - epidermis_mu_a_exp_align.min()) * 100
derm_mu_a_NRMSE = derm_mu_a_RMSE / (dermis_mu_a_exp_align.max() - dermis_mu_a_exp_align.min()) * 100
sub_mu_a_NRMSE = sub_mu_a_RMSE / (subcutaneous_mu_a_exp_align.max() - subcutaneous_mu_a_exp_align.min()) * 100
epi_mu_s_NRMSE = epi_mu_s_RMSE / (epidermis_mu_s_exp_align.max() - epidermis_mu_s_exp_align.min()) * 100
derm_mu_s_NRMSE = derm_mu_s_RMSE / (dermis_mu_s_exp_align.max() - dermis_mu_s_exp_align.min()) * 100
sub_mu_s_NRMSE = sub_mu_s_RMSE / (subcutaneous_mu_s_exp_align.max() - subcutaneous_mu_s_exp_align.min()) * 100

print(f"Epidermis Absorption NRMSE: {epi_mu_a_NRMSE} %")
print(f"Dermis Absorption NRMSE: {derm_mu_a_NRMSE} %")
print(f"Subcutaneous Absorption NRMSE: {sub_mu_a_NRMSE} %")
print(f"Epidermis Reduced Scattering NRMSE: {epi_mu_s_NRMSE} %")
print(f"Dermis Reduced Scattering NRMSE: {derm_mu_s_NRMSE} %")
print(f"Subcutaneous Reduced Scattering NRMSE: {sub_mu_s_NRMSE} %")


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
                  epidermis_mu_a_sim[lower_wavelength:upper_wavelength], 
                  label = "Analytical Models", 
                  color = "red")
axes_mu_a[0].text(lower_wavelength, 
                  epidermis_mu_a_sim[lower_wavelength], 
                  "M=2.5%", 
                  color = "red")

axes_mu_a[0].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  epidermis_mu_a_sim_Jolivot_SPT_II[lower_wavelength:upper_wavelength], 
                  color = "red")
axes_mu_a[0].text(lower_wavelength, 
                  epidermis_mu_a_sim_Jolivot_SPT_II[lower_wavelength], 
                  "SPT II:\nM=10.1%", 
                  color = "red")

axes_mu_a[0].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  epidermis_mu_a_sim_Jolivot_SPT_III[lower_wavelength:upper_wavelength], 
                  color = "red")
axes_mu_a[0].text(lower_wavelength, 
                  epidermis_mu_a_sim_Jolivot_SPT_III[lower_wavelength], 
                  "SPT III:\nM=18.7%", 
                  color = "red")

axes_mu_a[0].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  epidermis_mu_a_sim_Jolivot_SPT_IV[lower_wavelength:upper_wavelength], 
                  color = "red")
axes_mu_a[0].text(lower_wavelength, 
                  epidermis_mu_a_sim_Jolivot_SPT_IV[lower_wavelength], 
                  "SPT IV:\nM=26.5%", 
                  color = "red")

axes_mu_a[0].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  epidermis_mu_a_sim_Jolivot_SPT_V[lower_wavelength:upper_wavelength], 
                  color = "red")
axes_mu_a[0].text(lower_wavelength, 
                  epidermis_mu_a_sim_Jolivot_SPT_V[lower_wavelength], 
                  "SPT V:\nM=35.8%", 
                  color = "red")

axes_mu_a[0].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  epidermis_mu_a_sim_Jolivot_SPT_VI[lower_wavelength:upper_wavelength], 
                  color = "red")
axes_mu_a[0].text(lower_wavelength, 
                  epidermis_mu_a_sim_Jolivot_SPT_VI[lower_wavelength], 
                  "SPT VI:\nM=44.7%", 
                  color = "red")

axes_mu_a[0].set_title(r"$\mu_a$ epidermis", fontsize = 18)


axes_mu_a[1].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  dermis_mu_a_sim[lower_wavelength:upper_wavelength], 
                  color = "red")
axes_mu_a[1].text(lower_wavelength, 
                  dermis_mu_a_sim[lower_wavelength], 
                  "Composition: Jacques", 
                  color = "red")

axes_mu_a[1].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  lower_dermis_mu_a_sim_Jonasson[lower_wavelength:upper_wavelength], 
                  color = "red")
axes_mu_a[1].text(lower_wavelength, 
                  lower_dermis_mu_a_sim_Jonasson[lower_wavelength], 
                  "Composition: Jonasson", 
                  color = "red")

axes_mu_a[1].set_title(r"$\mu_a$ dermis",  fontsize = 18)


axes_mu_a[2].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  subcutaneous_mu_a_sim[lower_wavelength:upper_wavelength], 
                  color = "red")
axes_mu_a[2].text(lower_wavelength, 
                  subcutaneous_mu_a_sim[lower_wavelength], 
                  "Composition: Jacques", 
                  color = "red")

axes_mu_a[2].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], 
                  subcutaneous_mu_a_sim_Ganesan[lower_wavelength:upper_wavelength], 
                  color = "red")
axes_mu_a[2].text(lower_wavelength, 
                  subcutaneous_mu_a_sim_Ganesan[lower_wavelength], 
                  "Composition: Ganesan", 
                  color = "red")

axes_mu_a[2].set_title(r"$\mu_a$ subcutaneous", fontsize = 18)


# Experimental absorption spectra plots
axes_mu_a[0].plot(df_Shimojo_epi_mu_a.loc[lower_wavelength:upper_wavelength], 
                  label = "Shimojo 2020", 
                  color = "orange", 
                  linestyle = "--")
axes_mu_a[0].plot(df_Salomatina_epi_mu_a.loc[lower_wavelength:upper_wavelength], 
                  label = "Salomatina 2006", 
                  color = "green", 
                  linestyle = "--")


axes_mu_a[1].plot(df_Shimojo_derm_mu_a.loc[lower_wavelength:upper_wavelength], 
                  color = "orange", 
                  linestyle = "--")
axes_mu_a[1].plot(df_Jonasson_derm_mu_a.loc[lower_wavelength:upper_wavelength], 
                  color = "purple", 
                  label = "Jonasson 2023", 
                  linestyle = "--")
axes_mu_a[1].plot(df_Salomatina_derm_mu_a.loc[lower_wavelength:upper_wavelength], 
                  color = "green", 
                  linestyle = "--")


axes_mu_a[2].plot(df_Shimojo_sub_mu_a.loc[lower_wavelength:upper_wavelength], 
                  color = "orange", 
                  linestyle = "--")
axes_mu_a[2].plot(df_Ganesan_2016_mu_a.loc[lower_wavelength:upper_wavelength], 
                  color = "blue", 
                  label = "Ganesan 2016", 
                  linestyle = "--")
axes_mu_a[2].plot(df_Salomatina_sub_mu_a.loc[lower_wavelength:upper_wavelength], 
                  color = "green", 
                  linestyle = "--")

fig_mu_a.legend(fontsize = 10)




fig_mu_s, axes_mu_s = plt.subplots(1, 3, layout = "constrained")
fig_mu_s.get_layout_engine().set(h_pad = 0.1, w_pad = 0.1)
fig_mu_s.suptitle("Reduced Scattering Spectras for Layers in Three-Layer Tissue Model", fontsize = 18)

for ax in axes_mu_s:
    ax.set_xlabel("Wavelength (nm)", fontsize = 18)
    ax.set_ylabel("Reduced Scattering (1/mm)", fontsize = 18)
    ax.tick_params(axis = "both", labelsize = 18)

# Simulated reduced scattering coefficients
axes_mu_s[0].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], epidermis_mu_s_red_sim[lower_wavelength:upper_wavelength], label = "Jacques 2013 Model", color = "red")
axes_mu_s[0].set_title(r"$\mu_s^'$ epidermis", fontsize = 18)

axes_mu_s[1].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], dermis_mu_s_red_sim[lower_wavelength:upper_wavelength], color = "red")
axes_mu_s[1].set_title(r"$\mu^'_s$ dermis", fontsize = 18)

axes_mu_s[2].plot(norm_wavelength_nm[lower_wavelength:upper_wavelength], subcutaneous_mu_s_red_sim[lower_wavelength:upper_wavelength], color = "red")
axes_mu_s[2].set_title(r"$\mu^'_s$ subcutaneous tissue", fontsize = 18)

# Experimental reduced scattering coefficients
axes_mu_s[0].plot(df_Shimojo_epi_mu_s.loc[lower_wavelength:upper_wavelength], label = "Shimojo 2020",color = "orange", linestyle = "--")
axes_mu_s[0].plot(df_Salomatina_epi_mu_s.loc[lower_wavelength:upper_wavelength], label = "Salomatina 2006", color = "green", linestyle = "--")

axes_mu_s[1].plot(df_Shimojo_derm_mu_s.loc[lower_wavelength:upper_wavelength],color = "orange", linestyle = "--")
axes_mu_s[1].plot(df_Salomatina_derm_mu_s.loc[lower_wavelength:upper_wavelength], color = "green", linestyle = "--")

axes_mu_s[2].plot(df_Shimojo_sub_mu_s.loc[lower_wavelength:upper_wavelength], color = "orange", linestyle = "--")
axes_mu_s[2].plot(df_Ganesan_2016_mu_s.loc[lower_wavelength:upper_wavelength],color = "purple", label = "Ganesan 2016", linestyle = "--")
axes_mu_s[2].plot(df_Salomatina_sub_mu_s.loc[lower_wavelength:upper_wavelength], color = "green", linestyle = "--")

fig_mu_s.legend(fontsize = 10)




# Comment out to plot tissue spectras
# plt.close("all")

# =========================================
# MCeXtreme Simulation
# =========================================


sim_wavelength_nm = 1000

# Format: [mu_a, mu_s, g, n]
air_optic_params= [0.01, 0.01, 1, 1]
epidermis_optic_params = [epidermis_mu_a_sim[sim_wavelength_nm], epidermis_mu_s_sim[sim_wavelength_nm], epidermis_g, epidermis_n]
dermis_optic_params = [dermis_mu_a_sim[sim_wavelength_nm], dermis_mu_s_sim[sim_wavelength_nm], dermis_g, dermis_n]
subcutaneous_optic_params = [subcutaneous_mu_a_sim[sim_wavelength_nm], subcutaneous_mu_s_sim[sim_wavelength_nm], subcutaneous_g, subcutaneous_n]


# Define simulation values
prop = [
    air_optic_params,
    epidermis_optic_params,
    dermis_optic_params,
    subcutaneous_optic_params
]


Nx, Ny, Nz = 60, 60, 60

vol = np.zeros([Nx, Ny, Nz])
vol[:,:,0:3] = 1
vol[:,:,3:13] = 2
vol[:,:,13:60] = 3


sim_config = {
    "nphoton": 1e7, # 1e9: 57 seconds
    "vol": vol,
    "prop": prop,
    "srcpos": [30, 30, 0],
    "srcdir": [0, 0, 1],
    "detdir": [0, 0, 1, 1],
    "detpos": [[50, 30, 0, 10]],
    "tstart": 0,
    "tend": 10e-9,
    "tstep": 1e-9,
    "unitinmm": 0.1,
    "autopilot": 1,
    "gpuid": 1,
    "issrcfrom0": 1,
    "issavedet": 1,
    "issaveseed": 1,
    "savedetflag" : "dp"
}

# Run simulation


sim_out = pmcx.run(sim_config)

replay_config = sim_config.copy()
replay_config["outputtype"] = "jacobian"
replay_config["detphotons"] = sim_out["detp"]
replay_config["seed"] = sim_out["seeds"]
replay_config["debuglevel"] = "M"

replay_out = pmcx.run(replay_config)



jacobian_traj = replay_out["traj"]
jacobian_PMDF = replay_out["flux"]

PMDF = np.sum(jacobian_PMDF, axis = 3)
PMDF_epi = np.sum(PMDF[vol == 1])
PMDF_derm = np.sum(PMDF[vol == 2])
PMDF_sub = np.sum(PMDF[vol == 3])

PMDF_total = PMDF_epi + PMDF_derm + PMDF_sub

PMDF_epi_per = 100 * (PMDF_epi / PMDF_total)
PMDF_derm_per = 100 * (PMDF_derm / PMDF_total)
PMDF_sub_per = 100 * (PMDF_sub / PMDF_total)

print(f"{PMDF_epi_per:.2f} %")
print(f"{PMDF_derm_per:.2f} %")
print(f"{PMDF_sub_per:.2f} %")

# avg_path = pmcx.meanpath(sim_out["detp"], np.asarray(prop))
# print(avg_path)

flux = sim_out["flux"]

fig_flux, axes_flux = plt.subplots(1, 2, layout = "constrained")

# select one time gate
flux_1 = flux[:, :, :, 0]
flux_1_yslice = flux_1[:, 30, :]
axes_flux[0].set_title("0-1 ns timegate")
axes_flux[0].set_xlabel("x-axis (voxels)")
axes_flux[0].set_ylabel("z-axis (voxels)")
axes_flux[0].axhline(3, linestyle = "--", color = "black")
axes_flux[0].axhline(13, linestyle = "--", color = "black")
axes_flux[0].text(1, 1.5, "epidermis")
axes_flux[0].text(1, 8, "dermis")
axes_flux[0].text(1, 30, "subcutaneous")
axes_flux[0].arrow(sim_config["srcpos"][0], sim_config["srcpos"][2] - 1, 0, 2, head_width = 1, head_length = 1, color = "red", label = "source")
axes_flux[0].arrow(sim_config["detpos"][0][0], sim_config["detpos"][0][2] - 1, 0, 2, head_width = 1, head_length = 1, color = "blue", label = "detector" )

axes_flux[0].legend()
axes_flux[0].imshow(np.log10(np.abs(flux_1_yslice.T) + 1e-20))

jacobian_1 = jacobian_PMDF[:, :, :, 0]
jacobian_1_slice = jacobian_1[:, 30, :]

# jacobian_y = jacobian_y * (jacobian_y == 0)

axes_flux[1].imshow(np.log10(np.abs(jacobian_1_slice.T) + 1e-20))





























plt.show()