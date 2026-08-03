import pmcx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import layer
from sklearn.metrics import mean_squared_error


# From Jacques 2013
epidermis = {
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0.20,
    "fat" : 0,
    "melanin" : 0.025
}

dermis = {
    "blood" : 0.0020,
    "oxygen_saturation" : 0.39,
    "water" : 0.65,
    "fat" : 0,
    "melanin" : 0.001
}

subcutaneous = {
    "blood" : 0.00634,
    "oxygen_saturation" : 0.60,
    "water" : 0.2184,
    "fat" : 0.6433,
    "melanin" : 0.001
}

norm_wavelength_nm = np.arange(0, 2001, 1)

epidermis_mu_a_sim = layer.mu_a_layer(epidermis["blood"], epidermis["oxygen_saturation"], epidermis["water"], epidermis["fat"], epidermis["melanin"])
dermis_mu_a_sim = layer.mu_a_layer(dermis["blood"], dermis["oxygen_saturation"], dermis["water"], dermis["fat"], dermis["melanin"])
subcutaneous_mu_a_sim = layer.mu_a_layer(subcutaneous["blood"], subcutaneous["oxygen_saturation"], subcutaneous["water"], subcutaneous["fat"], subcutaneous["melanin"])

# For scattering coefficent
# epidermis_g = 0.9
# dermis_g = 0.9
# subcutaneous_g = 0.9

# For reduced scattering coefficent
epidermis_g = 0
dermis_g = 0
subcutaneous_g = 0

epidermis_mu_s_sim = layer.mu_s_layer("epidermis", epidermis_g)
dermis_mu_s_sim = layer.mu_s_layer("dermis", dermis_g)
subcutaneous_mu_s_sim = layer.mu_s_layer("subcutaneous", subcutaneous_g)

epidermis_n = 1.4
dermis_n = 1.4
subcutaneous_n = 1.4

# Obtain data for experimental absorption and scattering spectras for three-layer tissue model from Shimojo 2020

df_Shimojo_epi_mu_a = (
    pd.read_csv("Exp_Absorption/epidermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
epidermis_mu_a_exp = df_Shimojo_epi_mu_a["coefficient_mm^-1"]

df_Shimojo_derm_mu_a = (
    pd.read_csv("Exp_Absorption/dermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
dermis_mu_a_exp = df_Shimojo_derm_mu_a["coefficient_mm^-1"]

df_Shimojo_sub_mu_a = (
    pd.read_csv("Exp_Absorption/subcutaneous_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
subcutaneous_mu_a_exp = df_Shimojo_sub_mu_a["coefficient_mm^-1"]

df_Shimojo_epi_mu_s = (
    pd.read_csv("Exp_Scattering/epidermis_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
epidermis_mu_s_exp = df_Shimojo_epi_mu_s["coefficient_mm^-1"]

df_Shimojo_derm_mu_s = (
    pd.read_csv("Exp_Scattering/dermis_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
dermis_mu_s_exp = df_Shimojo_derm_mu_s["coefficient_mm^-1"]

df_Shimojo_sub_mu_s = (
    pd.read_csv("Exp_Scattering/subcutaneous_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
subcutaneous_mu_s_exp = df_Shimojo_sub_mu_s["coefficient_mm^-1"]



# Calculating RMSE
epidermis_mu_a_exp_align, epidermis_mu_a_sim_align = epidermis_mu_a_exp.align(epidermis_mu_a_sim[600:1060], join = "inner")
dermis_mu_a_exp_align, dermis_mu_a_sim_align = dermis_mu_a_exp.align(dermis_mu_a_sim[600:1060], join = "inner")
subcutaneous_mu_a_exp_align, subcutaneous_mu_a_sim_align = subcutaneous_mu_a_exp.align(subcutaneous_mu_a_sim[600:1060], join = "inner")
epidermis_mu_s_exp_align, epidermis_mu_s_sim_align = epidermis_mu_s_exp.align(epidermis_mu_s_sim[600:1060], join = "inner")
dermis_mu_s_exp_align, dermis_mu_s_sim_align = dermis_mu_s_exp.align(dermis_mu_s_sim[600:1060], join = "inner")
subcutaneous_mu_s_exp_align, subcutaneous_mu_s_sim_align = subcutaneous_mu_s_exp.align(subcutaneous_mu_s_sim[600:1060], join = "inner")

print(epidermis_mu_a_exp_align)
print("====")
print(epidermis_mu_a_sim_align)

epi_mu_a_RMSE = np.sqrt(mean_squared_error(epidermis_mu_a_exp_align, epidermis_mu_a_sim_align))
derm_mu_a_RMSE = np.sqrt(mean_squared_error(dermis_mu_a_exp_align, dermis_mu_a_sim_align))
sub_mu_a_RMSE = np.sqrt(mean_squared_error(subcutaneous_mu_a_exp_align, subcutaneous_mu_a_sim_align))
epi_mu_s_RMSE = np.sqrt(mean_squared_error(epidermis_mu_s_exp_align, epidermis_mu_s_sim_align))
derm_mu_s_RMSE = np.sqrt(mean_squared_error(dermis_mu_s_exp_align, dermis_mu_s_sim_align))
sub_mu_s_RMSE = np.sqrt(mean_squared_error(subcutaneous_mu_s_exp_align, subcutaneous_mu_s_sim_align))

print(f"{epi_mu_a_RMSE} (1/mm)")
print(f"{derm_mu_a_RMSE} (1/mm)")
print(f"{sub_mu_a_RMSE} (1/mm)")
print(f"{epi_mu_s_RMSE} (1/mm)")
print(f"{derm_mu_s_RMSE} (1/mm)")
print(f"{sub_mu_s_RMSE} (1/mm)")




# Note: This subcutaneous spectra is only valid within 429-1100 nm where reliable tissue constituent abosroption data is availible
fig_mu_a, axes_mu_a = plt.subplots(2, 3, layout = "constrained")
fig_mu_a.get_layout_engine().set(h_pad = 0.1, w_pad = 0.1)
fig_mu_a.suptitle("Absorption Spectras for Layers in Three-Layer Tissue Model")

# Simulated absorption spectra plots
axes_mu_a[0, 0].plot(norm_wavelength_nm[600:1060], epidermis_mu_a_sim[600:1060])
axes_mu_a[0, 0].set_xlabel("Wavelength (nm)")
axes_mu_a[0, 0].set_ylabel("Absorption (1/mm)")
axes_mu_a[0, 0].set_title(r"$\mu_a$ epidermis")

axes_mu_a[0, 1].plot(norm_wavelength_nm[600:1060], dermis_mu_a_sim[600:1060])
axes_mu_a[0, 1].set_xlabel("Wavelength (nm)")
axes_mu_a[0, 1].set_ylabel("Absorption (1/mm)")
axes_mu_a[0, 1].set_title(r"$\mu_a$ dermis")

axes_mu_a[0, 2].plot(norm_wavelength_nm[600:1060], subcutaneous_mu_a_sim[600:1060])
axes_mu_a[0, 2].set_xlabel("Wavelength (nm)")
axes_mu_a[0, 2].set_ylabel("Absorption (1/mm)")
axes_mu_a[0, 2].set_title(r"$\mu_a$ subcutaneous")


# Experimental absorption spectra plots
axes_mu_a[1, 0].plot(df_Shimojo_epi_mu_a.loc[600:1060])
axes_mu_a[1, 0].set_xlabel("Wavelength (nm)")
axes_mu_a[1, 0].set_ylabel("Absorption (1/mm)")
axes_mu_a[1, 0].set_title(r"$\mu_a$ epidermis")

axes_mu_a[1, 1].plot(df_Shimojo_derm_mu_a.loc[600:1060])
axes_mu_a[1, 1].set_xlabel("Wavelength (nm)")
axes_mu_a[1, 1].set_ylabel("Absorption (1/mm)")
axes_mu_a[1, 1].set_title(r"$\mu_a$ dermis")

axes_mu_a[1, 2].plot(df_Shimojo_sub_mu_a.loc[600:1060])
axes_mu_a[1, 2].set_xlabel("Wavelength (nm)")
axes_mu_a[1, 2].set_ylabel("Absorption (1/mm)")
axes_mu_a[1, 2].set_title(r"$\mu_a$ subcutaneous")


fig_mu_s, axes_mu_s = plt.subplots(2, 3, layout = "constrained")
fig_mu_s.get_layout_engine().set(h_pad = 0.1, w_pad = 0.1)
fig_mu_s.suptitle("Reduced Scattering Spectras for Layers in Three-Layer Tissue Model")

# Simulated reduced scattering coefficients
axes_mu_s[0, 0].plot(norm_wavelength_nm[600:1060], epidermis_mu_s_sim[600:1060])
axes_mu_s[0, 0].set_xlabel("Wavelength (nm)")
axes_mu_s[0, 0].set_ylabel("Reduced Scattering (1/mm)")
axes_mu_s[0, 0].set_title(r"$\mu_s$ epidermis")

axes_mu_s[0, 1].plot(norm_wavelength_nm[600:1060], dermis_mu_s_sim[600:1060])
axes_mu_s[0, 1].set_xlabel("Wavelength (nm)")
axes_mu_s[0, 1].set_ylabel("Reduced Scattering (1/mm)")
axes_mu_s[0, 1].set_title(r"$\mu^'_s$ dermis")

axes_mu_s[0, 2].plot(norm_wavelength_nm[600:1060], subcutaneous_mu_s_sim[600:1060])
axes_mu_s[0, 2].set_xlabel("Wavelength (nm)")
axes_mu_s[0, 2].set_ylabel("Reduced Scattering (1/mm)")
axes_mu_s[0, 2].set_title(r"$\mu^'_s$ subcutaneous tissue")

# Experimental reduced scattering coefficients
axes_mu_s[1, 0].plot(df_Shimojo_epi_mu_s.loc[600:1060])
axes_mu_s[1, 0].set_xlabel("Wavelength (nm)")
axes_mu_s[1, 0].set_ylabel("Reduced Scattering (1/mm)")
axes_mu_s[1, 0].set_title(r"$\mu^'_s$ epidermis")

axes_mu_s[1, 1].plot(df_Shimojo_derm_mu_s.loc[600:1060])
axes_mu_s[1, 1].set_xlabel("Wavelength (nm)")
axes_mu_s[1, 1].set_ylabel("Reduced Scattering (1/mm)")
axes_mu_s[1, 1].set_title(r"$\mu^'_s$ dermis")

axes_mu_s[1, 2].plot(df_Shimojo_sub_mu_s.loc[600:1060])
axes_mu_s[1, 2].set_xlabel("Wavelength (nm)")
axes_mu_s[1, 2].set_ylabel("Reduced Scattering (1/mm)")
axes_mu_s[1, 2].set_title(r"$\mu^'_s$ subcutaneous")




# Comment out to plot tissue spectras
# plt.close("all")


sim_wavelength_nm = 660

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


sim_config={
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
    "gpuid": 1
}

#Run simulation
sim_out = pmcx.run(sim_config)

flux = sim_out["flux"]

fig_flux, axes_flux = plt.subplots(1, 1, layout = "constrained")

# select one time gate
flux_1 = flux[:, :, :, 0]
flux_1_yslice = flux_1[:, 30, :]
axes_flux.set_title("0-1 ns timegate")
axes_flux.set_xlabel("x-axis (voxels)")
axes_flux.set_ylabel("z-axis (voxels)")
axes_flux.axhline(3, linestyle = "--", color = "black")
axes_flux.axhline(13, linestyle = "--", color = "black")
axes_flux.text(1, 1.5, "epidermis")
axes_flux.text(1, 8, "dermis")
axes_flux.text(1, 30, "subcutaneous")
axes_flux.arrow(sim_config["srcpos"][0], sim_config["srcpos"][2] - 1, 0, 2, head_width = 1, head_length = 1, color = "red", label = "source")
axes_flux.arrow(sim_config["detpos"][0][0], sim_config["detpos"][0][2] - 1, 0, 2, head_width = 1, head_length = 1, color = "blue", label = "detector" )

axes_flux.legend()
axes_flux.imshow(np.log10(flux_1_yslice.T))

plt.show()