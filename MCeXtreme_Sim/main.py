import pmcx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import layer




# From Jacques 2013
epidermis = {
    "blood" : 0,
    "blood oxygen" : 0,
    "water" : 0.20,
    "fat" : 0,
    "melanin" : 0.025
}

dermis = {
    "blood" : 0.0020, # was 0.0020
    "oxygen_saturation" : 0.39,
    "water" : 0.65,
    "fat" : 0,
    "melanin" : 0
}

subcutaneous = {
    "blood" : 0.00634, # was 0.00634
    "oxygen_saturation" : 0.60,
    "water" : 0.2184,
    "fat" : 0.6433,
    "melanin" : 0
}

norm_wavelength_nm = np.arange(0, 2001, 1)

epidermis_mu_a = layer.mu_a_layer(epidermis["blood"], epidermis["blood oxygen"], epidermis["water"], epidermis["fat"], epidermis["melanin"])
dermis_mu_a = layer.mu_a_layer(dermis["blood"], dermis["oxygen_saturation"], dermis["water"], dermis["fat"], dermis["melanin"])
subcutaneous_mu_a = layer.mu_a_layer(subcutaneous["blood"], subcutaneous["oxygen_saturation"], subcutaneous["water"], subcutaneous["fat"], subcutaneous["melanin"])

epidermis_g = 0.9
dermis_g = 0.9
subcutaneous_g = 0.9

epidermis_mu_s = layer.mu_s_layer("epidermis", epidermis_g)
dermis_mu_s = layer.mu_s_layer("dermis", dermis_g)
subcutaneous_mu_s = layer.mu_s_layer("subcutaneous", subcutaneous_g)

epidermis_n = 1.4
dermis_n = 1.4
subcutaneous_n = 1.4

# Obtain data for experimental absorption and scattering spectras for three-layer tissue model from Shimojo 2020

df_Shimojo_epi = (
    pd.read_csv("Exp_Absorption/epidermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
    )
Shimojo_wavelength_nm = df_Shimojo_epi.index
Shimojo_mu_a_epidermis = df_Shimojo_epi["coefficient_mm^-1"]

df_Shimojo_derm = (
    pd.read_csv("Exp_Absorption/dermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
    )
Shimojo_mu_a_dermis = df_Shimojo_derm["coefficient_mm^-1"]

df_Shimojo_sub = (
    pd.read_csv("Exp_Absorption/subcutaneous_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
    )
Shimojo_mu_a_subcutaneous = df_Shimojo_sub["coefficient_mm^-1"]

# Note: This subcutaneous spectra is only valid within 429-1100 nm where reliable tissue constituent abosroption data is availible
fig_mu_a, axes_mu_a = plt.subplots(2, 3, layout = "constrained")
fig_mu_a.get_layout_engine().set(h_pad = 0.1, w_pad = 0.1)
fig_mu_a.suptitle("Absorption Spectras for Layers in Three-Layer Tissue Model")

# Simulated absorption spectra plots
axes_mu_a[0, 0].plot(norm_wavelength_nm[429:1100], epidermis_mu_a[429:1100])
axes_mu_a[0, 0].set_xlabel("Wavelength (nm)")
axes_mu_a[0, 0].set_ylabel("Absorption (1/mm)")
axes_mu_a[0, 0].set_title(r"$\mu_a$ epidermis")

axes_mu_a[0, 1].plot(norm_wavelength_nm[429:1100], dermis_mu_a[429:1100])
axes_mu_a[0, 1].set_xlabel("Wavelength (nm)")
axes_mu_a[0, 1].set_ylabel("Absorption (1/mm)")
axes_mu_a[0, 1].set_title(r"$\mu_a$ dermis")

axes_mu_a[0, 2].plot(norm_wavelength_nm[429:1100], subcutaneous_mu_a[429:1100])
axes_mu_a[0, 2].set_xlabel("Wavelength (nm)")
axes_mu_a[0, 2].set_ylabel("Absorption (1/mm)")
axes_mu_a[0, 2].set_title(r"$\mu_a$ subcutaneous")


# Experimental absorption spectra plots
axes_mu_a[1, 0].plot(df_Shimojo_epi.loc[430:])
axes_mu_a[1, 0].set_xlabel("Wavelength (nm)")
axes_mu_a[1, 0].set_ylabel("Absorption (1/mm)")
axes_mu_a[1, 0].set_title(r"$\mu_a$ epidermis")

axes_mu_a[1, 1].plot(df_Shimojo_derm.loc[430:])
axes_mu_a[1, 1].set_xlabel("Wavelength (nm)")
axes_mu_a[1, 1].set_ylabel("Absorption (1/mm)")
axes_mu_a[1, 1].set_title(r"$\mu_a$ dermis")

axes_mu_a[1, 2].plot(df_Shimojo_sub.loc[430:])
axes_mu_a[1, 2].set_xlabel("Wavelength (nm)")
axes_mu_a[1, 2].set_ylabel("Absorption (1/mm)")
axes_mu_a[1, 2].set_title(r"$\mu_a$ subcutaneous")


fig_mu_s, axes_mu_s = plt.subplots(1, 3, layout = "constrained")
fig_mu_s.get_layout_engine().set(h_pad = 0.1, w_pad = 0.1)
fig_mu_s.suptitle("Scattering Spectras for Layers in Three-Layer Tissue Model, 429 - 1100nm")

axes_mu_s[0].plot(norm_wavelength_nm[429:1100], epidermis_mu_s[429:1100])
axes_mu_s[0].set_xlabel("Wavelength (nm)")
axes_mu_s[0].set_ylabel("Scattering (1/mm)")
axes_mu_s[0].set_title(r"$\mu_s$ epidermis")

axes_mu_s[1].plot(norm_wavelength_nm[429:1100], dermis_mu_s[429:1100])
axes_mu_s[1].set_xlabel("Wavelength (nm)")
axes_mu_s[1].set_ylabel("Scattering (1/mm)")
axes_mu_s[1].set_title(r"$\mu_s$ dermis")

axes_mu_s[2].plot(norm_wavelength_nm[429:1100], subcutaneous_mu_s[429:1100])
axes_mu_s[2].set_xlabel("Wavelength (nm)")
axes_mu_s[2].set_ylabel("Scattering (1/mm)")
axes_mu_s[2].set_title(r"$\mu_s$ subcutaneous tissue")

# Comment out to plot tissue spectras
# plt.close("all")


sim_wavelength_nm = 660

# Format: [mu_a, mu_s, g, n]
air_optic_params= [0.01,0.01, 1, 1]
epidermis_optic_params = [epidermis_mu_a[sim_wavelength_nm], epidermis_mu_s[sim_wavelength_nm], epidermis_g, epidermis_n]
dermis_optic_params = [dermis_mu_a[sim_wavelength_nm], dermis_mu_s[sim_wavelength_nm], dermis_g, dermis_n]
subcutaneous_optic_params = [subcutaneous_mu_a[sim_wavelength_nm], subcutaneous_mu_s[sim_wavelength_nm], subcutaneous_g, subcutaneous_n]


# Define simulation values
prop=[
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