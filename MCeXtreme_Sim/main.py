import pmcx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import layer




# From Jacques 2013
epidermis = {
    "blood" : 0,
    "oxygen" : 0,
    "water" : 0.20,
    "fat" : 0,
    "melanin" : 0.025
}

dermis = {
    "blood" : 0.0020,
    "oxygen" : 0.39,
    "water" : 0.65,
    "fat" : 0,
    "melanin" : 0
}

sub_tissue = {
    "blood" : 0.0054,
    "oxygen" : 0.76,
    "water" : 0.11,
    "fat" : 0.69,
    "melanin" : 0
}

norm_wavelength_nm = np.arange(200, 2001, 1)

epidermis_mu_a = layer.mu_a_layer(epidermis["blood"], epidermis["oxygen"], epidermis["water"], epidermis["fat"], epidermis["melanin"])
dermis_mu_a = layer.mu_a_layer(dermis["blood"], dermis["oxygen"], dermis["water"], dermis["fat"], dermis["melanin"])
sub_tissue_mu_a = layer.mu_a_layer(sub_tissue["blood"], sub_tissue["oxygen"], sub_tissue["water"], sub_tissue["fat"], sub_tissue["melanin"])


# Note: This spectra is only valid within 429-1100 nm where reliable tissue constituent abosroption data is availible
fig_lay, axes_lay = plt.subplots(1, 3, layout = "constrained")
fig_lay.get_layout_engine().set(h_pad = 0.1, w_pad = 0.1)
fig_lay.suptitle("Absorption Spectras for Layers in Three-Layer Tissue Model")

print(norm_wavelength_nm[229:900])
print(epidermis_mu_a[229:900])

axes_lay[0].plot(norm_wavelength_nm[229:900], epidermis_mu_a[229:900])
axes_lay[0].set_xlabel("Wavelength (nm)")
axes_lay[0].set_ylabel("Absorption (1/mm)")
axes_lay[0].set_title(r"$\mu_a$ epidermis")

axes_lay[1].plot(norm_wavelength_nm[229:900], dermis_mu_a[229:900])
axes_lay[1].set_xlabel("Wavelength (nm)")
axes_lay[1].set_ylabel("Absorption (1/mm)")
axes_lay[1].set_title(r"$\mu_a$ dermis")

axes_lay[2].plot(norm_wavelength_nm[229:900], sub_tissue_mu_a[229:900])
axes_lay[2].set_xlabel("Wavelength (nm)")
axes_lay[2].set_ylabel("Absorption (1/mm)")
axes_lay[2].set_title(r"$\mu_a$ subcutaneous tissue")

plt.show()

























# Format: [mu_a, mu_s, g, n]
air_optic_params= [0.01,0.01,0.01,0.01]
epidermis_optic_params = [0.01,0.01,0.01,0.01]
dermis_optic_params = [0.01,0.01,0.01,0.01]
sub_tissue_optic_params = [0.01,0.01,0.01,0.01]

#absoprtion spectra


#Define simulation values
prop=[
    air_optic_params,
    epidermis_optic_params,
    dermis_optic_params,
    sub_tissue_optic_params
]


Nx, Ny, Nz= 60, 60, 60

vol=np.zeros([Nx,Ny,Nz])
vol[:,:,0:3]=1
vol[:,:,3:20]=2
vol[:,:,20:60]=3


sim_config={
    "nphoton": 10000,
    "vol": vol,
    "prop": prop,
    "srcpos": [30, 30, 0],
    "srcdir": [0, 0, 1],
    "detdir": [0, 0, 1, 0.01],
    "detpos": [[32, 30, 0, 1]],
    "tstart": 0,
    "tend": 1,
    "tstep": 0.0001,
    "unitinmm": 0.1,
    "autopilot": 1,
    "gpuid": 1
}

#Run simulation
#sim_out=pmcx.run(sim_config)

#print(sim_out.keys())
#flux=sim_out["flux"]

#flo_plot=imshow(np.log10())