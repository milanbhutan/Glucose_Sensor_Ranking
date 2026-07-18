import pmcx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#aborption spectra extracted from Prahl


#Plot configuration
fig, axes = plt.subplots(2, 1, layout = "constrained")
fig.get_layout_engine().set(h_pad = 0.1, w_pad = 0.1)


# #Water absorption spectra from 382.5-750 nm, Pope 1993
# df_Pope_1993 = pd.read_csv("R_M_Pope_water_absorption_1993.csv", comment = "#")
# Pope_1993_wavelength_nm = df_Pope_1993["lambda (nm)"]
# Pope_1993_mu_a_water_cm = df_Pope_1993["absorption (1/cm)"]
# Pope_1993_mu_a_water_mm = Pope_1993_mu_a_water_cm / 10

# axes[0,0].plot(Pope_1993_wavelength_nm, Pope_1993_mu_a_water_mm)
# axes[0,0].set_xlabel("Wavelength (nm)")
# axes[0,0].set_xticks(np.arange(400,750,25), minor = True)
# axes[0,0].set_ylabel("Absorption (1/mm)")
# axes[0,0].set_title(r"$\mu_{a,water}$ from 400-750 nm" "\n R. M. Pope, 1993")


#Water absorption spectra from 667-2500 nm, Lou 1993
df_Lou_1993 = pd.read_csv("Water_Absorption_Data/L_Kou_water_absorption_1993.csv", comment = "#")
Lou_1993_wavelength_nm = df_Lou_1993["lambda (nm)"]
Lou_1993_mu_a_water_cm = df_Lou_1993["absorption (1/cm)"]
Lou_1993_mu_a_water_mm = Lou_1993_mu_a_water_cm / 10

axes[0].plot(Lou_1993_wavelength_nm, Lou_1993_mu_a_water_mm)
axes[0].set_xlabel("Wavelength (nm)")
axes[0].set_xticks(np.arange(750, 2500, 125), minor = True)
axes[0].set_ylabel("Absorption (1/mm)")
axes[0].set_title(r"$\mu_{a,water}$ from 667-2500nm" "\n L. Kou, 1993")


# #Water absorption spectra from 380-700 nm, Pope 1997
# df_Pope_1997 = pd.read_csv("R_M_Pope_water_absorption_1997.csv", comment = "#")
# Pope_1997_wavelength_nm = df_Pope_1997["lambda (nm)"]
# Pope_1997_mu_a_water_cm = df_Pope_1997["absorption (1/cm)"]
# Pope_1997_mu_a_water_mm = Pope_1997_mu_a_water_cm / 10;

# axes[1,0].plot(Pope_1997_wavelength_nm, Pope_1997_mu_a_water_mm);
# axes[1,0].set_xlabel("Wavelength (nm)")
# axes[1,0].set_xticks(np.arange(400,750,25), minor = True)
# axes[1,0].set_ylabel("Absorption (1/mm)")
# axes[1,0].set_title(r"$\mu_{a,water}$ from 380-700nm" "\n R. M. Pope, 1997");


#Water absorption from 171-2000nm, Fewell 2019
df_Fewell_2019 = pd.read_csv("Water_Absorption_Data/Fewell_water_absorption_modified_2019.csv", comment = "#")
Fewell_2019_wavelength_nm = df_Fewell_2019["wavelength (nm)"]
Fewell_2019_mu_a_water_m = df_Fewell_2019["pure water value"]
Fewell_2019_mu_a_water_mm = Fewell_2019_mu_a_water_m / 1000

axes[1].plot(Fewell_2019_wavelength_nm, Fewell_2019_mu_a_water_mm)
axes[1].set_xlabel("Wavelength (nm)")
axes[1].set_xticks(np.arange(500, 2000, 100), minor = True)
axes[1].set_ylabel("Absorption (1/mm)")
axes[1].set_title(r"$\mu_{a,water}$ from 400 to 2000 nm" "\nM. P. Fewell, 2019");




#Hemoglobin

#Deoxygninated Hemoglobin
















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
vol[:,:,0:3]=1;
vol[:,:,3:20]=2;
vol[:,:,20:60]=3;


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
#flux=sim_out["flux"];

#flo_plot=imshow(np.log10())