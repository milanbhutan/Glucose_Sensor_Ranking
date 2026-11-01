import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from paths import MCX_DIR

wavelength_bound_1 = 600
wavelength_bound_2 = 999
wavelength_bound_3 = 1001
wavelength_bound_4 = 1700

norm_wavelength_nm_np = np.arange(1, 2001, 1)
norm_wavelength_nm = pd.Series(data = norm_wavelength_nm_np, index = norm_wavelength_nm_np)

# Absorption Data
df_Iino_epi_mu_a = (
    pd.read_csv(MCX_DIR / "Exp_Absorption" / "iino_epidermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
Iino_epidermis_mu_a = df_Iino_epi_mu_a.loc[:, "epidermis_mu_a_mm^-1"]


df_Iino_derm_mu_a = (
    pd.read_csv(MCX_DIR / "Exp_Absorption" / "iino_dermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
Iino_dermis_mu_a = df_Iino_derm_mu_a.loc[:, "dermis_mu_a_mm^-1"]

df_Iino_sub_mu_a = (
    pd.read_csv(MCX_DIR / "Exp_Absorption" / "iino_subcutaneous_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
Iino_subcutaneous_mu_a = df_Iino_sub_mu_a.loc[:, "subcutaneous_mu_a_mm^-1"]

# Scattering Data
df_Iino_epi_red_mu_s = (
    pd.read_csv(MCX_DIR / "Exp_Scattering" / "iino_epidermis_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
Iino_epidermis_mu_red_s = df_Iino_epi_red_mu_s.loc[:, "epidermis_mu_s_prime_mm^-1"]

df_Iino_derm_red_mu_s = (
    pd.read_csv(MCX_DIR / "Exp_Scattering" / "iino_dermis_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
Iino_dermis_mu_red_s = df_Iino_derm_red_mu_s.loc[:, "dermis_mu_s_prime_mm^-1"]

df_Iino_sub_red_mu_s = (
    pd.read_csv(MCX_DIR / "Exp_Scattering" / "iino_subcutaneous_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
Iino_subcutaneous_mu_red_s = df_Iino_sub_red_mu_s.loc[:, "subcutaneous_mu_s_prime_mm^-1"]


# Water absorption spectra from 667-2500 nm, L. Kou 1993
df_L_Kou_1993 = pd.read_csv(MCX_DIR / "Water_Absorption_Data" / "L_Kou_water_absorption_1993.csv", comment = "#")
df_L_Kou_1993_wav_sorted = df_L_Kou_1993.sort_values("lambda (nm)")
L_Kou_1993_wavelength_nm = df_L_Kou_1993_wav_sorted["lambda (nm)"]
L_Kou_1993_mu_a_water_cm = df_L_Kou_1993_wav_sorted["absorption (1/cm)"]
L_Kou_1993_mu_a_water_mm = L_Kou_1993_mu_a_water_cm / 10

#Water absorption from 171-2000nm, Fewell 2019
df_Fewell_2019 = pd.read_csv(MCX_DIR / "Water_Absorption_Data" / "Fewell_water_absorption_modified_2019.csv", comment = "#")
Fewell_2019_wavelength_nm = df_Fewell_2019["wavelength (nm)"]
Fewell_2019_mu_a_water_m = df_Fewell_2019["pure water value"]
Fewell_2019_mu_a_water_mm = Fewell_2019_mu_a_water_m / 1000

# Hemoglobin and Deoxyginated Hemoglobin absorption spectra from 250-2000nm, Bosschart 2014
# Blank entries in data were removed beforehand
df_Bosschaart_2014 = pd.read_csv(MCX_DIR / "Blood_Absorption_Data" / "Bosschaart_blood_absorption_2014.csv", comment = "#")
Bosschaart_2014_wavelength_nm = df_Bosschaart_2014["wavelength_nm"]
Bosschaart_2014_mu_a_oxy_mm  = df_Bosschaart_2014["mu_a_oxy_mm-1"]
Bosschaart_2014_mu_a_deoxy_mm = df_Bosschaart_2014["mu_a_deoxy_mm-1"]

# Fat absorption spectra from Veen 2004
# Notes: 
# - Absorption coefficient is in (1/m)
# - Greater measurement uncertainty between 429-750 nm
# df_Veen_2004 = pd.read_csv("Fat_Absorption_Data/Veen_2004_Fat_Absorption_Spectra.csv", comment = "#")
# Veen_2004_wavelength_nm = df_Veen_2004["nm"]
# Veen_2004_mu_a_fat_m = df_Veen_2004["mu_a"]
# Veen_2004_mu_a_fat_mm = Veen_2004_mu_a_fat_m / 1000

df_Veen_Anderson = pd.read_csv(MCX_DIR / "Fat_Absorption_Data" / "fat_429_1700nm.csv", comment = "#")
Veen_Anderson_wavelength_nm = df_Veen_Anderson["nm"]
Veen_Anderson_mu_a_fat_mm = df_Veen_Anderson["mu_a"]



#SHOULD BE A SERIES
# melanosome absorption spectra from 200 - 2000 from Jacques 1998
# Need to find lower and upper wavelength range
Jacques_1998_wavelength_nm = np.arange(200, 2001, 1)
Jacques_1998_mu_a_melanosome_cm = (1.7*10**12)*(Jacques_1998_wavelength_nm)**-3.48
Jacques_1998_mu_a_melanosome_mm = Jacques_1998_mu_a_melanosome_cm / 10

if __name__ == "__main__":

    # Plot configuration
    fig, axes = plt.subplots(2, 3, layout = "constrained")
    fig.get_layout_engine().set(h_pad = 0.1, w_pad = 0.1)
    fig.suptitle("Tissue Constituent Absorption Spectra")

    axes[0, 0].plot(L_Kou_1993_wavelength_nm, L_Kou_1993_mu_a_water_mm)
    axes[0, 0].set_xlabel("Wavelength (nm)")
    axes[0, 0].set_ylabel("Absorption (1/mm)")
    axes[0, 0].set_title(r"$\mu_{a,water}$ from 667-2500nm" "\n L. Kou 1993")

    axes[1, 0].plot(Fewell_2019_wavelength_nm, Fewell_2019_mu_a_water_mm)
    axes[1, 0].set_xlabel("Wavelength (nm)")
    axes[1, 0].set_ylabel("Absorption (1/mm)")
    axes[1, 0].set_title(r"$\mu_{a,water}$ from 400 to 2000 nm" "\nM. P. Fewell 2019")

    axes[0, 1].plot(Bosschaart_2014_wavelength_nm, Bosschaart_2014_mu_a_oxy_mm)
    axes[0, 1].set_title(r"$\mu_{a, Hb_{O_2}}$ from 250-2000 nm""\nBosschaart 2014")
    axes[0, 1].set_xlabel("Wavelength (nm)")
    axes[0, 1].set_ylabel("Absorption (1/mm)")

    axes[1, 1].plot(Bosschaart_2014_wavelength_nm, Bosschaart_2014_mu_a_deoxy_mm)
    axes[1, 1].set_title(r"$\mu_{a,Hb}$ from 250-2000 nm""\nBosschaart 2014")
    axes[1, 1].set_xlabel("Wavelength (nm)")
    axes[1, 1].set_ylabel("Absorption (1/mm)")

    axes[0, 2].plot(Veen_Anderson_wavelength_nm, Veen_Anderson_mu_a_fat_mm)
    axes[0, 2].set_title(r"$\mu_{a, fat}$ for Veen 2004, 429-1098 nm""\nAnderson 2006, 1099-2600 nm")
    axes[0, 2].set_xlabel("Wavelength (nm)")
    axes[0, 2].set_ylabel("Absorption (1/mm)")

    axes[1, 2].plot(Jacques_1998_wavelength_nm, Jacques_1998_mu_a_melanosome_mm)
    axes[1, 2].set_title(r"$\mu_{a, melanosome}$ for 200 - 2000 nm""\nJacques 1998")
    axes[1, 2].set_xlabel("Wavelength (nm)")
    axes[1, 2].set_ylabel("Absorption (1/mm)")
    plt.show()






# interpolated tissue consituent absorption spectras
# the below spectras will now be interpolated within 200 - 2000 nm in 1nm increments
L_Kou_1993_interp_np = np.interp(norm_wavelength_nm, L_Kou_1993_wavelength_nm, L_Kou_1993_mu_a_water_mm, left = 0, right = 0)
Fewell_2019_interp_np = np.interp(norm_wavelength_nm, Fewell_2019_wavelength_nm, Fewell_2019_mu_a_water_mm, left = 0, right = 0)
Bosschaart_2014_interp_oxy_np = np.interp(norm_wavelength_nm, Bosschaart_2014_wavelength_nm, Bosschaart_2014_mu_a_oxy_mm, left = 0, right = 0)
Bosschaart_2014_interp_deoxy_np = np.interp(norm_wavelength_nm, Bosschaart_2014_wavelength_nm, Bosschaart_2014_mu_a_deoxy_mm, left = 0, right = 0)
Veen_Anderson_interp_np = np.interp(norm_wavelength_nm, Veen_Anderson_wavelength_nm, Veen_Anderson_mu_a_fat_mm, left = 0, right = 0)
Jacques_1998_interp_np = np.interp(norm_wavelength_nm, Jacques_1998_wavelength_nm, Jacques_1998_mu_a_melanosome_mm, left = 0, right = 0)

# Converting interpolated absorption values from NumPy to Series objects since we will need to keep track of the wavelength index for each absoprtion value
L_Kou_1993_interp_ser = pd.Series(data = L_Kou_1993_interp_np, index = norm_wavelength_nm)
Fewell_2019_interp_ser = pd.Series(data = Fewell_2019_interp_np, index = norm_wavelength_nm)
Bosschaart_2014_interp_oxy_ser = pd.Series(data = Bosschaart_2014_interp_oxy_np, index = norm_wavelength_nm)
Bosschaart_2014_interp_deoxy_ser = pd.Series(data = Bosschaart_2014_interp_deoxy_np, index = norm_wavelength_nm)
Veen_Anderson_interp_ser = pd.Series(data = Veen_Anderson_interp_np, index = norm_wavelength_nm)
Jacques_1998_interp_ser = pd.Series(data = Jacques_1998_interp_np, index = norm_wavelength_nm)