import matplotlib.pyplot as plt
import numpy as np
import pandas as pd



# Water absorption spectra from 667-2500 nm, L. Kou 1993
df_L_Kou_1993 = pd.read_csv("Water_Absorption_Data/L_Kou_water_absorption_1993.csv", comment = "#")
df_L_Kou_1993_wav_sorted = df_L_Kou_1993.sort_values("lambda (nm)")
L_Kou_1993_wavelength_nm = df_L_Kou_1993_wav_sorted["lambda (nm)"]
L_Kou_1993_mu_a_water_cm = df_L_Kou_1993_wav_sorted["absorption (1/cm)"]
L_Kou_1993_mu_a_water_mm = L_Kou_1993_mu_a_water_cm / 10

#Water absorption from 171-2000nm, Fewell 2019
df_Fewell_2019 = pd.read_csv("Water_Absorption_Data/Fewell_water_absorption_modified_2019.csv", comment = "#")
Fewell_2019_wavelength_nm = df_Fewell_2019["wavelength (nm)"]
Fewell_2019_mu_a_water_m = df_Fewell_2019["pure water value"]
Fewell_2019_mu_a_water_mm = Fewell_2019_mu_a_water_m / 1000

# Hemoglobin and Deoxyginated Hemoglobin absorption spectra from 250-2000nm, Bosschart 2014
# Blank entries in data were removed beforehand
df_Bosschaart_2014 = pd.read_csv("Blood_Absorption_Data/Bosschaart_blood_absorption_2014.csv", comment = "#")
Bosschaart_2014_wavelength_nm = df_Bosschaart_2014["wavelength_nm"]
Bosschaart_2014_mu_a_oxy_mm  = df_Bosschaart_2014["mu_a_oxy_mm-1"]
Bosschaart_2014_mu_a_deoxy_mm = df_Bosschaart_2014["mu_a_deoxy_mm-1"]

# Fat absorption spectra from Veen 2004
# Notes: 
# - Absorption coefficient is in (1/m)
# - Greater measurement uncertainty between 429-750 nm
df_Veen_2004 = pd.read_csv("Fat_Absorption_Data/Veen_2004_Fat_Absorption_Spectra.csv", comment = "#")
Veen_2004_wavelength_nm = df_Veen_2004["nm"]
Veen_2004_mu_a_fat_m = df_Veen_2004["mu_a"]
Veen_2004_mu_a_fat_mm = Veen_2004_mu_a_fat_m / 1000

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

    axes[0, 2].plot(Veen_2004_wavelength_nm, Veen_2004_mu_a_fat_mm)
    axes[0, 2].set_title(r"$\mu_{a, fat}$ for 429-1100 nm""\nVeen 2004")
    axes[0, 2].set_xlabel("Wavelength (nm)")
    axes[0, 2].set_ylabel("Absorption (1/mm)")

    axes[1, 2].plot(Jacques_1998_wavelength_nm, Jacques_1998_mu_a_melanosome_mm)
    axes[1, 2].set_title(r"$\mu_{a, melanosome}$ for 200 - 2000 nm""\nJacques 1998")
    axes[1, 2].set_xlabel("Wavelength (nm)")
    axes[1, 2].set_ylabel("Absorption (1/mm)")




norm_wavelength_nm = np.arange(0, 2001, 1)

# interpolated tissue consituent absorption spectras
# the below spectras will now be interpolated within 200 - 2000 nm in 1nm increments
L_Kou_1993_interp_np = np.interp(norm_wavelength_nm, L_Kou_1993_wavelength_nm, L_Kou_1993_mu_a_water_mm, left = 0, right = 0)
Fewell_2019_interp_np = np.interp(norm_wavelength_nm, Fewell_2019_wavelength_nm, Fewell_2019_mu_a_water_mm, left = 0, right = 0)
Bosschaart_2014_interp_oxy_np = np.interp(norm_wavelength_nm, Bosschaart_2014_wavelength_nm, Bosschaart_2014_mu_a_oxy_mm, left = 0, right = 0)
Bosschaart_2014_interp_deoxy_np = np.interp(norm_wavelength_nm, Bosschaart_2014_wavelength_nm, Bosschaart_2014_mu_a_deoxy_mm, left = 0, right = 0)
Veen_2004_interp_np = np.interp(norm_wavelength_nm, Veen_2004_wavelength_nm, Veen_2004_mu_a_fat_mm, left = 0, right = 0)
Jacques_1998_interp_np = np.interp(norm_wavelength_nm, Jacques_1998_wavelength_nm, Jacques_1998_mu_a_melanosome_mm, left = 0, right = 0)

# Converting interpolated absorption values from NumPy to Series objects since we will need to keep track of the wavelength index for each absoprtion value
L_Kou_1993_interp_ser = pd.Series(data = L_Kou_1993_interp_np, index = norm_wavelength_nm)
Fewell_2019_interp_ser = pd.Series(data = Fewell_2019_interp_np, index = norm_wavelength_nm)
Bosschaart_2014_interp_oxy_ser = pd.Series(data = Bosschaart_2014_interp_oxy_np, index = norm_wavelength_nm)
Bosschaart_2014_interp_deoxy_ser = pd.Series(data = Bosschaart_2014_interp_deoxy_np, index = norm_wavelength_nm)
Veen_2004_interp_ser = pd.Series(data = Veen_2004_interp_np, index = norm_wavelength_nm)
Jacques_1998_interp_ser = pd.Series(data = Jacques_1998_interp_np, index = norm_wavelength_nm)


# Plotting the interpolated Series data from 600 - 1100nm in 1nm increments
lower_wavelenghth_bound = 600
upper_wavelength_bound = 1100

# Using Jacques 2013 generic tissue model
# Per layer: (B * S * mu_a_HbO2) + (B * (1-S) * mu_a_Hb) + (W * mu_a_water) + (F * mu_a_fat) + (M * mu_a_melanosome) + baseline
# where B, S, W, F, and M represent blood, oxygen, water, fat, and melanosome compostion of tissue respectively

epidermis_Jacques = {
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0.20,
    "fat" : 0,
    "melanosome" : 0.025,
}

dermis_Jacques = {
    "blood" : 0.0020,
    "oxygen_saturation" : 0.39,
    "water" : 0.65,
    "fat" : 0,
    "melanosome" : 0,
}

subcutaneous_Jacques = {
    "blood" : 0.0054,
    "oxygen_saturation" : 0.76,
    "water" : 0.11,
    "fat" : 0.69,
    "melanosome" : 0,
}

subcutaneous_Ganesan = {
    "blood" : 0.00634,
    "oxygen_saturation" : 0.60,
    "water" : 0.2184,
    "fat" : 0.6433,
    "melanosome" : 0,
}

# Basline absorption effect
baseline_epi = (8.53 * np.exp(-(norm_wavelength_nm - 154)/66.2) + 0.0244)
baseline_epi_frac = (1 
                     - epidermis_Jacques["blood"] 
                     - epidermis_Jacques["fat"] 
                     - epidermis_Jacques["water"] 
                     - epidermis_Jacques["melanosome"]
                     ) 

baseline_derm = (1.682 * np.exp(-(norm_wavelength_nm - 400)/80.5) + 0.0244)
baseline_derm_frac = (1 
                     - dermis_Jacques["blood"] 
                     - dermis_Jacques["fat"] 
                     - dermis_Jacques["water"] 
                     - dermis_Jacques["melanosome"]
                     ) 

baseline_sub = (1.682 * np.exp(-(norm_wavelength_nm - 400)/80.5) + 0.0244)
baseline_sub_frac = (1 
                     - subcutaneous_Jacques["blood"] 
                     - subcutaneous_Jacques["fat"] 
                     - subcutaneous_Jacques["water"] 
                     - subcutaneous_Jacques["melanosome"]
                     ) 

# =========================================
# Weighted tissue constituent absorption spectra for per layer in three layer tissue model
# =========================================

# Epidermis
mu_a_water_epi = (
                 epidermis_Jacques["water"] 
                 * L_Kou_1993_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                 )
mu_a_oxy_Hb_epi = ( 
                  epidermis_Jacques["blood"] 
                  * epidermis_Jacques["oxygen_saturation"] 
                  * Bosschaart_2014_interp_oxy_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                  )
mu_a_deoxy_epi = (
                 epidermis_Jacques["blood"]
                 * (1 - epidermis_Jacques["oxygen_saturation"]) 
                 * Bosschaart_2014_interp_deoxy_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                 )
mu_a_fat_epi = (
               epidermis_Jacques["fat"] 
               * Veen_2004_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
               )
mu_a_melanosome_epi = (
                      epidermis_Jacques["melanosome"] 
                      * Jacques_1998_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                      )
mu_a_baseline_epi = (
                    baseline_epi_frac
                    * baseline_epi[lower_wavelenghth_bound : upper_wavelength_bound ]
                    )
# Dermis
mu_a_water_derm = (
                 dermis_Jacques["water"] 
                 * L_Kou_1993_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                 )
mu_a_oxy_Hb_derm = ( 
                  dermis_Jacques["blood"] 
                  * dermis_Jacques["oxygen_saturation"] 
                  * Bosschaart_2014_interp_oxy_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                  )
mu_a_deoxy_derm = (
                 dermis_Jacques["blood"]
                 * (1 - dermis_Jacques["oxygen_saturation"]) 
                 * Bosschaart_2014_interp_deoxy_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                 )
mu_a_fat_derm = (
               dermis_Jacques["fat"] 
               * Veen_2004_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
               )
mu_a_melanosome_derm = (
                      dermis_Jacques["melanosome"] 
                      * Jacques_1998_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                      )
mu_a_baseline_derm = (
                    baseline_derm_frac
                    * baseline_derm[lower_wavelenghth_bound : upper_wavelength_bound ]
                    )

# Subcutaneous
mu_a_water_sub = (
                 subcutaneous_Jacques["water"] 
                 * L_Kou_1993_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                 )
mu_a_oxy_Hb_sub = ( 
                  subcutaneous_Jacques["blood"] 
                  * subcutaneous_Jacques["oxygen_saturation"] 
                  * Bosschaart_2014_interp_oxy_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                  )
mu_a_deoxy_sub = (
                 subcutaneous_Jacques["blood"]
                 * (1 - subcutaneous_Jacques["oxygen_saturation"]) 
                 * Bosschaart_2014_interp_deoxy_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                 )
mu_a_fat_sub = (
               subcutaneous_Jacques["fat"] 
               * Veen_2004_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
               )
mu_a_melanosome_sub = (
                      subcutaneous_Jacques["melanosome"] 
                      * Jacques_1998_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                      )
mu_a_baseline_sub = (
                    baseline_sub_frac
                    * baseline_sub[lower_wavelenghth_bound : upper_wavelength_bound ]
                    )


fig_i, axes_i = plt.subplots(1,3, layout = "constrained")
fig_i.get_layout_engine().set(h_pad = 0.1, w_pad = 0.1)
fig_i.suptitle("Weighted Tissue Constituent Absorption Spectra from 600-1100 nm, Jacques 2013", fontsize = 18)

for ax in axes_i:
    ax.set_xlabel("Wavelength (nm)", fontsize = 18)
    ax.set_ylabel("Absorption (1/mm)", fontsize = 18)
    ax.tick_params(axis = "both", labelsize = 18)

# Epidermis
axes_i[0].set_title("Epidermis", fontsize = 18)
axes_i[0].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_water_epi, 
                    label = r"$\mu_{a,water}$" "\nL. Kou. 1993",
                    )
axes_i[0].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_oxy_Hb_epi, 
                    label = r"$\mu_{a, Hb_{O_2}}$""\nBosschaart, 2014",
                    )
axes_i[0].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_deoxy_epi, 
                    label = r"$\mu_{a, Hb}$""\nBosschaart, 2014",
                    )
axes_i[0].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_fat_epi, 
                    label = r"$\mu_{a, fat}$""\nVeen 2004",
                    )
axes_i[0].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_melanosome_epi, 
                    label = r"$\mu_{a, melanosome}$""\nJacques 1998",
                    )
axes_i[0].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_baseline_epi, 
                    label = "baseline""\nMoco 2018",
                    color = "pink",
                    )
# Dermis
axes_i[1].set_title("Dermis", fontsize = 18)
axes_i[1].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_water_derm, 
                    )
axes_i[1].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_oxy_Hb_derm, 
                    )
axes_i[1].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_deoxy_derm,
                    )
axes_i[1].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_fat_derm,
                    )
axes_i[1].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_melanosome_derm,
                    )
axes_i[1].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_baseline_derm,
                    color = "pink",
                    )

# Subcutaneous
axes_i[2].set_title("Subcutaneous", fontsize = 18)
axes_i[2].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_water_sub,
                    )
axes_i[2].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_oxy_Hb_sub,
                    )
axes_i[2].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_deoxy_sub,
                    )
axes_i[2].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_fat_sub,
                    )
axes_i[2].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_melanosome_sub,
                    )
axes_i[2].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_baseline_sub,
                    color = "pink",
                    )

fig_i.legend()



# Using Jacques 2013 generic tissue model
# Per layer: (B * S * mu_a_HbO2) + (B * (1-S) * mu_a_Hb) + (W * mu_a_water) + (F * mu_a_fat) + (M * mu_a_melanosome) + baseline
# where B, S, W, F, and M represent blood, oxygen, water, fat, and melanosome compostion of tissue respectively

def mu_a_epi (
    B: float,
    S: float,
    W: float,
    F: float,
    M: float,

) -> pd.Series:
    total_mu_a = ( 
    (B * S * Bosschaart_2014_interp_oxy_ser)
    + (B * (1-S) * Bosschaart_2014_interp_deoxy_ser)
    + (W * L_Kou_1993_interp_ser) 
    + (F * Veen_2004_interp_ser) 
    + (M * Jacques_1998_interp_ser)
    + ((1 - (B + W + F + M)) * ((8.53 * np.exp(-(norm_wavelength_nm - 154)/66.2)) + 0.0244)) # baseline tissue absorption
    )
    return total_mu_a

def mu_a_derm (
    B: float,
    S: float,
    W: float,
    F: float,
    M: float,

) -> pd.Series:
    total_mu_a = ( 
    (B * S * Bosschaart_2014_interp_oxy_ser)
    + (B * (1-S) * Bosschaart_2014_interp_deoxy_ser)
    + (W * L_Kou_1993_interp_ser) 
    + (F * Veen_2004_interp_ser) 
    + (M * Jacques_1998_interp_ser)
    + ((1 - (B + W + F + M)) * ((1.682 * np.exp(-(norm_wavelength_nm - 400)/80.5)) + 0.0244)) # baseline tissue absorption
    )
    return total_mu_a

def mu_a_sub (
    B: float,
    S: float,
    W: float,
    F: float,
    M: float,

) -> pd.Series:
    total_mu_a = ( 
    (B * S * Bosschaart_2014_interp_oxy_ser)
    + (B * (1-S) * Bosschaart_2014_interp_deoxy_ser)
    + (W * L_Kou_1993_interp_ser) 
    + (F * Veen_2004_interp_ser) 
    + (M * Jacques_1998_interp_ser)
    + ((1 - (B + W + F + M)) * ((1.682 * np.exp(-(norm_wavelength_nm - 400)/80.5)) + 0.0244)) # baseline tissue absorption
    )
    return total_mu_a

def mu_s_layer(
    lay: str, 
    g: float,
) -> pd.Series:
    match(lay):
        case("epidermis"):
            mu_s_rd = 6.87 * (norm_wavelength_nm[1:] / 500) ** -1.161
            mu_s = mu_s_rd/(1-g)
        case("dermis"):
            mu_s_rd = 4.53 * (norm_wavelength_nm[1:] / 500) ** -1.292
            mu_s = mu_s_rd/(1-g)
        case("subcutaneous"):
            mu_s_rd = 1.37 * (norm_wavelength_nm[1:] / 500) ** -0.385
            mu_s = mu_s_rd/(1-g)
        case _:
            mu_s = 0
            print("invalid layer")
    mu_s_ser = pd.Series(data = mu_s, index = norm_wavelength_nm[1:])
    return mu_s_ser
