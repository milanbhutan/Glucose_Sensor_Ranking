import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

derm_gamma = 0.5 # was 0.50
sub_gamma = 0.25

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
lower_wavelenghth_bound = 600 # original: 600
upper_wavelength_bound = 1100

# Using Jacques 2013 generic tissue model
# Per layer: (B * S * mu_a_HbO2) + (B * (1-S) * mu_a_Hb) + (W * mu_a_water) + (F * mu_a_fat) + (M * mu_a_melanosome) + baseline
# where B, S, W, F, and M represent blood, oxygen, water, fat, and melanosome compostion of tissue respectively


epi_Jacques = {
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0,
    "fat" : 0,
    "melanosome" : 0.025,
}

epi_Jolivot_low_mel ={
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0,
    "fat" : 0,
    "melanosome": 0.013
}

epi_Jolivot_SPT_II ={
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0,
    "fat" : 0,
    "melanosome": 0.101
}

epi_Jolivot_SPT_III ={
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0,
    "fat" : 0,
    "melanosome": 0.187
}

epi_Jolivot_SPT_IV ={
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0,
    "fat" : 0,
    "melanosome": 0.265
}

epi_Jolivot_SPT_V ={
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0,
    "fat" : 0,
    "melanosome": 0.358
}

epi_Jolivot_SPT_VI = {
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0,
    "fat" : 0,
    "melanosome": 0.447
}

derm_Jacques = {
    "blood" : 0.0020,  
    "oxygen_saturation" : 0.39,
    "water" : 0.65,
    "fat" : 0,
    "melanosome" : 0,
    "gamma" : derm_gamma
}

lower_derm_Jonasson = {
    "blood" : 0.0075,  
    "oxygen_saturation" : 0.50, 
    "water" : 0.0, # Jonasson neglects water
    "fat" : 0,
    "melanosome" : 0,
    "gamma" : derm_gamma,
}

sub_Jacques = {
    "blood" : 0.0054,
    "oxygen_saturation" : 0.76,
    "water" : 0.11,
    "fat" : 0.69,
    "melanosome" : 0,
    "gamma" : derm_gamma,
}

sub_Ganesan = {
    "blood" : 0.00634,
    "oxygen_saturation" : 0.60,
    "water" : 0.2184,
    "fat" : 0.6433,
    "melanosome" : 0,
    "gamma" : sub_gamma,
}



# Basline absorption effect
baseline_epi = (8.53 * np.exp(-(norm_wavelength_nm - 154)/66.2) + 0.0244)
baseline_epi_frac = (1 
                     - epi_Jacques["blood"] 
                     - epi_Jacques["fat"] 
                     - epi_Jacques["water"] 
                     - epi_Jacques["melanosome"]
                     ) 

baseline_derm = 0.5 * (derm_Jacques["water"] / 0.65 ) * (1.682 * np.exp(-(norm_wavelength_nm - 400)/80.5) + 0.0244)
baseline_derm_frac = (1 
                     - derm_Jacques["blood"] 
                     - derm_Jacques["fat"] 
                     - derm_Jacques["water"] 
                     - derm_Jacques["melanosome"]
                     ) 

baseline_sub = 0.25 * ((sub_Jacques["water"] / 0.65 )) * (1.682 * np.exp(-(norm_wavelength_nm - 400)/80.5) + 0.0244)
baseline_sub_frac = (1 
                     - sub_Jacques["blood"] 
                     - sub_Jacques["fat"] 
                     - sub_Jacques["water"] 
                     - sub_Jacques["melanosome"]
                     ) 

# =========================================
# Weighted tissue constituent absorption spectra for per layer in three layer tissue model
# =========================================

# Epidermis
mu_a_water_epi = (
                 epi_Jacques["water"] 
                 * Fewell_2019_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                 )
mu_a_oxy_Hb_epi = ( 
                  epi_Jacques["blood"] 
                  * epi_Jacques["oxygen_saturation"] 
                  * Bosschaart_2014_interp_oxy_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                  )
mu_a_deoxy_epi = (
                 epi_Jacques["blood"]
                 * (1 - epi_Jacques["oxygen_saturation"]) 
                 * Bosschaart_2014_interp_deoxy_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                 )
mu_a_fat_epi = (
               epi_Jacques["fat"] 
               * Veen_2004_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
               )
mu_a_melanosome_epi = (
                      epi_Jacques["melanosome"] 
                      * Jacques_1998_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                      )
mu_a_baseline_epi = (
                    baseline_epi_frac
                    * baseline_epi[lower_wavelenghth_bound : upper_wavelength_bound ]
                    )
# Dermis
mu_a_water_derm = (
                 derm_Jacques["water"] 
                 * Fewell_2019_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                 )
mu_a_oxy_Hb_derm = ( 
                  derm_Jacques["blood"] 
                  * derm_Jacques["oxygen_saturation"] 
                  * Bosschaart_2014_interp_oxy_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                  )
mu_a_deoxy_derm = (
                 derm_Jacques["blood"]
                 * (1 - derm_Jacques["oxygen_saturation"]) 
                 * Bosschaart_2014_interp_deoxy_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                 )
mu_a_fat_derm = (
               derm_Jacques["fat"] 
               * Veen_2004_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
               )
mu_a_melanosome_derm = (
                      derm_Jacques["melanosome"] 
                      * Jacques_1998_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                      )
mu_a_baseline_derm = (
                    baseline_derm_frac
                    * baseline_derm[lower_wavelenghth_bound : upper_wavelength_bound ]
                    )

# Subcutaneous
mu_a_water_sub = (
                 sub_Jacques["water"] 
                 * Fewell_2019_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                 )
mu_a_oxy_Hb_sub = ( 
                  sub_Jacques["blood"] 
                  * sub_Jacques["oxygen_saturation"] 
                  * Bosschaart_2014_interp_oxy_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                  )
mu_a_deoxy_sub = (
                 sub_Jacques["blood"]
                 * (1 - sub_Jacques["oxygen_saturation"]) 
                 * Bosschaart_2014_interp_deoxy_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                 )
mu_a_fat_sub = (
               sub_Jacques["fat"] 
               * Veen_2004_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
               )
mu_a_melanosome_sub = (
                      sub_Jacques["melanosome"] 
                      * Jacques_1998_interp_ser[lower_wavelenghth_bound : upper_wavelength_bound ]
                      )
mu_a_baseline_sub = (
                    baseline_sub_frac
                    * baseline_sub[lower_wavelenghth_bound : upper_wavelength_bound ]
                    )


fig_i, axes_i = plt.subplots(1,3, layout = "constrained")
fig_i.get_layout_engine().set(h_pad = 0.1, w_pad = 0.1)
fig_i.suptitle("Weighted Tissue Constituent Absorption Spectra from 600-1100 nm", fontsize = 18)

for ax in axes_i:
    ax.set_xlabel("Wavelength (nm)", fontsize = 18)
    ax.set_ylabel("Absorption (1/mm)", fontsize = 18)
    ax.tick_params(axis = "both", labelsize = 18)

# Epidermis
axes_i[0].set_title("Epidermis", fontsize = 18)
axes_i[0].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_water_epi, 
                    label = r"$\mu_{a,water}$" "\nFewell 2019",
                    )
axes_i[0].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_oxy_Hb_epi, 
                    label = r"$\mu_{a, Hb_{O_2}}$""\nBosschaart 2014",
                    )
axes_i[0].plot(norm_wavelength_nm[lower_wavelenghth_bound : upper_wavelength_bound ], 
                    mu_a_deoxy_epi, 
                    label = r"$\mu_{a, Hb}$""\nBosschaart 2014",
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


# ==========================================================
# Using Jacques 2013 generic tissue model
# Per layer: (B * S * mu_a_HbO2) + (B * (1-S) * mu_a_Hb) + (W * mu_a_water) + (F * mu_a_fat) + (M * mu_a_melanosome) + baseline
# where B, S, W, F, and M represent blood, oxygen, water, fat, and melanosome compostion of tissue respectively
# ==========================================================


def mu_a_epi (
    B: float,
    S: float,
    W: float,
    F: float,
    M: float,

) -> pd.Series:
    R = 1 - (B + W + F + M)
    if(R < 0):
        R = 0
    total_mu_a = ( 
    (B * S * Bosschaart_2014_interp_oxy_ser)
    + (B * (1-S) * Bosschaart_2014_interp_deoxy_ser)
    + (W * Fewell_2019_interp_ser) 
    + (F * Veen_2004_interp_ser) 
    + (M * Jacques_1998_interp_ser)
    + (R * ((8.53 * np.exp(-(norm_wavelength_nm - 154)/66.2)) + 0.0244)) # baseline tissue absorption
    )
    return total_mu_a

def mu_a_derm (
    B: float,
    S: float,
    W: float,
    F: float,
    M: float,
    gamma: float = derm_gamma

) -> pd.Series:
    R = 1 - (B + W + F + M)
    if(R < 0):
        R = 0
    # C_w = W / (1 - (B + F + M) ) # nested water fraction used by Moco 2018 for baseline
    total_mu_a = ( 
    (B * S * Bosschaart_2014_interp_oxy_ser)
    + (B * (1-S) * Bosschaart_2014_interp_deoxy_ser)
    + (W * Fewell_2019_interp_ser) 
    + (F * Veen_2004_interp_ser) 
    + (M * Jacques_1998_interp_ser)
    + (R * gamma  * ((1.682 * np.exp(-(norm_wavelength_nm - 400)/80.5)) + 0.0244)) # baseline tissue absorption
    )
    return total_mu_a

def mu_a_sub (
    B: float,
    S: float,
    W: float,
    F: float,
    M: float,
    gamma: float = sub_gamma,

) -> pd.Series:
    R = 1 - (B + W + F + M)
    if(R < 0):
        raise ValueError("Invalid Baseline Fraction")
    # C_w = W / (1 - (B + F + M) ) # nested water fraction used by Moco 2018 for baseline
    total_mu_a = ( 
    (B * S * Bosschaart_2014_interp_oxy_ser)
    + (B * (1-S) * Bosschaart_2014_interp_deoxy_ser)
    + (W * Fewell_2019_interp_ser) 
    + (F * Veen_2004_interp_ser) 
    + (M * Jacques_1998_interp_ser)
    + (R * gamma * ((1.682 * np.exp(-(norm_wavelength_nm - 400)/80.5)) + 0.0244)) # baseline tissue absorption
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
            mu_s_rd = 1.84 * (norm_wavelength_nm[1:] / 500) ** -0.672
            mu_s = mu_s_rd/(1-g)
        case _:
            mu_s = 0
            print("invalid layer")
    mu_s_ser = pd.Series(data = mu_s, index = norm_wavelength_nm[1:])
    return mu_s_ser