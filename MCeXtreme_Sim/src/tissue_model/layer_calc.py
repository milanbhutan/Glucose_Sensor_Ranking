import numpy as np
import pandas as pd
import tissue_model.exp_data as exp



# ==========================================================
# Using Jacques 2013 generic tissue model
# Per layer: (B * S * mu_a_HbO2) + (B * (1-S) * mu_a_Hb) + (W * mu_a_water) + (F * mu_a_fat) + (M * mu_a_melanosome) + baseline
# where B, S, W, F, and M represent blood, oxygen, water, fat, and melanosome compostion of tissue respectively
# ==========================================================

# Jacques 2013
n_dry = 1.5
n_water = 1.33

# van Gemert 1989
norm_wavelength_nm = np.arange(0, 2001, 1)
g_200_1000 = 0.62 + (0.29 * norm_wavelength_nm[0: 1000] * 10 ** -3)
g_1000_1700 = 0.9

g_200_1000_ser = pd.Series(data = g_200_1000, index = norm_wavelength_nm[0: 1000])
g_1000_1700_ser = pd.Series(data = g_1000_1700, index = norm_wavelength_nm[1000: 1701])

g = pd.concat([g_200_1000_ser, g_1000_1700_ser])

def epi(
    B: float,
    S: float,
    W: float,
    F: float,
    M: float,
    baseline: pd.Series,
) -> tuple[pd.Series, np.ndarray, np.ndarray, float, np.ndarray]:
    R = 1 - (B + W + F + M) # need R calculated here for uncertainty analysis
    
    if(R < 0):
        R = 0
    mu_a = ( 
    (B * S * exp.Bosschaart_2014_interp_oxy_ser)
    + (B * (1-S) * exp.Bosschaart_2014_interp_deoxy_ser)
    + (W * exp.Fewell_2019_interp_ser) 
    + (F * exp.Veen_Anderson_interp_ser) 
    + (M * exp.Jacques_1998_interp_ser)
    + (R * baseline) # baseline tissue absorption
    )

    # epidermis scattering
    mu_red_s = 6.87 * (norm_wavelength_nm[0: 1701] / 500) ** -1.161
    mu_s = mu_red_s / (1 - g)

    mu_red_s = pd.Series(data = mu_red_s, index = norm_wavelength_nm[0: 1701])
    mu_s = pd.Series(data = mu_s, index = norm_wavelength_nm[0: 1701])

    # Jacques 2013 refractive index
    n = ((1 - W) * n_dry) + W * n_water
    
    return mu_a, mu_s, g, n, mu_red_s

def derm(
    B: float,
    S: float,
    W: float,
    F: float,
    M: float,
    baseline: pd.Series,
) -> tuple[pd.Series, pd.Series, pd.Series, float, pd.Series]:
    R = 1 - (B + W + F + M) # need R calculated here for uncertainty analysis
    

    if(R < 0):
        R = 0
    mu_a = ( 
    (B * S * exp.Bosschaart_2014_interp_oxy_ser)
    + (B * (1-S) * exp.Bosschaart_2014_interp_deoxy_ser)
    + (W * exp.Fewell_2019_interp_ser) 
    + (F * exp.Veen_Anderson_interp_ser) 
    + (M * exp.Jacques_1998_interp_ser)
    + (R * baseline) # baseline tissue absorption
    )

    # dermis scattering
    mu_red_s = 4.53 * (norm_wavelength_nm[0: 1701] / 500) ** -1.292
    mu_s = mu_red_s / (1 - g)

    mu_red_s = pd.Series(data = mu_red_s, index = norm_wavelength_nm[0: 1701])
    mu_s = pd.Series(data = mu_s, index = norm_wavelength_nm[0: 1701])

    # Jacques 2013 refractive index
    n = ((1 - W) * n_dry) + W * n_water
    
    return mu_a, mu_s, g, n, mu_red_s


def sub(
    B: float,
    S: float,
    W: float,
    F: float,
    M: float,
    baseline: pd.Series,
) -> tuple[pd.Series, pd.Series, pd.Series, float, pd.Series]:
    R = 1 - (B + W + F + M) # need R calculated here for uncertainty analysis
    
    
    if(R < 0):
        R = 0
    mu_a = ( 
    (B * S * exp.Bosschaart_2014_interp_oxy_ser)
    + (B * (1-S) * exp.Bosschaart_2014_interp_deoxy_ser)
    + (W * exp.Fewell_2019_interp_ser) 
    + (F * exp.Veen_Anderson_interp_ser) 
    + (M * exp.Jacques_1998_interp_ser)
    + (R * baseline) # baseline tissue absorption
    )


    
    # subcutaneous scattering
    mu_red_s = 1.84 * (norm_wavelength_nm[0: 1701] / 500) ** -0.672
    mu_s = mu_red_s / (1 - g)

    mu_red_s = pd.Series(data = mu_red_s, index = norm_wavelength_nm[0: 1701])
    mu_s = pd.Series(data = mu_s, index = norm_wavelength_nm[0: 1701])

    # Jacques 2013 refractive index
    n = ((1 - W) * n_dry) + W * n_water
    
    return mu_a, mu_s, g, n, mu_red_s
