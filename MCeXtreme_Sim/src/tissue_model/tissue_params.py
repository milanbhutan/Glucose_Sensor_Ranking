import numpy as np
import pandas as pd
import tissue_model.layer_calc as calc

derm_gamma = 1 # was 0
sub_gamma = 0

lower_wavelength = 600 # original: 600
upper_wavelength = 1060
norm_wavelength_nm = np.arange(0, 2001, 1)

baseline_epi_np = (8.53 * np.exp(-(norm_wavelength_nm - 154)/66.2) + 0.0244)
baseline_derm_np = 7.84 * 10**7 * norm_wavelength_nm ** -3.255
baseline_sub_np = (1.682 * np.exp(-(norm_wavelength_nm - 400)/80.5) + 0.0244)

baseline_epi = pd.Series(
    baseline_epi_np,
    index = norm_wavelength_nm
)

baseline_derm = pd.Series(
    baseline_derm_np,
    index = norm_wavelength_nm
)

baseline_sub = pd.Series(
    baseline_sub_np,
    index = norm_wavelength_nm
)

derm_params = [
    "blood",
    "oxygen_saturation",
    "water",
]

sub_params = [
    "blood",
    "oxygen_saturation",
    "water",
    "fat",
]


epi_Jacques = {
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0,
    "fat" : 0,
    "melanosome" : 0.025,
}
epi_Jacques["baseline"] = (
    1 - 
    epi_Jacques["blood"] - 
    epi_Jacques["water"] - 
    epi_Jacques["fat"] - 
    epi_Jacques["melanosome"]
    )


epi_Jolivot_low_mel ={
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0,
    "fat" : 0,
    "melanosome": 0.013
}
epi_Jolivot_low_mel["baseline"] = (
    1 - 
    epi_Jolivot_low_mel["blood"] - 
    epi_Jolivot_low_mel["water"] - 
    epi_Jolivot_low_mel["fat"] - 
    epi_Jolivot_low_mel["melanosome"]
    )


epi_Jolivot_SPT_II ={
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0,
    "fat" : 0,
    "melanosome": 0.101
}
epi_Jolivot_SPT_II["baseline"] = (
    1 - 
    epi_Jolivot_SPT_II["blood"] - 
    epi_Jolivot_SPT_II["water"] - 
    epi_Jolivot_SPT_II["fat"] - 
    epi_Jolivot_SPT_II["melanosome"]
    )


epi_Jolivot_SPT_III ={
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0,
    "fat" : 0,
    "melanosome": 0.187
}
epi_Jolivot_SPT_III["baseline"] = (
    1 - 
    epi_Jolivot_SPT_III["blood"] - 
    epi_Jolivot_SPT_III["water"] - 
    epi_Jolivot_SPT_III["fat"] - 
    epi_Jolivot_SPT_III["melanosome"]
    )


epi_Jolivot_SPT_IV ={
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0,
    "fat" : 0,
    "melanosome": 0.265
}
epi_Jolivot_SPT_IV["baseline"] = (
    1 - 
    epi_Jolivot_SPT_IV["blood"] - 
    epi_Jolivot_SPT_IV["water"] - 
    epi_Jolivot_SPT_IV["fat"] - 
    epi_Jolivot_SPT_IV["melanosome"]
    )


epi_Jolivot_SPT_V ={
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0,
    "fat" : 0,
    "melanosome": 0.358
}
epi_Jolivot_SPT_V["baseline"] = (
    1 - 
    epi_Jolivot_SPT_V["blood"] - 
    epi_Jolivot_SPT_V["water"] - 
    epi_Jolivot_SPT_V["fat"] - 
    epi_Jolivot_SPT_V["melanosome"]
    )


epi_Jolivot_SPT_VI = {
    "blood" : 0,
    "oxygen_saturation" : 0,
    "water" : 0,
    "fat" : 0,
    "melanosome": 0.447
}
epi_Jolivot_SPT_VI["baseline"] = (
    1 - 
    epi_Jolivot_SPT_VI["blood"] - 
    epi_Jolivot_SPT_VI["water"] - 
    epi_Jolivot_SPT_VI["fat"] - 
    epi_Jolivot_SPT_VI["melanosome"]
    )


derm_Jacques = {
    "blood" : 0.0020,  
    "oxygen_saturation" : 0.39,
    "water" : 0.65,
    "fat" : 0,
    "melanosome" : 0,
    "gamma" : derm_gamma
}
derm_Jacques["baseline"] = (
    1 - 
    derm_Jacques["blood"] - 
    derm_Jacques["water"] - 
    derm_Jacques["fat"] - 
    derm_Jacques["melanosome"]
    )


lower_derm_Jonasson = {
    "blood" : 0.0075,  
    "oxygen_saturation" : 0.50, 
    "water" : 0.0, # Jonasson neglects water
    "fat" : 0,
    "melanosome" : 0,
    "gamma" : 0,
}
lower_derm_Jonasson["baseline"] = (
    1 - 
    lower_derm_Jonasson["blood"] - 
    lower_derm_Jonasson["water"] - 
    lower_derm_Jonasson["fat"] - 
    lower_derm_Jonasson["melanosome"]
    )


sub_Jacques = {
    "blood" : 0.0054,
    "oxygen_saturation" : 0.76,
    "water" : 0.11,
    "fat" : 0.69,
    "melanosome" : 0,
    "gamma" : derm_gamma,
}
sub_Jacques["baseline"] = (
    1 - 
    sub_Jacques["blood"] - 
    sub_Jacques["water"] - 
    sub_Jacques["fat"] - 
    sub_Jacques["melanosome"]
    )


sub_Ganesan = {
    "blood" : 0.00634,
    "oxygen_saturation" : 0.60,
    "water" : 0.2184,
    "fat" : 0.6433,
    "melanosome" : 0,
    "gamma" : sub_gamma,
}
sub_Ganesan["baseline"] = (
    1 - 
    sub_Ganesan["blood"] - 
    sub_Ganesan["water"] - 
    sub_Ganesan["fat"] - 
    sub_Ganesan["melanosome"]
    )


derm_upper_bound = {
    "blood" : 0.015,
    "oxygen_saturation" : 0.65,    
    "water" : 0.75,
    "gamma" : 1 # original: 0.8

}

derm_lower_bound = {
    "blood" : 0.002,
    "oxygen_saturation" : 0.39,    
    "water" : 0.65,
    "gamma" : 0.1 # original: 0.2

}

sub_upper_bound = {
    "blood" : 0.00813,
    "oxygen_saturation" : 0.67,    
    "water" : 0.304,
    "fat" : 0.702,
    "gamma" : 0, # original: 0.3

}

sub_lower_bound = {
    "blood" : 0.00398,
    "oxygen_saturation" : 0.493,    
    "water" : 0.145,
    "fat" : 0.599, 
    "gamma" : 0, #original: 0.2

}

epi_mu_a_sim, epi_mu_s_sim, epi_g_sim, epi_n_sim, epi_mu_s_red_sim = calc.epi(
    epi_Jacques["blood"], 
    epi_Jacques["oxygen_saturation"], 
    epi_Jacques["water"], 
    epi_Jacques["fat"], 
    epi_Jacques["melanosome"],
    baseline_epi,
    )

epi_mu_a_sim_Jolivot_SPT_II = calc.epi(
    epi_Jolivot_SPT_II["blood"], 
    epi_Jolivot_SPT_II["oxygen_saturation"], 
    epi_Jolivot_SPT_II["water"], 
    epi_Jolivot_SPT_II["fat"], 
    epi_Jolivot_SPT_II["melanosome"],
    baseline_epi,
    )[0]

epi_mu_a_sim_Jolivot_SPT_III = calc.epi(
    epi_Jolivot_SPT_III["blood"], 
    epi_Jolivot_SPT_III["oxygen_saturation"], 
    epi_Jolivot_SPT_III["water"], 
    epi_Jolivot_SPT_III["fat"], 
    epi_Jolivot_SPT_III["melanosome"],
    baseline_epi,
    )[0]

epi_mu_a_sim_Jolivot_SPT_IV = calc.epi(
    epi_Jolivot_SPT_IV["blood"], 
    epi_Jolivot_SPT_IV["oxygen_saturation"], 
    epi_Jolivot_SPT_IV["water"], 
    epi_Jolivot_SPT_IV["fat"], 
    epi_Jolivot_SPT_IV["melanosome"],
    baseline_epi
    )[0]

epi_mu_a_sim_Jolivot_SPT_V = calc.epi(
    epi_Jolivot_SPT_V["blood"], 
    epi_Jolivot_SPT_V["oxygen_saturation"], 
    epi_Jolivot_SPT_V["water"], 
    epi_Jolivot_SPT_V["fat"], 
    epi_Jolivot_SPT_V["melanosome"],
    baseline_epi
    )[0]

epi_mu_a_sim_Jolivot_SPT_VI = calc.epi(
    epi_Jolivot_SPT_VI["blood"], 
    epi_Jolivot_SPT_VI["oxygen_saturation"], 
    epi_Jolivot_SPT_VI["water"], 
    epi_Jolivot_SPT_VI["fat"], 
    epi_Jolivot_SPT_VI["melanosome"],
    baseline_epi
    )[0]

derm_mu_a_sim, derm_mu_s_sim, derm_g_sim, derm_n_sim, derm_mu_s_red_sim = calc.derm(
    derm_Jacques["blood"], 
    derm_Jacques["oxygen_saturation"], 
    derm_Jacques["water"], 
    derm_Jacques["fat"], 
    derm_Jacques["melanosome"],
    derm_Jacques["gamma"],
    baseline_derm,
    )

Jonasson_derm_mu_a_sim = calc.derm(
    lower_derm_Jonasson["blood"], 
    lower_derm_Jonasson["oxygen_saturation"], 
    lower_derm_Jonasson["water"], 
    lower_derm_Jonasson["fat"], 
    lower_derm_Jonasson["melanosome"],
    lower_derm_Jonasson["gamma"],
    baseline_derm,
    )[0]

sub_mu_a_sim, sub_mu_s_sim, sub_g_sim, sub_n_sim, sub_mu_s_red_sim = calc.sub(
    sub_Ganesan["blood"], 
    sub_Ganesan["oxygen_saturation"],
    sub_Ganesan["water"],
    sub_Ganesan["fat"],
    sub_Ganesan["melanosome"],
    baseline_sub
    )