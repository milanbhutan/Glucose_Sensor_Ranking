import numpy as np
import pandas as pd
import tissue_model.tissue_params as params
import tissue_model.layer as layer

rng = np.random.default_rng(45)

epi_mu_a_sim_set = []
derm_mu_a_sim_set = []
sub_mu_a_sim_set = []

for i in range(200):
    derm_blood = rng.uniform(params.derm_comp_lower["blood"], params.derm_comp_upper["blood"])
    derm_oxygen_saturation = rng.uniform(params.derm_comp_lower["oxygen_saturation"], params.derm_comp_upper["oxygen_saturation"])
    derm_water = rng.uniform(params.derm_comp_lower["water"], params.derm_comp_upper["water"])
    derm_gamma = rng.uniform(params.derm_comp_lower["gamma"], params.derm_comp_upper["gamma"])
    if(1 - (derm_blood + derm_water) > 0):
        derm_mu_a_sim_set.append(
            layer.mu_a_derm(
                B = derm_blood, 
                S = derm_oxygen_saturation, 
                W = derm_water, 
                F = 0, 
                M = 0,
                gamma = derm_gamma,
        )
        )

for i in range(200):
    sub_blood = rng.uniform(params.sub_comp_lower["blood"], params.sub_comp_upper["blood"])
    sub_oxygen_saturation = rng.uniform(params.sub_comp_lower["oxygen_saturation"], params.sub_comp_upper["oxygen_saturation"])
    sub_water = rng.uniform(params.sub_comp_lower["water"], params.sub_comp_upper["water"])
    sub_fat = rng.uniform(params.sub_comp_lower["fat"], params.sub_comp_upper["fat"])
    sub_gamma = rng.uniform(params.sub_comp_lower["gamma"], params.sub_comp_upper["gamma"])
    if(1 - (sub_blood + sub_water + sub_fat) > 0): # Remove samples with a negative baseline fraction
        sub_mu_a_sim_set.append(
            layer.mu_a_sub(
                B = sub_blood, 
                S = sub_oxygen_saturation, 
                W = sub_water, 
                F = sub_fat, 
                M = 0,
                gamma = sub_gamma
        )
        )

epi_mu_a_sim_df = []
derm_mu_a_sim_df = pd.concat(derm_mu_a_sim_set, axis = 1)
sub_mu_a_sim_df = pd.concat(sub_mu_a_sim_set, axis = 1)