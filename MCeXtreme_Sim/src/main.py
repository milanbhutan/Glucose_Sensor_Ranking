import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pysr import PySRRegressor as SR
import physo
import torch
from scipy.optimize import least_squares
from pathlib import Path

import tissue_model.layer_calc as calc
import tissue_model.tissue_params as params
import analysis.error_analysis as err
import analysis.uncertainty_analysis as uncertain
import photon_sim.monte_carlo as mc



plt.rcParams["text.usetex"] = False

csv_path = Path("outputs/PySR_Best_eq.csv")

PySR_runs = range(20)





# Obtain data for experimental absorption and scattering spectras for three-layer tissue model from Shimojo 2020

# Absorption Data
df_Wan_epi_mu_a = (
    pd.read_csv("Exp_Absorption/wan_epidermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Shimojo_epi_mu_a = (
    pd.read_csv("Exp_Absorption/shimojo_epidermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)


df_Simpson_derm_mu_a = (
    pd.read_csv("Exp_Absorption/simpson_dermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
mu_a_derm_Simpson = df_Simpson_derm_mu_a.loc[params.lower_wavelength:params.upper_wavelength, "caucasian_dermis_mu_a_mm^-1"]

df_Shimojo_derm_mu_a = (
    pd.read_csv("Exp_Absorption/shimojo_dermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Jonasson_derm_mu_a = (
    pd.read_csv("Exp_Absorption/jonasson_lower_dermis_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Ganesan_sub_mu_a = (
    pd.read_csv("Exp_Absorption/ganesan_subcutaneous_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
mu_a_sub_Ganesan = df_Ganesan_sub_mu_a.loc[params.lower_wavelength:params.upper_wavelength, "mu_a_mm^-1"]

df_Simpson_sub_mu_a = (
    pd.read_csv("Exp_Absorption/simpson_subcutaneous_absorption.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)
mu_a_sub_Simpson = df_Simpson_sub_mu_a.loc[params.lower_wavelength:params.upper_wavelength, "caucasian_subdermis_mu_a_mm^-1"]

# Scattering Data
df_Shimojo_epi_red_mu_s = (
    pd.read_csv("Exp_Scattering/shimojo_epidermis_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Shimojo_derm_red_mu_s = (
    pd.read_csv("Exp_Scattering/shimojo_dermis_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Simpson_derm_red_mu_s = (
    pd.read_csv("Exp_Scattering/simpson_dermis_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Shimojo_sub_red_mu_s = (
    pd.read_csv("Exp_Scattering/shimojo_subcutaneous_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Ganesan_sub_red_mu_s = (
    pd.read_csv("Exp_Scattering/ganesan_subcutaneous_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

df_Simpson_sub_red_mu_s = (
    pd.read_csv("Exp_Scattering/simpson_subcutaneous_reduced_scattering.csv", comment = "#")
    .set_index("wavelength_nm")
    .sort_index()
)

# ===============================================
# Plotting tissue spectra from 600 - 1100nm
# ===============================================

# Note: This subcutaneous spectra is only valid within 429-1100 nm where reliable tissue constituent abosroption data is availible
fig_mu_a, axes_mu_a = plt.subplots(1, 3, layout = "constrained")
fig_mu_a.get_layout_engine().set(h_pad = 0.1, w_pad = 0.1)
fig_mu_a.suptitle("Absorption Spectras for Layers in Three-Layer Tissue Model", fontsize = 18)

for ax in axes_mu_a:
    ax.set_xlabel("Wavelength (nm)", fontsize = 18)
    ax.set_ylabel("Absorption (1/mm)", fontsize = 18)
    ax.tick_params(axis = "both", labelsize = 18)


# Simulated absorption spectra plots
axes_mu_a[0].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
                  params.epi_mu_a_sim[params.lower_wavelength:params.upper_wavelength], 
                  label = "Analytical Model", 
                  color = "red")
axes_mu_a[0].text(600, 
                  0.2, 
                  "M=2.5%", 
                  color = "red",
                  fontsize = 10)

axes_mu_a[0].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
                  params.epi_mu_a_sim_Jolivot_SPT_II[params.lower_wavelength:params.upper_wavelength], 
                  color = "red")


axes_mu_a[0].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
                  params.epi_mu_a_sim_Jolivot_SPT_III[params.lower_wavelength:params.upper_wavelength], 
                  color = "red")

axes_mu_a[0].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
                  params.epi_mu_a_sim_Jolivot_SPT_IV[params.lower_wavelength:params.upper_wavelength], 
                  color = "red")

axes_mu_a[0].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
                  params.epi_mu_a_sim_Jolivot_SPT_V[params.lower_wavelength:params.upper_wavelength], 
                  color = "red")

axes_mu_a[0].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
                  params.epi_mu_a_sim_Jolivot_SPT_VI[params.lower_wavelength:params.upper_wavelength], 
                  color = "red")
axes_mu_a[0].text(713.5, 
                  9.11, 
                  "M=44.7%", 
                  color = "red",
                  fontsize = 10)

axes_mu_a[0].set_title(r"$\mu_a$ epidermis", fontsize = 18)


axes_mu_a[1].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
                  params.derm_mu_a_sim[params.lower_wavelength:params.upper_wavelength], 
                  color = "red",
                  label = "Analytical Model")

axes_mu_a[1].plot(params.norm_wavelength_nm[params.lower_wavelength:850], 
                  params.Jonasson_derm_mu_a_sim[params.lower_wavelength:850], 
                  color = "red",
                  )
axes_mu_a[1].text(850, 
                  0, 
                  "Simplified dermis\nabsorption model\n", 
                  color = "red",
                  fontsize = 10)

axes_mu_a[1].plot(
    uncertain.derm_mu_a_sim_df[params.lower_wavelength:params.upper_wavelength], 
    color = "red",
    alpha = 0.02
    )
# axes_mu_a[1].plot(
#     uncertain.Jonasson_derm_mu_a_sim_df[params.lower_wavelength:params.upper_wavelength], 
#     color = "red",
#     alpha = 0.02
#     )
axes_mu_a[1].set_title(r"$\mu_a$ dermis",  fontsize = 18)



axes_mu_a[2].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
                  params.sub_mu_a_sim[params.lower_wavelength:params.upper_wavelength], 
                  color = "red",
                  label = "Analytical Model")

axes_mu_a[2].plot(
    uncertain.sub_mu_a_sim_df[params.lower_wavelength:params.upper_wavelength], 
    color = "red",
    alpha = 0.02
    )

axes_mu_a[2].set_title(r"$\mu_a$ subcutaneous", fontsize = 18)


# Experimental absorption spectra plots
axes_mu_a[0].plot(df_Wan_epi_mu_a.loc[params.lower_wavelength:params.upper_wavelength], 
                  color = "orange", 
                  label = "Wan 1981: Dark Skin",
                  linestyle = "--")
axes_mu_a[0].plot(df_Shimojo_epi_mu_a.loc[params.lower_wavelength:params.upper_wavelength], 
                  color = "black", 
                  label = "Shimojo 2020: Asian Skin",
                  linestyle = "--")

axes_mu_a[1].plot(df_Simpson_derm_mu_a.loc[params.lower_wavelength:params.upper_wavelength], 
                  color = "green",
                  label = "Simpson 1998",
                  linestyle = "--")

axes_mu_a[1].plot(df_Jonasson_derm_mu_a.loc[params.lower_wavelength:params.upper_wavelength], 
                  color = "purple",
                  label = "Jonasson 2023",
                  linestyle = "--")

# axes_mu_a[1].plot(df_Shimojo_derm_mu_a.loc[params.lower_wavelength:params.upper_wavelength], 
#                   color = "black",
#                   linestyle = "--")

axes_mu_a[2].plot(df_Ganesan_sub_mu_a.loc[params.lower_wavelength:params.upper_wavelength], 
                  color = "black", 
                  label = "Ganesan 2016", 
                  linestyle = "--")
axes_mu_a[2].plot(df_Simpson_sub_mu_a.loc[params.lower_wavelength:params.upper_wavelength], 
                  color = "green", 
                  label = "Simpson 1998", 
                  linestyle = "--")
axes_mu_a[2].scatter(x = 650, 
                     y = 0.006,
                     color = "pink",
                     label = "Weber 2009",
                     s = 100)

axes_mu_a[0].legend(fontsize = 10)
axes_mu_a[1].legend(fontsize = 10)
axes_mu_a[2].legend(fontsize = 10)

fig_mu_s, axes_mu_s = plt.subplots(1, 3, layout = "constrained")
fig_mu_s.get_layout_engine().set(h_pad = 0.1, w_pad = 0.1)
fig_mu_s.suptitle("Reduced Scattering Spectras for Layers in Three-Layer Tissue Model", fontsize = 18)

for ax in axes_mu_s:
    ax.set_xlabel("Wavelength (nm)", fontsize = 18)
    ax.set_ylabel("Reduced Scattering (1/mm)", fontsize = 18)
    ax.tick_params(axis = "both", labelsize = 18)

# Simulated reduced scattering coefficients
axes_mu_s[0].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
                  params.epi_mu_s_red_sim[params.lower_wavelength:params.upper_wavelength], 
                  label = "Jacques 2013 Model", 
                  color = "red")
axes_mu_s[0].set_title(r"$\mu_s^'$ epidermis", fontsize = 18)



axes_mu_s[1].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
                  params.derm_mu_s_red_sim[params.lower_wavelength:params.upper_wavelength], 
                  label = "Jacques 2013 Model",
                  color = "red")
axes_mu_s[1].set_title(r"$\mu^'_s$ dermis", fontsize = 18)

axes_mu_s[2].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
                  params.sub_mu_s_red_sim[params.lower_wavelength:params.upper_wavelength],
                  label = "Jacques 2013 Model", 
                  color = "red")
axes_mu_s[2].set_title(r"$\mu^'_s$ subcutaneous tissue", fontsize = 18)

# Experimental reduced scattering coefficients
axes_mu_s[0].plot(df_Shimojo_epi_red_mu_s.loc[params.lower_wavelength:params.upper_wavelength], 
                  color = "cyan", 
                  label = "Shimojo 2020", 
                  linestyle = "--")

axes_mu_s[1].plot(df_Shimojo_derm_red_mu_s.loc[params.lower_wavelength:params.upper_wavelength], 
                  color = "cyan",
                  label = "Shimojo 2020", 
                  linestyle = "--")
axes_mu_s[1].plot(df_Simpson_derm_red_mu_s.loc[params.lower_wavelength:params.upper_wavelength], 
                  color = "green", 
                  label = "Simpson 1998", 
                  linestyle = "--")

axes_mu_s[2].plot(df_Shimojo_sub_red_mu_s.loc[params.lower_wavelength:params.upper_wavelength], 
                  color = "cyan",
                  label = "Shimojo 2020", 
                  linestyle = "--")
axes_mu_s[2].plot(df_Simpson_sub_red_mu_s.loc[params.lower_wavelength:params.upper_wavelength], 
                  color = "green",
                  label = "Simpson 1998",
                  linestyle = "--")
axes_mu_s[2].plot(df_Ganesan_sub_red_mu_s.loc[params.lower_wavelength:params.upper_wavelength], 
                  color = "purple", 
                  label = "Ganesan 2016", 
                  linestyle = "--")

axes_mu_s[0].legend(fontsize = 10)
axes_mu_s[1].legend(fontsize = 10)
axes_mu_s[2].legend(fontsize = 10)


# print(f"Maximum mu_a_derm_660nm: {max(mu_a_derm_660)}")
# print(f"Minimum mu_a_derm_660nm: {min(mu_a_derm_660)}")
# print(f"Maximum mu_a_derm_940nm: {max(mu_a_derm_940)}")
# print(f"Minimum mu_a_derm_940nm: {min(mu_a_derm_940)}")
# print(f"Maximum mu_a_sub_650nm: {max(mu_a_sub_650)}")
# print(f"Minimum mu_a_sub_650nm: {min(mu_a_sub_650)}")
# print(f"Maximum mu_a_sub_RMSE: {max(mu_a_sub_RMSEs)}")
# print(f"Minimum mu_a_sub_RMSE: {min(mu_a_sub_RMSEs)}")

# # Find tissue compostion which matches Ganesan 2016
# result_Ganesan = least_squares(
#     err.residual_sub, 
#     [params.sub_Ganesan["blood"],
#      params.sub_Ganesan["oxygen_saturation"],
#      params.sub_Ganesan["water"],
#      params.sub_Ganesan["fat"],
#      ],
#     args = (mu_a_sub_Ganesan, calc.mu_a_sub,),
#     bounds =([params.sub_lower_bound[p] for p in params.sub_params],
#              [params.sub_upper_bound[p] for p in params.sub_params])
#     )


# mu_a_sub_Ganesan_fitted = calc.mu_a_sub(
#     B = result_Ganesan.x[0], 
#     S = result_Ganesan.x[1],
#     W = result_Ganesan.x[2],
#     F = result_Ganesan.x[3],
#     M = 0.0,
#     )

# axes_mu_a[2].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
#                   mu_a_sub_Ganesan_fitted[params.lower_wavelength:params.upper_wavelength],
#                   label = "Model Prediction",
#                   color = "brown")

# # Find tissue compostion which matches Simpson 1998, subcutaneous layer
# # Simpson 1998 compostion not known so using Ganesan 2016 composition as an initial guess
# result_Simpson_sub = least_squares(
#     err.residual_sub, 
#     [params.sub_Ganesan["blood"],
#      params.sub_Ganesan["oxygen_saturation"],
#      params.sub_Ganesan["water"],
#      params.sub_Ganesan["fat"],
#      0.25
#      ],
#     args = (mu_a_sub_Simpson, calc.mu_a_sub,),
#     bounds =([params.sub_lower_bound[p] for p in params.sub_params],
#              [params.sub_upper_bound[p] for p in params.sub_params])
#     )

# mu_a_sub_Simpson_fitted = calc.mu_a_sub(
#     B = result_Simpson_sub.x[0], 
#     S = result_Simpson_sub.x[1],
#     W = result_Simpson_sub.x[2],
#     F = result_Simpson_sub.x[3],
#     M = 0.0,
#     gamma = result_Simpson_sub.x[4]
#     )

# axes_mu_a[2].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
#                   mu_a_sub_Simpson_fitted[params.lower_wavelength:params.upper_wavelength],
#                   color = "brown")


# # Find tissue compostion which matches Simpson 1998, dermis layer
# # Simpson 1998 compostion not know so using Jacques 2013 composition as an initial guess
# result_Simpson_derm = least_squares(
#     err.residual_derm, 
#     [params.derm_Jacques["blood"],
#      params.derm_Jacques["oxygen_saturation"],
#      params.derm_Jacques["water"],
#      params.derm_Jacques["gamma"],
#      ],
#     args = (mu_a_derm_Simpson, calc.mu_a_derm,),
#     bounds =([params.derm_lower_bound[p] for p in params.derm_params],
#              [params.derm_upper_bound[p] for p in params.derm_params])
#     )

# mu_a_derm_Simpson_fitted = calc.mu_a_derm(
#     B = result_Simpson_derm.x[0], 
#     S = result_Simpson_derm.x[1],
#     W = result_Simpson_derm.x[2],
#     F = 0.0,
#     M = 0.0,
#     gamma = result_Simpson_derm.x[3]
#     )

# axes_mu_a[1].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
#                   mu_a_derm_Simpson_fitted[params.lower_wavelength:params.upper_wavelength],
#                   color = "brown")

#fig_mu_a.legend(fontsize = 10) was here but moved up


# 930 nm lipid absorption feature
axes_mu_a[2].vlines(
    x = 930,
    ymin = 0.01,
    ymax = 0.020,
    linestyle=":",
    linewidth = 1.5,
    color = "black"
)

axes_mu_a[2].text(
    900,
    0.020,
    "930 nm\nC-H overtone",
    ha = "center",
    va = "bottom",
    fontsize = 10
)

# 970 nm water absorption feature
axes_mu_a[1].vlines(
    x = 970,
    ymin = 0.027,
    ymax = 0.05,
    linestyle=":",
    linewidth = 1.5,
    color = "black"
)

axes_mu_a[1].text(
    970,
    0.05,
    "970 nm\nH$_2$O O-H overtone",
    ha = "center",
    va = "bottom",
    fontsize = 10
)

axes_mu_a[2].vlines(
    x = 970,
    ymin = 0.01,
    ymax = 0.025,
    linestyle=":",
    linewidth = 1.5,
    color = "black"
)

axes_mu_a[2].text(
    970,
    0.025,
    "970 nm\nH$_2$O O-H overtone",
    ha = "center",
    va = "bottom",
    fontsize = 10
)

# Deoxy Hb absorption
axes_mu_a[1].vlines(
    x = 760,
    ymin = 0.002,
    ymax = 0.022,
    linestyle=":",
    linewidth=1.5,
    color="black"
)

axes_mu_a[1].text(
    760,
    0.022,
    "760 nm\nHb Electron Transfer",
    ha = "center",
    va = "bottom",
    fontsize = 10
)

axes_mu_a[2].vlines(
    x = 760,
    ymin = 0.0025,
    ymax = 0.014,
    linestyle=":",
    linewidth=1.5,
    color="black"
)

axes_mu_a[2].text(
    760,
    0.014,
    "760 nm\nHb Electron Transfer",
    ha = "center",
    va = "bottom",
    fontsize = 10
)

# # Symbolic Regression using Ganesan 2016
# mu_a_sub_Ganesan_align, sub_mu_a_sim_align = mu_a_sub_Ganesan.align(params.sub_mu_a_sim[params.lower_wavelength:params.upper_wavelength], join = "inner")
# mu_a_sub_Ganesan_align, baseline_sub_align = mu_a_sub_Ganesan.align(params.baseline_sub[params.lower_wavelength:params.upper_wavelength], join = "inner")

# alpha = (mu_a_sub_Ganesan_align -  sub_mu_a_sim_align) / (params.sub_Ganesan["baseline"] * baseline_sub_align)
# beta = (mu_a_sub_Ganesan_align -  sub_mu_a_sim_align) / params.sub_Ganesan["baseline"]

# fig_scale, axes_scale = plt.subplots(1, 2)

# fig_scale.suptitle("Baseline discrepancy between analytical and experimental\n for $\mu_{a,subcutaneous}$ , Ganesan 2016")

# axes_scale[0].plot(alpha.index, 
#                   alpha.loc[params.lower_wavelength:params.upper_wavelength], 
#                   color = "black")
# axes_scale[1].plot(beta.index, 
#                   beta.loc[params.lower_wavelength:params.upper_wavelength], 
#                   color = "black")


# fig_symb, axes_symb = plt.subplots(1, 2)

# fig_symb.suptitle(r"Physically Implausible $\alpha(\lambda)$ Expressions from PySR")

# Will test later
# expression, logs = physo.SR(
#     binary_operators = ["*", "+", "-"],
#     unary_operators = ["exp"]
# )

# for i in PySR_runs:
#     model = SR(
#         binary_operators = ["*", "+", "-", "/"],
#         unary_operators = ["exp"],
#         elementwise_loss = """
# loss(prediction, target) = prediction >= 0 ? (target - prediction) ^ 2 : Inf
# """
#     )
#     model.fit(alpha.index.to_numpy().reshape(-1, 1), alpha.loc[params.lower_wavelength:params.upper_wavelength])

#     print(model)

#     best = model.get_best()
#     best.to_frame().T.to_csv(
#         csv_path,
#         mode = "a",
#         header = not csv_path.exists(),
#         index = False
        
#         )


# correction_1 = np.exp(params.norm_wavelength_nm * 0.008014623) * -0.00020990851

# print(f"calc input: {params.baseline_sub}")

# cor_sub_mu_a_sim = calc.sub(params.sub_Ganesan["blood"], 
#                                                params.sub_Ganesan["oxygen_saturation"],
#                                                params.sub_Ganesan["water"],
#                                                params.sub_Ganesan["fat"],
#                                                params.sub_Ganesan["melanosome"],
#                                                params.baseline_sub,
#                                                gamma = correction_1
#                                                )[0]



# axes_symb[0].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
#                   params.sub_mu_a_sim[params.lower_wavelength:params.upper_wavelength], 
#                   color = "red",
#                   label = "Analytical Model",)
# axes_symb[0].plot(df_Ganesan_sub_mu_a.loc[params.lower_wavelength:params.upper_wavelength], 
#                   color = "black", 
#                   label = "Ganesan 2016", 
#                   linestyle = "--")
# axes_symb[0].set_title(r"$\alpha(\lambda) = 0$")

# axes_symb[1].plot(params.norm_wavelength_nm[params.lower_wavelength:params.upper_wavelength], 
#                   cor_sub_mu_a_sim[params.lower_wavelength:params.upper_wavelength], 
#                   color = "red")
# axes_symb[1].plot(df_Ganesan_sub_mu_a.loc[params.lower_wavelength:params.upper_wavelength], 
#                   color = "black", 
#                   linestyle = "--")
# axes_symb[1].set_title(r"$\alpha(\lambda) < 0$")


# fig_symb.legend(fontsize = 10)



# print(
#     f"===Tissue Composition Prediction of Ganesan 2016===\n"
#     f"Blood: {result_Ganesan.x[0]:.4f}\n"
#     f"Oxygen Saturation: {result_Ganesan.x[1]:.4f}\n"
#     f"Water: {result_Ganesan.x[2]:.4f}\n"
#     f"Fat: {result_Ganesan.x[3]:.4f}\n"
#     f"Melanin: 0.0\n"
#     )

# print(
#     f"===Tissue Composition Prediction of Simpson 1998, Subcutaneous===\n"
#     f"Blood: {result_Simpson_sub.x[0]:.4f}\n"
#     f"Oxygen Saturation: {result_Simpson_sub.x[1]:.4f}\n"
#     f"Water: {result_Simpson_sub.x[2]:.4f}\n"
#     f"Fat: {result_Simpson_sub.x[3]:.4f}\n"
#     f"Melanin: 0.0\n"
#     f"Gamma: {result_Simpson_sub.x[4]:.4f}\n"
#     )

# print(
#     f"===Tissue Composition Prediction of Simpson 1998, Dermis===\n"
#     f"Blood: {result_Simpson_derm.x[0]:.4f}\n"
#     f"Oxygen Saturation: {result_Simpson_derm.x[1]:.4f}\n"
#     f"Water: {result_Simpson_derm.x[2]:.4f}\n"
#     f"Melanin: 0.0\n"
#     f"Gamma: {result_Simpson_derm.x[3]:.4f}\n"
#     )


# uncomment to run monte carlo simulation


sim_wavelength_nm = 1050
canidate2_wavelength_nm = 1050

mu_s_values = np.arange(5, 14, 1)
A = []
for mu_s in mu_s_values:
    val_A = mc.run_simulation(
            [
            params.epi_mu_a_sim[sim_wavelength_nm], 
            params.epi_mu_s_sim[sim_wavelength_nm], 
            params.epi_g_sim[sim_wavelength_nm], 
            params.epi_n_sim ],

            [
            params.derm_mu_a_sim[sim_wavelength_nm], 
            mu_s, 
            params.derm_g_sim[sim_wavelength_nm], 
            params.derm_n_sim],

            [
            params.sub_mu_a_sim[sim_wavelength_nm], 
            params.sub_mu_s_sim[sim_wavelength_nm], 
            params.sub_g_sim[sim_wavelength_nm], 
            params.sub_n_sim]
    )
    A.append(val_A)




# np.t
# np.corrcoef
plt.figure()
plt.scatter(mu_s_values, A)
plt.xlabel("$\mu_s$ $(mm^-1)$")
plt.ylabel("$A$ $(a.u.)$")

m, b = np.polyfit(mu_s_values, A, 1)
print(m)

plt.show()