import pandas as pd
import matplotlib.pyplot as plt
import tissue_model.layer_calc as calc
import tissue_model.tissue_params as params
from paths import MCX_DIR
from paths import sim_results_DIR

mean_path_epi_df = pd.read_csv(MCX_DIR / sim_results_DIR / "mean_metrics"/ "Metric_Mean_mean_path_epi.csv", index_col = "Wavelength (nm)")
mean_path_derm_df = pd.read_csv(MCX_DIR / sim_results_DIR / "mean_metrics"/ "Metric_Mean_mean_path_derm.csv", index_col = "Wavelength (nm)")
mean_path_sub_df = pd.read_csv(MCX_DIR / sim_results_DIR / "mean_metrics"/ "Metric_Mean_mean_path_sub.csv", index_col = "Wavelength (nm)")
dA_dmu_a_df = pd.read_csv(MCX_DIR / sim_results_DIR / "mean_metrics"/ "Metric_Mean_dA_dmu_a.csv", index_col = "Wavelength (nm)")
dA_dmu_s_df = pd.read_csv(MCX_DIR / sim_results_DIR / "mean_metrics"/ "Metric_Mean_dA_dmu_s.csv", index_col = "Wavelength (nm)")
dA_dCg_mu_a_df = pd.read_csv(MCX_DIR / sim_results_DIR / "mean_metrics"/ "Metric_Mean_dA_dCg_mu_a.csv", index_col = "Wavelength (nm)")
dA_dCg_mu_s_df = pd.read_csv(MCX_DIR / sim_results_DIR / "mean_metrics"/ "Metric_Mean_dA_dCg_mu_s.csv", index_col = "Wavelength (nm)")
dA_dCg_df = pd.read_csv(MCX_DIR / sim_results_DIR / "mean_metrics"/ "Metric_Mean_dA_dCg_total.csv", index_col = "Wavelength (nm)")
beta_derm_df = 100 * pd.read_csv(MCX_DIR / sim_results_DIR / "mean_metrics"/ "Metric_Mean_beta_derm.csv", index_col = "Wavelength (nm)")
power_dB_df = pd.read_csv(MCX_DIR / sim_results_DIR / "mean_metrics"/ "Metric_Mean_power_dB.csv", index_col = "Wavelength (nm)")

#Calculate Water Interference per layer
dW = 0.05
epi_dmu_a_dW = (
    - params.epi_mu_a_sim 
    + calc.epi(
        params.epi_Jacques["blood"], 
        params.epi_Jacques["oxygen_saturation"], 
        params.epi_Jacques["water"] + dW, 
        params.epi_Jacques["fat"], 
        params.epi_Jacques["melanosome"],
        params.baseline_epi,
    )[0]) / dW

derm_dmu_a_dW = (
    - params.derm_mu_a_sim 
    + calc.derm(
        params.derm_Jacques["blood"], 
        params.derm_Jacques["oxygen_saturation"], 
        params.derm_Jacques["water"] + dW, 
        params.derm_Jacques["fat"], 
        params.derm_Jacques["melanosome"],
        params.baseline_derm,
    )[0]) / dW

sub_dmu_a_dW = (
    - params.sub_mu_a_sim 
    + calc.sub(
        params.sub_Ganesan["blood"], 
        params.sub_Ganesan["oxygen_saturation"], 
        params.sub_Ganesan["water"] + dW, 
        params.sub_Ganesan["fat"], 
        params.sub_Ganesan["melanosome"],
        0, # No baseline applied by Lapitan 2024
    )[0]) / dW

# chain rule: dA/dmu_a * dmu_a/dW = dA/dW
epi_dmu_a_dW = epi_dmu_a_dW.reindex(
    mean_path_epi_df.index
)

derm_dmu_a_dW = derm_dmu_a_dW.reindex(
    mean_path_derm_df.index
)

sub_dmu_a_dW = sub_dmu_a_dW.reindex(
    mean_path_sub_df.index
)

epi_dA_dW = mean_path_epi_df.mul(
    epi_dmu_a_dW,
    axis = 0)
derm_dA_dW = mean_path_derm_df.mul(
    derm_dmu_a_dW,
    axis = 0)
sub_dA_dW = mean_path_sub_df.mul(
    sub_dmu_a_dW,
    axis = 0)

dA_dW_df = epi_dA_dW + derm_dA_dW + sub_dA_dW




fig_gluc, axes_gluc = plt.subplots(2, 3, layout = "constrained")

dA_dmu_a_df.plot(ax = axes_gluc[0, 0])
axes_gluc[0, 0].set_ylabel("$\delta A/\delta\mu_a$", fontsize = 18)
axes_gluc[0, 0].set_xlabel("")
axes_gluc[0, 0].tick_params(axis = "x", labelbottom = False)

dA_dmu_s_df.plot(ax = axes_gluc[1, 0], legend = False)
axes_gluc[1, 0].set_ylabel("$\delta A/\delta\mu_s$", fontsize = 18)
axes_gluc[1, 0].set_xlabel("Wavelength (nm)", fontsize = 18)

dA_dCg_mu_a_df.plot(ax = axes_gluc[0, 1], legend = False)
axes_gluc[0, 1].set_ylabel("$(\delta A/\delta\mu_a)*(d\mu_a /dC_g)$", fontsize = 18)
axes_gluc[0, 1].set_xlabel("")
axes_gluc[0, 1].tick_params(axis = "x", labelbottom = False)


dA_dCg_mu_s_df.plot(ax = axes_gluc[1, 1], legend = False)
axes_gluc[1, 1].set_ylabel("$(\delta A/\delta\mu_s)*(d\mu_s /dC_g)$", fontsize = 18)
axes_gluc[1, 1].set_xlabel("Wavelength (nm)", fontsize = 18)


dA_dCg_df.plot(ax = axes_gluc[1, 2], legend = False)
axes_gluc[1, 2].set_ylabel("$dA/dC_g$", fontsize = 18)
axes_gluc[1, 2].set_xlabel("Wavelength (nm)", fontsize = 18)

fig_gluc.delaxes(axes_gluc[0, 2])

# keys = ["dA_dmu_a","dA_dmu_s","dA_dCg_mu_a","dA_dCg_mu_s","dA_dCg_total"]
# for ax, key in zip(axes_gluc.flat, keys):
#     ax.errorbar(
#         x = mean_df[key].index,
#         y = mean_df[key]["SDS_1.7_mm"],
#         width = stderr_df[key],
#         color = "blue"

#     )

for ax in axes_gluc.flat:
    ax.tick_params(axis = "both", labelsize = 18)
    ax.grid(True)


def normalize(df):
    minimum = df.min().min()
    maximum = df.max().max()


    return (df - minimum) / (maximum - minimum)



glucose_score = normalize(abs(dA_dCg_df))
dermal_score = normalize(beta_derm_df)
power_score = normalize(power_dB_df)
water_score = normalize(abs(dA_dW_df))

fig_WPI, ax_WPI = plt.subplots(2, 4, layout = "constrained")
fig_WPI.suptitle("WPI Inputs", fontsize = 18)

dA_dCg_df.plot(ax = ax_WPI[0, 0])
ax_WPI[0 ,0].set_title("Raw Input 1: Glucose Sensitivity")
ax_WPI[0, 0].set_xlabel("")
ax_WPI[0, 0].set_ylabel("$dA/dC_g$", fontsize = 18)
ax_WPI[0, 0].tick_params(axis = "both", labelsize = 18)
ax_WPI[0, 0].tick_params(axis = "x", labelbottom = False)

beta_derm_df.plot(ax = ax_WPI[0, 1], legend = False)
ax_WPI[0, 1].set_title("Raw Input 2: Dermal Aborption Sensitivity")
ax_WPI[0, 1].set_xlabel("")
ax_WPI[0, 1].set_ylabel("Path Fraction (%)", fontsize = 18)
ax_WPI[0, 1].tick_params(axis = "both", labelsize = 18)
ax_WPI[0, 1].tick_params(axis = "x", labelbottom = False)

power_dB_df.plot(ax = ax_WPI[0, 2], legend = False)
ax_WPI[0, 2].set_title("Raw Input 3: Detected Optical Power")
ax_WPI[0, 2].set_xlabel("")
ax_WPI[0, 2].set_ylabel("Power (dB)", fontsize = 18)
ax_WPI[0, 2].tick_params(axis = "both", labelsize = 18)
ax_WPI[0, 2].tick_params(axis = "x", labelbottom = False)

dA_dW_df.plot(ax = ax_WPI[0, 3], legend = False)
ax_WPI[0, 3].set_title("Raw Input 4: Water Fraction Sensitivity")
ax_WPI[0, 3].set_xlabel("")
ax_WPI[0, 3].set_ylabel("$dA/dW$", fontsize = 18)
ax_WPI[0, 3].tick_params(axis = "x", labelbottom = False)


glucose_score.plot(ax = ax_WPI[1, 0], legend = False)
ax_WPI[1, 0].set_title("Normalized Input 1: Glucose Sensitivity")
ax_WPI[1, 0].set_xlabel("Wavelength (nm)", fontsize = 18)
ax_WPI[1, 0].set_ylabel("Glucose Sensitivity Score", fontsize = 18)
ax_WPI[1, 0].tick_params(axis = "both", labelsize = 18)

dermal_score.plot(ax = ax_WPI[1, 1], legend = False)
ax_WPI[1, 1].set_title("Normalized Input 2: Dermal Aborption Sensitivity")
ax_WPI[1, 1].set_xlabel("Wavelength (nm)", fontsize = 18)
ax_WPI[1, 1].set_ylabel("Dermal Absorption Score", fontsize = 18)
ax_WPI[1, 1].tick_params(axis = "both", labelsize = 18)

power_score.plot(ax = ax_WPI[1, 2], legend = False)
ax_WPI[1, 2].set_title("Normalized Input 3: Detected Optical Power")
ax_WPI[1, 2].set_xlabel("Wavelength (nm)", fontsize = 18)
ax_WPI[1, 2].set_ylabel("Power Score", fontsize = 18)
ax_WPI[1, 2].tick_params(axis = "both", labelsize = 18)

water_score.plot(ax = ax_WPI[1, 3], legend = False)
ax_WPI[1, 3].set_title("Normalized Input 4: Water Fraction Sensitivity")
ax_WPI[1, 3].set_xlabel("Wavelength (nm)", fontsize = 18)
ax_WPI[1, 3].set_ylabel("Water Interference Score", fontsize = 18)
ax_WPI[1, 3].tick_params(axis = "both", labelsize = 18)

for ax in ax_WPI.flat:
    ax.grid(True)

# Calculate Wavelength Performance Index
WPI = (
    0.25 * glucose_score
    + 0.25 * dermal_score
    + 0.25 * power_score
    - 0.25 * water_score
)



# Plot WPI
WPI.plot(figsize=(10, 6))

plt.xlabel("Wavelength (nm)", fontsize = 18)
plt.ylabel("WPI", fontsize = 18)
plt.title("Wavelength Performance Index", fontsize = 18)
plt.tick_params(axis = "both", labelsize = 18)
plt.grid(True)
plt.show()