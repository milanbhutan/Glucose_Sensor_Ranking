import pandas as pd
import matplotlib.pyplot as plt
import tissue_model.layer_calc as calc
import tissue_model.tissue_params as params
import numpy as np
from paths import MCX_DIR
from paths import sim_results_DIR

P_LED = 5e-3          # W
responsivity = 0.75   # A/W
I_dark = 50e-9        # A
R = 10_000            # ohms
T = 298.15            
bandwidth = 10        # Hz

q = 1.602176634e-19
kB = 1.380649e-23

mean_path_epi_df = pd.read_csv(MCX_DIR / sim_results_DIR / "mean_metrics"/ "Metric_Mean_mean_path_epi.csv", index_col = "Wavelength (nm)")
mean_path_derm_df = pd.read_csv(MCX_DIR / sim_results_DIR / "mean_metrics"/ "Metric_Mean_mean_path_derm.csv", index_col = "Wavelength (nm)")
mean_path_sub_df = pd.read_csv(MCX_DIR / sim_results_DIR / "mean_metrics"/ "Metric_Mean_mean_path_sub.csv", index_col = "Wavelength (nm)")

derm_dA_dmu_s = pd.read_csv(MCX_DIR / sim_results_DIR / "mean_metrics"/ "Metric_Mean_derm_dA_dmu_s.csv", index_col = "Wavelength (nm)")


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

# chain rule: dA/dW = dA/dmu_a * dmu_a/dW
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

# Calculate Blood Interference per layer, following Ge 2025
# 50% increase in blood volume so multiply original blood fraction by 1.5
B_multiplier = 1.50

derm_dmu_a_dB = (
    - params.derm_mu_a_sim 
    + calc.derm(
        params.derm_Jacques["blood"] * B_multiplier, 
        params.derm_Jacques["oxygen_saturation"], 
        params.derm_Jacques["water"], 
        params.derm_Jacques["fat"], 
        params.derm_Jacques["melanosome"],
        params.baseline_derm,
    )[0]) / (0.5 * params.derm_Jacques["blood"])




# From Ge 2025: dmu_s/dB = mu_s / B
derm_dmu_s_dB = params.derm_mu_s_sim / params.derm_Jacques["blood"]

plt.plot(params.derm_mu_s_sim)

# chain rule: dA/dB = (dA/dmu_a * dmu_a/dB) + (dA/dmu_s * dmu_s/dB)
derm_dmu_a_dB = derm_dmu_a_dB.reindex(
    mean_path_derm_df.index
)


derm_dmu_s_dB = derm_dmu_s_dB.reindex(
    derm_dA_dmu_s.index
)

derm_dA_dB_dmu_a = mean_path_derm_df.mul(
    derm_dmu_a_dB,
    axis = 0)

derm_dA_dB_dmu_s = derm_dA_dmu_s.mul(
    derm_dmu_s_dB,
    axis = 0)

dA_dB_df = derm_dA_dB_dmu_a + derm_dA_dB_dmu_s

P_det = P_LED * 10**(power_dB_df / 10)
I_det = responsivity * P_det

i_thermal = np.sqrt(4 * kB * T * bandwidth / R)
i_dark = np.sqrt(2 * q * I_dark * bandwidth)
i_shot = np.sqrt(2 * q * I_det * bandwidth)

i_noise = np.sqrt(i_thermal**2 + i_dark**2 + i_shot**2)

delta_Cg = 1.0  # mmol/L

S_g = dA_dCg_df

#i_signal = I_det * np.expm1(-S_g * delta_Cg)


SNR_G = abs(I_det) / i_noise



# ============== PLOTTING ==============



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
blood_score = normalize(abs(dA_dB_df))

fig_WPI, ax_WPI = plt.subplots(2, 5, layout = "constrained")
fig_WPI.suptitle("WPI Inputs", fontsize = 18)

dA_dCg_df.plot(ax = ax_WPI[0, 0])
ax_WPI[0 ,0].set_title("Raw Input 1: Glucose Sensitivity")
ax_WPI[0, 0].set_xlabel("")
ax_WPI[0, 0].set_ylabel("$dA/dC_g$", fontsize = 18)
ax_WPI[0, 0].tick_params(axis = "both", labelsize = 18)
ax_WPI[0, 0].tick_params(axis = "x", labelbottom = False)

beta_derm_df.plot(ax = ax_WPI[0, 1], legend = False)
ax_WPI[0, 1].set_title("Raw Input 2: Dermal Absorption Sensitivity")
ax_WPI[0, 1].set_xlabel("")
ax_WPI[0, 1].set_ylabel("Dermal Path Fraction (%)", fontsize = 18)
ax_WPI[0, 1].tick_params(axis = "both", labelsize = 18)
ax_WPI[0, 1].tick_params(axis = "x", labelbottom = False)

power_dB_df.plot(ax = ax_WPI[0, 2], legend = False)
ax_WPI[0, 2].set_title("Raw Input 3: Detected Optical Power")
ax_WPI[0, 2].set_xlabel("")
ax_WPI[0, 2].set_ylabel("Power (dB)", fontsize = 18)
ax_WPI[0, 2].tick_params(axis = "both", labelsize = 18)
ax_WPI[0, 2].tick_params(axis = "x", labelbottom = False)

dA_dW_df.plot(ax = ax_WPI[0, 3], legend = False)
ax_WPI[0, 3].set_title("Raw Input 4: Water Interference")
ax_WPI[0, 3].set_xlabel("")
ax_WPI[0, 3].set_ylabel("$dA/dW$", fontsize = 18)
ax_WPI[0, 3].tick_params(axis = "x", labelbottom = False)

dA_dB_df.plot(ax = ax_WPI[0, 4], legend = False)
ax_WPI[0, 4].set_title("Raw Input 5: Blood Interference")
ax_WPI[0, 4].set_xlabel("")
ax_WPI[0, 4].set_ylabel("$dA/dB$", fontsize = 18)
ax_WPI[0, 4].tick_params(axis = "x", labelbottom = False)


glucose_score.plot(ax = ax_WPI[1, 0], legend = False)
ax_WPI[1, 0].set_title("Normalized Input 1:\nGlucose Sensitivity")
ax_WPI[1, 0].set_xlabel("Wavelength (nm)", fontsize = 18)
ax_WPI[1, 0].set_ylabel("Glucose Sensitivity Score", fontsize = 18)
ax_WPI[1, 0].tick_params(axis = "both", labelsize = 18)

dermal_score.plot(ax = ax_WPI[1, 1], legend = False)
ax_WPI[1, 1].set_title("Normalized Input 2:\nDermal Aborption Sensitivity")
ax_WPI[1, 1].set_xlabel("Wavelength (nm)", fontsize = 18)
ax_WPI[1, 1].set_ylabel("Dermal Absorption Score", fontsize = 18)
ax_WPI[1, 1].tick_params(axis = "both", labelsize = 18)

power_score.plot(ax = ax_WPI[1, 2], legend = False)
ax_WPI[1, 2].set_title("Normalized Input 3:\nDetected Optical Power")
ax_WPI[1, 2].set_xlabel("Wavelength (nm)", fontsize = 18)
ax_WPI[1, 2].set_ylabel("Power Score", fontsize = 18)
ax_WPI[1, 2].tick_params(axis = "both", labelsize = 18)

water_score.plot(ax = ax_WPI[1, 3], legend = False)
ax_WPI[1, 3].set_title("Normalized Input 4:\nWater Interference")
ax_WPI[1, 3].set_xlabel("Wavelength (nm)", fontsize = 18)
ax_WPI[1, 3].set_ylabel("Water Interference Score", fontsize = 18)
ax_WPI[1, 3].tick_params(axis = "both", labelsize = 18)

blood_score.plot(ax = ax_WPI[1, 4], legend = False)
ax_WPI[1, 4].set_title("Normalized Input 5:\nBlood Interference")
ax_WPI[1, 4].set_xlabel("Wavelength (nm)", fontsize = 18)
ax_WPI[1, 4].set_ylabel("Blood Interference Score", fontsize = 18)
ax_WPI[1, 4].tick_params(axis = "both", labelsize = 18)

dA_dB_df.plot(legend = False)
plt.title("Raw Input 5: Blood Interference")
plt.xlabel("Wavelength (nm)", fontsize = 18)
plt.ylabel("$dA/dB$", fontsize = 18)
plt.grid(visible = True)

SNR_G.plot()
plt.title("Baseline Signal to Noise Ratio")
plt.xlabel("Wavelength (nm)", fontsize = 18)
plt.ylabel("$SNR$", fontsize = 18)
plt.grid(visible = True)


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