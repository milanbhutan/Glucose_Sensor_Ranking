import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paths import MCX_DIR
from paths import sim_results_DIR
import tissue_model.tissue_params as params
import photon_sim.monte_carlo as mc


mean_metrics_path = MCX_DIR / sim_results_DIR / "mean_metrics"
mean_metrics_path.mkdir(parents = True, exist_ok = True)

error_metrics_path = MCX_DIR / sim_results_DIR / "error_metrics"
error_metrics_path.mkdir(parents = True, exist_ok = True)

wavelength_range = range(600, 1701, 25) # was 10
batches = 2
sds = ["SDS_1.7mm", "SDS_2.0mm", "SDS_2.3mm", "SDS_2.6mm", "SDS_2.9mm"]

# Glucose Effects
dmu_s_dCg = -0.0029

df_dmu_a_dCg_600_999 = (
    pd.read_csv(MCX_DIR / "Glucose_Absorption_Data" / "haque2021_linear_glucose_absorption_600_999nm.csv")
    .set_index("wavelength_nm")
    .sort_index()
)

df_dmu_a_dCg_1000_1700 = (
    pd.read_csv(MCX_DIR / "Glucose_Absorption_Data" / "ge2024_dmu_a_dCg.csv")
    .set_index("wavelength_nm")
    .sort_index()
)

df_dmu_a_dCg = pd.concat([df_dmu_a_dCg_600_999, df_dmu_a_dCg_1000_1700])

dmu_a_dCg = df_dmu_a_dCg.loc[:, "dmu_a_dCg_mm_inv_per_mmol_L"]






plt.figure()
plt.plot(df_dmu_a_dCg)
plt.title("$d\mu_a/dC_g$ from Haque 2021 (600 - 999nm)\nand Ge 2024 (1000 - 1700nm)", fontsize = 18)
plt.xlabel("Wavelength (nm)", fontsize = 18)
plt.xticks(fontsize = 18)
plt.ylabel("$d\mu_a/dC_g$", fontsize = 18)
plt.yticks(fontsize = 18)



batch_metrics = {
    "dA_dCg_total": [None] * batches,
    "dA_dmu_a": [None] * batches,
    "dA_dmu_s": [None] * batches,
    "dA_dCg_mu_a": [None] * batches,
    "dA_dCg_mu_s": [None] * batches,
    "mean_path_epi": [None] * batches,
    "mean_path_derm": [None] * batches,
    "mean_path_sub": [None] * batches,
    "epi_dA_dmu_s": [None] * batches,
    "derm_dA_dmu_s": [None] * batches,
    "sub_dA_dmu_s": [None] * batches,
    "beta_epi": [None] * batches,
    "beta_derm": [None] * batches,
    "beta_sub": [None] * batches,
    "power_dB": [None] * batches,
}

mean_df = dict.fromkeys(batch_metrics, None)
stderr_df = dict.fromkeys(batch_metrics, None)

for batch in range(batches):
    # Create dataframes for one batch
    dA_dCg_total = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )
    dA_dmu_a = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )
    dA_dmu_s = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )
    dA_dCg_mu_a = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )
    dA_dCg_mu_s = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )


    mean_path_epi = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )
    mean_path_derm = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )
    mean_path_sub = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )


    epi_dA_dmu_s = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )
    derm_dA_dmu_s = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )
    sub_dA_dmu_s = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )


    beta_epi = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )
    beta_derm = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )
    beta_sub = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )


    power_dB = pd.DataFrame(
        np.nan,
        index = wavelength_range,
        columns = sds
    )
    for sim_wavelength_nm in wavelength_range:
        print(f"Simulating {sim_wavelength_nm} in batch {batch}")
        #Following Ge 2024 procedure for computing dA/dCg
        baseline = mc.run_simulation(
            [
            params.glass_mu_a_sim,
            params.glass_mu_s_sim,
            params.glass_g_sim,
            params.glass_n_sim
            ],

            [
            params.epi_mu_a_sim[sim_wavelength_nm], 
            params.epi_mu_s_sim[sim_wavelength_nm], 
            params.epi_g_sim[sim_wavelength_nm], 
            params.epi_n_sim],

            [
            params.derm_mu_a_sim[sim_wavelength_nm], 
            params.derm_mu_s_sim[sim_wavelength_nm], 
            params.derm_g_sim[sim_wavelength_nm], 
            params.derm_n_sim],

            [
            params.sub_mu_a_sim[sim_wavelength_nm], 
            params.sub_mu_s_sim[sim_wavelength_nm], 
            params.sub_g_sim[sim_wavelength_nm], 
            params.sub_n_sim]
        )


        dA_dmu_a.loc[sim_wavelength_nm] = (
            baseline["mean_path"]["derm_mean_path"]
            +
            baseline["mean_path"]["sub_mean_path"]
            #(base_dmu_a["absorbance"] - baseline["absorbance"])
            )

        dA_dmu_s.loc[sim_wavelength_nm] = (
            baseline["dA_dmu_s"]["derm_dA_dmu_s"]
            +
            baseline["dA_dmu_s"]["sub_dA_dmu_s"]
            #(base_dmu_s["absorbance"] - baseline["absorbance"])
            )

        dA_dCg_mu_a.loc[sim_wavelength_nm] = (
            baseline["mean_path"]["derm_mean_path"] * dmu_a_dCg.loc[sim_wavelength_nm]
            +
            baseline["mean_path"]["sub_mean_path"] * dmu_a_dCg.loc[sim_wavelength_nm]
            #(base_dmu_a["absorbance"] - baseline["absorbance"]) * dmu_a_dCg.loc[sim_wavelength_nm]
            )

        dA_dCg_mu_s.loc[sim_wavelength_nm] = (
            baseline["dA_dmu_s"]["derm_dA_dmu_s"] * dmu_s_dCg * params.derm_mu_s_sim[sim_wavelength_nm]
            +
            baseline["dA_dmu_s"]["sub_dA_dmu_s"] * dmu_s_dCg * params.sub_mu_s_sim[sim_wavelength_nm]
            #(base_dmu_s["absorbance"] - baseline["absorbance"]) * dmu_s_dCg * params.derm_mu_s_sim[sim_wavelength_nm] 

            )

        dA_dCg_total.loc[sim_wavelength_nm] = (
            dA_dCg_mu_a.loc[sim_wavelength_nm] + dA_dCg_mu_s.loc[sim_wavelength_nm]
            )

        mean_path_epi.loc[sim_wavelength_nm] = baseline["mean_path"]["epi_mean_path"]
        mean_path_derm.loc[sim_wavelength_nm] = baseline["mean_path"]["derm_mean_path"]
        mean_path_sub.loc[sim_wavelength_nm] = baseline["mean_path"]["sub_mean_path"]

        epi_dA_dmu_s.loc[sim_wavelength_nm] = baseline["dA_dmu_s"]["epi_dA_dmu_s"]
        derm_dA_dmu_s.loc[sim_wavelength_nm] = baseline["dA_dmu_s"]["derm_dA_dmu_s"]
        sub_dA_dmu_s.loc[sim_wavelength_nm] = baseline["dA_dmu_s"]["sub_dA_dmu_s"]
        
        beta_epi.loc[sim_wavelength_nm] = baseline["sensitivity"]["epi_sensitivity"]
        beta_derm.loc[sim_wavelength_nm] = baseline["sensitivity"]["derm_sensitivity"]
        beta_sub.loc[sim_wavelength_nm] = baseline["sensitivity"]["sub_sensitivity"]

        power_dB.loc[sim_wavelength_nm] = baseline["power_dB"]
        

    #Add wavelength sweep results as a batch entry
    batch_metrics["dA_dCg_total"][batch] = dA_dCg_total
    batch_metrics["dA_dmu_a"][batch] = dA_dmu_a
    batch_metrics["dA_dmu_s"][batch] = dA_dmu_s
    batch_metrics["dA_dCg_mu_a"][batch] = dA_dCg_mu_a
    batch_metrics["dA_dCg_mu_s"][batch] = dA_dCg_mu_s
    batch_metrics["mean_path_epi"][batch] = mean_path_epi
    batch_metrics["mean_path_derm"][batch] = mean_path_derm
    batch_metrics["mean_path_sub"][batch] = mean_path_sub
    batch_metrics["epi_dA_dmu_s"][batch] = epi_dA_dmu_s
    batch_metrics["derm_dA_dmu_s"][batch] = derm_dA_dmu_s
    batch_metrics["sub_dA_dmu_s"][batch] = sub_dA_dmu_s
    batch_metrics["beta_epi"][batch] = beta_epi
    batch_metrics["beta_derm"][batch] = beta_derm
    batch_metrics["beta_sub"][batch] = beta_sub
    batch_metrics["power_dB"][batch] = power_dB

for key, value in batch_metrics.items():
    mean_df[key] = pd.DataFrame(
        data = np.average(value, axis = 0),
        index = wavelength_range,
        columns = sds
    )
    stderr_df[key] = pd.DataFrame(
        data = np.std(value, axis = 0, ddof = 1) / np.sqrt(batches),
        index = wavelength_range,
        columns = sds
    )


for metric, df in mean_df.items():
    df.to_csv(mean_metrics_path / f"Metric_Mean_{metric}.csv", index = True, index_label = "Wavelength (nm)")

for metric, df in stderr_df.items():
    df.to_csv(error_metrics_path / f"Metric_Error_{metric}.csv", index = True, index_label = "Wavelength (nm)")


plt.show()