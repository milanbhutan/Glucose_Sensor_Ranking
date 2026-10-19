import pmcx
import numpy as np
import sys
from pathlib import Path
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parents[1]))

import tissue_model.tissue_params as params


air_optic_params = [
    0.01, 
    0.01, 
    1, 
    1.4]

Nx, Ny, Nz = params.vox_x, params.vox_y, params.vox_z

vol = np.zeros((Nx, Ny, Nz), dtype = np.uint8)
vol[:,
    :,
    0 : params.bot_epi_depth_vox
    ] = 1

vol[:,
    :,
    params.bot_epi_depth_vox : params.bot_derm_depth_vox] = 2

vol[:,
    :,
    params.bot_derm_depth_vox : params.bot_sub_depth_vox] = 3


def run_simulation(glass_optic_params, epi_optic_params, derm_optic_params, sub_optic_params):    
    # =========================================
    # MCeXtreme Simulation
    # =========================================

    # Format: [mu_a, mu_s, g, n]


    # Define simulation values
    prop = [
        air_optic_params,
        glass_optic_params,
        epi_optic_params,
        derm_optic_params,
        sub_optic_params
    ]


    sim_config = {
        "nphoton": 1e8, # 1e9: 57 seconds
        "vol": vol,
        "prop": prop,
        "srcpos": [params.src_x_vox, params.src_y_vox, params.src_z_vox],
        "srcdir": [0, 0, 1],
        "detpos": [[params.det1_x_vox, params.det1_y_vox, params.det1_z_vox, params.det1_r_vox], 
                   [params.det2_x_vox, params.det2_y_vox, params.det2_z_vox, params.det2_r_vox],
                   [params.det3_x_vox, params.det3_y_vox, params.det3_z_vox, params.det3_r_vox],
                   [params.det4_x_vox, params.det4_y_vox, params.det4_z_vox, params.det4_r_vox],
                   [params.det5_x_vox, params.det5_y_vox, params.det5_z_vox, params.det5_r_vox]],
        "tstart": 0,
        "tend": 10e-9,
        "tstep": 10e-9,
        "unitinmm": params.vox_size_mm,
        "autopilot": 1,
        "gpuid": 1,
        "issrcfrom0": 1,
        "issavedet": 1,
        "issaveseed": 0,
        "issave2pt": 0,
        "maxdetphoton": 1e7,
        "savedetflag" : "dp"
    }


    # Run simulation
    sim_out = pmcx.mcxlab(sim_config)

    detp = sim_out["detp"]
    photon_weights = pmcx.detweight(detp)

    mask_1 = detp["detid"] == 1
    mask_2 = detp["detid"] == 2
    mask_3 = detp["detid"] == 3
    mask_4 = detp["detid"] == 4
    mask_5 = detp["detid"] == 5

    R_1 = photon_weights[mask_1].sum() / sim_config["nphoton"]
    R_2 = photon_weights[mask_2].sum() / sim_config["nphoton"]
    R_3 = photon_weights[mask_3].sum() / sim_config["nphoton"]
    R_4 = photon_weights[mask_4].sum() / sim_config["nphoton"]
    R_5 = photon_weights[mask_5].sum() / sim_config["nphoton"]

    A_1 = -np.log(R_1)
    A_2 = -np.log(R_2)
    A_3 = -np.log(R_3)
    A_4 = -np.log(R_4)
    A_5 = -np.log(R_5)

    # replay_config = sim_config.copy()
    # replay_config["outputtype"] = "jacobian"
    # replay_config["detphotons"] = sim_out["detp"]
    # replay_config["seed"] = sim_out["seeds"]

    # replay_out = pmcx.run(replay_config)

    # jacobian_PMDF = replay_out["flux"]



    # flux = sim_out["flux"]

    # fig_flux, axes_flux = plt.subplots(1, 2, layout = "constrained")

    # # select one time gate
    # flux_1 = flux[:, :, :, 0]
    # flux_1_yslice = flux_1[:, 30, :]
    # axes_flux[0].set_title("Fluence Map")
    # axes_flux[0].set_xlabel("x-axis (voxels)")
    # axes_flux[0].set_ylabel("z-axis (voxels)")
    # axes_flux[0].axhline(params.bot_epi_depth_vox, linestyle = "--", color = "black")
    # axes_flux[0].axhline(params.bot_derm_depth_vox, linestyle = "--", color = "black")
    # axes_flux[0].text(1, params.bot_epi_depth_vox / 2, "epidermis")
    # axes_flux[0].text(1, params.bot_derm_depth_vox / 2, "dermis")
    # axes_flux[0].text(1, params.bot_sub_depth_vox / 2, "subcutaneous")
    # axes_flux[0].arrow(sim_config["srcpos"][0], sim_config["srcpos"][2] - 1, 0, 2, head_width = 1, head_length = 1, color = "red", label = "source")
    # axes_flux[0].arrow(sim_config["detpos"][0][0], sim_config["detpos"][0][2] - 1, 0, 2, head_width = 1, head_length = 1, color = "blue", label = "detector" )

    # axes_flux[0].legend()
    # axes_flux[0].imshow(np.log10(np.abs(flux_1_yslice.T) + 1e-20))

    # jacobian_1 = np.sum(jacobian_PMDF, axis = 3)
    # jacobian_1_slice = jacobian_1[:, 30, :]
    # axes_flux[1].set_title("Spatial Absorption Sensitivity")
    # axes_flux[1].set_xlabel("x-axis (voxels)")
    # axes_flux[1].set_ylabel("z-axis (voxels)")
    # axes_flux[1].axhline(3, linestyle = "--", color = "black")
    # axes_flux[1].axhline(13, linestyle = "--", color = "black")
    # axes_flux[1].text(1, 1.5, "epidermis")
    # axes_flux[1].text(1, 8, "dermis")
    # axes_flux[1].text(1, 30, "subcutaneous")
    # axes_flux[1].arrow(sim_config["srcpos"][0], sim_config["srcpos"][2] - 1, 0, 2, head_width = 1, head_length = 1, color = "red", label = "source")
    # axes_flux[1].arrow(sim_config["detpos"][0][0], sim_config["detpos"][0][2] - 1, 0, 2, head_width = 1, head_length = 1, color = "blue", label = "detector" )

    # axes_flux[1].legend()
    # image = axes_flux[1].imshow(np.log10(np.abs(jacobian_1_slice.T) + 1e-20))

    # # fig_flux.colorbar(
    # #     image,
    # #     ax = axes_flux[0],
    # # )

    # fig_flux.colorbar(
    #     image,
    #     ax = axes_flux[1],
    # )

    # PMDF = np.sum(jacobian_PMDF, axis = 3)
    # PMDF_epi = np.sum(PMDF[vol == 1])
    # PMDF_derm = np.sum(PMDF[vol == 2])
    # PMDF_sub = np.sum(PMDF[vol == 3])

    # PMDF_total = PMDF_epi + PMDF_derm + PMDF_sub

    # PMDF_epi_per = 100 * (PMDF_epi / PMDF_total)
    # PMDF_derm_per = 100 * (PMDF_derm / PMDF_total)
    # PMDF_sub_per = 100 * (PMDF_sub / PMDF_total)

    # print(f"ES(%): {PMDF_epi_per:.2f} %")
    # print(f"DS(%): {PMDF_derm_per:.2f} %")
    # print(f"SS(%): {PMDF_sub_per:.2f} %")

    A_array = np.array([A_1, A_2, A_3, A_4, A_5])
    return A_array