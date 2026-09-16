# import pmcx


# # Comment out to plot tissue spectras
# # plt.close("all")

# # =========================================
# # MCeXtreme Simulation
# # =========================================

# sim_wavelength_nm = 600

# # Format: [mu_a, mu_s, g, n]
# air_optic_params= [0.01, 0.01, 1, 1]
# epi_optic_params = [epi_mu_a_sim[sim_wavelength_nm], epi_mu_s_sim[sim_wavelength_nm], epi_g, epi_n]
# derm_optic_params = [derm_mu_a_sim[sim_wavelength_nm], derm_mu_s_sim[sim_wavelength_nm], derm_g, derm_n]
# sub_optic_params = [sub_mu_a_sim[sim_wavelength_nm], sub_mu_s_sim[sim_wavelength_nm], sub_g, sub_n]


# # Define simulation values
# prop = [
#     air_optic_params,
#     epi_optic_params,
#     derm_optic_params,
#     sub_optic_params
# ]


# Nx, Ny, Nz = 160, 160, 160

# vol = np.zeros((Nx, Ny, Nz), dtype = np.uint8)
# vol[:,:,0:3] = 1
# vol[:,:,3:13] = 2
# vol[:,:,13:160] = 3


# sim_config = {
#     "nphoton": 1e7, # 1e9: 57 seconds
#     "vol": vol,
#     "prop": prop,
#     "srcpos": [40, 30, 0],
#     "srcdir": [0, 0, 1],
#     "detdir": [0, 0, 1, 20],
#     "detpos": [[120, 30, 0, 20]],
#     "tstart": 0,
#     "tend": 10e-9,
#     "tstep": 1e-9,
#     "unitinmm": 0.1,
#     "autopilot": 1,
#     "gpuid": 1,
#     "issrcfrom0": 1,
#     "issavedet": 1,
#     "issaveseed": 1,
#     "savedetflag" : "dp"
# }

# # Run simulation


# sim_out = pmcx.run(sim_config)

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
# axes_flux[0].axhline(3, linestyle = "--", color = "black")
# axes_flux[0].axhline(13, linestyle = "--", color = "black")
# axes_flux[0].text(1, 1.5, "epidermis")
# axes_flux[0].text(1, 8, "dermis")
# axes_flux[0].text(1, 30, "subcutaneous")
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
