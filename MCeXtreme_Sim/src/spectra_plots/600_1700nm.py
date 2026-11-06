"""Plot the adopted three-layer optical model, marking modeled scattering tails."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from paths import MCX_DIR
import tissue_model.tissue_params as params


wavelength_min, wavelength_max = 600, 1700
absorption_spectra = [
    params.epi_mu_a_sim,
    params.derm_mu_a_sim,
    params.sub_mu_a_sim,
]
layer_info = [
    ("Epidermis", "epidermis", 1598),
    ("Dermis", "dermis", 1584),
    ("Subcutaneous", "subcutaneous_fat", 1585),
]

fig, axes = plt.subplots(
    2, 3, figsize=(16, 9), layout="constrained", sharex="col"
)
fig.suptitle("Optical Properties of Three-Layer Tissue Model", fontsize=20)

for i, (layer, slug, last_source_nm) in enumerate(layer_info):
    # Top row: your existing constituent-based absorption model.
    mu_a = absorption_spectra[i].loc[wavelength_min:wavelength_max]
    axes[0, i].plot(
        mu_a.index, mu_a.values,
        color="red", linewidth=1.8, label="Constituent-Based Model"
    )
    axes[0, i].set_title(layer, fontsize=18)
    axes[0, i].set_ylabel(r"Absorption $\mu_a$ (mm$^{-1}$)", fontsize=14)

    # Bottom row: the revised Salomatina CSV, split at its actual data limit.
    filename = f"salomatina_2006_{slug}_reduced_scattering_extended_1700nm.csv"
    data = pd.read_csv(
        MCX_DIR / "Exp_Scattering" / filename,
        comment="#",
    ).set_index("wavelength_nm").sort_index()
    reduced_scat = data["mu_s_prime_mm^-1"].loc[wavelength_min:wavelength_max]

    measured = reduced_scat.loc[:last_source_nm]
    modeled = reduced_scat.loc[last_source_nm:wavelength_max]

    axes[1, i].plot(
        measured.index, measured.values,
        color="black", linewidth=1.8, label="Salomatina 2006 (Source)"
    )
    axes[1, i].plot(
        modeled.index, modeled.values,
        color="black", linewidth=2.0, linestyle="--",
        label="Extrapolated (Model)"
    )
    axes[1, i].set_ylabel(
        r"Reduced Scattering $\mu_s'$ (mm$^{-1}$)", fontsize=14
    )
    axes[1, i].set_xlabel("Wavelength (nm)", fontsize=14)

    for ax in axes[:, i]:
        ax.set_xlim(wavelength_min, wavelength_max)
        ax.grid(True, alpha=0.25)
        ax.tick_params(axis="both", labelsize=12)
        ax.legend(fontsize=9)

for ax in axes[0, :]:
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelbottom=False)

plt.show()
