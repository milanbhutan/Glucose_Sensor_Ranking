import pmcx
import matplotlib.pyplot as plt
import numpy as np


# Format: [mu_a, mu_s, g, n]
air_optic_params= [0.01,0.01,0.01,0.01]
epidermis_optic_params = [0.01,0.01,0.01,0.01]
dermis_optic_params = [0.01,0.01,0.01,0.01]
sub_tissue_optic_params = [0.01,0.01,0.01,0.01]

#absoprtion spectra


#Define simulation values
prop=[
    air_optic_params,
    epidermis_optic_params,
    dermis_optic_params,
    sub_tissue_optic_params
]


Nx, Ny, Nz= 60, 60, 60

vol=np.zeros((Nx,Ny,Nz),dtype=np.uint8)
vol[:,:,0:3]=1
vol[:,:,3:20]=2
vol[:,:,20:]=3

sim_config={
    "nphoton": 10000,
    "vol": vol,
    "prop": prop,
    "srcpos": [30, 30, 0],
    "srcdir": [0, 0, 1],
    "detdir": [0, 0, 1, 0.01],
    "detpos": [[32, 30, 0, 1]],
    "tstart": 0,
    "tend": 1,
    "tstep": 0.0001,
    "unitinmm": 0.1,
    "autopilot": 1,
    "gpuid": 1
}

#Run simulation
sim_out=pmcx.run(sim_config)

print(sim_out.keys())
flux=sim_out["flux"];

#flo_plot=imshow(np.log10())