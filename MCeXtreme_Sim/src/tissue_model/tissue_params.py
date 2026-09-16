derm_params = [
    "blood",
    "oxygen_saturation",
    "water",
    "gamma"
]

sub_params = [
    "blood",
    "oxygen_saturation",
    "water",
    "fat",
    "gamma",
]

derm_comp_upper = {
    "blood" : 0.015,
    "oxygen_saturation" : 0.65,    
    "water" : 0.75,
    "gamma" : 1 # original: 0.8

}

derm_comp_lower = {
    "blood" : 0.002,
    "oxygen_saturation" : 0.39,    
    "water" : 0.65,
    "gamma" : 0.1 # original: 0.2

}

sub_comp_upper = {
    "blood" : 0.010,
    "oxygen_saturation" : 0.80,    
    "water" : 0.30,
    "fat" : 0.80,
    "gamma" : 1, # original: 0.3

}

sub_comp_lower = {
    "blood" : 0.002,
    "oxygen_saturation" : 0.50,    
    "water" : 0.10,
    "fat" : 0.3, #original: 0.55
    "gamma" : 0.1, #original: 0.2

}