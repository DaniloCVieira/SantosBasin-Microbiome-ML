# Microbial Composition and Environmental Properties in the Santos Basin: Data and Code Repository

This repository contains the resources and data associated with the manuscript **"A Machine Learning Approach Elucidates Spatial Patterns of Environmental Properties Driving Microbial Composition Over Santos Basin, South Atlantic."** 

The study investigates the structure and diversity of microbial communities in the Santos Basin (SB), Brazil's largest marine sedimentary basin and an ecologically and industrially significant region. By integrating amplicon sequencing data and flow cytometry quantitative cell counts with environmental parameters, we provide insights into microbial diversity and the environmental gradients shaping these communities.


 ├── family_data.csv        # Microbial family abundance data
    ├── Cito_data.csv          # Cytometry data for microbial cell counts
    ├── Coords_Depth_data.csv  # Coordinates and depth information
    ├── environ_data.csv       # Environmental variable dataset
    ├── depth_layers.csv       # Depth layer metadata
    ├── args_map3D.rds         # Arguments for 3D maps visualization
    ├── 01-Script_ML_supplementary.R       # Main machine learning analysis script
    ├── 02-Script_3D_maps_supplementary.R  # Script for generating 3D visualizations
    ├── Auxiliar_functions.R               # Functions for data analysis
