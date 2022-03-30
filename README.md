# NSM-SA

Code utilizing the standard analysis (SA) for Nanofluidic Scattering Microscopy (NSM), described in "B. Špačková et al.: Label-Free Nanofluidic Scattering Microscopy of Size and Mass of Single Diffusing Molecules and Nanoparticles". 

The code performs 3 main tasks:
1. Image processing including image stabilization and background subtraction - generation of a kymograph from raw data  
2. Particle tracking - finding particle trajectories within the kymograph
3. Characterization of the trajectories in terms of particles' integrated optical contrast (iOC), diffusivity (D), molecular mass (MW), and hydrodynamic radius (HR). 

## Installation
Install Matlab and download the repository from git. The code can be run directly from the supplied scripts. 

## Instructions for Demo 

1. Sample data

Sample data are provided in data folder. The names of the files are comprised of ExperimentTimeStamp + Idenfier.

ExperimentTimeStamp:

- Ferritin in Channel I: 13-10-20_14-01-50
tbc...

Indetifier:

- M - raw data
- C - processed image data
- D - particle trajectories data

2. Process the raw data

- Change the required field (ExperimentTimeStamp) in the saveKymograph.m to specify the name of the file to be analyzed. 
- A serie of settings is described in the heading of saveKymograph.m. We recommend to use the default values (optimized for the collected data).
- Run saveKymograph.m
- Results of the image processing is saved as ExperimentTimeStamp_C.m
- Results of the particle tracking is saved as ExperimentTimeStamp_D.m

3. Plot the kymograph

- Load the ExperimentTimeStamp_C.m and ExperimentTimeStamp_D.m 
- Run plots/plotKymograph.m

4. Plot the scatter plot and histograms of iOC/MW and D/HR corresponding to the found particles
- Load ExperimentTimeStamp_D.m 
- Run plots/plotParticles.m

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

