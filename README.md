# NSM-SA

Code utilizing the standard analysis (SA) for Nanofluidic Scattering Microscopy (NSM), described in [1] "B. Špačková et al.: Label-Free Nanofluidic Scattering Microscopy of Size and Mass of Single Diffusing Molecules and Nanoparticles". 

The code performs 3 main tasks:
1. Image processing including image stabilization and background subtraction - generation of a kymograph from raw data  
2. Particle tracking - finding particle trajectories within the kymograph
3. Characterization of the trajectories in terms of particles' integrated optical contrast (iOC), diffusivity (D), molecular mass (MW), and hydrodynamic radius (HR). 

## Installation
Install Matlab and download the repository from git. The code can be run directly from the supplied scripts. 

## Instructions for Demo 
To reproduce a subset of data presented in [1] we supply example data and instructions how to use to code.

**1. Sample data**

Sample data are provided in _data_ folder. The names of the files are comprised of [_ExperimentTimeStamp_]_[_Identifier_].

Indetifier:

- M - raw data
- C - processed image data
- D - particle trajectories data

The whole sample dataset can be downloaded from https://chalmers-my.sharepoint.com/:f:/g/personal/hmoberg_chalmers_se/ElpFotTfPNNInU2YghmGN3IBs5cPPgMxjsdPwtsXserlYA?e=53WB9a
Place the folder "Demo Data" in directory NSM-SA/data/.
This data corresponds to the artificial dataset used for Fig. S14.

**2. Processing the raw data**

- Change the required field (ExperimentTimeStamp) in the saveKymograph.m to specify the name of the file to be analyzed. 
- A serie of settings is described in the heading of saveKymograph.m. We recommend to use the default values        (optimized for the collected data).
- Run _saveKymograph.m_
- Results of the image processing is saved as _ExperimentTimeStamp_C.m_
- Results of the particle tracking is saved as _ExperimentTimeStamp_D.m_

**3. Plotting the kymograph**

- Load the _ExperimentTimeStamp_C.m_ and _ExperimentTimeStamp_D.m_ 
- Run _plots/plotKymograph.m_

**4. Plotting the scatter plot and histograms of iOC/MW and D/HR corresponding to the found particles**
- Load _ExperimentTimeStamp__D.m 
- Run _plots/plotParticles.m_

