# Master Thesis MATLAB Programs

This repository contains two MATLAB programs developed as part of my master's thesis. 
These programs are designed to process experimental images for specific analyses. 
Below you will find instructions on how to use each script, along with the required data structure.

## Programs

### 1. Cropping Script

The `Cropping_script.m` script is designed for image cropping and preprocessing. It utilizes a `Greenmask.m` function to isolate areas of interest in the images.

#### Usage

- Run the `Cropping_script.m` script directly in MATLAB.
- Before processing each experiment, update the thresholding values in the `Greenmask.m` function. This can be done using the ColorThresholder app in Matlab.
- The `Greenmask.m` function is stored as a separate file.

### 2. StrainMeasuring Script

The `StrainMeasuring` script is used to process images to find strains on the test specimens.
It is structured in blocks, with the first block dedicated to loading and processing images, and subsequent blocks focusing on smoothing and post-processing.
The use of the subsequent block is optional. For all data presented in the final report, the post-processing was performed using Mcalibration.

#### Usage

- Run the `StrainMeasuring` script directly in MATLAB.
- The script uses an inline `Greenmask.m` function for initial processing. Before processing each experiment, update the thresholding values in the `Greenmask.m` function.

## Data Structure

For both programs, data should be organized in a specific manner:

- Create a master folder for all your experiments.
- Inside the master folder, create separate folders for each experiment. Name this folder how you like.
- Within each experiment's folder, organize the images into a folder named `Photos`. Ensure that when sorted alphabetically the order is correct.
- In the experiment folder there must be a text file name named `specimen_information.tex` in each experiment folder, containing the dimensions of the test specimen.

## Requirements

- MATLAB (Version recommended: 2023b).
# Toolboxes
- Image processing toolbox
- Curve Fitting Toolbox
- Computer Vision Toolbox

# Matlab Add-ons
The program uses a Matlab addon to display the progress of the programs. 
downloaded here: https://se.mathworks.com/matlabcentral/fileexchange/56871-command-line-progress-bar-waitbar




## Contact

For questions or further information, please contact Asger Riis Vienberg at avienberg@hotmail.com
