


# Description of cortical activity underyling rhythmic and discrete movements. 

The following code creates the first figure of the manuscript "Rhythmic and discrete arm movements arise from 
the same solution in Primary Motor Cortex but not in the Supplementary Motor Area". 
The code on this repository is for code review purposes only.

# How to run this code:

1. Download the data of the "Two-target cycling task" from Russo et al. 2020 from its [repository](https://data.mendeley.com/datasets/tfcwp8bp5j/1)
2. Download the code of this repository
3. Place the dataset in a folder separate from the code and save its path
4. Open the file `Main_reproduce_figures.m`
5. Replace the path in the variable "dataset_path" for the path where the dataset is.
6. Run `Main_reproduce_figures.m`

# Expected output

The code creates a file for each animal and cortical area recorded and store them in a folder called "Output files". Each file contains the behavioural and neural data necessary for later analyses. 

The main output of this code is the first figure of the manuscript. The figure should look like the figure below

<img width="1131" height="878" alt="image" src="https://github.com/user-attachments/assets/69e8c4ad-67e2-4ba0-bea7-92a251fd2a3d" />


  DEPENDENCIES
  ------------
 jPCA
 dPCA
 Statistics and Machine Learning Toolbox (for PCA, SUBSPACE)
