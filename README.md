# Description of cortical activity underyling rhythmic and discrete movements. 

The following code creates the first figure of the manuscript "Rhythmic and discrete arm movements arise from 
the same solution in Primary Motor Cortex but not in the Supplementary Motor Area". 
The code on this repository is for code review purposes only.

# How to run this code:
## Dependencies
This code performs jPCA and dPCA. The toolboxes requiered for these analyses can be downladed here:
-  dPCA: [dPCA toolbox](https://github.com/machenslab/dPCA)
-  jPCA: [dPCA toolbox](https://churchland.zuckermaninstitute.columbia.edu/content/code)

1. Download the data of the "Two-target cycling task" from Russo et al. 2020 from its [repository](https://data.mendeley.com/datasets/tfcwp8bp5j/1)
2. Download the code of this repository
3. Place the dataset in a folder separate from the code and save its path
4. Download the jPCA and dPCA toolboxes 
5. Open the file `Main_reproduce_figures.m`
6. Replace the path in the variable `dataset_path` for the path where the dataset is.
7. Replace the paths in the variables  `jPCA_path` and `dPCA_path` for the paths where the dPCA toolbox and this code are
8. Run `Main_reproduce_figures.m`

# Expected output

The code creates a file for each animal and cortical area recorded and store them in a folder called "Output files". Each file contains the behavioural and neural data necessary for later analyses. 

This script generates:

- 3 figures corresponding to the main figures of the text.

- 7 supplementary figures that were used to compute complementary statistics.



