# Motor cortical dynamics driving rhythmic and discrete movements. 

The following code creates the figures of the manuscript "Primary and Supplementary Motor cortex implement parallel
solutions for the control of rhythmic and discrete arm movements". 

> Primary and Supplementary Motor cortex implement parallel solutions for the control of rhythmic and discrete arm movements
> 
> Andrea Colins Rodriguez, Romulo Fuentes, Mark D. Humphries
> 
>  [DOI: XX](XX)



# How to run this code:
## Dependencies
This code performs jPCA and dPCA. The toolboxes required for these analyses can be downloaded here:
-  dPCA: [dPCA toolbox](https://github.com/machenslab/dPCA)
-  jPCA: [dPCA toolbox](https://churchland.zuckermaninstitute.columbia.edu/content/code)

## Running the analysis
1. Download the data of the "Two-target cycling task" from Russo et al. 2020 from its [repository](https://data.mendeley.com/datasets/tfcwp8bp5j/1)
2. Download the code of this repository
3. Place the dataset in a folder separate from the code and save its path
4. Download the jPCA and dPCA toolboxes 
5. Open the file `Main_reproduce_figures.m`
6. Replace the path in the variable `dataset_path` with the path where the dataset is.
7. Replace the paths in the variables  `jPCA_path` and `dPCA_path` with the paths where the dPCA toolbox and this code are
8. Run `Main_reproduce_figures.m`

# Expected output

The code creates a file for each animal and cortical area recorded and stores them in a folder called "Output files". Each file contains the behavioural and neural data necessary for later analyses. 

This script generates:

- 3 figures corresponding to the main figures of the text.

- 7 supplementary figures that were used to compute complementary statistics.



