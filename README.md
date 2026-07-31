# Motor cortical dynamics driving rhythmic and discrete movements. 

The following code creates the figures of the manuscript "Primary and Supplementary Motor cortex implement parallel
solutions for the control of rhythmic and discrete arm movements". 

> Primary and Supplementary Motor cortex implement parallel solutions for the control of rhythmic and discrete arm movements
> 
> Andrea Colins Rodriguez*, Romulo Fuentes, Mark D. Humphries*
> 
>  [DOI: https://doi.org/10.64898/2026.02.16.706176](https://doi.org/10.64898/2026.02.16.706176)

# Hardware requirements
This package requires only a standard computer. Examples in Jupyter Notebooks can be executed in Google Colab. 

# Software requirements

- Matlab 2023 or 2024: To reproduce the figures of the paper.
- Python (and Tensorflow): To train RNNs.
  
# Data
All data used for the analyses were published by Russo et al. (2020) and can be accessed from their [repository](https://data.mendeley.com/datasets/tfcwp8bp5j/1)

# How to run this code:
## Dependencies
This code performs jPCA and dPCA on Matlab. The toolboxes required for these analyses can be downloaded here:
-  dPCA: [dPCA toolbox](https://github.com/machenslab/dPCA)
-  jPCA: [dPCA toolbox](https://churchland.zuckermaninstitute.columbia.edu/content/code)

## Reproducing figures of the papers
1. Download the data of the "Two-target cycling task" from Russo et al. 2020 from its [repository](https://data.mendeley.com/datasets/tfcwp8bp5j/1)
2. Download the code of this repository
3. Place the dataset in a folder separate from the code and save its path
4. Download the jPCA and dPCA toolboxes 
5. Open the file `Main_reproduce_figures.m` in Matlab
6. Replace the path in the variable `dataset_path` with the path where the dataset is.
7. Replace the paths in the variables  `jPCA_path` and `dPCA_path` with the paths where the dPCA toolbox and this code are
8. Run `Main_reproduce_figures.m`

   Note: By default, the code performs the analyses to reproduce the three figures of the paper. To add the supplementary figures and videos, please set the value of the variable `plot_supp_figs.do_plot` as 1. 

# Expected output

The code creates a file for each animal and cortical area recorded and stores them in a folder called "Output files". Each file contains the behavioural and neural data necessary for later analyses. 

This script generates:

- 3 figures corresponding to the main figures of the text.

- 7 supplementary figures that were used to compute complementary statistics.
  
- 2 supplementary videos

# License
This project is covered under the MIT License.

