# AugmentedDesign

## Background
In an augmented experimental design, a set of checks are replicated across blocks, while the test genotypes are evaluated only once or a few times. 
This is common in plant breeding when limited experimental resources are available for field trials.
The main challenge is to separate variation due to block and genotype effects and then estimate heritability.

## Model
A possible model for analyzing an augmented design is 

$y_{ij}=\mu+B_j+G_i+\varepsilon_{ij}$

where 
$y_{ij}$: observed phenotype of genotype $i$ in block $j$.  
$\mu$: overall mean.  
$B_j$: block effect  
$G_i$: genetic effect of genotype $i$  
$\varepsilon_{ij}$: residual error.  

This model assumes that both block and genotype effects are random that is  
$B_j\sim N(0,\sigma_B^2)$  
$G_i\sim N(0,\sigma_G^2)$  
$\varepsilon_{ij}\sim N(0,\sigma_e^2)$  


## Results
### Estimation of broad sense heritability
| Trait          | N_Genos | Mean_BLUE | SD_BLUE | CV_BLUE_pct | Sigma2_G   | Sigma2_E    | H2_Aug | H2_Cullis | H2_AugmentedRCBD | H2_sommer |
|----------------|---------|-----------|---------|-------------|------------|-------------|--------|-----------|------------------|-----------|
| TotCH          | 179     | 3         | 10.46   | 348.74      | 16.1623    | 80.4727     | 0.167  | 0.184     | NA               | 0.167     |
| Fv_Fm          | 179     | 0         | 0.06    | 1425.35     | 0          | 2.00E-04    | 0.149  | 0.164     | 0.499            | 0.15      |
| Sm             | 179     | 2.12      | 82.27   | 3888.6      | 659.0161   | 4698.1263   | 0.123  | 0.139     | 0.185            | 0.123     |
| Pi_Abs         | 179     | 0.02      | 1.97    | 8056.2      | 0.3366     | 3.2483      | 0.094  | 0.108     | 0.255            | 0.094     |
| plant_weight_g | 178     | 4.53      | 251.62  | 5559.89     | 14579.1013 | 46855.7174  | 0.237  | 0.253     | 0.243            | 0.237     |
| fruit_weight_g | 177     | 7.11      | 411.36  | 5781.66     | 38987.6187 | 120560.5735 | 0.244  | 0.258     | 0.391            | 0.244     |
