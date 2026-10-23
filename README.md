# AugmentedDesign

## Background
In an augmented experimental design, a set of checks are replicated across blocks, while the test genotypes are evaluated only once or a few times. 
This is common in plant breeding when limited experimental resources are available for field trials.
The main challenge is to separate variation due to block and genotype effects and then estimate heritability.

## Model
Models for analyzing an experiment augmented design can might described as follows: 

$y_{ij}=\mu+B_j+G_i+\varepsilon_{ij}$ **(Model 1)**

where 
$y_{ij}$: observed phenotype of genotype $i$ in block $j$.  
$\mu$: overall mean.  
$B_j$: block effect  
$G_i$: genetic effect of genotype $i$  
$\varepsilon_{ij}$: residual error.  

This model assumes that both block and genotype effects are random that is:  
$B_j\sim N(0,\sigma_B^2)$  
$G_i\sim N(0,\sigma_G^2)$  
$\varepsilon_{ij}\sim N(0,\sigma_e^2)$  

The aforementioned model was used to 1) estimate broad sense heritability and 2) compute the adjusted means of genotypes in a panel of hybrids tested across different site-by-season combinations.
## Results


### Estimation of broad sense heritability 
Variance components of **Model 1** fitted using lme4, sommer and augmentedRCBD R packages were used to estimate broad sense heritability **(Table 1)**.  

**Table 1**: Estimation of broad sense heritability in hybrids tested in 2026  
| Trait          | *H*<sup>2</sup> Augmented | *H*<sup>2</sup> Cullis | *H*<sup>2</sup> AugmentedRCBD | *H*<sup>2</sup> sommer |
|----------------|-------------------------|----------------------|-----------------------------|----------------------|
| TotCH          | 0.167      | 0.184         | 0.137                | 0.167         |
| Fv_Fm          | 0.149      | 0.164         | 0.509                | 0.15          |
| Sm             | 0.123      | 0.139         | 0.193                | 0.123         |
| Pi_Abs         | 0.094      | 0.108         | 0.293                | 0.094         |
| plant_weight_g | 0.237      | 0.253         | 0.268                | 0.237         |
| fruit_weight_g | 0.244      | 0.258         | 0.395                | 0.244         |

*H*<sup>2</sup> Augmented : broad sense heritability estimated from **Model 1** fitted with **lme4** package using standard formula.  
*H*<sup>2</sup> Cullis: broad sense heritability estimated from **Model 1** fitted with **lme4** package using Cullis' formula.  
*H*<sup>2</sup> AugmentedRCBD: broad sense heritability estimated from **Model 1** fitted with **AugmentedRCBD** package (see manual for more info regarding the formula used for computing this parameter).  
*H*<sup>2</sup> sommer: broad sense heritability estimated from **Model 1** fitted with **sommer** package using standard formula.  

To double check the results reported in **Table 1**, *H*<sup>2</sup> was also computed using a specific R package (heritable), whch allows also the computation of different estimates of *H*<sup>2</sup> (**Table 2**)

**Table 2**: Estimation of broad sense heritability in hybrids tested in 2026 using heritable package.
|Trait         |*H*<sup>2</sup> Cullis| *H*<sup>2</sup> Oakey|*H*<sup>2</sup> Delta| *H*<sup>2</sup> Piepho|
|--------------|------|-----|-----|------|
|TotCH         |0.184 |0.184|0.184|0.139 |
|Fv_Fm         |0.165 |0.165|0.165|0.16  |
|Sm            |0.139 |0.139|0.139|0.115 |
|Pi_Abs        |0.109 |0.109|0.109|0.079 |
|plant_weight_g|0.254 |0.254|0.254|0.228 |
|fruit_weight_g|0.259 |0.259|0.259|0.226 |

*H*<sup>2</sup> Cullis:  broad sense heritability estimated from **Model 1** fitted with **lme4** package using Cullis' formula.  
*H*<sup>2</sup>  Oakey:  broad sense heritability estimated from **Model 1** fitted with **lme4** package using  Oakey's formula.  
*H*<sup>2</sup> Delta:  broad sense heritability estimated from **Model 1** fitted with **lme4** package using Delta's formula.  
*H*<sup>2</sup> Piepho:  broad sense heritability estimated from **Model 1** fitted with **lme4** package using Piepho's formula.  





