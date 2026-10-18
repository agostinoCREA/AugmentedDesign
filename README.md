# AugmentedDesign

## Background
In an augmented experimental design, a set of checks are replicated across blocks, while the test genotypes are evaluated only once or a few times. 
This is common in plant breeding when limited experimental resources are available for field trials.
The main challenge is to separate variation due to block and genotype effects and then estimate heritability.

## Model
A possible model for analyzing an augmented design is 

$y_{ij}=\mu+B_j+G_i+\varepsilon_{ij}$

where 
$y_{ij}$: observed phenotype of genotype $i4 in block $j$.  
$\(\mu\): overall mean.  
$\(B_j\): block effect  
$\(G_i\): genetic effect of genotype \(i\)  
$\(\varepsilon_{ij}\): residual error.  
