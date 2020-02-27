# R package: GLMMMiRKAT

Title: A distance-based kernel association test based on the generalized linear mixed model

Version: 1.2

Date: 2020-02-27

Author: Hyunwook Koh

Maintainer: Hyunwook Koh <hkoh7@jhu.edu>

Description: This sofware package provides facilities for GLMM-MiRKAT/aGLMM-MiRKAT which tests the association betweem microbial community composition and a host trait of interest for correlated (e.g., family-based or longitudinal) microbiome studies. For the host trait of interest, continuous (e.g., body mass index), binary (e.g., disease status, treatment/placebo) or Poisson (e.g., number of tumors/treatments) traits can be handled. 

NeedsCompilation: No

Depends: R(>= 3.4.1)

Imports: CompQuadForm, dirmult, ecodist, GUniFrac, lme4, MASS, Matrix, permute, phyloseq

License: GPL-2

URL: https://github.com/hk1785/GLMM-MiRKAT

## Reference

* Koh H, Li Y, Zhan X, Chen J, Zhao N. (2019) A distance-based kernel association test based on the generalized linear mixed model for correlated microbiome studies. Front. Genet. 458(10), 1-14.
* DOI: https://doi.org/10.3389/fgene.2019.00458

## Installation

GLMMMiRKAT
```
library(devtools)
install_github("hk1785/GLMM-MiRKAT", force=T)
```
Please make suere if you have the most recent software version.
```
library(GLMMMiRKAT)
sessionInfo()
```

## Prerequites

phyloseq
```
source("https://bioconductor.org/biocLite.R")
biocLite("phyloseq")
```
CompQuadForm
```
install.packages("CompQuadForm")
```
devtools
```
install.packages("devtools")
```
ecodist
```
install.packages("ecodist")
```
GUniFrac
```
install.packages("GUniFrac")
```
lme4
```
install.packages("lme4")
```
MASS
```
install.packages("MASS")
```
Matrix
```
install.packages("Matrix")
```
permute
```
install.packages("permute")
```

## Data format

```
library(phyloseq)
URL: https://joey711.github.io/phyloseq/
```

## Troubleshooting Tips

If you have any problems for using this R package, please report in Issues (https://github.com/hk1785/GLMM-MiRKAT/issues) or email Hyunwook Koh (hkoh@jhu.edu).

* Tip 1. Depending on your pre-installed R libraries, this R package can require you to install additional R packages such as "gh", "usethis", "cli", etc using the command: install.packages("package_name").
* Tip 2. Please make sure if you have the most recent package version.

---------------------------------------------------------------------------------------------------------------------------------------

# Manual
This R package includes two core functions, Kernels and GLMMMiRKAT. Please find the details below.

## :mag: Kernels

### Description
This function creates kernel matrices based on ecological distance measures.

### Usage
```
Kernels(otu.tab, tree)
```

### Arguments
* _otu.tab_ - A matrix of OTU count table. (1. Rows are samples and columns are OTUs. 2. Monotone/singletone OTUs need to be removed.)
* _tree_ - A rooted phylogenetic tree.

### References
* Koh H, Li Y, Zhan X, Chen J, Zhao N. (2019) A distance-based kernel association test based on the generalized linear mixed model for correlated microbiome studies. Front. Genet. 458(10), 1-14.
* Jaccard P. (1912) The distribution of the flora in the alpine zone. New Phytol. 11, 37-50.
* Lozupone CA Knight R. (2005) UniFrac: a new phylogenetic method for comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-35.
* Lozupone CA et al. (2007) Quantitative and qualitative β diversity measures lead to different insights into factors that structure microbial communities. Appl. Environ. Microbiol. 73, 1576–85.
* Chen J et al. (2012) Associating microbiome composition with environmental covariates using generalized UniFrac distances. Bioinformatics 28, 2106–13. 

### Example
Import requisite R packages
```
library(GLMMMiRKAT)
library(CompQuadForm)
library(dirmult)
library(ecodist)
library(GUniFrac)
library(lme4)
library(MASS)
library(Matrix)
library(permute)
library(phyloseq)
```
Import example microbiome data
```
data(nor.phy)
```
Rarefy the microbiome data using phyloseq::rarefy_even_depth to control differing total reads per sample (recommended).
```
set.seed(100)
nor.phy <- rarefy_even_depth(nor.phy, rngseed=TRUE)

otu.tab <- otu_table(nor.phy)
tree <- phy_tree(nor.phy)
```
Create kernel matrices
```
Kernels(otu.tab, tree)
```


## :mag: GLMMMiRKAT

### Description
This function tests the association betweem microbial community composition and a host trait of interest based on correlated (e.g., family-based or longitudinal) microbiome studies. For the host trait of interest, continuous (e.g., body mass index), binary (e.g., disease status, treatment/placebo) or Poisson (e.g., number of tumors/treatments) traits can be handled. 

### Usage
```
GLMMMiRKAT(y, cov = NULL, id, time.pt = NULL, Ks, model, slope = FALSE, n.perm = 5000)
```

### Arguments
* _y_ - A numeric vector of continuous (e.g., body mass index), binary (e.g., disease status, treatment/placebo) or Poisson (e.g., number of tumors/treatments) traits.
* _cov_ - A data.frame for covariate (e.g., age, gender) adjustment(s). Default is NULL for no covariate adjustment.
* _id_ - A vector of cluster (e.g., family, longitudinal cluster) IDs.
* _time.pt_ - A vector of the time points for the longitudinal studies. 'time.pt' is not required for the random intercept model. Default is time.pt = NULL. 
* _Ks_ - A list of the kernel matrices created using GLMMMiRKAT::Kernels (?Kernels).
* _model_ - "gaussian" for continuous, "binomial" for binary or "poisson" for Poisson traits.
* _slope_ - An indicator to include random slopes in the model (slope = TRUE) or not (slope = FALSE). 'slope = FALSE' is for the random intercept model. 'slope = TRUE' is for the random slope model. For the random slope model (slope=TRUE), 'time.pt' is required. 
* _n.perm_ - A number of permutations. Default is 5000. 

### Values
_$ItembyItem_ - A vector of the estimated p-values for the item-by-item GLMM-MiRKAT analyses

_$aGLMMMiRKAT_ - The estimated p-value for the adaptive GLMM-MiRKAT (aGLMM-MiRKAT) analysis

### References
* Koh H, Li Y, Zhan X, Chen J, Zhao N. (2019) A distance-based kernel association test based on the generalized linear mixed model for correlated microbiome studies. Front. Genet. 458(10), 1-14.

### Example 1. Continuous (e.g., body mass index) traits

Import requisite R packages
```
library(GLMMMiRKAT)
library(CompQuadForm)
library(dirmult)
library(ecodist)
library(GUniFrac)
library(lme4)
library(MASS)
library(Matrix)
library(permute)
library(phyloseq)
```
Import example microbiome data with continuous traits
```
data(nor.phy)
```
Rarefy the microbiome data using phyloseq::rarefy_even_depth to control differing total reads per sample (recommended).
```
set.seed(100)
nor.phy <- rarefy_even_depth(nor.phy, rngseed=TRUE)
```		
Extract microbiome and meta information
```
otu.tab <- otu_table(nor.phy)
tree <- phy_tree(nor.phy)
meta <- sample_data(nor.phy)
y <- meta$y
id <- meta$id
x1 <- meta$x1
x2 <- meta$x2
covs <- as.data.frame(cbind(x1,x2))
covs[,2] <- as.factor(covs[,2]) 
# 'as.factor()' is needed for categorical covariates.
```		
Create kernel matrices
```
Ks <- Kernels(otu.tab, tree)
```		
Run GLMM-MiRKAT
```
GLMMMiRKAT(y, cov=covs, id=id, Ks=Ks, model="gaussian")
```

### Example 2. Binary (e.g., disease status, treatment/placebo) traits
		
Import requisite R packages
```
library(GLMMMiRKAT)
library(CompQuadForm)
library(dirmult)
library(ecodist)
library(GUniFrac)
library(lme4)
library(MASS)
library(Matrix)
library(permute)
library(phyloseq)
```
Import example microbiome data with binary traits
```
data(bin.phy)
```
Rarefy the microbiome data using phyloseq::rarefy_even_depth to control differing total reads per sample (recommended).
```
set.seed(100)
bin.phy <- rarefy_even_depth(bin.phy, rngseed=TRUE)
```		
Extract microbiome and meta information
```
otu.tab <- otu_table(bin.phy)
tree <- phy_tree(bin.phy)
meta <- sample_data(bin.phy)
y <- meta$y
id <- meta$id
x1 <- meta$x1
x2 <- meta$x2
covs <- as.data.frame(cbind(x1,x2))
covs[,2] <- as.factor(covs[,2]) 
# 'as.factor()' is needed for categorical covariates.
```		
Create kernel matrices
```
Ks <- Kernels(otu.tab, tree)
```		
Run GLMM-MiRKAT
```
GLMMMiRKAT(y, cov=covs, id=id, Ks=Ks, model="binomial")
```

### Example 3. Poisson (e.g., number of tumors/treatments) traits

Import requisite R packages
```
library(GLMMMiRKAT)
library(CompQuadForm)
library(dirmult)
library(ecodist)
library(GUniFrac)
library(lme4)
library(MASS)
library(Matrix)
library(permute)
library(phyloseq)
```
Import example microbiome data with Poisson traits
```
data(pos.phy)
```
Rarefy the microbiome data using phyloseq::rarefy_even_depth to control differing total reads per sample (recommended).
```
set.seed(100)
pos.phy <- rarefy_even_depth(pos.phy, rngseed=TRUE)
```		
Extract microbiome and meta information
```
otu.tab <- otu_table(pos.phy)
tree <- phy_tree(pos.phy)
meta <- sample_data(pos.phy)
y <- meta$y
id <- meta$id
x1 <- meta$x1
x2 <- meta$x2
covs <- as.data.frame(cbind(x1,x2))
covs[,2] <- as.factor(covs[,2]) 
# 'as.factor()' is needed for categorical covariates.
```		
Create kernel matrices
```
Ks <- Kernels(otu.tab, tree)
```		
Run GLMM-MiRKAT
```
GLMMMiRKAT(y, cov=covs, id=id, Ks=Ks, model="poisson")
```
