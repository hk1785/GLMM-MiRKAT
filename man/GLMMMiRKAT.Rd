\name{GLMMMiRKAT}
\alias{GLMMMiRKAT}
\title{
The microbiome regression-based kernel association test based on the generalized linear mixed model
}
\description{
This function tests the association betweem microbial community (e.g., bacterial kingdom) composition and a host trait of interest (e.g., health/disease status, medical intervention, behavioral/environmental factor) based on correlated (e.g., family-based or longitudinal) microbiome studies. For the host trait of interest, Gaussian (e.g., body mass index), Binomial (e.g., disease status, treatment/placebo) or Poisson (e.g., number of tumors/treatments) traits can be handled.
}
\usage{
GLMMMiRKAT(y, cov = NULL, id, time.pt = NULL, Ks, model, slope = FALSE, n.perm = 5000)
}
\arguments{
  \item{y}{
A numeric vector of Gaussian (e.g., body mass index), Binomial (e.g., disease status, treatment/placebo) or Poisson (e.g., number of tumors/treatments) traits.
}
  \item{cov}{
A data.frame for covariate (e.g., age, gender) adjustment(s). Default is cov = NULL for no covariate adjustment.
}
  \item{id}{
A vector of cluster (e.g., family, subject including repeated measurements) IDs.
}
  \item{time.pt}{
A vector of time points for the longitudinal studies. 'time.pt' is not required (i.e., 'time.pt = NULL') for the random intercept model. Default is time.pt = NULL. 
}
  \item{Ks}{
A list of the kernel matrices created using GLMMMiRKAT::Kernels (?Kernels).
}
  \item{model}{
"gaussian" for Gaussian outcomes, "binomial" for Binomial outcomes, "poisson" for Poisson outcomes.
}
  \item{slope}{
An indicator to include random slopes in the model (slope = TRUE) or not (slope = FALSE). 'slope = FALSE' is for the random intercept model. 'slope = TRUE' is for the random slope model. For the random slope model (slope = TRUE), 'time.pt' is required. 
}
  \item{n.perm}{
A number of permutations. Default is 5000.  
}
}
\value{
ItembyItem: A vector of the estimated p-values for the item-by-item GLMM-MiRKAT analyses

aGLMMMiRKAT: The estimated p-value for the adaptive GLMM-MiRKAT (aGLMM-MiRKAT) analysis
}
\references{
Koh H, Li Y, Zhan X, Chen J, Zhao N. (2019) A distance-based kernel association test based on the generalized linear mixed model for correlated microbiome studies. Front. Genet. 458(10), 1-14.
}
\author{
Hyunwook Koh
}

\examples{
# Import requisite R packages
require(CompQuadForm)
require(dirmult)
require(ecodist)
require(GUniFrac)
require(lme4)
require(MASS)
require(Matrix)
require(permute)
require(phyloseq)

## Example 1. Gaussian (e.g., body mass index) traits
		
# Import example microbiome data with Gaussian traits
data(nor.phy)

# Rarefy the microbiome data using phyloseq::rarefy_even_depth 
# to control differing total reads per sample (recommended).
set.seed(100)
nor.phy <- rarefy_even_depth(nor.phy, rngseed=TRUE)
		
# Extract microbiome and meta information
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

# Create kernel matrices
Ks <- Kernels(otu.tab, tree)
		
# Run GLMM-MiRKAT
GLMMMiRKAT(y, cov=covs, id=id, Ks=Ks, model="gaussian")
	
## Example 2. Binomial (e.g., disease status, treatment/placebo) traits
		
# Import microbiome data with binomial traits
data(bin.phy)

# Rarefy the microbiome data using phyloseq::rarefy_even_depth 
# to control differing total reads per sample (recommended).
set.seed(100)
bin.phy <- rarefy_even_depth(bin.phy, rngseed=TRUE)

# Extract microbiome and meta information
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
		
# Create kernel matrices
Ks <- Kernels(otu.tab, tree)
		
# Run GLMM-MiRKAT
GLMMMiRKAT(y, cov=covs, id=id, Ks=Ks, model="binomial")

## Example 3. Poisson (e.g., number of tumors/treatments) traits
		
# Import microbiome data with Poisson traits
data(pos.phy)

# Rarefy the microbiome data using phyloseq::rarefy_even_depth 
# to control differing total reads per sample (recommended).
set.seed(100)
pos.phy <- rarefy_even_depth(pos.phy, rngseed=TRUE)

# Extract microbiome and meta information
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

# Create kernel matrices
Ks <- Kernels(otu.tab, tree)
	
# Run GLMM-MiRKAT
GLMMMiRKAT(y, cov=covs, id=id, Ks=Ks, model="poisson")

}