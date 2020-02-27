\name{Kernels}
\alias{Kernels}
\title{
Create kernel matrices based on ecological distance measures
}
\description{
This function creates kernel matrices based on ecological distance measures.
}
\usage{
Kernels(otu.tab, tree)
}
\arguments{
  \item{otu.tab}{
A matrix of the OTU count table. (1. Rows are samples and columns are OTUs. 2. Monotone/singletone OTUs need to be removed.)
}
  \item{tree}{
A rooted phylogenetic tree.
}
}
\value{
Jaccard: For the use of Jaccard dissimilarity (Jaccard, 1912)

BC: For the use of Bray-Curtis dissimilarity (Bray and Curtis, 1957)

U.UniFrac: For the use of unweighted UniFrac distance (Lozupone and Knight, 2005)

G.UniFrac.0.5: For the use of generalized UniFrac distance(0.5) (Chen et al., 2012)

W.UniFrac: For the use of weighted UniFrac distance (Lozupone et al., 2007)

}
\references{
Koh H, Li Y, Zhan X, Chen J, Zhao N. (2019) A distance-based kernel association test based on the generalized linear mixed model for correlated microbiome studies. Front. Genet. 458(10), 1-14.

Jaccard P. (1912) The distribution of the flora in the alpine zone. New Phytol. 11, 37-50.

Lozupone CA, Knight R. (2005). UniFrac: a new phylogenetic method for comparing microbial communities. Appl. Environ. Microbiol. 71, 8228-35.

Lozupone CA, et al. (2007). Quantitative and qualitative β diversity measures lead to different insights into factors that structure microbial communities. Appl. Environ. Microbiol. 73, 1576–85.

Chen J, et al. (2012). Associating microbiome composition with environmental covariates using generalized UniFrac distances. Bioinformatics 28, 2106–13. 
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

# Import example microbiome data
data(nor.phy)

# Rarefy the microbiome data using phyloseq::rarefy_even_depth 
# to control differing total reads per sample (recommended).
set.seed(100)
nor.phy <- rarefy_even_depth(nor.phy, rngseed=TRUE)

otu.tab <- otu_table(nor.phy)
tree <- phy_tree(nor.phy)

# Create kernel matrices
Kernels(otu.tab, tree)
	
}