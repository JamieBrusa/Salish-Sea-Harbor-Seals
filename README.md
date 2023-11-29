# Salish-Sea-Harbor-Seals

#### *In Review*  

#### Jamie L. Brusa, Scott F. Pearson, Martin G. Raphael, Beth Gardner 

##### Please contact the first author for questions about the code or data: Jamie Brusa (jlbwcc@gmail.com)

_______________________________________________________________________________________

## Abstract

Understanding the role environmental characteristics play in driving variation in species occurrence and abundance is important not only for basic ecological knowledge but also for informing management decisions. In particular, information on marine mammal spatial distributions can help to manage human activities, species of concern, and species of commercial importance in a dynamic system. To identify environmental drivers associated with variation in densities of harbor seals (Phoca vitulina) in the Salish Sea, Washington, USA, we analyzed 20 years of boat-based survey data and environmental covariates with a hierarchical distance sampling model. We included spatial, temporal, and spatiotemporal environmental covariates in our model and produced fine-scale predictive maps displaying these parameters.  We found that spatial covariates were the strongest predictors for harbor seal densities in the Salish Sea.  Harbor seals were more abundant closer to major river mouths, near shore, in shallower waters, and in areas with more haul-out sites.  Additionally, harbor seal density varied with shoreline type.  Changes in predicted spatial use of the Salish Sea varied but with little difference between the breeding/molting and non-breeding/non-molting seasons.  Our results revealed patterns of greater harbor seal density in the Salish Sea, which are particularly important for conservation planning.

### Required Packages
dplyr
stringr
runjags
rjags
coda
ggmcmc
here

### How to Use this Repository
The model code can be used to estimate fine-scale densities of harbor seals using distance sampling data from vessel-based surveys. We provide code for the hierarchical distance sampling model described the manuscript to use with our data, which are provided in the data folder. The model code can be adapted to fit data for future surveys in the Salish Sea or data for other studies using a similar protocol.
