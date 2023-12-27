library(dplyr)
library(ggplot2)
library(stringr)
library(runjags)
library(rjags)
library(coda)
library(ggmcmc)
library(here)

#Bring in data and do some clean up (work-around to have large dataset on Github)
setwd(here('Data', 'Partial data'))
seals.9 <- list.files(path = here('Data', 'Partial data'), pattern = "09")
seals.9 <- do.call("rbind", lapply(seals.9, FUN=function(files){read.csv(files)}))
colnames(seals.9)[2] <- "StratPSUSeg"
seals.11 <- list.files(path = here('Data', 'Partial data'), pattern = "11")
seals.11 <- do.call("rbind", lapply(seals.11, FUN=function(files){read.csv(files)}))
colnames(seals.11)[2] <- "mean_Length"
seals.14 <- list.files(path = here('Data', 'Partial data'), pattern = "14")
seals.14 <- do.call("rbind", lapply(seals.14, FUN=function(files){read.csv(files)}))
colnames(seals.14)[2] <- "MeanDepth"
seals.19 <- list.files(path = here('Data', 'Partial data'), pattern = "19")
seals.19 <- do.call("rbind", lapply(seals.19, FUN=function(files){read.csv(files)}))
colnames(seals.19)[2] <- "StratPSUSegYJ"
seals.33 <- list.files(path = here('Data', 'Partial data'), pattern = "33")
seals.33 <- do.call("rbind", lapply(seals.33, FUN=function(files){read.csv(files)}))
colnames(seals.33)[2] <- "RM.dist"
seals.34 <- list.files(path = here('Data', 'Partial data'), pattern = "34")
seals.34 <- do.call("rbind", lapply(seals.34, FUN=function(files){read.csv(files)}))
colnames(seals.34)[2] <- "SST"
seals.35 <- list.files(path = here('Data', 'Partial data'), pattern = "35")
seals.35 <- do.call("rbind", lapply(seals.35, FUN=function(files){read.csv(files)}))
colnames(seals.35)[2] <- "chlorophyll"
seals.36 <- list.files(path = here('Data', 'Partial data'), pattern = "36")
seals.36 <- do.call("rbind", lapply(seals.36, FUN=function(files){read.csv(files)}))
colnames(seals.36)[2] <- "Salinity"
seals.37 <- list.files(path = here('Data', 'Partial data'), pattern = "37")
seals.37 <- do.call("rbind", lapply(seals.37, FUN=function(files){read.csv(files)}))
colnames(seals.37)[2] <- "area"
seals.parts <- cbind(seals.9, seals.11, seals.14, seals.19, seals.33, seals.34, seals.35, seals.36, seals.37)
seals.parts[, "X" == names(seals.parts)] <- NULL
setwd(here('Data', 'SS Harbor seal data'))
path <- here('Data', 'SS Harbor seal data')
multiseal <- list.files(path = path, pattern = "csv")
multiseal_df <- do.call("cbind",lapply(multiseal,FUN=function(files){ read.csv(files)}))
multiseal_df[, "X" == names(multiseal_df)] <- NULL
colnames(multiseal_df) <- c("YYYY", "ESIcode", "NPGO", "Upwelling", "Spring.Transition.JD", "Fall.Transition.JD", 
                            "dclass", "MM", "shore7", "shore5", "shore9A", "shore1A", "shore2A", "shore6A", "shore4", 
                            "shore6D", "shore8A", "numeric_pairs", "N", "HCount", "Perp_Dist", "Stratum", "PSU",
                            "Segment", "Beaufort")
multiseal_df <- cbind(multiseal_df, seals.parts)

#Figure out how many sites
nSites <- length(unique(multiseal_df$StratPSUSeg)) #number of segments/zigzags

strip.width <- 275 #95% quantile for distances from 2001-2018
int.w <- 27 #binwidth for distances histogram
dist.breaks <- seq(0, strip.width, by=int.w) #setting up bins
nG<-length(dist.breaks)-1 #number of bins/rectangles, 12 for seals and 1 for no seals
#Give each rectangle an identification number
multiseal_df$dclass<-with(multiseal_df,ifelse(Perp_Dist<=dist.breaks[1],1,ifelse(Perp_Dist > dist.breaks[1]
                                                                                 & Perp_Dist <= dist.breaks[2], 2,
                                                                                 ifelse(Perp_Dist > dist.breaks[2]
                                                                                 & Perp_Dist <= dist.breaks[3], 3,
                                                                                 ifelse(Perp_Dist > dist.breaks[3]
                                                                                 & Perp_Dist <= dist.breaks[4], 4,
                                                                                 ifelse(Perp_Dist > dist.breaks[4]
                                                                                 & Perp_Dist <= dist.breaks[5], 5,
                                                                                 ifelse(Perp_Dist > dist.breaks[5]
                                                                                 & Perp_Dist <= dist.breaks[6], 6,
                                                                                 ifelse(Perp_Dist > dist.breaks[6]
                                                                                 & Perp_Dist <= dist.breaks[7], 7,
                                                                                 ifelse(Perp_Dist > dist.breaks[7]
                                                                                 & Perp_Dist <= dist.breaks[8], 8,
                                                                                 ifelse(Perp_Dist > dist.breaks[8]
                                                                                 & Perp_Dist <= dist.breaks[9], 9,
                                                                                 ifelse(Perp_Dist > dist.breaks[9]
                                                                                 & Perp_Dist <= dist.breaks[10], 10,
                                                                                 ifelse(Perp_Dist > dist.breaks[10]
                                                                                 & Perp_Dist <= dist.breaks[11], 11,
                                                                                 12))))))))))))


##negative binomial dispersion parameter in the Poisson-Gamma formulation
r <- 1


###########Full model##############
#Make an array with 1 species x 129 transects x 12 distance classes and populate with N
transect <- unique(multiseal_df$StratPSUSegYJ)
nSites <- length(unique(multiseal_df$StratPSUSegYJ))
distclass <- unique(multiseal_df$dclass)
distclass <- sort(distclass)

seals <- array(0, c("1", length(transect), length(distclass)-1))

for(t in 1:nG){
  for(c in 1:(nSites)){
    tcall <- distclass[t]
    ccall <- transect[c]
    vec <- multiseal_df$N[multiseal_df$StratPSUSegYJ==ccall & multiseal_df$dclass==tcall]
    seals[1,c,t] <- ifelse(length(vec) == 0, rep(0,nSites), vec)
  }}

N <- multiseal_df %>%
  group_by(StratPSUSegYJ) %>%
  summarise(sumN = round(sum(N),0))

N <- N$sumN

nind <- sum(seals)

seals.sum<-apply(seals,1:2, sum)


site<-dclass<-NULL

for(j in 1:nSites){
  for (k in 1:nG){
    if (seals[1,j,k]==0) next
    site<-c(site, rep(j, seals[1,j,k]))
    dclass<-c(dclass, rep(k, seals[1,j,k]))
    
  }}

Year <- multiseal_df %>% group_by(StratPSUSegYJ) %>%
  summarise(Year = mean(YYYY))

Year <- Year$Year - 2000

multiseal_df$StratPSU <- paste(multiseal_df$Stratum, multiseal_df$PSU)

SP<- unique(multiseal_df$StratPSU)
dumb <- data.frame(SP = SP,
                   numbers = 1:length(SP))

multiseal_df$numeric_SP <- 0
for(i in dumb$SP){
  multiseal_df$numeric_SP[multiseal_df$StratPSU == i] <- match(i, dumb$SP)
}

StratPSU <- multiseal_df %>% group_by(StratPSUSegYJ) %>%
  summarise(numeric_SP = max(numeric_SP))

StratPSU <- StratPSU$numeric_SP

pair <- multiseal_df %>% group_by((StratPSUSegYJ)) %>%
  summarise(Pair = max(numeric_pairs))

pair <- pair$Pair


dclass[dclass == 1] <- 2  #assign all the 1's to be 2s 
dclass=dclass-1           #then subtract 1 from all to start at category 1

# #Need to have 1 obs per segment Year JD
multiseal_df$BSS.1 <- with(multiseal_df, ifelse(multiseal_df$Beaufort == 1, 1, 0))
BSS.1 <- multiseal_df %>% group_by(StratPSUSegYJ) %>% summarise(BSS.1 = round(mean(BSS.1)))
multiseal_df$BSS.2 <- with(multiseal_df, ifelse(multiseal_df$Beaufort == 2, 1, 0))
BSS.2 <- multiseal_df %>% group_by(StratPSUSegYJ) %>% summarise(BSS.2 = round(mean(BSS.2)))
multiseal_df$BSS.3 <- with(multiseal_df, ifelse(multiseal_df$Beaufort == 3 | multiseal_df$Beaufort == 4, 1, 0))
BSS.3 <- multiseal_df %>% group_by(StratPSUSegYJ) %>% summarise(BSS.3 = round(mean(BSS.3)))
BSS.1 = BSS.1$BSS.1
BSS.2 = BSS.2$BSS.2
BSS.3 = BSS.3$BSS.3


haul.count <- multiseal_df %>% group_by(StratPSUSegYJ) %>%
  summarise(haul.count = mean(HCount))
haul.count <- haul.count$haul.count
#Center/scale 
meanh <- mean(haul.count)
h2sd <- 2*sd(haul.count)
for(i in 1:length(haul.count)){
  haul.count[i] <- (haul.count[i] - meanh)/h2sd
}

#Set up categorical shoretype covariate with ESI code 8C as the intercept
multiseal_df$shore7 <- with(multiseal_df, ifelse(multiseal_df$ESIcode == "7", 1, 0))
shore7 <- multiseal_df %>% group_by(StratPSUSegYJ) %>% summarise(shore7 = mean(shore7))
multiseal_df$shore5 <- with(multiseal_df, ifelse(multiseal_df$ESIcode == "5", 1, 0))
shore5 <- multiseal_df %>% group_by(StratPSUSegYJ) %>% summarise(shore5 = mean(shore5))
multiseal_df$shore9A <- with(multiseal_df, ifelse(multiseal_df$ESIcode == "9A", 1, 0))
shore9A <- multiseal_df %>% group_by(StratPSUSegYJ) %>% summarise(shore9A = mean(shore9A))
multiseal_df$shore1A <- with(multiseal_df, ifelse(multiseal_df$ESIcode == "1A", 1, 0))
shore1A <- multiseal_df %>% group_by(StratPSUSegYJ) %>% summarise(shore1A = mean(shore1A))
multiseal_df$shore2A <- with(multiseal_df, ifelse(multiseal_df$ESIcode == "2A", 1, 0))
shore2A <- multiseal_df %>% group_by(StratPSUSegYJ) %>% summarise(shore2A = mean(shore2A))
multiseal_df$shore6A <- with(multiseal_df, ifelse(multiseal_df$ESIcode == "6A", 1, 0))
shore6A <- multiseal_df %>% group_by(StratPSUSegYJ) %>% summarise(shore6A = mean(shore6A))
multiseal_df$shore4 <- with(multiseal_df, ifelse(multiseal_df$ESIcode == "4", 1, 0))
shore4 <- multiseal_df %>% group_by(StratPSUSegYJ) %>% summarise(shore4 = mean(shore4))
multiseal_df$shore6D <- with(multiseal_df, ifelse(multiseal_df$ESIcode == "6D", 1, 0))
shore6D <- multiseal_df %>% group_by(StratPSUSegYJ) %>% summarise(shore6D = mean(shore6D))
multiseal_df$shore8A <- with(multiseal_df, ifelse(multiseal_df$ESIcode == "8A", 1, 0))
shore8A <- multiseal_df %>% group_by(StratPSUSegYJ) %>% summarise(shore8A = mean(shore8A))
shore7 = shore7$shore7
shore5 = shore5$shore5
shore9A = shore9A$shore9A
shore1A = shore1A$shore1A
shore2A = shore2A$shore2A 
shore6A = shore6A$shore6A
shore4 = shore4$shore4
shore6D = shore6D$shore6D 
shore8A = shore8A$shore8A

multiseal_df$offshore <- with(multiseal_df, ifelse(str_detect(multiseal_df$Segment, "Z"), 1, 0))
Offshore <- multiseal_df %>% group_by(StratPSUSegYJ) %>% summarise(offshore = mean(offshore))
Offshore <- Offshore$offshore

River <- multiseal_df %>% group_by(StratPSUSegYJ) %>%
  summarise(River = mean(RM.dist))
River <- River$River
#Center/scale 
meanRM <- mean(River)
Rsd <- 2*sd(River)
for(i in 1:length(River)){
  River[i] <- (River[i] - meanRM)/Rsd
}

NPGO <- multiseal_df %>% group_by(StratPSUSegYJ) %>%
  summarise(NPGO = mean(NPGO))
NPGO <- NPGO$NPGO
#Center/scale 
meanN <- mean(NPGO)
N2sd <- 2*sd(NPGO)
for(i in 1:length(NPGO)){
  NPGO[i] <- (NPGO[i] - meanN)/N2sd
}

MeanDepth <- multiseal_df %>% group_by(StratPSUSegYJ) %>%
  summarise(MeanDepth = mean(MeanDepth))
MeanDepth <- MeanDepth$MeanDepth
#Center/scale 
meanD <- mean(MeanDepth)
D2sd <- 2*sd(MeanDepth)
for(i in 1:length(MeanDepth)){
  MeanDepth[i] <- (MeanDepth[i] - meanD)/D2sd
}

Upwell <- multiseal_df %>% group_by(StratPSUSegYJ) %>%
  summarise(Upwell = mean(Upwelling))
Upwell <- Upwell$Upwell
#Center/scale 
meanU <- mean(Upwell)
U2sd <- 2*sd(Upwell)
for(i in 1:length(Upwell)){
  Upwell[i] <- (Upwell[i] - meanU)/U2sd
}

ST <- multiseal_df %>% group_by(StratPSUSegYJ) %>%
  summarise(ST = mean(Spring.Transition.JD))
ST <- ST$ST
#Center/scale 
meanST <- mean(ST)
ST2sd <- 2*sd(ST)
for(i in 1:length(ST)){
  ST[i] <- (ST[i] - meanST)/ST2sd
}

FT <- multiseal_df %>% group_by(StratPSUSegYJ) %>%
  summarise(FT = mean(Fall.Transition.JD))
FT <- FT$FT
#Center/scale 
meanFT <- mean(FT)
FT2sd <- 2*sd(FT)
for(i in 1:length(FT)){
  FT[i] <- (FT[i] - meanFT)/FT2sd
}

sst <- multiseal_df %>% group_by(StratPSUSegYJ) %>%
  summarise(sst = mean(SST))
sst <- sst$sst
#Center/scale 
meansst <- mean(sst)
sst2sd <- 2*sd(sst)
for(i in 1:length(sst)){
  sst[i] <- (sst[i] - meansst)/sst2sd
}

chl <- multiseal_df %>% group_by(StratPSUSegYJ) %>%
  summarise(chl = mean(chlorophyll))
chl <- chl$chl
#Center/scale 
meanchl <- mean(chl)
chl2sd <- 2*sd(chl)
for(i in 1:length(chl)){
  chl[i] <- (chl[i] - meanchl)/chl2sd
}

sal <- multiseal_df %>% group_by(StratPSUSegYJ) %>%
  summarise(sal = mean(Salinity))
sal <- sal$sal
#Center/scale 
meansal <- mean(sal)
sal2sd <- 2*sd(sal)
for(i in 1:length(sal)){
  sal[i] <- (sal[i] - meansal)/sal2sd
}

area <- multiseal_df %>% group_by(StratPSUSegYJ) %>%
  summarise(area = mean(area/(5000*550)))
area <- area$area

#Set up breeding + molting/non-breeding (non-molting) seasons 
multiseal_df$breedmolt <- ifelse(multiseal_df$MM == 6 | multiseal_df$MM == 7 | multiseal_df$MM == 8 |
                                   multiseal_df$MM == 9, 1, 0)
breedmolt <- multiseal_df %>% group_by(StratPSUSegYJ) %>%
  summarise(breedmolt = mean(breedmolt))
breedmolt <- breedmolt$breedmolt


dataCovs<-list(nG=nG, xg=dist.breaks[-1]-13.5, nsites=nSites, Year = Year, pair = pair,
               pi=rep(1/(length(dist.breaks)-1), length(dist.breaks)-1), 
               nind=nind, dclass=dclass, y=as.vector(t(seals.sum)), sst = sst,
               haul.count = haul.count, NPGO = NPGO, MeanDepth = MeanDepth, River = River,
               Upwell = Upwell, ST = ST, FT = FT, shore7 = shore7, shore5 = shore5,
               shore9A = shore9A, shore1A = shore1A, shore2A = shore2A, shore6A = shore6A,
               shore4 = shore4, shore6D = shore6D, shore8A = shore8A, chl = chl, sal = sal,
               BSS.1 = BSS.1, BSS.2 = BSS.2, BSS.3=BSS.3, breedmolt = breedmolt, Offshore = Offshore,
               area = area, site = site, SP = StratPSU)


### initial values for N
N.in<-t(seals.sum)+3

initsCovs<-function(){list(N=as.vector(N.in), alpha = runif(1,1, 3), sigma0 = runif(1, 0, 7),
                           beta.bss.1=rnorm(1, 0, 0.1),
                           beta.bss.2=rnorm(1, 0, 0.1),
                           beta.bss.3=rnorm(1, 0, 0.1),
                           beta.shore7=rnorm(1, 0, 0.1),
                           beta.shore5=rnorm(1, 0, 0.1),
                           beta.shore9A=rnorm(1, 0, 0.1),
                           beta.shore1A=rnorm(1, 0, 0.1),
                           beta.shore2A=rnorm(1, 0, 0.1),
                           beta.shore6A=rnorm(1, 0, 0.1),
                           beta.shore4=rnorm(1, 0, 0.1),
                           beta.shore6D=rnorm(1, 0, 0.1),
                           beta.shore8A=rnorm(1, 0, 0.1),
                           beta.river=rnorm(1, 0, 0.1),
                           beta.NPGO=rnorm(1,0, 0.1),
                           beta.Depth=rnorm(1,0, 0.1),
                           beta.upwell=rnorm(1,0, 0.1),
                           beta.ST=rnorm(1,0, 0.1),
                           beta.FT=rnorm(1,0, 0.1),
                           beta.sst=rnorm(1, 0, 0.1),
                           beta.chl=rnorm(1, 0, 0.1),
                           beta.sal=rnorm(1, 0, 0.1),
                           beta.BM = rnorm(1,0, 0.1),
                           beta.offshore = rnorm(1, 0, 0.1))}

params.Covs<-c('alpha', 'sigma0', 'Bp.N', 'Nseals', 'beta.NPGO', 'beta.Depth', 'beta.ho',
               'beta.upwell', 'sigma.eps.year', 'beta.shore7', 'beta.shore5', 'sigma.eps.pair',
               'beta.shore9A', 'beta.shore1A', 'beta.shore2A', 'beta.shore6A', 'beta.FT', 'beta.ST',
               'beta.shore4', 'beta.shore6D', 'beta.shore8A', 'beta.river', 'Bp.N',
               'beta.sst', 'beta.chl', 'beta.sal', 'beta.BM', 'beta.bss.1', 'beta.bss.2',
               'beta.bss.3', 'Bp.Obs', 'r.N', 'eps.year', 'beta.offshore', 'sigma.eps.psu', 'eps.PSU')

setwd(here())
modelFileCovs='NegBinomBreedMolt.txt'


covs.mod<-run.jags(model = modelFileCovs,
                   monitor = params.Covs,
                   data = dataCovs,
                   n.chains = 3,
                   adapt = 2000,
                   burnin = 10000,
                   sample = 50000,
                   inits = initsCovs,
                   thin = 1)


#Diagnostics
covs.mod_mcmc <- as.mcmc.list(covs.mod)
covs.mod_ggs <- ggs(covs.mod_mcmc)
ggs_geweke(covs.mod_ggs)
ggs_Rhat(covs.mod_ggs)

#Can swap covariates below for more investigation
ggs_traceplot(covs.mod_ggs, c("beta.ho"))
ggs_traceplot(covs.mod_ggs, c("alpha"))
ggs_traceplot(covs.mod_ggs, c("N"))
ggs_traceplot(covs.mod_ggs, c("sigma"))

ggs_density(covs.mod_ggs, c("beta.b"))
ggs_density(covs.mod_ggs, c("alpha"))
ggs_density(covs.mod_ggs, c("N"))
ggs_density(covs.mod_ggs, c("sigma"))






