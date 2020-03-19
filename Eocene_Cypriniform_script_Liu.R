### step 0 install R and geomorph package
## JL uses a Mac, this is for install R and geomorph package on mac

## download  R, or upgrade R to current version
## geomorph dependency rgl require Xquartz. Download and install XQuartz (X11) http://xquartz.macosforge.org/trac if it is not already on your Mac. It will be installed in the Utilities folder. This program must be running every time you use geomorph (required by rgl). launch Xquartz to install geomorph


## install geomorph, using commends below, or use the menu "Packages&Data", check "Install dependenceis"

install.packages("geomorph", dependencies = TRUE)
##if warning and error shows, check the library location

library(geomorph)
library(RColorBrewer) #clolor panel for colorblind friendly colors
display.brewer.all(colorblindFriendly=TRUE) #use this to choose color palettes


getwd()
#check the work directory

setwd("/Users/juanliu/Box Sync/R_projects/Eocene_Cypriniform_Body_Shape")
#complete directory in the brackets and quotation marks. drag folder from Finder to Terminal to print directory path

### step 1 import data (tps files of each species)

Amyzon_aggregatum <- readland.tps("A aggregatum ALL.tps", specID = "ID", negNA = FALSE, readcurves = FALSE, warnmsg = TRUE)
# the .tps file"A aggregatum ALL.tps" is in working directory. use species name  Amyzon_aggregatum to define the landmark file

Amyzon_gosiutense <- readland.tps("A gosiutense_ALL.tps", specID = "ID", negNA = FALSE, readcurves = FALSE, warnmsg = TRUE)
# Note: if the file name contains blank, use space in file name; if use _ in file name, type in _;  space and underscore "_" are not interchangeble when import file using their file name in the work directory. 

UDG_brevippine <- readland.tps("A. brevippine.tps", specID = "ID", negNA = FALSE, readcurves = FALSE, warnmsg = TRUE)
# Note: the file name of A. brevippine contains dot "." and space. use exact spell when import file. 

Amyzon_commune <- readland.tps("A.commune ALL.tps", specID = "ID", negNA = FALSE, readcurves = FALSE, warnmsg = TRUE)

Jianghanichthys_hubeiensis <- readland.tps("HuiBei_LJ_1.tps", specID = "ID", negNA = FALSE, readcurves = FALSE, warnmsg = TRUE)

Undes_cyprinid_Hunan <- readland.tps("hunan_cyprinid_201806_TM_Scale.tps", specID = "ID", negNA = FALSE, readcurves = FALSE, warnmsg = TRUE)

Amyzon_hunanense <- readland.tps("Hunan_JJ.TPS", specID = "ID", negNA = FALSE, readcurves = FALSE, warnmsg = TRUE)


Jianghanichthys_sanshuiensis <- readland.tps("J_sanshuiensis_TM_Scale.tps", specID = "ID", negNA = FALSE, readcurves = FALSE, warnmsg = TRUE)

Amyzon_kishenehnicum <- readland.tps("Montana Amyzon.TPS", specID = "ID", negNA = FALSE, readcurves = FALSE, warnmsg = TRUE)

Tianshanicus_liui <- readland.tps("Tianshanichthys201806_TM.TPS", specID = "ID", negNA = FALSE, readcurves = FALSE, warnmsg = TRUE)


### Step 2, combine and grouping sub-dataset
library(abind)

All_Eoc_Cypr <- abind(Amyzon_aggregatum, Amyzon_gosiutense, UDG_brevippine, Amyzon_commune, Jianghanichthys_hubeiensis, Undes_cyprinid_Hunan, Amyzon_hunanense, Jianghanichthys_sanshuiensis, Amyzon_kishenehnicum, Tianshanicus_liui)
##not run, for checking data
#plotAllSpecimens(All_Eoc_Cypr)


## grouping/classify samples
classifier <- read.csv("All_Eoc_Cyp_Classifier.csv", header=T, row.names=1)
# create .csv file containing specimen ID, species, family, and distribution. import. 

#is.factor(classifier$Species)
#[1] FALSE
# check if Species is a factor, the anwer is False here, because in the csv, species was spelled as "Speceies". correct it, and re-run:
##not run, for checking data
#is.factor(classifier$Species)
#[1] TRUE
##not run, for checking data
#is.factor(classifier$ID)
#[1] FALSE
# cheick if ID is a facor 
##not run, for checking data
#is.factor(classifier$Distribution)
#[1] TRUE
#check if Distribution is factor
##not run, for checking data
#is.factor(classifier$Family)
#[1] TRUE
#check if Family a factor


### Step 3, sumperimposition and Procrustes
gpa_All_Eoc_Cypr <- gpagen(All_Eoc_Cypr)
#general Proscrustes analysis on all Eocene cypriniforms. 
# should always perform gpa on the entire dataset.

gpa_All_Eoc_Cypr$coords
# display coordinates

gpa_All_Eoc_Cypr$Csize
# display centroid size of all samples
##not run, for checking data
#plotAllSpecimens(gpa_All_Eoc_Cypr$coords)
# visulization of the GPA




### Step 4, PCA analysis and visulization

# PCA on Procrusted coordinates of all Eocene cypriniforms.
plotTangentSpace(gpa_All_Eoc_Cypr$coords) 
gp<-classifier$Species
##To change colors of groups
## col.gp <- rainbow(length(levels(gp))) 
col.gp <- brewer.pal(length(levels(gp)), "RdBu")
names(col.gp) <- levels(gp)
col.gp <- col.gp[match(gp, names(col.gp))] # col.gp must NOT be a factor
PCA_species <- plotTangentSpace(gpa_All_Eoc_Cypr$coords, groups = col.gp, legend = TRUE)
# PCA with grouping by species 

# PCA with grouping by families
gp.fam<-classifier$Family
##col.gp.fam <- rainbow(length(levels(gp.fam))) 
col.gp.fam <- brewer.pal(length(levels(gp.fam)), "RdBu")
names(col.gp.fam) <- levels(gp.fam)
col.gp.fam <- col.gp.fam[match(gp.fam, names(col.gp.fam))] 
PCA_family <- plotTangentSpace(gpa_All_Eoc_Cypr$coords, groups = col.gp.fam, legend = TRUE)

# PCA, grouping by interaction of family and distribution
DisFam <- interaction(classifier$Distribution, classifier$Family)
col.gp.DisFam <- brewer.pal(length(levels(DisFam)), "RdBu")
names(col.gp.DisFam) <- levels(DisFam)
col.gp.DisFam <- col.gp.DisFam[match(DisFam, names(col.gp.DisFam))] 
plotTangentSpace(gpa_All_Eoc_Cypr$coords, groups = col.gp.DisFam, legend = TRUE)


## pairwise comparison on shape between species, families, and geographical distribution
fit.sp <- procD.lm(coords ~ classifier$Species, data = gpa_All_Eoc_Cypr, iter = 999, print.progress = FALSE)
summary(fit.sp)
fit.fam <- procD.lm(coords ~ classifier$Family, data = gpa_All_Eoc_Cypr, iter = 999, print.progress = FALSE)
summary(fit.fam)
fit.dist <- procD.lm(coords ~ classifier$Distribution, data = gpa_All_Eoc_Cypr, iter = 999, print.progress = FALSE)
summary(fit.dist)

PW.sp <- pairwise(fit.sp, groups = gp, covariate = NULL) 
# Pairwise distances between means, summarized two ways (replaces advanced.procD.lm): 
summary(PW.sp, test.type = "dist", confidence = 0.95, stat.table = TRUE)

PW.fam <- pairwise(fit.fam, groups = classifier$Family, covariate = NULL) 
summary(PW.fam, test.type = "dist", confidence = 0.95, stat.table = TRUE)
PW.dist <- pairwise(fit.dist, groups = classifier$Distribution, covariate = NULL) 
summary(PW.dist, test.type = "dist", confidence = 0.95, stat.table = TRUE)



##PCA on consensus shape of each species
EC.gdf <- geomorph.data.frame(gpa_All_Eoc_Cypr,  taxon= classifier$taxon)
new.coords <- coords.subset(A = EC.gdf $coords, group = EC.gdf $taxon)
names(new.coords)
spmean <- lapply(new.coords, mshape)
spmean <- simplify2array(spmean, higher=TRUE)
EC.msPC <- gm.prcomp(spmean)
summary(EC.msPC)
gp<-classifier$Species
##col.gp <- rainbow(length(dimnames(spmean)[[3]])) 
col.gp <- brewer.pal(length(levels(gp)), "RdBu")
   names(col.gp) <- dimnames(spmean)[[3]]
pca.msEC <- plotTangentSpace(spmean, groups = col.gp, legend = TRUE)


# PCA on largest specimen of each species
# GM of largest specimen comparison. average or consensus shape may not represent the species because of juvinile/adult ratio varies across locality, use shape of largest individule to represent the adult shape.
# bigAA <- A_aggregatum_12, A_gosiutense_6, UDG_brevippine_1, A_commune_2, J_hubeiensis_18, UD_Cyprn_hunan_2, A_hunanense_3, J_sanshuiensis_1, A_kishenehnicum_5, T_liui_1
BigEC<- readland.tps("bigEC.TPS", specID = "ID", negNA = FALSE, readcurves = FALSE, warnmsg = TRUE)
gpa.BigEC <-gpagen(BigEC)
##col.gp <- rainbow(length(dimnames(spmean)[[3]])) 
col.gp <- brewer.pal(length(levels(gp)), "RdBu")
   names(col.gp) <- dimnames(spmean)[[3]]
pca.BigEC <- plotTangentSpace (gpa.BigEC$coords, legend = TRUE, group= col.gp)


### step 5, size-shape covariation, simple allometry

# shape-size global morphosapce of all species
fit.simp <- procD.lm(coords ~ log(Csize), data = gpa_All_Eoc_Cypr, iter = 999, print.progress = FALSE)

fit.comm <- procD.lm(coords ~ log(Csize)+classifier$Species, data = gpa_All_Eoc_Cypr, iter = 999, print.progress = FALSE)

fit.uniq <- procD.lm(coords ~ log(Csize)*classifier$Species, data = gpa_All_Eoc_Cypr, iter = 999, print.progress = FALSE)
summary(fit.simp)
summary(fit.comm)
summary(fit.uniq)


RSS <- function(fit) 
sum(diag(resid(fit)%*%t(resid(fit))))
RSS(fit.simp)
RSS(fit.comm)
RSS(fit.uniq)

PW.sp.sim <- pairwise(fit.simp, groups = gp, covariate = NULL) 
# Pairwise distances between means, summarized two ways (replaces advanced.procD.lm): 
summary(PW.sp.sim, test.type = "dist", confidence = 0.95, stat.table = TRUE)

fit.fam.simp <- procD.lm(coords ~ log(Csize)+classifier$Family, data = gpa_All_Eoc_Cypr, iter = 999, print.progress = FALSE)
PW.fam.simp <- pairwise(fit.fam.simp, groups = classifier$Family, covariate = NULL) 
summary(PW.fam.simp, test.type = "dist", confidence = 0.95, stat.table = TRUE)
fit.dist.simp <- procD.lm(coords ~ log(Csize)+classifier$Distribution, data = gpa_All_Eoc_Cypr, iter = 999, print.progress = FALSE)
PW.dist.simp <- pairwise(fit.dist.simp, groups = classifier$Distribution, covariate = NULL) 
summary(PW.dist.simp, test.type = "dist", confidence = 0.95, stat.table = TRUE)


#linear regression, species specific allometry. 
gp<-classifier$Species
##col.gp <- rainbow(length(levels(gp))) 
col.gp <- brewer.pal(length(levels(gp)), "RdBu")
names(col.gp) <- levels(gp)
col.gp <- col.gp[match(gp, names(col.gp))] # col.gp must NOT be a factor
pc.plot.ss <- plotAllometry(fit.uniq, size = gpa_All_Eoc_Cypr$Csize, logsz = TRUE, method = "size.shape", pch = 19, col = col.gp)
pc.plot.CAC <- plotAllometry(fit.uniq, size = gpa_All_Eoc_Cypr$Csize, logsz = TRUE, method = "CAC", pch = 19, col = col.gp)
plotAllometry(fit.uniq, size = gpa_All_Eoc_Cypr$Csize, logsz = TRUE, method = "CAC", pch = 19, col = col.gp)


summary(pc.plot.ss$size.shape.PCA)
summary(pc.plot.CAC$size.shape.PCA)

# method="size.shape" indication simple allometry (size and shape covariation), no matter what type of fit (defined above) won't affect the plot. 

#plot by family
gp.fam<-classifier$Family
##col.gp.fam <- rainbow(length(levels(gp.fam))) 
col.gp.fam <- brewer.pal(length(levels(gp.fam)), "RdBu")
names(col.gp.fam) <- levels(gp.fam)
col.gp.fam <- col.gp.fam[match(gp.fam, names(col.gp.fam))] 
pc.plot.fam <- plotAllometry(fit, size = gpa_All_Eoc_Cypr$Csize, logsz = TRUE, method = "size.shape", pch = 19, col = col.gp.fam)
summary(pc.plot.fam$size.shape.PCA)

#plot by distributionVSfamily
gp.fam<-classifier$Family
##col.gp.fam <- rainbow(length(levels(gp.fam))) 
DisFam <- interaction(classifier$Distribution, classifier$Family)
col.gp.DisFam <- brewer.pal(length(levels(DisFam)), "RdBu")
names(col.gp.DisFam) <- levels(DisFam)
col.gp.DisFam <- col.gp.DisFam[match(DisFam, names(col.gp.DisFam))] 
pc.plot.fam <- plotAllometry(fit, size = gpa_All_Eoc_Cypr$Csize, logsz = TRUE, method = "size.shape", pch = 19, col = col.gp.DisFam)
summary(pc.plot.fam$size.shape.PCA)




### Step 6 allometry of Amyzon
Amyzon <- abind(Amyzon_aggregatum, Amyzon_gosiutense, Amyzon_commune, Amyzon_hunanense, Amyzon_kishenehnicum)
gpa_Amyzon <- gpagen(Amyzon)
##not run, for checking data
#gpa_Amyzon $coords
##not run, for checking data
#gpa_Amyzon $Csize
Am_classifier <- read.csv("Am_Classifier.csv", header=T, row.names=1)
#PCA of Amyzon
Am.gp<-Am_classifier$Species
##To change colors of groups
## col.gp <- rainbow(length(levels(gp))) 
col.Am.gp <- brewer.pal(length(levels(Am.gp)), "RdBu")
names(col.Am.gp) <- levels(Am.gp)
col.Am.gp <- col.Am.gp[match(Am.gp, names(col.Am.gp))] # col.gp must NOT be a factor
PCA.Am <- plotTangentSpace(gpa_Amyzon$coords, groups = col.Am.gp, legend = TRUE)
Am_fit1 <- procD.lm(coords ~ log(Csize)+ Am_classifier$Species, data = gpa_Amyzon, iter = 999, print.progress = FALSE)
summary(Am_fit1)

Am_fit2 <- procD.lm(coords ~ log(Csize) * Am_classifier$Species, data = gpa_Amyzon, iter = 999, print.progress = FALSE)
summary(Am_fit2)

RSS <- function(fit) 
sum(diag(resid(fit)%*%t(resid(fit))))
RSS(Am_fit1)
RSS(Am_fit2)



agp <- Am_classifier$Species
#col.agp <- rainbow(length(levels(agp))) 
col.agp <- brewer.pal(length(levels(agp)), "RdBu")
names(col.agp) <- levels(agp)
col.agp <- col.agp[match(agp, names(col.agp))] # col.gp must NOT be a factor
Am_cac <- plotAllometry (Am_fit1, size = gpa_Amyzon$Csize, logsz = TRUE, method = "CAC", pch = 19, col = col.agp)
Am_ss <- plotAllometry (Am_fit1, size = gpa_Amyzon$Csize, logsz = TRUE, method = "size.shape", pch = 19, col = col.agp)
##Am_ss2 <- plotAllometry (Am_fit2, size = gpa_Amyzon$Csize, logsz = TRUE, method = "size.shape", pch = 19, col = col.agp) # not run, size.shape generate same plot, no matter which fit


# pca on residules of regression
Am_fit1.resid <- Am_fit1$residuals
Am.fit1.resid.pca <- prcomp(Am_fit1.resid)
summary(Am.fit1.resid.pca)
plot(Am.fit1.resid.pca$x[,1:2], col=col.agp, pch=19, cex=1.5) #use pch=c(1:5) for different symbols
legend(0.07, -0.035, legend=unique(names(col.agp)), col=unique(col.agp), pch=19, cex=0.8)



### Step 7 phylo-morph
library(ape)

# import tree, alternative function read.tree()
tree <- read.nexus(file="EoCy_tree.nex")
## is.phylo(tree) #zjt: this is a function in phylosim package, which is not installed

EC.gdf <- geomorph.data.frame(gpa_All_Eoc_Cypr,  taxon= classifier$taxon)
new.coords <- coords.subset(A = EC.gdf $coords, group = EC.gdf $taxon)
names(new.coords)
spmean <- lapply(new.coords, mshape)
spmean <- simplify2array(spmean, higher=TRUE)
#test for phylogenetic signal
#physignal(A=spmean , phy= EC.gdf $tree, iter = 999, seed = NULL, print.progress = TRUE) #"EC.dgf$tree" is undefined
tree <- compute.brlen(tree, 1) #set uniform branch length, as it was previously undefined during import.
physignal(A=spmean, phy= tree, iter = 999, seed = NULL, print.progress = TRUE)

#Plot phylomorphospace using mean shape
EC.PC.phylo <- gm.prcomp(spmean, phy=tree)
#zjt plot using spmean and resoling multichotomies
bi.tree <- multi2di(tree, random=F)
plotGMPhyloMorphoSpace(bi.tree, spmean)
 #EC.pgls <- procD.pgls(coords ~ Csize, phy = tree, data = gpa_All_Eoc_Cypr, iter = 999) #this calls a dataset that is different from tree; the "gpa_All_Eoc_Cypr" dataset is all specimens, not mean. PGLS needs match between tree and dataset specimen sample size
EC.gdf.cs <- EC.gdf$Csize
names(EC.gdf.cs) <- classifier$taxon
spmean.cs <- aggregate(x = EC.gdf.cs,                # Specify data column
          by = list(names(EC.gdf.cs)),              # Specify group indicator
          FUN = mean) 
#check if mean csize and mean shape entries are identically ordered
dimnames(spmean)[[3]] == spmean.cs$Group.1


EC.pgls <- procD.pgls(spmean ~ spmean.cs$x, phy=tree, iter=999)
anova(EC.pgls) 

#ANOVA wih phylo, not working
#EC.pgls.cs <- procD.pgls(spmean$coords~ Csize+spmean.cs$Group.1, data= EC.PC.phylo, phy=tree, iter=999)
#procD.pgls(spmean ~ spmean.cs$x+Csize, phy=tree, iter=999)
#summary(EC.pgls.cs)anova(EC.pgls) 
## summary(EC.pgls) return the same result as anova(EC.pgls)

##phylo-morph using largest specimen
big.classifier <- read.csv("bigEC_Classifier.CSV", header=T, row.names=1)
#bigEC.gdf <- geomorph.data.frame(gpa.BigEC,  taxon= big.classifier $taxon)
big.classifier$name <- rownames(big.classifier)

BigEC<- readland.tps("bigEC.TPS", specID = "ID", negNA = FALSE, readcurves = FALSE, warnmsg = TRUE)

dimnames(BigEC)[[3]] <- trimws(dimnames(BigEC)[[3]])

ordered.big.classifier <- big.classifier[match(dimnames(BigEC)[[3]], big.classifier$name),]

dimnames(BigEC)[[3]] <- big.classifier$taxon

tree <- compute.brlen(tree, 1) #set uniform branch length, as it was previously undefined during import.

physignal(A= BigEC, phy= tree, iter = 999, seed = NULL, print.progress = TRUE)

gpa.BigEC <-gpagen(BigEC)
bigEC.PC.phylo <- gm.prcomp(gpa.BigEC$coords, phy=tree)


#ANOVA wih phylo
bigEC.pgls <- procD.pgls(coords ~ big.classifier$Species, phy=tree, data= gpa.BigEC, iter=999)
summary(bigEC.pgls)
#physi signal in size
bigEC.pgls.reg <- procD.pgls(coords ~ Csize, phy=tree, data= gpa.BigEC, iter=999)
plotGMPhyloMorphoSpace(bi.tree, gpa.BigEC$coords)
summary(bigEC.pgls.reg)

#ANOVA wih phylo by family
bigEC.pgls.fam <- procD.pgls(coords ~ big.classifier$Family, phy=tree, data= gpa.BigEC, iter=999)
summary(bigEC.pgls.fam)
bigEC.pgls.fam.reg <- procD.pgls(coords ~ Csize+big.classifier$Family, phy=tree, data= gpa.BigEC, iter=999)
summary(bigEC.pgls.fam.reg)
