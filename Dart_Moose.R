library(dartR)
library(adegenet)
library(magrittr)

####Bring in DArT Data
#2 SNP data
MooseDArT <- gl.read.dart(filename="New_All_SNP_singlerow_Chase_191003.csv", covfilename = "covfile.csv")

glMoose_clean <- 
  #Quality check
  gl.filter.repavg(MooseDArT,t=1) %>% 
  #1 SNP per locus
  gl.filter.secondaries(method = "best") %>%    
  #Loci with missing data
  gl.filter.callrate(method = "loc", threshold = 0.9) %>%   
  #Remove individuals with >10% missing data
  gl.filter.callrate(method = "ind", threshold = 0.9) %>%
  #Remove maf
  gl.filter.maf(threshold = 0.05) %>%
  #HWE filter
  gl.filter.hwe()

#####Visualization
pc <- gl.pcoa(glMoose_clean)
moosePCA <- 
  gl.pcoa.plot(pc, glMoose_clean, labels = "legend", ellipse = T) +
  theme_classic()
moosePCA

library(plotly)
ggplotly()

#GENETIC diVersity
library(diveRsity)
basicStats(infile = 'diveRsity_subspecies.txt', outfile = 'test', fis_ci = T,
           ar_ci = T, fis_boots = 999, ar_boots = 999, 
           mc_reps = 9999, rarefaction = TRUE, ar_alpha = 0.05, 
           fis_alpha = 0.05)
############ANALYSES##############
#####FSTATS
library(StAMPP) #you may need to install the package
#Fst
MooseFst <-stamppFst(glMoose_clean, nboots=999, percent=95, nclusters=8)
MooseFst

Moose_NMT <- gl.recode.ind(glMoose_clean, "ind_recode.csv")
glMoose_clean$ind.names

Moose_sub <- gl.recode.pop(glMoose_clean, "poprecode.csv")
subFst <-stamppFst(Moose_sub, nboots=999, percent=95, nclusters=8)
subFst


#Input for Structure
gl2structure(glMoose_clean,outfile="Moose_Structure_1L.str", outpath = "C:/Users/Chase Smith/Desktop/Moose DArT")

###TESS3
moosegi <- gl2gi(glMoose_clean)
geno2 <- glMoose_clean
    geno2 <- glMoose_clean
    geno <- as.matrix(geno2)
    sample <- row.names(geno)
    pop.names <- pop(geno2)
    ploidy <- ploidy(geno2)
    geno = geno
    geno[is.na(geno)] = NaN
    format <- vector(length = length(geno[, 1]))
    format[1:length(geno[, 1])] = "genlight"
    pops <- unique(pop.names)
    pop.num <- vector(length = length(geno[, 1]))
    for (i in 1:length(geno[, 1])) {
      pop.num[i] = which(pop.names[i] == pops)
    }
    genoLHS <- as.data.frame(cbind(sample, pop.names, pop.num, 
                                   ploidy, format))
    geno <- cbind(genoLHS, geno)
    geno[, 2] = as.character(pop.names)
    geno[, 4] = as.numeric(as.character(geno[, 4]))
    row.names(geno) = NULL

coordinates <- moosegi$other$latlong
coordinates <- as.matrix(coordinates)
coordinates <- coordinates[,c(2,1)]
genotype <- geno[, -c(1:5)]

library(tess3r)
# Running the tess3 function
tess3.obj <- tess3(genotype, XProba = NULL, coordinates, K = 1:7, ploidy = 2, method = "projected.ls")

# Plot error
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# Retrieve the Q-matrix for K = 3 clusters
q.matrix2 <- qmatrix(tess3.obj, K = 2)
q.matrix3 <- qmatrix(tess3.obj, K = 3)

## STRUCTURE-like barplot for the Q-matrix
barplot(q.matrix2, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix", col.palette = CreatePalette(color.vector = c("purple", "blue"))) -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)

barplot(q.matrix3, border = NA, space = 0,
        xlab = "Individuals", ylab = "Ancestry proportions",
        main = "Ancestry matrix", col.palette = CreatePalette(color.vector = c("tomato", "orange", "lightblue"))) -> bp
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4)