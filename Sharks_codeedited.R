#Analysis of Shark Samples
#Among treatments what are the broad changes in diversity and richness of microbial communities and composition?

#Picante and vegan libraries used in this analysis
library(picante)
library(vegan)


#OTU picking file, converted from biom fromat as output by qiime, picking from the current silva database
OTUS_SILVA <- read.csv(file.choose(), header=T, row.names=1)


##############DATA FORMATTING#############################################################################
#transform the data
#transpose the data set
tdata2 = t(OTUS_SILVA)
#Format transposed data as a dataframe
tdata3 = as.data.frame(tdata2)
class(tdata3)

#Import meta data for the entire sample set
Real_Meta <- read.csv(file.choose(), header=T, row.names=1)
All.Real_Meta <- read.csv(file.choose(), header=T, row.names=1)

#Format meta data as a data frame
as.data.frame(Real_Meta)
#Use row names in tdata3 as are found in the Real_Meta file
tdata3<-tdata3[rownames(Real_Meta),]
#ensure that all row names in tdata3 and Real_Meta match
all.equal(rownames(tdata3),rownames(Real_Meta))


# check total abundance in each sample
apply(tdata3, 2, sum)

# Turn percent cover to relative abundance by dividing each value by sample total abundance

#Relative abundance
tdata3.ra <- decostand(tdata3, method = "total")

#Presence/Absence of new taxa
tdata3.pa <- decostand(tdata3, method = "pa")
# check total abundance in each sample
apply(tdata3.ra, 1, sum)


AllMeta<- read.csv(file.choose(), header=T, row.names=1)
#Use row names in tdata3 as are found in the Real_Meta file
tdata3W<-tdata3[rownames(AllMeta),]
#ensure that all row names in tdata3 and Real_Meta match
all.equal(rownames(tdata3W),rownames(AllMeta))


###########TEETH##############################################################################################

#Set up Teeth OTUS
#Import meta data file with only the meta data for the teeth samples
TeethMeta<- read.csv(file.choose(), header=T, row.names=1)

#create a subset of just the teeth sample IDs and the OTUS which coorelate 

Teeth.sample.ids <- intersect(rownames(TeethMeta), rownames(tdata3))
Teeth.sub <- tdata3[Teeth.sample.ids,]
meta.sub <- TeethMeta[Teeth.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(Teeth.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}



# Turn percent cover to relative abundance by dividing each value by sample total abundance
#Relative Abundance
teethdata2.ra <- decostand(Teeth.sub, method = "total")
#Presence/Absence of new taxa
teethdata2.pa <- decostand(Teeth.sub, method = "pa")
# check total abundance in each sample
apply(teethdata2.ra, 1, sum)



########ENVIRONMENTAL#####################################################################################

#Set up Water OTUS
#Import meta data file with only the meta data for the water samples
WaterMeta<- read.csv(file.choose(), header=T, row.names=1)

#create a subset of just the water sample IDs and the OTUS which coorelate 
Water.sample.ids <- intersect(rownames(WaterMeta), rownames(tdata3))
Water.sub <- tdata3[Water.sample.ids,]
meta.sub <- Real_Meta[Water.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(Water.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}


# Turn percent cover to relative abundance by dividing each value by sample total abundance
#Relative Abundance
Waterdata.ra <- decostand(Water.sub, method = "total")
#Presence/Absence of new taxa
Waterdata.pa <- decostand(Water.sub, method = "pa")
# check total abundance in each sample
apply(Waterdata.ra, 1, sum)


#Overlap between teeth and water
Water.otus = which(colSums(Waterdata.ra) > .01)
TWaterData = teethdata2.ra[,Water.otus]
rowSums(TWaterData)
#Gives mean percentage of overlap
mean(rowSums(TWaterData))

#Pull only OTUS specific to Teeth
Teeth.Spec<-setdiff(colnames(teethdata2.ra), colnames(TWaterData))
TeethDataSansWater<- teethdata2.ra[,Teeth.Spec]
mean(rowSums(TeethDataSansWater))

#Analysis of variance using the distance matrices (Adonis)
adonis(TeethDataSansWater~TeethMeta$Species)
adonis(TWaterData~TeethMeta$Species)
adonis(teethdata2.ra~TeethMeta$Species)

richnessW<-specnumber(tdata3W)
richnessWaov<-aov(richnessW~SampleType, data=AllMeta)
summary(richnessWaov)

boxplot(richnessW~AllMeta$SampleType, col=c("green", "blue"))


Diversity1W<-diversity(tdata3W, index="shannon")
watershanaov<-aov(Diversity1W~SampleType,data=AllMeta)  
summary(watershanaov)  

boxplot(Diversity1W~AllMeta$SampleType, col=c("green", "blue"))


Diversity2W<-diversity(tdata3W, index="invsimpson")
watersimpaov<-aov(Diversity2W~SampleType,data=AllMeta)  
summary(watersimpaov)  

boxplot(Diversity2W~AllMeta$SampleType, col=c("green", "blue"))

tdata3W.ra <- decostand(tdata3W, method = "total")
tdata3W.pa <- decostand(tdata3W, method = "pa")
# calculate Bray-Curtis distance among samples
tdata3W.bc.dist <- vegdist(tdata3W.ra, method = "bray")
tdata3W.pa.bc.dist <- vegdist(tdata3W.pa, method = "bray")

# cluster communities using average-linkage algorithm
tdata3W.bc.clust <- hclust(tdata3W.bc.dist, method = "average")
tdata3W.pa.bc.clust <- hclust(tdata3W.pa.bc.dist, method = "average")


SharkWater_NMDS=metaMDS(tdata3W.bc.dist,k=2)
stressplot(SharkWater_NMDS)
plot(SharkWater_NMDS)
ordiplot(SharkWater_NMDS, display = "sites", type = "text")

# ordination plots are highly customizable set up the plotting area but
# don't plot anything yet
mds.fig <- ordiplot(SharkWater_NMDS, type = "none", legend("bottomright",legend=as.character(paste(" ",unique(AllMeta$SampleType))), 
                                                           pch=19, col=c("green", "blue", "red", "blue", "black"), cex = 0.8))
# plot just the samples, colour by species, pch=19 means plot a circle
points(mds.fig, "sites", pch = 19, col = "green", select = AllMeta$SampleType=="Shark")
points(mds.fig, "sites", pch = 19, col = "blue", select = AllMeta$SampleType=="Water")

#Check the species richness to determine if any samples need to be excluded
richnessk<-specnumber(TeethDataSansWater)
richnessk
richnessA<-specnumber(TWaterData)
richnessA
richnessB<-specnumber(teethdata2.ra)
richnessB
##################SHARKS(ALL)#################################################################################

#Set up Shark OTUS
#Import meta data file with only the meta data for the shark samples, excluding those with low sample size 
#and species richness
SharkMeta<- read.csv(file.choose(), header=T, row.names=1)
#create a subset of just the Shark sample IDs and the OTUS which coorelate 
Shark.sample.ids <- intersect(rownames(All.Real_Meta), rownames(tdata3))
Shark.sub <- tdata3[Shark.sample.ids,]
Sharkmeta.sub <- All.Real_Meta[Shark.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(Shark.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}


# Turn percent cover to relative abundance by dividing each value by sample total abundance
#Relative Abundance
Sharkdata.ra <- decostand(Shark.sub, method = "total")
#Presence/Absence of new taxa
Sharkdata.pa <- decostand(Shark.sub, method = "pa")

#Overlap between Environmental and Shark Samples
Water.otus = which(colSums(Waterdata.ra) > .01)
Shark.otus = which(colSums(Sharkdata.ra) > .01)

WaterOnlyData = Shark.sub[,Water.otus]

Shark.Spec<-setdiff(colnames(Shark.sub), colnames(WaterOnlyData))
SharkDataSansWater<- Sharkdata.ra[,Shark.Spec]

rowSums(SharkDataSansWater)
#Gives mean percentage of overlap
mean(rowSums(SharkDataSansWater))

#Check the species richness of the samples without the overlap of the water
RichnessShark<-specnumber(SharkDataSansWater)
RichnessShark
mean(RichnessShark)
sd(RichnessShark)

SharkWaterMeta<- read.csv(file.choose(), header=T, row.names=1)



adonis(tdata3W~AllMeta$SampleType)
adonis(SharkDataSansWater~All.Real_Meta$Gender , select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~Gender*Location, data=All.Real_Meta, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~Gender*Species, data=All.Real_Meta, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~Gender*Month, data=All.Real_Meta, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~Gender*Habitat, data=All.Real_Meta, select=All.Real_Meta$SampleType=="Shark")

adonis(SharkDataSansWater~All.Real_Meta$Location, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Location*Gender, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Location*Species, data=All.Real_Meta, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Location*Month, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Location*Habitat, select=All.Real_Meta$SampleType=="Shark")

adonis(SharkDataSansWater~All.Real_Meta$Species, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Species*Gender, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Species*Location, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Species*Month, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Species*Habitat, select=All.Real_Meta$SampleType=="Shark")


adonis(SharkDataSansWater~All.Real_Meta$Month, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Month*Gender, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Month*Location, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Month*Species, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Month*Habitat, select=All.Real_Meta$SampleType=="Shark")

adonis(SharkDataSansWater~All.Real_Meta$Habitat, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Habitat*Gender, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Habitat*Location, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Habitat*Species, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~All.Real_Meta$Habitat*Month, select=All.Real_Meta$SampleType=="Shark")


adonis(SharkDataSansWater~Gender*Species*Location*Month, data= All.Real_Meta, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~Gender*Species*Location*Gender, data= All.Real_Meta, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~Gender*Species*Location*Habitat, data= All.Real_Meta, select=All.Real_Meta$SampleType=="Shark")

adonis(SharkDataSansWater~Gender*Species*Location*Month*Habitat, data= All.Real_Meta, select=All.Real_Meta$SampleType=="Shark")
adonis(SharkDataSansWater~Gender+Species+Location+Month+Habitat, data= All.Real_Meta, select=All.Real_Meta$SampleType=="Shark")

#Rarefraction curves to determine if samples need to be excluded
rarecurve(Shark.sub,label = FALSE)

###########Phylogenetic Diversity################################################################################


phy <- read.tree(file.choose())
class(phy)


### root
print("tree was rooted to the longest edge before tips were removed")
longest.edge <- which.max(phy$edge.length)
longest.edge
long.nodes <- phy$edge[longest.edge,]
long.nodes
new.outgroup <- long.nodes[long.nodes < Ntip(phy)]
new.outgroup
phy.rooted <- root(phy, outgroup=new.outgroup, resolve.root=T)   ### one rooting method
phy.rooted


common.samps.rows<-intersect(rownames(Sharkdata.ra),rownames(tdata3))
common.samps.cols<-intersect(colnames(SharkDataSansWater),colnames(tdata3))
tdata3.mod<-tdata3.ra[common.samps.rows,common.samps.cols]


# Turn percent cover to relative abundance by dividing each value by sample
# total abundance
tdata3.mod.ra <- decostand(tdata3.mod, method = "total")
tdata3.mod.pa <- decostand(tdata3.mod, method = "pa")

#Remove singletons
SharksMinusSing<- which(colSums(tdata3.mod.pa)>1)
SharksMinusSing
Shark.phy<- tdata3.mod[,SharksMinusSing]
Shark.phy
Shark.phy.ra <- decostand(tdata3.mod, method = "total")

combined <- match.phylo.comm(phy.rooted, Shark.phy)

phy2 <- combined$phy
comm <- combined$comm

#######################WHOLE DATA SET STATS##########################################################################
# Turn percent cover to relative abundance by dividing each value by sample
# total abundance
tdata3.ra <- decostand(tdata3, method = "total")
tdata3.pa <- decostand(tdata3, method = "pa")
# calculate Bray-Curtis distance among samples
tdata3.bc.dist <- vegdist(tdata3.ra, method = "bray")
tdata3.pa.bc.dist <- vegdist(tdata3.pa, method = "bray")

# cluster communities using average-linkage algorithm
tdata3.bc.clust <- hclust(tdata3.bc.dist, method = "average")
tdata3.pa.bc.clust <- hclust(tdata3.pa.bc.dist, method = "average")

# plot cluster diagram
plot(tdata3.bc.clust, ylab = "Bray-Curtis dissimilarity")
plot(tdata3.pa.bc.clust, ylab = "Bray-Curtis dissimilarity")




# calculate Bray-Curtis distance among samples
tdata3.mod.bc.dist <- vegdist(Shark.phy, method = "bray")
tdata3.mod.pa.bc.dist <- vegdist(Shark.phy.pa, method = "bray")

# cluster communities using average-linkage algorithm
tdata3.mod.bc.clust <- hclust(tdata3.mod.bc.dist, method = "average")
tdata3.mod.pa.bc.clust <- hclust(tdata3.mod.pa.bc.dist, method = "average")




# plot cluster diagram
plot(tdata3.mod.bc.clust, ylab = "Bray-Curtis dissimilarity")
plot(tdata3.mod.pa.bc.clust, ylab = "Bray-Curtis dissimilarity")


###Faith's
#Determine faiths phylogenetic diversity
Sharks.pd <- pd(comm, phy2)
head(Sharks.pd)
#Plot faiths PD
boxplot(Sharks.pd$PD ~ SharkMeta$Species, xlab = "Species", ylab = "Faith's PD")
boxplot(Sharks.pd$PD ~ SharkMeta$Location, xlab = "Location", ylab = "Faith's PD")
plot(Sharks.pd$PD ~ Sharks.pd$SR, xlab = "Species richness", ylab = "Faith's PD")

# convert phylogenety to a distance matrix
phy.dist <- cophenetic(phy2)

comm.bc.dist <- vegdist(comm, method = "bray")

# calculate phylogenetic Mean Nearest Taxon Distance beta diversity
comm.mntd.dist <- comdistnt(comm, phy.dist, abundance.weighted = TRUE)

# calculate Mantel correlation for taxonomic Bray-Curtis vs. phylogenetic MNTD diversity
mantel(comm.bc.dist, comm.mntd.dist)

#Mean MNTD???#################





# ANOVA
# richness
richness<-specnumber(Shark.phy)

data.aov<-aov(richness~Species,data=SharkMeta)
summary(data.aov)
mean(richness)
sd(richness)
boxplot(richness~All.Real_Meta$Species, col=c("gray", "red", "green", "blue", "cyan"))
boxplot(richness~All.Real_Meta$Location, col=c("gray", "red", "green", "blue", "cyan"))
boxplot(richness~Species+Location, data=SharkMeta, col=c("gray", "red", "green", "blue", "cyan"), las=2)



data.aov2<-aov(richness~Location,data=SharkMeta)
summary(data.aov2)


data.aov3<-aov(richness~Species*Location,data=SharkMeta)
summary(data.aov3)

richness2<-specnumber(tdata3)
data.aov4<-aov(richness2~Type,data=Real_Meta)
summary(data.aov4)

# diversity
#Shannon
divers<-diversity(Shark.phy,index="shannon")
tdata.aov<-aov(divers~Species,data=SharkMeta)
summary(tdata.aov)

boxplot(divers~SharkMeta$Species, col=c("gray", "red", "green", "blue", "cyan"))
boxplot(divers~SharkMeta$Location, col=c("gray", "red", "green", "blue", "cyan"))
boxplot(divers~Species+Location, data=SharkMeta, col=c("gray", "red", "green", "blue", "cyan"), las=2)



tdata.aov2<-aov(divers~Location,data=SharkMeta)
summary(tdata.aov2)

tdata.aov3<-aov(divers~Species*Location,data=SharkMeta)
summary(tdata.aov3)

divers2<-diversity(tdata3,index="shannon")
tdata.aov4<-aov(divers2~Type,data=Real_Meta)
summary(tdata.aov4)


#Inverse Simpson
diversi<-diversity(Shark.phy,index="invsimpson")
t1data.aov<-aov(diversi~Species,data=SharkMeta)
summary(t1data.aov)

boxplot(diversi~SharkMeta$Species, col=c("gray", "red", "green", "blue", "cyan"))
boxplot(diversi~SharkMeta$Location, col=c("gray", "red", "green", "blue", "cyan"))
boxplot(diversi~Species+Location, data=SharkMeta, col=c("gray", "red", "green", "blue", "cyan"), las=2)


t2data.aov2<-aov(diversi~Location,data=SharkMeta)
summary(t2data.aov2)

t3data.aov3<-aov(diversi~Species*Location,data=SharkMeta)
summary(t3data.aov3)

diversi2<-diversity(tdata3,index="invsimpson")
t1data.aov4<-aov(diversi2~Type,data=Real_Meta)
summary(t1data.aov4)

# Tukey
#Validate (post-hoc) analyses of varince
TukeyHSD(data.aov)
TukeyHSD(data.aov2)
TukeyHSD(data.aov3)

TukeyHSD(tdata.aov)
TukeyHSD(tdata.aov2)
TukeyHSD(tdata.aov3)

TukeyHSD(t1data.aov)
TukeyHSD(t2data.aov2)
TukeyHSD(t3data.aov3)

# Taxonomic (Bray-Curtis) dissimilarity explained
adonis(tdata3.mod.bc.dist~Species*Location, data=SharkMeta)
adonis(tdata3.mod.bc.dist~Location, data=SharkMeta)
adonis(tdata3.mod.bc.dist~Species, data=SharkMeta)
adonis(tdata3.bc.dist~Type, data=Real_Meta)

Shark_NMDS=metaMDS(tdata3.mod.bc.dist,k=2)
stressplot(Shark_NMDS)
plot(Shark_NMDS)
ordiplot(Shark_NMDS, display = "sites", type = "text")

# ordination plots are highly customizable set up the plotting area but
# don't plot anything yet
mds.fig <- ordiplot(Shark_NMDS, type = "none")
# plot just the samples, colour by species, pch=19 means plot a circle
points(mds.fig, "sites", pch = 19, col = "green", select = SharkMeta$Species+Location=="Nurse")
points(mds.fig, "sites", pch = 19, col = "cyan", select = SharkMeta$Species=="Tiger")
points(mds.fig, "sites", pch = 19, col = "red", select = SharkMeta$Species=="Lemon")
points(mds.fig, "sites", pch = 19, col = "blue", select = SharkMeta$Species=="Sandbar")
points(mds.fig, "sites", pch = 19, col = "black", select = SharkMeta$Species=="Caribbean Reef")
points(mds.fig, "sites", pch = 19, select = SharkMeta$Location=="Teeth")

Sharkall_NMDS=metaMDS(tdata3.mod.bc.dist,k=2)
stressplot(Sharkall_NMDS)
plot(Sharkall_NMDS)
ordiplot(Sharkall_NMDS, display = "sites", type = "text")

# ordination plots are highly customizable set up the plotting area but
# don't plot anything yet
mds.fig <- ordiplot(Sharkall_NMDS, type = "none")
# plot just the samples, colour by species, pch=19 means plot a circle
points(mds.fig, "sites", pch = 19, col = "green", select = SharkMeta$Species=="Nurse")
points(mds.fig, "sites", pch = 19, col = "cyan", select = SharkMeta$Species=="Tiger")
points(mds.fig, "sites", pch = 19, col = "red", select = SharkMeta$Species=="Lemon")
points(mds.fig, "sites", pch = 19, col = "blue", select = SharkMeta$Species=="Sandbar")
points(mds.fig, "sites", pch = 19, col = "black", select = SharkMeta$Species=="Caribbean Reef")
points(mds.fig, "sites", pch = 0, col = "black", select = SharkMeta$Location=="Gills")
points(mds.fig, "sites", pch = 2, col = "black", select = SharkMeta$Location=="Teeth")
points(mds.fig, "sites", pch = 5, col = "black", select = SharkMeta$Location=="Skin")
points(mds.fig, "sites", pch = 8, col = "black", select = SharkMeta$Location=="Cloaca")


# ordination plots are highly customizable set up the plotting area but
# don't plot anything yet
mds.fig3 <- ordiplot(Shark_NMDS, type = "none", legend("bottomright",legend=as.character(paste(" ",unique(SharkMeta$Location))), 
                                                       pch=19, col=c("green", "cyan", "red", "blue"), cex = 0.8))

# plot just the samples, colour by species, pch=19 means plot a circle
points(mds.fig3, "sites", pch = 19, col = "green", select = SharkMeta$Location=="Gills")
points(mds.fig3, "sites", pch = 19, col = "cyan", select = SharkMeta$Location=="Teeth")
points(mds.fig3, "sites", pch = 19, col = "red", select = SharkMeta$Location=="Skin")
points(mds.fig3, "sites", pch = 19, col = "blue", select = SharkMeta$Location=="Cloaca")



adonis(tdata3.mod.pa.bc.dist~Species*Location, data=SharkMeta)
adonis(tdata3.mod.pa.bc.dist~Species, data=SharkMeta)
adonis(tdata3.mod.pa.bc.dist~Location, data=SharkMeta)
adonis(tdata3.pa.bc.dist~Type, data=Real_Meta)
library(RVAideMemoire)
pairwise.perm.manova(Shark.phy, SharkMeta$Species, nperm=20)


#Test beta dispersion of samples based on relative abundance and BC distance
sharkbeta<-betadisper(tdata3.mod.bc.dist, SharkMeta$Species, type = c("median","centroid"), bias.adjust = FALSE,
                      sqrt.dist = FALSE, add = FALSE)
sharkbeta

anova(sharkbeta)
permutest(sharkbeta)

labs <- paste("Dimension", 1:4, "(", 
              round(100*sharkbeta$eig / sum(sharkbeta$eig), 2), "%)")

plot(sharkbeta, cex=1, pch=15:20,
     main="Shark Species: MDS coordinates", cex.lab=1.25,
     xlab=labs[1], ylab=labs[2],
     hull=FALSE, ellipse=TRUE, conf=.68, lwd=-1)

boxplot(sharkbeta, xlab="Species", notch=FALSE, col=c("gray", "red", "green", "blue", "cyan"))

#Test beta dispersion of samples based on Presence/Absence of new taxa and BC distance
sharkbetapa<-betadisper(tdata3.mod.pa.bc.dist, SharkMeta$Species, type = c("median","centroid"), bias.adjust = FALSE,
                        sqrt.dist = FALSE, add = FALSE)
sharkbetapa

anova(sharkbetapa)
permutest(sharkbetapa)

labspa <- paste("Dimension", 1:4, "(", round(100*sharkbetapa$eig / sum(sharkbetapa$eig), 2), "%)")

plot(sharkbetapa, cex=1, pch=15:20,
     main="Shark Species: MDS coordinates", cex.lab=1.25,
     xlab=labspa[1], ylab=labspa[2],
     hull=FALSE, ellipse=TRUE, conf=.68, lwd=-1)

boxplot(sharkbetapa, xlab="Species", notch=FALSE, col=c("gray", "red", "green", "blue", "cyan"))



#NMDS analysis: species
Shark_NMDS=metaMDS(tdata3.mod,k=2)
stressplot(Shark_NMDS)
dev.new()
plot(Shark_NMDS)
ordiplot(Shark_NMDS, display = "sites", type = "text")

# ordination plots are highly customizable set up the plotting area but
# don't plot anything yet
mds.fig <- ordiplot(Shark_NMDS, type = "none")
# plot just the samples, colour by species, pch=19 means plot a circle
points(mds.fig, "sites", pch = 19, col = "green", select = SharkMeta$Species=="Nurse")
points(mds.fig, "sites", pch = 19, col = "cyan", select = SharkMeta$Species=="Tiger")
points(mds.fig, "sites", pch = 19, col = "red", select = SharkMeta$Species=="Lemon")
points(mds.fig, "sites", pch = 19, col = "blue", select = SharkMeta$Species=="Sandbar")
points(mds.fig, "sites", pch = 19, col = "black", select = SharkMeta$Species=="Caribbean Reef")


##Plot PA matrix 
Shark_NMDS.pa=metaMDS(tdata3.mod.pa,k=2)
stressplot(Shark_NMDS.pa)
plot(Shark_NMDS.pa)
ordiplot(Shark_NMDS.pa, display = "sites", type = "text")

# ordination plots are highly customizable set up the plotting area but
# don't plot anything yet
mds.fig.pa <- ordiplot(Shark_NMDS.pa, type = "none")
# plot just the samples, colour by species, pch=19 means plot a circle
points(mds.fig.pa, "sites", pch = 19, col = "green", select = SharkMeta$Species=="Nurse")
points(mds.fig.pa, "sites", pch = 19, col = "cyan", select = SharkMeta$Species=="Tiger")
points(mds.fig.pa, "sites", pch = 19, col = "red", select = SharkMeta$Species=="Lemon")
points(mds.fig.pa, "sites", pch = 19, col = "blue", select = SharkMeta$Species=="Sandbar")
points(mds.fig.pa, "sites", pch = 19, col = "black", select = SharkMeta$Species=="Caribbean Reef")

#NMDS analysis: Location
#Relative Abundance
# ordination plots are highly customizable set up the plotting area but
# don't plot anything yet
dev.new()
mds.fig2 <- ordiplot(Shark_NMDS, type = "none")
# plot just the samples, colour by species, pch=19 means plot a circle
points(mds.fig2, "sites", pch = 19, col = "green", select = SharkMeta$Location=="Gills")
points(mds.fig2, "sites", pch = 19, col = "cyan", select = SharkMeta$Location=="Teeth")
points(mds.fig2, "sites", pch = 19, col = "red", select = SharkMeta$Location=="Skin")
points(mds.fig2, "sites", pch = 19, col = "blue", select = SharkMeta$Location=="Cloaca")

#Presence/Absence
mds.fig2pa <- ordiplot(Shark_NMDS.pa, type = "none")
# plot just the samples, colour by species, pch=19 means plot a circle
points(mds.fig2pa, "sites", pch = 19, col = "green", select = SharkMeta$Location=="Gills")
points(mds.fig2pa, "sites", pch = 19, col = "cyan", select = SharkMeta$Location=="Teeth")
points(mds.fig2pa, "sites", pch = 19, col = "red", select = SharkMeta$Location=="Skin")
points(mds.fig2pa, "sites", pch = 19, col = "blue", select = SharkMeta$Location=="Cloaca")



adonis(SharkDataSansWater~Species+Location, data=SharkMeta)


# simper
sink("simper.SharkInteraction.csv")
SharkWatersimp<-simper(Shark.phy, Real_Meta$Species+Location, permutations = 499)
summary(SharkWatersimp)
sink()

sink("simper.spec.csv")
specsimp<-simper(Shark.phy, SharkMeta$Species, permutations = 499)
summary(specsimp)
sink()

sink("simper.location.csv")
locsimp<-simper(Shark.phy, SharkMeta$Location, permutations = 499)
summary(locsimp)
sink()





##########RA Curves




library(goeveg)
Tig_Meta <- read.csv(file.choose(), header=T, row.names=1)
as.data.frame(Tig_Meta)
Tigdata3<-Shark.phy[rownames(Tig_Meta),]
sorted<-Tigdata3[order(colSums(-Tigdata3))]
TMinusSing<-sorted[1:20]
racurve(TMinusSing, main = "Rank Abundance:Tiger Sharks", nlab = 5, ylog = FALSE,
        frequency = FALSE)
Tigcolm<-colMeans(sorted)
Tigcolm
Tigsdmatrix<-colSds(as.matrix(sorted))
Tigsdmatrix



Nur_Meta <- read.csv(file.choose(), header=T, row.names=1)
as.data.frame(Nur_Meta)
Nurdata3<-Shark.phy[rownames(Nur_Meta),]
Nsorted<-Nurdata3[order(colSums(-Nurdata3))]
NMinusSing<-Nsorted[1:20]
racurve(NMinusSing, main = "Rank Abundance:Nurse Sharks", nlab = 3, ylog = FALSE,
        frequency = FALSE)
Nurcolm<-colMeans(Nsorted)
Nurcolm
Nursdmatrix<-colSds(as.matrix(Nsorted))
Nursdmatrix



Lem_Meta <- read.csv(file.choose(), header=T, row.names=1)
as.data.frame(Lem_Meta)
Lemdata3<-Shark.phy[rownames(Lem_Meta),]
Lsorted<-Lemdata3[order(colSums(-Lemdata3))]
LMinusSing<-Lsorted[1:20]
racurve(LMinusSing, main = "Rank Abundance:Lemon Sharks", nlab = 3, ylog = FALSE,
        frequency = FALSE)
Lemcolm<-colMeans(Lsorted)
Lemcolm
Lemsdmatrix<-colSds(as.matrix(Lsorted))
Lemsdmatrix

OTUtaxa_Meta <- read.csv(file.choose(), header=T, row.names=1)


SB_Meta <- read.csv(file.choose(), header=T, row.names=1)
as.data.frame(SB_Meta)
SBdata3<-Shark.phy[rownames(SB_Meta),]
SBsorted<-SBdata3[order(colSums(-SBdata3))]
SBMinusSing<-SBsorted[1:20]
racurve(SBMinusSing, main = "Rank Abundance:Sandbar Sharks", nlab = 3, ylog = FALSE,
        frequency = FALSE)
Sandcolm<-colMeans(SBsorted)
Sandcolm
Sandsdmatrix<-colSds(as.matrix(SBsorted))
Sandsdmatrix



CR_Meta <- read.csv(file.choose(), header=T, row.names=1)
as.data.frame(CR_Meta)
CRdata3<-Shark.phy[rownames(CR_Meta),]
CRsorted<-CRdata3[order(colSums(-CRdata3))]
CRMinusSing<-CRsorted[1:20]
racurve(CRMinusSing, main = "Rank Abundance:Caribbean Reef Sharks", nlab = 3, ylog = FALSE,
        frequency = FALSE)
Carcolm<-colMeans(CRsorted)
Carcolm
Carsdmatrix<-colSds(as.matrix(CRsorted))
Carsdmatrix



#NMDS Relative Abundance
Sharkwater_NMDS=metaMDS(teethwaterdata3,k=3)
dev.new()
ordiplot(Sharkwater_NMDS, display = "sites", type = "text")
mds.fig22<- ordiplot(Sharkwater_NMDS, type = "none")
# plot just the samples, colour by type, pch=19 means plot a circle
points(mds.fig22, "sites", pch = 19, col = "black", select = SharkWaterMeta$SampleType=="Shark")
points(mds.fig22, "sites", pch = 19, col = "blue", select = SharkWaterMeta$SampleType=="Environment")


dev.new()
mds.fig <- ordiplot(Teethwater_NMDS, type = "none")
# plot just the samples, colour by type, pch=19 means plot a circle
points(mds.fig, "sites", pch = 19, col = "black", select = TeethWaterMeta$Location=="Teeth")
points(mds.fig, "sites", pch = 19, col = "blue", select = TeethWaterMeta$Location=="Water")
