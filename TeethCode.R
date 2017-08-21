#Analysis of Shark Teeth Samples
#Among treatments what are the broad changes in diversity and richness of microbial communities and composition?
library(picante)
library(vegan)
OTUS_SILVA <- read.csv(file.choose(), header=T, row.names=1)


##############DATA FORMATTING#############################################################################

#Entire data set
#transform the data
#transpose the data set
tdata2 = t(OTUS_SILVA)
#Format transposed data as a dataframe
tdata3 = as.data.frame(tdata2)
class(tdata3)

#Import meta data for the entire sample set
Real_Meta <- read.csv(file.choose(), header=T, row.names=1)

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


#Teeth only
#transform the data
teethdata2 = t(OTUS_SILVA)
teethdata3 = as.data.frame(teethdata2)
class(teethdata3)

TeethMeta <- read.csv(file.choose(), header=T, row.names=1)

as.data.frame(TeethMeta)
teethdata3<-teethdata3[rownames(TeethMeta),]
all.equal(rownames(teethdata3),rownames(TeethMeta))


#Check to see if any OTUs are found ONLY in teeth and nowhere else
Teeth.Spec<-setdiff(colnames(teethdata2.ra), colnames(tdata3))
TeethDataSansAll<- teethdata2.ra[,Teeth.Spec]
mean(rowSums(TeethDataSansAll))

# check total abundance in each sample
apply(teethdata3, 1, sum)

# Turn percent cover to relative abundance by dividing each value by sample
# total abundance
teethdata3.ra <- decostand(teethdata3, method = "total")

teethdata3.pa <- decostand(teethdata3, method = "pa")
# check total abundance in each sample
apply(teethdata3.ra, 1, sum)


########ENVIRONMENTAL#####################################################################################

#Set up Water OTUS

WaterMeta<- read.csv(file.choose(), header=T, row.names=1)

TWater.sample.ids <- intersect(rownames(WaterMeta), rownames(tdata3))
TWater.sub <- tdata3[TWater.sample.ids,]
Tmeta.sub <- Real_Meta[TWater.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(TWater.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}


# Turn percent cover to relative abundance by dividing each value by sample
# total abundance
TWaterdata.ra <- decostand(TWater.sub, method = "total")

TWaterdata.pa <- decostand(TWater.sub, method = "pa")

###########TEETH##############################################################################################

#Set up Teeth OTUS
Teeth.sample.ids <- intersect(rownames(TeethMeta), rownames(teethdata3))
Teeth.sub <- teethdata3[Teeth.sample.ids,]
Tmeta.sub <- TeethMeta[Teeth.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(Teeth.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}


# Turn percent cover to relative abundance by dividing each value by sample
# total abundance
teethdata2.ra <- decostand(Teeth.sub, method = "total")

teethdata2.pa <- decostand(Teeth.sub, method = "pa")

# check total abundance in each sample
apply(teethdata2.ra, 1, sum)

#Overlap between teeth and water
Water.otus = which(colSums(TWaterdata.ra) > .01)
TWaterData = teethdata2.ra[,Water.otus]
rowSums(TWaterData)
#Gives mean percentage of overlap
mean(rowSums(TWaterData))
 

#Pull only OTUS specific to Teeth
Teeth.Spec<-setdiff(colnames(teethdata2.ra), colnames(TWaterData))
TeethDataSansWater<- teethdata2.ra[,Teeth.Spec]
mean(rowSums(TeethDataSansWater))

adonis(TeethDataSansWater~TeethMeta$Species)
adonis(TWaterData~TeethMeta$Species)
adonis(teethdata3~TeethMeta$Species)


richness<-specnumber(TeethDataSansWater)
richness
richnessA<-specnumber(TWaterData)
richnessA
richnessB<-specnumber(teethdata3)
richnessB


###########Phylogenetic Diversity###############################################################################


teethphy <- read.tree(file.choose())
class(teethphy)


### root
print("tree was rooted to the longest edge before tips were removed")
longest.edge <- which.max(teethphy$edge.length)
longest.edge
long.nodes <- teethphy$edge[longest.edge,]
long.nodes
new.outgroup <- long.nodes[long.nodes < Ntip(teethphy)]
new.outgroup
Teeth.phy <- root(teethphy, outgroup=new.outgroup, resolve.root=T)   ### one rooting method
Teeth.phy.pa<-decostand(Teeth.phy, method = "pa")

teeth.common.samps.rows<-intersect(rownames(teethdata3.ra),rownames(teethdata3))
teeth.common.samps.cols<-intersect(colnames(TeethDataSansWater),colnames(teethdata3))
teethdata3.mod<-teethdata3.ra[teeth.common.samps.rows,teeth.common.samps.cols]



# Turn percent cover to relative abundance by dividing each value by sample
# total abundance
teethdata3.mod.ra <- decostand(teethdata3.mod, method = "total")
teethdata3.mod.pa <- decostand(teethdata3.mod, method = "pa")

#Remove singletons
teethMinusSing<- which(colSums(teethdata3.mod)>.1)
teethMinusSing
teeth.phy<- teethdata3.mod[,teethMinusSing]
teeth.phy
teeth.phy.pa <- decostand(teethdata3.mod, method = "pa")

combined <- match.phylo.comm(phy.rooted, Shark.phy)

phy2 <- combined$phy
comm <- combined$comm


###Faith's

Teeth.pd <- pd(comm, Teeth.phy)
head(Teeth.pd)


teeth.phy.dist <- cophenetic(Teeth.phy)


# calculate phylogenetic MNTD beta diversity
comm.mntd.dist <- comdistnt(comm, teeth.phy.dist, abundance.weighted = TRUE)

# calculate Mantel correlation for taxonomic Bray-Curtis vs. phylogenetic
# MNTD diversity
mantel(comm.bc.dist, comm.mntd.dist)

# calculate Mantel correlation for taxonomic Bray-Curtis vs. trait MNTD
# diversity
mantel(comm.bc.dist, comm.mntd.traits.dist)



#############TEETH STATS########################################################################################

# calculate Bray-Curtis distance among samples
teethdata3.bc.dist <- vegdist(teethdata3.mod.ra, method = "bray")
teethdata3.pa.bc.dist <- vegdist(teethdata3.mod.pa, method = "bray")

# cluster communities using average-linkage algorithm
teethdata3.bc.clust <- hclust(teethdata3.bc.dist, method = "average")
teethdata3.pa.bc.clust <- hclust(teethdata3.pa.bc.dist, method = "average")




# plot cluster diagram
plot(teethdata3.bc.clust, ylab = "Bray-Curtis dissimilarity")
plot(teethdata3.pa.bc.clust, ylab = "Bray-Curtis dissimilarity")


# ANOVA
# richness
TeethRichness<-specnumber(TeethDataSansWater)
Teeth.aov<-aov(TeethRichness~Species,data=TeethMeta)
summary(Teeth.aov)  
boxplot(TeethRichness~TeethMeta$Species, col=c("gray", "red", "green", "blue", "cyan"), las=2)




# diversity
Teethdivers<-diversity(TeethDataSansWater,index="shannon")
Teeth.aov2<-aov(Teethdivers~Species,data=TeethMeta)
summary(Teeth.aov2)
boxplot(Teethdivers~TeethMeta$Species, col=c("gray", "red", "green", "blue", "cyan"), las=2)



Teethdiversi<-diversity(TeethDataSansWater,index="invsimpson")
Teeth.aov3<-aov(Teethdiversi~Species,data=TeethMeta)
summary(Teeth.aov3)
boxplot(Teethdiversi~TeethMeta$Species, col=c("gray", "red", "green", "blue", "cyan"), las=2)


# Tukey
TukeyHSD(Teeth.aov)
TukeyHSD(Teeth.aov2)
TukeyHSD(Teeth.aov3)
#Changing dominant structure across respective treatments

# Taxonomic (Bray-Curtis) dissimilarity explained
adonis(teethdata3.bc.dist~Species, data=TeethMeta)
adonis(teethdata3.pa.bc.dist~Species, data=TeethMeta)

library(RVAideMemoire)
pairwise.perm.manova(teeth.phy, TeethMeta$Species, nperm=49)




# simper
sink("simper.teeth.csv")
teethsimp<-simper(TeethDataSansWater, TeethMeta$Species, permutations = 499)
summary(teethsimp)
sink()

#############TEETH TO WATER########################################


TeethWaterMeta <- read.csv(file.choose(), header=T, row.names=1)

as.data.frame(TeethWaterMeta)
teethwaterdata3<-tdata3[rownames(TeethWaterMeta),]
all.equal(rownames(teethwaterdata3),rownames(TeethWaterMeta))
# check total abundance in each sample
apply(teethwaterdata3, 1, sum)

# Turn percent cover to relative abundance by dividing each value by sample
# total abundance
teethwaterdata3.ra <- decostand(teethwaterdata3, method = "total")

teethwaterdata3.pa <- decostand(teethwaterdata3, method = "pa")
# check total abundance in each sample
apply(teethwaterdata3.ra, 1, sum)



# calculate Bray-Curtis distance among samples
teethwaterdata3.bc.dist <- vegdist(teethwaterdata3.ra, method = "bray")
teethwaterdata3.pa.bc.dist <- vegdist(teethwaterdata3.pa, method = "bray")

# cluster communities using average-linkage algorithm
teethwaterdata3.bc.clust <- hclust(teethwaterdata3.bc.dist, method = "average")
teethwaterdata3.pa.bc.clust <- hclust(teethwaterdata3.pa.bc.dist, method = "average")




# plot cluster diagram
plot(teethwaterdata3.bc.clust, ylab = "Bray-Curtis dissimilarity")
plot(teethwaterdata3.pa.bc.clust, ylab = "Bray-Curtis dissimilarity")


# ANOVA
# richness
TeethWaterRichness<-specnumber(teethwaterdata3)
TeethWater.aov<-aov(TeethWaterRichness~Species,data=TeethWaterMeta)
summary(TeethWater.aov)  


# diversity
TeethWaterdivers<-diversity(teethwaterdata3,index="shannon")
TeethWater.aov2<-aov(Teethwaterdivers~Species,data=TeethWaterMeta)
summary(TeethWater.aov2)



TeethWaterdiversi<-diversity(teethwaterdata3,index="invsimpson")
TeethWater.aov3<-aov(TeethWaterdiversi~Species,data=TeethWaterMeta)
summary(TeethWater.aov3)


# Tukey
TukeyHSD(TeethWater.aov)
TukeyHSD(TeethWater.aov2)
TukeyHSD(TeethWater.aov3)
#Changing dominant structure across respective treatments

# Taxonomic (Bray-Curtis) dissimilarity explained
adonis(teethwaterdata3.bc.dist~Species, data=TeethWaterMeta)
adonis(teethwaterdata3.pa.bc.dist~Species, data=TeethWaterMeta)
adonis(teethwaterdata3.bc.dist~Species*Location, data=TeethWaterMeta)




sink("simper.TeethWater.csv")
TeethWatersimp<-simper(teethwaterdata3.ra, TeethWaterMeta$Species, permutations = 499)
summary(TeethWatersimp)
sink()



##################################Data visualization####################
######DATA VISUALIZATION######
Teeth_NMDS=metaMDS(teethdata3,k=2)
stressplot(Teeth_NMDS)
plot(Teeth_NMDS)
ordiplot(Teeth_NMDS, display = "sites", type = "text")

# ordination plots are highly customizable set up the plotting area but
# don't plot anything yet

dev.new()

mds.fig <- ordiplot(Teeth_NMDS, type = "none")
# plot just the samples, colour by species, pch=19 means plot a circle
points(mds.fig, "sites", pch = 19, col = "green", select = TeethMeta$Species=="Nurse")
points(mds.fig, "sites", pch = 19, col = "cyan", select = TeethMeta$Species=="Tiger")
points(mds.fig, "sites", pch = 19, col = "red", select = TeethMeta$Species=="Lemon")
points(mds.fig, "sites", pch = 19, col = "blue", select = TeethMeta$Species=="Sandbar")
points(mds.fig, "sites", pch = 19, col = "black", select = TeethMeta$Species=="Caribbean Reef")

adonis(teethwaterdata3~TeethWaterMeta$Location)
#NMDS Relative Abundance
Teethwater_NMDS=metaMDS(teethwaterdata3,k=2)
dev.new()

ordiplot(Teethwater_NMDS, display = "sites", type = "text")
mds.fig <- ordiplot(Teethwater_NMDS, type = "none")

# plot just the samples, colour by type, pch=19 means plot a circle
points(mds.fig, "sites", pch = 19, col = "black", select = TeethWaterMeta$Location=="Teeth")
points(mds.fig, "sites", pch = 19, col = "blue", select = TeethWaterMeta$Location=="Water")

##NMDS Presence/Absence
Teethwater_NMDSpa=metaMDS(teethwaterdata3.pa,k=2)
ordiplot(Teethwater_NMDSpa, display = "sites", type = "text")

mds.fig.pa <- ordiplot(Teethwater_NMDSpa, type = "none")
# plot just the samples, colour by type, pch=19 means plot a circle
points(mds.fig.pa, "sites", pch = 19, col = "black", select = TeethWaterMeta$Location=="Teeth")
points(mds.fig.pa, "sites", pch = 19, col = "blue", select = TeethWaterMeta$Location=="Water")


library(gplots)
comm100<-teethwaterdata3.ra
comm100.sqrt<-sqrt(comm100)
dev.new()
heatmap.2(as.matrix(comm100.sqrt), col=colorRampPalette(c("floralwhite", "moccasin","navajowhite2","firebrick1", "firebrick", "firebrick4","black"))(100),dendrogram = "row",Rowv=TRUE, Colv=FALSE, keysize = 1, cexCol=.7, trace="none", margins = c(5, 12))

#Heat Map
library(gplots)
library(vegan)
library(RColorBrewer)
                   
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)

maxab <- apply(teethwaterdata3.ra, 2, max)
head(maxab)
n1 <- names(which(maxab < 0.04))
data.prop.1 <- teethwaterdata3.ra[, -which(names(teethwaterdata3.ra) %in% n1)]
heatmap(as.matrix(data.prop.1), Rowv = NA, Colv = NA, col = scaleyellowred, margins = c(10, 2))
data.dist <- vegdist(teethdata3.mod.ra, method = "bray")
row.clus <- hclust(data.dist, "aver")
heatmap(as.matrix(data.prop.1), Rowv = as.dendrogram(row.clus), Colv = NA, col = scaleyellowred, margins = c(10, 3))
data.dist.g <- vegdist(t(data.prop.1), method = "bray")
col.clus <- hclust(data.dist.g, "aver")
heatmap(as.matrix(data.prop.1), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, margins = c(10, 3))
var1 <- round(runif(n = 33, min = 1, max = 2))  # this randomly samples from a uniform distribution and rounds the result to an integer value
var1 <- replace(var1, which(var1 == 1), "deepskyblue")
var1 <- replace(var1, which(var1 == 2), "magenta")
cbind(row.names(data.prop), var1)
heatmap.2(as.matrix(data.prop.1),Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, RowSideColors = var1) 
# this puts in the annotation for the samples margins = c(10, 3))
heatmap.2(as.matrix(data.prop.1), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, RowSideColors = var1, margins = c(11, 5), trace = "none", density.info = "none", xlab = "genera", ylab = "Samples", main = "Heatmap example", lhei = c(2, 8))
# this makes the colour-key legend a little thinner
plot(annHeatmap2(as.matrix(data.prop.1),
col = colorRampPalette(c("lightyellow", "red"), space = "rgb")(51), breaks = 50,
dendrogram = list(Row = list(dendro = as.dendrogram(row.clus)), Col = list(dendro = as.dendrogram(col.clus))), legend = 3,
labels = list(Col = list(nrow = 12))))

#Rank Abundance Curves
library(goeveg)
library(matrixStats)
Teethsorted<-teethdata3.mod.ra[order(colSums(-teethdata3.mod.ra))]
SteethMinusSing<-Teethsorted[1:20]
colSums(SteethMinusSing)
Tcolm<-colMeans(SteethMinusSing)
Tcolm
Tsdmatrix<-colSds(as.matrix(SteethMinusSing))
Tsdmatrix






CRT_Meta <- read.csv(file.choose(), header=T, row.names=1)
as.data.frame(CRT_Meta)
crteethdata3<-teethdata3.ra[rownames(CRT_Meta),]
Crteethsorted<-crteethdata3[order(colSums(-crteethdata3))]
CrteethMinusSing<-Crteethsorted[1:20]
racurve(CrteethMinusSing, main = "Rank Abundance:Caribbean Reef Sharks", nlab = 10, ylog = FALSE,
        frequency = FALSE)


LT_Meta <- read.csv(file.choose(), header=T, row.names=1)
as.data.frame(LT_Meta)
Lteethdata3<-teethdata3.ra[rownames(LT_Meta),]
Lteethsorted<-Lteethdata3[order(colSums(-Lteethdata3))]
LteethMinusSing<-Lteethsorted[1:20]
racurve(LteethMinusSing, main = "Rank Abundance:Lemon Sharks", nlab = 10, ylog = FALSE,
        frequency = FALSE)

NT_Meta <- read.csv(file.choose(), header=T, row.names=1)
as.data.frame(NT_Meta)
Nteethdata3<-teethdata3.ra[rownames(NT_Meta),]
Nteethsorted<-Nteethdata3[order(colSums(-Nteethdata3))]
NteethMinusSing<-Nteethsorted[1:20]
racurve(NteethMinusSing, main = "Rank Abundance:Nurse Sharks", nlab = 8, ylog = FALSE,
        frequency = FALSE)

SBT_Meta <- read.csv(file.choose(), header=T, row.names=1)
as.data.frame(SBT_Meta)
SBteethdata3<-teethdata3.ra[rownames(SBT_Meta),]
SBteethsorted<-SBteethdata3[order(colSums(-SBteethdata3))]
SBteethMinusSing<-SBteethsorted[1:20]
racurve(SBteethMinusSing, main = "Rank Abundance:Sandbar Sharks", nlab = 4, ylog = FALSE,
        frequency = FALSE)


TT_Meta <- read.csv(file.choose(), header=T, row.names=1)
as.data.frame(TT_Meta)
Tteethdata3<-teethdata3.ra[rownames(TT_Meta),]
Tteethsorted<-Tteethdata3[order(colSums(-Tteethdata3))]
TteethMinusSing<-Tteethsorted[1:20]
racurve(TteethMinusSing, main = "Rank Abundance:Tiger Sharks", nlab = 6, ylog = FALSE,
        frequency = FALSE)

