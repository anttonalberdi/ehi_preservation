
setwd("/Users/anttonalberdi/github/ehi_preservation/")

#######
# Load libraries
#######

library(hilldiv)
library(ape)
library(phytools)
library(vegan)
library(nlme)
library(usedist)

#######
# Load data
#######

counts <- read.csv("data/ASV_counts.csv",row.names=1)
metadata <- read.csv("data/metadata.csv")
taxonomy <- read.csv("data/ASV_taxa.txt",row.names=1)
tree <- read.tree("data/ASVs.tre")
tree <- force.ultrametric(tree,method="extend")

#######
# Sanity check
#######

#Remove samples with <1000 reads
counts <- depth_filt(counts, 1000)
metadata <- metadata[metadata$Sample %in% colnames(counts),]
metadata$Species <- as.factor(metadata$Species)

#######
# Initial visualisation
#######

#Buffer test
qvalue=0
divtest <- div_test(counts,qvalue=qvalue,hierarchy=metadata[,c("Sample","Buffer")])
divtest <- div_test(counts,qvalue=qvalue,tree=tree,hierarchy=metadata[,c("Sample","Buffer")])
div_test_plot(divtest)

#Species test
divtest <- div_test(counts,qvalue=1,hierarchy=metadata[,c("Sample","Species")])
div_test_plot(divtest)

#Time test
divtest <- div_test(counts,qvalue=1,hierarchy=metadata[,c("Sample","Time")])
div_test_plot(divtest)

#Protocol test
divtest <- div_test(counts,qvalue=1,hierarchy=metadata[,c("Sample","Protocol")])
div_test_plot(divtest)

#Dissimilarity NMDS
qvalue=1
pairdis_q0 <- pair_dis(counts,qvalue=1)

pdf("results/nmds_q1.pdf",width=8,height=6)
dis_nmds(pairdis_q0$L1_CqN,hierarchy=metadata[,c("Sample","Species")],colour=c("#e776ae","#76aee7","#76e7a3","#d3e776","#e79676"),labels="sample")
dev.off()

#######
# Diversity differences
#######

hill_q0 <- hill_div(counts,qvalue=0)
hill_q1 <- hill_div(counts,qvalue=1)
hill_q1phy <- hill_div(counts,qvalue=1,tree=tree)
hill_table <- cbind(q0=hill_q0,q1=hill_q1,q1phy=hill_q1phy)

hill_metadata <- merge(hill_table,metadata,by.x="row.names",by.y="Sample")
rownames(hill_metadata) <- hill_metadata[,1]
hill_metadata <- hill_metadata[,-1]

summary(lme(q0~Buffer*Time,random = ~1|Species, data=hill_metadata))
summary(lme(q1~Buffer*Time,random = ~1|Species, data=hill_metadata))
summary(lme(q1phy~Buffer*Time,random = ~1|Species, data=hill_metadata))

#######
# Compositional differences
#######
#Compute pairwise dissimilarities
#q0
pairdis_q0 <- pair_dis(counts, qvalue=0)
distance <- as.dist(pairdis_q0$L1_CqN)

#q1
pairdis_q1 <- pair_dis(counts, qvalue=1)
distance <- as.dist(pairdis_q0$L1_CqN)

#q1phy (overnight)
pairdis_q1phy <- pair_dis(counts, qvalue=1, tree=tree)
distance <- as.dist(pairdis_q1phy$L1_CqN)

#Prepare
drex_samples<-metadata[metadata$Protocol == "D","Sample"]
zymo_samples<-metadata[metadata$Protocol == "Z","Sample"]
distance_drex <- dist_subset(distance, drex_samples)
distance_zymo <- dist_subset(distance, zymo_samples)

adonis(distance_drex ~ Buffer*Time*Species, metadata[metadata$Protocol == "D",], permutations = 999 )
adonis(distance_zymo ~ Buffer*Time*Species, metadata[metadata$Protocol == "Z",], permutations = 999 )
