###############
# Script for Sparcc analysis: start from csv with Species counts and end with adjacency matrix saved to csv 


library(phyloseq)
library(devtools)
library(SpiecEasi)
library(Matrix)
library(igraph)

###############


# Load Metadata
metadata<- read.csv("MVCO_metadata_AllSeqeunceSamples.csv")
metadata$sampledates <- as.Date(metadata$sampledates, "%d-%b-%Y")

row.names(metadata) <- metadata$filenames;


# Load in count data for Species assingments after cut to confidence threshhold 

readfile <- read.csv("chrono_species_gt80.csv")
inputmat <- data.matrix(readfile[,9:133] )


taxalist = as.list(readfile$Species)

## Here is crux of function 
sparccout = sparcc(t(inputmat)); 

#save.image("sparccout_revised.RData")


###############  Start Here if you want to just work from sparcc output ###############
#load('sparccout.RData')

# Load in csv with taxa broken into different levels & column for which taxa are small phytoplankton 

taxafile <- read.csv("chrono_species_gt80.csv")


# cut down network to only chlorophytes and strong positive correlations 

sparcc.graph <- abs(sparccout$Cor) >= 0.3; 
diag(sparcc.graph) <- 0; 

ig.sparcc <- adj2igraph(sparcc.graph)

degs <- degree(ig.sparcc)

indtouse = as.logical(taxafile$sum_is_chlorophyta)&degs>1
sparcc_sub.graph <- subgraph(ig.sparcc, indtouse)


#now we want to convert sparcc_sub.graph into a matrix

 sparcc_sub_mat <- as_adj(sparcc_sub.graph) #matrix.type = "adjacency")
 
 #sparcc_sub_mat.names
 
 
 # save adjacency matrix 
 write.matrix(sparcc_sub_mat,file="SparCCMatrix.csv")

 
 # save taxa names for adjacency matrix 
 B <- taxafile$Species[indtouse]
 write.csv(B, 'SparCCtaxa.csv')
 
 
 
 