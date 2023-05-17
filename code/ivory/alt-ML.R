
#### libraries ####

library(vegan)
library(ape)
library(tidyr)
library(dplyr)
library(ggplot2)

sessionInfo()

#### read meta data ####

metadataPSMatchedDPQCFiltered <- read.csv('../data/Fig4_plasma/Metadata-Plasma-Filtered-For-Analysis.csv', row.names=1)

#### find data files ####

decontaminationFolder = "../results/data/ivory/decontamination"
inputPattern = "_vsnm.csv"
infiles = dir(decontaminationFolder, pattern=inputPattern, full.names = T, recursive = T)

# create directory for output
predictionFolder = "../results/data/ivory/prediction"
suppressWarnings(dir.create(predictionFolder))

## raed data --- later do this in loop
infile = infiles[3]
message("Reading data from file: ", infile)
data = read.csv(infile, row.names = 1)

#### train ####

# pick one ---later do this as a loop through all samples.
test1 = sample(row.names(data), size=1)
message("Leave out sample ", test1, " as the test.")

train = filter(data, rownames(data) != test1)

categories = c("Control", "SKCM")

hasCat = metadataPSMatchedDPQCFiltered$disease_type_consol %in% categories
nameHasCat = row.names(metadataPSMatchedDPQCFiltered)[which(hasCat)]
# metaCat = metadataPSMatchedDPQCFiltered[nameHasCat, ]
# datasub = data[nameHasCat, ]

##### pcoa #####
dd = vegan::vegdist(train)
# dd = vegan::vegdist(data)

res = ape::pcoa(dd)

# bioplot(res)

axes = res$vectors

pcoaDF = merge(axes, metadataPSMatchedDPQCFiltered, by=0)
row.names(pcoaDF) = pcoaDF$Row.names
pcoaDF = pcoaDF %>% select(-Row.names)

plot1 = ggplot2::ggplot(data=pcoaDF[nameHasCat,]) +
    geom_point(mapping = aes(x=Axis.1, y=Axis.2, color = disease_type)) +
    ggtitle("PCoA")
#plot1
                        
plot2 = ggplot2::ggplot(data=pcoaDF[nameHasCat,]) +
    geom_point(mapping = aes(x=Axis.3, y=Axis.4, color = disease_type_consol)) +
    ggtitle("PCoA")
#plot2

axisSet = colnames(axes)
if (length(axisSet) > 12) axisSet = axisSet[1:12]


##### pcoa boxplot #####

d4 = pcoaDF %>% 
    filter(disease_type_consol %in% categories) %>% 
    select(axisSet, disease_type_consol) %>% 
    gather("value", key="Axis", -disease_type_consol) %>%
    mutate(axisNum = as.numeric(gsub("Axis.", "", Axis))) %>%
    arrange(axisNum) 

d5 = d4 %>%
    arrange(axisNum) %>%
    mutate(Axis = factor(Axis, levels=unique(Axis))) %>%
    select(-axisNum)

# plot3 = ggplot2::ggplot(data=d4) +
#     geom_boxplot(mapping = aes(x=Axis, y=value, color = disease_type_consol)) +
#     ggtitle("PCoA values split by sample type")
# plot3

plot4 = ggplot2::ggplot(data=d5) +
    geom_boxplot(mapping = aes(y=value, x=disease_type_consol, color = disease_type_consol)) +
    ggtitle("PCoA values split by sample type") +
    facet_wrap(~Axis, nrow=3)
#plot4

##### pcoa choose best #####
# Choose the single best axis for separating the data.
axisPs = sapply(axisSet, function(ax) {
    t.test(data=d4 %>% filter(Axis == ax), 
           value ~ disease_type_consol)$p.value
})
bestAxis = names(axisPs)[which(axisPs == min(axisPs))]
message("Axis ", bestAxis, " is the best at separating ", paste(categories, collapse=" from "))

# highlight the best axis in the plot
plot6 = plot4 +
    geom_rect(data = subset(d5, Axis == bestAxis), 
              fill = NA, colour = "red", 
              xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf)
plot6

# define the model,
# To the values for this 'axis' use the vector

##### pcoa test sample #####

# ruh-roh

#### prcomp ####

res = stats::prcomp(train)

axes = res$x

pcDF = merge(axes, metadataPSMatchedDPQCFiltered, by=0)
row.names(pcDF) = pcDF$Row.names
pcDF = pcDF %>% select(-Row.names)

#### classic biplots
#
# plot1 = ggplot2::ggplot(data=pcDF[nameHasCat,]) +
#     geom_point(mapping = aes(x=PC1, y=PC2, color = disease_type)) +
#     ggtitle("prcomp")
# plot1
# 
# plot2 = ggplot2::ggplot(data=pcDF[nameHasCat,]) +
#     geom_point(mapping = aes(x=PC3, y=PC4, color = disease_type)) +
#     ggtitle("prcomp")
# plot2

axisSet = colnames(axes)
if (length(axisSet) > 12) axisSet = axisSet[1:12]

d4 = pcDF %>% 
    filter(disease_type_consol %in% categories) %>% 
    select(axisSet, disease_type_consol) %>% 
    gather("value", key="PC", -disease_type_consol) %>%
    mutate(axisNum = as.numeric(gsub("PC", "", PC))) %>%
    arrange(axisNum) 

d5 = d4 %>%
    arrange(axisNum) %>%
    mutate(PC = factor(PC, levels=unique(PC))) %>%
    select(-axisNum)

# plot3 = ggplot2::ggplot(data=d4) +
#     geom_boxplot(mapping = aes(x=Axis, y=value, color = disease_type_consol)) +
#     ggtitle("PCoA values split by sample type")
# plot3

plot4 = ggplot2::ggplot(data=d5) +
    geom_boxplot(mapping = aes(y=value, x=disease_type_consol, color = disease_type_consol)) +
    ggtitle("prcomp values split by sample type") +
    facet_wrap(~PC, nrow=3)
#plot4

##### prcomp choose best #####
# Choose the single best axis for separating the data.
axisPs = sapply(axisSet, function(ax) {
    t.test(data=d4 %>% filter(PC == ax), 
           value ~ disease_type_consol)$p.value
})
bestAxis = names(axisPs)[which(axisPs == min(axisPs))]
message("Axis ", bestAxis, " is the best at separating ", paste(categories, collapse=" from "))

# highlight the best axis in the plot
plot6 = plot4 +
    geom_rect(data = subset(d5, PC == bestAxis), 
              fill = NA, colour = "red", 
              xmin = -Inf,xmax = Inf, ymin = -Inf,ymax = Inf)
plot6
plotFile = sub(inputPattern, "_prcomp-boxplot.png", infile)
ggsave(plotFile, plot=plot6)

##### prcomp test sample #####

###### sanity check ######

# pc1 rotation values (1 per taxon)
pc1.rot = res$rotation[,"PC1"]
# pc1 values (1 per sample)
pc1.vals = res$x[,"PC1"]

# sample name
train1 = row.names(train)[1]
# taxa counts for that sample
train1.dat = unlist(data[train1,])

# single value, the pc1 value for the first train sample
pc1.vals[train1]

# ....how to translate those vectors into that one value...
tp = train1.dat * pc1.rot


# from https://rdrr.io/github/nchlis/pca.utils/man/project_pca.html
# install.packages("remotes")
# remotes::install_github("nchlis/pca.utils")
pro.check = pca.utils::project_pca(Xnew = data[train1,], pc = res)
pro.check[,"PC1"] == pc1.vals[train1]


###### actual test ######

pro.test = pca.utils::project_pca(Xnew = data[test1,], pc = res)
pro.test.val = pro.test[,bestAxis]


claims = sapply(categories, function(cat){
    catVals = pcDF %>% 
        filter(disease_type_consol == cat) %>%
        select(bestAxis) %>% 
        unlist()
    t.test(catVals, pro.test.val, var.equal = T)$p.value
})
claims = claims[order(claims, decreasing = T)]

# which category had the highest p? that's our prediction.
predict = names(claims)[1]
confidenceP = claims[2]


message("Done!")
