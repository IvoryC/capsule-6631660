
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
infiles = dir(decontaminationFolder, pattern="vsnm.csv", full.names = T, recursive = T)

# create directory for output
predictionFolder = "../results/data/ivory/prediction"
suppressWarnings(dir.create(predictionFolder))

#### pcoa ####

## raed data --- later do this in loop
data = read.csv(infiles[3], row.names = 1)

# pick one ---later do this as a loop through all samples.
test1 = sample(row.names(data), size=1)
message("Leave out sample ", test1, " as the test.")

train = filter(data, rownames(data) != test1)

categories = c("Control", "SKCM")


hasCat = metadataPSMatchedDPQCFiltered$disease_type_consol %in% categories
nameHasCat = row.names(metadataPSMatchedDPQCFiltered)[which(hasCat)]
# metaCat = metadataPSMatchedDPQCFiltered[nameHasCat, ]
# datasub = data[nameHasCat, ]

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


##### boxplot #####

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

# plot5 = plot4 +
#     ggpubr::stat_compare_means(aes(y=value, x=disease_type_consol),
#                                label = "p.signif", 
#                                method = "t.test",
#                                ref.group = "0.5")
# plot5


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


#### test sample ####




message("Done!")
