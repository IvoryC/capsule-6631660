
#### libraries ####

library(vegan)
library(ape)
library(tidyr)
library(dplyr)
library(ggplot2)

# install.packages("remotes")
# remotes::install_github("nchlis/pca.utils")
# library("pca.utils")
#
# not sure why I have problems loading this library. This is the function we use from it:
project_pca <- function (Xnew = NULL, pc = NULL) 
{
    string = "Plesae cite github('nchlis/pca.utils') see: https://rdrr.io/github/nchlis/pca.utils/man/project_pca.html"
    message(string); print(string)
    return(scale(Xnew, pc$center, pc$scale) %*% pc$rotation)
}


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

# interest point... is it better to train on only the categories of interest? or on all data?
# datasub = data[nameHasCat, ]

##### prcomp #####

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

# reset the factor levels so that the plots appear in order
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

#### test sample ####

###### sanity check ######

# from https://rdrr.io/github/nchlis/pca.utils/man/project_pca.html
# install.packages("remotes")
# remotes::install_github("nchlis/pca.utils")
#
# example: pca.utils::project_pca(Xnew = NULL, pc = NULL)
pro.check = project_pca(Xnew = data[train1,], pc = res)
pro.check[,"PC1"] == res$x[,"PC1"][train1]
# TRUE #--yay!

### the project_pca function is very straight forward:
# function (Xnew = NULL, pc = NULL) 
# {
#     return(scale(Xnew, pc$center, pc$scale) %*% pc$rotation)
# }

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

realCategory = metadataPSMatchedDPQCFiltered[test1,"disease_type_consol"]
isCorrect = predict == realCategory
correctOrNot = ifelse(isCorrect, "CORRECT", paste0("WRONG (is actually ", realCategory, ")"))

message("We predict (with confidence p=", round(confidenceP,3),") that sample [", test1, "] belongs in category [", predict, "], which is ", correctOrNot, "." )


message("Done!")
