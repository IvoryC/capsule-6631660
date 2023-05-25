
# summarize performance

#### libraries ####

library(ggplot2)
library(ggalluvial)
library(dplyr)



#### find data ####


predictionFolder = "../results--devRuns/data/ivory/prediction"
inputPattern = "prediction-summary.txt"
infiles = dir(predictionFolder, pattern=inputPattern, full.names = T, recursive = T)
message("Found ", length(infiles), " files of predictions.")


infiles = infiles[grep("trial_1", infiles)]
message("Found ", length(infiles), " files of predictions that are 'trial_1'.")

resultsDir = "../results/data/ivory/accuracySummary"
dir.create(resultsDir, recursive = TRUE)


#### read and record accuracy ####





rawFile = "../results--devRuns/data/ivory/prediction/raw/not_decontaminated_leave-1-out-prediction-summary.txt"
predNone = read.delim2(rawFile, comment.char = "#")

acc = data.frame(numBlanks=0, accuracy=sum(predNone$isCorrect)/nrow(predNone), numCorrect=sum(predNone$isCorrect))

predNone = predNone[,c("sampleID",  "prediction")]
names(predNone)[2] = "predictionRaw"
predNone$numBlanks = 0

preds = list()
#preds = lapply(infiles, function(infile) {
for (infile in infiles){
    nBlanks = strsplit(infile, split="SCRuB_noPlate_")[[1]][2]
    nBlanks = strsplit(nBlanks, "-blanks")[[1]][1]
    nBlanks = as.numeric(nBlanks)
    df = read.delim2(infile, comment.char = "#")
    df$numBlanks = nBlanks
    acc = rbind(acc, data.frame(numBlanks=nBlanks, 
                                accuracy=sum(df$isCorrect)/nrow(df), 
                                numCorrect=sum(df$isCorrect)))
    #return(df)
    preds[[infile]] = df
}
#})
names(preds) = infiles

acc = acc %>% mutate_at(names(acc), as.numeric)

#### scatter plot ####

s1 = ggplot(acc) +
    geom_hline(yintercept=0, color="white", size=2) +
    geom_hline(yintercept=1, color="white", size=2) +
    geom_point(aes(x=numBlanks, y=accuracy)) +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2), minor_breaks=seq(0,1,.1)) +
    scale_x_continuous(labels=unique(acc$numBlanks), breaks = unique(acc$numBlanks)) +
    xlab("number of blanks for decontamination") +
    ylab("accuracy") +
    ggtitle(paste(unique(predNone$predictionRaw), collapse=" vs ")) +
    theme(axis.title = element_text(size = 20)) +
    theme(plot.title = element_text(size = 25)) +
    theme(axis.text = element_text(size = 15))
s1

imgFile1 = file.path(resultsDir, "accuracy-scatterplot.png")
message("Saving image to: ", imgFile1)
ggsave(imgFile1, s1, device = png, units="in", width=6, height = 5)


#### alluvial ####

predAll = preds[[3]][,c("sampleID",  "actual", "prediction")]
names(predAll)[3] = "prediction90"
predHalf = preds[[2]][,c("sampleID",  "prediction")]
names(predHalf)[2] = "prediction45"
predFew = preds[[1]][,c("sampleID",  "prediction")]
names(predFew)[2] = "prediction23"

m0 = merge(predNone, predAll, by="sampleID")
m1 = merge(m0, predHalf, by="sampleID")
m2 = merge(m1, predFew, by="sampleID")


p2 = m2 %>% 
    group_by(predictionRaw, actual, prediction90, prediction45, prediction23) %>% 
    summarise(Freq = n())

## helpful pages
# https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html
# https://felixfan.github.io/ggplot2-remove-grid-background-margin/

a1 = ggplot(p2,
       aes(y = Freq, 
           axis1=predictionRaw, axis2 = actual, axis3 = prediction90, axis4=prediction45, axis5=prediction23,
           #axis1 = actual, axis2 = prediction90, axis3=prediction45, axis4=prediction23,
           fill = actual)) +
    theme_bw() +
    theme(panel.border = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_flow() +
    geom_stratum(alpha = .5) +
    geom_text(stat = "stratum", size = 3, aes(label = after_stat(stratum))) + 
    theme(axis.text=element_text(size=12)) +
    scale_x_continuous(breaks = 1:5, labels = c("raw", "truth", "90 blanks", "45 blanks", "23 blanks")) + 
    ggtitle("prediction of leave-1-out samples")
a1    

imgFile2 = file.path(resultsDir, "accuracy-alluvial.png")
message("Saving image to: ", imgFile2)
ggsave(imgFile2, a1, device = png, units="in", width=9, height = 7)


message("Done!")
