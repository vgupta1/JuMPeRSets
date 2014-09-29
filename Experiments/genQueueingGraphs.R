####
# Queuing Plots
####
setwd("/Users/vishalgupta/Documents/Research/DataDriven Uncertainty SEts/JuMPeRSets/Experiments/")
library(ggplot2)
library(reshape)
library(plyr)
library(extrafont)

#Experiment 2
#N Varies
#Redo this with more points between 1000 and 2000
#Consider dropping the CS values
dat2 = read.csv(file = "exp2.csv", header=TRUE)
dat2.melt = melt(dat2, id.vars = c("iRun", "N", "n"), variable_name="Method")
font = "CM Roman"
font = "Times New Roman"

# add text labels or each run
#Change everything to the tex font and size...
g <- ggplot(aes(x=N, y=value, color=Method, group=Method), 
               data=subset(dat2.melt, Method != "Kingman" & N<2000)) +
  stat_summary(fun.data="mean_cl_boot", geom="errorbar", width=20) +
  theme_bw(base_size=10) + 
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.position="none") +
  ylab("Waiting Time")

dat.means = ddply(dat2.melt, ~N + Method, summarize, value=mean(value))
dat2.labels = subset(dat.means, N==1000)
dat2.labels$Method = revalue(dat2.labels$Method, 
                             c(CSBPD="CS3", FBBPD="FB3"))
g<- g + geom_text(aes(x=N + 100, y=value, label=Method), 
                  data=subset(dat2.labels, Method != "Kingman"), 
                  color="black", 
                  family=font, 
                  size=3) + xlim(250, 1150)

ggsave("/Users/vishalgupta/Documents/Research/DataDriven Uncertainty SEts/Tex Documents/OR_Submission_v2/Figures/queue_exp2.pdf", 
       g, width=3, height=2, units="in")

###############################
# Experiment 1
###############################
dat = read.csv(file = "exp1.csv")
dat.melt = melt(dat, id.vars = c("iRun", "n"), variable_name = "Method")
ggplot(aes(x=n, y=value, group=Method, color=Method), data=dat.melt) + 
  stat_summary(fun.data="mean_cl_boot", geom="line") + 
  scale_x_log10() + ylim(0, 100)


##############################
# Experiment 3
##############################
dat3 = read.csv(file = "exp3.csv", header=TRUE)
quants = read.csv(file="simQuants.csv", header=FALSE)
dat3 = data.frame(dat3, quants)
names(dat3)[10] = "Quants"

dat3.melt = melt(dat3, id.vars = c("n", "epsbar"), variable_name = "Method")
g<- ggplot(aes(y=1-epsbar, x=value, color=Method), 
       data=subset(dat3.melt, Method %in% c("FB3", "Quants_", "CS3_", "Kingman_"))) + 
  theme_bw(base_size=10) + 
  theme(legend.position="none", 
        text=element_text(family=font))+
  geom_line() + geom_point() +
  xlab("Probability") + 
  ylab("Waiting Time")

ggsave("exp3_3.pdf", g)  

ggsave("../../Tex Documents/OR_Submission_v2/Figures/queueCDF.pdf", g, 
       height=2, width=3, units="in")
