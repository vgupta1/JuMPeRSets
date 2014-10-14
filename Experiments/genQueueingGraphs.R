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

##VG drop the CS values which aren't useul
dat2.melt = subset(dat2.melt, Method %in% c("FB", "FB2", "FBBPD", "Kingman"))


font = "CM Roman"
font = "Times New Roman"

# # add text labels or each run
# #Change everything to the tex font and size...
# g <- ggplot(aes(x=N, y=value, color=Method, group=Method), 
#                data=subset(dat2.melt, N<2000)) +
#   stat_summary(fun.data="mean_cl_boot", geom="errorbar", width=20) +
#   theme_bw(base_size=10) + 
#   theme(legend.title=element_blank(), 
#         text=element_text(family=font), 
#         legend.position="none") +
#   ylab("Waiting Time")

dat.summary = ddply(dat2.melt, ~N + Method, summarise, 
                    mean=mean(value), std=sd(value), 
                    ub=quantile(value, .9), lb=quantile(value, .1))

#most obnoxious way to get hte last kingman value
king_approx = subset(dat.summary, Method=="Kingman" & N == max(dat.summary$N))$mean

pd <- position_dodge(30)
dat2.labels = subset(dat.summary, N==1000)
dat2.labels[dat2.labels$Method == "Kingman", "mean"] = king_approx -3
dat2.labels[dat2.labels$Method == "Kingman", "N"] = 550

dat2.labels$Method = revalue(dat2.labels$Method, 
                             c(CSBPD="CS3", FBBPD="FB3"))

g<- ggplot(aes(x=N, y=mean, color=Method, group=Method), 
       data=subset(dat.summary, Method!="Kingman" & N < 2000)) + 
  geom_line(linetype="dashed",position=pd) + geom_point(position=pd) +
  geom_errorbar(aes(ymin=lb, ymax=ub), width=50, position=pd)+
  theme_bw(base_size=10) + 
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.position="none") +
  ylab("Waiting Time") + 
  geom_hline(yintercept=king_approx, linetype="dashed")


g<- g + geom_text(aes(x=N + 50, y=mean, label=Method), 
                  data=dat2.labels, 
                  color="black", 
                  family=font, 
                  size=3,
                  hjust=0) + 
    xlim(250, 1125)

ggsave("/Users/vishalgupta/Documents/Research/DataDriven Uncertainty SEts/Tex Documents/OR_Submission_v2/Figures/queue_exp2.pdf", 
       g, width=3, height=2, units="in")

# for N = 1000 std of Kingman = 597.022518, quants(31, 103)
# for N = 10000 std of Kingman = 8.7, quants(46, 67)

#include this table
# N  Method     mean       std       ub       lb
# 41 10000      FB 34.57965 0.4182423 35.17161 34.01688
# 42 10000     FB2 25.77690 0.2972324 26.18773 25.37503
# 43 10000   FBBPD 14.36647 1.2456163 15.45071 13.54583
# 44 10000 Kingman 55.06982 8.6522245 67.43304 46.03020

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
#limit to the relevant subset of methods
dat3.melt = subset(dat3.melt, Method %in% c("FB3", "FB2", "FB", "Quants", "Kingman"))
dat3.melt$Method = factor(dat3.melt$Method, c("FB", "FB2", "FB3", "Kingman", "Quants"))
dat3.melt$Method = revalue(dat3.melt$Method, c(Quants="Empirical", Kingman="King"))

#label each one with the appropriate thing
dat3.labels=subset(dat3.melt, abs(epsbar-0.51551724) < 1e-6)

g<- ggplot(aes(y=1-epsbar, x=value, color=Method, linetype=Method), 
       data=dat3.melt) + 
  theme_bw(base_size=10) + 
  theme(legend.position="none", 
        text=element_text(family=font))+
  geom_line() + geom_point() +
  ylab("Probability") + 
  xlab("Waiting Time") + xlim(0, 50)  
  
g<- g + geom_text(aes(x=value + 3, y=1-epsbar+.07, label=Method), 
              data=dat3.labels,
                  color="black", 
                  family=font, 
                  size=3,
                  hjust=0) 

ggsave("../../Tex Documents/OR_Submission_v2/Figures/queueCDF.pdf", g, 
       height=2, width=3.5, units="in")
