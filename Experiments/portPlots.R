#####
# New Portfolio Experiments
#####
setwd("/Users/vishalgupta/Documents/Research/DataDriven Uncertainty SEts/JuMPeRSets/Experiments/")
library(ggplot2)
library(plyr)
library(dplyr)
library(extrafont)
library(reshape)

cv_dat = read.csv("portCrossVal500.csv")
csnp_dat = read.csv("CSNP500.csv") #Add in the vals for CS no boundary
cv_dat = rbind(cv_dat, csnp_dat)
rm(csnp_dat)

#RENAME and reorder the methods for convenience
cv_dat$Method = revalue(cv_dat$Method, c(CS="CS2", CSNSP="CS", CC="CM"))
cv_dat$Method = factor(cv_dat$Method, levels=c("M", "LCX", "CS", "DY", "CM", "CS2"))

cv.sum = cv_dat %>% group_by(Method, delta) %>% 
                    summarise(avg = mean(score), 
                              min = min(score), 
                              max=max(score))

#Plot the various Methods
g <- ggplot(aes(x=Method, y=avg), 
       data=filter(cv.sum, Method %in% c("CM", "LCX", "M", "CS"))) + 
  geom_pointrange(aes(y=avg, ymin=min, ymax=max), shape=15) +
  geom_point(aes(x=Method, y=score), color="blue", alpha=.5, 
             data=filter(cv_dat, Method %in% c("CM", "LCX", "M", "CS"))) +
  theme_bw(base_size = 10) + 
  theme(text=element_text(family=font)) + 
  xlab("") + ylab("Return")

ggsave("/Users/vishalgupta/Documents/Research/DataDriven Uncertainty SEts/Tex Documents/OR_Submission_v2/Figures/portCrossVal500.pdf", 
       g, width=3.25, height=2, units="in")


#### Graph the in-sample portolio holdings
library(reshape)
ports = read.csv("all_ports500.csv", header=TRUE)
ports$Method = revalue(ports$Method, c(CC="CM", CSNP="CS"))
  

#### 
# True Performance of those portfolios
####
out_perfs = read.csv("out_perfs500.csv", header=TRUE)
out_perfs = melt(out_perfs, variable_name = "Method")
out_perfs$Method=revalue(out_perfs$Method, 
                         c(CC="CM", CSNP="CS"))
out_perfs$Method = factor(out_perfs$Method, 
                          levels=c("M", "LCX", "CS", "CM"))
out_perfs %>% group_by(Method) %>%
            summarize(avg=mean(value), std=sd(value), sharpe=avg/std, var = quantile(value, .10))

#####
### Portfolio holdings or different Data Sets
out.ports = read.csv("outPorts500.csv")
out.ports = filter(out.ports, iRun <= 100)

out.ports$Method = revalue(out.ports$Method, c(CC="CM", CSNP="CS"))
out.ports$Method = factor(out.ports$Method, 
                          levels=c("M", "LCX", "CS", "CM"))

out.ports %>% group_by(Method) %>% 
              summarize(meanVar = mean(VaR), 
                        VaR10 = quantile(VaR, .1), 
                        Var90=quantile(VaR, .9))

#plot the portfolio holdings by asset
ports.x = melt(select(out.ports, -c(zIn, Avg, Std, VaR)),
            id.vars = c("iRun", "Method"), 
            variable_name = "Asset")

ports.avg = ports.x %>% group_by(Method, Asset) %>% 
                  summarise(avg=mean(value), q10=quantile(value, .1), q90=quantile(value, .9))

pd = position_dodge(.5)
g3<- ggplot(aes(x=Asset, y=avg, color=Method, group=Method), 
       data=filter(ports.avg, Method!="M")) +
  geom_line(aes(linetype=Method), position=pd) +
  geom_pointrange(aes(ymin=q10, ymax=q90, shape=Method), position=pd) +
  theme_bw(base_size=10) + 
  theme(legend.title = element_blank(), 
        legend.direction="horizontal", 
        legend.position=c(.63, .85), 
        text=element_text(family=font)) + 
  ylab("Holding (%)") + xlab("") + ylim(0, .22)

ggsave("/Users/vishalgupta/Documents/Research/DataDriven Uncertainty SEts/Tex Documents/OR_Submission_v2/Figures/outPorts500.pdf", 
       g3, width=3.25, height=2, units="in")


### Plot the out-sample performance
g4 <- ggplot(aes(x=VaR, group=Method, color=Method), 
       data=filter(out.ports, Method != "M")) + 
  geom_density(aes(linetype=Method)) + theme_bw(base_size=10) + 
  theme(legend.title=element_blank(), 
        text=element_text(family=font), 
        legend.position=c(.2, .65)) + 
  xlab("Worst Case Return") + ylab("")

ggsave("/Users/vishalgupta/Documents/Research/DataDriven Uncertainty SEts/Tex Documents/OR_Submission_v2/Figures/outVar500.pdf", 
       g4, width=3.25, height=2, units="in")


##### Create a unified table for the 500 Results zVals
#Avg out of sample (avg over data runs)
sum500 = out.ports %>% group_by(Method) %>%
                      summarise(Avg_z = mean(VaR))
#The In sample behavior
sum500 = inner_join(sum500, ports %>% group_by(Method) %>% select(zIn))

#The CrossVal Behavior
sum500 = inner_join(sum500, select(cv.sum, avg))
names(sum500)[4] = "CV"

#Out sample of the specific portfolio
t = out_perfs %>% group_by(Method) %>%
  summarize(avg=mean(value), std=sd(value), sharpe=avg/std, var = quantile(value, .10))
sum500 = inner_join(sum500, select(t, Method, var))
names(sum500)[5] = "z_out"

t = sum500 %>% mutate(zIn=round(zIn, 3), CV=round(CV, 3), zOut=round(z_out, 3), zAvg=round(Avg_z, 3)) %>% 
  select(-c(Avg_z, z_out))
write.csv(t, "/Users/vishalgupta/Documents/Research/DataDriven Uncertainty SEts/Tex Documents/OR_Submission_v2/Tables/portSummary.csv")
