# Making plots for the suppFcns from output

library(ggplot2)
library(extrafont)
setwd("/Users/vishalgupta/Documents/Research/DataDriven Uncertainty SEts/JuMPeRSets/Experiments/")

### Data for the outer octagon
A = data.frame(matrix(c(0, sqrt(2), 
            1,          1,  
      sqrt(2),   0,
      1,        -1,
      0,        -sqrt(2),
      -1,       -1,
      -sqrt(2),  0,
      -1,        1), nrow=8, byrow=TRUE))
A$probs = c(".25", ".50", ".01", ".01", ".01", ".03", ".06", ".13")

font = "Times New Roman"
font = "CM Roman"

datChiSq = read.csv("ChiSq.csv")
datG = read.csv("G.csv")
datExact = read.csv("Exact.csv")


#create the vertices with annotated probabilities
g <- ggplot(aes(x=X1, y=X2), data = A) + 
  geom_point(size=3, shape=18, color="blue") +
  theme_bw(base_size=12)  +   
  theme(text=element_text(family=font)) +
  xlab("") + ylab("") + 
  geom_text(aes(label=probs), hjust=-.2, family=font, 
            size=4) + 
  xlim(-1.6, 1.6) + ylim(-1.6, 1.6)

#add the exact cvar points
g<- g + geom_polygon(aes(x=u1, y=u2), data=datExact, fill="blue", alpha=.2)

#add the dotted points for the ChiSq
# ggplot(aes(x=u1, y=u2, group=N), 
#                data=datChiSq) +
#   geom_path(aes(linetype=factor(N)))

#merge the two data sets
datChiSq$Method = "ChiSq"
datG$Method = "G"
dat2 = rbind(datChiSq, datG)
dat2$Method = factor(dat2$Method)
dat2$Index = factor(paste(dat2$Method, dat2$N, sep="_"))

g1a <- g + geom_path(aes(x=u1, y=u2, group=Index, color=Method, linetype=factor(N)), 
                     data=subset(dat2, N %in% c(100, 1000, 10000)))

g1a<- g1a + scale_linetype_manual(breaks=c("100","1000", "10000"), 
                          values=c(5,1, 3), name="N") + 
      scale_color_manual(breaks=c("ChiSq", "G"), 
                         values=c("grey", "black"), 
                         name="", 
                         labels=c("Chi Sq.", "G"))
  
g1a <- g1a + annotate("text", x=-.5, y=.5, label="U^{CVAR}", family=font, size=4, parse=TRUE)

width = 4
height = width / 5.15 * 4.18 - .25
ggsave("../../Tex Documents/OR_Submission_v2/Figures/disSets1.pdf", 
       g1a, height=height, width=width)

ggsave("disc_sets.pdf", g1a)


########### 
# Similar set of plots for the UI and UFB sets
###########
datFB = read.csv(file ="UFB_100.csv")
datFB$Method = "FB"
datUI = read.csv(file="UICuts_100.csv")
datUI$Method = "UI"
dat_raw = read.csv("dat_100.csv", header=FALSE)

datFBInf = read.csv(file ="FBCuts_Inf.csv")
datFBInf$Method = "FB"
datUIInf = read.csv(file="UICuts_Inf.csv")
datUIInf$Method = "UI"

#massaging some data to improve plots
datFB = rbind(datFB, datFB[1, ])
datUI = rbind(datUI, datUI[1, ])
datUI$zstar[50] = 1
datUI$u2[50] = 1
datUI$zstar[150] = -1
datUI$u2[150] = -1
datFBInf = rbind(datFBInf, datFBInf[1, ])
datUIInf = rbind(datUIInf, datUIInf[1, ])
datUIInf$zstar[150] = -1
datUIInf$u2[150] = -1

dat = rbind(datFB, datUI, datFBInf, datUIInf)
dat$ID = paste(dat$Method, dat$N, sep="_")

##VG To do:  
#Repair the top portio of the graph and bottom portions to get linkage
#trip the top of the plot to get more bang for buck...
#Figure out a better color/labeling scheme..
#Add a dotted lien for the unit square
#use heavier lines to distinguish between two sets

g<- ggplot(aes(x=u1, y=u2), 
           data=subset(dat, u2>=-1)) + 
  geom_path(aes(group=ID, size=Method, color=factor(N), linetype=Method)) + 
  theme_bw(base_size=10) + 
  theme(text=element_text(family=font), 
        legend.title=element_blank(), 
        legend.position="right") +
  xlab(expression(u[1])) + ylab(expression(u[2])) + 
  geom_point(aes(x=V1, y=V2), data=dat_raw, color="blue", alpha=.5) + 
  scale_linetype_discrete(breaks=c("FB", "UI"), 
                       labels=c(expression(U^{FB}), expression(U^I))) +
  scale_size_manual(values=c(.8, 1.5), guide=FALSE) + 
  scale_color_manual(breaks=c("100", "1e+06"), 
                     values=c("black", "blue"), 
                     labels=c("100", expression(infinity)))

#Add a bounding box
box  = data.frame(matrix(c(1,          1,  
                          1,        -1,
                          -1,       -1,
                          -1,        1, 
                         1,          1), nrow=5, byrow=TRUE))
g<- g + geom_path(aes(x=X1, y=X2), data=box, linetype="dotted", color="black")
g<- g + xlim(-1, 1.4) + theme(legend.position=c(.9, .5))


ggsave("../../Tex Documents/OR_Submission_v2/Figures/UFB_UI_Comp.pdf", 
       g, width=.65 * 6.5, height=2.2)


