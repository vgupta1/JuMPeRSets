#VG Fix the experiment to distinguish between pstar and phat
#http://blog.revolutionanalytics.com/2012/09/how-to-use-your-favorite-fonts-in-r-charts.html
#Run that to get the font to look same as latex document
#Add the anotations for the true probabilities

# Making plots for the discrete sets from output

library(ggplot2)
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

# g1 <- g + geom_path(aes(x=u1, y=u2, group=factor(N)), 
#        data=subset(dat2, N %in% c(100, 1000, 10000)), 
#        color="red", linetype=5) 
  
#add the paths for the G
# g2 <- g1 + geom_path(aes(x=u1, y=u2, group=N), 
#                data=subset(datG, N %in% c(100, 1000, 10000)), 
#                color="blue", linetype=1) + 
#   scale_linetype_manual(breaks=c("100","1000", "10000"), 
#                         values=c(5,1, 3), name="N")
         
# g2 

ggsave("disc_sets.pdf", g1a)



ggplot(aes(x=u1, y=u2), data=datChiSq) + 
  geom_point(aes(x=X1, y=X2), color="red", data=data.frame(A)) +  
  geom_path(linetype="dotted") + 
  geom_path(data=datG) + 
    theme_bw(base_size=12) + 
  xlab("") + ylab("")



