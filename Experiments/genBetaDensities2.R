library(ggplot2)
library(extrafont)

x_grid = seq(from=-1, to=1, by=.01)
y1 = dbeta( .5 * (x_grid + 1 ) , 4, 4)
y2 = dbeta( .5 * (x_grid + 1 ), .4, 2)
font = "Times New Roman"
font = "CM Roman"

dat = data.frame( x = x_grid, U1 = y1, U2 = y2)

g<- ggplot(aes(x=x_grid, y=U1), data=dat) + geom_line(color="red") + 
  geom_line(aes(y=U2), color="blue") + 
  theme_minimal(base_size=10) + 
  theme(text=element_text(family=font)) +
  xlab("") + ylab("") + 
  ylim(0, 4)

g1<- g + annotate("text", x=.27, y=2.1, label="u[1]", parse=TRUE, family=font) + 
    annotate("text", x=-.66, y=2.1, label="u[2]", parse=TRUE, family=font) +
  ylim(0, 2.5)

ggsave("../../Tex Documents/OR_Submission_v2/Figures/betaDensity.pdf", 
       g1, width=.35 * 6.5, height=2.2)

