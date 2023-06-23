# Sterility Seed count for Skyfall and Howley - ggplots for sterility seed count. 

library(readxl)
library(readr)
library(ggplot2)
library (readxl)

getwd()

skyfall_sterility <- read.csv("Sterility seed count - R import list.csv")
skyfall_sterility

rm(skyfall_sterility)

skyfall_sterility <- read.csv("extended sterility counts.csv")
head(skyfall_sterility)
Summary_data <- tapply(skyfall_sterility$No.Seed.set, skyfall_sterility$Tray, summary)
Summary_data
sd <- tapply(skyfall_sterility$No.Seed.set, skyfall_sterility$Tray, sd)
sd
mean <- tapply(skyfall_sterility$No.Seed.set, skyfall_sterility$Tray, mean)
mean
df <- data.frame(Mean = c(3.7380952, 3.0666667, 28.3725490, 5.5238095, 0.4255319, 27.9795918), sd=c(6.953011,5.557141,7.185989,8.12175,1.514347,5.673365), Category =c("Tray 1", "Tray 2", "Tray 3", "Tray 4", "Tray 5", "Tray 6"), Insert = c(0.0, 0.1, 0.3, 0.5,0.7,0.9))

p <-ggplot(df, aes(x=Category, y= Mean,))+
  geom_point(aes(colour = category)+
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2, position=position_dodge(0.05))
p


df_2 <- data.frame(Mean = c(3.7380952, 3.0666667, 28.3725490, 5.5238095, 0.4255319, 27.9795918), sd=c(6.953011,5.557141,7.185989,8.12175,1.514347,5.673365), Treatment=as.factor (c("90% F.D", "90% F.D","Water Treated", "80% F.D","80% F.D", "Water Treated")), Category=c("Tray 1", "Tray 2", "Tray 3", "Tray 4", "Tray 5", "Tray 6"),  Insert = c(0.0, 0.1, 0.3, 0.5,0.7,0.9))

head(df_2)

P2 <-ggplot(df_2, aes(x=Category, y=Mean, fill=Treatment))+
  geom_point(aes(color= Treatment, shape=Treatment))+
  ylim(-5,50)+
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2, position=position_dodge(0.05))

P2
  

p3<- P2 +theme(plot.background = element_rect(fill="white"),panel.background = element_rect(fill="white"),axis.line.x = element_line(color="Black"), axis.line.y = element_line("black"))+
  scale_colour_manual(breaks = c("80% F.D","90% F.D","Water Treated"), values=c("red","Green","blue"))
p3

p4 <- p3 +labs(x="Tray number", y="Mean number of seed set")

p4

p5 <- p4 + theme(axis.title.x = element_text(size= 12, face = "bold", vjust = -1.5),
           axis.title.y= element_text(size = 12, face = "bold", vjust =2.5))

p5
p6 <- p5 +theme(axis.text.x = element_text(angle=60, hjust=1))
p6

p7 <- p6 + ggtitle("RGT Skyfall")

p7

p8 <- p7 + theme(plot.title = element_text(size = 16, face="bold"))
p8

p9<- p8 + theme(plot.title = element_text(hjust=0.4))
p9

#sterility analysis for Howley trays 

getwd()

Howley_sterility <- read.csv("Howley sterility counts.csv")
Howley_sterility
Howley_summary <- tapply(Howley_sterility$No.Seed.set, Howley_sterility$Tray, summary)
Howley_summary

sd_Howley <- tapply(Howley_sterility$No.Seed.set, Howley_sterility$Tray, sd)
mean_howley <- tapply(Howley_sterility$No.Seed.set, Howley_sterility$Tray, mean)

mean_howley


df_2_howley <- data.frame(Mean = c(4.642857, 5.187500, 36.724138, 1.777778, 28.675676, 1.826087 ), sd=c(8.079303,9.623350,8.263070,5.990999,8.158343,6.124854), Treatment=as.factor (c("60% F.D", "60% F.D","Water Treated", "70% F.D","Water Treated", "70% F.D")), Category=c("Tray 1", "Tray 2", "Tray 3", "Tray 4", "Tray 5", "Tray 6"),  Insert = c(0.0, 0.1, 0.3, 0.5,0.7,0.9))

Howley_p <-ggplot(df_2_howley, aes(x=Category, y=Mean, fill=Treatment))+
  geom_point(aes(color= Treatment, shape=Treatment))+
  ylim(-5,50)+
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2, position=position_dodge(0.05))

Howley_p

Howley_p2<- Howley_p +theme(plot.background = element_rect(fill="white"),panel.background = element_rect(fill="white"),axis.line.x = element_line(color="Black"), axis.line.y = element_line("black"))

Howley_p2

Howley_p3 <- Howley_p2 +labs(x="Tray number", y="")
Howley_p3

Howley_p4 <- Howley_p3 + theme(axis.title.x = element_text(size= 12, face = "bold", vjust = -1.5),
                 axis.title.y= element_text(size = 12, face = "bold", vjust =2.5))

HOwley_p5 <- Howley_p4 +theme(axis.text.x = element_text(angle=60, hjust=1))

Howley_p6 <- HOwley_p5 + ggtitle("RGT Howley")
Howley_p7<- Howley_p6 + theme(plot.title = element_text(size = 16, face="bold"))
Howley_p8<- Howley_p7 + theme(plot.title = element_text(hjust=0.4))

install.packages("gridExtra")
library("gridExtra")

grid.arrange (p9,Howley_p8, ncol=2)
