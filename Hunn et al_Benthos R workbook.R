rm(list=ls())
setwd()
getwd() 

library(tidyverse)
library(MASS) 
library(plyr)
library(modelr)
library(heplots)
library(ggthemes)
library(ggpubr)
library(arm)
library(effects)
library(RColorBrewer)
library(car)
Benth <- read.table("Benthic data summary.txt", header=TRUE) 
summary(Benth)

Benth$BLOCK <- as.factor(Benth$Block)
Benth$Temperature <- as.factor(Benth$Temp)
Benth$CO2 <- as.factor(Benth$Carbon_Dioxide)
Benth$FLOW <- as.factor(Benth$Flow)
Benth$SEDIMENT <- as.factor(Benth$Sediment)

#summary
par(mfrow=c(1,2))
hist(Benth$Total_abundance)
boxplot(Benth$Total_abundance)
leveneTest(Total_abundance ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)
hist(Benth$Taxon_richness)
boxplot(Benth$Taxon_richness)
leveneTest(Taxon_richness ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)
hist(Benth$EPT_total_abundance)
boxplot(Benth$EPT_total_abundance)
leveneTest(EPT_total_abundance ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)
hist(Benth$EPT_taxon_richness)
boxplot(Benth$EPT_taxon_richness)
leveneTest(EPT_taxon_richness ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)
#all good

hist(Benth$Body_size)
shapiro.test(Benth$Body_size)
leveneTest(Body_size ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)

hist(Benth$Orthoclad_body_size)
shapiro.test(Benth$Orthoclad_body_size)
leveneTest(Orthoclad_body_size ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)
#good

hist(Benth$Simpsons)
shapiro.test(Benth$Simpsons)
leveneTest(Simpsons ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)
#negatively skewed..
#square transform
sqSimpsons<-(Benth$Simpsons)^2
hist(sqSimpsons)
shapiro.test(sqSimpsons)
#not much improvement
logSimpsons<-log(Benth$Simpsons+1)
hist(logSimpsons)
#no
#exponential
expSimpsons <- exp(Benth$Simpsons)
hist(expSimpsons)
shapiro.test(expSimpsons)
#no
#arcsine
arcSimpsons<- asin(Benth$Simpsons)
hist(arcSimpsons)
shapiro.test(arcSimpsons) #definitely improved but not passing test
Benth <- Benth %>%
  mutate(arcSimpsons = asin(Simpsons))
         
leveneTest(arcSimpsons ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)
shapiro.test(Benth$arcSimpsons)



#individual taxa
hist(Benth$Cladocera_abundance) #OK
boxplot(Benth$Cladocera_abundance)
leveneTest(Cladocera_abundance ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth) #failed

hist(Benth$Orthocladiinae_abundance) #OK
boxplot(Benth$Orthocladiinae_abundance)
leveneTest(Orthocladiinae_abundance ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)

hist(Benth$Annelida_abundance) #OK
boxplot(Benth$Annelida_abundance)
leveneTest(Annelida_abundance ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)


hist(Benth$Potamopyrgus_abundance)#OK
boxplot(Benth$Potamopyrgus_abundance)
shapiro.test(Benth$Potamopyrgus_abundance)
leveneTest(Potamopyrgus_abundance ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)


hist(Benth$Copepoda_abundance) #positive skew
boxplot(Benth$Copepoda_abundance)
leveneTest(Copepoda_abundance ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)

hist(Benth$Nematoda_abundance) #skewed 
boxplot(Benth$Nematoda_abundance)


hist(Benth$Ostracoda_abundance) #OK
boxplot(Benth$Ostracoda_abundance)
leveneTest(Ostracoda_abundance ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)

hist(Benth$Tanypodinae_abundance) #OK
boxplot(Benth$Tanypodinae_abundance)
shapiro.test(Benth$Tanypodinae_abundance) #fail
leveneTest(Tanypodinae_abundance ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)

hist(Benth$Chironominae_abundance) #OK
boxplot(Benth$Chironominae_abundance)
shapiro.test(Benth$Chironominae_abundance)
leveneTest(Chironominae_abundance ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)

hist(Benth$Pycnocentrodes_abundance) #skewed 
boxplot(Benth$Pycnocentrodes_abundance)
shapiro.test(Benth$Pycnocentrodes_abundance)
leveneTest(Pycnocentrodes_abundance ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)

hist(Benth$Oxyethira_abundance) #skewed
boxplot(Benth$Oxyethira_abundance)
hist(Benth$Hydrobiosidae_abundance) #skewed 
boxplot(Benth$Hydrobiosidae_abundance)

hist(Benth$Leptophlebiidae_abundance) #skewed 
boxplot(Benth$Leptophlebiidae_abundance)

#Transformations - sqrt for all taxa for manova

Benth <- Benth %>%
  mutate(sqCladocera = sqrt(Cladocera_abundance)) 
leveneTest(sqCladocera ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)
shapiro.test(Benth$sqCladocera) #GOOD

Benth <- Benth %>%
  mutate(sqOrthoclad = sqrt(Orthocladiinae_abundance)) 
hist(Benth$sqOrthoclad)
leveneTest(sqOrthoclad ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth) #GOOD

Benth <- Benth %>%
  mutate(sqAnnelida = sqrt(Annelida_abundance)) 
leveneTest(sqAnnelida ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)#GOOD

Benth <- Benth %>%
  mutate(sqPotamopyrgus = sqrt(Potamopyrgus_abundance)) 
leveneTest(sqPotamopyrgus ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth) #GOOD

Benth <- Benth %>%
  mutate(sqCopepoda = sqrt(Copepoda_abundance)) 
leveneTest(sqCopepoda ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)
shapiro.test(Benth$sqCopepoda) 

Benth <- Benth %>%
  mutate(sqNematoda = sqrt(Nematoda_abundance)) 
leveneTest(sqNematoda ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)

Benth <- Benth %>%
  mutate(sqTanypodinae = sqrt(Tanypodinae_abundance)) 
shapiro.test(Benth$sqTanypodinae)

Benth <- Benth %>%
  mutate(sqChironominae = sqrt(Chironominae_abundance)) 
hist(Benth$Chironominae) 

Benth <- Benth %>%
  mutate(sqPycnocentrodes = sqrt(Pycnocentrodes_abundance)) 
hist(Benth$sqPycnocentrodes) 

Benth <- Benth %>%
  mutate(sqOxyethira = sqrt(Oxyethira_abundance)) 
leveneTest(sqOxyethira ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)

Benth <- Benth %>%
  mutate(sqHydrobiosidae = sqrt(Hydrobiosidae_abundance)) 
leveneTest(sqHydrobiosidae ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)

Benth <- Benth %>%
  mutate(sqLeptophlebiidae = sqrt(Leptophlebiidae_abundance)) 
leveneTest(sqLeptophlebiidae ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)

Benth <- Benth %>%
  mutate(sqOstracoda = sqrt(Ostracoda_abundance)) 
leveneTest(sqOstracoda ~ Temperature*FLOW*CO2*SEDIMENT, data= Benth)
hist(Benth$sqOstracoda)



library(ggplot2); library(plyr)
se <- function(x) sd(x)/sqrt(length(x))
summaryAbundanceTemp <- ddply(Benth, .(Temperature ), 
                      summarise, mean = mean(Total_abundance), se = se(Total_abundance))
ggplot(summaryAbundanceTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryAbundanceCarbon <- ddply(Benth, .(CO2 ), 
                           summarise, mean = mean(Total_abundance), se = se(Total_abundance))
ggplot(summaryAbundanceCarbon, aes(x = CO2, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryAbundanceFlow <- ddply(Benth, .(FLOW ), 
                             summarise, mean = mean(Total_abundance), se = se(Total_abundance))
ggplot(summaryAbundanceFlow, aes(x = FLOW, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryAbundanceSediment <- ddply(Benth, .(SEDIMENT ), 
                           summarise, mean = mean(Total_abundance), se = se(Total_abundance))
ggplot(summaryAbundanceSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

#Total abundance anova

library(car)
library(carData)
library(heplots)

totalabundanceanova <- aov(Total_abundance ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(totalabundanceanova)
etasq(totalabundanceanova)
TukeyHSD(totalabundanceanova)

#Taxon richness anova

taxonrichanova <- aov(Taxon_richness ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(taxonrichanova)
etasq(taxonrichanova)


#temp & sediment - significant effects:

summaryTaxrichnessTemp <- ddply(Benth, .(Temperature ), 
                              summarise, mean = mean(Taxon_richness), se = se(Taxon_richness))
ggplot(summaryTaxrichnessTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryTaxrichnessSediment <- ddply(Benth, .(SEDIMENT ), 
                                  summarise, mean = mean(Taxon_richness), se = se(Taxon_richness))
ggplot(summaryTaxrichnessSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

#CO2 - almost significant

summaryTaxrichCO2 <- ddply(Benth, .(CO2 ), 
                              summarise, mean = mean(Taxon_richness), se = se(Taxon_richness))
ggplot(summaryTaxrichCO2, aes(x = CO2, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

#EPT abundance anova

EPTabundanceanova <- aov(EPT_total_abundance ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(EPTabundanceanova)
etasq(EPTabundanceanova)

#temp & sediment significant again 

summaryEPTabundanceTemp <- ddply(Benth, .(Temperature ), 
                                summarise, mean = mean(EPT_total_abundance), se = se(EPT_total_abundance))
ggplot(summaryEPTabundanceTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryEPTabundanceSediment <- ddply(Benth, .(SEDIMENT ), 
                                 summarise, mean = mean(EPT_total_abundance), se = se(EPT_total_abundance))
ggplot(summaryEPTabundanceSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)


#EPT taxon richness anova

EPTtaxanova<- aov(EPT_taxon_richness ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(EPTtaxanova)
etasq(EPTtaxanova)

#temp & sediment again 
summaryEPTTaxrichnessTemp <- ddply(Benth, .(Temperature ), 
                                summarise, mean = mean(EPT_taxon_richness), se = se(EPT_taxon_richness))
ggplot(summaryTaxrichnessTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryEPTTaxrichSediment <- ddply(Benth, .(SEDIMENT ), 
                               summarise, mean = mean(EPT_taxon_richness), se = se(EPT_taxon_richness))
ggplot(summaryEPTTaxrichSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)


#Body size anova

bodysizeanova <- aov(Body_size ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(bodysizeanova)
etasq(bodysizeanova)

summaryBodysizeTemp <- ddply(Benth, .(Temperature ), 
                              summarise, mean = mean(Body_size), se = se(Body_size))
ggplot(summaryBodysizeTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryBodysizeCarbon <- ddply(Benth, .(CO2 ), 
                                summarise, mean = mean(Body_size), se = se(Body_size))
ggplot(summaryBodysizeCarbon, aes(x = CO2, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryBodysizeFlow <- ddply(Benth, .(FLOW ), 
                              summarise, mean = mean(Body_size), se = se(Body_size))
ggplot(summaryBodysizeFlow, aes(x = FLOW, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryBodysizeSediment <- ddply(Benth, .(SEDIMENT ), 
                                  summarise, mean = mean(Body_size), se = se(Body_size))
ggplot(summaryBodysizeSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)


#Orthoclad size anova

orthosizeanova <- aov(Orthoclad_body_size ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(orthosizeanova)
etasq(orthosizeanova)

summaryOrthosizeCarbon <- ddply(Benth, .(CO2 ), 
                               summarise, mean = mean(Orthoclad_body_size), se = se(Orthoclad_body_size))
ggplot(summaryOrthosizeCarbon, aes(x = CO2, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)


interaction.plot(Benth$Temperature, Benth$CO2, Benth$Orthoclad_body_size, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="Temp",
                 ylab="Orthoclad body size",
                 trace.label= "CO2",
                 main="Orthoclad body size",
                 xpd = TRUE)

interaction.plot(Benth$Temperature, Benth$SEDIMENT, Benth$Orthoclad_body_size, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="Temp",
                 ylab="Orthoclad body size",
                 trace.label= "Sediment",
                 main="Orthoclad body size",
                 xpd = TRUE)

#Simpsons anova
Simpsonsanova <- aov(arcSimpsons ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(Simpsonsanova)
etasq(Simpsonsanova)

Simpanova2 <-aov(Simpsons~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(Simpanova2)

summarySimpsTemp <- ddply(Benth, .(Temperature ), 
                             summarise, mean = mean(Simpsons), se = se(Simpsons))
ggplot(summarySimpsTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summarySimpsSediment <- ddply(Benth, .(SEDIMENT ), 
                                 summarise, mean = mean(Simpsons), se = se(Simpsons))
ggplot(summarySimpsSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)


#interactions plot for total abundance CO2*SEDIMENT

par(mfrow=c(1,1))

interaction.plot(Benth$SEDIMENT, Benth$CO2, Benth$Total_abundance, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="Sediment",
                 ylab="Total abundance",
                 trace.label= "CO2",
                 main="Total benthic invertebrates",
                 xpd = TRUE)

interaction.plot(Benth$CO2, Benth$SEDIMENT, Benth$Total_abundance, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="CO2",
                 ylab="Total abundance",
                 trace.label= "Sediment",
                 main="Total benthic invertebrates",
                 xpd = TRUE)

#interaction plot for taxon richness CO2*TEMP (almost signif)

interaction.plot(Benth$Temperature, Benth$CO2, Benth$Taxon_richness, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="Temp",
                 ylab="Taxon richness",
                 trace.label= "CO2",
                 main="Taxon richness",
                 xpd = TRUE)

#interaction plot for taxon richness FLOW*CO2 (almost signif)

interaction.plot(Benth$FLOW, Benth$CO2, Benth$Taxon_richness, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="Flow",
                 ylab="Taxon richness",
                 trace.label= "CO2",
                 main="Taxon richness",
                 xpd = TRUE)

#interaction plot for body size FLOW*SEDIMENT

interaction.plot(Benth$SEDIMENT, Benth$FLOW, Benth$Body_size, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="Sediment",
                 ylab="Body size",
                 trace.label= "Flow",
                 main="Body size",
                 xpd = TRUE)

#Four-way bar plots

#Total_abundance


topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Total_abundance),
        se = se(Total_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Total_abundance),
        se = se(Total_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Total_abundance),
        se = se(Total_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Total_abundance),
        se = se(Total_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 2500) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 2500) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

topright

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 2500) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 2500) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomright

Total_abundance_plot <- ggarrange(topleft  + font("title", size = 12),
          topright + font("title", size = 12),
          bottomleft + font("title", size = 12),
          bottomright+ font("title", size = 12),
          ncol = 2, nrow = 2,                       
          common.legend = TRUE,
          legend = "none") 
Total_abundance_plot

Total_abundance_plot1<-  annotate_figure(Total_abundance_plot,
                                              top = text_grob("Sediment added", size = 12),
                                               left = text_grob("Individuals per channel", size = 12, rot = 90),
                                         right = text_grob("Flow treatment", size = 12, rot=270),
                                              bottom = text_grob("CO2 treatment", size = 12))

Total_abundance_plot1



#individual taxa anovas

#Cladocera

Cladoceranova<- aov(Cladocera_abundance ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(Cladoceranova)
etasq(Cladoceranova)

summaryCladoceraTemp <- ddply(Benth, .(Temperature ), 
                                summarise, mean = mean(Cladocera_abundance), se = se(Cladocera_abundance))
ggplot(summaryCladoceraTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryCladoceraFlow <- ddply(Benth, .(FLOW ), 
                              summarise, mean = mean(Cladocera_abundance), se = se(Cladocera_abundance))
ggplot(summaryCladoceraFlow, aes(x = FLOW, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryCladoceraSediment <- ddply(Benth, .(SEDIMENT ), 
                              summarise, mean = mean(Cladocera_abundance), se = se(Cladocera_abundance))
ggplot(summaryCladoceraSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

interaction.plot(Benth$Temperature, Benth$SEDIMENT, Benth$Cladocera_abundance, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="Temperature",
                 ylab="Cladocera abundance",
                 trace.label= "Sediment",
                 main="Cladocera",
                 xpd = TRUE)

#SQRT Cladocera

sqCladoceranova<- aov(sqCladocera ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqCladoceranova)
etasq(sqCladoceranova)

#Orthoclad

Orthoanova<- aov(Orthocladiinae_abundance ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(Orthoanova)
etasq(Orthoanova)

sqOrthanova<- aov(sqOrthoclad ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqOrthanova)
etasq(sqOrthanova)

summaryOrthoTemp <- ddply(Benth, .(Temperature ), 
                                  summarise, mean = mean(Orthocladiinae_abundance), se = se(Orthocladiinae_abundance))
ggplot(summaryOrthoTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryOrthoFlow <- ddply(Benth, .(FLOW ), 
                                  summarise, mean = mean(Orthocladiinae_abundance), se = se(Orthocladiinae_abundance))
ggplot(summaryOrthoFlow, aes(x = FLOW, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryOrthoSediment <- ddply(Benth, .(SEDIMENT ), 
                                  summarise, mean = mean(Orthocladiinae_abundance), se = se(Orthocladiinae_abundance))
ggplot(summaryOrthoSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryOrthoCO2 <- ddply(Benth, .(CO2 ), 
                                  summarise, mean = mean(Orthocladiinae_abundance), se = se(Orthocladiinae_abundance))
ggplot(summaryOrthoCO2, aes(x = CO2, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

interaction.plot(Benth$FLOW, Benth$SEDIMENT, Benth$sqOrthoclad, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="FLOW",
                 ylab="Orthoclad abundance",
                 trace.label= "Sediment",
                 main="Cladocera",
                 xpd = TRUE)

#Annelida

sqAnnanova<- aov(sqAnnelida ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqAnnanova)
etasq(sqAnnanova)


summaryAnnFlow <- ddply(Benth, .(FLOW ), 
                          summarise, mean = mean(Annelida_abundance), se = se(Annelida_abundance))
ggplot(summaryAnnFlow, aes(x = FLOW, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryAnnSediment <- ddply(Benth, .(SEDIMENT ), 
                              summarise, mean = mean(Annelida_abundance), se = se(Annelida_abundance))
ggplot(summaryAnnSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

interaction.plot(Benth$SEDIMENT, Benth$CO2, Benth$sqAnnelida, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="Sediment",
                 ylab="Annelida abundance",
                 trace.label= "CO2",
                 main="Annelida",
                 xpd = TRUE)


#interactions - temp*sediment and CO2*temp*sediment

interaction.plot(Benth$SEDIMENT, Benth$Temperature, Benth$sqAnnelida, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="Sediment",
                 ylab="Annelida abundance",
                 trace.label= "TEMP",
                 main="Annelida",
                 xpd = TRUE)


#Potamopyrgus

sqPotamanova<- aov(sqPotamopyrgus ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqPotamanova)
etasq(sqPotamanova)

summaryPotamFlow <- ddply(Benth, .(FLOW ), 
                          summarise, mean = mean(Potamopyrgus_abundance), se = se(Potamopyrgus_abundance))
ggplot(summaryPotamFlow, aes(x = FLOW, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryPotamCO2 <- ddply(Benth, .(CO2 ), 
                         summarise, mean = mean(Potamopyrgus_abundance), se = se(Potamopyrgus_abundance))
ggplot(summaryPotamCO2, aes(x = CO2, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

interaction.plot(Benth$CO2, Benth$SEDIMENT, Benth$Potamopyrgus_abundance, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="CO2",
                 ylab="Annelida abundance",
                 trace.label= "Sediment",
                 main="Annelida",
                 xpd = TRUE)




interaction.plot(Benth$CO2, Benth$Temperature, Benth$Potamopyrgus_abundance, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="CO2",
                 ylab="Potamopyrgus abundance",
                 trace.label= "Temp",
                 main="Potamopyrgus",
                 xpd = TRUE)

#Copepoda 


sqCopanova<- aov(sqCopepoda ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqCopanova)
etasq(sqCopanova)

summaryCopepodCO2 <- ddply(Benth, .(CO2 ), 
                         summarise, mean = mean(copelog), se = se(copelog))
ggplot(summaryCopepodCO2, aes(x = CO2, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryCopepodSediment <- ddply(Benth, .(SEDIMENT), 
                           summarise, mean = mean(copelog), se = se(copelog))
ggplot(summaryCopepodSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

interaction.plot(Benth$CO2, Benth$SEDIMENT, Benth$Copepoda_abundance, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="CO2",
                 ylab="Copepoda abundance",
                 trace.label= "Sediment",
                 main="Copepoda",
                 xpd = TRUE)

#Nematoda 

Nemanova<- aov(sqNematoda ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(Nemanova)
etasq(Nemanova)

summaryNemaTemp <- ddply(Benth, .(Temperature ), 
                            summarise, mean = mean(sqNematoda), se = se(sqNematoda))
ggplot(summaryNemaTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryNemaFlow <- ddply(Benth, .(FLOW ), 
                            summarise, mean = mean(sqNematoda), se = se(sqNematoda))
ggplot(summaryNemaFlow, aes(x = FLOW, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryNemaCO2 <- ddply(Benth, .(CO2 ), 
                           summarise, mean = mean(sqNematoda), se = se(sqNematoda))
ggplot(summaryNemaCO2, aes(x = CO2, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryNemaSediment <- ddply(Benth, .(SEDIMENT ), 
                                summarise, mean = mean(sqNematoda), se = se(sqNematoda))
ggplot(summaryNemaSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

interaction.plot(Benth$Temperature, Benth$SEDIMENT, Benth$sqNematoda, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="Temperature",
                 ylab="Annelida abundance",
                 trace.label= "Sediment",
                 main="Nematoda",
                 xpd = TRUE)

interaction.plot(Benth$SEDIMENT, Benth$Temperature, Benth$sqNematoda, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="Sediment",
                 ylab="Annelida abundance",
                 trace.label= "Temp",
                 main="Nematoda",
                 xpd = TRUE)

#Ostracoda

sqOstranova<- aov(sqOstracoda ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqOstranova)
etasq(sqOstranova)

summaryOstraTemp <- ddply(Benth, .(Temperature ), 
                          summarise, mean = mean(Ostracoda_abundance), se = se(Ostracoda_abundance))
ggplot(summaryOstraTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryOstraFlow <- ddply(Benth, .(FLOW ), 
                          summarise, mean = mean(Ostracoda_abundance), se = se(Ostracoda_abundance))
ggplot(summaryOstraFlow, aes(x = FLOW, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

#Tanypodinae

sqTanypodanova<- aov(sqTanypodinae ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqTanypodanova)
etasq(sqTanypodanova)


summaryTanypodTemp <- ddply(Benth, .(Temperature ), 
                          summarise, mean = mean(Tanypodinae_abundance), se = se(Tanypodinae_abundance))
ggplot(summaryTanypodTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryTanypodFlow <- ddply(Benth, .(FLOW ), 
                          summarise, mean = mean(Tanypodinae_abundance), se = se(Tanypodinae_abundance))
ggplot(summaryTanypodFlow, aes(x = FLOW, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryTanypodSediment <- ddply(Benth, .(SEDIMENT ), 
                              summarise, mean = mean(Tanypodinae_abundance), se = se(Tanypodinae_abundance))
ggplot(summaryTanypodSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

#Chironominae - sqrt

sqChironominanova<- aov(sqChironominae ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqChironominanova)
etasq(sqChironominanova)

summaryChironominTemp <- ddply(Benth, .(Temperature ), 
                            summarise, mean = mean(sqChironominae), se = se(sqChironominae))
ggplot(summaryChironominTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryChironominFlow <- ddply(Benth, .(FLOW ), 
                            summarise, mean = mean(sqChironominae), se = se(sqChironominae))
ggplot(summaryChironominFlow, aes(x = FLOW, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryChironominCO2 <- ddply(Benth, .(CO2 ), 
                             summarise, mean = mean(sqChironominae), se = se(sqChironominae))
ggplot(summaryChironominCO2, aes(x = CO2, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryChironominSediment <- ddply(Benth, .(SEDIMENT ), 
                                summarise, mean = mean(sqChironominae), se = se(sqChironominae))
ggplot(summaryChironominSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

#Pycnocentrodes 

sqPycnoanova<- aov(sqPycnocentrodes ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqPycnoanova)
etasq(sqPycnoanova)

summaryPycnoTemp <- ddply(Benth, .(Temperature ), 
                             summarise, mean = mean(sqPycnocentrodes), se = se(sqPycnocentrodes))
ggplot(summaryPycnoTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)


#Oxyethira

sqOxyethiranova<- aov(sqOxyethira ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqOxyethiranova)
etasq(sqOxyethiranova)

summaryOxyTemp <- ddply(Benth, .(Temperature ), 
                          summarise, mean = mean(sqOxyethira), se = se(sqOxyethira))
ggplot(summaryOxyTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryOxyFlow <- ddply(Benth, .(FLOW ), 
                        summarise, mean = mean(sqOxyethira), se = se(sqOxyethira))
ggplot(summaryOxyFlow, aes(x = FLOW, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryOxySediment <- ddply(Benth, .(SEDIMENT ), 
                        summarise, mean = mean(sqOxyethira), se = se(sqOxyethira))
ggplot(summaryOxySediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)


#Hydrobiosidae

sqHydrobiosidaenova<- aov(sqHydrobiosidae ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqHydrobiosidaenova)
etasq(sqHydrobiosidaenova)

summaryHydrobioTemp <- ddply(Benth, .(Temperature ), 
                        summarise, mean = mean(sqHydrobiosidae), se = se(sqHydrobiosidae))
ggplot(summaryHydrobioTemp, aes(x = Temperature, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryHydrobioSediment <- ddply(Benth, .(SEDIMENT ), 
                            summarise, mean = mean(sqHydrobiosidae), se = se(sqHydrobiosidae))
ggplot(summaryHydrobioSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

interaction.plot(Benth$FLOW, Benth$SEDIMENT, Benth$sqHydrobiosidae, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="Flow",
                 ylab="Hydrobiosidae abundance",
                 trace.label= "Sediment",
                 main="Hydrobiosidae",
                 xpd = TRUE)

#Leptophlebiidae

sqLeptophlebiidaeanova<- aov(sqLeptophlebiidae ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqLeptophlebiidaeanova)
etasq(sqLeptophlebiidaeanova)

summaryLeptoCO2 <- ddply(Benth, .(CO2 ), 
                             summarise, mean = mean(sqLeptophlebiidae), se = se(sqLeptophlebiidae))
ggplot(summaryLeptoCO2, aes(x = CO2, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

summaryLeptoSediment <- ddply(Benth, .(SEDIMENT ), 
                                 summarise, mean = mean(sqLeptophlebiidae), se = se(sqLeptophlebiidae))
ggplot(summaryLeptoSediment, aes(x = SEDIMENT, y= mean)) +
  geom_bar(position="dodge", stat="identity", fill="skyblue")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),position=position_dodge(0.9), width=0.2)

interaction.plot(Benth$CO2, Benth$FLOW, Benth$Leptophlebiidae_abundance, type="b", col=c(1:3),
                 leg.bty="o", leg.bg="white", lwd=2, pch=c(18,20,22),
                 xlab="CO2",
                 ylab="Leptophlebiidae abundance",
                 trace.label= "FLOW",
                 main="Leptophlebiidae",
                 xpd = TRUE)


#3-way interactions - 4-way bar plots
#Cladocera 

topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Cladocera_abundance),
        se = se(Cladocera_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Cladocera_abundance),
        se = se(Cladocera_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Cladocera_abundance),
        se = se(Cladocera_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Cladocera_abundance),
        se = se(Cladocera_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 1000) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=16),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 1000) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=16),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 1000) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=16),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 1000) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=16),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Cladocera_plot <- ggarrange(topleft  + font("title", size = 12),
                                  topright + font("title", size = 12),
                                  bottomleft + font("title", size = 12),
                                  bottomright+ font("title", size = 12),
                                  ncol = 2, nrow = 2,                       
                                  common.legend = TRUE,
                                  legend = "none") 

Cladocera_plot1<-  annotate_figure(Cladocera_plot,
                                         top = text_grob("Sediment added", size = 12),
                                         left = text_grob("Individuals per channel", size = 12, rot = 90),
                                         right = text_grob("Flow treatment", size = 12, rot=270),
                                         bottom = text_grob("CO2 treatment", size = 12))

Cladocera_plot1


#Orthocladiinae 4 way plot

topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Orthocladiinae_abundance),
        se = se(Orthocladiinae_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Orthocladiinae_abundance),
        se = se(Orthocladiinae_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Orthocladiinae_abundance),
        se = se(Orthocladiinae_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Orthocladiinae_abundance),
        se = se(Orthocladiinae_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 650) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 650) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 650) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 650) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Orthocladiinae_plot <- ggarrange(topleft  + font("title", size = 12),
                            topright + font("title", size = 12),
                            bottomleft + font("title", size = 12),
                            bottomright+ font("title", size = 12),
                            ncol = 2, nrow = 2,                       
                            common.legend = TRUE,
                            legend = "none") 

Orthocladiinae_plot1<-  annotate_figure(Orthocladiinae_plot,
                                   top = text_grob("Sediment added", size = 12),
                                   left = text_grob("Individuals per channel", size = 12, rot = 90),
                                   right = text_grob("Flow treatment", size = 12, rot=270),
                                   bottom = text_grob("CO2 treatment", size = 12))

Orthocladiinae_plot1

#Annelida 4 way plot

topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Annelida_abundance),
        se = se(Annelida_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Annelida_abundance),
        se = se(Annelida_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Annelida_abundance),
        se = se(Annelida_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Annelida_abundance),
        se = se(Annelida_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 600) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 600) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 600) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 600) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Annelida_plot <- ggarrange(topleft  + font("title", size = 12),
                                 topright + font("title", size = 12),
                                 bottomleft + font("title", size = 12),
                                 bottomright+ font("title", size = 12),
                                 ncol = 2, nrow = 2,                       
                                 common.legend = TRUE,
                                 legend = "none") 

Annelida_plot1<-  annotate_figure(Annelida_plot,
                                        top = text_grob("Sediment added", size = 12),
                                        left = text_grob("Individuals per channel", size = 12, rot = 90),
                                        right = text_grob("Flow treatment", size = 12, rot=270),
                                        bottom = text_grob("CO2 treatment", size = 12))

Annelida_plot1


#Potamopyrgus

topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Potamopyrgus_abundance),
        se = se(Potamopyrgus_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Potamopyrgus_abundance),
        se = se(Potamopyrgus_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Potamopyrgus_abundance),
        se = se(Potamopyrgus_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Potamopyrgus_abundance),
        se = se(Potamopyrgus_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 400) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 300) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 300) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 300) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Potamopyrgus_plot <- ggarrange(topleft  + font("title", size = 12),
                              topright + font("title", size = 12),
                              bottomleft + font("title", size = 12),
                              bottomright+ font("title", size = 12),
                              ncol = 2, nrow = 2,                       
                              common.legend = TRUE,
                              legend = "none") 

Potamopyrgus_plot1<-  annotate_figure(Potamopyrgus_plot,
                                     top = text_grob("Sediment added", size = 12),
                                     left = text_grob("Individuals per channel", size = 12, rot = 90),
                                     right = text_grob("Flow treatment", size = 12, rot=270),
                                     bottom = text_grob("CO2 treatment", size = 12))

Potamopyrgus_plot1

#copepoda

topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Copepoda_abundance),
        se = se(Copepoda_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Copepoda_abundance),
        se = se(Copepoda_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Copepoda_abundance),
        se = se(Copepoda_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Copepoda_abundance),
        se = se(Copepoda_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 400) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 400) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 400) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 400) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Copepoda_plot <- ggarrange(topleft  + font("title", size = 12),
                               topright + font("title", size = 12),
                               bottomleft + font("title", size = 12),
                               bottomright+ font("title", size = 12),
                               ncol = 2, nrow = 2,                       
                               common.legend = TRUE,
                               legend = "none") 

Copepoda_plot1<-  annotate_figure(Copepoda_plot,
                                      top = text_grob("Sediment added", size = 12),
                                      left = text_grob("Individuals per channel", size = 12, rot = 90),
                                      right = text_grob("Flow treatment", size = 12, rot=270),
                                      bottom = text_grob("CO2 treatment", size = 12))

Copepoda_plot1

#Nematoda

topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Nematoda_abundance),
        se = se(Nematoda_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Nematoda_abundance),
        se = se(Nematoda_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Nematoda_abundance),
        se = se(Nematoda_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Nematoda_abundance),
        se = se(Nematoda_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 250) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 250) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 250) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 250) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Nematoda_plot <- ggarrange(topleft  + font("title", size = 12),
                           topright + font("title", size = 12),
                           bottomleft + font("title", size = 12),
                           bottomright+ font("title", size = 12),
                           ncol = 2, nrow = 2,                       
                           common.legend = TRUE,
                           legend = "none") 

Nematoda_plot1<-  annotate_figure(Nematoda_plot,
                                  top = text_grob("Sediment added", size = 12),
                                  left = text_grob("Individuals per channel", size = 12, rot = 90),
                                  right = text_grob("Flow treatment", size = 12, rot=270),
                                  bottom = text_grob("CO2 treatment", size = 12))

Nematoda_plot1

#Ostracoda

topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Ostracoda_abundance),
        se = se(Ostracoda_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Ostracoda_abundance),
        se = se(Ostracoda_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Ostracoda_abundance),
        se = se(Ostracoda_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Ostracoda_abundance),
        se = se(Ostracoda_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 125) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 125) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 125) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 125) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Ostracoda_plot <- ggarrange(topleft  + font("title", size = 12),
                           topright + font("title", size = 12),
                           bottomleft + font("title", size = 12),
                           bottomright+ font("title", size = 12),
                           ncol = 2, nrow = 2,                       
                           common.legend = TRUE,
                           legend = "none") 

Ostracoda_plot1<-  annotate_figure(Ostracoda_plot,
                                  top = text_grob("Sediment added", size = 12),
                                  left = text_grob("Individuals per channel", size = 12, rot = 90),
                                  right = text_grob("Flow treatment", size = 12, rot=270),
                                  bottom = text_grob("CO2 treatment", size = 12))

Ostracoda_plot1

#Tanypodinae


topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Tanypodinae_abundance),
        se = se(Tanypodinae_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Tanypodinae_abundance),
        se = se(Tanypodinae_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Tanypodinae_abundance),
        se = se(Tanypodinae_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Tanypodinae_abundance),
        se = se(Tanypodinae_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 100) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 100) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 100) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 100) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Tanypod_plot <- ggarrange(topleft  + font("title", size = 12),
                            topright + font("title", size = 12),
                            bottomleft + font("title", size = 12),
                            bottomright+ font("title", size = 12),
                            ncol = 2, nrow = 2,                       
                            common.legend = TRUE,
                            legend = "none") 

Tanypod_plot1<-  annotate_figure(Tanypod_plot,
                                   top = text_grob("Sediment added", size = 12),
                                   left = text_grob("Individuals per channel", size = 12, rot = 90),
                                   right = text_grob("Flow treatment", size = 12, rot=270),
                                   bottom = text_grob("CO2 treatment", size = 12))

Tanypod_plot1

#Chironominae



topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Chironominae_abundance),
        se = se(Chironominae_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Chironominae_abundance),
        se = se(Chironominae_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Chironominae_abundance),
        se = se(Chironominae_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Chironominae_abundance),
        se = se(Chironominae_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 120) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 120) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 120) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 120) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Chironomin_plot <- ggarrange(topleft  + font("title", size = 12),
                          topright + font("title", size = 12),
                          bottomleft + font("title", size = 12),
                          bottomright+ font("title", size = 12),
                          ncol = 2, nrow = 2,                       
                          common.legend = TRUE,
                          legend = "none") 

Chironomin_plot1<-  annotate_figure(Chironomin_plot,
                                 top = text_grob("Sediment added", size = 12),
                                 left = text_grob("Individuals per channel", size = 12, rot = 90),
                                 right = text_grob("Flow treatment", size = 12, rot=270),
                                 bottom = text_grob("CO2 treatment", size = 12))

Chironomin_plot1

#Pycnocentrodes

topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Pycnocentrodes_abundance),
        se = se(Pycnocentrodes_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Pycnocentrodes_abundance),
        se = se(Pycnocentrodes_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Pycnocentrodes_abundance),
        se = se(Pycnocentrodes_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Pycnocentrodes_abundance),
        se = se(Pycnocentrodes_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 40) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 40) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 40) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 40) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Pycnocentrodes_plot <- ggarrange(topleft  + font("title", size = 12),
                           topright + font("title", size = 12),
                           bottomleft + font("title", size = 12),
                           bottomright+ font("title", size = 12),
                           ncol = 2, nrow = 2,                       
                           common.legend = TRUE,
                           legend = "none") 

Pycnocentrodes_plot1<-  annotate_figure(Pycnocentrodes_plot,
                                  top = text_grob("Sediment added", size = 12),
                                  left = text_grob("Individuals per channel", size = 12, rot = 90),
                                  right = text_grob("Flow treatment", size = 12, rot=270),
                                  bottom = text_grob("CO2 treatment", size = 12))

Pycnocentrodes_plot1

# Oxyethira


topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Oxyethira_abundance),
        se = se(Oxyethira_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Oxyethira_abundance),
        se = se(Oxyethira_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Oxyethira_abundance),
        se = se(Oxyethira_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Oxyethira_abundance),
        se = se(Oxyethira_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 40) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 40) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 40) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 40) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Oxyethira_plot <- ggarrange(topleft  + font("title", size = 12),
                                 topright + font("title", size = 12),
                                 bottomleft + font("title", size = 12),
                                 bottomright+ font("title", size = 12),
                                 ncol = 2, nrow = 2,                       
                                 common.legend = TRUE,
                                 legend = "none") 

Oxyethira_plot1<-  annotate_figure(Oxyethira_plot,
                                        top = text_grob("Sediment added", size = 12),
                                        left = text_grob("Individuals per channel", size = 12, rot = 90),
                                        right = text_grob("Flow treatment", size = 12, rot=270),
                                        bottom = text_grob("CO2 treatment", size = 12))

Oxyethira_plot1

#Hydrobiosidae

topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Hydrobiosidae_abundance),
        se = se(Hydrobiosidae_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Hydrobiosidae_abundance),
        se = se(Hydrobiosidae_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Hydrobiosidae_abundance),
        se = se(Hydrobiosidae_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Hydrobiosidae_abundance),
        se = se(Hydrobiosidae_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 30) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 30) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 30) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 30) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Hydrobiosidae_plot <- ggarrange(topleft  + font("title", size = 12),
                            topright + font("title", size = 12),
                            bottomleft + font("title", size = 12),
                            bottomright+ font("title", size = 12),
                            ncol = 2, nrow = 2,                       
                            common.legend = TRUE,
                            legend = "none") 

Hydrobiosidae_plot1<-  annotate_figure(Hydrobiosidae_plot,
                                   top = text_grob("Sediment added", size = 12),
                                   left = text_grob("Individuals per channel", size = 12, rot = 90),
                                   right = text_grob("Flow treatment", size = 12, rot=270),
                                   bottom = text_grob("CO2 treatment", size = 12))

Hydrobiosidae_plot1

#Leptophlebiidae

topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Leptophlebiidae_abundance),
        se = se(Leptophlebiidae_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Leptophlebiidae_abundance),
        se = se(Leptophlebiidae_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Leptophlebiidae_abundance),
        se = se(Leptophlebiidae_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Leptophlebiidae_abundance),
        se = se(Leptophlebiidae_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 30) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 30) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 30) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 30) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Leptophlebiidae_plot <- ggarrange(topleft  + font("title", size = 12),
                                topright + font("title", size = 12),
                                bottomleft + font("title", size = 12),
                                bottomright+ font("title", size = 12),
                                ncol = 2, nrow = 2,                       
                                common.legend = TRUE,
                                legend = "none") 

Leptophlebiidae_plot1<-  annotate_figure(Leptophlebiidae_plot,
                                       top = text_grob("Sediment added", size = 12),
                                       left = text_grob("Individuals per channel", size = 12, rot = 90),
                                       right = text_grob("Flow treatment", size = 12, rot=270),
                                       bottom = text_grob("CO2 treatment", size = 12))

Leptophlebiidae_plot1

#MANOVA

community_manova <- manova(cbind(sqCladocera, sqOrthoclad, sqAnnelida, sqPotamopyrgus, sqCopepoda, sqNematoda, sqOstracoda, sqTanypodinae, sqChironominae, sqPycnocentrodes, sqOxyethira, sqHydrobiosidae, sqLeptophlebiidae) ~ Temperature*CO2*FLOW*SEDIMENT+BLOCK, data = Benth)
summary(community_manova)
community_manova
etasq(community_manova, anova=TRUE, partial=TRUE)
summary.aov(community_manova)


#all four way plots
#tax richness


topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Taxon_richness),
        se = se(Taxon_richness))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Taxon_richness),
        se = se(Taxon_richness))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Taxon_richness),
        se = se(Taxon_richness))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Taxon_richness),
        se = se(Taxon_richness))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 20) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 20) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 20) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 20) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

TaxRich_plot <- ggarrange(topleft  + font("title", size = 12),
                                 topright + font("title", size = 12),
                                 bottomleft + font("title", size = 12),
                                 bottomright+ font("title", size = 12),
                                 ncol = 2, nrow = 2,                       
                                 common.legend = TRUE,
                          legend='none') 

Taxrich_plot1<-  annotate_figure(TaxRich_plot,
                                        top = text_grob("Sediment added", size = 12),
                                        left = text_grob("Taxa per channel", size = 12, rot = 90),
                                        right = text_grob("Flow treatment", size = 12, rot=270),
                                        bottom = text_grob("CO2 treatment", size = 12))

Taxrich_plot1

#EPT abundance


topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(EPT_total_abundance),
        se = se(EPT_total_abundance))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(EPT_total_abundance),
        se = se(EPT_total_abundance))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(EPT_total_abundance),
        se = se(EPT_total_abundance))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(EPT_total_abundance),
        se = se(EPT_total_abundance))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 100) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 100) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 100) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 100) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

EPTabundance_plot <- ggarrange(topleft  + font("title", size = 12),
                          topright + font("title", size = 12),
                          bottomleft + font("title", size = 12),
                          bottomright+ font("title", size = 12),
                          ncol = 2, nrow = 2,                       
                          common.legend = TRUE,
                          legend='none') 

EPTabundance_plot1<-  annotate_figure(EPTabundance_plot,
                                 top = text_grob("Sediment added", size = 12),
                                 left = text_grob("EPT individuals per channel", size = 12, rot = 90),
                                 right = text_grob("Flow treatment", size = 12, rot=270),
                                 bottom = text_grob("CO2 treatment", size = 12))
EPTabundance_plot1

#EPT tax richness plot


topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(EPT_taxon_richness),
        se = se(EPT_taxon_richness))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(EPT_taxon_richness),
        se = se(EPT_taxon_richness))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(EPT_taxon_richness),
        se = se(EPT_taxon_richness))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(EPT_taxon_richness),
        se = se(EPT_taxon_richness))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 8) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 8) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 8) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 8) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

EPTtaxrich_plot <- ggarrange(topleft  + font("title", size = 12),
                               topright + font("title", size = 12),
                               bottomleft + font("title", size = 12),
                               bottomright+ font("title", size = 12),
                               ncol = 2, nrow = 2,                       
                               common.legend = TRUE,
                               legend='none') 

EPTtax_plot1<-  annotate_figure(EPTtaxrich_plot,
                                      top = text_grob("Sediment added", size = 12),
                                      left = text_grob("EPT taxa per channel", size = 12, rot = 90),
                                      right = text_grob("Flow treatment", size = 12, rot=270),
                                      bottom = text_grob("CO2 treatment", size = 12))
EPTtax_plot1


#simpsons D plot


topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Simpsons),
        se = se(Simpsons))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Simpsons),
        se = se(Simpsons))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Simpsons),
        se = se(Simpsons))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Simpsons),
        se = se(Simpsons))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0.0, 1) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 1) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 1) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 1) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Simpsons_plot <- ggarrange(topleft  + font("title", size = 12),
                               topright + font("title", size = 12),
                               bottomleft + font("title", size = 12),
                               bottomright+ font("title", size = 12),
                               ncol = 2, nrow = 2,                       
                               common.legend = TRUE,
                               legend='none') 

Simpsons_plot1<-  annotate_figure(Simpsons_plot,
                                      top = text_grob("Sediment added", size = 12),
                                      left = text_grob("Simpson's Diversity Index", size = 12, rot = 90),
                                      right = text_grob("Flow treatment", size = 12, rot=270),
                                      bottom = text_grob("CO2 treatment", size = 12))

Simpsons_plot1

#Body size plot


topleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Body_size),
        se = se(Body_size))

topright_subset <- Benth %>% 
  filter(Sediment=='Yes',
         Flow=='Constant') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Body_size),
        se = se(Body_size))

bottomleft_subset <- Benth %>% 
  filter(Sediment=='No',
         Flow=='Variable') %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Body_size),
        se = se(Body_size))

bottomright_subset <- Benth %>% 
  filter(Sediment=="Yes",
         Flow=="Variable") %>%  
  ddply(c("CO2", "Temperature"),
        summarise,
        mean = mean(Body_size),
        se = se(Body_size))

topleft <- ggplot(topleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0.0, 2.5) +
  labs(title = 'No', x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))

topleft

topright <- ggplot(topright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 2.5) +
  labs(title = 'Yes', x = NULL, y = NULL, tag= 'Constant')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(plot.title= element_text(size = 14, hjust=0.5),
        strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12,angle=-90),
        plot.tag.position=c(1.0, 0.5))

bottomleft <- ggplot(bottomleft_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 2.5) +
  labs(title = NULL, x = NULL, y = NULL)+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))


bottomright <- ggplot(bottomright_subset, aes(x = CO2, y = mean, fill = Temperature)) +
  geom_bar(position = "dodge", stat = "identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(0.9), width=0.2) + 
  ylim(0, 2.5) +
  labs(title = NULL, x = NULL, y = NULL, tag = 'Variable')+
  scale_fill_manual(values = c("lightblue2", "tan1")) + theme_classic() +
  theme(strip.text = element_text(size = 10),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        plot.tag=element_text(size=12, angle=-90),
        plot.tag.position=c(1.0, 0.5))

Body_plot <- ggarrange(topleft  + font("title", size = 12),
                           topright + font("title", size = 12),
                           bottomleft + font("title", size = 12),
                           bottomright+ font("title", size = 12),
                           ncol = 2, nrow = 2,                       
                           common.legend = TRUE,
                           legend='none') 

Body_plot1<-  annotate_figure(Body_plot,
                                  top = text_grob("Sediment added", size = 12),
                                  left = text_grob("Invertebrate body length (mm)", size = 12, rot = 90),
                                  right = text_grob("Flow treatment", size = 12, rot=270),
                                  bottom = text_grob("CO2 treatment", size = 12))

Body_plot1


#NMDS

rm(list=ls())
install.packages("vegan")

library(vegan)
set.seed(2)

communitydata<-read.table("NMDSmatrix.txt", header = TRUE)
NMDSgroups<-read.table("NMDS groups.txt", header = TRUE)

community_relative <-decostand(communitydata, method = "total")

community_distmat <-  vegdist(community_relative, method = "bray", autotransform = TRUE)

community_distmat <- 
  as.matrix(community_distmat, labels = T)
write.table(community_distmat, "community_distmat_distmat.txt")

community_distmat_NMDS <-
  metaMDS(community_distmat,
          distance = "bray",
          k = 2,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)
sppscores(community_distmat_NMDS)<-communitydata

goodness(community_distmat_NMDS) 

stressplot(community_distmat_NMDS) 



# ggplot

site.scrs <- as.data.frame(scores(community_distmat_NMDS, display = "sites")) 
site.scrs <- cbind(site.scrs, Site = rownames(site.scrs)) 

site.scrs <- cbind(site.scrs, CO2 = NMDSgroups$CO2)
site.scrs <- cbind(site.scrs, Temperature = NMDSgroups$Temp)
site.scrs <- cbind(site.scrs, Flow = NMDSgroups$Flow)
site.scrs <- cbind(site.scrs, Sediment = NMDSgroups$Sediment)

head(site.scrs)

species.scores <- as.data.frame(scores(community_distmat_NMDS, "species"))  
head(species.scores)

spp.scrs<-cbind(species.scores, Species=rownames(species.scores))

head(spp.scrs)

#CO2 plot

hull_datacarbon <- 
  site.scrs %>%
  drop_na() %>%
  group_by(CO2) %>% 
  slice(chull(NMDS1, NMDS2))


nmds.plotcarbon <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, colour = factor(site.scrs$CO2)), size = 1)+ 
  coord_fixed()+
  theme_classic()+ 
  scale_colour_manual(values = c("violet","springgreen3"))+
  scale_fill_manual(values = c("violet","springgreen3"))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 0.5, linetype = "solid"))+
  labs(colour = "CO2")+ # add legend labels for Management and Landuse
  theme(legend.position = "right", legend.text = element_text(size = 9), legend.title = element_text(size = 10), axis.text = element_text(size = 10))+
  geom_polygon(data = hull_datacarbon,
               aes(fill = CO2,
                   colour = CO2),
               alpha = 0.3,
               show.legend = TRUE)+
  ggrepel::geom_text_repel(data = spp.scrs, aes(x=NMDS1, y=NMDS2,label = Species), cex = 2.7, direction = "both", segment.size = 0.25)

nmds.plotcarbon + labs(title = "CO2") 

#TEMP

hull_datatemp <- 
  site.scrs %>%
  drop_na() %>%
  group_by(Temperature) %>% 
  slice(chull(NMDS1, NMDS2))


nmds.plotTemp <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(site.scrs$Temperature)), size = 1)+ 
  coord_fixed()+
  theme_classic()+ 
  scale_colour_manual(values = c("skyblue", "tan1"))+
  scale_fill_manual(values = c("skyblue", "tan1"))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 0.5, linetype = "solid"))+
  labs(colour = "Temperature")+ 
  theme(legend.position = "right", legend.text = element_text(size = 9), legend.title = element_text(size = 10), axis.text = element_text(size = 10))+
  geom_polygon(data = hull_datatemp,
               aes(fill = Temperature,
                   colour = Temperature),
               alpha = 0.3,
               show.legend = TRUE)+
  ggrepel::geom_text_repel(data = spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 2.7, direction = "both", segment.size = 0.25)

nmds.plotTemp + labs(title = "Temperature") 

#FLOW

hull_dataflow <- 
  site.scrs %>%
  drop_na() %>%
  group_by(Flow) %>% 
  slice(chull(NMDS1, NMDS2))


nmds.plotFlow <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ 
  geom_point(aes(NMDS1, NMDS2, colour = factor(site.scrs$Flow)), size = 1)+ 
  coord_fixed()+
  theme_classic()+ 
  scale_colour_manual(values = c("goldenrod2","mediumpurple1"))+
  scale_fill_manual(values = c("goldenrod2","mediumpurple1"))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 0.5, linetype = "solid"))+
  labs(colour = "Flow")+ # add legend labels for Management and Landuse
  theme(legend.position = "right", legend.text = element_text(size = 9), legend.title = element_text(size = 10), axis.text = element_text(size = 10))+
  geom_polygon(data = hull_dataflow,
               aes(fill = Flow,
                   colour = Flow),
               alpha = 0.3,
               show.legend = TRUE)+
  ggrepel::geom_text_repel(data = spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 2.7, direction = "both", segment.size = 0.25)

nmds.plotFlow + labs(title = "Flow") 

#SEDIMENT

hull_dataSediment <- 
  site.scrs %>%
  drop_na() %>%
  group_by(Sediment) %>% 
  slice(chull(NMDS1, NMDS2))


nmds.plotSediment <- ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  geom_point(aes(NMDS1, NMDS2, colour = factor(site.scrs$Sediment)), size = 1)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic()+ 
  scale_colour_manual(values = c("yellow3","grey"))+
  scale_fill_manual(values = c("yellow3","grey"))+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 0.5, linetype = "solid"))+
  labs(colour = "Sediment")+ # add legend labels for Management and Landuse
  theme(legend.position = "right", legend.text = element_text(size = 9), legend.title = element_text(size = 10), axis.text = element_text(size = 10))+
  geom_polygon(data = hull_dataSediment,
               aes(fill = Sediment,
                   colour = Sediment),
               alpha = 0.3,
               show.legend = TRUE)+
  ggrepel::geom_text_repel(data = spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 2.7, direction = "both", segment.size = 0.25)#add labels for species, use ggrepel::ge


nmds.plotSediment + labs(title = "Sediment") #displays plot 


#interaction ggplots
library(dplyr)
detach(package:plyr)

Benth %>% 
  group_by(SEDIMENT, CO2) %>% 
  summarise(mean = mean(Total_abundance), se = se(Total_abundance)) -> totabundanceint
limits = aes(ymax = mean + se, ymin=mean - se)


totabundanceint %>% 
  ggplot() +
  aes(x = CO2, y = mean, color = SEDIMENT) +
  geom_line(aes(group = SEDIMENT)) +
  geom_errorbar(limits, width=0.1)+
  geom_point()+
  ylab('Individuals per channel')+
  labs(colour = "Sediment")+
  scale_colour_manual(values = c("goldenrod","grey2"))+
  theme_classic()

#other way around

totabundanceint %>% 
  ggplot() +
  aes(x = SEDIMENT, y = mean, color = CO2) +
  geom_line(aes(group = CO2)) +
  geom_errorbar(limits, width=0.1)+
  geom_point()+
  ylab('Individuals per channel')+
  xlab('Sediment')+
  labs(colour = "CO2")+
  scale_colour_manual(values = c("violet","springgreen3"))+
  theme_classic()


#body size flow x sediment

Benth %>% 
  group_by(SEDIMENT, FLOW) %>% 
  summarise(mean = mean(Body_size), se = se(Body_size)) -> bodysizeint
limits = aes(ymax = mean + se, ymin=mean - se)


bodysizeint %>% 
  ggplot() +
  aes(x = FLOW, y = mean, color = SEDIMENT) +
  geom_line(aes(group = SEDIMENT)) +
  geom_errorbar(limits, width=0.1)+
  geom_point()+
  ylab('Invertebrate body length (mm)')+
  labs(colour = "Sediment")+
  scale_colour_manual(values = c("goldenrod","grey2"))+
  theme_classic()

#other way around

bodysizeint %>% 
  ggplot() +
  aes(x = SEDIMENT, y = mean, color = FLOW) +
  geom_line(aes(group = FLOW)) +
  geom_errorbar(limits, width=0.1)+
  geom_point()+
  ylab('Invertebrate body length (mm)')+
  xlab('Sediment')+
  labs(colour = "Flow")+
  scale_colour_manual(values = c("goldenrod2","mediumpurple1"))+
  theme_classic()

bodysizeanova <- aov(Body_size ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(bodysizeanova)
etasq(bodysizeanova)
TukeyHSD(bodysizeanova)

#Orthoclad body size CO2 x Temp

Benth %>% 
  group_by(Temperature, CO2) %>% 
  summarise(mean = mean(Orthoclad_body_size), se = se(Orthoclad_body_size)) -> orthobodysizeint
limits = aes(ymax = mean + se, ymin=mean - se)

orthobodysizeint %>% 
  ggplot() +
  aes(x = Temperature, y = mean, color = CO2) +
  geom_line(aes(group = CO2)) +
  geom_errorbar(limits, width=0.1)+
  geom_point()+
  ylab('Body length (mm)')+
  xlab('Temperature')+
  labs(colour = "CO2")+
  scale_colour_manual(values = c("violet","springgreen3"))+
  theme_classic()

#Orthoclad body size Sediment x Temp

Benth %>% 
  group_by(Temperature, SEDIMENT) %>% 
  summarise(mean = mean(Orthoclad_body_size), se = se(Orthoclad_body_size)) -> orthobodysizeint2
limits = aes(ymax = mean + se, ymin=mean - se)

orthobodysizeint2 %>% 
  ggplot() +
  aes(x = SEDIMENT, y = mean, color = Temperature) +
  geom_line(aes(group = Temperature)) +
  geom_errorbar(limits, width=0.1)+
  geom_point()+
  ylab('Body length (mm)')+
  xlab('Sediment')+
  labs(colour = "Temperature")+
  scale_colour_manual(values = c("skyblue3","tan1"))+
  theme_classic()

orthosizeanova <- aov(Orthoclad_body_size ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(orthosizeanova)
etasq(orthosizeanova)
TukeyHSD(orthosizeanova)




#Annelida CO2 x Sediment

Benth %>% 
  group_by(SEDIMENT, CO2) %>% 
  summarise(mean = mean(Annelida_abundance), se = se(Annelida_abundance)) -> oligint
limits = aes(ymax = mean + se, ymin=mean - se)


oligint %>% 
  ggplot() +
  aes(x = SEDIMENT, y = mean, color = CO2) +
  geom_line(aes(group = CO2)) +
  geom_errorbar(limits, width=0.1)+
  geom_point()+
  ylab('Individuals per channel')+
  labs(colour = "CO2")+
  xlab('Sediment')+
  scale_colour_manual(values = c("violet","springgreen3"))+
  theme_classic()

oligint %>% 
  ggplot() +
  aes(x = CO2, y = mean, color = SEDIMENT) +
  geom_line(aes(group = SEDIMENT)) +
  geom_errorbar(limits, width=0.1)+
  geom_point()+
  ylab('Individuals per channel')+
  labs(colour = "SEDIMENT")+
  xlab('CO2')+
  scale_colour_manual(values = c("goldenrod","grey2"))+
  theme_classic()

sqAnnanova<- aov(sqAnnelida ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqAnnanova)
etasq(sqAnnanova)

TukeyHSD(sqAnnanova)

#Lepto interaction CO2 x Flow

Benth %>% 
  group_by(FLOW, CO2) %>% 
  summarise(mean = mean(Leptophlebiidae_abundance), se = se(Leptophlebiidae_abundance)) -> Leptoint
limits = aes(ymax = mean + se, ymin=mean - se)



Leptoint %>% 
  ggplot() +
  aes(x = CO2, y = mean, color = FLOW) +
  geom_line(aes(group = FLOW)) +
  geom_errorbar(limits, width=0.1)+
  geom_point()+
  ylab('Individuals per channel')+
  labs(colour = "FLOW")+
  xlab('CO2')+
  scale_colour_manual(values = c("goldenrod2","mediumpurple1"))+
  theme_classic()

Leptoint %>% 
  ggplot() +
  aes(x = FLOW, y = mean, color = CO2) +
  geom_line(aes(group = CO2)) +
  geom_errorbar(limits, width=0.1)+
  geom_point()+
  ylab('Individuals per channel')+
  labs(colour = "CO2")+
  xlab('Flow')+
  scale_colour_manual(values = c("violet","springgreen3"))+
  theme_classic()

sqLeptophlebiidaeanova<- aov(sqLeptophlebiidae ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqLeptophlebiidaeanova)
etasq(sqLeptophlebiidaeanova)
TukeyHSD(sqLeptophlebiidaeanova)

#Copepoda 2 way CO2 x Sediment

Benth %>% 
  group_by(SEDIMENT, CO2) %>% 
  summarise(mean = mean(Copepoda_abundance), se = se(Copepoda_abundance)) -> copeint
limits = aes(ymax = mean + se, ymin=mean - se)

copeint

copeint %>% 
  ggplot() +
  aes(x = SEDIMENT, y = mean, color = CO2) +
  geom_line(aes(group = CO2)) +
  geom_errorbar(limits, width=0.1)+
  geom_point()+
  ylab('Individuals per channel')+
  labs(colour = "CO2")+
  xlab('Sediment')+
  scale_colour_manual(values = c("violet","springgreen3"))+
  theme_classic()

sqCopanova<- aov(sqCopepoda ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqCopanova)
etasq(sqCopanova)
TukeyHSD(sqCopanova)


#Oxyethira 3 way Temp x CO2 x Flow

sqOxyethiranova<- aov(sqOxyethira ~ Temperature*CO2*FLOW*SEDIMENT + BLOCK, data = Benth)
summary(sqOxyethiranova)
etasq(sqOxyethiranova)
TukeyHSD(sqOxyethiranova)


groups <- Benth %>% 
  ddply(c("FLOW", "Temperature", "CO2"),
        summarise,
        mean = mean(Oxyethira_abundance),
        se = se(Oxyethira_abundance))


groups<-c(group1,group2)


oxyplot<-ggplot(groups)+
  aes(x = FLOW, y = mean, color = Temperature) +
  geom_line(aes(group = Temperature)) +
  ylim(0, 30) +
  geom_errorbar(limits, width=0.1)+
  geom_point()+
  labs(colour = "Temperature", x= "Flow", y= "Individuals per channel")+
  theme(plot.title = element_text(hjust = -0.5))+
  scale_colour_manual(values = c("skyblue3","tan1"))+
  facet_grid(.~ CO2)+
  theme_classic()
oxyplot

oxy_3way <- ggarrange(left,
                      right,
                       ncol = 2, nrow = 1,                       
                       common.legend = TRUE,
                      hjust = -0.5,
                      legend='right') 

oxy_3way

oxy_3way<-  annotate_figure(oxy_3way,
                              top = text_grob("Sediment added", size = 12),
                              left = text_grob("Invertebrate body length (mm)", size = 12, rot = 90),
                              right = text_grob("Flow treatment", size = 12, rot=270),
                              bottom = text_grob("CO2 treatment", size = 12))





