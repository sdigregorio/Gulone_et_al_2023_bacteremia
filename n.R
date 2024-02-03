#"Analysis performed in [Gulone et al 2023](https://www.liebertpub.com/doi/10.1089/mdr.2023.0219)"

#loading libraries
library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)
library(ggsci)
library(scales)
library(gridExtra)



#MRSA isolates susceptibility profile

BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(Period)%>%
  count(Bact_isolated)

BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(Period)%>%
  count(MDR)%>%
  mutate(pct=n/sum(n)*100)

MDR_pos = c(27,26,10)
MDR_n = c(44,44,39)
prop.test(MDR_pos,MDR_n)
prop.trend.test(MDR_pos,MDR_n)


#Graph ERY resistance phenotype by period
ery <- BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(ERY_Phenotype)%>%
  ggplot(aes(Period, fill=ERY_Phenotype)) +
  geom_bar(position="fill")+
  theme_minimal()+
  theme(title = element_text(size=7, face="bold"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        legend.text = element_text(size = 8),
        legend.title = element_blank())+
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values = c("I"="#f2931f",
                               "R"="#E64B35FF", 
                               "S"="#4DBBD5FF"))+
  labs(y="Relative frequencies",
       title="Erythromycin")

#Graph CLIN resistance phenotype by period
clin<-BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(CLIN_Phenotype)%>%
  ggplot(aes(Period, fill=CLIN_Phenotype)) +
  geom_bar(position="fill")+
  theme_minimal()+
  theme(title = element_text(size=7, face="bold"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        legend.text = element_text(size = 8),
        legend.title = element_blank())+
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values = c("I"="#f2931f",
                               "R"="#E64B35FF", 
                               "S"="#4DBBD5FF"))+
  labs(y="Relative frequencies",
       title="Clindamycin")

#Graph GEN resistance phenotype by period
gen<-BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(GEN_Phenotype)%>%
  ggplot(aes(Period, fill=GEN_Phenotype)) +
  geom_bar(position="fill")+
  theme_minimal()+
  theme(title = element_text(size=7, face="bold"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        legend.text = element_text(size = 8),
        legend.title = element_blank())+
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values = c("I"="#f2931f",
                               "R"="#E64B35FF", 
                               "S"="#4DBBD5FF"))+
  labs(y="Relative frequencies",
       title="Gentamycin")


#Graph tet (tet+min) resistance phenotype by period
tet<-BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(Tetraclyclines_phenotype)%>%
  ggplot(aes(Period, fill=Tetraclyclines_phenotype)) +
  geom_bar(position="fill")+
  theme_minimal()+
  theme(title = element_text(size=7, face="bold"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        legend.text = element_text(size = 8),
        legend.title = element_blank())+
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values = c("I"="#f2931f",
                               "R"="#E64B35FF", 
                               "S"="#4DBBD5FF"))+
  labs(y="Relative frequencies",
       title="Tetracyclines")


#Graph RIF resistance phenotype by period
rif<-BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(RIF_Phenotype)%>%
  ggplot(aes(Period, fill=RIF_Phenotype)) +
  geom_bar(position="fill")+
  theme_minimal()+
  theme(title = element_text(size=7, face="bold"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        legend.text = element_text(size = 8),
        legend.title = element_blank())+
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values = c("I"="#f2931f",
                               "R"="#E64B35FF", 
                               "S"="#4DBBD5FF"))+
  labs(y="Relative frequencies",
       title="Rifampin")


#Graph fluoroquinolones (CIP+LEV) resistance phenotype by period
fq<-BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(Fluorquinolones_phenotype)%>%
  ggplot(aes(Period, fill=Fluorquinolones_phenotype)) +
  geom_bar(position="fill")+
  theme_minimal()+
  theme(title = element_text(size=7, face="bold"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        legend.text = element_text(size = 8),
        legend.title = element_blank())+
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values = c("I"="#f2931f",
                               "R"="#E64B35FF", 
                               "S"="#4DBBD5FF"))+
  labs(y="Relative frequencies",
       title="Fluorquinolones")


#Graph tms resistance phenotype by period
tms<-BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(TMS_Phenotype)%>%
  ggplot(aes(Period, fill=TMS_Phenotype)) +
  geom_bar(position="fill")+
  theme_minimal()+
  theme(title = element_text(size=7, face="bold"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        legend.text = element_text(size = 8),
        legend.title = element_blank())+
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values = c("I"="#f2931f",
                               "R"="#E64B35FF", 
                               "S"="#4DBBD5FF"))+
  labs(y="Relative frequencies",
       title="Trimethoprim-sulfamethoxazole")

#Graph VAN resistance phenotype by period
van<-BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(VAN_Phenotype)%>%
  ggplot(aes(Period, fill=VAN_Phenotype)) +
  geom_bar(position="fill")+
  theme_minimal()+
  theme(title = element_text(size=7, face="bold"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        legend.text = element_text(size = 8),
        legend.title = element_blank())+
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values = c("I"="#f2931f",
                               "R"="#E64B35FF", 
                               "S"="#4DBBD5FF"))+
  labs(y="Relative frequencies",
       title="Vancomycin")

#Graph LZD resistance phenotype by period
lzd<-BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(LZD_Phenotype)%>%
  ggplot(aes(Period, fill=LZD_Phenotype)) +
  geom_bar(position="fill")+
  theme_minimal()+
  theme(title = element_text(size=7, face="bold"),
        axis.text=element_text(size=7),
        axis.title=element_text(size=7,face="bold"),
        legend.text = element_text(size = 8),
        legend.title = element_blank())+
  scale_y_continuous(labels=scales::percent) +
  scale_fill_manual(values = c("I"="#f2931f",
                               "R"="#E64B35FF", 
                               "S"="#4DBBD5FF"))+
  labs(y="Relative frequencies",
       title="Linezolid")

#join graphs into one
grid.arrange(ery, clin, fq, gen,  rif, tet, tms, van, lzd, 
             nrow = 3,
             ncol = 3)



#Count Erythromycin resistance phenotype by period
BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(Period)%>%
  count(ERY_Phenotype)%>%
  mutate(pct=n/sum(n)*100)

#define ERY R positives and negatives and perform prop test
ERY_R_pos = c(26,20,10)
ERY_R_n = c(44,44,39)
prop.test(ERY_R_pos,ERY_R_n)
prop.trend.test(ERY_R_pos,ERY_R_n)

#Count clindamycin resistance phenotype by period
BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(Period)%>%
  count(CLIN_Phenotype)%>%
  mutate(pct=n/sum(n)*100)

#define CLIN R positives and negatives and perform prop test
CLIN_R_pos = c(20,17,3)
CLIN_R_n = c(44,44,39)
prop.test(CLIN_R_pos,CLIN_R_n)
prop.trend.test(CLIN_R_pos,CLIN_R_n)

#Count fluoroquinolones resistance phenotype by period
BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
   group_by(Period)%>%
   count(Fluorquinolones_phenotype)

#define FQ R positives and negatives and perform prop test
FQ_R_pos = c(26,25,8)
FQ_R_n = c(44,44,39)
prop.test(FQ_R_pos,FQ_R_n)
prop.trend.test(FQ_R_pos,FQ_R_n)



#Molecular characterization
#Graph SCCmec by period
SCCmec_period<-BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(SCCmec)%>%
  ggplot(aes(Period, fill=SCCmec)) +
  geom_bar(position="fill")+
  theme_minimal()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))+
  scale_fill_jama()+
  scale_y_continuous(labels=scales::percent) +
  labs(y="Relative frequencies")

SCCmec_period

#Count SCCmec by period
BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(Period)%>%
  count(SCCmec)%>%
  mutate(pct=n/sum(n)*100)

#Perform prop test
SCCmecIV_pos = c(17,26,34)
SCCmecIV_n = c(44,44,39)
prop.test(SCCmecIV_pos,SCCmecIV_n)
prop.trend.test(SCCmecIV_pos,SCCmecIV_n)

#Graph MDR by SCCmec type
SCCmec_A<-BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(MDR)%>%
  ggplot(aes(SCCmec, fill=MDR)) +
  geom_bar(position="fill")+
  theme_minimal()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))+
  scale_fill_lancet()+
  scale_y_continuous(labels=scales::percent) +
  labs(y="Relative frequencies")
SCCmec_A


#Count SCCmecIV isolates with PVL
BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  filter(SCCmec=="IV")%>%
  group_by(Period)%>%
  count(PVL)%>%
  mutate(pct=n/sum(n)*100)

#Perform prop test
pvl_pos = c(13,13,16)
pvl_n = c(17,26,34)
prop.test(pvl_pos,pvl_n)
prop.trend.test(pvl_pos,pvl_n)


#Graph Genotype by period
Genotype_period<-BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(Genotype)%>%
  ggplot(aes(Period, fill=Genotype)) +
  geom_bar(position="fill")+
  theme_minimal()+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))+
  scale_fill_simpsons()+
  scale_y_continuous(labels=scales::percent) +
  labs(y="Relative frequencies",
       title="B")

Genotype_period


#Graph spa type by period and Genotype
spageno2<-BACTERIEMIAS_SAMR_2009_2016_MLS.mod%>%
  group_by(spa_type)%>%
  filter(Period!="P1")
spageno2$Genotype <- factor(spageno2$Genotype,
                              levels = c("ND-MRSA-I",
                                         "ST100-MRSA-V",
                                         "ST238-MRSA-III",
                                         "ST72-MRSA-I",
                                         "ND-MRSA-NT",
                                         "ST5-MRSA-I",
                                         "ND-MRSA-IV",
                                         "ST100-MRSA-IV",
                                         "ST72-MRSA-IV",
                                         "ST97-MRSA-IV",
                                         "ST5-MRSA-IV",
                                         "ST8-MRSA-IV",
                                         "ST30-MRSA-IV"
                              ))

spageno3<-  ggplot(spageno2, aes(x=Genotype, fill=spa_type)) +
  geom_bar(position="fill",
           width = 0.8)+
  theme_minimal()+
  theme(axis.text = element_text(size = 8),
        axis.title = element_text(size = 9))+
  scale_fill_manual(values = c("ND"="#D3D3D3",
                               "Not in DB"="#DDA0DD",
                               "Novel-t008 related"="#A9A9A9FF",
                               "Novel-t019 related"="#F0E685FF",
                               "Novel-t11608 related"="#466983FF",
                               "Novel-t1346 related"="#7FFFD4",
                               "Novel-t149 related"="#D595A7FF",
                               "Novel-t2724 related"="#802268FF",
                               "Novel-t512 related"="#6BD76BFF",
                               "Novel-t942 related"="#5DB1DDFF",
                               "t002"="#924822FF",
                               "t008"="#837B8DFF",
                               "t019"="#e25322",
                               "t024"="#FFC0CB",
                               "t037"="#7A65A5FF",
                               "t067"="#E4AF69FF",
                               "t11644"="#3B1B53FF",
                               "t121"="#CDDEB7FF",
                               "t148"="#612A79FF",
                               "t149"="#AE1F63FF",
                               "t168"="#fe935a",
                               "t189"="#5A655EFF",
                               "t2051"="#CC9900FF",
                               "t215"="#99CC00FF",
                               "t2168"="#749B58FF",
                               "t311"="#303386",
                               "t330"="#FBEF3A",
                               "t359"="#33CC00FF",
                               "t5160"="#008000",
                               "t521"="#00CC99FF",
                               "t7892"="#0099CCFF",
                               "t9645"="#0A47FFFF"))+
  scale_y_continuous(labels=scales::percent) +
  labs(y="Relative frequencies",
       fill="spa type")
                                        

spageno3 +facet_wrap(~Period)+coord_flip()+
  theme(strip.text.x = element_text(size = 9))


