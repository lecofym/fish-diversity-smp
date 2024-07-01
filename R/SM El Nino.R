###El Nino Supplementary material####
# S1 Text: Table C.
simpson$transect <- rownames(simpson)
test <- simpson[simpson$transect %in% rownames(PNH1), 5]
SM_Nino<- cbind(PNH1, test)
SM_Nino<- SM_Nino %>% rename(simpson = test)
names(SM_Nino)
SM_Nino$Year<- as.factor(SM_Nino$Year)

NP <- subset(SM_Nino, SM_Nino$Human_use_level == "Non-protected")
MPA <- subset(SM_Nino, SM_Nino$Human_use_level == "MPA")

#NP El Nino effect
NP.lm.S <- lm(Nb_sp ~ Year,data= NP)
anova(NP.lm.S)
summary(NP.lm.S)
S.year.boxplot<- ggplot(NP,aes(x=Year,y=Nb_sp, fill= Human_use_level))+ geom_boxplot()+ ylab("S")+ theme_bw()+ theme(panel.grid=element_blank(), panel.background=element_rect(size=2,color="black"), axis.title.x=element_blank(), axis.title.y=element_text(size=14),axis.text.x=element_blank(), axis.text.y=element_text(size=12), legend.position="none")
S.year.boxplot

NP.lm.d <- lm(log(Tot_weight,base=2) ~ Year,data= NP)
anova(NP.lm.d)
summary(NP.lm.d)
d.year.boxplot<- ggplot(NP,aes(x=Year,y=log(Tot_weight,base=2), fill= Human_use_level))+ geom_boxplot()+ ylab("d (log2)")+ theme_bw()+ theme(panel.grid=element_blank(), panel.background=element_rect(size=2,color="black"), axis.title.x=element_blank(), axis.title.y=element_text(size=14),axis.text.x=element_blank(), axis.text.y=element_text(size=12), legend.position="none")
d.year.boxplot

NP.lm.D <- lm(simpson ~ Year,data= NP)
anova(NP.lm.D)
summary(NP.lm.D)
D.year.boxplot<- ggplot(NP,aes(x=Year,y=simpson, fill= Human_use_level))+ geom_boxplot()+ ylab("D")+ theme_bw()+ theme(panel.grid=element_blank(), panel.background=element_rect(size=2,color="black"), axis.title.x=element_blank(), axis.title.y=element_text(size=14),axis.text.x=element_blank(), axis.text.y=element_text(size=12), legend.position="none")
D.year.boxplot

NP.lm.FR <- lm(FRic ~ Year, data= NP)
anova(NP.lm.FR)
summary(NP.lm.FR)
FR.year.boxplot<- ggplot(NP,aes(x=Year,y=FRic, fill= Human_use_level))+ geom_boxplot()+ ylab("FRic")+ theme_bw()+ theme(panel.grid=element_blank(), panel.background=element_rect(size=2,color="black"), axis.title.x=element_blank(), axis.title.y=element_text(size=14),axis.text.x=element_blank(), axis.text.y=element_text(size=12), legend.position="none")
FR.year.boxplot

NP.lm.FD <-lm(FDiv ~ Year, data= NP)
anova(NP.lm.FD)
summary(NP.lm.FD)
FD.year.boxplot<- ggplot(NP,aes(x=Year,y=FDiv, fill= Human_use_level))+ geom_boxplot()+ ylab("FDiv")+ theme_bw()+ theme(panel.grid=element_blank(), panel.background=element_rect(size=2,color="black"), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=12, angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=12), legend.position="none")
FD.year.boxplot

NP.lm.FO <- lm(FOri ~ Year, data= NP)
anova(NP.lm.FO)
summary(NP.lm.FO)
FO.year.boxplot<- ggplot(NP,aes(x=Year,y=FOri, fill= Human_use_level))+ geom_boxplot()+ ylab("FOri")+ theme_bw()+ theme(panel.grid=element_blank(), panel.background=element_rect(size=2,color="black"), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=12,angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=12), legend.position="none")
FO.year.boxplot

NP.Nino.graph <- S.year.boxplot + d.year.boxplot + D.year.boxplot + FR.year.boxplot + FD.year.boxplot + FO.year.boxplot + plot_layout(nrow = 3, guides = "collect")

ggsave('Figure D.tiff', plot = Spatial.graph,
       width = 2250, height = 2000, units = 'px', dpi = 320,
       compression = "lzw")

ggsave('Figure D.jpeg', plot = NP.Nino.graph,
       width = 2250, height = 2000, units = 'px', dpi = 320)


#MPA El Nino effect
MPA.lm.S <- lm(Nb_sp ~ Year,data= MPA)
anova(MPA.lm.S)
summary(MPA.lm.S)
MPA.S.year.boxplot<- ggplot(MPA,aes(x=Year,y=Nb_sp, fill= Human_use_level))+ geom_boxplot()+ ylab("S")+ theme_bw() + scale_fill_manual(values = cpalette) + theme(panel.grid=element_blank(), panel.background=element_rect(size=2,color="black"), axis.title.x=element_blank(), axis.title.y=element_text(size=14),axis.text.x=element_blank(), axis.text.y=element_text(size=12), legend.position="none")
MPA.S.year.boxplot

MPA.lm.d <- lm(log(Tot_weight,base=2) ~ Year,data= MPA)
anova(MPA.lm.d)
summary(MPA.lm.d)
MPA.d.year.boxplot<- ggplot(MPA,aes(x=Year,y=log(Tot_weight,base=2), fill= Human_use_level))+ geom_boxplot()+ ylab("d (log2)")+ theme_bw() + scale_fill_manual(values = cpalette) + theme(panel.grid=element_blank(), panel.background=element_rect(size=2,color="black"), axis.title.x=element_blank(), axis.title.y=element_text(size=14),axis.text.x=element_blank(), axis.text.y=element_text(size=12), legend.position="none")
MPA.d.year.boxplot

MPA.lm.D <- lm(simpson ~ Year,data= MPA)
anova(MPA.lm.D)
summary(MPA.lm.D)
MPA.D.year.boxplot<- ggplot(MPA,aes(x=Year,y=simpson, fill= Human_use_level))+ geom_boxplot()+ ylab("D")+ theme_bw()+ scale_fill_manual(values = cpalette) + theme(panel.grid=element_blank(), panel.background=element_rect(size=2,color="black"), axis.title.x=element_blank(), axis.title.y=element_text(size=14),axis.text.x=element_blank(), axis.text.y=element_text(size=12), legend.position="none")
MPA.D.year.boxplot

MPA.lm.FR <- lm(FRic ~ Year, data= MPA)
anova(MPA.lm.FR)
summary(MPA.lm.FR)
MPA.FR.year.boxplot<- ggplot(MPA,aes(x=Year,y=FRic, fill= Human_use_level))+ geom_boxplot()+ ylab("FRic")+ theme_bw()+ scale_fill_manual(values = cpalette) + theme(panel.grid=element_blank(), panel.background=element_rect(size=2,color="black"), axis.title.x=element_blank(), axis.title.y=element_text(size=14),axis.text.x=element_blank(), axis.text.y=element_text(size=12), legend.position="none")
MPA.FR.year.boxplot

MPA.lm.FD <-lm(FDiv ~ Year, data= MPA)
anova(MPA.lm.FD)
summary(MPA.lm.FD)
MPA.FD.year.boxplot<- ggplot(MPA,aes(x=Year,y=FDiv, fill= Human_use_level))+ geom_boxplot()+ ylab("FDiv")+ theme_bw()+ scale_fill_manual(values = cpalette) + theme(panel.grid=element_blank(), panel.background=element_rect(size=2,color="black"), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=12, angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=12), legend.position="none")
MPA.FD.year.boxplot

MPA.lm.FO <- lm(FOri ~ Year, data= MPA)
anova(MPA.lm.FO)
summary(MPA.lm.FO)
MPA.FO.year.boxplot<- ggplot(MPA,aes(x=Year,y=FOri, fill= Human_use_level))+ geom_boxplot()+ ylab("FOri")+ theme_bw()+ scale_fill_manual(values = cpalette) + theme(panel.grid=element_blank(), panel.background=element_rect(size=2,color="black"), axis.title.x=element_text(size=14), axis.title.y=element_text(size=14),axis.text.x=element_text(size=12,angle = 90, vjust = 0.5, hjust=1), axis.text.y=element_text(size=12), legend.position="none")
MPA.FO.year.boxplot

MPA.Nino.graph <- MPA.S.year.boxplot + MPA.d.year.boxplot + MPA.D.year.boxplot + MPA.FR.year.boxplot + MPA.FD.year.boxplot + MPA.FO.year.boxplot + plot_layout(nrow = 3, guides = "collect")

ggsave('Figure C.tiff', plot = MPA.Nino.graph,
       width = 2250, height = 2000, units = 'px', dpi = 320,
       compression = "lzw")

ggsave('Figure C.jpeg', plot = MPA.Nino.graph,
       width = 2250, height = 2000, units = 'px', dpi = 320)
