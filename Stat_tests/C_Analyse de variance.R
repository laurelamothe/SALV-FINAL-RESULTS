library(car)

NomFichier = read.csv("donnees-pl_qui2_an variance_reg lin_corr.csv", sep=";")

############ fa?on 1 - ANOVA lot ###################

mod = aov(pl~lact , data=NomFichier) 
summary (mod)
anova(mod)
############ fa?on 2 - ANOVA lot*lact ###################
mod<-lm(as.numeric(NomFichier$Vexpliquee)~as.factor(NomFichier$Vexplicative1)*as.factor(NomFichier$Vexplicative2)) #OK
Anova(mod, type=c(2))
summary (mod)


########## fa?on julien breda##############
# faire d'abord CTRL+C de la plage de données dans le fichier excel à importer
NomFichier<-read.table(file("clipboard"),sep="\t",h=T, dec=".")#si fichier xls
lm1<-lm(NomFichier$ps ~ NomFichier$lot)
anova(lm1)


####################################################
########## exemple autonome #######################
############ Importation de données   #############
#chosir le fichier donnees-pl.csv#
pl<-read.table(choose.files(), h=T, sep=";")
pl$lact=as.factor(pl$lact)


############ façon 1 - ANOVA lot ###################
aov.pl = aov(pl~lot , data=pl) 
summary (aov.pl) 


############ façon 2 - ANOVA lot*lact ###################
mod<-lm(as.numeric(pl$pl)~as.factor(pl$lot)*as.factor(pl$lact))
Anova(mod, type=c(2))
summary (mod)
###################################################

