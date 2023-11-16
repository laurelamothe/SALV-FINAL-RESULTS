library(car)
library(agricolae)
letter_strain = c("A", "B", "C", "D")
list_parameters = c("r", "K", "Kr", "μ", "ϵ", "ϕ", "η", "β")

out=data.frame(param = rep(0,length(list_parameters)*length(letter_strain)), 
               strain = rep(0,length(list_parameters)*length(letter_strain)), 
               anova_p_value = rep(1,length(list_parameters)*length(letter_strain)))
c=0
## ANOVA1
for (param in 1:length(list_parameters)) {
  for (strain in 1:length(letter_strain)) {
      
x=read.csv(paste("./data/basal_model_result_", list_parameters[param], "_" ,letter_strain[strain], "_", "2023-11-09.csv",sep=""))
samp=runif(n=1000, min=1, max=8000)
y=data.frame(salinity = rep(c("S1", "S2", "S3", "S4", "S5", "S6"), each = 1000),
             value = c(x[samp,1], x[samp,2], x[samp,3], x[samp,4], x[samp,5], x[samp,6]))

#print(paste("SHAPIRO", list_parameters[param], letter_strain[strain]))
#print(shapiro.test(x[samp,1])$p.value)
#print(paste("BARTLETT", list_parameters[param], letter_strain[strain]))
#print(bartlett.test(value ~ salinity, data=y)$p.value)
#leveneTest(value ~ salinity, data=y)

#print(paste("ANOVA", list_parameters[param], letter_strain[strain]))
#print(oneway.test(value ~ salinity, data = y, var.equal = FALSE)$p.value)
print(c(list_parameters[param], letter_strain[strain]))
  c = c + 1
out[c,] = c(list_parameters[param], letter_strain[strain],oneway.test(value ~ salinity, data = y, var.equal = FALSE)$p.value )
if (oneway.test(value ~ salinity, data = y, var.equal = FALSE)$p.value < 0.05) {
  #print(paste(list_parameters[param], letter_strain[strain], "significatif"))
  #print(oneway.test(value ~ salinity, data = y, var.equal = FALSE))
}
else {
 # print(paste(list_parameters[param], letter_strain[strain], " NON significatif"))
  
}
  }
}

write.csv(out, "C://Users//33625//OneDrive//PROJET_SAL//SALV-FINAL-RESULTS//Stat_tests//anova.csv" , row.names=FALSE)

##ANOVA 2
NomFichier = read.csv("./data/basal_model_result_Kr_A_2023-11-09.csv", sep=",")

df = data.frame(sal = rep(c("S5", "S10", "S15", "S25", "S35", "S40"), each = 8000), 
value = c(NomFichier$S5,NomFichier$S10,NomFichier$S15,NomFichier$S25,NomFichier$S35,NomFichier$S40))

mod = aov(value~sal , data=df) 
anova(mod)

for (param in 1:length(list_parameters)) {
  for (strain in 1:length(letter_strain)) {
    NomFichier=read.csv(paste("./data/basal_model_result_", list_parameters[param], "_" ,letter_strain[strain], "_", "2023-11-09.csv",sep=""))
    
    df = data.frame(sal = rep(c("S5", "S10", "S15", "S25", "S35", "S40"), each = 8000), 
                    value = c(NomFichier$S5,NomFichier$S10,NomFichier$S15,NomFichier$S25,NomFichier$S35,NomFichier$S40))
    
    mod = aov(value~sal , data=df) 
    print (paste(list_parameters[param],letter_strain[strain], ":" ))
    print(anova(mod)$Pr)

  }}



## TEST COMPARAISON DES MOYENNES

out2 = data.frame(param = rep(0,length(list_parameters)*length(letter_strain)),  
                  strain = rep(0,length(list_parameters)*length(letter_strain)), 
                  S1=rep(0,length(list_parameters)*length(letter_strain)), 
                  S2=rep(0,length(list_parameters)*length(letter_strain)),
                  S3=rep(0,length(list_parameters)*length(letter_strain)),
                  S4=rep(0,length(list_parameters)*length(letter_strain)),
                  S5=rep(0,length(list_parameters)*length(letter_strain)),
                  S6=rep(0,length(list_parameters)*length(letter_strain)))
c = 0
for (param in 1:length(list_parameters)) {
  for (strain in 1:length(letter_strain)) {
    NomFichier=read.csv(paste("./data/basal_model_result_", list_parameters[param], "_" ,letter_strain[strain], "_", "2023-11-09.csv",sep=""))
    df = data.frame(sal = rep(c("S5", "S10", "S15", "S25", "S35", "S40"), each = 8000), 
                    value = c(NomFichier$S5,NomFichier$S10,NomFichier$S15,NomFichier$S25,NomFichier$S35,NomFichier$S40))
    res.aov=aov(df$value~as.factor(df$sal))   #faire l'analyse de variance pl=f(lact)
    d<-df.residual(res.aov)				#garder le df du modèle
    MSerror<-deviance(res.aov)/d				#garder le MSerror du modele
    #car df et MSerror vont servir à la comparaison des moyennes
    hsd.res=HSD.test(df$value,df$sal,d,MSerror, group=TRUE, main="value=f(sal)")
    print (paste(list_parameters[param],letter_strain[strain], ":" ))
    print(hsd.res$groups)
    #comparaison des moyennes
    c = c + 1
    out2[c,] = c(list_parameters[param],letter_strain[strain], t(hsd.res$groups)[2,])
}}
write.csv(out2, "C://Users//33625//OneDrive//PROJET_SAL//SALV-FINAL-RESULTS//Stat_tests//comparaison_moyennes.csv", row.names=FALSE)




## TEST COMPARAISON DES MOYENNES
pl = read.csv("donnees-pl_qui2_an variance_reg lin_corr.csv", sep=";")
res.aov=aov(pl$pl~as.factor(pl$lact))   #faire l'analyse de variance pl=f(lact)
df<-df.residual(res.aov)				#garder le df du modèle
MSerror<-deviance(res.aov)/df				#garder le MSerror du modele
#car df et MSerror vont servir à la comparaison des moyennes
hsd.res=HSD.test(pl$pl,pl$lact,df,MSerror, group=TRUE, main="pl=f(lact)")
#comparaison des moyennes
