setwd("C:/Users/pablo/OneDrive - unizar.es/Documentos/simul_houssem")

library(MASS)
library(MBESS)
library(pedigree)
library(AlphaPart)
library(dplyr)
library(corrplot)
library(ggplot2)
library(patchwork)

#varcovar<-matrix(0,8,8)
varE<-c(0.127,0.244,0.259,0.349)
varR<-c(0.137,0.190,0.159,0.334) #PPE PME CPE CME PPR PMR CPR CMR

corE<-matrix(0,4,4)
corE[lower.tri(corE)]<-c(0.52,0.53,0.4,0.21,0.32,0.4)
#corE[upper.tri(corE)]<-corE[lower.tri(corE)]

corR<-matrix(0,4,4)
corR[lower.tri(corR)]<-c(0.09,0.17,0.04,0.46,0.41,0.47)
#corR[upper.tri(corR)]<-corR[lower.tri(corR)]

diag(corE)<-1
diag(corR)<-1

varcovarE<-cor2cov(corE,sqrt(varE))
varcovarR<-cor2cov(corR,sqrt(varR))

varcovar0<-matrix(0,4,4)

varcovartot<-rbind(cbind(varcovarE,varcovar0),cbind(varcovar0,varcovarR))
varcovartot<-matrix(Matrix::forceSymmetric(varcovartot,uplo="L"),8,8)

medias<-c(rep(0,8))

geneaEE<-data.frame(V1=c(1:1000),V2=0,V3=0,V4=c(1:2),Raza="EE")
geneaRR<-data.frame(V1=c(1001:2000),V2=0,V3=0,V4=c(1:2),Raza="RR")
geneaER<-data.frame(V1=c(2001:2500),V2=0,V3=0,V4=c(2),Raza="ER")
geneaRE<-data.frame(V1=c(2501:3000),V2=0,V3=0,V4=c(2),Raza="RE")

geneaER$V2<-sample(geneaEE$V1[geneaEE$V4==1],nrow(geneaER),replace = T)
geneaER$V3<-sample(geneaRR$V1[geneaRR$V4==2],nrow(geneaER),replace = T)

geneaRE$V3<-sample(geneaEE$V1[geneaEE$V4==2],nrow(geneaER),replace = T)
geneaRE$V2<-sample(geneaRR$V1[geneaRR$V4==1],nrow(geneaER),replace = T)

geneatot<-rbind(geneaEE,geneaRR,geneaER,geneaRE)

TBVE<-as.data.frame(mvrnorm(nrow(geneaEE)*2,medias,varcovartot))
TBVR<-as.data.frame(mvrnorm(nrow(geneaRR)*2,medias,varcovartot))
TBVtot<-rbind(TBVE,TBVR)

var(TBVE)
var(TBVR)
corrplot.mixed(cor(TBVE))
corrplot.mixed(cor(TBVR))
corrplot.mixed(cor(TBVtot))

TBVER<-as.data.frame(matrix(0,nrow = nrow(geneaER)*2,ncol = 8))
TBVRE<-as.data.frame(matrix(0,nrow = nrow(geneaRE)*2,ncol = 8))

filas<-seq_len(nrow(TBVER)) %% 2
TBVER[filas==1,]<-TBVtot[(geneaER$V2*2)-1,]*0.5+TBVtot[geneaER$V2*2,]*0.5
TBVER[filas==0,]<-TBVtot[(geneaER$V3*2)-1,]*0.5+TBVtot[geneaER$V3*2,]*0.5


filas<-seq_len(nrow(TBVRE)) %% 2
TBVRE[filas==1,]<-TBVtot[(geneaRE$V2*2)-1,]*0.5+TBVtot[geneaRE$V2*2,]*0.5
TBVRE[filas==0,]<-TBVtot[(geneaRE$V3*2)-1,]*0.5+TBVtot[geneaRE$V3*2,]*0.5

TBVER<-TBVER+(mvrnorm(nrow(TBVER),medias,varcovartot/2))#les sumo el muestreo mendeliano porque aunque no son reproductores, 
TBVRE<-TBVRE+(mvrnorm(nrow(TBVRE),medias,varcovartot/2))#pueden variar la varianza de la población

corrplot.mixed(cor(TBVER))
corrplot.mixed(cor(TBVRE))

TBVtot<-rbind(TBVtot,TBVER,TBVRE)

corrplot.mixed(cor(TBVtot))

var(TBVER)
var(TBVRE)

#generaciones sin selección

for (Generation in 1:20) {
  
  print(Generation)
  
  newgeneaEE<-data.frame(V1=c(1:1000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="EE")
  newgeneaRR<-data.frame(V1=c(1001:2000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="RR")
  newgeneaER<-data.frame(V1=c(2001:2500)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="ER")
  newgeneaRE<-data.frame(V1=c(2501:3000)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="RE")
  
  #escoger los reproductores EE
  machosEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==1], #TBV machos EE
                       V2=TBVtot$V1[(geneaEE$V1[geneaEE$V4==1]*2)-1]+
                         TBVtot$V1[(geneaEE$V1[geneaEE$V4==1]*2)])
  machosEE<-machosEE$V1[order(machosEE$V2,decreasing = T)]
  
  hembrasEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==2], #TBV hembras EE
                        V2=TBVtot$V2[(geneaEE$V1[geneaEE$V4==2]*2)]+
                          TBVtot$V2[(geneaEE$V1[geneaEE$V4==2]*2)-1])
  hembrasEE<-hembrasEE$V1[order(hembrasEE$V2,decreasing = T)]
  
  newgeneaEE$V2<-sample(machosEE,nrow(newgeneaEE),replace = T) #cojo los 100 mejores machos de 500
  newgeneaEE$V3<-sample(hembrasEE,nrow(newgeneaEE),replace = T) #cojo las 100 mejores hembras de 500
  
  
  #escoger los reproductores RR
  machosRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==1], #TBV machos RR
                       V2=TBVtot$V5[geneaRR$V1[geneaRR$V4==1]]+
                         TBVtot$V6[geneaRR$V1[geneaRR$V4==1]])
  machosRR<-machosRR$V1[order(machosRR$V2,decreasing = T)]
  
  hembrasRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==2], #TBV hembras RR
                        V2=TBVtot$V5[geneaRR$V1[geneaRR$V4==2]]+
                          TBVtot$V6[geneaRR$V1[geneaRR$V4==2]])
  hembrasRR<-hembrasRR$V1[order(hembrasRR$V2,decreasing = T)]
  
  newgeneaRR$V2<-sample(machosRR,nrow(newgeneaRR),replace = T) #cojo los 100 mejores machos de 500
  newgeneaRR$V3<-sample(hembrasRR,nrow(newgeneaRR),replace = T) #cojo las 100 mejores hembras de 500
  
  #animales cruzados
  
  newgeneaER$V2<-sample(newgeneaEE$V1[newgeneaEE$V4==1],nrow(newgeneaER),replace = T)
  newgeneaER$V3<-sample(newgeneaRR$V1[newgeneaRR$V4==2],nrow(newgeneaER),replace = T)
  
  newgeneaRE$V2<-sample(newgeneaRR$V1[newgeneaRR$V4==2],nrow(newgeneaER),replace = T)
  newgeneaRE$V3<-sample(newgeneaEE$V1[newgeneaEE$V4==1],nrow(newgeneaER),replace = T)
  
  geneaEE<-newgeneaEE
  geneaRR<-newgeneaRR
  geneaER<-newgeneaER
  geneaRE<-newgeneaRE
  
  #genea completa
  
  geneatot<-rbind(geneatot,newgeneaEE,newgeneaRR,newgeneaER,newgeneaRE)
  
  TBVEnew<-as.data.frame(matrix(0,nrow = nrow(geneaEE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVEnew)) %% 2
  TBVEnew[filas==1,]<-TBVtot[(geneaEE$V2*2)-1,]*0.5+TBVtot[geneaEE$V2*2,]*0.5
  TBVEnew[filas==0,]<-TBVtot[(geneaEE$V3*2)-1,]*0.5+TBVtot[geneaEE$V3*2,]*0.5
  TBVEnew<-TBVEnew+(mvrnorm(nrow(geneaEE)*2,medias,varcovartot/2))
  
  TBVRnew<-as.data.frame(matrix(0,nrow = nrow(geneaRR)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRnew)) %% 2
  TBVRnew[filas==1,]<-TBVtot[(geneaRR$V2*2)-1,]*0.5+TBVtot[geneaRR$V2*2,]*0.5
  TBVRnew[filas==0,]<-TBVtot[(geneaRR$V3*2)-1,]*0.5+TBVtot[geneaRR$V3*2,]*0.5
  TBVRnew<-TBVRnew+(mvrnorm(nrow(geneaRR)*2,medias,varcovartot/2))
  
  TBVE<-TBVEnew
  TBVR<-TBVRnew
  
  TBVtot<-rbind(TBVtot,TBVE,TBVR)
  
  TBVER<-as.data.frame(matrix(0,nrow = nrow(geneaER)*2,ncol = 8))
  filas<-seq_len(nrow(TBVER)) %% 2
  TBVER[filas==1,]<-TBVtot[(geneaER$V2*2)-1,]*0.5+TBVtot[geneaER$V2*2,]*0.5
  TBVER[filas==0,]<-TBVtot[(geneaER$V3*2)-1,]*0.5+TBVtot[geneaER$V3*2,]*0.5
  TBVER<-TBVER+(mvrnorm(nrow(geneaER)*2,medias,varcovartot/2))
  
  TBVRE<-as.data.frame(matrix(0,nrow = nrow(geneaRE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRE)) %% 2
  TBVRE[filas==1,]<-TBVtot[(geneaRE$V2*2)-1,]*0.5+TBVtot[geneaRE$V2*2,]*0.5
  TBVRE[filas==0,]<-TBVtot[(geneaRE$V3*2)-1,]*0.5+TBVtot[geneaRE$V3*2,]*0.5
  TBVRE<-TBVRE+(mvrnorm(nrow(geneaRE)*2,medias,varcovartot/2))
  
  TBVtot<-rbind(TBVtot,TBVER,TBVRE)
  
}

var(TBVtot)
corrplot.mixed(cor(TBVtot))
TBVtot2<-TBVtot
colnames(TBVtot2)<-as.character(round(diag(varcovartot),2))
corrplot.mixed(cor(TBVtot2))

#datos
geneaEEtot<-geneatot[geneatot$Raza=="EE",]
geneaRRtot<-geneatot[geneatot$Raza=="RR",]
geneaERtot<-geneatot[geneatot$Raza=="ER",]
geneaREtot<-geneatot[geneatot$Raza=="RE",]

datosEEmachos<-data.frame(V1=(geneaEEtot$V1[geneaEEtot$V4==1]*2)-1,V2=geneaEEtot$V1[geneaEEtot$V4==1]*2,V3=1,V4=10,V5=0,V6=0,V7=0)
datosEEmachos$V4<-TBVtot$V1[datosEEmachos$V1]+TBVtot$V2[datosEEmachos$V2]

datosEEhembras<-data.frame(V1=(geneaEEtot$V1[geneaEEtot$V4==2]*2)-1,V2=geneaEEtot$V1[geneaEEtot$V4==2]*2,V3=1,V4=10,V5=0,V6=0,V7=0)
datosEEhembras$V4<-TBVtot$V1[datosEEhembras$V1]+TBVtot$V2[datosEEhembras$V2]

datosRRmachos<-data.frame(V1=(geneaRRtot$V1[geneaRRtot$V4==1]*2)-1,V2=geneaRRtot$V1[geneaRRtot$V4==1]*2,V3=1,V4=0,V5=10,V6=0,V7=0)
datosRRmachos$V5<-TBVtot$V5[datosRRmachos$V1]+TBVtot$V6[datosRRmachos$V2]

datosRRhembras<-data.frame(V1=(geneaRRtot$V1[geneaRRtot$V4==2]*2)-1,V2=geneaRRtot$V1[geneaRRtot$V4==2]*2,V3=1,V4=0,V5=10,V6=0,V7=0)
datosRRhembras$V5<-TBVtot$V5[datosRRhembras$V1]+TBVtot$V6[datosRRhembras$V2]

datosER<-data.frame(V1=(geneaERtot$V1*2)-1,V2=geneaERtot$V1*2,V3=1,V4=0,V5=0,V6=10,V7=0) #solo hay hembras
datosER$V6<-TBVtot$V3[datosER$V1]+TBVtot$V8[datosER$V2]

datosRE<-data.frame(V1=(geneaREtot$V1*2)-1,V2=geneaREtot$V1*2,V3=1,V4=0,V5=0,V6=0,V7=10) #solo hay hembras
datosRE$V7<-TBVtot$V7[datosRE$V1]+TBVtot$V4[datosRE$V2]

datostot<-rbind(datosEEmachos,datosEEhembras,datosRRmachos,datosRRhembras,datosER,datosRE)

write.table(datostot,"datos.txt",col.names = F,row.names = F,quote = F)


geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(c(1:21),each=1000))
geneaEEprint$TBVmatER<-TBVtot$V4[(geneaEEprint$V1*2)-1]+TBVtot$V4[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(c(1:21),each=1000))
geneaRRprint$TBVmatER<-TBVtot$V8[(geneaRRprint$V1*2)-1]+TBVtot$V8[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatER",colBy = "Gen")
plot(AlphaPart,ylim = c(-1,10),color = c(2,3))

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(c(1:21),each=1000))
geneaEEprint$TBVmatRE<-TBVtot$V3[(geneaEEprint$V1*2)-1]+TBVtot$V3[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(c(1:21),each=1000))
geneaRRprint$TBVmatRE<-TBVtot$V7[(geneaRRprint$V1*2)-1]+TBVtot$V7[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatRE",colBy = "Gen")
plot(AlphaPart,ylim = c(-1,10),color = c(2,3))


round(varcovartot,4)

save.image("pobbase.RData")

#BLUP

BLUP<-function(genea,datos){
  genea$V4<-calcInbreeding(genea)
  write.table(genea,"genea.txt",quote = F,row.names = F,col.names = F)
  system(paste('echo ',nrow(genea), '|./gam'))
  write.table(datos,'datos.txt',quote = F,row.names = F,col.names = F)
  #system('ulimit -s unlimited')
  system(paste('echo ',nrow(genea)*2, '|./doparam'))
  system('blupf90+ param.txt')
}

GIBBS<-function(genea,datos){
  genea$V4<-calcInbreeding(genea)
  write.table(genea,"genea.txt",quote = F,row.names = F,col.names = F)
  system(paste('echo ',nrow(genea), '|./gam'))
  write.table(datos,'datos.txt',quote = F,row.names = F,col.names = F)
  #system('ulimit -s unlimited')
  system('gibbsf90+ param.txt --samples 50000 --burnin 20000 --interval 1')
}

#####
#EBV#
#####

#seleccionar por madres puras#

load("pobbase.RData")

for (Generation in 21:30) {
  BLUP(geneatot,datostot)
  solutions<-read.table('solutions',skip=1)
  
  EBVtot<-data.frame(V1=solutions$V4[solutions$V1==1 & solutions$V2==2],
                     V2=solutions$V4[solutions$V1==1 & solutions$V2==3],
                     V3=solutions$V4[solutions$V1==3 & solutions$V2==2],
                     V4=solutions$V4[solutions$V1==4 & solutions$V2==3],
                     V5=solutions$V4[solutions$V1==2 & solutions$V2==2],
                     V6=solutions$V4[solutions$V1==2 & solutions$V2==3],
                     V7=solutions$V4[solutions$V1==4 & solutions$V2==2],
                     V8=solutions$V4[solutions$V1==3 & solutions$V2==3])
  
  newgeneaEE<-data.frame(V1=c(1:1000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="EE")
  newgeneaRR<-data.frame(V1=c(1001:2000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="RR")
  newgeneaER<-data.frame(V1=c(2001:2500)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="ER")
  newgeneaRE<-data.frame(V1=c(2501:3000)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="RE")
  
  #escoger los reproductores EE
  machosEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==1]) 
  machosEE$V2<-EBVtot$V2[(machosEE$V1)*2-1]+EBVtot$V2[(machosEE$V1)*2]
  machosEE<-machosEE$V1[order(machosEE$V2,decreasing = T)]
  
  hembrasEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==2])
  hembrasEE$V2<-EBVtot$V2[(hembrasEE$V1)*2-1]+EBVtot$V2[(hembrasEE$V1)*2]
  hembrasEE<-hembrasEE$V1[order(hembrasEE$V2,decreasing = T)]
  
  newgeneaEE$V2<-sample(machosEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo los 100 mejores machos de 500
  newgeneaEE$V3<-sample(hembrasEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo las 100 mejores hembras de 500
  
  
  #escoger los reproductores RR
  machosRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==1]) 
  machosRR$V2<-EBVtot$V6[(machosRR$V1)*2-1]+EBVtot$V6[(machosRR$V1)*2]
  machosRR<-machosRR$V1[order(machosRR$V2,decreasing = T)]
  
  hembrasRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==2])
  hembrasRR$V2<-EBVtot$V6[(hembrasRR$V1)*2-1]+EBVtot$V6[(hembrasRR$V1)*2]
  hembrasRR<-hembrasRR$V1[order(hembrasRR$V2,decreasing = T)]
  
  newgeneaRR$V2<-sample(machosRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo los 100 mejores machos de 500
  newgeneaRR$V3<-sample(hembrasRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo las 100 mejores hembras de 500
  
  #animales cruzados
  
  newgeneaER$V2<-sample(newgeneaEE$V1[newgeneaEE$V4==1],nrow(newgeneaER),replace = T)
  newgeneaER$V3<-sample(newgeneaRR$V1[newgeneaRR$V4==2],nrow(newgeneaER),replace = T)
  
  newgeneaRE$V3<-sample(newgeneaEE$V1[newgeneaEE$V4==2],nrow(newgeneaER),replace = T)
  newgeneaRE$V2<-sample(newgeneaRR$V1[newgeneaRR$V4==1],nrow(newgeneaER),replace = T)
  
  geneaEE<-newgeneaEE
  geneaRR<-newgeneaRR
  geneaER<-newgeneaER
  geneaRE<-newgeneaRE
  
  #genea completa
  
  geneatot<-rbind(geneatot,newgeneaEE,newgeneaRR,newgeneaER,newgeneaRE)
  
  TBVEnew<-as.data.frame(matrix(0,nrow = nrow(geneaEE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVEnew)) %% 2
  TBVEnew[filas==1,]<-TBVtot[(geneaEE$V2*2)-1,]*0.5+TBVtot[geneaEE$V2*2,]*0.5
  TBVEnew[filas==0,]<-TBVtot[(geneaEE$V3*2)-1,]*0.5+TBVtot[geneaEE$V3*2,]*0.5
  TBVEnew<-TBVEnew+(mvrnorm(nrow(geneaEE)*2,medias,varcovartot/2))
  
  TBVRnew<-as.data.frame(matrix(0,nrow = nrow(geneaRR)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRnew)) %% 2
  TBVRnew[filas==1,]<-TBVtot[(geneaRR$V2*2)-1,]*0.5+TBVtot[geneaRR$V2*2,]*0.5
  TBVRnew[filas==0,]<-TBVtot[(geneaRR$V3*2)-1,]*0.5+TBVtot[geneaRR$V3*2,]*0.5
  TBVRnew<-TBVRnew+(mvrnorm(nrow(geneaRR)*2,medias,varcovartot/2))
  
  TBVE<-TBVEnew
  TBVR<-TBVRnew
  
  TBVtot<-rbind(TBVtot,TBVE,TBVR)
  
  TBVER<-as.data.frame(matrix(0,nrow = nrow(geneaER)*2,ncol = 8))
  filas<-seq_len(nrow(TBVER)) %% 2
  TBVER[filas==1,]<-TBVtot[(geneaER$V2*2)-1,]*0.5+TBVtot[geneaER$V2*2,]*0.5
  TBVER[filas==0,]<-TBVtot[(geneaER$V3*2)-1,]*0.5+TBVtot[geneaER$V3*2,]*0.5
  TBVER<-TBVER+(mvrnorm(nrow(geneaER)*2,medias,varcovartot/2))
  
  TBVRE<-as.data.frame(matrix(0,nrow = nrow(geneaRE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRE)) %% 2
  TBVRE[filas==1,]<-TBVtot[(geneaRE$V2*2)-1,]*0.5+TBVtot[geneaRE$V2*2,]*0.5
  TBVRE[filas==0,]<-TBVtot[(geneaRE$V3*2)-1,]*0.5+TBVtot[geneaRE$V3*2,]*0.5
  TBVRE<-TBVRE+(mvrnorm(nrow(geneaRE)*2,medias,varcovartot/2))
  
  TBVtot<-rbind(TBVtot,TBVER,TBVRE)
  
  #datos
  geneaEEtot<-geneatot[geneatot$Raza=="EE",]
  geneaRRtot<-geneatot[geneatot$Raza=="RR",]
  geneaERtot<-geneatot[geneatot$Raza=="ER",]
  geneaREtot<-geneatot[geneatot$Raza=="RE",]
  
  datosEEmachos<-data.frame(V1=(geneaEEtot$V1[geneaEEtot$V4==1]*2)-1,V2=geneaEEtot$V1[geneaEEtot$V4==1]*2,V3=1,V4=10,V5=0,V6=0,V7=0)
  datosEEmachos$V4<-TBVtot$V1[datosEEmachos$V1]+TBVtot$V2[datosEEmachos$V2]
  
  datosEEhembras<-data.frame(V1=(geneaEEtot$V1[geneaEEtot$V4==2]*2)-1,V2=geneaEEtot$V1[geneaEEtot$V4==2]*2,V3=1,V4=10,V5=0,V6=0,V7=0)
  datosEEhembras$V4<-TBVtot$V1[datosEEhembras$V1]+TBVtot$V2[datosEEhembras$V2]
  
  datosRRmachos<-data.frame(V1=(geneaRRtot$V1[geneaRRtot$V4==1]*2)-1,V2=geneaRRtot$V1[geneaRRtot$V4==1]*2,V3=1,V4=0,V5=10,V6=0,V7=0)
  datosRRmachos$V5<-TBVtot$V5[datosRRmachos$V1]+TBVtot$V6[datosRRmachos$V2]
  
  datosRRhembras<-data.frame(V1=(geneaRRtot$V1[geneaRRtot$V4==2]*2)-1,V2=geneaRRtot$V1[geneaRRtot$V4==2]*2,V3=1,V4=0,V5=10,V6=0,V7=0)
  datosRRhembras$V5<-TBVtot$V5[datosRRhembras$V1]+TBVtot$V6[datosRRhembras$V2]
  
  datosER<-data.frame(V1=(geneaERtot$V1*2)-1,V2=geneaERtot$V1*2,V3=1,V4=0,V5=0,V6=10,V7=0) #solo hay hembras
  datosER$V6<-TBVtot$V3[datosER$V1]+TBVtot$V8[datosER$V2]
  
  datosRE<-data.frame(V1=(geneaREtot$V1*2)-1,V2=geneaREtot$V1*2,V3=1,V4=0,V5=0,V6=0,V7=10) #solo hay hembras
  datosRE$V7<-TBVtot$V7[datosRE$V1]+TBVtot$V4[datosRE$V2]
  
  datostot<-rbind(datosEEmachos,datosEEhembras,datosRRmachos,datosRRhembras,datosER,datosRE)
  
}

save.image("selEBVmadrespuras.R")
load("selEBVmadrespuras.R")

geneaEEtot<-geneatot[geneatot$Raza=="EE",]
geneaRRtot<-geneatot[geneatot$Raza=="RR",]
geneaERtot<-geneatot[geneatot$Raza=="ER",]
geneaREtot<-geneatot[geneatot$Raza=="RE",]

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaEEprint$TBVmatER<-TBVtot$V4[(geneaEEprint$V1*2)-1]+TBVtot$V4[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaRRprint$TBVmatER<-TBVtot$V8[(geneaRRprint$V1*2)-1]+TBVtot$V8[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatER",colBy = "Gen")
plot(AlphaPart,ylim = c(-1,10),color = c(2,3))

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaEEprint$TBVmatRE<-TBVtot$V3[(geneaEEprint$V1*2)-1]+TBVtot$V3[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaRRprint$TBVmatRE<-TBVtot$V7[(geneaRRprint$V1*2)-1]+TBVtot$V7[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatRE",colBy = "Gen")
plot(AlphaPart,ylim = c(-1,10),color = c(2,3))


#seleccionar por madres cruzadas#

load("pobbase.RData")

for (Generation in 21:30) {
  BLUP(geneatot,datostot)
  solutions<-read.table('solutions',skip=1)
  
  EBVtot<-data.frame(V1=solutions$V4[solutions$V1==1 & solutions$V2==2],
                    V2=solutions$V4[solutions$V1==1 & solutions$V2==3],
                    V3=solutions$V4[solutions$V1==3 & solutions$V2==2],
                    V4=solutions$V4[solutions$V1==4 & solutions$V2==3],
                    V5=solutions$V4[solutions$V1==2 & solutions$V2==2],
                    V6=solutions$V4[solutions$V1==2 & solutions$V2==3],
                    V7=solutions$V4[solutions$V1==4 & solutions$V2==2],
                    V8=solutions$V4[solutions$V1==3 & solutions$V2==3])
  
  newgeneaEE<-data.frame(V1=c(1:1000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="EE")
  newgeneaRR<-data.frame(V1=c(1001:2000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="RR")
  newgeneaER<-data.frame(V1=c(2001:2500)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="ER")
  newgeneaRE<-data.frame(V1=c(2501:3000)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="RE")
  
  #escoger los reproductores EE
  machosEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==1]) 
  machosEE$V2<-EBVtot$V4[(machosEE$V1)*2-1]+EBVtot$V4[(machosEE$V1)*2]
  machosEE<-machosEE$V1[order(machosEE$V2,decreasing = T)]
  
  hembrasEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==2])
  hembrasEE$V2<-EBVtot$V4[(hembrasEE$V1)*2-1]+EBVtot$V4[(hembrasEE$V1)*2]
  hembrasEE<-hembrasEE$V1[order(hembrasEE$V2,decreasing = T)]
  
  newgeneaEE$V2<-sample(machosEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo los 100 mejores machos de 500
  newgeneaEE$V3<-sample(hembrasEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo las 100 mejores hembras de 500
  
  
  #escoger los reproductores RR
  machosRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==1]) 
  machosRR$V2<-EBVtot$V8[(machosRR$V1)*2-1]+EBVtot$V8[(machosRR$V1)*2]
  machosRR<-machosRR$V1[order(machosRR$V2,decreasing = T)]
  
  hembrasRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==2])
  hembrasRR$V2<-EBVtot$V8[(hembrasRR$V1)*2-1]+EBVtot$V8[(hembrasRR$V1)*2]
  hembrasRR<-hembrasRR$V1[order(hembrasRR$V2,decreasing = T)]
  
  newgeneaRR$V2<-sample(machosRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo los 100 mejores machos de 500
  newgeneaRR$V3<-sample(hembrasRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo las 100 mejores hembras de 500
  
  #animales cruzados
  
  newgeneaER$V2<-sample(newgeneaEE$V1[newgeneaEE$V4==1],nrow(newgeneaER),replace = T)
  newgeneaER$V3<-sample(newgeneaRR$V1[newgeneaRR$V4==2],nrow(newgeneaER),replace = T)
  
  newgeneaRE$V3<-sample(newgeneaEE$V1[newgeneaEE$V4==2],nrow(newgeneaER),replace = T)
  newgeneaRE$V2<-sample(newgeneaRR$V1[newgeneaRR$V4==1],nrow(newgeneaER),replace = T)
  
  geneaEE<-newgeneaEE
  geneaRR<-newgeneaRR
  geneaER<-newgeneaER
  geneaRE<-newgeneaRE
  
  #genea completa
  
  geneatot<-rbind(geneatot,newgeneaEE,newgeneaRR,newgeneaER,newgeneaRE)
  
  TBVEnew<-as.data.frame(matrix(0,nrow = nrow(geneaEE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVEnew)) %% 2
  TBVEnew[filas==1,]<-TBVtot[(geneaEE$V2*2)-1,]*0.5+TBVtot[geneaEE$V2*2,]*0.5
  TBVEnew[filas==0,]<-TBVtot[(geneaEE$V3*2)-1,]*0.5+TBVtot[geneaEE$V3*2,]*0.5
  TBVEnew<-TBVEnew+(mvrnorm(nrow(geneaEE)*2,medias,varcovartot/2))
  
  TBVRnew<-as.data.frame(matrix(0,nrow = nrow(geneaRR)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRnew)) %% 2
  TBVRnew[filas==1,]<-TBVtot[(geneaRR$V2*2)-1,]*0.5+TBVtot[geneaRR$V2*2,]*0.5
  TBVRnew[filas==0,]<-TBVtot[(geneaRR$V3*2)-1,]*0.5+TBVtot[geneaRR$V3*2,]*0.5
  TBVRnew<-TBVRnew+(mvrnorm(nrow(geneaRR)*2,medias,varcovartot/2))
  
  TBVE<-TBVEnew
  TBVR<-TBVRnew
  
  TBVtot<-rbind(TBVtot,TBVE,TBVR)
  
  TBVER<-as.data.frame(matrix(0,nrow = nrow(geneaER)*2,ncol = 8))
  filas<-seq_len(nrow(TBVER)) %% 2
  TBVER[filas==1,]<-TBVtot[(geneaER$V2*2)-1,]*0.5+TBVtot[geneaER$V2*2,]*0.5
  TBVER[filas==0,]<-TBVtot[(geneaER$V3*2)-1,]*0.5+TBVtot[geneaER$V3*2,]*0.5
  TBVER<-TBVER+(mvrnorm(nrow(geneaER)*2,medias,varcovartot/2))
  
  TBVRE<-as.data.frame(matrix(0,nrow = nrow(geneaRE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRE)) %% 2
  TBVRE[filas==1,]<-TBVtot[(geneaRE$V2*2)-1,]*0.5+TBVtot[geneaRE$V2*2,]*0.5
  TBVRE[filas==0,]<-TBVtot[(geneaRE$V3*2)-1,]*0.5+TBVtot[geneaRE$V3*2,]*0.5
  TBVRE<-TBVRE+(mvrnorm(nrow(geneaRE)*2,medias,varcovartot/2))
  
  TBVtot<-rbind(TBVtot,TBVER,TBVRE)
  
  #datos
  geneaEEtot<-geneatot[geneatot$Raza=="EE",]
  geneaRRtot<-geneatot[geneatot$Raza=="RR",]
  geneaERtot<-geneatot[geneatot$Raza=="ER",]
  geneaREtot<-geneatot[geneatot$Raza=="RE",]
  
  datosEEmachos<-data.frame(V1=(geneaEEtot$V1[geneaEEtot$V4==1]*2)-1,V2=geneaEEtot$V1[geneaEEtot$V4==1]*2,V3=1,V4=10,V5=0,V6=0,V7=0)
  datosEEmachos$V4<-TBVtot$V1[datosEEmachos$V1]+TBVtot$V2[datosEEmachos$V2]
  
  datosEEhembras<-data.frame(V1=(geneaEEtot$V1[geneaEEtot$V4==2]*2)-1,V2=geneaEEtot$V1[geneaEEtot$V4==2]*2,V3=1,V4=10,V5=0,V6=0,V7=0)
  datosEEhembras$V4<-TBVtot$V1[datosEEhembras$V1]+TBVtot$V2[datosEEhembras$V2]
  
  datosRRmachos<-data.frame(V1=(geneaRRtot$V1[geneaRRtot$V4==1]*2)-1,V2=geneaRRtot$V1[geneaRRtot$V4==1]*2,V3=1,V4=0,V5=10,V6=0,V7=0)
  datosRRmachos$V5<-TBVtot$V5[datosRRmachos$V1]+TBVtot$V6[datosRRmachos$V2]
  
  datosRRhembras<-data.frame(V1=(geneaRRtot$V1[geneaRRtot$V4==2]*2)-1,V2=geneaRRtot$V1[geneaRRtot$V4==2]*2,V3=1,V4=0,V5=10,V6=0,V7=0)
  datosRRhembras$V5<-TBVtot$V5[datosRRhembras$V1]+TBVtot$V6[datosRRhembras$V2]
  
  datosER<-data.frame(V1=(geneaERtot$V1*2)-1,V2=geneaERtot$V1*2,V3=1,V4=0,V5=0,V6=10,V7=0) #solo hay hembras
  datosER$V6<-TBVtot$V3[datosER$V1]+TBVtot$V8[datosER$V2]
  
  datosRE<-data.frame(V1=(geneaREtot$V1*2)-1,V2=geneaREtot$V1*2,V3=1,V4=0,V5=0,V6=0,V7=10) #solo hay hembras
  datosRE$V7<-TBVtot$V7[datosRE$V1]+TBVtot$V4[datosRE$V2]
  
  datostot<-rbind(datosEEmachos,datosEEhembras,datosRRmachos,datosRRhembras,datosER,datosRE)
  
}

save.image("selEBVmadrescruce.R")
load("selEBVmadrescruce.R")

geneaEEtot<-geneatot[geneatot$Raza=="EE",]
geneaRRtot<-geneatot[geneatot$Raza=="RR",]
geneaERtot<-geneatot[geneatot$Raza=="ER",]
geneaREtot<-geneatot[geneatot$Raza=="RE",]

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaEEprint$TBVmatER<-TBVtot$V4[(geneaEEprint$V1*2)-1]+TBVtot$V4[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaRRprint$TBVmatER<-TBVtot$V8[(geneaRRprint$V1*2)-1]+TBVtot$V8[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatER",colBy = "Gen")
plot(AlphaPart,ylim = c(-1,10),color = c(2,3))

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(c(1:41),each=1000))
geneaEEprint$TBVmatRE<-TBVtot$V3[(geneaEEprint$V1*2)-1]+TBVtot$V3[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(c(1:41),each=1000))
geneaRRprint$TBVmatRE<-TBVtot$V7[(geneaRRprint$V1*2)-1]+TBVtot$V7[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatRE",colBy = "Gen")
plot(AlphaPart,ylim = c(-1,10),color = c(2,3))

#seleccionar por madres cruzadas no info RE#

load("pobbase.RData")

datosRE<-data.frame(V1=(geneaREtot$V1*2)-1,V2=geneaREtot$V1*2,V3=1,V4=0,V5=0,V6=0,V7=0) #solo hay hembras
#datosRE$V7<-TBVtot$V8[datosRE$V1]+TBVtot$V4[datosRE$V2]

datostot<-rbind(datosEEmachos,datosEEhembras,datosRRmachos,datosRRhembras,datosER,datosRE)

for (Generation in 21:30) {
  BLUP(geneatot,datostot)
  solutions<-read.table('solutions',skip=1)
  
  EBVtot<-data.frame(V1=solutions$V4[solutions$V1==1 & solutions$V2==2],
                     V2=solutions$V4[solutions$V1==1 & solutions$V2==3],
                     V3=solutions$V4[solutions$V1==3 & solutions$V2==2],
                     V4=solutions$V4[solutions$V1==4 & solutions$V2==3],
                     V5=solutions$V4[solutions$V1==2 & solutions$V2==2],
                     V6=solutions$V4[solutions$V1==2 & solutions$V2==3],
                     V7=solutions$V4[solutions$V1==4 & solutions$V2==2],
                     V8=solutions$V4[solutions$V1==3 & solutions$V2==3])
  
  newgeneaEE<-data.frame(V1=c(1:1000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="EE")
  newgeneaRR<-data.frame(V1=c(1001:2000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="RR")
  newgeneaER<-data.frame(V1=c(2001:2500)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="ER")
  newgeneaRE<-data.frame(V1=c(2501:3000)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="RE")
  
  #escoger los reproductores EE
  machosEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==1]) 
  machosEE$V2<-EBVtot$V4[(machosEE$V1)*2-1]+EBVtot$V4[(machosEE$V1)*2]
  machosEE<-machosEE$V1[order(machosEE$V2,decreasing = T)]
  
  hembrasEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==2])
  hembrasEE$V2<-EBVtot$V4[(hembrasEE$V1)*2-1]+EBVtot$V4[(hembrasEE$V1)*2]
  hembrasEE<-hembrasEE$V1[order(hembrasEE$V2,decreasing = T)]
  
  newgeneaEE$V2<-sample(machosEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo los 100 mejores machos de 500
  newgeneaEE$V3<-sample(hembrasEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo las 100 mejores hembras de 500
  
  
  #escoger los reproductores RR
  machosRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==1]) 
  machosRR$V2<-EBVtot$V8[(machosRR$V1)*2-1]+EBVtot$V8[(machosRR$V1)*2]
  machosRR<-machosRR$V1[order(machosRR$V2,decreasing = T)]
  
  hembrasRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==2])
  hembrasRR$V2<-EBVtot$V8[(hembrasRR$V1)*2-1]+EBVtot$V8[(hembrasRR$V1)*2]
  hembrasRR<-hembrasRR$V1[order(hembrasRR$V2,decreasing = T)]
  
  newgeneaRR$V2<-sample(machosRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo los 100 mejores machos de 500
  newgeneaRR$V3<-sample(hembrasRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo las 100 mejores hembras de 500
  
  #animales cruzados
  
  newgeneaER$V2<-sample(newgeneaEE$V1[newgeneaEE$V4==1],nrow(newgeneaER),replace = T)
  newgeneaER$V3<-sample(newgeneaRR$V1[newgeneaRR$V4==2],nrow(newgeneaER),replace = T)
  
  newgeneaRE$V3<-sample(newgeneaEE$V1[newgeneaEE$V4==2],nrow(newgeneaER),replace = T)
  newgeneaRE$V2<-sample(newgeneaRR$V1[newgeneaRR$V4==1],nrow(newgeneaER),replace = T)
  
  geneaEE<-newgeneaEE
  geneaRR<-newgeneaRR
  geneaER<-newgeneaER
  geneaRE<-newgeneaRE
  
  #genea completa
  
  geneatot<-rbind(geneatot,newgeneaEE,newgeneaRR,newgeneaER,newgeneaRE)
  
  TBVEnew<-as.data.frame(matrix(0,nrow = nrow(geneaEE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVEnew)) %% 2
  TBVEnew[filas==1,]<-TBVtot[(geneaEE$V2*2)-1,]*0.5+TBVtot[geneaEE$V2*2,]*0.5
  TBVEnew[filas==0,]<-TBVtot[(geneaEE$V3*2)-1,]*0.5+TBVtot[geneaEE$V3*2,]*0.5
  TBVEnew<-TBVEnew+(mvrnorm(nrow(geneaEE)*2,medias,varcovartot/2))
  
  TBVRnew<-as.data.frame(matrix(0,nrow = nrow(geneaRR)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRnew)) %% 2
  TBVRnew[filas==1,]<-TBVtot[(geneaRR$V2*2)-1,]*0.5+TBVtot[geneaRR$V2*2,]*0.5
  TBVRnew[filas==0,]<-TBVtot[(geneaRR$V3*2)-1,]*0.5+TBVtot[geneaRR$V3*2,]*0.5
  TBVRnew<-TBVRnew+(mvrnorm(nrow(geneaRR)*2,medias,varcovartot/2))
  
  TBVE<-TBVEnew
  TBVR<-TBVRnew
  
  TBVtot<-rbind(TBVtot,TBVE,TBVR)
  
  TBVER<-as.data.frame(matrix(0,nrow = nrow(geneaER)*2,ncol = 8))
  filas<-seq_len(nrow(TBVER)) %% 2
  TBVER[filas==1,]<-TBVtot[(geneaER$V2*2)-1,]*0.5+TBVtot[geneaER$V2*2,]*0.5
  TBVER[filas==0,]<-TBVtot[(geneaER$V3*2)-1,]*0.5+TBVtot[geneaER$V3*2,]*0.5
  TBVER<-TBVER+(mvrnorm(nrow(geneaER)*2,medias,varcovartot/2))
  
  TBVRE<-as.data.frame(matrix(0,nrow = nrow(geneaRE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRE)) %% 2
  TBVRE[filas==1,]<-TBVtot[(geneaRE$V2*2)-1,]*0.5+TBVtot[geneaRE$V2*2,]*0.5
  TBVRE[filas==0,]<-TBVtot[(geneaRE$V3*2)-1,]*0.5+TBVtot[geneaRE$V3*2,]*0.5
  TBVRE<-TBVRE+(mvrnorm(nrow(geneaRE)*2,medias,varcovartot/2))
  
  TBVtot<-rbind(TBVtot,TBVER,TBVRE)
  
  #datos
  geneaEEtot<-geneatot[geneatot$Raza=="EE",]
  geneaRRtot<-geneatot[geneatot$Raza=="RR",]
  geneaERtot<-geneatot[geneatot$Raza=="ER",]
  geneaREtot<-geneatot[geneatot$Raza=="RE",]
  
  datosEEmachos<-data.frame(V1=(geneaEEtot$V1[geneaEEtot$V4==1]*2)-1,V2=geneaEEtot$V1[geneaEEtot$V4==1]*2,V3=1,V4=10,V5=0,V6=0,V7=0)
  datosEEmachos$V4<-TBVtot$V1[datosEEmachos$V1]+TBVtot$V2[datosEEmachos$V2]
  
  datosEEhembras<-data.frame(V1=(geneaEEtot$V1[geneaEEtot$V4==2]*2)-1,V2=geneaEEtot$V1[geneaEEtot$V4==2]*2,V3=1,V4=10,V5=0,V6=0,V7=0)
  datosEEhembras$V4<-TBVtot$V1[datosEEhembras$V1]+TBVtot$V2[datosEEhembras$V2]
  
  datosRRmachos<-data.frame(V1=(geneaRRtot$V1[geneaRRtot$V4==1]*2)-1,V2=geneaRRtot$V1[geneaRRtot$V4==1]*2,V3=1,V4=0,V5=10,V6=0,V7=0)
  datosRRmachos$V5<-TBVtot$V5[datosRRmachos$V1]+TBVtot$V6[datosRRmachos$V2]
  
  datosRRhembras<-data.frame(V1=(geneaRRtot$V1[geneaRRtot$V4==2]*2)-1,V2=geneaRRtot$V1[geneaRRtot$V4==2]*2,V3=1,V4=0,V5=10,V6=0,V7=0)
  datosRRhembras$V5<-TBVtot$V5[datosRRhembras$V1]+TBVtot$V6[datosRRhembras$V2]
  
  datosER<-data.frame(V1=(geneaERtot$V1*2)-1,V2=geneaERtot$V1*2,V3=1,V4=0,V5=0,V6=10,V7=0) #solo hay hembras
  datosER$V6<-TBVtot$V3[datosER$V1]+TBVtot$V8[datosER$V2]
  
  datosRE<-data.frame(V1=(geneaREtot$V1*2)-1,V2=geneaREtot$V1*2,V3=1,V4=0,V5=0,V6=0,V7=10) #solo hay hembras
  #datosRE$V7<-TBVtot$V7[datosRE$V1]+TBVtot$V4[datosRE$V2]
  
  datostot<-rbind(datosEEmachos,datosEEhembras,datosRRmachos,datosRRhembras,datosER,datosRE)
  
}

save.image("selEBVmadrescrucenoinfoRE.R")
load("selEBVmadrescrucenoinfoRE.R")

geneaEEtot<-geneatot[geneatot$Raza=="EE",]
geneaRRtot<-geneatot[geneatot$Raza=="RR",]
geneaERtot<-geneatot[geneatot$Raza=="ER",]
geneaREtot<-geneatot[geneatot$Raza=="RE",]

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaEEprint$TBVmatER<-TBVtot$V4[(geneaEEprint$V1*2)-1]+TBVtot$V4[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaRRprint$TBVmatER<-TBVtot$V8[(geneaRRprint$V1*2)-1]+TBVtot$V8[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatER",colBy = "Gen")
plot(AlphaPart,ylim = c(-1,10),color = c(2,3))

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(c(1:41),each=1000))
geneaEEprint$TBVmatRE<-TBVtot$V3[(geneaEEprint$V1*2)-1]+TBVtot$V3[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(c(1:41),each=1000))
geneaRRprint$TBVmatRE<-TBVtot$V7[(geneaRRprint$V1*2)-1]+TBVtot$V7[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatRE",colBy = "Gen")
plot(AlphaPart,ylim = c(-1,10),color = c(2,3))

#seleccionar por EE todo y RR madre#

load("pobbase.RData")

for (Generation in 21:30) {
  BLUP(geneatot,datostot)
  solutions<-read.table('solutions',skip=1)
  
  EBVtot<-data.frame(V1=solutions$V4[solutions$V1==1 & solutions$V2==2],
                     V2=solutions$V4[solutions$V1==1 & solutions$V2==3],
                     V3=solutions$V4[solutions$V1==3 & solutions$V2==2],
                     V4=solutions$V4[solutions$V1==4 & solutions$V2==3],
                     V5=solutions$V4[solutions$V1==2 & solutions$V2==2],
                     V6=solutions$V4[solutions$V1==2 & solutions$V2==3],
                     V7=solutions$V4[solutions$V1==4 & solutions$V2==2],
                     V8=solutions$V4[solutions$V1==3 & solutions$V2==3])
  
  newgeneaEE<-data.frame(V1=c(1:1000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="EE")
  newgeneaRR<-data.frame(V1=c(1001:2000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="RR")
  newgeneaER<-data.frame(V1=c(2001:2500)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="ER")
  newgeneaRE<-data.frame(V1=c(2501:3000)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="RE")
  
  #escoger los reproductores EE
  machosEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==1]) 
  machosEE$V2<-EBVtot$V1[(machosEE$V1)*2-1]+EBVtot$V1[(machosEE$V1)*2]+EBVtot$V2[(machosEE$V1)*2-1]+EBVtot$V2[(machosEE$V1)*2]
  machosEE<-machosEE$V1[order(machosEE$V2,decreasing = T)]
  
  hembrasEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==2])
  hembrasEE$V2<-EBVtot$V1[(hembrasEE$V1)*2-1]+EBVtot$V1[(hembrasEE$V1)*2]+EBVtot$V2[(hembrasEE$V1)*2-1]+EBVtot$V2[(hembrasEE$V1)*2]
  hembrasEE<-hembrasEE$V1[order(hembrasEE$V2,decreasing = T)]
  
  newgeneaEE$V2<-sample(machosEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo los 100 mejores machos de 500
  newgeneaEE$V3<-sample(hembrasEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo las 100 mejores hembras de 500
  
  
  #escoger los reproductores RR
  machosRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==1]) 
  machosRR$V2<-EBVtot$V6[(machosRR$V1)*2-1]+EBVtot$V6[(machosRR$V1)*2]
  machosRR<-machosRR$V1[order(machosRR$V2,decreasing = T)]
  
  hembrasRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==2])
  hembrasRR$V2<-EBVtot$V6[(hembrasRR$V1)*2-1]+EBVtot$V6[(hembrasRR$V1)*2]
  hembrasRR<-hembrasRR$V1[order(hembrasRR$V2,decreasing = T)]
  
  newgeneaRR$V2<-sample(machosRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo los 100 mejores machos de 500
  newgeneaRR$V3<-sample(hembrasRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo las 100 mejores hembras de 500
  
  #animales cruzados
  
  newgeneaER$V2<-sample(newgeneaEE$V1[newgeneaEE$V4==1],nrow(newgeneaER),replace = T)
  newgeneaER$V3<-sample(newgeneaRR$V1[newgeneaRR$V4==2],nrow(newgeneaER),replace = T)
  
  newgeneaRE$V3<-sample(newgeneaEE$V1[newgeneaEE$V4==2],nrow(newgeneaER),replace = T)
  newgeneaRE$V2<-sample(newgeneaRR$V1[newgeneaRR$V4==1],nrow(newgeneaER),replace = T)
  
  geneaEE<-newgeneaEE
  geneaRR<-newgeneaRR
  geneaER<-newgeneaER
  geneaRE<-newgeneaRE
  
  #genea completa
  
  geneatot<-rbind(geneatot,newgeneaEE,newgeneaRR,newgeneaER,newgeneaRE)
  
  TBVEnew<-as.data.frame(matrix(0,nrow = nrow(geneaEE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVEnew)) %% 2
  TBVEnew[filas==1,]<-TBVtot[(geneaEE$V2*2)-1,]*0.5+TBVtot[geneaEE$V2*2,]*0.5
  TBVEnew[filas==0,]<-TBVtot[(geneaEE$V3*2)-1,]*0.5+TBVtot[geneaEE$V3*2,]*0.5
  TBVEnew<-TBVEnew+(mvrnorm(nrow(geneaEE)*2,medias,varcovartot/2))
  
  TBVRnew<-as.data.frame(matrix(0,nrow = nrow(geneaRR)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRnew)) %% 2
  TBVRnew[filas==1,]<-TBVtot[(geneaRR$V2*2)-1,]*0.5+TBVtot[geneaRR$V2*2,]*0.5
  TBVRnew[filas==0,]<-TBVtot[(geneaRR$V3*2)-1,]*0.5+TBVtot[geneaRR$V3*2,]*0.5
  TBVRnew<-TBVRnew+(mvrnorm(nrow(geneaRR)*2,medias,varcovartot/2))
  
  TBVE<-TBVEnew
  TBVR<-TBVRnew
  
  TBVtot<-rbind(TBVtot,TBVE,TBVR)
  
  TBVER<-as.data.frame(matrix(0,nrow = nrow(geneaER)*2,ncol = 8))
  filas<-seq_len(nrow(TBVER)) %% 2
  TBVER[filas==1,]<-TBVtot[(geneaER$V2*2)-1,]*0.5+TBVtot[geneaER$V2*2,]*0.5
  TBVER[filas==0,]<-TBVtot[(geneaER$V3*2)-1,]*0.5+TBVtot[geneaER$V3*2,]*0.5
  TBVER<-TBVER+(mvrnorm(nrow(geneaER)*2,medias,varcovartot/2))
  
  TBVRE<-as.data.frame(matrix(0,nrow = nrow(geneaRE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRE)) %% 2
  TBVRE[filas==1,]<-TBVtot[(geneaRE$V2*2)-1,]*0.5+TBVtot[geneaRE$V2*2,]*0.5
  TBVRE[filas==0,]<-TBVtot[(geneaRE$V3*2)-1,]*0.5+TBVtot[geneaRE$V3*2,]*0.5
  TBVRE<-TBVRE+(mvrnorm(nrow(geneaRE)*2,medias,varcovartot/2))
  
  TBVtot<-rbind(TBVtot,TBVER,TBVRE)
  
  #datos
  geneaEEtot<-geneatot[geneatot$Raza=="EE",]
  geneaRRtot<-geneatot[geneatot$Raza=="RR",]
  geneaERtot<-geneatot[geneatot$Raza=="ER",]
  geneaREtot<-geneatot[geneatot$Raza=="RE",]
  
  datosEEmachos<-data.frame(V1=(geneaEEtot$V1[geneaEEtot$V4==1]*2)-1,V2=geneaEEtot$V1[geneaEEtot$V4==1]*2,V3=1,V4=10,V5=0,V6=0,V7=0)
  datosEEmachos$V4<-TBVtot$V1[datosEEmachos$V1]+TBVtot$V2[datosEEmachos$V2]
  
  datosEEhembras<-data.frame(V1=(geneaEEtot$V1[geneaEEtot$V4==2]*2)-1,V2=geneaEEtot$V1[geneaEEtot$V4==2]*2,V3=1,V4=10,V5=0,V6=0,V7=0)
  datosEEhembras$V4<-TBVtot$V1[datosEEhembras$V1]+TBVtot$V2[datosEEhembras$V2]
  
  datosRRmachos<-data.frame(V1=(geneaRRtot$V1[geneaRRtot$V4==1]*2)-1,V2=geneaRRtot$V1[geneaRRtot$V4==1]*2,V3=1,V4=0,V5=10,V6=0,V7=0)
  datosRRmachos$V5<-TBVtot$V5[datosRRmachos$V1]+TBVtot$V6[datosRRmachos$V2]
  
  datosRRhembras<-data.frame(V1=(geneaRRtot$V1[geneaRRtot$V4==2]*2)-1,V2=geneaRRtot$V1[geneaRRtot$V4==2]*2,V3=1,V4=0,V5=10,V6=0,V7=0)
  datosRRhembras$V5<-TBVtot$V5[datosRRhembras$V1]+TBVtot$V6[datosRRhembras$V2]
  
  datosER<-data.frame(V1=(geneaERtot$V1*2)-1,V2=geneaERtot$V1*2,V3=1,V4=0,V5=0,V6=10,V7=0) #solo hay hembras
  datosER$V6<-TBVtot$V3[datosER$V1]+TBVtot$V8[datosER$V2]
  
  datosRE<-data.frame(V1=(geneaREtot$V1*2)-1,V2=geneaREtot$V1*2,V3=1,V4=0,V5=0,V6=0,V7=10) #solo hay hembras
  datosRE$V7<-TBVtot$V7[datosRE$V1]+TBVtot$V4[datosRE$V2]
  
  datostot<-rbind(datosEEmachos,datosEEhembras,datosRRmachos,datosRRhembras,datosER,datosRE)
  
}

save.image("selEBVtodoEEmadresRR.R")
load("selEBVtodoEEmadresRR.R")

geneaEEtot<-geneatot[geneatot$Raza=="EE",]
geneaRRtot<-geneatot[geneatot$Raza=="RR",]
geneaERtot<-geneatot[geneatot$Raza=="ER",]
geneaREtot<-geneatot[geneatot$Raza=="RE",]

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaEEprint$TBVmatER<-TBVtot$V4[(geneaEEprint$V1*2)-1]+TBVtot$V4[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaRRprint$TBVmatER<-TBVtot$V8[(geneaRRprint$V1*2)-1]+TBVtot$V8[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatER",colBy = "Gen")
plot(AlphaPart,ylim = c(-1,10),color = c(2,3))

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaEEprint$TBVmatRE<-TBVtot$V3[(geneaEEprint$V1*2)-1]+TBVtot$V3[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaRRprint$TBVmatRE<-TBVtot$V7[(geneaRRprint$V1*2)-1]+TBVtot$V7[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatRE",colBy = "Gen")
plot(AlphaPart,ylim = c(-1,10),color = c(2,3))

#seleccionar por madres cruzadas sin datos cruce#

load("pobbase.RData")

datosER<-data.frame(V1=(geneaERtot$V1*2)-1,V2=geneaERtot$V1*2,V3=1,V4=0,V5=0,V6=0,V7=0) #solo hay hembras
#datosER$V6<-TBVtot$V4[datosER$V1]+TBVtot$V8[datosER$V2]

datosRE<-data.frame(V1=(geneaREtot$V1*2)-1,V2=geneaREtot$V1*2,V3=1,V4=0,V5=0,V6=0,V7=0) #solo hay hembras
#datosRE$V7<-TBVtot$V8[datosRE$V1]+TBVtot$V4[datosRE$V2]

datostot<-rbind(datosEEmachos,datosEEhembras,datosRRmachos,datosRRhembras,datosER,datosRE)

for (Generation in 21:30) {
  BLUP(geneatot,datostot)
  solutions<-read.table('solutions',skip=1)
  
  EBVtot<-data.frame(V1=solutions$V4[solutions$V1==1 & solutions$V2==2],
                     V2=solutions$V4[solutions$V1==1 & solutions$V2==3],
                     V3=solutions$V4[solutions$V1==3 & solutions$V2==2],
                     V4=solutions$V4[solutions$V1==4 & solutions$V2==3],
                     V5=solutions$V4[solutions$V1==2 & solutions$V2==2],
                     V6=solutions$V4[solutions$V1==2 & solutions$V2==3],
                     V7=solutions$V4[solutions$V1==4 & solutions$V2==2],
                     V8=solutions$V4[solutions$V1==3 & solutions$V2==3])
  
  newgeneaEE<-data.frame(V1=c(1:1000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="EE")
  newgeneaRR<-data.frame(V1=c(1001:2000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="RR")
  newgeneaER<-data.frame(V1=c(2001:2500)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="ER")
  newgeneaRE<-data.frame(V1=c(2501:3000)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="RE")
  
  #escoger los reproductores EE
  machosEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==1]) 
  machosEE$V2<-EBVtot$V4[(machosEE$V1)*2-1]+EBVtot$V4[(machosEE$V1)*2]
  machosEE<-machosEE$V1[order(machosEE$V2,decreasing = T)]
  
  hembrasEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==2])
  hembrasEE$V2<-EBVtot$V4[(hembrasEE$V1)*2-1]+EBVtot$V4[(hembrasEE$V1)*2]
  hembrasEE<-hembrasEE$V1[order(hembrasEE$V2,decreasing = T)]
  
  newgeneaEE$V2<-sample(machosEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo los 100 mejores machos de 500
  newgeneaEE$V3<-sample(hembrasEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo las 100 mejores hembras de 500
  
  
  #escoger los reproductores RR
  machosRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==1]) 
  machosRR$V2<-EBVtot$V8[(machosRR$V1)*2-1]+EBVtot$V8[(machosRR$V1)*2]
  machosRR<-machosRR$V1[order(machosRR$V2,decreasing = T)]
  
  hembrasRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==2])
  hembrasRR$V2<-EBVtot$V8[(hembrasRR$V1)*2-1]+EBVtot$V8[(hembrasRR$V1)*2]
  hembrasRR<-hembrasRR$V1[order(hembrasRR$V2,decreasing = T)]
  
  newgeneaRR$V2<-sample(machosRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo los 100 mejores machos de 500
  newgeneaRR$V3<-sample(hembrasRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo las 100 mejores hembras de 500
  
  #animales cruzados
  
  newgeneaER$V2<-sample(newgeneaEE$V1[newgeneaEE$V4==1],nrow(newgeneaER),replace = T)
  newgeneaER$V3<-sample(newgeneaRR$V1[newgeneaRR$V4==2],nrow(newgeneaER),replace = T)
  
  newgeneaRE$V3<-sample(newgeneaEE$V1[newgeneaEE$V4==2],nrow(newgeneaER),replace = T)
  newgeneaRE$V2<-sample(newgeneaRR$V1[newgeneaRR$V4==1],nrow(newgeneaER),replace = T)
  
  geneaEE<-newgeneaEE
  geneaRR<-newgeneaRR
  geneaER<-newgeneaER
  geneaRE<-newgeneaRE
  
  #genea completa
  
  geneatot<-rbind(geneatot,newgeneaEE,newgeneaRR,newgeneaER,newgeneaRE)
  
  TBVEnew<-as.data.frame(matrix(0,nrow = nrow(geneaEE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVEnew)) %% 2
  TBVEnew[filas==1,]<-TBVtot[(geneaEE$V2*2)-1,]*0.5+TBVtot[geneaEE$V2*2,]*0.5
  TBVEnew[filas==0,]<-TBVtot[(geneaEE$V3*2)-1,]*0.5+TBVtot[geneaEE$V3*2,]*0.5
  TBVEnew<-TBVEnew+(mvrnorm(nrow(geneaEE)*2,medias,varcovartot/2))
  
  TBVRnew<-as.data.frame(matrix(0,nrow = nrow(geneaRR)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRnew)) %% 2
  TBVRnew[filas==1,]<-TBVtot[(geneaRR$V2*2)-1,]*0.5+TBVtot[geneaRR$V2*2,]*0.5
  TBVRnew[filas==0,]<-TBVtot[(geneaRR$V3*2)-1,]*0.5+TBVtot[geneaRR$V3*2,]*0.5
  TBVRnew<-TBVRnew+(mvrnorm(nrow(geneaRR)*2,medias,varcovartot/2))
  
  TBVE<-TBVEnew
  TBVR<-TBVRnew
  
  TBVtot<-rbind(TBVtot,TBVE,TBVR)
  
  TBVER<-as.data.frame(matrix(0,nrow = nrow(geneaER)*2,ncol = 8))
  filas<-seq_len(nrow(TBVER)) %% 2
  TBVER[filas==1,]<-TBVtot[(geneaER$V2*2)-1,]*0.5+TBVtot[geneaER$V2*2,]*0.5
  TBVER[filas==0,]<-TBVtot[(geneaER$V3*2)-1,]*0.5+TBVtot[geneaER$V3*2,]*0.5
  TBVER<-TBVER+(mvrnorm(nrow(geneaER)*2,medias,varcovartot/2))
  
  TBVRE<-as.data.frame(matrix(0,nrow = nrow(geneaRE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRE)) %% 2
  TBVRE[filas==1,]<-TBVtot[(geneaRE$V2*2)-1,]*0.5+TBVtot[geneaRE$V2*2,]*0.5
  TBVRE[filas==0,]<-TBVtot[(geneaRE$V3*2)-1,]*0.5+TBVtot[geneaRE$V3*2,]*0.5
  TBVRE<-TBVRE+(mvrnorm(nrow(geneaRE)*2,medias,varcovartot/2))
  
  TBVtot<-rbind(TBVtot,TBVER,TBVRE)
  
  #datos
  geneaEEtot<-geneatot[geneatot$Raza=="EE",]
  geneaRRtot<-geneatot[geneatot$Raza=="RR",]
  geneaERtot<-geneatot[geneatot$Raza=="ER",]
  geneaREtot<-geneatot[geneatot$Raza=="RE",]
  
  datosEEmachos<-data.frame(V1=(geneaEEtot$V1[geneaEEtot$V4==1]*2)-1,V2=geneaEEtot$V1[geneaEEtot$V4==1]*2,V3=1,V4=10,V5=0,V6=0,V7=0)
  datosEEmachos$V4<-TBVtot$V1[datosEEmachos$V1]+TBVtot$V2[datosEEmachos$V2]
  
  datosEEhembras<-data.frame(V1=(geneaEEtot$V1[geneaEEtot$V4==2]*2)-1,V2=geneaEEtot$V1[geneaEEtot$V4==2]*2,V3=1,V4=10,V5=0,V6=0,V7=0)
  datosEEhembras$V4<-TBVtot$V1[datosEEhembras$V1]+TBVtot$V2[datosEEhembras$V2]
  
  datosRRmachos<-data.frame(V1=(geneaRRtot$V1[geneaRRtot$V4==1]*2)-1,V2=geneaRRtot$V1[geneaRRtot$V4==1]*2,V3=1,V4=0,V5=10,V6=0,V7=0)
  datosRRmachos$V5<-TBVtot$V5[datosRRmachos$V1]+TBVtot$V6[datosRRmachos$V2]
  
  datosRRhembras<-data.frame(V1=(geneaRRtot$V1[geneaRRtot$V4==2]*2)-1,V2=geneaRRtot$V1[geneaRRtot$V4==2]*2,V3=1,V4=0,V5=10,V6=0,V7=0)
  datosRRhembras$V5<-TBVtot$V5[datosRRhembras$V1]+TBVtot$V6[datosRRhembras$V2]
  
  datosER<-data.frame(V1=(geneaERtot$V1*2)-1,V2=geneaERtot$V1*2,V3=1,V4=0,V5=0,V6=10,V7=0) #solo hay hembras
  #datosER$V6<-TBVtot$V3[datosER$V1]+TBVtot$V8[datosER$V2]
  
  datosRE<-data.frame(V1=(geneaREtot$V1*2)-1,V2=geneaREtot$V1*2,V3=1,V4=0,V5=0,V6=0,V7=10) #solo hay hembras
  #datosRE$V7<-TBVtot$V7[datosRE$V1]+TBVtot$V4[datosRE$V2]
  
  datostot<-rbind(datosEEmachos,datosEEhembras,datosRRmachos,datosRRhembras,datosER,datosRE)
  
}

save.image("selEBVmadrescrucesindato.R")
load("selEBVmadrescrucesindato.R")

geneaEEtot<-geneatot[geneatot$Raza=="EE",]
geneaRRtot<-geneatot[geneatot$Raza=="RR",]
geneaERtot<-geneatot[geneatot$Raza=="ER",]
geneaREtot<-geneatot[geneatot$Raza=="RE",]

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaEEprint$TBVmatER<-TBVtot$V4[(geneaEEprint$V1*2)-1]+TBVtot$V4[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaRRprint$TBVmatER<-TBVtot$V8[(geneaRRprint$V1*2)-1]+TBVtot$V8[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatER",colBy = "Gen")
plot(AlphaPart,ylim = c(-1,10),color = c(2,3))

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaEEprint$TBVmatRE<-TBVtot$V3[(geneaEEprint$V1*2)-1]+TBVtot$V3[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaRRprint$TBVmatRE<-TBVtot$V7[(geneaRRprint$V1*2)-1]+TBVtot$V7[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatRE",colBy = "Gen")
plot(AlphaPart,ylim = c(-1,10),color = c(2,3))

#seleccionar media cruce#

load("pobbase.RData")

for (Generation in 21:30) {
  BLUP(geneatot,datostot)
  solutions<-read.table('solutions',skip=1)
  
  EBVtot<-data.frame(V1=solutions$V4[solutions$V1==1 & solutions$V2==2],
                     V2=solutions$V4[solutions$V1==1 & solutions$V2==3],
                     V3=solutions$V4[solutions$V1==3 & solutions$V2==2],
                     V4=solutions$V4[solutions$V1==4 & solutions$V2==3],
                     V5=solutions$V4[solutions$V1==2 & solutions$V2==2],
                     V6=solutions$V4[solutions$V1==2 & solutions$V2==3],
                     V7=solutions$V4[solutions$V1==4 & solutions$V2==2],
                     V8=solutions$V4[solutions$V1==3 & solutions$V2==3])
  
  newgeneaEE<-data.frame(V1=c(1:1000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="EE")
  newgeneaRR<-data.frame(V1=c(1001:2000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="RR")
  newgeneaER<-data.frame(V1=c(2001:2500)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="ER")
  newgeneaRE<-data.frame(V1=c(2501:3000)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="RE")
  
  #escoger los reproductores EE
  machosEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==1]) 
  machosEE$V2<-EBVtot$V3[(machosEE$V1)*2-1]+EBVtot$V3[(machosEE$V1)*2]+EBVtot$V4[(machosEE$V1)*2-1]+EBVtot$V4[(machosEE$V1)*2]
  machosEE<-machosEE$V1[order(machosEE$V2,decreasing = T)]
  
  hembrasEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==2])
  hembrasEE$V2<-EBVtot$V3[(hembrasEE$V1)*2-1]+EBVtot$V3[(hembrasEE$V1)*2]+EBVtot$V4[(hembrasEE$V1)*2-1]+EBVtot$V4[(hembrasEE$V1)*2]
  hembrasEE<-hembrasEE$V1[order(hembrasEE$V2,decreasing = T)]
  
  newgeneaEE$V2<-sample(machosEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo los 100 mejores machos de 500
  newgeneaEE$V3<-sample(hembrasEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo las 100 mejores hembras de 500
  
  
  #escoger los reproductores RR
  machosRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==1]) 
  machosRR$V2<-EBVtot$V7[(machosRR$V1)*2-1]+EBVtot$V7[(machosRR$V1)*2]+EBVtot$V8[(machosRR$V1)*2-1]+EBVtot$V8[(machosRR$V1)*2]
  machosRR<-machosRR$V1[order(machosRR$V2,decreasing = T)]
  
  hembrasRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==2])
  hembrasRR$V2<-EBVtot$V7[(hembrasRR$V1)*2-1]+EBVtot$V7[(hembrasRR$V1)*2]+EBVtot$V8[(hembrasRR$V1)*2-1]+EBVtot$V8[(hembrasRR$V1)*2]
  hembrasRR<-hembrasRR$V1[order(hembrasRR$V2,decreasing = T)]
  
  newgeneaRR$V2<-sample(machosRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo los 100 mejores machos de 500
  newgeneaRR$V3<-sample(hembrasRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo las 100 mejores hembras de 500
  
  #animales cruzados
  
  newgeneaER$V2<-sample(newgeneaEE$V1[newgeneaEE$V4==1],nrow(newgeneaER),replace = T)
  newgeneaER$V3<-sample(newgeneaRR$V1[newgeneaRR$V4==2],nrow(newgeneaER),replace = T)
  
  newgeneaRE$V3<-sample(newgeneaEE$V1[newgeneaEE$V4==2],nrow(newgeneaER),replace = T)
  newgeneaRE$V2<-sample(newgeneaRR$V1[newgeneaRR$V4==1],nrow(newgeneaER),replace = T)
  
  geneaEE<-newgeneaEE
  geneaRR<-newgeneaRR
  geneaER<-newgeneaER
  geneaRE<-newgeneaRE
  
  #genea completa
  
  geneatot<-rbind(geneatot,newgeneaEE,newgeneaRR,newgeneaER,newgeneaRE)
  
  TBVEnew<-as.data.frame(matrix(0,nrow = nrow(geneaEE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVEnew)) %% 2
  TBVEnew[filas==1,]<-TBVtot[(geneaEE$V2*2)-1,]*0.5+TBVtot[geneaEE$V2*2,]*0.5
  TBVEnew[filas==0,]<-TBVtot[(geneaEE$V3*2)-1,]*0.5+TBVtot[geneaEE$V3*2,]*0.5
  TBVEnew<-TBVEnew+(mvrnorm(nrow(geneaEE)*2,medias,varcovartot/2))
  
  TBVRnew<-as.data.frame(matrix(0,nrow = nrow(geneaRR)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRnew)) %% 2
  TBVRnew[filas==1,]<-TBVtot[(geneaRR$V2*2)-1,]*0.5+TBVtot[geneaRR$V2*2,]*0.5
  TBVRnew[filas==0,]<-TBVtot[(geneaRR$V3*2)-1,]*0.5+TBVtot[geneaRR$V3*2,]*0.5
  TBVRnew<-TBVRnew+(mvrnorm(nrow(geneaRR)*2,medias,varcovartot/2))
  
  TBVE<-TBVEnew
  TBVR<-TBVRnew
  
  TBVtot<-rbind(TBVtot,TBVE,TBVR)
  
  TBVER<-as.data.frame(matrix(0,nrow = nrow(geneaER)*2,ncol = 8))
  filas<-seq_len(nrow(TBVER)) %% 2
  TBVER[filas==1,]<-TBVtot[(geneaER$V2*2)-1,]*0.5+TBVtot[geneaER$V2*2,]*0.5
  TBVER[filas==0,]<-TBVtot[(geneaER$V3*2)-1,]*0.5+TBVtot[geneaER$V3*2,]*0.5
  TBVER<-TBVER+(mvrnorm(nrow(geneaER)*2,medias,varcovartot/2))
  
  TBVRE<-as.data.frame(matrix(0,nrow = nrow(geneaRE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRE)) %% 2
  TBVRE[filas==1,]<-TBVtot[(geneaRE$V2*2)-1,]*0.5+TBVtot[geneaRE$V2*2,]*0.5
  TBVRE[filas==0,]<-TBVtot[(geneaRE$V3*2)-1,]*0.5+TBVtot[geneaRE$V3*2,]*0.5
  TBVRE<-TBVRE+(mvrnorm(nrow(geneaRE)*2,medias,varcovartot/2))
  
  TBVtot<-rbind(TBVtot,TBVER,TBVRE)
  
  #datos
  geneaEEtot<-geneatot[geneatot$Raza=="EE",]
  geneaRRtot<-geneatot[geneatot$Raza=="RR",]
  geneaERtot<-geneatot[geneatot$Raza=="ER",]
  geneaREtot<-geneatot[geneatot$Raza=="RE",]
  
  datosEEmachos<-data.frame(V1=(geneaEEtot$V1[geneaEEtot$V4==1]*2)-1,V2=geneaEEtot$V1[geneaEEtot$V4==1]*2,V3=1,V4=10,V5=0,V6=0,V7=0)
  datosEEmachos$V4<-TBVtot$V1[datosEEmachos$V1]+TBVtot$V2[datosEEmachos$V2]
  
  datosEEhembras<-data.frame(V1=(geneaEEtot$V1[geneaEEtot$V4==2]*2)-1,V2=geneaEEtot$V1[geneaEEtot$V4==2]*2,V3=1,V4=10,V5=0,V6=0,V7=0)
  datosEEhembras$V4<-TBVtot$V1[datosEEhembras$V1]+TBVtot$V2[datosEEhembras$V2]
  
  datosRRmachos<-data.frame(V1=(geneaRRtot$V1[geneaRRtot$V4==1]*2)-1,V2=geneaRRtot$V1[geneaRRtot$V4==1]*2,V3=1,V4=0,V5=10,V6=0,V7=0)
  datosRRmachos$V5<-TBVtot$V5[datosRRmachos$V1]+TBVtot$V6[datosRRmachos$V2]
  
  datosRRhembras<-data.frame(V1=(geneaRRtot$V1[geneaRRtot$V4==2]*2)-1,V2=geneaRRtot$V1[geneaRRtot$V4==2]*2,V3=1,V4=0,V5=10,V6=0,V7=0)
  datosRRhembras$V5<-TBVtot$V5[datosRRhembras$V1]+TBVtot$V6[datosRRhembras$V2]
  
  datosER<-data.frame(V1=(geneaERtot$V1*2)-1,V2=geneaERtot$V1*2,V3=1,V4=0,V5=0,V6=10,V7=0) #solo hay hembras
  datosER$V6<-TBVtot$V3[datosER$V1]+TBVtot$V8[datosER$V2]
  
  datosRE<-data.frame(V1=(geneaREtot$V1*2)-1,V2=geneaREtot$V1*2,V3=1,V4=0,V5=0,V6=0,V7=10) #solo hay hembras
  datosRE$V7<-TBVtot$V7[datosRE$V1]+TBVtot$V4[datosRE$V2]
  
  datostot<-rbind(datosEEmachos,datosEEhembras,datosRRmachos,datosRRhembras,datosER,datosRE)
  
}

save.image("selEBVmediacruce.R")
load("selEBVmediacruce.R")

geneaEEtot<-geneatot[geneatot$Raza=="EE",]
geneaRRtot<-geneatot[geneatot$Raza=="RR",]
geneaERtot<-geneatot[geneatot$Raza=="ER",]
geneaREtot<-geneatot[geneatot$Raza=="RE",]

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaEEprint$TBVmatER<-TBVtot$V4[(geneaEEprint$V1*2)-1]+TBVtot$V4[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaRRprint$TBVmatER<-TBVtot$V8[(geneaRRprint$V1*2)-1]+TBVtot$V8[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatER",colBy = "Gen")
plot(AlphaPart,ylim = c(-1,10),color = c(2,3))

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaEEprint$TBVmatRE<-TBVtot$V3[(geneaEEprint$V1*2)-1]+TBVtot$V3[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaRRprint$TBVmatRE<-TBVtot$V7[(geneaRRprint$V1*2)-1]+TBVtot$V7[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatRE",colBy = "Gen")
plot(AlphaPart,ylim = c(-1,10),color = c(2,3))

#####
#TBV#
#####

#Selección por raza
resultsEE1<-matrix(ncol=12,nrow=100)
resultsRR1<-matrix(ncol=12,nrow=100)
resultsSUM1<-matrix(ncol=12,nrow=100)

resultsEE2<-matrix(ncol=12,nrow=100)
resultsRR2<-matrix(ncol=12,nrow=100)
resultsSUM2<-matrix(ncol=12,nrow=100)

prodER<-matrix(ncol=1,nrow=100)
prodRE<-matrix(ncol=1,nrow=100)


for (rep in 1:100) {
  

load("pobbase.RData")
  set.seed(Sys.time())

for (Generation in 21:30) {

  newgeneaEE<-data.frame(V1=c(1:1000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="EE")
  newgeneaRR<-data.frame(V1=c(1001:2000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="RR")
  newgeneaER<-data.frame(V1=c(2001:2500)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="ER")
  newgeneaRE<-data.frame(V1=c(2501:3000)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="RE")
  
  #escoger los reproductores EE
  machosEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==1]) 
  machosEE$V2<-TBVtot$V1[(machosEE$V1)*2-1]+TBVtot$V2[(machosEE$V1)*2]
  machosEE<-machosEE$V1[order(machosEE$V2,decreasing = T)]
  
  hembrasEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==2])
  hembrasEE$V2<-TBVtot$V1[(hembrasEE$V1)*2-1]+TBVtot$V2[(hembrasEE$V1)*2]
  hembrasEE<-hembrasEE$V1[order(hembrasEE$V2,decreasing = T)]
  
  newgeneaEE$V2<-sample(machosEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo los 100 mejores machos de 500
  newgeneaEE$V3<-sample(hembrasEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo las 100 mejores hembras de 500
  
  
  #escoger los reproductores RR
  machosRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==1]) 
  machosRR$V2<-TBVtot$V1[(machosRR$V1)*2-1]+TBVtot$V2[(machosRR$V1)*2]
  machosRR<-machosRR$V1[order(machosRR$V2,decreasing = T)]
  
  hembrasRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==2])
  hembrasRR$V2<-TBVtot$V1[(hembrasRR$V1)*2-1]+TBVtot$V2[(hembrasRR$V1)*2]
  hembrasRR<-hembrasRR$V1[order(hembrasRR$V2,decreasing = T)]
  
  newgeneaRR$V2<-sample(machosRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo los 100 mejores machos de 500
  newgeneaRR$V3<-sample(hembrasRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo las 100 mejores hembras de 500
  
  #animales cruzados
  
  newgeneaER$V2<-sample(newgeneaEE$V1[newgeneaEE$V4==1],nrow(newgeneaER),replace = T)
  newgeneaER$V3<-sample(newgeneaRR$V1[newgeneaRR$V4==2],nrow(newgeneaER),replace = T)
  
  newgeneaRE$V3<-sample(newgeneaEE$V1[newgeneaEE$V4==2],nrow(newgeneaER),replace = T)
  newgeneaRE$V2<-sample(newgeneaRR$V1[newgeneaRR$V4==1],nrow(newgeneaER),replace = T)
  
  geneaEE<-newgeneaEE
  geneaRR<-newgeneaRR
  geneaER<-newgeneaER
  geneaRE<-newgeneaRE
  
  #genea completa
  
  geneatot<-rbind(geneatot,newgeneaEE,newgeneaRR,newgeneaER,newgeneaRE)
  
  TBVEnew<-as.data.frame(matrix(0,nrow = nrow(geneaEE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVEnew)) %% 2
  TBVEnew[filas==1,]<-TBVtot[(geneaEE$V2*2)-1,]*0.5+TBVtot[geneaEE$V2*2,]*0.5
  TBVEnew[filas==0,]<-TBVtot[(geneaEE$V3*2)-1,]*0.5+TBVtot[geneaEE$V3*2,]*0.5
  TBVEnew<-TBVEnew+(mvrnorm(nrow(geneaEE)*2,medias,varcovartot/2))
  
  TBVRnew<-as.data.frame(matrix(0,nrow = nrow(geneaRR)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRnew)) %% 2
  TBVRnew[filas==1,]<-TBVtot[(geneaRR$V2*2)-1,]*0.5+TBVtot[geneaRR$V2*2,]*0.5
  TBVRnew[filas==0,]<-TBVtot[(geneaRR$V3*2)-1,]*0.5+TBVtot[geneaRR$V3*2,]*0.5
  TBVRnew<-TBVRnew+(mvrnorm(nrow(geneaRR)*2,medias,varcovartot/2))
  
  TBVE<-TBVEnew
  TBVR<-TBVRnew
  
  TBVtot<-rbind(TBVtot,TBVE,TBVR)
  
  TBVER<-as.data.frame(matrix(0,nrow = nrow(geneaER)*2,ncol = 8))
  filas<-seq_len(nrow(TBVER)) %% 2
  TBVER[filas==1,]<-TBVtot[(geneaER$V2*2)-1,]*0.5+TBVtot[geneaER$V2*2,]*0.5
  TBVER[filas==0,]<-TBVtot[(geneaER$V3*2)-1,]*0.5+TBVtot[geneaER$V3*2,]*0.5
  TBVER<-TBVER+(mvrnorm(nrow(geneaER)*2,medias,varcovartot/2))
  
  TBVRE<-as.data.frame(matrix(0,nrow = nrow(geneaRE)*2,ncol = 8))
  filas<-seq_len(nrow(TBVRE)) %% 2
  TBVRE[filas==1,]<-TBVtot[(geneaRE$V2*2)-1,]*0.5+TBVtot[geneaRE$V2*2,]*0.5
  TBVRE[filas==0,]<-TBVtot[(geneaRE$V3*2)-1,]*0.5+TBVtot[geneaRE$V3*2,]*0.5
  TBVRE<-TBVRE+(mvrnorm(nrow(geneaRE)*2,medias,varcovartot/2))
  
  TBVtot<-rbind(TBVtot,TBVER,TBVRE)
  
}

geneaEEtot<-geneatot[geneatot$Raza=="EE",]
geneaRRtot<-geneatot[geneatot$Raza=="RR",]

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaEEprint$TBVmatER<-TBVtot$V4[(geneaEEprint$V1*2)-1]+TBVtot$V4[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaRRprint$TBVmatER<-TBVtot$V8[(geneaRRprint$V1*2)-1]+TBVtot$V8[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatER",colBy = "Gen")
resultsEE1[rep,]<-AlphaPart[["TBVmatER"]][["EE"]]
resultsRR1[rep,]<-AlphaPart[["TBVmatER"]][["RR"]]
resultsSUM1[rep,]<-AlphaPart[["TBVmatER"]][["Sum"]]
#plot(AlphaPart,ylim = c(-1,10),color = c(2,3))

geneaEEprint<-geneaEEtot
geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaEEprint$TBVmatRE<-TBVtot$V3[(geneaEEprint$V1*2)-1]+TBVtot$V3[(geneaEEprint$V1*2)] 

geneaRRprint<-geneaRRtot
geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
geneaRRprint$TBVmatRE<-TBVtot$V7[(geneaRRprint$V1*2)-1]+TBVtot$V7[(geneaRRprint$V1*2)] 

geneatotprint<-rbind(geneaEEprint,geneaRRprint)

AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatRE",colBy = "Gen")
resultsEE2[rep,]<-AlphaPart[["TBVmatRE"]][["EE"]]
resultsRR2[rep,]<-AlphaPart[["TBVmatRE"]][["RR"]]
resultsSUM2[rep,]<-AlphaPart[["TBVmatRE"]][["Sum"]]
#plot(AlphaPart,ylim = c(-1,10),color = c(2,3))

prodER[rep]<-TBVER$V3+TBVER$V8
prodRE[rep]<-TBVRE$V4+TBVRE$V7
}

resumen1<-rbind(data.frame(V1=20:31,V2=apply(resultsEE1,2,mean),
                    V3=apply(resultsEE1,2,sd),V4="E"),
                data.frame(V1=20:31,V2=apply(resultsRR1,2,mean),
                    V3=apply(resultsRR1,2,sd),V4="R"),
                data.frame(V1=20:31,V2=apply(resultsSUM1,2,mean),
                V3=apply(resultsSUM1,2,sd),V4="SUM"))


resumen2<-rbind(data.frame(V1=20:31,V2=apply(resultsEE2,2,mean),
                           V3=apply(resultsEE2,2,sd),V4="E"),
                data.frame(V1=20:31,V2=apply(resultsRR2,2,mean),
                           V3=apply(resultsRR2,2,sd),V4="R"),
                data.frame(V1=20:31,V2=apply(resultsSUM2,2,mean),
                           V3=apply(resultsSUM2,2,sd),V4="SUM"))

ERplot<-ggplot(resumen1,aes(x=V1, y=V2,fill=V4,color=V4)) +
  geom_line(show.legend = F)+
  geom_ribbon(aes(x=V1,ymin=V2-V3*2,ymax=V2+V3*2),alpha=.15)+
  theme(legend.background = element_rect(fill = "lightblue",colour = 1))+
  lims(y=c(-1,10.5))+
  labs(x="Generacion", y="TBV")

REplot<-ggplot(resumen1,aes(x=V1, y=V2,fill=V4,color=V4)) +
  geom_line(show.legend = F)+
  geom_ribbon(aes(x=V1,ymin=V2-V3*2,ymax=V2+V3*2),alpha=.15)+
  theme(legend.background = element_rect(fill = "lightblue",colour = 1))+
  lims(y=c(-1,10.5))+
  labs(x="Generacion", y="TBV")

(ERplot + REplot + plot_layout(guides="collect"))

prod<-rbind(data.frame(V1="ER",V2=prodER),
            data.frame(V1="RE",V2=prodRE))
prodplot<-ggplot(prod)+
  geom_violin(aes(x=V1,y=V2,fill=V1),show.legend = F)+
  labs(x="Cruce", y="TBV")+
  lims(y=c(-1,6))

prodplot


#Selección por raza EE y hembra RR

resultsEE1<-matrix(ncol=12,nrow=100)
resultsRR1<-matrix(ncol=12,nrow=100)
resultsSUM1<-matrix(ncol=12,nrow=100)

resultsEE2<-matrix(ncol=12,nrow=100)
resultsRR2<-matrix(ncol=12,nrow=100)
resultsSUM2<-matrix(ncol=12,nrow=100)

prodER<-matrix(ncol=1,nrow=100)
prodRE<-matrix(ncol=1,nrow=100)


for (rep in 1:100) {
  
  
  load("pobbase.RData")
  set.seed(Sys.time())
  
  for (Generation in 21:30) {
    
    newgeneaEE<-data.frame(V1=c(1:1000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="EE")
    newgeneaRR<-data.frame(V1=c(1001:2000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="RR")
    newgeneaER<-data.frame(V1=c(2001:2500)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="ER")
    newgeneaRE<-data.frame(V1=c(2501:3000)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="RE")
    
    #escoger los reproductores EE
    machosEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==1]) 
    machosEE$V2<-TBVtot$V1[(machosEE$V1)*2-1]+TBVtot$V2[(machosEE$V1)*2]
    machosEE<-machosEE$V1[order(machosEE$V2,decreasing = T)]
    
    hembrasEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==2])
    hembrasEE$V2<-TBVtot$V1[(hembrasEE$V1)*2-1]+TBVtot$V2[(hembrasEE$V1)*2]
    hembrasEE<-hembrasEE$V1[order(hembrasEE$V2,decreasing = T)]
    
    newgeneaEE$V2<-sample(machosEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo los 100 mejores machos de 500
    newgeneaEE$V3<-sample(hembrasEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo las 100 mejores hembras de 500
    
    
    #escoger los reproductores RR
    machosRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==1]) 
    machosRR$V2<-TBVtot$V2[(machosRR$V1)*2-1]+TBVtot$V2[(machosRR$V1)*2]
    machosRR<-machosRR$V1[order(machosRR$V2,decreasing = T)]
    
    hembrasRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==2])
    hembrasRR$V2<-TBVtot$V2[(hembrasRR$V1)*2-1]+TBVtot$V2[(hembrasRR$V1)*2]
    hembrasRR<-hembrasRR$V1[order(hembrasRR$V2,decreasing = T)]
    
    newgeneaRR$V2<-sample(machosRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo los 100 mejores machos de 500
    newgeneaRR$V3<-sample(hembrasRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo las 100 mejores hembras de 500
    
    #animales cruzados
    
    newgeneaER$V2<-sample(newgeneaEE$V1[newgeneaEE$V4==1],nrow(newgeneaER),replace = T)
    newgeneaER$V3<-sample(newgeneaRR$V1[newgeneaRR$V4==2],nrow(newgeneaER),replace = T)
    
    newgeneaRE$V3<-sample(newgeneaEE$V1[newgeneaEE$V4==2],nrow(newgeneaER),replace = T)
    newgeneaRE$V2<-sample(newgeneaRR$V1[newgeneaRR$V4==1],nrow(newgeneaER),replace = T)
    
    geneaEE<-newgeneaEE
    geneaRR<-newgeneaRR
    geneaER<-newgeneaER
    geneaRE<-newgeneaRE
    
    #genea completa
    
    geneatot<-rbind(geneatot,newgeneaEE,newgeneaRR,newgeneaER,newgeneaRE)
    
    TBVEnew<-as.data.frame(matrix(0,nrow = nrow(geneaEE)*2,ncol = 8))
    filas<-seq_len(nrow(TBVEnew)) %% 2
    TBVEnew[filas==1,]<-TBVtot[(geneaEE$V2*2)-1,]*0.5+TBVtot[geneaEE$V2*2,]*0.5
    TBVEnew[filas==0,]<-TBVtot[(geneaEE$V3*2)-1,]*0.5+TBVtot[geneaEE$V3*2,]*0.5
    TBVEnew<-TBVEnew+(mvrnorm(nrow(geneaEE)*2,medias,varcovartot/2))
    
    TBVRnew<-as.data.frame(matrix(0,nrow = nrow(geneaRR)*2,ncol = 8))
    filas<-seq_len(nrow(TBVRnew)) %% 2
    TBVRnew[filas==1,]<-TBVtot[(geneaRR$V2*2)-1,]*0.5+TBVtot[geneaRR$V2*2,]*0.5
    TBVRnew[filas==0,]<-TBVtot[(geneaRR$V3*2)-1,]*0.5+TBVtot[geneaRR$V3*2,]*0.5
    TBVRnew<-TBVRnew+(mvrnorm(nrow(geneaRR)*2,medias,varcovartot/2))
    
    TBVE<-TBVEnew
    TBVR<-TBVRnew
    
    TBVtot<-rbind(TBVtot,TBVE,TBVR)
    
    TBVER<-as.data.frame(matrix(0,nrow = nrow(geneaER)*2,ncol = 8))
    filas<-seq_len(nrow(TBVER)) %% 2
    TBVER[filas==1,]<-TBVtot[(geneaER$V2*2)-1,]*0.5+TBVtot[geneaER$V2*2,]*0.5
    TBVER[filas==0,]<-TBVtot[(geneaER$V3*2)-1,]*0.5+TBVtot[geneaER$V3*2,]*0.5
    TBVER<-TBVER+(mvrnorm(nrow(geneaER)*2,medias,varcovartot/2))
    
    TBVRE<-as.data.frame(matrix(0,nrow = nrow(geneaRE)*2,ncol = 8))
    filas<-seq_len(nrow(TBVRE)) %% 2
    TBVRE[filas==1,]<-TBVtot[(geneaRE$V2*2)-1,]*0.5+TBVtot[geneaRE$V2*2,]*0.5
    TBVRE[filas==0,]<-TBVtot[(geneaRE$V3*2)-1,]*0.5+TBVtot[geneaRE$V3*2,]*0.5
    TBVRE<-TBVRE+(mvrnorm(nrow(geneaRE)*2,medias,varcovartot/2))
    
    TBVtot<-rbind(TBVtot,TBVER,TBVRE)
    
  }
  
  geneaEEtot<-geneatot[geneatot$Raza=="EE",]
  geneaRRtot<-geneatot[geneatot$Raza=="RR",]
  
  geneaEEprint<-geneaEEtot
  geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaEEprint$TBVmatER<-TBVtot$V4[(geneaEEprint$V1*2)-1]+TBVtot$V4[(geneaEEprint$V1*2)] 
  
  geneaRRprint<-geneaRRtot
  geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaRRprint$TBVmatER<-TBVtot$V8[(geneaRRprint$V1*2)-1]+TBVtot$V8[(geneaRRprint$V1*2)] 
  
  geneatotprint<-rbind(geneaEEprint,geneaRRprint)
  
  AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatER",colBy = "Gen")
  resultsEE1[rep,]<-AlphaPart[["TBVmatER"]][["EE"]]
  resultsRR1[rep,]<-AlphaPart[["TBVmatER"]][["RR"]]
  resultsSUM1[rep,]<-AlphaPart[["TBVmatER"]][["Sum"]]
  #plot(AlphaPart,ylim = c(-1,10),color = c(2,3))
  
  geneaEEprint<-geneaEEtot
  geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaEEprint$TBVmatRE<-TBVtot$V3[(geneaEEprint$V1*2)-1]+TBVtot$V3[(geneaEEprint$V1*2)] 
  
  geneaRRprint<-geneaRRtot
  geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaRRprint$TBVmatRE<-TBVtot$V7[(geneaRRprint$V1*2)-1]+TBVtot$V7[(geneaRRprint$V1*2)] 
  
  geneatotprint<-rbind(geneaEEprint,geneaRRprint)
  
  AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatRE",colBy = "Gen")
  resultsEE2[rep,]<-AlphaPart[["TBVmatRE"]][["EE"]]
  resultsRR2[rep,]<-AlphaPart[["TBVmatRE"]][["RR"]]
  resultsSUM2[rep,]<-AlphaPart[["TBVmatRE"]][["Sum"]]
  #plot(AlphaPart,ylim = c(-1,10),color = c(2,3))
  
  prodER[rep]<-TBVER$V3+TBVER$V8
  prodRE[rep]<-TBVRE$V4+TBVRE$V7
}

resumen1<-rbind(data.frame(V1=20:31,V2=apply(resultsEE1,2,mean),
                           V3=apply(resultsEE1,2,sd),V4="E"),
                data.frame(V1=20:31,V2=apply(resultsRR1,2,mean),
                           V3=apply(resultsRR1,2,sd),V4="R"),
                data.frame(V1=20:31,V2=apply(resultsSUM1,2,mean),
                           V3=apply(resultsSUM1,2,sd),V4="SUM"))


resumen2<-rbind(data.frame(V1=20:31,V2=apply(resultsEE2,2,mean),
                           V3=apply(resultsEE2,2,sd),V4="E"),
                data.frame(V1=20:31,V2=apply(resultsRR2,2,mean),
                           V3=apply(resultsRR2,2,sd),V4="R"),
                data.frame(V1=20:31,V2=apply(resultsSUM2,2,mean),
                           V3=apply(resultsSUM2,2,sd),V4="SUM"))

ERplot<-ggplot(resumen1,aes(x=V1, y=V2,fill=V4,color=V4)) +
  geom_line(show.legend = F)+
  geom_ribbon(aes(x=V1,ymin=V2-V3*2,ymax=V2+V3*2),alpha=.15)+
  theme(legend.background = element_rect(fill = "lightblue",colour = 1))+
  lims(y=c(-1,10.5))+
  labs(x="Generacion", y="TBV")

REplot<-ggplot(resumen1,aes(x=V1, y=V2,fill=V4,color=V4)) +
  geom_line(show.legend = F)+
  geom_ribbon(aes(x=V1,ymin=V2-V3*2,ymax=V2+V3*2),alpha=.15)+
  theme(legend.background = element_rect(fill = "lightblue",colour = 1))+
  lims(y=c(-1,10.5))+
  labs(x="Generacion", y="TBV")

(ERplot + REplot + plot_layout(guides="collect"))

prod<-rbind(data.frame(V1="ER",V2=prodER),
            data.frame(V1="RE",V2=prodRE))
prodplot<-ggplot(prod)+
  geom_violin(aes(x=V1,y=V2,fill=V1),show.legend = F)+
  labs(x="Cruce", y="TBV")+
  lims(y=c(-1,6))

prodplot

#Selección por cruce

resultsEE1<-matrix(ncol=12,nrow=100)
resultsRR1<-matrix(ncol=12,nrow=100)
resultsSUM1<-matrix(ncol=12,nrow=100)

resultsEE2<-matrix(ncol=12,nrow=100)
resultsRR2<-matrix(ncol=12,nrow=100)
resultsSUM2<-matrix(ncol=12,nrow=100)

prodER<-matrix(ncol=1,nrow=100)
prodRE<-matrix(ncol=1,nrow=100)


for (rep in 1:100) {
  
  
  load("pobbase.RData")
  set.seed(Sys.time())
  
  for (Generation in 21:30) {
    
    newgeneaEE<-data.frame(V1=c(1:1000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="EE")
    newgeneaRR<-data.frame(V1=c(1001:2000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="RR")
    newgeneaER<-data.frame(V1=c(2001:2500)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="ER")
    newgeneaRE<-data.frame(V1=c(2501:3000)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="RE")
    
    #escoger los reproductores EE
    machosEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==1]) 
    machosEE$V2<-TBVtot$V3[(machosEE$V1)*2-1]+TBVtot$V3[(machosEE$V1)*2]+TBVtot$V4[(machosEE$V1)*2-1]+TBVtot$V4[(machosEE$V1)*2]
    machosEE<-machosEE$V1[order(machosEE$V2,decreasing = T)]
    
    hembrasEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==2])
    hembrasEE$V2<-TBVtot$V3[(hembrasEE$V1)*2-1]+TBVtot$V3[(hembrasEE$V1)*2]+TBVtot$V4[(hembrasEE$V1)*2-1]+TBVtot$V4[(hembrasEE$V1)*2]
    hembrasEE<-hembrasEE$V1[order(hembrasEE$V2,decreasing = T)]
    
    newgeneaEE$V2<-sample(machosEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo los 100 mejores machos de 500
    newgeneaEE$V3<-sample(hembrasEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo las 100 mejores hembras de 500
    
    
    #escoger los reproductores RR
    machosRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==1]) 
    machosRR$V2<-TBVtot$V7[(machosRR$V1)*2-1]+TBVtot$V7[(machosRR$V1)*2]+TBVtot$V8[(machosRR$V1)*2-1]+TBVtot$V8[(machosRR$V1)*2]
    machosRR<-machosRR$V1[order(machosRR$V2,decreasing = T)]
    
    hembrasRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==2])
    hembrasRR$V2<-TBVtot$V7[(hembrasRR$V1)*2-1]+TBVtot$V7[(hembrasRR$V1)*2]+TBVtot$V8[(hembrasRR$V1)*2-1]+TBVtot$V8[(hembrasRR$V1)*2]
    hembrasRR<-hembrasRR$V1[order(hembrasRR$V2,decreasing = T)]
    
    newgeneaRR$V2<-sample(machosRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo los 100 mejores machos de 500
    newgeneaRR$V3<-sample(hembrasRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo las 100 mejores hembras de 500
    
    #animales cruzados
    
    newgeneaER$V2<-sample(newgeneaEE$V1[newgeneaEE$V4==1],nrow(newgeneaER),replace = T)
    newgeneaER$V3<-sample(newgeneaRR$V1[newgeneaRR$V4==2],nrow(newgeneaER),replace = T)
    
    newgeneaRE$V3<-sample(newgeneaEE$V1[newgeneaEE$V4==2],nrow(newgeneaER),replace = T)
    newgeneaRE$V2<-sample(newgeneaRR$V1[newgeneaRR$V4==1],nrow(newgeneaER),replace = T)
    
    geneaEE<-newgeneaEE
    geneaRR<-newgeneaRR
    geneaER<-newgeneaER
    geneaRE<-newgeneaRE
    
    #genea completa
    
    geneatot<-rbind(geneatot,newgeneaEE,newgeneaRR,newgeneaER,newgeneaRE)
    
    TBVEnew<-as.data.frame(matrix(0,nrow = nrow(geneaEE)*2,ncol = 8))
    filas<-seq_len(nrow(TBVEnew)) %% 2
    TBVEnew[filas==1,]<-TBVtot[(geneaEE$V2*2)-1,]*0.5+TBVtot[geneaEE$V2*2,]*0.5
    TBVEnew[filas==0,]<-TBVtot[(geneaEE$V3*2)-1,]*0.5+TBVtot[geneaEE$V3*2,]*0.5
    TBVEnew<-TBVEnew+(mvrnorm(nrow(geneaEE)*2,medias,varcovartot/2))
    
    TBVRnew<-as.data.frame(matrix(0,nrow = nrow(geneaRR)*2,ncol = 8))
    filas<-seq_len(nrow(TBVRnew)) %% 2
    TBVRnew[filas==1,]<-TBVtot[(geneaRR$V2*2)-1,]*0.5+TBVtot[geneaRR$V2*2,]*0.5
    TBVRnew[filas==0,]<-TBVtot[(geneaRR$V3*2)-1,]*0.5+TBVtot[geneaRR$V3*2,]*0.5
    TBVRnew<-TBVRnew+(mvrnorm(nrow(geneaRR)*2,medias,varcovartot/2))
    
    TBVE<-TBVEnew
    TBVR<-TBVRnew
    
    TBVtot<-rbind(TBVtot,TBVE,TBVR)
    
    TBVER<-as.data.frame(matrix(0,nrow = nrow(geneaER)*2,ncol = 8))
    filas<-seq_len(nrow(TBVER)) %% 2
    TBVER[filas==1,]<-TBVtot[(geneaER$V2*2)-1,]*0.5+TBVtot[geneaER$V2*2,]*0.5
    TBVER[filas==0,]<-TBVtot[(geneaER$V3*2)-1,]*0.5+TBVtot[geneaER$V3*2,]*0.5
    TBVER<-TBVER+(mvrnorm(nrow(geneaER)*2,medias,varcovartot/2))
    
    TBVRE<-as.data.frame(matrix(0,nrow = nrow(geneaRE)*2,ncol = 8))
    filas<-seq_len(nrow(TBVRE)) %% 2
    TBVRE[filas==1,]<-TBVtot[(geneaRE$V2*2)-1,]*0.5+TBVtot[geneaRE$V2*2,]*0.5
    TBVRE[filas==0,]<-TBVtot[(geneaRE$V3*2)-1,]*0.5+TBVtot[geneaRE$V3*2,]*0.5
    TBVRE<-TBVRE+(mvrnorm(nrow(geneaRE)*2,medias,varcovartot/2))
    
    TBVtot<-rbind(TBVtot,TBVER,TBVRE)
    
  }
  
  geneaEEtot<-geneatot[geneatot$Raza=="EE",]
  geneaRRtot<-geneatot[geneatot$Raza=="RR",]
  
  geneaEEprint<-geneaEEtot
  geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaEEprint$TBVmatER<-TBVtot$V4[(geneaEEprint$V1*2)-1]+TBVtot$V4[(geneaEEprint$V1*2)] 
  
  geneaRRprint<-geneaRRtot
  geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaRRprint$TBVmatER<-TBVtot$V8[(geneaRRprint$V1*2)-1]+TBVtot$V8[(geneaRRprint$V1*2)] 
  
  geneatotprint<-rbind(geneaEEprint,geneaRRprint)
  
  AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatER",colBy = "Gen")
  resultsEE1[rep,]<-AlphaPart[["TBVmatER"]][["EE"]]
  resultsRR1[rep,]<-AlphaPart[["TBVmatER"]][["RR"]]
  resultsSUM1[rep,]<-AlphaPart[["TBVmatER"]][["Sum"]]
  #plot(AlphaPart,ylim = c(-1,10),color = c(2,3))
  
  geneaEEprint<-geneaEEtot
  geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaEEprint$TBVmatRE<-TBVtot$V3[(geneaEEprint$V1*2)-1]+TBVtot$V3[(geneaEEprint$V1*2)] 
  
  geneaRRprint<-geneaRRtot
  geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaRRprint$TBVmatRE<-TBVtot$V7[(geneaRRprint$V1*2)-1]+TBVtot$V7[(geneaRRprint$V1*2)] 
  
  geneatotprint<-rbind(geneaEEprint,geneaRRprint)
  
  AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatRE",colBy = "Gen")
  resultsEE2[rep,]<-AlphaPart[["TBVmatRE"]][["EE"]]
  resultsRR2[rep,]<-AlphaPart[["TBVmatRE"]][["RR"]]
  resultsSUM2[rep,]<-AlphaPart[["TBVmatRE"]][["Sum"]]
  #plot(AlphaPart,ylim = c(-1,10),color = c(2,3))
  
  prodER[rep]<-TBVER$V3+TBVER$V8
  prodRE[rep]<-TBVRE$V4+TBVRE$V7
}

resumen1<-rbind(data.frame(V1=20:31,V2=apply(resultsEE1,2,mean),
                           V3=apply(resultsEE1,2,sd),V4="E"),
                data.frame(V1=20:31,V2=apply(resultsRR1,2,mean),
                           V3=apply(resultsRR1,2,sd),V4="R"),
                data.frame(V1=20:31,V2=apply(resultsSUM1,2,mean),
                           V3=apply(resultsSUM1,2,sd),V4="SUM"))


resumen2<-rbind(data.frame(V1=20:31,V2=apply(resultsEE2,2,mean),
                           V3=apply(resultsEE2,2,sd),V4="E"),
                data.frame(V1=20:31,V2=apply(resultsRR2,2,mean),
                           V3=apply(resultsRR2,2,sd),V4="R"),
                data.frame(V1=20:31,V2=apply(resultsSUM2,2,mean),
                           V3=apply(resultsSUM2,2,sd),V4="SUM"))

ERplot<-ggplot(resumen1,aes(x=V1, y=V2,fill=V4,color=V4)) +
  geom_line(show.legend = F)+
  geom_ribbon(aes(x=V1,ymin=V2-V3*2,ymax=V2+V3*2),alpha=.15)+
  theme(legend.background = element_rect(fill = "lightblue",colour = 1))+
  lims(y=c(-1,10.5))+
  labs(x="Generacion", y="TBV")

REplot<-ggplot(resumen1,aes(x=V1, y=V2,fill=V4,color=V4)) +
  geom_line(show.legend = F)+
  geom_ribbon(aes(x=V1,ymin=V2-V3*2,ymax=V2+V3*2),alpha=.15)+
  theme(legend.background = element_rect(fill = "lightblue",colour = 1))+
  lims(y=c(-1,10.5))+
  labs(x="Generacion", y="TBV")

(ERplot + REplot + plot_layout(guides="collect"))

prod<-rbind(data.frame(V1="ER",V2=prodER),
            data.frame(V1="RE",V2=prodRE))
prodplot<-ggplot(prod)+
  geom_violin(aes(x=V1,y=V2,fill=V1),show.legend = F)+
  labs(x="Cruce", y="TBV")+
  lims(y=c(-1,6))

prodplot

#Selección por hembras cruce

resultsEE1<-matrix(ncol=12,nrow=100)
resultsRR1<-matrix(ncol=12,nrow=100)
resultsSUM1<-matrix(ncol=12,nrow=100)

resultsEE2<-matrix(ncol=12,nrow=100)
resultsRR2<-matrix(ncol=12,nrow=100)
resultsSUM2<-matrix(ncol=12,nrow=100)

prodER<-matrix(ncol=1,nrow=100)
prodRE<-matrix(ncol=1,nrow=100)


for (rep in 1:100) {
  
  
  load("pobbase.RData")
  set.seed(Sys.time())
  
  for (Generation in 21:30) {
    
    newgeneaEE<-data.frame(V1=c(1:1000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="EE")
    newgeneaRR<-data.frame(V1=c(1001:2000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="RR")
    newgeneaER<-data.frame(V1=c(2001:2500)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="ER")
    newgeneaRE<-data.frame(V1=c(2501:3000)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="RE")
    
    #escoger los reproductores EE
    machosEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==1]) 
    machosEE$V2<-TBVtot$V4[(machosEE$V1)*2-1]+TBVtot$V4[(machosEE$V1)*2]+TBVtot$V4[(machosEE$V1)*2-1]+TBVtot$V4[(machosEE$V1)*2]
    machosEE<-machosEE$V1[order(machosEE$V2,decreasing = T)]
    
    hembrasEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==2])
    hembrasEE$V2<-TBVtot$V1[(hembrasEE$V1)*2-1]+TBVtot$V4[(hembrasEE$V1)*2]+TBVtot$V4[(hembrasEE$V1)*2-1]+TBVtot$V4[(hembrasEE$V1)*2]
    hembrasEE<-hembrasEE$V1[order(hembrasEE$V2,decreasing = T)]
    
    newgeneaEE$V2<-sample(machosEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo los 100 mejores machos de 500
    newgeneaEE$V3<-sample(hembrasEE[c(1:100)],nrow(newgeneaEE),replace = T) #cojo las 100 mejores hembras de 500
    
    
    #escoger los reproductores RR
    machosRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==1]) 
    machosRR$V2<-TBVtot$V8[(machosRR$V1)*2-1]+TBVtot$V8[(machosRR$V1)*2]+TBVtot$V8[(machosRR$V1)*2-1]+TBVtot$V8[(machosRR$V1)*2]
    machosRR<-machosRR$V1[order(machosRR$V2,decreasing = T)]
    
    hembrasRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==2])
    hembrasRR$V2<-TBVtot$V8[(hembrasRR$V1)*2-1]+TBVtot$V8[(hembrasRR$V1)*2]+TBVtot$V8[(hembrasRR$V1)*2-1]+TBVtot$V8[(hembrasRR$V1)*2]
    hembrasRR<-hembrasRR$V1[order(hembrasRR$V2,decreasing = T)]
    
    newgeneaRR$V2<-sample(machosRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo los 100 mejores machos de 500
    newgeneaRR$V3<-sample(hembrasRR[c(1:100)],nrow(newgeneaRR),replace = T) #cojo las 100 mejores hembras de 500
    
    #animales cruzados
    
    newgeneaER$V2<-sample(newgeneaEE$V1[newgeneaEE$V4==1],nrow(newgeneaER),replace = T)
    newgeneaER$V3<-sample(newgeneaRR$V1[newgeneaRR$V4==2],nrow(newgeneaER),replace = T)
    
    newgeneaRE$V3<-sample(newgeneaEE$V1[newgeneaEE$V4==2],nrow(newgeneaER),replace = T)
    newgeneaRE$V2<-sample(newgeneaRR$V1[newgeneaRR$V4==1],nrow(newgeneaER),replace = T)
    
    geneaEE<-newgeneaEE
    geneaRR<-newgeneaRR
    geneaER<-newgeneaER
    geneaRE<-newgeneaRE
    
    #genea completa
    
    geneatot<-rbind(geneatot,newgeneaEE,newgeneaRR,newgeneaER,newgeneaRE)
    
    TBVEnew<-as.data.frame(matrix(0,nrow = nrow(geneaEE)*2,ncol = 8))
    filas<-seq_len(nrow(TBVEnew)) %% 2
    TBVEnew[filas==1,]<-TBVtot[(geneaEE$V2*2)-1,]*0.5+TBVtot[geneaEE$V2*2,]*0.5
    TBVEnew[filas==0,]<-TBVtot[(geneaEE$V3*2)-1,]*0.5+TBVtot[geneaEE$V3*2,]*0.5
    TBVEnew<-TBVEnew+(mvrnorm(nrow(geneaEE)*2,medias,varcovartot/2))
    
    TBVRnew<-as.data.frame(matrix(0,nrow = nrow(geneaRR)*2,ncol = 8))
    filas<-seq_len(nrow(TBVRnew)) %% 2
    TBVRnew[filas==1,]<-TBVtot[(geneaRR$V2*2)-1,]*0.5+TBVtot[geneaRR$V2*2,]*0.5
    TBVRnew[filas==0,]<-TBVtot[(geneaRR$V3*2)-1,]*0.5+TBVtot[geneaRR$V3*2,]*0.5
    TBVRnew<-TBVRnew+(mvrnorm(nrow(geneaRR)*2,medias,varcovartot/2))
    
    TBVE<-TBVEnew
    TBVR<-TBVRnew
    
    TBVtot<-rbind(TBVtot,TBVE,TBVR)
    
    TBVER<-as.data.frame(matrix(0,nrow = nrow(geneaER)*2,ncol = 8))
    filas<-seq_len(nrow(TBVER)) %% 2
    TBVER[filas==1,]<-TBVtot[(geneaER$V2*2)-1,]*0.5+TBVtot[geneaER$V2*2,]*0.5
    TBVER[filas==0,]<-TBVtot[(geneaER$V3*2)-1,]*0.5+TBVtot[geneaER$V3*2,]*0.5
    TBVER<-TBVER+(mvrnorm(nrow(geneaER)*2,medias,varcovartot/2))
    
    TBVRE<-as.data.frame(matrix(0,nrow = nrow(geneaRE)*2,ncol = 8))
    filas<-seq_len(nrow(TBVRE)) %% 2
    TBVRE[filas==1,]<-TBVtot[(geneaRE$V2*2)-1,]*0.5+TBVtot[geneaRE$V2*2,]*0.5
    TBVRE[filas==0,]<-TBVtot[(geneaRE$V3*2)-1,]*0.5+TBVtot[geneaRE$V3*2,]*0.5
    TBVRE<-TBVRE+(mvrnorm(nrow(geneaRE)*2,medias,varcovartot/2))
    
    TBVtot<-rbind(TBVtot,TBVER,TBVRE)
    
  }
  
  geneaEEtot<-geneatot[geneatot$Raza=="EE",]
  geneaRRtot<-geneatot[geneatot$Raza=="RR",]
  
  geneaEEprint<-geneaEEtot
  geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaEEprint$TBVmatER<-TBVtot$V4[(geneaEEprint$V1*2)-1]+TBVtot$V4[(geneaEEprint$V1*2)] 
  
  geneaRRprint<-geneaRRtot
  geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaRRprint$TBVmatER<-TBVtot$V8[(geneaRRprint$V1*2)-1]+TBVtot$V8[(geneaRRprint$V1*2)] 
  
  geneatotprint<-rbind(geneaEEprint,geneaRRprint)
  
  AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatER",colBy = "Gen")
  resultsEE1[rep,]<-AlphaPart[["TBVmatER"]][["EE"]]
  resultsRR1[rep,]<-AlphaPart[["TBVmatER"]][["RR"]]
  resultsSUM1[rep,]<-AlphaPart[["TBVmatER"]][["Sum"]]
  #plot(AlphaPart,ylim = c(-1,10),color = c(2,3))
  
  geneaEEprint<-geneaEEtot
  geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaEEprint$TBVmatRE<-TBVtot$V3[(geneaEEprint$V1*2)-1]+TBVtot$V3[(geneaEEprint$V1*2)] 
  
  geneaRRprint<-geneaRRtot
  geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaRRprint$TBVmatRE<-TBVtot$V7[(geneaRRprint$V1*2)-1]+TBVtot$V7[(geneaRRprint$V1*2)] 
  
  geneatotprint<-rbind(geneaEEprint,geneaRRprint)
  
  AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatRE",colBy = "Gen")
  resultsEE2[rep,]<-AlphaPart[["TBVmatRE"]][["EE"]]
  resultsRR2[rep,]<-AlphaPart[["TBVmatRE"]][["RR"]]
  resultsSUM2[rep,]<-AlphaPart[["TBVmatRE"]][["Sum"]]
  #plot(AlphaPart,ylim = c(-1,10),color = c(2,3))
  
  prodER[rep]<-TBVER$V3+TBVER$V8
  prodRE[rep]<-TBVRE$V4+TBVRE$V7
}

resumen1<-rbind(data.frame(V1=20:31,V2=apply(resultsEE1,2,mean),
                           V3=apply(resultsEE1,2,sd),V4="E"),
                data.frame(V1=20:31,V2=apply(resultsRR1,2,mean),
                           V3=apply(resultsRR1,2,sd),V4="R"),
                data.frame(V1=20:31,V2=apply(resultsSUM1,2,mean),
                           V3=apply(resultsSUM1,2,sd),V4="SUM"))


resumen2<-rbind(data.frame(V1=20:31,V2=apply(resultsEE2,2,mean),
                           V3=apply(resultsEE2,2,sd),V4="E"),
                data.frame(V1=20:31,V2=apply(resultsRR2,2,mean),
                           V3=apply(resultsRR2,2,sd),V4="R"),
                data.frame(V1=20:31,V2=apply(resultsSUM2,2,mean),
                           V3=apply(resultsSUM2,2,sd),V4="SUM"))

ERplot<-ggplot(resumen1,aes(x=V1, y=V2,fill=V4,color=V4)) +
  geom_line(show.legend = F)+
  geom_ribbon(aes(x=V1,ymin=V2-V3*2,ymax=V2+V3*2),alpha=.15)+
  theme(legend.background = element_rect(fill = "lightblue",colour = 1))+
  lims(y=c(-1,10.5))+
  labs(x="Generacion", y="TBV")

REplot<-ggplot(resumen1,aes(x=V1, y=V2,fill=V4,color=V4)) +
  geom_line(show.legend = F)+
  geom_ribbon(aes(x=V1,ymin=V2-V3*2,ymax=V2+V3*2),alpha=.15)+
  theme(legend.background = element_rect(fill = "lightblue",colour = 1))+
  lims(y=c(-1,10.5))+
  labs(x="Generacion", y="TBV")

(ERplot + REplot + plot_layout(guides="collect"))

prod<-rbind(data.frame(V1="ER",V2=prodER),
            data.frame(V1="RE",V2=prodRE))
prodplot<-ggplot(prod)+
  geom_violin(aes(x=V1,y=V2,fill=V1),show.legend = F)+
  labs(x="Cruce", y="TBV")+
  lims(y=c(-1,6))

prodplot

#Selección por linea macho y linea hembra

resultsEE1<-matrix(ncol=12,nrow=100)
resultsRR1<-matrix(ncol=12,nrow=100)
resultsSUM1<-matrix(ncol=12,nrow=100)

resultsEE2<-matrix(ncol=12,nrow=100)
resultsRR2<-matrix(ncol=12,nrow=100)
resultsSUM2<-matrix(ncol=12,nrow=100)

prodER<-matrix(ncol=1,nrow=100)
prodRE<-matrix(ncol=1,nrow=100)


for (rep in 1:100) {
  
  
  load("pobbase.RData")
  set.seed(Sys.time())
  
  for (Generation in 21:30) {
    
    newgeneaEE<-data.frame(V1=c(1:1000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="EE")
    newgeneaRR<-data.frame(V1=c(1001:2000)+3000*Generation,V2=0,V3=0,V4=c(1:2),Raza="RR")
    newgeneaER<-data.frame(V1=c(2001:2500)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="ER")
    newgeneaRE<-data.frame(V1=c(2501:3000)+3000*Generation,V2=0,V3=0,V4=c(2),Raza="RE")
    
    #escoger los reproductores EE
    machosEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==1]) 
    machosEE$V2<-TBVtot$V3[(machosEE$V1)*2-1]+TBVtot$V3[(machosEE$V1)*2]
    machosEE1<-machosEE$V1[order(machosEE$V2,decreasing = T)]
    machosEE$V2<-TBVtot$V4[(machosEE$V1)*2-1]+TBVtot$V4[(machosEE$V1)*2]
    machosEE2<-machosEE$V1[order(machosEE$V2,decreasing = T)]
    
    hembrasEE<-data.frame(V1=geneaEE$V1[geneaEE$V4==2])
    hembrasEE$V2<-TBVtot$V3[(hembrasEE$V1)*2-1]+TBVtot$V3[(hembrasEE$V1)*2]
    hembrasEE1<-hembrasEE$V1[order(hembrasEE$V2,decreasing = T)]
    hembrasEE$V2<-TBVtot$V4[(hembrasEE$V1)*2-1]+TBVtot$V4[(hembrasEE$V1)*2]
    hembrasEE2<-hembrasEE$V1[order(hembrasEE$V2,decreasing = T)]
    
    newgeneaEE$V2<-c(sample(machosEE1[c(1:100)],nrow(newgeneaEE)/2,replace = T),sample(machosEE2[c(1:100)],nrow(newgeneaEE)/2,replace = T)) #cojo los 100 mejores machos de 500
    newgeneaEE$V3<-c(sample(hembrasEE1[c(1:100)],nrow(newgeneaEE)/2,replace = T),sample(hembrasEE1[c(1:100)],nrow(newgeneaEE)/2,replace = T)) #cojo las 100 mejores hembras de 500
    
    
    #escoger los reproductores RR
    machosRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==1]) 
    machosRR$V2<-TBVtot$V7[(machosRR$V1)*2-1]+TBVtot$V7[(machosRR$V1)*2]
    machosRR1<-machosRR$V1[order(machosRR$V2,decreasing = T)]
    machosRR$V2<-TBVtot$V8[(machosRR$V1)*2-1]+TBVtot$V8[(machosRR$V1)*2]
    machosRR2<-machosRR$V1[order(machosRR$V2,decreasing = T)]
    
    hembrasRR<-data.frame(V1=geneaRR$V1[geneaRR$V4==2])
    hembrasRR$V2<-TBVtot$V7[(hembrasRR$V1)*2-1]+TBVtot$V7[(hembrasRR$V1)*2]
    hembrasRR1<-hembrasRR$V1[order(hembrasRR$V2,decreasing = T)]
    hembrasRR$V2<-TBVtot$V8[(hembrasRR$V1)*2-1]+TBVtot$V8[(hembrasRR$V1)*2]
    hembrasRR2<-hembrasRR$V1[order(hembrasRR$V2,decreasing = T)]
    
    newgeneaRR$V2<-c(sample(machosRR1[c(1:100)],nrow(newgeneaRR)/2,replace = T),sample(machosRR2[c(1:100)],nrow(newgeneaRR)/2,replace = T)) #cojo los 100 mejores machos de 500
    newgeneaRR$V3<-c(sample(hembrasRR1[c(1:100)],nrow(newgeneaRR)/2,replace = T),sample(hembrasRR2[c(1:100)],nrow(newgeneaRR)/2,replace = T)) #cojo las 100 mejores hembras de 500
    
    #animales cruzados
    
    newgeneaER$V2<-sample(newgeneaEE$V1[newgeneaEE$V4==1][c(1:250)],nrow(newgeneaER),replace = T)
    newgeneaER$V3<-sample(newgeneaRR$V1[newgeneaRR$V4==2][c(251:500)],nrow(newgeneaER),replace = T)
    
    newgeneaRE$V3<-sample(newgeneaEE$V1[newgeneaEE$V4==2][c(251:500)],nrow(newgeneaER),replace = T)
    newgeneaRE$V2<-sample(newgeneaRR$V1[newgeneaRR$V4==1][c(1:250)],nrow(newgeneaER),replace = T)
    
    geneaEE<-newgeneaEE
    geneaRR<-newgeneaRR
    geneaER<-newgeneaER
    geneaRE<-newgeneaRE
    
    #genea completa
    
    geneatot<-rbind(geneatot,newgeneaEE,newgeneaRR,newgeneaER,newgeneaRE)
    
    TBVEnew<-as.data.frame(matrix(0,nrow = nrow(geneaEE)*2,ncol = 8))
    filas<-seq_len(nrow(TBVEnew)) %% 2
    TBVEnew[filas==1,]<-TBVtot[(geneaEE$V2*2)-1,]*0.5+TBVtot[geneaEE$V2*2,]*0.5
    TBVEnew[filas==0,]<-TBVtot[(geneaEE$V3*2)-1,]*0.5+TBVtot[geneaEE$V3*2,]*0.5
    TBVEnew<-TBVEnew+(mvrnorm(nrow(geneaEE)*2,medias,varcovartot/2))
    
    TBVRnew<-as.data.frame(matrix(0,nrow = nrow(geneaRR)*2,ncol = 8))
    filas<-seq_len(nrow(TBVRnew)) %% 2
    TBVRnew[filas==1,]<-TBVtot[(geneaRR$V2*2)-1,]*0.5+TBVtot[geneaRR$V2*2,]*0.5
    TBVRnew[filas==0,]<-TBVtot[(geneaRR$V3*2)-1,]*0.5+TBVtot[geneaRR$V3*2,]*0.5
    TBVRnew<-TBVRnew+(mvrnorm(nrow(geneaRR)*2,medias,varcovartot/2))
    
    TBVE<-TBVEnew
    TBVR<-TBVRnew
    
    TBVtot<-rbind(TBVtot,TBVE,TBVR)
    
    TBVER<-as.data.frame(matrix(0,nrow = nrow(geneaER)*2,ncol = 8))
    filas<-seq_len(nrow(TBVER)) %% 2
    TBVER[filas==1,]<-TBVtot[(geneaER$V2*2)-1,]*0.5+TBVtot[geneaER$V2*2,]*0.5
    TBVER[filas==0,]<-TBVtot[(geneaER$V3*2)-1,]*0.5+TBVtot[geneaER$V3*2,]*0.5
    TBVER<-TBVER+(mvrnorm(nrow(geneaER)*2,medias,varcovartot/2))
    
    TBVRE<-as.data.frame(matrix(0,nrow = nrow(geneaRE)*2,ncol = 8))
    filas<-seq_len(nrow(TBVRE)) %% 2
    TBVRE[filas==1,]<-TBVtot[(geneaRE$V2*2)-1,]*0.5+TBVtot[geneaRE$V2*2,]*0.5
    TBVRE[filas==0,]<-TBVtot[(geneaRE$V3*2)-1,]*0.5+TBVtot[geneaRE$V3*2,]*0.5
    TBVRE<-TBVRE+(mvrnorm(nrow(geneaRE)*2,medias,varcovartot/2))
    
    TBVtot<-rbind(TBVtot,TBVER,TBVRE)
    
  }
  
  geneaEEtot<-geneatot[geneatot$Raza=="EE",]
  geneaRRtot<-geneatot[geneatot$Raza=="RR",]
  
  geneaEEprint<-geneaEEtot
  geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaEEprint$TBVmatER<-TBVtot$V4[(geneaEEprint$V1*2)-1]+TBVtot$V4[(geneaEEprint$V1*2)] 
  
  geneaRRprint<-geneaRRtot
  geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaRRprint$TBVmatER<-TBVtot$V8[(geneaRRprint$V1*2)-1]+TBVtot$V8[(geneaRRprint$V1*2)] 
  
  geneatotprint<-rbind(geneaEEprint,geneaRRprint)
  
  AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatER",colBy = "Gen")
  resultsEE1[rep,]<-AlphaPart[["TBVmatER"]][["EE"]]
  resultsRR1[rep,]<-AlphaPart[["TBVmatER"]][["RR"]]
  resultsSUM1[rep,]<-AlphaPart[["TBVmatER"]][["Sum"]]
  #plot(AlphaPart,ylim = c(-1,10),color = c(2,3))
  
  geneaEEprint<-geneaEEtot
  geneaEEprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaEEprint$TBVmatRE<-TBVtot$V3[(geneaEEprint$V1*2)-1]+TBVtot$V3[(geneaEEprint$V1*2)] 
  
  geneaRRprint<-geneaRRtot
  geneaRRprint$Gen<-c(rep(20,times=20000),rep(c(21:31),each=1000))
  geneaRRprint$TBVmatRE<-TBVtot$V7[(geneaRRprint$V1*2)-1]+TBVtot$V7[(geneaRRprint$V1*2)] 
  
  geneatotprint<-rbind(geneaEEprint,geneaRRprint)
  
  AlphaPart<-AlphaPart(geneatotprint,colId = "V1",colFid = "V2",colMid = "V3",colPath = "Raza",colBV = "TBVmatRE",colBy = "Gen")
  resultsEE2[rep,]<-AlphaPart[["TBVmatRE"]][["EE"]]
  resultsRR2[rep,]<-AlphaPart[["TBVmatRE"]][["RR"]]
  resultsSUM2[rep,]<-AlphaPart[["TBVmatRE"]][["Sum"]]
  #plot(AlphaPart,ylim = c(-1,10),color = c(2,3))
  
  prodER[rep]<-TBVER$V3+TBVER$V8
  prodRE[rep]<-TBVRE$V4+TBVRE$V7
}

resumen1<-rbind(data.frame(V1=20:31,V2=apply(resultsEE1,2,mean),
                           V3=apply(resultsEE1,2,sd),V4="E"),
                data.frame(V1=20:31,V2=apply(resultsRR1,2,mean),
                           V3=apply(resultsRR1,2,sd),V4="R"),
                data.frame(V1=20:31,V2=apply(resultsSUM1,2,mean),
                           V3=apply(resultsSUM1,2,sd),V4="SUM"))


resumen2<-rbind(data.frame(V1=20:31,V2=apply(resultsEE2,2,mean),
                           V3=apply(resultsEE2,2,sd),V4="E"),
                data.frame(V1=20:31,V2=apply(resultsRR2,2,mean),
                           V3=apply(resultsRR2,2,sd),V4="R"),
                data.frame(V1=20:31,V2=apply(resultsSUM2,2,mean),
                           V3=apply(resultsSUM2,2,sd),V4="SUM"))

ERplot<-ggplot(resumen1,aes(x=V1, y=V2,fill=V4,color=V4)) +
  geom_line(show.legend = F)+
  geom_ribbon(aes(x=V1,ymin=V2-V3*2,ymax=V2+V3*2),alpha=.15)+
  theme(legend.background = element_rect(fill = "lightblue",colour = 1))+
  lims(y=c(-1,10.5))+
  labs(x="Generacion", y="TBV")

REplot<-ggplot(resumen1,aes(x=V1, y=V2,fill=V4,color=V4)) +
  geom_line(show.legend = F)+
  geom_ribbon(aes(x=V1,ymin=V2-V3*2,ymax=V2+V3*2),alpha=.15)+
  theme(legend.background = element_rect(fill = "lightblue",colour = 1))+
  lims(y=c(-1,10.5))+
  labs(x="Generacion", y="TBV")

(ERplot + REplot + plot_layout(guides="collect"))

prod<-rbind(data.frame(V1="ER",V2=prodER),
            data.frame(V1="RE",V2=prodRE))
prodplot<-ggplot(prod)+
  geom_violin(aes(x=V1,y=V2,fill=V1),show.legend = F)+
  labs(x="Cruce", y="TBV")+
  lims(y=c(-1,6))

prodplot

#NEW VARCOVARTOT

newvarcovartot<-matrix(0,8,8)
newa<-c(1,5,3,7,2,6,8,4)
for (a in 1:8) {
  for (b in 1:8) {
   newvarcovartot[a,b]<-varcovartot[newa[a],newa[b]] 
  }
  
  
}
round(newvarcovartot,4)
