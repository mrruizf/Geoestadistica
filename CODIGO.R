############################INGENIERIA CATSTRAL Y GEODESIA############################
####################################GEOESTADISTICA####################################

library(gstat)
library(geoR)
library(sgeostat)
library(geospt)
library(scatterplot3d)
library(ggplot2)
library(car)
library(xtable)
library(stargazer)
library(sp)
library(spTest)
library(nortest)
library(intamap)
library(spatial)
library(ggmap)
library(maptools)
library(leaflet)
library(sf)
library(rgdal)
library(dplyr)
library(RGeostats)
library(maps)
library(mapview)


      #####################################
      ##CARGAR BASE DE DATOS PARA CROACIA##
      #####################################

load("/BBDD/Croacia/IDSTA.ov.rda") #leer la base de datos
muestra<-as.data.frame(IDSTA.ov)
datos=muestra[c("HRdem","HRdsea","HRtwi","Lon","Lat","LST2008_09_29")]

      ##########################
      ##ESPACIALIZAR LOS DATOS##
      ##########################

coords <- cbind(datos$Lon,datos$Lat)
SP = SpatialPoints(coords)
SPDF2 <- SpatialPointsDataFrame(coords, datos)
##DAR UN SISTEMA DE REFERENCIA GEOGRAFICO
proj4string(SPDF2)
proj4string(SPDF2)<-CRS("+proj=longlat +datum=WGS84")
##REPROYECTAR A COORDENADAS PLANAS
utm33 <- "+proj=utm +zone=33 +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
zona_estudio <- spTransform(SPDF2,CRS(utm33))
summary(zona_estudio)
class(zona_estudio)

      ###################################################################
      ##DATOS EXPORTADOS EN ARCHIVO SHAPEFILE LIBRERIA RGDAL (opcional)##
      ###################################################################
writeOGR(zona_estudio, "/BBDD/Croacia/croacia.shp",layer="temp", driver="ESRI Shapefile")
shape1<-readOGR("/BBDD/Croacia/croacia.shp")
#VISUALIZACI?N DEL SHAPE
plot(shape1, border="BROWN",axes=T)
proj4string(shape1)

###################################################
####################CONTINUANDO###################
###################################################
croacia<-data.frame(zona_estudio)
names(croacia)<-c("HRdem","HRdsea","HRtwi","lon","lat","temp","E","N")
#estadisticas de la temperatura
(summary(croacia$temp))

#eliminar NA de la variable temperatura
croacia1<-na.omit(croacia,cols=temp)
summary(croacia1$temp)

      #########################################
      ##Mapa Area de estudio libreria leaflet##
      #########################################
pal <- colorNumeric(
  palette = "YlOrRd",
  domain = croacia1$temp
)

leaflet(data=croacia1) %>% addTiles() %>%addProviderTiles(providers$Esri.NatGeoWorldMap)%>%
addCircleMarkers(lng = croacia1$lon, lat = croacia1$lat,radius = 3,color = ~pal(temp), stroke = FALSE, fillOpacity = 1)%>%
addLegend("bottomright", pal = pal ,values = ~temp,title = "Temperatura ?C",opacity = 1) 

      #############################
      ##Estadisticas descriptivas##
      #############################

xtable(summary(croacia1))
##Anal?sis de normalidad
##m?todos graficos
qplot(croacia1$temp, geom = 'blank') +   
  geom_line(aes(y = ..density.., colour = 'Empirica'), stat = 'density') +  
  stat_function(fun = dnorm, aes(colour = 'Normal')) +                       
  geom_histogram(aes(y = ..density..), alpha = 0.4) +                        
  scale_colour_manual(name = 'Densidad', values = c('red', 'blue')) + 
  theme(legend.position = c(0.8, 0.8)) 
##m?todos Estad?sticos
shapiro.test(croacia1$temp)
ks.test(as.numeric(scale(sort(croacia1$temp))),pnorm)

      ##################################################
      ##Anamorfosis Gaussiana de la libreria RGeostats##
      ##################################################

croacia.rgdb <- db.create(croacia1,ndim=2,autoname=F)
croacia.herm <- anam.fit(croacia.rgdb,name="temp",type="gaus")
croacia.hermtrans <- anam.z2y(croacia.rgdb,names="temp",anam=croacia.herm)
temp.trans <- croacia.hermtrans@items$Gaussian.temp
##m?todos graficos para la variable transfomanda
qplot(temp.trans, geom = 'blank') +   
  geom_line(aes(y = ..density.., colour = 'Empirica'), stat = 'density') +  
  stat_function(fun = dnorm, aes(colour = 'Normal')) +                       
  geom_histogram(aes(y = ..density..), alpha = 0.4) +                        
  scale_colour_manual(name = 'Densidad', values = c('red', 'blue')) + 
  theme(legend.position = c(0.8, 0.8)) 
##m?todos Estad?sticos para la variable transfomanda
shapiro.test(temp.trans)
ks.test(as.numeric(scale(sort(temp.trans))),pnorm)

      ####################################################
      ##An?lisis de tendencia y elecci?n del mejor modelo##
      ####################################################

modelo_lineal<-lm(temp.trans~HRdem+HRdsea+HRtwi+N+E,data=croacia1)
summary(modelo_lineal)
MejorModelo<-stepAIC(modelo_lineal,direction = "both")
summary(MejorModelo)
#Test de Normalidad a los residuos
residuales<-residuals(MejorModelo)
##graficos
qplot(residuales, geom = 'blank') +   
  geom_line(aes(y = ..density.., colour = 'Empirica'), stat = 'density') +  
  stat_function(fun = dnorm, aes(colour = 'Normal')) +                       
  geom_histogram(aes(y = ..density..), alpha = 0.4) +                        
  scale_colour_manual(name = 'Densidad', values = c('red', 'blue')) + 
  theme(legend.position = c(0.8, 0.8)) 
##estadisticos
shapiro.test(residuales)
ks.test(as.numeric(scale(sort(residuales))),pnorm)

      ###########################
      ##An?lisis de anisotrop?a##
      ###########################

##primero se genera la geodata para realizar los variogramas
croacia2<-data.frame(croacia1,temp.trans)
geodata.croacia<-as.geodata(croacia2,coords.col = 7:8,data.col =10)
Dmax<-sqrt((max(croacia2[, 7])-min(croacia2[, 7]))^2+(max(croacia2[, 8])-min(croacia2[, 8]))^2)
Dusar<-Dmax/2
x11()
plot(variog4(geodata.croacia,max.dist = Dusar), xlab = "Distancia (m)", ylab = "Semivariograma estimado", legend = F)
legend(0,2.4, legend = c(expression(0*degree), expression(45*degree), expression(90*degree), expression(135*degree)), col = 1:4, lty = 1:4)
title("Semivariogramas experimentales (estimador clasico)")
##m?todos "estadisticos" para la anisotropia
estimateAnisotropy(IDSTA.ov,"temp")
##source("/BBDD/rose.R")  -->NO ENTENDI COMO USARLA
##rose(variog4(geodata.croacia),300000)

        #######################################################
        ##Modelamiento del semivariograma para Kriging Simple##
        #######################################################
mu<-mean(croacia2$temp.trans)
error<-croacia2$temp.trans-mu
croacia2<-data.frame(croacia2,error)
Point <- point(croacia2, x = "E", y = "N") 
Pair <- pair(Point,num.lags=25,maxdist=Dusar)
EstimaVariogramaks<- est.variograms(Point,Pair,"error",trim=0.1) 

##Ajuste de los Semivariogramas teoricos a los experimentales
##CLASICO
Exp.ml.ks<-likfit(geodata = geodata.croacia,nugget = 0.05,ini = c(1,60000),fix.nug=T) 
Sph.ml.ks<-likfit(geodata = geodata.croacia,nugget = 0.01, ini = c(2,65000),cov.model="sph",fix.nug=T) 
Matern.ml.ks<-likfit(geodata = geodata.croacia,nugget = 0.2, ini = c(2,67000),cov.model="mat",kappa=0.66,fix.kappa = T,fix.nug=T) 
Cir.ml.ks<-likfit(geodata = geodata.croacia,nugget = 0.091, ini = c(1.9,72000),cov.model="cir",fix.nug=T) 
Pow.ml.ks<-likfit(geodata = geodata.croacia, nugget = 0.1,ini = c(2.2,52000),cov.model="powered.exponential",kappa=1.15,fix.nug=T) 
x11()
plot(EstimaVariogramaks$bins,EstimaVariogramaks$classic,lty=1,ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental Clasico",xlab="Distancia (m)", ylab="Semivarianza") 
lines(Exp.ml.ks,max.dist=Dusar,lwd=3,col="lawngreen") 
lines(Sph.ml.ks,max.dist=Dusar,lwd=3,col="tan1") 
lines(Matern.ml.ks,max.dist=Dusar,lwd=3,col="tomato") 
lines(Cir.ml.ks,max.dist=Dusar,lwd=3,col="mediumorchid") 
lines(Pow.ml.ks,max.dist=Dusar,lwd=3,col="yellow") 
legend(locator(1),c('Exponencial','Esferico','Matern','Cir','Pow'),col=c("lawngreen","tan1","tomato","mediumorchid", "yellow"), lty=c(1,1,1,1,1,1,1))

##ROBUSTO
Exp.mlr.ks<-likfit(geodata = geodata.croacia,nugget = 0.02,ini = c(1.5,55000),fix.nug=T) 
Sph.mlr.ks<-likfit(geodata = geodata.croacia,nugget = 0.01, ini = c(2,70000),cov.model="sph",fix.nug=T) 
Matern.mlr.ks<-likfit(geodata = geodata.croacia,nugget = 0.31, ini = c(1.4,58000),cov.model="mat",kappa=2,fix.nug=T) 
Cir.mlr.ks<-likfit(geodata = geodata.croacia,nugget = 0.01, ini = c(1,58000),cov.model="cir",fix.nug=T) 
Pow.mlr.ks<-likfit(geodata = geodata.croacia, nugget = 0.1,ini = c(0.7,30000),cov.model="powered.exponential",kappa=1.75,fix.nug=T) 
x11()
plot(EstimaVariogramaks$bins,EstimaVariogramaks$robust,lty=1,ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental Robusto",xlab="Distancia (m)", ylab="Semivarianza") 
lines(Exp.mlr.ks,max.dist=Dusar,lwd=3,col="lawngreen") 
lines(Sph.mlr.ks,max.dist=Dusar,lwd=3,col="tan1") 
lines(Matern.mlr.ks,max.dist=Dusar,lwd=3,col="tomato") 
lines(Cir.mlr.ks,max.dist=Dusar,lwd=3,col="mediumorchid") 
lines(Pow.mlr.ks,max.dist=Dusar,lwd=3,col="yellow") 
legend(locator(1),c('Exponencial','Esferico','Matern','Cir','Pow'),col=c("lawngreen","tan1","tomato","mediumorchid", "yellow"), lty=c(1,1,1,1,1,1,1))

##MEDIANA
Exp.mlm.ks<-likfit(geodata = geodata.croacia,nugget = 0.12,ini = c(2.5,90000),fix.nug=T) 
Sph.mlm.ks<-likfit(geodata = geodata.croacia,nugget = 0.08, ini = c(1.3,90000),cov.model="sph",fix.nug=T) 
Matern.mlm.ks<-likfit(geodata = geodata.croacia,nugget = 0.18, ini = c(1.9,47000),cov.model="mat",kappa=0.8,fix.kappa = T,fix.nug=T)
Cir.mlm.ks<-likfit(geodata = geodata.croacia,nugget = 0.01, ini = c(1,60000),cov.model="cir",fix.nug=T) 
Pow.mlm.ks<-likfit(geodata = geodata.croacia, nugget = 0.1,ini = c(0.7,28200),cov.model="powered.exponential",kappa=1.75,fix.nug=T) 
x11()
plot(EstimaVariogramaks$bins,EstimaVariogramaks$med,lty=1,ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental a partir de la Mediana",xlab="Distancia (m)", ylab="Semivarianza") 
lines(Exp.mlm.ks,max.dist=Dusar,lwd=3,col="lawngreen") 
lines(Sph.mlm.ks,max.dist=Dusar,lwd=3,col="tan1") 
lines(Matern.mlm.ks,max.dist=Dusar,lwd=3,col="tomato") 
lines(Cir.mlm.ks,max.dist=Dusar,lwd=3,col="mediumorchid") 
lines(Pow.mlm.ks,max.dist=Dusar,lwd=3,col="yellow") 
legend(locator(1),c('Exponencial','Esferico','Matern','Cir','Pow'),col=c("lawngreen","tan1","tomato","mediumorchid", "yellow"), lty=c(1,1,1,1,1,1,1))

##MEDIA RECORTADA
Exp.mlmr.ks<-likfit(geodata = geodata.croacia,nugget = 0.09,ini = c(1.8,60000),fix.nug=T) 
Sph.mlmr.ks<-likfit(geodata = geodata.croacia,nugget = 0.04, ini = c(1.6,80000),cov.model="sph",fix.nug=T) 
Matern.mlmr.ks<-likfit(geodata = geodata.croacia,nugget = 0.17, ini = c(1.9,45000),cov.model="mat",kappa=0.8,fix.kappa = T,fix.nug=T)
Cir.mlmr.ks<-likfit(geodata = geodata.croacia,nugget = 0.1, ini = c(1,75000),cov.model="cir",fix.nug=T) 
Pow.mlmr.ks<-likfit(geodata = geodata.croacia, nugget = 0.1,ini = c(0.7,28200),cov.model="powered.exponential",kappa=1.72,fix.nug=T) 
x11()
plot(EstimaVariogramaks$bins,EstimaVariogramaks$trimmed.mean,lty=1,ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental media recortada",xlab="Distancia (m)", ylab="Semivarianza") 
lines(Exp.mlmr.ks,max.dist=Dusar,lwd=3,col="lawngreen") 
lines(Sph.mlmr.ks,max.dist=Dusar,lwd=3,col="tan1") 
lines(Matern.mlmr.ks,max.dist=Dusar,lwd=3,col="tomato") 
lines(Cir.mlmr.ks,max.dist=Dusar,lwd=3,col="mediumorchid") 
lines(Pow.mlmr.ks,max.dist=Dusar,lwd=3,col="yellow") 
legend(locator(1),c('Exponencial','Esferico','Matern','Cir','Pow'),col=c("lawngreen","tan1","tomato","mediumorchid", "yellow"), lty=c(1,1,1,1,1,1,1))

##AIC modelos teoricos estimaci?n clasica.
s2<-matrix(c(Exp.ml.ks$AIC,Sph.ml.ks$AIC,Matern.ml.ks$AIC,Cir.ml.ks$AIC,Pow.ml.ks$AIC))
rownames(s2)<-c("Modelo Exponencial","Modelo esf?rico","Modelo Matern","Modelo circular","Modelo de potencia")
colnames(s2)<-c("AIC")
s2
xtable(t(s2))

#Clasico, OLS WLS RML ML

geodata.croacia1<-as.geodata(croacia2,coords.col = 7:8,data.col =11)

variograma.clasico.ks<-variog(geodata.croacia1,trend="cte",max.dist=Dusar, option = "cloud",estimator.type="modulus")
matern.ks.wls<-variofit(vario = variograma.clasico.ks, nugget = 0.2,ini = c(2.110284, 67000),kappa = 0.69,fix.nugget= T,weights="npairs",cov.model = "matern")
matern.ks.ml<-likfit(geodata = geodata.croacia, nugget = 0.2,ini = c(2.110284, 67000),kappa = 0.68,fix.nugget= T)
matern.ks.rml<-likfit(geodata = geodata.croacia,nugget = 0.2,ini = c(2.110284, 67000),kappa = 0.64,fix.nugget= T,method='RML')

x11()##Eleccion del m?todo de estimaci?n por medios gr?ficos
plot(EstimaVariogramaks$bins,EstimaVariogramaks$classic,lty=1,ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental Cl?sico",xlab="Distancia (m)", ylab="Semivarianza") 
lines(matern.ks.ml,max.dist=Dusar,lwd=2,col="navyblue")
lines(matern.ks.rml,max.dist=Dusar,lwd=3,lty=1, col="yellow")
lines(matern.ks.wls,max.dist=Dusar,lwd=2,lty=1,col="coral")
legend(locator(1),c('ML','RML','WLS'),lty=c(1,1,1),col=c("navyblue","yellow","coral"))

      ##########################################################
      ##Modelamiento del semivariograma para Kriging ordinario##
      ##########################################################

EstimaVariograma<- est.variograms(Point,Pair,"temp.trans",trim=0.1) 
##Estimador Clasico
plot(EstimaVariograma$bins,EstimaVariograma$classic,lty=1,ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental Clasico",xlab="Distancia (m)", ylab="Semivarianza") 
lines(EstimaVariograma$bins,EstimaVariograma$classic, col=2)
##Estimador Robusto, pepita 1.0, rango 125000, meseta 2.0 aprox
plot(EstimaVariograma$bins,EstimaVariograma$robust,lty=1,ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental Robusto",xlab="Distancia (m)", ylab="Semivarianza") 
lines(EstimaVariograma$bins,EstimaVariograma$robust, col=2)
##Estimador Mediana, pepita 0.9, rango 125000, meseta 2.5 aprox
plot(EstimaVariograma$bins,EstimaVariograma$med,lty=1,ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental Mediana",xlab="Distancia (m)", ylab="Semivarianza") 
lines(EstimaVariograma$bins,EstimaVariograma$med, col=2)
##Estimador Media recortada, pepita 0.8, rango 150000, meseta 2.0 aprox
plot(EstimaVariograma$bins,EstimaVariograma$trimmed.mean,lty=1,ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental Media recortada",xlab="Distancia (m)", ylab="Semivarianza") 
lines(EstimaVariograma$bins,EstimaVariograma$trimmed.mean, col=2)

    ##Ajuste de los Semivariogramas teoricos a los experimentales
##CLASICO
Exp.ml<-likfit(geodata = geodata.croacia,nugget = 0.05,ini = c(2,55000),fix.nug=T) 
Sph.ml<-likfit(geodata = geodata.croacia,nugget = 0.01, ini = c(1.4,65000),cov.model="sph",fix.nug=T) 
Matern.ml<-likfit(geodata = geodata.croacia,nugget = 0.3, ini = c(2,67000),cov.model="mat",kappa=1.5,fix.nug=T) 
Cir.ml<-likfit(geodata = geodata.croacia,nugget = 0.01, ini = c(2.2,52000),cov.model="cir",fix.nug=T) 
Pow.ml<-likfit(geodata = geodata.croacia, nugget = 0.1,ini = c(2,27500),cov.model="powered.exponential",kappa=1.75,fix.nug=T) 
x11()
  plot(EstimaVariograma$bins,EstimaVariograma$classic,lty=1,ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental Clasico",xlab="Distancia (m)", ylab="Semivarianza") 
lines(Exp.ml,max.dist=Dusar,lwd=3,col="lawngreen") 
lines(Sph.ml,max.dist=Dusar,lwd=3,col="tan1") 
lines(Matern.ml,max.dist=Dusar,lwd=3,col="tomato") 
lines(Cir.ml,max.dist=Dusar,lwd=3,col="mediumorchid") 
lines(Pow.ml,max.dist=Dusar,lwd=3,col="yellow") 
legend(locator(1),c('Exponencial','Esferico','Matern','Cir','Pow'),col=c("lawngreen","tan1","tomato","mediumorchid", "yellow"), lty=c(1,1,1,1,1,1,1))

##ROBUSTO
Exp.mlr<-likfit(geodata = geodata.croacia,nugget = 0.05,ini = c(1,60000),fix.nug=T) 
Sph.mlr<-likfit(geodata = geodata.croacia,nugget = 0.01, ini = c(0.9,70000),cov.model="sph",fix.nug=T) 
Matern.mlr<-likfit(geodata = geodata.croacia,nugget = 0.3, ini = c(1,66000),cov.model="mat",kappa=1.5,fix.nug=T) 
Cir.mlr<-likfit(geodata = geodata.croacia,nugget = 0.01, ini = c(1,58000),cov.model="cir",fix.nug=T) 
Pow.mlr<-likfit(geodata = geodata.croacia, nugget = 0.1,ini = c(0.7,30000),cov.model="powered.exponential",kappa=1.75,fix.nug=T) 
x11()
plot(EstimaVariograma$bins,EstimaVariograma$robust,lty=1,ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental Robusto",xlab="Distancia (m)", ylab="Semivarianza") 
lines(Exp.mlr,max.dist=Dusar,lwd=3,col="lawngreen") 
lines(Sph.mlr,max.dist=Dusar,lwd=3,col="tan1") 
lines(Matern.mlr,max.dist=Dusar,lwd=3,col="tomato") 
lines(Cir.mlr,max.dist=Dusar,lwd=3,col="mediumorchid") 
lines(Pow.mlr,max.dist=Dusar,lwd=3,col="yellow") 
legend(locator(1),c('Exponencial','Esferico','Matern','Cir','Pow'),col=c("lawngreen","tan1","tomato","mediumorchid", "yellow"), lty=c(1,1,1,1,1,1,1))

##MEDIANA
Exp.mlm<-likfit(geodata = geodata.croacia,nugget = 0.05,ini = c(0.7,55000),fix.nug=T) 
Sph.mlm<-likfit(geodata = geodata.croacia,nugget = 0.01, ini = c(0.9,72000),cov.model="sph",fix.nug=T) 
Matern.mlm<-likfit(geodata = geodata.croacia,nugget = 0.3, ini = c(0.5,62000),cov.model="mat",kappa=1.5,fix.nug=T) 
Cir.mlm<-likfit(geodata = geodata.croacia,nugget = 0.01, ini = c(1,60000),cov.model="cir",fix.nug=T) 
Pow.mlm<-likfit(geodata = geodata.croacia, nugget = 0.1,ini = c(0.7,28200),cov.model="powered.exponential",kappa=1.75,fix.nug=T) 
x11()
plot(EstimaVariograma$bins,EstimaVariograma$med,lty=1,ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental a partir de la Mediana",xlab="Distancia (m)", ylab="Semivarianza") 
lines(Exp.mlm,max.dist=Dusar,lwd=3,col="lawngreen") 
lines(Sph.mlm,max.dist=Dusar,lwd=3,col="tan1") 
lines(Matern.mlm,max.dist=Dusar,lwd=3,col="tomato") 
lines(Cir.mlm,max.dist=Dusar,lwd=3,col="mediumorchid") 
lines(Pow.mlm,max.dist=Dusar,lwd=3,col="yellow") 
legend(locator(1),c('Exponencial','Esferico','Matern','Cir','Pow'),col=c("lawngreen","tan1","tomato","mediumorchid", "yellow"), lty=c(1,1,1,1,1,1,1))

##MEDIA RECORTADA
Exp.mlmr<-likfit(geodata = geodata.croacia,nugget = 0.02,ini = c(1,50000),fix.nug=T) 
Sph.mlmr<-likfit(geodata = geodata.croacia,nugget = 0.02, ini = c(1.5,70000),cov.model="sph",fix.nug=T) 
Matern.mlmr<-likfit(geodata = geodata.croacia,nugget = 0.28, ini = c(0.5,62000),cov.model="mat",kappa=1.5,fix.nug=T) 
Cir.mlmr<-likfit(geodata = geodata.croacia,nugget = 0.1, ini = c(1,75000),cov.model="cir",fix.nug=T) 
Pow.mlmr<-likfit(geodata = geodata.croacia, nugget = 0.11,ini = c(0.7,28200),cov.model="powered.exponential",kappa=1.75,fix.nug=T) 
x11()
plot(EstimaVariograma$bins,EstimaVariograma$trimmed.mean,lty=1,ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental media recortada",xlab="Distancia (m)", ylab="Semivarianza") 
lines(Exp.mlmr,max.dist=Dusar,lwd=3,col="lawngreen") 
lines(Sph.mlmr,max.dist=Dusar,lwd=3,col="tan1") 
lines(Matern.mlmr,max.dist=Dusar,lwd=3,col="tomato") 
lines(Cir.mlmr,max.dist=Dusar,lwd=3,col="mediumorchid") 
lines(Pow.mlmr,max.dist=Dusar,lwd=3,col="yellow") 
legend(locator(1),c('Exponencial','Esferico','Matern','Cir','Pow'),col=c("lawngreen","tan1","tomato","mediumorchid", "yellow"), lty=c(1,1,1,1,1,1,1))

##AIC modelos teoricos estimaci?n robusta.
s1<-matrix(c(Exp.mlr$AIC,Sph.mlr$AIC,Matern.mlr$AIC,Cir.mlr$AIC,Pow.mlr$AIC))
rownames(s1)<-c("Modelo Exponencial","Modelo esf?rico","Modelo Matern","Modelo circular","Modelo de potencia")
colnames(s1)<-c("AIC")
s1
xtable(t(s1))

#Robusto, OLS WLS RML ML

variograma.robusto<-variog(geodata.croacia,trend="cte",max.dist=Dusar, option = "cloud",estimator.type="modulus")
exp.wls<-variofit(vario = variograma.robusto, nugget = 0.93,ini = c(1.5, 60000),fix.nugget= T,weights="npairs",cov.model = "exponential")
exp.ml<-likfit(geodata = geodata.croacia, nugget = 0.05,ini = c(1.5, 51000),fix.nugget= T)
exp.rml<-likfit(geodata = geodata.croacia,nugget = 0.1,ini = c(1.6,65000),fix.nugget= T,method='RML')

x11()##Eleccion del m?todo de estimaci?n por medios gr?ficos
plot(EstimaVariograma$bins,EstimaVariograma$robust,lty=1,ylim=c(0,3), col =1,pch= 16, main = "Modelo experimental Robusto",xlab="Distancia (m)", ylab="Semivarianza") 
lines(exp.ml,max.dist=Dusar,lwd=2,col="navyblue")
lines(exp.rml,max.dist=Dusar,lwd=3,lty=1, col="yellow")
lines(exp.wls,max.dist=Dusar,lwd=2,lty=1,col="coral")
legend(locator(1),c('ML','RML','WLS'),lty=c(1,1,1),col=c("navyblue","yellow","coral"))
##residuales
#resi.Exp.ML <- resid(exp.ml, spatial = FALSE)
#resi.Exp.RML <- resid(exp.rml, spatial = FALSE)
#resi.Exp.WLS <- sum((EstimaVariograma$robust-exp.wls$gamma)^2)/7
#print(data.frame(resi.Exp.ML,resi.Exp.RML,resi.Exp.WLS))

      ##################################################
      #####M?todos de interpolaci?n geoestad?sticos#####
      ##################################################

##SE CARGA EL SHAPE DE CROACIA
shape<-readOGR("/BBDD/Croacia/CroatiaL.shp")
mapView(shape)##cambiar a open street map para visualizar el shape
#TRANSFORMACION A COORDS. UTM 33N
proj4string(shape)
shape <- spTransform(shape,CRS(utm33))
# Rejilla discreta de predicci?n//se utilizan valores max y min de las estaciones con incertidumbre
xx <- seq((min(croacia2[, 7])-15000), (max(croacia2[, 7])+15000), l = 95)##455.000m
yy <- seq((min(croacia2[, 8])-15000), (max(croacia2[, 8])+15000), l = 95)##432.000m
pred.grid <- expand.grid(x = xx, y = yy) 
### VISUALIZACION DE LA GRILLA DE COORDS., LAS ESTACIONES METEREOLOGICAS Y CROACIA
x11()
plot(shape, border="BROWN",axes=T)
points(geodata.croacia$coords, pch = 20)
points(pred.grid, pch = 3, cex = 0.2)
locations.inside(pred.grid, shape)
############################################################################
##kriging ordinario para la variable transformada con un vecindario global##
############################################################################
ko.global<- krige.conv(geodata.croacia, loc = pred.grid, krige = krige.control(obj.m = exp.wls))
coords <- cbind(pred.grid$x,pred.grid$y)
SP1 = SpatialPoints(coords)
SPDF3 <- SpatialPointsDataFrame(coords, pred.grid)
################################################
## antitransformacion usando RGeostats
temp.bts = cbind(coordinates(SPDF3),ko.global$predict)
temp.bts.db = db.create(temp.bts,autoname = T)
tempdb = anam.y2z(temp.bts.db,names="V3",anam = croacia.herm)
#Prediction map
ko.global$predict<- tempdb@items$Raw.V3
###########################################################################
##kriging ordinario para la variable transformada con un vecindario local##
###########################################################################
args(ksline)
ko.local<- ksline(geodata.croacia ,locations = pred.grid,cov.model ="exp", cov.pars = c(1.5,60000),nugget =0.93 ,m0 = "ok",nwin=30 )
#antitransformaci?n
temp.bts1 = cbind(coordinates(SPDF3),ko.local$predict)
temp.bts.db1 = db.create(temp.bts1,autoname = T)
tempdb1 = anam.y2z(temp.bts.db1,names="V3",anam = croacia.herm)
#Prediction map
ko.local$predict<- tempdb1@items$Raw.V3
#############################################
##mapa de predicciones sobre la grilla
x11()
oldpar <- par(mfrow = c(1, 2))
image(ko.global,col = heat.colors(12)) #superficie de predicci?n
map(shape, add = T, col = "sienna4",lwd =3)
title("Predicciones Kriging ordinario con vecindario global")
points(geodata.croacia$coords, pch=20) #a?adir posiciones datos
contour(ko.global,add=T, nlev=10) #a?adir gr?fico de contorno
legend(locator(1), c("14","15","16","17","18","19"), bty = "o", title = "Temperatura ?C", fill = heat.colors(12), xpd = NA)

image(ko.local,col = heat.colors(12))
map(shape, add = T, col = "sienna4",lwd =3)
title("Predicciones Kriging ordinario con vecindario local")
points(geodata.croacia$coords, pch=20) #a?adir posiciones datos
contour(ko.local,add=T, nlev=10) #a?adir gr?fico de contorno
legend(locator(1), c("14","15","16","17","18","19"), bty = "o", title = "Temperatura ?C", fill = heat.colors(12), xpd = NA)

X11()
oldpar <- par(mfrow = c(1, 2))
image(ko.global, val = ko.global$krige.var,col = heat.colors(12)) #superficie de varianzas
map(shape, add = T, col = "springgreen4",lwd =3)
title("Superficie de varianzas (ko) con vecindario global")
points(geodata.croacia$coords, pch=20)
contour(ko.global,val=sqrt(ko.global$krige.var),add=T, nlev=10)
legend(locator(1), c("1","1.1","1.2","1.3","1.4","1.5"), bty = "o", title = "Desviaci?n estandar", fill = heat.colors(12), xpd = NA)

image(ko.local, val = ko.local$krige.var,col = heat.colors(12)) #superficie de varianzas
map(shape, add = T, col = "springgreen4",lwd =3)
title("Superficie de varianzas (ko) con vecindario local")
points(geodata.croacia$coords, pch=20)
contour(ko.local,val=sqrt(ko.local$krige.var),add=T, nlev=10)
legend(locator(1), c("1","1.1","1.2","1.3","1.4","1.5"), bty = "o", title = "Desviaci?n estandar", fill = heat.colors(12), xpd = NA)

##otras opciones para graficar
contour(ko.wls,filled = TRUE)
library(plot3D)
persp3D(xx, yy, matrix(ko.wls$predict, nrow = length(xx)), theta=-60, phi=40)
library(npsp)
spersp(xx, yy, ko.wls$predict, theta=-60, phi=40)

############################################################
###kriging simple para los residuos con vecindario global###
############################################################
ks.global<- krige.conv(geodata.croacia1, loc = pred.grid, krige = krige.control(type.krige = "sk",obj.m = matern.ks.rml,beta = mu))
## antitransformacion usando RGeostats
temp.bts.ks = cbind(coordinates(SPDF3),ks.global$predict)
temp.bts.db.ks = db.create(temp.bts.ks,autoname = T)
tempdb.ks = anam.y2z(temp.bts.db.ks,names="V3",anam = croacia.herm)
#Prediction map
ks.global$predict<- tempdb.ks@items$Raw.V3
###########################################################
###kriging simple para los residuos con vecindario local###
###########################################################
ks.local<- ksline(geodata.croacia ,locations = pred.grid,cov.model ="matern", cov.pars = c(2.062649,66999.999955),nugget =0.2,kappa = 0.64 ,m0 = "av",nwin=30 )
#antitransformaci?n
temp.bts.ks1 = cbind(coordinates(SPDF3),ks.local$predict)
temp.bts.db.ks1 = db.create(temp.bts.ks1,autoname = T)
tempdb.ks1 = anam.y2z(temp.bts.db.ks1,names="V3",anam = croacia.herm)
#Prediction map
ks.local$predict<- tempdb.ks1@items$Raw.V3
############################################
####mapa de predicciones sobre la grilla####
############################################
x11()
oldpar <- par(mfrow = c(1, 2))
image(ks.global,col = heat.colors(12)) #superficie de predicci?n
map(shape, add = T, col = "sienna4",lwd =3)
title("Predicciones Kriging simple con vencindario global")
points(geodata.croacia$coords, pch=20) #a?adir posiciones datos
contour(ks.global,add=T, nlev=5) #a?adir gr?fico de contorno
legend(locator(1), c("14","15","16","17","18","19"), bty = "o", title = "Temperatura ?C", fill = heat.colors(12), xpd = NA)

image(ks.local,col = heat.colors(12))
map(shape, add = T, col = "sienna4",lwd =3)
title("Predicciones Kriging simple con vecindario local")
points(geodata.croacia$coords, pch=20) #a?adir posiciones datos
contour(ks.local,add=T, nlev=5) #a?adir gr?fico de contorno
legend(locator(1), c("14","15","16","17","18","19"), bty = "o", title = "Temperatura ?C", fill = heat.colors(12), xpd = NA)

x11()
oldpar <- par(mfrow = c(1, 2))
image(ks.global, val = ks.global$krige.var) #superficie de varianzas
map(shape, add = T, col = "springgreen4",lwd =3)
title("Superficie de varianzas (ks) con vecindario global")
points(geodata.croacia$coords, pch=20)
contour(ks.global,val=sqrt(ks.global$krige.var),add=T, nlev=8)
legend(locator(1), c("0.6","0.8","1","1.2","1.4","1.6"), bty = "o", title = "Desviaci?n estandar", fill = heat.colors(12), xpd = NA)


image(ko.local, val = ko.local$krige.var,col = heat.colors(12)) #superficie de varianzas
map(shape, add = T, col = "springgreen4",lwd =3)
title("Superficie de varianzas (ks) con vecindario local")
points(geodata.croacia$coords, pch=20)
contour(ko.local,val=sqrt(ko.local$krige.var),add=T, nlev=8)
legend(locator(1), c("1","1.1","1.2","1.3","1.4","1.5"), bty = "o", title = "Desviaci?n estandar", fill = heat.colors(12), xpd = NA)


##otras opciones para graficar
contour(ks.local,filled = TRUE)
library(plot3D)
persp3D(xx, yy, matrix(ks.local$predict, nrow = length(xx)), theta=-60, phi=40)
library(npsp)
spersp(xx, yy, ks.local$predict, theta=-60, phi=40)
######################
##VALIDACI?N CRUZADA##
######################

EXPO<-vgm(psill = 1.5, "Exp", range = 60000,nugget = 0.93)
MAT<-vgm(psill = 2.063, "Mat", range = 67000,nugget = 0.2,kappa = 0.64)
KO.exp.cv.z <- krige.cv(temp.trans~1, ~E+N, croacia2, EXPO, nmin=0, nmax=30)    
KS.mat.cv.z <- krige.cv(temp.trans~1, ~E+N, croacia2,MAT, beta=mu, nmin=0, nmax=154)

resultados.cv.z <- rbind(criterio.cv(KO.exp.cv.z), criterio.cv(KS.mat.cv.z))
rownames(resultados.cv.z) <- c("KO.esf.cv.z", "KS.esf.cv.z")
resultados.cv.z
xtable(resultados.cv.z)

        ##################################################
        #####M?todos de interpolaci?n deterministicos#####
        ##################################################

#####################################
#####DISTANCIA INVERSA PONDERADA#####
#####################################

#P optimizado para el IDW
p.optimo <- function(p, formula, locations, data, newdata, nmax, nmin, maxdist, var.reg){
  idw.pred <- as.data.frame(matrix(NA,nrow= nrow(data), ncol=4))
  colnames(idw.pred) <- c("x","y","var1.pred","var1.var")
  for(i in 1:(nrow(data))){
    idw.pred[i,] <- idw(formula, locations, data[-i,], newdata[i,], nmax, nmin, maxdist, idp=p)
  } 
  RMSPE <-  sqrt(sum((idw.pred$var1.pred-var.reg)^2)/nrow(data))
  RMSPE
}
temp<-croacia2$temp

P <- optimize(p.optimo, c(0,10), formula=temp~1, locations=~E+N, data=croacia2, newdata=croacia2, nmax=10, nmin=10, maxdist=Inf, var.reg=temp)
cat("Par?metro ?ptimo IDW: ", "\n", "P       = ", P$minimum, "\n", "RMSPE   = ", P$objective, "\n")

p1 <- p.optimo(p=1, formula=temp~1, locations=~E+N, data=croacia2, newdata=croacia2, nmax=10, nmin=10, maxdist= Inf, var.reg=temp)
p2 <- P$objective
p3 <- p.optimo(p=2, formula=temp~1, locations=~E+N, data=croacia2, newdata=croacia2, nmax=10, nmin=10, maxdist= Inf, var.reg=temp)
p4 <- p.optimo(p=3, formula=temp~1, locations=~E+N, data=croacia2, newdata=croacia2, nmax=10, nmin=10, maxdist= Inf, var.reg=temp)
p5 <- p.optimo(p=4, formula=temp~1, locations=~E+N, data=croacia2, newdata=croacia2, nmax=10, nmin=10, maxdist= Inf, var.reg=temp)
p6 <- p.optimo(p=5, formula=temp~1, locations=~E+N, data=croacia2, newdata=croacia2, nmax=10, nmin=10, maxdist= Inf, var.reg=temp)

RMSPE <- matrix(c(p1, p2, p3, p4, p5, p6))
rownames(RMSPE)<-c("P1","P2","P3","P4","P5","P6")
colnames(RMSPE)<-c("RMSPE")
xtable(t(RMSPE))
x11()
plot(c(1,P$minimum,2,3,4,5),RMSPE, main="Gr?fico de optimizaci?n del Par?metro (P)\n Distancia Inversa Ponderada", ylab="RMSPE", xlab="p ?ptimo =1.328093", type="b",col=2)

#IDW
knowndt = data.frame(croacia2$E,croacia2$N,croacia2$temp)
names(knowndt)<-c("x","y","temp")
coordinates(knowndt)  <- ~ x+y
idw.temp <- idw(temp~1, knowndt, SPDF3, maxdist = Inf, idp=1.328093)
##ESTAD?STICOS DE LA PREDICCI?N IDW
estadisticosIDW<-matrix(c(min(idw.temp$var1.pred),mean(idw.temp$var1.pred),max(idw.temp$var1.pred)))
rownames(estadisticosIDW)<-c("min","media","max")
colnames(estadisticosIDW)<-c("IDW")
xtable(t(estadisticosIDW))
##se realiza un cambio de los valores predecidos en IDW con krige para poder graficar
#ya que como NO se uso la grilla suministrada porque estaba mal proyectada (utm 32N)
idwplot<- krige.conv(geodata.croacia, loc = pred.grid, krige = krige.control(obj.m = exp.wls))##krige sin importancia, solo para cambiar los valores de IDW
idwplot$predict<- idw.temp$var1.pred ##CAMBIO DE VALORES DE IDW CON KRIGE MOMENTANEO

x11() ## GENERACION DE MAPAS IDW
image(idwplot,col = heat.colors(30)) #superficie de predicci?n
map(shape, add = T, col = "sienna4",lwd =3)
title("Predicciones Distancia inversa ponderada con P?ptimo")
points(geodata.croacia$coords, pch=20, nlev=21) #a?adir posiciones datos
contour(idwplot,add=T, nlev=10) #a?adir gr?fico de contorno
legend(locator(1), c("14","15","16","17","18","19"), bty = "o", title = "Temperatura ?C", fill = heat.colors(12), xpd = NA)

#############################################################
###################FUNCIONES DE BASE RADIAL##################
#############################################################
croacia2
coordinates(croacia2) <- ~E+N
spsample()
##OPTIMIZACI?N DE LOS PAR?METROS ETA Y RHO, SOLO SE DEBE CAMBIAR func="ID DE LA FUNCION"
args(graph.rbf)
a<-graph.rbf(temp~1, croacia2, eta.opt=T, rho.opt=F, n.neigh=30,np = 154, func="TPS",
          eta.dmax=2, iter=600,P.T = T)
##SELECCI?N DE LA FUNCI?N DE BASE RADIAL A USAR RMSPE MAS BAJO.
RMSPE.M<-rbf.cv(temp~1, croacia2, eta=1e-05, rho=0, n.neigh=30, func="M")
RMSPE.IM<-rbf.cv(temp~1, croacia2, eta=100, rho=0, n.neigh=30, func="IM")
RMSPE.ST<-rbf.cv(temp~1, croacia2, eta=0.013, rho=0, n.neigh=30, func="ST")
RMSPE.CRS<-rbf.cv(temp~1, croacia2, eta=0.009025519, rho=0, n.neigh=30, func="CRS")
RMSPE.TPS<-rbf.cv(temp~1, croacia2, eta=1.43, rho=0, n.neigh=30, func="TPS")
RMSPE.EXP<-rbf.cv(temp~1, croacia2, eta=0.03601307, rho=0, n.neigh=30, func="EXPON")
RMSPE.GAU<-rbf.cv(temp~1, croacia2, eta=0.01, rho=0, n.neigh=30, func="GAU")

RMSPE.rbf <- matrix(c(RMSPE.M, RMSPE.IM, RMSPE.ST, RMSPE.CRS, RMSPE.TPS,RMSPE.EXP,RMSPE.GAU))
rownames(RMSPE.rbf)<-c("MULTICUADRATICA","INVERSA MULTICUADRATICA","SPLINE CON TENSI?N","COMPLETAMENTE REGULARIZADA SPLINE","SPLINE CAPA DELGADA","EXPONENCIAL","GAUSSIANA")
colnames(RMSPE.rbf)<-c("RMSPE")
RMSPE.rbf
xtable((RMSPE.rbf))

# PREDICCI?N EN LA GRILLA DE PUNTOS GENERADA, CON LA FUNCION DE BASE RADIAL SPLINE CON TENSI?N
pred.rbf <- rbf(temp~1, croacia2, eta=0.009025519, rho=0, newdata= pred.grid,
                n.neigh=30, func="CRS")

##ESTAD?STICOS DE LA PREDICCI?N RBF
estadisticosRBF<-matrix(c(min(pred.rbf$var1.pred),mean(pred.rbf$var1.pred),max(pred.rbf$var1.pred)))
rownames(estadisticosRBF)<-c("min","media","max")
colnames(estadisticosRBF)<-c("RBF")
xtable(t(estadisticosRBF))

##se realiza un cambio de los valores predecidos en RBF con krige para poder graficar
#ya que como NO se uso la grilla suministrada porque estaba mal proyectada (utm 32N)
RBF.M<- krige.conv(geodata.croacia, loc = pred.grid, krige = krige.control(obj.m = exp.wls))##krige sin importancia, solo para cambiar los valores de IDW
RBF.M$predict<- pred.rbf$var1.pred ##CAMBIO DE VALORES DE RBF CON KRIGE MOMENTANEO

x11() ## GENERACION DE MAPAS RBF
image(RBF.M) #superficie de predicci?n
map(shape, add = T, col = "sienna4",lwd =3)
title("Predicciones Spline completamente regularizada CRS")
points(geodata.croacia$coords, pch=20, nlev=21) #a?adir posiciones datos
contour(RBF.M,add=T, nlev=10) #a?adir gr?fico de contorno
legend(locator(1), c("14","15","16","17","18","19"), bty = "o", title = "Temperatura ?C", fill = heat.colors(12), xpd = NA)

#############################################################################
################################DISE?O DE RED################################
#############################################################################

assign("network.design",
       function(formula, model, npoint.x, npoint.y, npoints, boundary=NULL, mu, type="regular", ...){
         if (is.null(boundary)) {
           grid.p<-expand.grid(x=seq(min(x),max(x),(max(x)-min(x))/npoint.x), y=seq(min(y),max(y),(max(y)-min(y))/npoint.y))
           plot(grid.p,pch=19,cex=0.5)
           grid.p$z <- grid.p$x
         }
         else if (is.null(boundary)==FALSE) {
           df.pts<-spsample(boundary, n=npoints, type=type)
           plot(boundary,axes=T,border="Blue",main="Dise?o de Red basado en ASEPE")
           points(df.pts,pch=19,cex=0.5)
           grid.p <- data.frame(df.pts,df.pts@coords[,1])
           names(grid.p) <- c("x","y","z")
         }
         K = krige.cv(formula, ~x+y, grid.p, model, ...)
         ASEPE <- mean((K[,2])^0.5)             
         ASEPE
       }
)

borde<-readOGR("/BBDD/Croacia/Croatia.shp")
x <- 0:10
y<- 0:10

NDP1 <- network.design(x~x+y,model=EXPO, npoints=100, boundary=borde, nmax=6, type="stratified")
NDP2 <- network.design(x~1, model=EXPO, npoints=100, boundary=borde, nmax=6, type="stratified")
NDP3 <- network.design(x~1, model=EXPO, npoints=100, boundary=borde, nmax=8, type="stratified")
NDP4 <- network.design(log(x)~x+ y + x*y + I(x^2)+I(y^2),model=EXPO, npoints=100, boundary=borde, nmax=12, type="stratified")
x11()
Networks.P <- rbind(NDP1,NDP2,NDP3,NDP4)
colnames(Networks.P) <- c("ASEPE")

xtable(Networks.P)
