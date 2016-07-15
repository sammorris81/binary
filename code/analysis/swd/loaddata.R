#Krishna Pacifici
#July 15, 2013
#get data formatted
library(fields)
dat<-read.csv("SWD.csv")

get_a_var<-function(dat,index){
  D<-dat[,index]
  D<-matrix(D,200,200,byrow=T)
  D<-t(D[200:1,])
}

Y1      <- get_a_var(dat,10)
Y2      <- get_a_var(dat,9)
tci     <- get_a_var(dat,3)
solar_s <- get_a_var(dat,4)   
slope   <- get_a_var(dat,5)
road_di <- get_a_var(dat,6)
elev    <- get_a_var(dat,7)


X1 <-(tci-mean(tci))/sd(as.vector(tci))
X2 <-(solar_s-mean(solar_s))/sd(as.vector(solar_s))
X3 <-(slope-mean(slope))/sd(as.vector(slope))
X4 <-(road_di-mean(road_di))/sd(as.vector(road_di))
X5 <-(elev-mean(elev))/sd(as.vector(elev))


rm(get_a_var,dat,tci,solar_s,slope,road_di,elev)

