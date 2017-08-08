##AUTHOR: Thomas Opitz (BioSP, INRA, Avignon, France; thomas AT opitz AT inra POINT fr)
## R code file accompanying the manuscript 
## "Point process-based modeling of multiple debris flow landslides using INLA: an application to the 2009 Messina disaster"
## Due to confidentiality reasons, no data are provided (the code cannot be executed "as is").


#provide folder paths (here for Unix-like operating systems) and some constants  ####
DATA="~/research/data/"
WORK="~/research/implementations/landslides/"
OUT="~/research/implementations/landslides/results/"

#the following constant has been used to rescale coordinate values (in m) for avoiding numerical instabilities
scalecoor=10^4

#extract significant effects from INLA fit 
extractSignificantEffects=function(results.df){
  print(results.df[sign(results.df$`0.025quant`)==sign(results.df$`0.975quant`),])
}

#load and preprocess data ####
setwd(DATA)

#load data: here given in matrix format (one line per pixel)
mydata=read.delim("DataMatrix.txt",header=TRUE)
#names(mydata) gives the following variable names:
#number of landslide events observed at each pixel (integer): status 
#"X"          "Y"          "status"          "Dist2Fault"          "NDVI"       "RSP"     "Slope"  
#"SPI"        "TWI"        "idsu"       "Aspect"          "Landform"          "Landuse"   
#"Plan_Cur"   "Prof_Cur"   "Lithology"  "DEM"  
#Description of variables:
#coordinates in m: X,Y
#numerical covariates: Dist2Fault,NDVI,RSP, Slope,SPI,TWI,Plan_Cur,Prof_Cur,DEM
#categorical covariates: Landform,Landuse,Lithology, Aspect (16 classes)
#slope unit identifier: idsu

#rescale some of the numerical covariates to have mean 0 and variance 1
mydatascale=mydata
vars2scale=setdiff(names(mydatascale),c("X","Y", "status", "idsu","Aspect","Landform","Landuse","Lithology"))
mydatascale[,vars2scale]=apply(mydatascale[,vars2scale],2,scale)
#create data frame for use with inla(...)
covar.inla=mydatascale[,c("Aspect", "DEM","Dist2Fault","Lithology", "Landuse", "Landform","NDVI","Plan_Cur","Prof_Cur","RSP","Slope","SPI","TWI","idsu")]
covar.inla=cbind(intercept=1,covar.inla)
#regroup a number of Lithology categories with very few occurrences into a single class with identifier 0
doReplace=covar.inla$Lithology %in% c(8,14,15,16,17,18,20,21,22,23,24,29)
covar.inla$Lithology[doReplace]=0
#use rm(mydatascale,mydata) if you want to remove those data frames for memory reasons 

#for the numerical covariates where nonlinear effects are checked, we discretize observations into 20 categories using functionality from R-INLA:
library(INLA) #see www.r-inla.org for installation
covar.inla$DEMlevels=inla.group(covar.inla$DEM,n=20)
covar.inla$Dist2Faultlevels=inla.group(covar.inla$Dist2Fault,n=20)
covar.inla$NDVIlevels=inla.group(covar.inla$NDVI,n=20)
covar.inla$proflevels=inla.group(covar.inla$Prof_Cur,n=20)
covar.inla$planlevels=inla.group(covar.inla$Plan_Cur,n=20)
covar.inla$Slopelevels=inla.group(covar.inla$Slope,n=20)
covar.inla$SPIlevels=inla.group(covar.inla$SPI,n=20)
covar.inla$TWIlevels=inla.group(covar.inla$TWI,n=20)


#prepare vector of Poisson responses for INLA
#we have to rescale the Poisson mean by an multiplicate "offset" according to the area of each grid cell
#(Note: the original grid cells are 15mx15m)
E.poisson=rep((15/scalecoor)^2,nrow(covar.inla))

#define response vector
y.poisson=covar.inla$status

#to avoid problems with intercept values that are far from being O(1), we will include a fixed offset, named "offs", into the linear (gaussian) predictor
offs=log(sum(y.poisson))


#fitting with INLA ####

#put all data elements into a "data stack"
stack=inla.stack(data=list(y=y.poisson,e=E.poisson),A=list(1),effects=list(covar.inla))

#fit model with fixed effects and nonlinear Aspect:
#create model formula: 
form=y~-1+intercept+DEM+Dist2Fault+NDVI+Plan_Cur+Prof_Cur+Slope+SPI+TWI+#fixed effects
  f(Landform,model="iid",hyper=list(theta=list(initial=log(10^2),fixed=TRUE)),constr=TRUE)+
  f(Landuse,model="iid",hyper=list(prec=list(initial=log(10^2),fixed=TRUE)),constr=TRUE)+
  f(Lithology,model="iid",hyper=list(prec=list(initial=log(10^2),fixed=TRUE)),constr=TRUE)+
  f(Aspect,model="rw1",cyclic=TRUE,constr=TRUE,hyper=list(theta=list(initial=log(5^2),fixed=TRUE)))
#run INLA (ATTENTION: this can be very time and memory consuming depending on the size of your data set and the fitted model)
#(for our data, it has been run on machines with >= 32GB of memory)
#we modify the prior of the fixed effect coefficients and of the intercept in control.fixed(...)
#the multiplicative offset must be indicated through the "E" argument
#num.threads allows to fix the (maximum) number of threads to be executed in parallel
#(each thread has memory requirement, hence fixing an upper bound may be useful to avoid crashing)
fit=inla(form, 
         data=inla.stack.data(stack), 
         family="poisson",
         offset=offs,
         control.fixed=list(prec=2,prec.intercept=1,mean.intercept=-2),
         E=inla.stack.data(stack)$e,
         verbose=TRUE,
         num.threads=2)
#show summary of fitted model
summary(fit)

#fit model with fixed and nonlinear effects (Model 2 in paper):
form=y~-1+intercept+DEM+Dist2Fault+NDVI+Plan_Cur+Prof_Cur+Slope+SPI+TWI+
  f(Landuse,model="iid",hyper=list(prec=list(initial=log(10^2),fixed=TRUE)),constr=TRUE)+
  f(Landform,model="iid",hyper=list(prec=list(initial=log(10^2),fixed=TRUE)),constr=TRUE)+
  f(Lithology,model="iid",hyper=list(prec=list(initial=log(10^2),fixed=TRUE)),constr=TRUE)+
  f(Aspect,model="rw1",cyclic=TRUE,constr=TRUE,hyper=list(theta=list(initial=log(5^2),fixed=TRUE)))+
  f(DEMlevels,model="rw1",hyper=list(prec=list(initial=log(5^2),fixed=TRUE)),constr=TRUE)+
  f(Dist2Faultlevels,model="rw1",hyper=list(prec=list(initial=log(5^2),fixed=TRUE)),constr=TRUE)+
  f(NDVIlevels,model="rw1",hyper=list(prec=list(initial=log(5^2),fixed=TRUE)),constr=TRUE)+
  f(planlevels,model="rw1",hyper=list(prec=list(initial=log(5^2),fixed=TRUE)),constr=TRUE)+
  f(proflevels,model="rw1",hyper=list(prec=list(initial=log(5^2),fixed=TRUE)),constr=TRUE)+
  f(Slopelevels,model="rw1",hyper=list(prec=list(initial=log(5^2),fixed=TRUE)),constr=TRUE)+ 
  f(SPIlevels,model="rw1",hyper=list(prec=list(initial=log(5^2),fixed=TRUE)),constr=TRUE)+ 
  f(TWIlevels,model="rw1",hyper=list(prec=list(initial=log(5^2),fixed=TRUE)),constr=TRUE)
#run INLA 
fit=inla(form, 
         data=inla.stack.data(stack), 
         family="poisson",
         offset=offs,
         control.fixed=list(prec=2,prec.intercept=1,mean.intercept=-2),
         E=inla.stack.data(stack)$e,
         verbose=TRUE,
         num.threads=2)
summary(fit)

#fit Model with partially nonlinear covariate effects and spatial effect (Model 3 in paper):
#the file adjgraph.txt contains the graph structure of the slope units
#in our case, it has been produced beforehand using the following commands:
#library(spdep)
#counties <- readShapeSpatial("slu.shp")
#nb2INLA(paste0(DATA,"adjgraph.txt"),poly2nb(counties,queen=F,row.names=counties$label))
form=y~-1+intercept+DEM+Dist2Fault+NDVI+Plan_Cur+Prof_Cur+Slope+SPI+TWI+
  f(Landuse,model="iid",hyper=list(prec=list(initial=log(10^2),fixed=TRUE)),constr=TRUE)+
  f(Landform,model="iid",hyper=list(prec=list(initial=log(10^2),fixed=TRUE)),constr=TRUE)+
  f(Lithology,model="iid",hyper=list(prec=list(initial=log(10^2),fixed=TRUE)),constr=TRUE)+
  f(Aspect,model="rw1",cyclic=TRUE,constr=TRUE,hyper=list(theta=list(initial=log(5^2),fixed=TRUE)))+
  f(DEMlevels,model="rw1",hyper=list(prec=list(initial=log(5^2),fixed=TRUE)),constr=TRUE)+
  f(Slopelevels,model="rw1",hyper=list(prec=list(initial=log(5^2),fixed=TRUE)),constr=TRUE)+f(Dist2Faultlevels,model="rw1",hyper=list(prec=list(initial=log(5^2),fixed=TRUE)),constr=TRUE)+
  f(idsu,model="besag",graph=paste0(DATA,"adjgraph.txt"),hyper=list(theta=list(initial=log(1),fixed=FALSE,prior="loggamma",param=c(.25,.25))))
#run INLA 
fit=inla(form, 
         data=inla.stack.data(stack), 
         family="poisson",
         offset=offs,
         control.compute=list(dic=TRUE,cpo=TRUE,waic=TRUE,mlik=TRUE,hyperpar=TRUE,return.marginals=TRUE),
         control.fixed=list(prec=2,prec.intercept=1,mean.intercept=-2),
         E=inla.stack.data(stack)$e,
         verbose=TRUE,
         num.threads=2)
#show summary of fitted model
summary(fit)




