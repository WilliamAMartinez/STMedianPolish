\name{VRes}
\alias{VRes}
\docType{data}
\title{
Precomputed Variogram for residuals of monthly precipitation
}
\description{
Precomputed Variogram for residuals of monthly precipitation \code{\link{Metadb}}. For this, the 'variogram' [package "gstat"] function is used. 
}
\usage{data(VRes)}
\format{
  The format is: 'StVariogram' 'data.frame'
}
\examples{
## Not run:
library(zoo)
library(gstat)
library(sp)
library(spacetime)
data(Metadb)
x<-matrix(0,1,121)
for(i in 1:121){
  x[,i] <- 2000 + (seq(0, 120)/12)[i]
}
x<-as.Date (as.yearmon(x), frac = 1)
time = as.POSIXct(x, tz = "GMT")
MPST<-ConstructMPst(sqrt(0.5+Metadb[,-c(1:4)]),time,pts=Metadb[,2:4],Delta=c(7,6,5))
residuals<-removetrendMPst(MPST,eps=0.01, maxiter=2)
rain.loc<-Metadb[,c("Station","East","North","Height")]
coordinates(rain.loc) = ~East+North+Height
proj4string(rain.loc) = CRS(proj4string(DemMeta))
rain_residual = stConstruct(data.frame(Res=residuals[,7]), space = list(values = 1),
                            time, SpatialObj = rain.loc,interval=TRUE)
#Empirical variogram
VRes = variogram(values ~ 1, rain_residual, cutoff=90000,tlags=0:24,width=2650,
                 assumeRegular=TRUE, na.omit=TRUE)			   
plot(VRes)
## End(Not run)		   
#Fit variogram
data(VRes)
FitPar_st = function(p, gfn, v, trace = FALSE, ...) {
  mod = gfn(v$spacelag, v$timelag,p,  ...)
  resid = v$gamma - mod
  if (trace)
    print(c(p, MSE = mean(resid^2)))
  mean(resid^2)
}
ModSpatial = function(h,p){p[2]*(1-exp(-h/p[3]))}
ModTemporal = function(u,p){p[4]*(1-exp(-u/p[5]))+ p[6]*(1-cos(pi*u/180))+
							p[7]*abs(sin(pi*u/180))}
VariogST=function(h,u,p)
		{ModTemporal(u,p)+ModSpatial(h,p)+p[8]*ModTemporal(u,p)*ModSpatial(h,p)}
#Parametros Iniciales
p1<-c(2,14.5,13900,5.9,29,1.55,3.7,-0.07)
pars.st = optim(p1, FitPar_st, method =  "BFGS",
                gfn = VariogST, v = VRes, hessian=TRUE)
fit_Variog_ST<-VRes
fit_Variog_ST$gamma<-VariogST(VRes$spacelag, VRes$timelag, pars.st$par)
plot(fit_Variog_ST)
}
\keyword{datasets}
\references{
Martínez, W. A., Melo, C. E., & Melo, O. O. (2017). \emph{Median Polish Kriging for space--time analysis of precipitation} Spatial Statistics, 19, 1-20. \href{http://www.sciencedirect.com/science/article/pii/S2211675316301336}{[link]}

Pebesma, E.J., 2004. \emph{Multivariable geostatistics in S: the gstat package}. Computers & Geosciences, 30: 683-691 \href{https://cran.r-project.org/web/packages/gstat/index.html}{[link]}
}




