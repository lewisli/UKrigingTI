library(RGeoS,lib.loc = '/bibli/R64')
library(MASS)
library(fields)
source('../cubic.s')

# data 
ref = read.csv('trend1.dat')#,skip = 1)
path = '~/Dropbox/SlideingScheme/Data/'

	
# grid

nx = 100
ny = 200	
x = seq(1,nx)
y = seq(1,ny)
g = expand.grid(x,y)
N = nx*ny

# figure
	
#image.plot(x,y,matrix(dat$Z,nx))
	
# sampling

number = c(25,50,100,150,200)

kmapo = matrix(0,nx,ny)
kmapg = kmapo
kmapgv = kmapo
vmapo = kmapo
vmapg = kmapo
vmapgv = kmapo
	
sig0 = 0.05
rg = 25
ray = 40
Ns = 100
mse_ols = matrix(0,Ns,length(number))
mse_gls = mse_ols
mse_gls_vario = mse_gls
	
dats = read.table(paste(path,'hardData/DS-SimpleTrend/100/hData3',sep=''),skip=6)
dats = dats[,-3]
names(dats) = c('x','y','Z')
	
	rego = lm(Z~x+y,data = dats) # OLS regression
	dist = sqrt(outer(dats[,1],dats[,1],'-')^2+outer(dats[,2],dats[,2],'-')^2) # distance matrix
	COVdat_t= matrix(sig0*expo(dist,25),n) # true covariance matrix 
	regg = lm.gls(Z~x+y,data = dats,W=COVdat) # GLS regression
	invCOVdat_t = solve(COVdat)
	# dual kriging
        uo = invCOVdat_t%*%(rego$residuals)
	X = as.matrix(dats[,1:2])
	xx1x = solve(t(X)%*%X)%*%t(X)
	c1x = invCOVdat_t%*%X
	cx_ols = xx1x%*%COVdat_t%*%t(xx1x)
	# dual kriging
        ug = invCOVdat_t%*%(regg$residuals)
	cx_gls = chol2inv(chol(t(X)%*%invCOVdat_t%*%X))
	# variogram estimation from residuals
	dbdat = db.create(x1 = dats[,1],x2 = dats[,2],z1 = rego$res,ndim=2)
	gammadat = vario.calc(dbdat)
	modeldat = model.auto(gammadat,struct =melem.name(c(1,2)),npairpt=T,draw=F)
	COVdat = matrix(model.eval(modeldat,dist,as.cov=T),n)
	reg = lm.gls(Z~x+y,data = dats,W=COVdat) # GLS regression
	# variogram estimation from better residuals
	dbdat = db.create(x1 = dats[,1],x2 = dats[,2],z1 = reg$res,ndim=2)
	gammadat = vario.calc(dbdat)
	modeldat = model.auto(gammadat,struct =melem.name(c(1,2)),npairpt=T,draw=F)
	COVdat = matrix(model.eval(modeldat,dist,as.cov=T),n)
	invCOVdat_e = chol2inv(chol(COVdat))
	# dual kriging
        ugv = invCOVdat_e%*%(reg$res)
	cx_glse = chol2inv(chol(t(X)%*%invCOVdat_e%*%X))
	sig0e = model.eval(modeldat,0,as.cov=T)
	for (i in 1:(nx*ny)){
		h=sqrt((dats[,1]-g[i,1])^2+(dats[,2]-g[i,2])^2)
		c0 = model.eval(modeldat,h,as.cov=T)
		kmapo[i]=predict(rego,newdata = dat[i,])+t(c0)%*%uo
		cc = invCOVdat_t%*%c0
		vmapo[i]=sig0 - t(c0)%*%invCOVdat_t%*%c0 + as.matrix(g[i,]-t(c0)%*%c1x)%*%cx_ols%*%t(as.matrix(g[i,]-t(c0)%*%c1x))
		kmapg[i]=regg$coefficients%*%c(1,dat$x[i],dat$y[i])+t(c0)%*%ug
		vmapg[i]=sig0 - t(c0)%*%invCOVdat_t%*%c0 + as.matrix(g[i,]-t(c0)%*%c1x)%*%cx_gls%*%t(as.matrix(g[i,]-t(c0)%*%c1x))
		c0 = model.eval(modeldat,h,as.cov=T)
		kmapgv[i]=reg$coefficients%*%c(1,dat$x[i],dat$y[i])+t(c0)%*%ugv
		vmapgv[i]=sig0e - t(c0)%*%invCOVdat_t%*%c0 + as.matrix(g[i,]-t(c0)%*%c1x)%*%cx_glse%*%t(as.matrix(g[i,]-t(c0)%*%c1x))
}		


x11()
image.plot(x,y,kmapo,col=heat.colors(24))
#points(dats[,1],dats[,2],pch='x')
#dev.copy2eps(file='kriging_map_ols.eps')
image.plot(x,y,kmapg,col=heat.colors(24))
#points(dats[,1],dats[,2],pch='x')
#dev.copy2eps(file='kriging_map_gls.eps')
image.plot(x,y,kmapgv,col=heat.colors(24))
#points(dats[,1],dats[,2],pch='x')
#dev.copy2eps(file='kriging_map_gls_vario.eps')
	
write.csv(kmapo,'kmap_ols.csv',row.names = F)
write.csv(vmapo,'vmap_ols.csv',row.names = F)

write.csv(kmapg,'kmap_gls.csv',row.names = F)
write.csv(vmapg,'vmap_gls.csv',row.names = F)

write.csv(kmapgv,'kmap_glse.csv',row.names = F)
write.csv(vmapgv,'vmap_glse.csv',row.names = F)

	