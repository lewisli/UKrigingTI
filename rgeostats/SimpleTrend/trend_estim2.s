library(RGeoS)
library(MASS)
source('../cubic.s')

# data 
ref = read.csv('trend1.dat')#,skip = 1)

# grid

nx = 100
ny = 200	
x = seq(1,nx)
y = seq(1,ny)
g = expand.grid(x,y)
N = nx*ny
dat = data.frame(cbind(g,ref[,1],ref[,2]))
names(dat) = c('x','y','trend','Z')

# figure
	
#image(x,y,matrix(dat$Z,nx))
	
# sampling

number = c(25,50,100,150,200)

kmapo = matrix(0,nx,ny)
kmapg = kmapo
#kmapgv = kmapo

sig0 = 0.05
rg = 25
ray = 40
Ns = 100
mse_ols = matrix(0,Ns,length(number))
mse_gls = mse_ols
mse_gls_vario = mse_gls
	
for (n in number){
for (k in 1:Ns){
	sam = sample(1:N,n) # sampling
	dats = dat[sam,] # sampling data
	rego = lm(Z~x+y,data = dats) # OLS regression
	dist = sqrt(outer(dats[,1],dats[,1],'-')^2+outer(dats[,2],dats[,2],'-')^2) # distance matrix
	COVdat = matrix(sig0*expo(dist,25),n) # true covariance matrix 
	regg = lm.gls(Z~x+y,data = dats,W=COVdat) # GLS regression
	invCOVdat = chol2inv(chol(COVdat))
	# dual kriging
        uo = invCOVdat%*%(rego$residuals)
	# dual kriging
        ug = invCOVdat%*%(regg$residuals)
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
	invCOVdat = chol2inv(chol(COVdat))
	# dual kriging
        ugv = invCOVdat%*%(reg$res)
	for (i in 1:(nx*ny)){
		h=sqrt((dats[,1]-g[i,1])^2+(dats[,2]-g[i,2])^2)
		c0 = sig0*expo(h,rg)
		kmapo[i]=predict(rego,newdata = dat[i,])+t(c0)%*%uo
		kmapg[i]=regg$coefficients%*%c(1,dat$x[i],dat$y[i])+t(c0)%*%ug
		kmapgv[i]=reg$coefficients%*%c(1,dat$x[i],dat$y[i])+t(c0)%*%ugv
}
mse_ols[k,which(number==n)] = 1/(N-n)*sum((dat$Z[-sam]-kmapo[-sam])^2)
mse_gls[k,which(number==n)] = 1/(N-n)*sum((dat$Z[-sam]-kmapg[-sam])^2)
mse_gls_vario[k,which(number==n)] = 1/(N-n)*sum((dat$Z[-sam]-kmapgv[-sam])^2)
print(c(n,k))
}
}

x11()
image(x,y,kmapo,col=heat.colors(24))
points(dats[,1],dats[,2],pch='x')
dev.copy2eps(file='kriging_map_ols.eps')
image(x,y,kmapg,col=heat.colors(24))
points(dats[,1],dats[,2],pch='x')
dev.copy2eps(file='kriging_map_gls.eps')
image(x,y,kmapgv,col=heat.colors(24))
points(dats[,1],dats[,2],pch='x')
dev.copy2eps(file='kriging_map_gls_vario.eps')
	
mse = rbind(apply(mse_ols,2,mean),apply(mse_gls,2,mean),apply(mse_gls_vario,2,mean))
mse=mse/sig0
mse = rbind(number,mse)
	write.csv(round(mse,2),'mse2.csv',row.names = F)

vmse = rbind(apply(mse_ols/sig0,2,var),apply(mse_gls/sig0,2,var),apply(mse_gls_vario/sig0,2,var))
vmse = rbind(number,vmse)
write.csv(round(100*vmse,2),'vmse2.csv',row.names = F)

write.csv(mse_ols,'rmse_ols2.csv')
write.csv(mse_gls,'rmse_gls2.csv')
write.csv(mse_gls_vario,'rmse_gls_vario2.csv')
