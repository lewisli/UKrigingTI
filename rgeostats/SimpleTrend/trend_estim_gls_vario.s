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

kmap = matrix(0,nx,ny)
varmap = kmap

sig0 = 0.05
rg = 25
ray = 40
Ns = 100
rmse = matrix(0,Ns,length(number))
	
for (n in number){
for (k in 1:Ns){
	sam = sample(1:N,n)
	dats = dat[sam,]
	reg = lm(Z~trend,data = dats) # first OLS regression
	dbdat = db.create(x1 = dats[,1],x2 = dats[,2],z1 = reg$res,ndim=2)
	gammadat = vario.calc(dbdat)
	modeldat = model.auto(gammadat,struct = melem.name(c(1,2)),npairpt=T,draw=F)
	dist = sqrt(outer(dats[,1],dats[,1],'-')^2+outer(dats[,2],dats[,2],'-')^2)
	COVdat = matrix(model.eval(modeldat,dist,as.cov=T),n)
	reg = lm.gls(Z~trend,data = dats,W=COVdat) # GLS regression
	dbdat = db.create(x1 = dats[,1],x2 = dats[,2],z1 = reg$res,ndim=2)
	gammadat = vario.calc(dbdat)
	modeldat = model.auto(gammadat,struct = melem.name(c(1,2)),npairpt=T,draw=F)
	COVdat = matrix(model.eval(modeldat,dist,as.cov=T),n)
	dats = cbind(dats,reg$fitted)
	invCOVdat = chol2inv(chol(COVdat))
        u = solve(COVdat,reg$res)
	for (i in 1:(nx*ny)){
		h=sqrt((dats[,1]-g[i,1])^2+(dats[,2]-g[i,2])^2)
		c0 = model.eval(modeldat,h,as.cov=T)
		kmap[i]=reg$coefficients%*%c(1,dat$trend[i])+t(c0)%*%u
		#varmap[i]=(sig0-u%*%c0)
	}
rmse[k,which(number==n)] = 1/(N-n)*sum((dat$Z[-sam]-kmap[-sam])^2)
print(k)
}
}
