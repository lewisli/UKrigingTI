library(RGeostats)
library(fields)
source('../cubic.s')

# data 
DataPath="/Volumes/Data/Data/UK-TI/"
ref = scan(paste(DataPath,'reference/Stationary',sep=''),skip = 3)

plot.new()
frame()
image(matrix(ref,150))
	
number = c(25)
file = 'hardData'

# Simulation Grid Dimensions
nx = 150
ny = 150
x = seq(1,nx)
y = seq(1,ny)
g = expand.grid(x,y)

# Allocate output variance and kriged estimate map
kmap = matrix(0,nx,ny)
vmap = matrix(0,nx,ny)
kmapg = kmap
vmapg = vmap
dat = cbind(g,ref)
N = nx*ny
	
# Parameters
rg = 25
ray = 40
Ns = 100
mse_mv = matrix(0,Ns,length(number))
mse_g = mse_mv

dats = read.table(paste(DataPath,'hardData/Stationary/100/hData2',sep=''),skip=6)
dats = dats[,-3]
dbdat = db.create(x1 = dats[,1],x2 = dats[,2],z1 = dats[,3],ndim=2)
gammadat = vario.calc(dbdat,lag=10,nlag=10)
modeldat = model.auto(gammadat,struct =melem.name(c(1,3)),npairpt=T,draw=F)
sig0 = model.eval(modeldat,0,as.cov=T)
dist = sqrt(outer(dats[,1],dats[,1],'-')^2+outer(dats[,2],dats[,2],'-')^2)
n = Ns

COVdat = matrix(model.eval(modeldat,dist,as.cov=T),n)
COVdat = cbind(COVdat,rep(1,n))
COVdat = rbind(COVdat,c(rep(1,n),0))
COVinv = solve(COVdat)
ug = COVinv%*%c(dats[,3],0)

## continuous moving neighborhood
al = ray - 10
be = ray + 10
for (i in 1:(nx*ny))
{
	h=sqrt((dats[,1]-g[i,1])^2+(dats[,2]-g[i,2])^2)
	c0 = c(model.eval(modeldat,h,as.cov=T),1)
	kmapg[i] = t(c0)%*%ug
	vmapg[i] = sig0 - t(c0)%*%COVinv%*%c0
	vois = which(h<be)

	if (length(vois)>0)
	{
		zeros = c(rep(0,length(vois)),1)
		c0 = c0[c(vois,n+1)]
		W = sig0*diag(c(4*(h[vois]-al)^2/(be-h[vois])^2*(h[vois]>al),0))
		if (dim(W)[1]==0){W=0}
		COVinvdatv = solve(COVdat[c(vois,n+1),c(vois,n+1)]+W)
		u=t(c0)%*%COVinvdatv
		kmap[i]=u[-length(u)]%*%dats[vois,3]
		vmap[i]= sig0 - t(c0)%*%COVinvdatv%*%c0
	}
}

# image.plot(1:nx,1:ny,matrix(kmap,nx))
# image.plot(1:nx,1:ny,matrix(vmap,nx))
	
# write.csv(kmap,'kmap_cok.csv',row.names = F)
# write.csv(vmap,'vmap_cok.csv',row.names = F)

# image.plot(1:nx,1:ny,matrix(kmapg,nx))
# image.plot(1:nx,1:ny,matrix(vmapg,nx))
	
# write.csv(kmapg,'kmap_g.csv',row.names = F)
# write.csv(vmapg,'vmap_g.csv',row.names = F)


