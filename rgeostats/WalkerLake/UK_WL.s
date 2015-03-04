library(fields)
library(RGeoS)
library(plotrix)
source('cubic.s')
source('dico.s')

set.seed(3009)
	
#nDEM = scan('data/nDEM',skip = 3)
#nDEM = matrix(nDEM,286)

WL = scan('WalkerLake',skip  = 3)

nx = 260
x = seq(1,nx)
ny = 300
y = seq(1,ny)
g = expand.grid(x,y)
zl =  range(WL)
image.plot(x,y,matrix(WL,260))
# smoothing
#wl_mav = smooth.2d(log(1+WL[,1]),g,theta=2,nrow= 260,ncol=300)$z
#zl = range(wl_mav)
#image.plot(x,y,wl_mav)

# sampling
n = 1000	
sam = sample(1:length(WL),n)
points(g[sam,])
# data
dat = WL[sam]
#dat = wl_mav[sam]

# Thin plate spline regression
tps_trend = Tps(g[sam,],dat,lambda = 0.005)
summary(tps_trend)
set.panel(2,2)
plot(tps_trend)
x11()
surface(tps_trend)

rmap = predict(tps_trend,g)
image.plot(x,y,matrix(rmap,nx),zlim=zl)


	
#image.plot(x,y,rmap-wl_mav)
#hist(as.vector(rmap - wl_mav))
#var(as.vector(rmap - wl_mav))

	
# analyse des résidus
dbglm = db.create(x1 = g[sam,1],x2 = g[sam,2],z1 = tps_trend$res )
gamma = vario.calc(dbglm,lag = 5, nlag = 20)
mod = model.auto(gamma,struct =melem.name(c(1,2)),npairpt=T)

print(mod)
	
d =  sqrt(outer(g[sam,1],g[sam,1],'-')^2+outer(g[sam,2],g[sam,2],'-')^2)
C = matrix(model.eval(mod,d,as.cov=T),n)
u = solve(C,tps_trend$res)
kmap = 0*rmap
for (i in 1:(nx*ny)){
	h=sqrt((g[sam,1]-g[i,1])^2+(g[sam,2]-g[i,2])^2)
	c0 = model.eval(mod,h,as.cov=T)
	kmap[i]=rmap[i]+t(c0)%*%u
}

image.plot(x,y,matrix(kmap,nx),zlim = zl)
print(1/(nx*ny)*sum((kmap-WL)^2))
#mse = 1/(nx*ny)*sum((kmap-wl_mav)^2)


###############################################################
## comparison mn-kriging
dbdat = db.create(x1 = g[sam,1],x2 = g[sam,2],z1 = dat,ndim=2)
gammadat = vario.calc(dbdat,lag = 10, nlag = 20)
modeldat = model.auto(gammadat,struct =melem.name(c(1,3)),npairpt=T)
sig0 = model.eval(modeldat,0,as.cov=T)

COVdat = matrix(model.eval(modeldat,d,as.cov=T),n)
COVdat = cbind(COVdat,rep(1,n))
COVdat = rbind(COVdat,c(rep(1,n),0))


print(modeldat)

## continuous moving neighborhood
ray = 70
al = ray - 10
be = ray + 10
kmap_vgc = 0*kmap
mmap_vgc = kmap
varmap_vgc = kmap

	
for (i in 1:(nx*ny)){
h=sqrt((g[sam,1]-g[i,1])^2+(g[sam,2]-g[i,2])^2)
vois = which(h<be)
if (length(vois)>0){
zeros = c(rep(0,length(vois)),1)
c0 = c(model.eval(modeldat,h[vois],as.cov=T),1)
W = sig0*diag(c(4*(h[vois]-al)^2/(be-h[vois])^2*(h[vois]>al),0))
if (dim(W)[1]==0){W=0}
#u=t(c0)%*%solve(COVdat[c(vois,n+1),c(vois,n+1)]+W)
invC = solve(COVdat[c(vois,n+1),c(vois,n+1)]+W)
u=t(c0)%*%invC
kmap_vgc[i]=u[-length(u)]%*%dat[vois]
mmap_vgc[i] = (zeros%*%invC)[-length(u)]%*%dat[vois]
varmap_vgc[i]=(sig0-u%*%c0)
}
print(i)
}

image.plot(x,y,matrix(kmap_vgc,nx),zlim=zl)
print(1/(nx*ny)*sum((kmap_vgc-WL)^2))

image.plot(x,y,mmap_vgc,zlim=zl)
points(g[sam,],col=tim.colors(64)[cut(dat,64,labels=F)])
var(dat - mmap_vgc[sam])

write.table(kmap_vgc,'kmap_cok_1000.dat',row.names=F)
write.table(rmap,'rmap_1000.dat',row.names=F)
write.table(kmap,'kmap_uk_1000.dat',row.names=F)
write.table(cbind(g[sam,],dat),'samples_1000.dat',row.names=F)
