library(RGeoS)
source('../cubic.s')

path = '/home/thomas/Dropbox/SlideingScheme/GPUTIKriging/MPSKriging/data/GaussianExample/'
# data 
ref = scan(paste(path,'ref',sep=''),skip = 3)
image(matrix(ref,150))
	
number = c(25,50,100,150,200)
file = 'hardData'

# Grid
nx = 150
x = seq(1,nx)
ny = 150
y = seq(1,ny)
g = expand.grid(x,y)
kmap = matrix(0,nx,ny)
kmapg = kmap
dat = cbind(g,ref)
N = 150^2
	
# Parameters
sig0 = 1
rg = 25
ray = 40
Ns = 100
mse_mv = matrix(0,Ns,length(number))
mse_g = mse_mv

	
for (n in number){
for (k in 1:Ns){
	sam = sample(1:N,n)
	dats = dat[sam,]
	dbdat = db.create(x1 = dats[,1],x2 = dats[,2],z1 = dats[,3],ndim=2)
	dist = sqrt(outer(dats[,1],dats[,1],'-')^2+outer(dats[,2],dats[,2],'-')^2)
	COVdat = matrix(sig0*sferic(dist,rg),n)
	COVdat = cbind(COVdat,rep(1,n))
	COVdat = rbind(COVdat,c(rep(1,n),0))
	#invCOVdat = chol2inv(chol(COVdat))
	ug = solve(COVdat,c(dats[,3],0))
	## continuous moving neighborhood
	al = ray - 10
	be = ray + 10
	for (i in 1:(nx*ny)){
		h=sqrt((dats[,1]-g[i,1])^2+(dats[,2]-g[i,2])^2)
		c0 = c(sig0*sferic(h,rg),1)
		kmapg[i] = t(c0)%*%ug
		vois = which(h<be)
		if (length(vois)>0){
			zeros = c(rep(0,length(vois)),1)
			c0 = c0[c(vois,n+1)]
			W = sig0*diag(c(4*(h[vois]-al)^2/(be-h[vois])^2*(h[vois]>al),0))
			if (dim(W)[1]==0){W=0}
			u=t(c0)%*%solve(COVdat[c(vois,n+1),c(vois,n+1)]+W)
			kmap[i]=u[-length(u)]%*%dats[vois,3]
		}
	}
mse_mv[k,which(number==n)] = 1/(N-n)*sum((dat[-sam,3]-kmap[-sam])^2)
mse_g[k,which(number==n)] = 1/(N-n)*sum((dat[-sam,3]-kmapg[-sam])^2)
print(k)
}
}

mse = rbind(apply(mse_mv,2,mean),apply(mse_g,2,mean))
mse = rbind(number,mse)
write.csv(round(mse,2),'mse2t.csv',row.names = F)

	
vmse = rbind(apply(mse_mv,2,var),apply(mse_g,2,var))
vmse = rbind(number,vmse)
write.csv(round(100*vmse,2),'vmse2t.csv',row.names = F)

