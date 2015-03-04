matern=function(h,a,nu=2)
{h=2*sqrt(nu)*h/a
sol=vector("numeric",length(h))
i=(h==0)
sol=(1/(gamma(nu)*2**(nu-1))*(h)**nu*besselK(h, nu))
sol[i]=1
sol
}
	
storkey=function(h,a)
{h=h/a
(2*(1-h)*(1+(cos(2*pi*h))/2)+3/(2*pi)*sin(2*pi*h))/3*as.numeric(h<=1)
}
	
penta=function(h,a)
{h=h/a
(1-22/3*h^2+33*h^4-77/2*h^5+33/2*h^7-11/2*h^9+5/6*h^11)*as.numeric(h<=1)
}

cubic=function(h,a)
{h=h/a
(1-7*h^2+35/4*h**3-7/2*h**5+3/4*h**7)*as.numeric(h<=1)
}	

sferic=function(h,a)
{h=h/a
(1-3/2*h+1/2*h^3)*as.numeric(h<=1)
}

gauss=function(h,a)
{h=h/a
(exp(-3*h^2))
}

sincard=function(h,a)
{h=h/a
c=sin(h)/h
c[is.nan(c)]=1
c
}

expo=function(h,a)
{h=h/a
(exp(-3*h))
}

wend1=function(h,a)
{h=h/a
((1-10*h**2+20*h**3-15*h**4+4*h**5)*as.numeric(h<=1))
}
	
wend2=function(h,a)
{h=h/a
((35/3*h**8-64*h**7+140*h**6-448/3*h**5+70*h**4-28/3*h**2+1)*as.numeric(h<=1))
}

tap=function(h,mod1,a1,mod2,a2)
{
(sapply(h,mod1,a1)*sapply(h,mod2,a2))
}


inchol=function(mod,a,x,y,n,eta){
G=NULL
P=1:n
i=1
diagG=rep(1,n)
crit=sum(diagG)
while(sum(diagG[i:n])>eta*n & i<n){
	G=cbind(G,rep(0,n))
	j=match(max(diagG[i:n]),diagG[i:n])+i-1
	P[c(i,j)]=P[c(j,i)]
	G[c(i,j),1:i]=G[c(j,i),1:i]
	G[i,i]=sqrt(diagG[j])
	if (i<n){
		Kcol=apply(sqrt(outer(x[P[(i+1):n]],x[P[i]],"-")^2+outer(y[P[(i+1):n]],y[P[i]],"-")^2),1,mod,a)
		if(i>1){
			G[(i+1):n,i]=1/G[i,i]*(Kcol-G[(i+1):n,1:(i-1)]%*%as.matrix(G[i,1:(i-1)],i-1,1))
			diagG[(i+1):n]=1-apply(matrix(G[(i+1):n,1:i]^2,(n-i),i),1,sum)
		}		
		else{
			G[(i+1):n,i]=1/G[i,i]*Kcol
			diagG[(i+1):n]=1-G[(i+1):n]^2
		}
	}
	crit=c(crit,sum(diagG[i:n]))
	i=i+1
}
if(i==n){
	G=cbind(G,rep(0,n))
	G[n,n]=sqrt(diagG[n])
	i=i+1
	crit=c(crit,0)
}	
return(list(G[order(P),],P,i-1,crit))
}


covspam=function(x,y,mod,a){
A=spam(0,n,n)
mh=NULL
for (i in 1:(n-1)){
	h=sqrt(outer(x[(i+1):n],x[i],"-")^2+outer(y[(i+1):n],y[i],"-")^2)
	mh=rbind(mh,c(i,min(h)))
	id=which(h<a)
	if(length(id)>0) A[i,((i+1):n)[id]]=sapply(h[id],mod,a)
}
(A+t(A)+diag.spam(n))
}

covtap=function(x,y,mod,a,mod1,eta){
A=spam(0,n,n)
mh=NULL
for (i in 1:(n-1)){
	h=sqrt(outer(x[(i+1):n],x[i],"-")^2+outer(y[(i+1):n],y[i],"-")^2)
	mh=rbind(mh,c(i,min(h)))
	id=which(h<a*eta)
	if(length(id)>0) A[i,((i+1):n)[id]]=sapply(h[id],tap,mod,a,mod1,eta*a)
}
(A+t(A)+diag.spam(n))
}

varioh=function(Z,hmax){
	d=dim(Z)
	varioh=0
	for(i in 1:hmax){
		varioh=c(varioh,1/2*mean((Z[,1:(d[2]-i)]-Z[,(i+1):d[2]])^2))
	}
}

object.sizes <- function()
 {
     return(rev(sort(sapply(ls(envir=.GlobalEnv), function (object.name)
         object.size(get(object.name))))))
 }
