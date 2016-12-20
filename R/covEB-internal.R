.AicW<-function(S,z){
	np<-z[lower.tri(z)]
	p<-sum(np!=0)
	v<-ncol(z)
	aic<-2*p-2*dwishart(S,v,z,TRUE)
	return(aic)
}

.AicN<-function(x,z){
	np<-z[lower.tri(z)]
	p<-sum(np!=0)
	#assume x is rows are observations, columns are variables
	muest<-colMeans(x)
	aic<-2*p-sum(2*dmvnorm(x,muest,z,TRUE))
	return(aic)
}


.SelThresh<-function(S,data=NULL,dist=c("N","IW")){

	if(!is.positive.definite(S)){
		S<-nearPD(S)$mat
		S<-as.matrix(S)
	}
	dist<-match.arg(dist)
	Sc<-cov2cor(S)
	Sc<-as.matrix(Sc)
	sl<-Sc[lower.tri(Sc)]
	sl<-abs(sl)
	lamrange<-seq(median(sl),max(sl)-0.1,by=0.01)
	outm<-matrix(0,nrow=length(lamrange),ncol=2)
	outm[,1]<-lamrange
	for(i in 1:length(lamrange)){
		ztemp<-Sc
		ztemp[abs(ztemp)<lamrange[i]]<-0
		gz<-graph.adjacency(ztemp,mode="upper",weighted=TRUE)
		zclust<-clusters(gz)
		zm<-zclust$membership
		nc<-max(zm)
		z<-diag(0,nrow=nrow(S),ncol=ncol(S))
		for(j in 1:nc){
			sel<-which(zm==j)
			z[sel,sel]<-S[sel,sel]
		}
		if(dist=="N"){
			if(is.null(data)){
				stop('Data needs to be provided to use MVN likelihood')
			}
			aico<-.AicN(data,z)
		}else{
			aico<-.AicW(S,z)}
			
		outm[i,2]<-aico
		
	}
	wm<-min(outm[,2],na.rm=TRUE)
	swm<-which(outm[,2]==wm)
	
	return(outm[swm[1],1])	
}
.EBWishsingle <-
function(S,z,gamma,n,happrox){
#S needs to be the covariance matrix not the correlation matrix.

r=cov2cor(S)
r=r[lower.tri(r)]

vmat<-matrix(diag(S),nrow=nrow(S),ncol=nrow(S))
vmat2<-matrix(diag(S),nrow=nrow(S),ncol=nrow(S),byrow=TRUE)

multiplier<-sqrt(vmat*vmat2)
zcov<-z*multiplier
diag(z)<-diag(S)
diag(zcov)<-diag(S)
#so now zcov is made up of diagonal sample variances and estimated correlations (converted to estimated covariances using the multiplier) These are fixed or mixed depending on what for z is passed from initial function. z is sample variances and gamma estimates for correlations.

#now move onto calculation of the d.f (lambda the second parameter of the IW):


#hpprox- hypergeometric approximations:
if(happrox){
	rsq<-1-((n-2)/(n-1))*(1-r^2)*hyperg_2F1(1,1,(n+1)/2,(1-r^2))
}else{
	rsq<-r^2
}
rsqm<-matrix(0,nrow=nrow(S),ncol=nrow(S))
rsqm[lower.tri(rsqm)]<-rsq
rst<-t(rsqm)
rsqm[upper.tri(rsqm)]<-rst[upper.tri(rst)]

if(happrox){
	hg<-r*hyperg_2F1(0.5,0.5,(n-1)/2,1-r^2)
}else{
	#use sample correlations
	hg<-r
}
hgm<-matrix(0,nrow=nrow(S),ncol=nrow(S))
hgm[lower.tri(hgm)]<-hg
hgt<-t(hgm)
hgm[upper.tri(hgm)]<-hgt[upper.tri(hgt)]

#replaced gamma with z as no longer a flat prior
knmat<-(rsqm-2*hgm*z-z^2)
denom<-(1-z^2)^2
knoffdiag<-knmat[lower.tri(knmat)]


ksq<-sum(knoffdiag)/(sum(denom[lower.tri(denom)]))
ifelse(ksq<=0,lambda<-n,lambda<-(1/ksq)-3)
if(lambda<1){lambda=1}

unsmoothsigma<-(lambda*zcov+S*(n-1))/(lambda+n)

#the smoothing modification to the variances (diagonal entries)
z1<-matrix(diag(zcov),nrow=nrow(S),ncol=ncol(S))
z2<-matrix(diag(zcov),nrow=nrow(S),ncol=ncol(S),byrow=TRUE)
multiplier2<-sqrt(vmat*vmat2/z1*z2)
smoothsigma<-unsmoothsigma
diag(smoothsigma)<-diag(unsmoothsigma)*diag(multiplier2)

return(list(smoothsigma=smoothsigma,unsmoothsigma=unsmoothsigma,ksq=ksq,lambda=lambda,zcov=zcov))
}

.EBWishartc <-
function(S,cutoff=0.5,n){
	
#n=ncol(S)
r=cov2cor(S)
r=r[lower.tri(S)]
#inc<-r!=0
#change 2/6/15:
#hg<-r*hyperg_2F1(0.5,0.5,(n-1)/2,(1-r^2))
#hgm<-matrix(hg,nrow=nrow(r))
#gamma<-mean(hgm[lower.tri(hgm)])

#r<-abs(r)
#change 2/6/15:
hg<-r
#hg<-r*hyperg_2F1(0.5,0.5,(n-1)/2,(1-r^2))
gamma<-mean(hg)
#for a flat prior, with sample variances on the diagonal and all off diagonal elements equal:
z<-matrix(gamma,nrow=nrow(S),ncol=nrow(S))
diag(z)<-1
vmat<-matrix(diag(S),nrow=nrow(S),ncol=nrow(S))
vmat2<-matrix(diag(S),nrow=nrow(S),ncol=nrow(S),byrow=TRUE)
multiplier<-sqrt(vmat*vmat2)
z<-z*multiplier
diag(z)<-diag(S)
#so now z is made up of diagonal sample variances and estimated correlations (converted to estimated covariances using the multiplier)
#check with the cutoff for z independent or not:
if(gamma<cutoff){
	z<-matrix(0,nrow=nrow(S),ncol=nrow(S))
	diag(z)<-diag(S)
	gamma<-0
}
#now move onto calculation of the d.f (lambda the second parameter of the IW):
#change 2/6/15:
rsq<-r^2
#rsq<-1-((n-2)/(n-1))*(1-r^2)*hyperg_2F1(1,1,(n+1)/2,(1-r^2))

test<-hyperg_2F1(1,1,(n+1)/2,(1-r^2))
test2<-1-r^2
#rsqm<-matrix(rsq,nrow=nrow(r))
#knmat<-rsqm-2*hgm*gamma-gamma^2
knma<-rsq-2*hg*gamma-gamma^2
#knoffdiag<-knmat[lower.tri(knmat)]
#ksq<-sum(2*knoffdiag)/n*(n-1)*((1-gamma^2)^2)
ksq<-sum(2*knma/n*(n-1)*((1-gamma^2)^2))

ifelse(ksq<=0,lambda<-n,lambda<-(1/ksq)-3)
if(lambda<1){lambda=1}

unsmoothsigma<-(lambda*z+S*(n-1))/(lambda+n)

#the smoothing modification to the variances (diagonal entries)
z1<-matrix(diag(z),nrow=nrow(S),ncol=ncol(S))
z2<-matrix(diag(z),nrow=nrow(S),ncol=ncol(S),byrow=TRUE)
multiplier2<-sqrt(vmat*vmat2/z1*z2)
smoothsigma<-unsmoothsigma
diag(smoothsigma)<-diag(unsmoothsigma)*diag(multiplier2)
return(list(smoothsigma=unsmoothsigma,unsmoothsigma=unsmoothsigma,ksq=ksq,lambda=lambda))
}
.EBEMWishart <-
function(S,n){
#S needs to be the covariance matrix not the correlation matrix.
#print(S)	
#n=ncol(S)
n=n
k<-S[lower.tri(S)]
r=cov2cor(S)
r=r[lower.tri(S)]
#change 26/5/15:
#r<-abs(r)
#change 26/5/15:
hg<-r
#hg<-r*hyperg_2F1(0.5,0.5,(n-1)/2,(1-r^2))

hgm<-matrix(0,nrow=nrow(S),ncol=nrow(S))
hgm[lower.tri(hgm)]<-hg
hgt<-t(hgm)
hgm[upper.tri(hgm)]<-hgt[upper.tri(hgt)]


#use as initial estimate of gamma
gamma<-mean(hg[abs(k)>0.001])
#gamma2<-mean(hg[k==0.001])

#oldll<-0
#newll<-100

#p1<-sum(k>0.001)/length(k)
#p2<-1-p1
#Need to create a mixture to give a weighted mean of the original individual gamma values. 
#while(abs(oldll-newll)>0.001){
#oldll<-newll
#nonzero<-(r-gamma)*(sqrt((l-2)/(1-(r-gamma)^2)))


#tvalsnz<-dt(nonzero,l-2)

#testvals<-(r-gamma2)*(sqrt((l-2)/(1-(r-gamma2)^2)))
#tvalzero<-dt(testvals,l-2)

#T1<-p1*tvalsnz/(p1*tvalsnz+p2*tvalzero)
#T2<-p2*tvalzero/(p1*tvalsnz+p2*tvalzero)

#p1<-mean(T1)
#p2<-1-p1

#gamma<-sum(T1*hg)/sum(T1)
#gamma2<-sum(T2*hg)/sum(T2)
#newll<-sum(log(T1*tvalsnz)+log(T2*tvalzero))
#}

#for a flat prior, with sample variances on the diagonal and all off diagonal elements equal:
z1<-matrix(gamma,nrow=nrow(S),ncol=nrow(S))

#IDnonzero<-T1>0.5

IDnonzero<-abs(k)>0.001
zeromat<-matrix(0,nrow=nrow(S),ncol=nrow(S))
zeromat[lower.tri(zeromat)]<-IDnonzero
zmt<-t(zeromat)
zeromat[upper.tri(zeromat)]<-zmt[upper.tri(zmt)]
z<-matrix(0,nrow=nrow(S),ncol=nrow(S))
z[as.logical(zeromat)]<-gamma

#diag(z)<-1
vmat<-matrix(diag(S),nrow=nrow(S),ncol=nrow(S))
vmat2<-matrix(diag(S),nrow=nrow(S),ncol=nrow(S),byrow=TRUE)

multiplier<-sqrt(vmat*vmat2)
z<-z*multiplier
diag(z)<-diag(S)
#so now z is made up of diagonal sample variances and estimated correlations (converted to estimated covariances using the multiplier)
#want to create a mixture value for gamma depending on whether more likely to be non-zero or not. So a gamma value close to zero will get closer to zero. One larger will still be large.
#first calc non-zero values for testing:
#testvals<-S[lower.tri(S)]

#now move onto calculation of the d.f (lambda the second parameter of the IW):
#r<-rorig
#change 26/5/15:
rsq<-r^2
#rsq<-1-((n-2)/(n-1))*(1-r^2)*hyperg_2F1(1,1,(n+1)/2,(1-r^2))
test<-hyperg_2F1(1,1,(n+1)/2,(1-r^2))
test2<-1-r^2
rsqm<-matrix(0,nrow=nrow(S),ncol=nrow(S))
rsqm[lower.tri(rsqm)]<-rsq
rst<-t(rsqm)
rsqm[upper.tri(rsqm)]<-rst[upper.tri(rst)]
#replaced gamma with z as no longer a flat prior
knmat<-rsqm-2*hgm*z-z^2
#knmat<-rsq-2*hg*gamma-gamma^2
knoffdiag<-knmat[lower.tri(knmat)]
knoffdiag<-knoffdiag
ksq<-sum(2*knoffdiag)/n*(n-1)*((1-gamma^2)^2)
#ksq<-sum(2*knmat)/n*(n-1)*((1-gamma^2)^2)
ifelse(ksq<=0,lambda<-n,lambda<-(1/ksq)-3)
if(lambda<1){lambda=1}

unsmoothsigma<-(lambda*z+S*(n-1))/(lambda+n)

#the smoothing modification to the variances (diagonal entries)
z1<-matrix(diag(z),nrow=nrow(S),ncol=ncol(S))
z2<-matrix(diag(z),nrow=nrow(S),ncol=ncol(S),byrow=TRUE)
multiplier2<-sqrt(vmat*vmat2/z1*z2)
smoothsigma<-unsmoothsigma
diag(smoothsigma)<-diag(unsmoothsigma)*diag(multiplier2)
return(list(smoothsigma=unsmoothsigma,unsmoothsigma=unsmoothsigma,ksq=ksq,lambda=lambda))
}
