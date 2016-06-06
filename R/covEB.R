covEB <-
function(Covmat,delta=0.1,shift=0.05,cutoff=NULL,startlambda=0.5){
	#check data and parameter inputs are in valid range:
	eigen<-eigen(Covmat)$values
	check<-all(eigen>0)
	if(!check){
		cat('Warning: Covariance matrix is not positive semi definite')
	}
	if(delta<0|delta>1){
		stop('Delta must be between 0 and 1')
	}
	if(shift<0|shift>1){
		stop('Shift parameter must be between 0 and 1')
	}
	if(startlambda>1){
		stop('Starting lambda value must be less than 1')
	}
	origmat<-Covmat
	Cormat<-cov2cor(Covmat)
	Cormat<-abs(Cormat)
	testvals<-Cormat[lower.tri(Cormat)]
	#minl<-min(testvals)
	maxl<-max(testvals)
	#s=minl
	s=startlambda
	covlist<-list()
	i=1

	resultsmat<-matrix(0,nrow=nrow(Cormat),ncol=ncol(Cormat))
	
	countmat<-matrix(0,nrow=nrow(Cormat),ncol=ncol(Cormat))
	while(s<(maxl-delta)){
		lseq<-seq(from=s,to=maxl,by=delta)
		s=s+shift
		if(length(lseq)>1){
		for(i in 1:(length(lseq)-1)){
			temp<-Cormat
			
			temp[temp>lseq[i+1]]<-0
			temp[temp<lseq[i]]<-0
			diag(temp)<-0
			tempclust<-clusters(graph.adjacency(temp,mode="upper",weighted=TRUE))
			mem<-tempclust$membership
			nocl<-max(mem)
			reslist<-list(length=nocl)
			uncon<-c()
			for(i in 1:nocl){
				w<-which(mem==i)
				#m<-Covmat[w,w]
				check<-temp[w,w]
				
				
				if(length(w)==1){uncon=c(uncon,w)}else{
					
					m<-Covmat[w,w]
					m[check==0]<-0.001
					
					diag(m)<-diag(Covmat[w,w])
					
					}
				if(sum(abs(check))>0&&length(w)>1){
					if(is.null(cutoff)){
						reslist[[i]]<-.EBEMWishart(m)
					}else{
						reslist[[i]]<-.EBWishartc(m,cutoff=cutoff)
					}
					outmat<-cov2cor(reslist[[i]]$smoothsigma)
					outmat[check==0]<-0
					#replace w,w with check!=0
					
					resultsmat[w,w]<-resultsmat[w,w]+outmat
					m[outmat!=0]<-1
					m[outmat==0]<-0
					countmat[w,w]<-countmat[w,w]+m
				}	
			
				
			}
			
			#need to identify the off-block diagonal elements
			
			
		}
		i=i+1
		}else{
			s<-maxl
		}
	}
	
	
	finalmat<-resultsmat/countmat
	
	#need to check don't have any NaN entries, if so replace with sample values for the moment
		#finalmat is a correlation matrix (result from Wishart functions is a covariance matrix)
	finalmat[is.na(finalmat)]<-Cormat[is.na(finalmat)]
	sf<-sign(finalmat)
	sc<-sign(origmat)
	change<-sf==sc
	finalmat[!change]<-(-1*finalmat[!change])
	rownames(finalmat)<-rownames(Covmat)
	colnames(finalmat)<-colnames(Covmat)
	return(finalmat)
	
}
