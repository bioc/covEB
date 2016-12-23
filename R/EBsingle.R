EBsingle<-function(Covmat,startlambda=0.5,n,happrox=FALSE,groups=NULL,aic=FALSE,data=NULL){
	if(aic&&is.null(data)){
		stop('Data needs to be provided to use AIC')
	}
	#check data and parameter inputs are in valid range:
	eigen<-eigen(Covmat)$values
	check<-all(eigen>0)
	#if(!check){
	#	stop('Covariance matrix is not positive definite')
	#}
	
	if(startlambda>1){
		stop('Starting lambda value must be less than 1')
	}
	origmat<-Covmat
	Cormat<-cov2cor(Covmat)
	Cormat<-abs(Cormat)

	if(is.null(colnames(Covmat))){
		cnames<-paste("V",1:ncol(Covmat),sep="")
		rownames(Covmat)<-cnames
		colnames(Covmat)<-cnames	
	}
	reslist<-list()
	
	if(aic){
		startlambda<-.selThresh(Covmat,data)
	}
	if(is.null(groups)){
			temp<-Cormat
			
			temp[abs(temp)<startlambda]<-0
			diag(temp)<-0
			tempclust<-clusters(graph.adjacency(temp,mode="upper",weighted=TRUE))
			mem<-tempclust$membership
			nocl<-max(mem)
			
			fullmat<-matrix(0,nrow=nrow(Covmat),ncol=ncol(Covmat))
			diag(fullmat)<-1
			rownames(fullmat)<-rownames(Covmat)
			colnames(fullmat)<-colnames(Covmat)
		
				z<-matrix(1,nrow=length(mem),ncol=length(mem))
				rownames(z)<-rownames(Covmat)
				colnames(z)<-colnames(Covmat)
					for(i in 1:nocl){
						w<-names(which(mem==i))
						
						#set entries outside block to zero:
						z[w,!(colnames(z)%in%w)]<-0
						z[!(colnames(z)%in%w),w]<-0
						
						#use as initial estimate of gamma
						ct<-Cormat[w,w]
						if(happrox){
							rijs<-ct[lower.tri(ct)]
							rhoijs<-rijs*hyperg_2F1(0.5,0.5,(n-1)/2,1-rijs^2)
							gamma<-mean(rhoijs)
						}else{
							
							gamma<-mean(ct[lower.tri(ct)])
						}
							
						z[w,w]<-gamma


						
					}
					
				}else{
					#user has passed a set of pre-defined groupings.
					if(!is.list(groups)){
						stop('Input groups needs to be a list of groupings')
					}
					nocl<-length(groups)
					allnames<-unlist(groups)
					if(length(allnames)!=nrow(Covmat)){
						stop('Please provide a group for each variable in Covmat. Dimensions do not match.')
					}
					fullmat<-matrix(0,nrow=nrow(Covmat),ncol=ncol(Covmat))
					diag(fullmat)<-1
					rownames(fullmat)<-allnames
					colnames(fullmat)<-allnames
		
				z<-matrix(1,nrow=length(allnames),ncol=length(allnames))
				rownames(z)<-allnames
				colnames(z)<-allnames
					for(i in 1:length(groups)){
						w<-groups[[i]]
						
						#set entries outside block to zero:
						z[w,!(colnames(z)%in%w)]<-0
						z[!(colnames(z)%in%w),w]<-0
						
						#use as initial estimate of gamma
						ct<-Cormat[w,w]
						if(happrox){
							rijs<-ct[lower.tri(ct)]
							rhoijs<-rijs*hyperg_2F1(0.5,0.5,(n-1)/2,1-rijs^2)
							gamma<-mean(rhoijs)
						}else{
							
							gamma<-mean(ct[lower.tri(ct)])
						}
							
						z[w,w]<-gamma


						
					}

				}
				
				
				#now run algorithm on all genes with block diagonal z prior:
				diag(z)<-1
				reslist<-.EBWishsingle(S=Covmat,z=z,gamma=gamma,n,happrox)
				finalmat<-cov2cor(reslist$unsmoothsigma)
			
		

#need to check don't have any NaN entries, if so replace with sample values for the moment
	#finalmat is a correlation matrix (result from Wishart functions is a covariance matrix)
	finalmat[is.na(finalmat)]<-Cormat[is.na(finalmat)]
	sf<-sign(finalmat)
	sc<-sign(origmat)
	change<-sf==sc
	finalmat[!change]<-(-1*finalmat[!change])
	return(finalmat)
	
}
