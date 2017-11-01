#Equations to compute mean and variance of the L statistics
#
# Lee (2004). A generalized significance testing method for global measures
# of spatial association: an extension of the Mantel test. 
# Environment and Planning A 2004, volume 36, pages 1687 - 1703
#

Pmatrix<-function(W)
{
	n<-nrow(W)
	II<-matrix(rep(1, n), ncol=1)
	P<-t(W)%*%W
	P<-P/((t(II)%*%(P%*%II))[1,1])

	return(P)
}

Qmatrix<-function(x,y)
{
	n<-length(x)
	Q<- outer(x-mean(x), y-mean(y)) /(((n-1)/n)*sd(x)*sd(y))

	return(Q)
}


#Some quantities on the matrices, from Lee (2004), table 1, page 1690
PQquant<-function(M)
{
        II<-matrix(rep(1, nrow(M)), ncol=1)

        aux<-list(
                Foff0=as.numeric(t(II)%*%(M%*%II)-sum(diag(M))),
                Fon0=sum(diag(M)),
                Foff1=sum(diag(t(M)%*%M))-sum(diag(M)^2),
                Fon1=sum(diag(M)^2),
                Foff2 = sum((M%*%II-diag(M))^2),
                Fall2 = as.numeric(sum((M%*%II)^2)) #as.numeric(t(II)%*%(t(M)%*%M)%*%II)
        )

        return(aux)
}




#Values used to compute the expectation
EGamma<-function(x, y, W, P=NULL, Q=NULL)
{
	n<-length(x)

	if(is.null(P))
		P<-PQquant(Pmatrix(W))
	if(is.null(Q))
		Q<-PQquant(Qmatrix(x,y))

	xx<-list(
		EGammaoff = P$Foff0*Q$Foff0/(n*(n-1)),
		EGammaon  = P$Fon0*Q$Fon0/n
	)

	return(xx)
}

#Values used to compute the variance
VarGamma<-function(x,y,W, P=NULL, Q=NULL)
{
	n<-length(x)

	if(is.null(P))
		P<-PQquant(Pmatrix(W))
	
	if(is.null(Q))
		Q<-PQquant(Qmatrix(x,y))

	EG<-EGamma(x,y,W, P=P, Q=Q)
	


#varGammaoff
varGammaoff<-2*P$Foff1*Q$Foff1/(n*(n-1))
varGammaoff<-varGammaoff+4*(P$Foff2-P$Foff1)*(Q$Foff2-Q$Foff1)/(n*(n-1)*(n-2))

varGammaoff<-varGammaoff+((P$Foff0^2+2*P$Foff1-4*P$Foff2)*(Q$Foff0^2+2*Q$Foff1-4*Q$Foff2))/(n*(n-1)*(n-2)*(n-3))

varGammaoff<-varGammaoff-EG$EGammaoff^2

#VarGammaon
varGammaon<- P$Fon1*Q$Fon1/n

varGammaon<- varGammaon + (P$Fon0^2-P$Fon1)*(Q$Fon0^2-Q$Fon1)/(n*(n-1))

varGammaon<- varGammaon - EG$EGammaon^2


#Covarianza gammaon-off
varGammaonoff<- (P$Fall2-P$Fon1-P$Foff2)*(Q$Fall2-Q$Fon1-Q$Foff2)/(2*n*(n-1))

varGammaonoff<-varGammaonoff+(P$Fon0*P$Foff0-(P$Fall2-P$Fon1-P$Foff2))*( Q$Fon0*Q$Foff0-(Q$Fall2-Q$Fon1-Q$Foff2)  )/(n*(n-1)*(n-2))

varGammaonoff<-varGammaonoff - EG$EGammaoff*EG$EGammaon


	return(list(varGammaon=varGammaon,varGammaoff=varGammaoff,
		varGammaonoff=varGammaonoff))

}

