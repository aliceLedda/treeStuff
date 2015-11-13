
library(ape)

source("progs/coalescentExpFunc1.R")
source("progs/coalescentNonExpFunc1.R")
source("progs/coalescentBothFuncNested1.R")



############# case 1. : Non expanding population Only ########################################### 
NofSim<-100
testExp<-c(0.001,1,2,3,4,5,10) 		#give the exponential rate to test against
times0<-		#give the times at which the leaves appear

times0<-rep(0,10)+0.3*c(1:10)
theta<-1

AnalysisNon<-function(testExp,times0,NofSim,theta){
Results<-matrix(ncol=4,nrow=NofSim,data=0)
for (i in 1:NofSim){
			
	tree<-MakeTreeNon(timesNon=times0)
	Like<-c()

	LikeNon<-log(LikeliHood(tree,theta))
	Like<-LikeNon
	for (j in 1:length(testExp)){
		alpha<-testExp[j]
		LikeExp<-(log(LikeliHoodExp(tree,theta,alpha)))
		
		Like<-c(Like,LikeExp)
	}
	if(anyNA(Like)){
	i<-i-1}
	else{
	Max<-which(Like==max(Like))
	testExpTmp<-c(0,testExp)
	Ll<-(-2*(Like[Max]-Like[1]))
	Results[i,]<-c(testExpTmp[(Max)], Like[1], Like[Max],Ll)
	}
	}
return(Results)
}
	
############# case 2. : Expanding population Only ########################################### 
NofSim<-100
testExp<- 		#give the exonential rate to test against
times0<-		#give the times at which the leaves appear
alpha<-
AnalysisExp<-function(testExp,times0,alpha,NofSim,theta){
Results<-matrix(ncol=5,nrow=NofSim,data=0)
for (i in 1:NofSim){
	
	tree<-MakeTreeExp(alpha=alpha,timesExp=times0)
	Like<-c()

	LikeNon<-(log(LikeliHood(tree,theta)))
	Like<-LikeNon
	for (j in 1:length(testExp)){
		beta<-testExp[j]
		LikeExp<-(log(LikeliHoodExp(tree,theta,beta)))
		
		Like<-c(Like,LikeExp)
	}
	if(anyNA(Like)){
	i<-i-1}
	else{
	Max<-which(Like==max(Like))
	testExpTmp<-c(0,testExp)
	Alpha<-which(testExpTmp==alpha)
	Ll<-(-2*(Like[Max]-Like[Alpha]))
	Results[i,]<-c(alpha,testExpTmp[Max], Like[Alpha], Like[Max],Ll)
	}
	}
return(Results)
}

