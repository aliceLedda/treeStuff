library(ape)
library(phangorn)
source("NonExpFunc1.R")
source("ExpFunc2.R")
source("coalescentBothFuncNested1.R")
#source("progs/AnalysisFunct.R")
source("like0.5Sum.R")

############# function usefulNodes ##################################
##########we will use only the nodes that have at least 4 (or maxChild) children ##############
########## if the expansion happens in nodes that have less than 4 children is very difficult to detect  
####### but the nodes that have fewer children are more numerous in some kinds of trees
###### (in the caterpillar tree it is not true)

AvailableNodes<-function(tree, maxChilds){
if (is.na(maxChilds)){maxChilds=4}
allnodes<-unique(tree$edge[,1])
nodes<-c()
for (i in 1:length(allnodes)){
	if ((length(unlist(Descendants(tree,allnodes[i])))>maxChilds)&&(length(unlist(Descendants(tree,allnodes[i])))<length(tree$tip.label))){
		nodes<-c(nodes,allnodes[i])
		}
	}
return(nodes)
}

########### function extractNode  ############
ChooseNode<-function(nodes,branch){
	l<-length(nodes)
	one<- sample(1:l,1)  
	node<-nodes[one]
	return(node)
}

########### function ChooseRate ############
ChooseRate<-function(alpha){
	if(sqrt(alpha)>5){
	newrate<-rnorm(1,mean=alpha,sd=sqrt(alpha))}
	else{newrate<-abs(rnorm(1,mean=alpha,sd=5))}
	return(newrate)
}





########## LogLikeliHoodBoth #####################

LogLikeliHoodBothO<-function(ETree,NTree,theta,alpha, branch){
	Lnon<-LogLikeliHood(NTree, theta)
	Lexp<-LogLikeliHoodExp(ETree,theta,alpha)
	LBoth<-Lnon+Lexp
	return(list(likelihood=LBoth,rate=alpha,branch=branch))
	}

########## LogLikeliHoodBoth #####################

LogLikeliHoodBoth<-function(tree,theta,alpha, branch){
	ETree<-ExpTree(tree, branch)
	
	write.tree(ETree, file="tmp2.tre")
	ETree<-read.tree("tmp2.tre")

	NTree<-NonExpTree(tree,branch)
	
	write.tree(NTree, file="tmp1.tre")
	NTree<-read.tree("tmp1.tre")

	Lnon<-LogLikeliHood(NTree, theta)
	Lexp<-LogLikeliHoodExp(ETree,theta,alpha)
	LBoth<-Lnon+Lexp
	return(list(likelihood=LBoth,rate=alpha,branch=branch))
	}


############ NewState #########################


State<-function(tree, branch, theta, alpha){
    ETree<-ExpTree(tree, branch)
	write.tree(ETree, file="tmp2.tre")
	ETree<-read.tree("tmp2.tre")

	NTree<-NonExpTree(tree,branch)
	write.tree(NTree, file="tmp1.tre")
	NTree<-read.tree("tmp1.tre")

	LikeNew<-LogLikeliHoodBoth(tree,theta,alpha,branch)
	return(LikeNew)
}

########## Random Accept ############
RandomAccept<-function(Old, New){
	p<-New$likelihood/Old$likelihood
	q<-runif(1, min=0,max=1)
	if(p>q){
		Accepted<-New
	}
	else{
		Accepted<-Old
	}
	return(Accepted)
}
	

#############  MCMC  #################

MCMC<-function(tree,trials,burnon){
	nodes<-AvailableNodes(tree,4)
	branch<-ChooseNode(nodes)
	alpha<-ChooseRate(20)
	Oldstate<-State(tree, branch, theta, alpha)
	i<-1
	while(i<burnon){
		Newstate<-State(tree, branch,theta, alpha)
		if(Newstate$likelihood>Oldstate$likelihood){
			Oldstate<-Newstate
		}
		else{
			Oldstate<-RandomAccept(Oldstate,Newstate)
		}

		i<-i+1
	}
	check<-matrix(nrow=(trials-burnon),ncol=9)
	colnames(check)<-c("old.likelihood","old.branch","old.rate","proposed.likelihood","proposed.branch","proposed.rate","acccepted.likelihood","accepted.branch","accepted.rate")
	i<-1
	while(i<(trials-burnon)){
		branch<-ChooseNode(nodes)
		alpha<-ChooseRate(Oldstate$rate)
		Newstate<-State(tree, branch,theta, alpha)
		check[i,1]<-Newstate$likelihood
		check[i,2]<-Newstate$branch
		check[i,3]<-Newstate$rate
		check[i,4]<-Oldstate$likelihood
		check[i,5]<-Oldstate$branch
		check[i,6]<-Oldstate$rate

		if(is.na(Newstate$likelihood)){
			check[i,7]<-0
			check[i,8]<-0
			check[i,9]<-0
	
		}
		else{
		if(Newstate$likelihood>Oldstate$likelihood){
			Oldstate<-Newstate
		}
		else{
			Oldstate<-RandomAccept(Oldstate,Newstate)
		}
	check[i,7]<-Oldstate$likelihood
	check[i,8]<-Oldstate$branch
	check[i,9]<-Oldstate$rate
	}
	i<-i+1
	}
return(check)
}
