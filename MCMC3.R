library(ape)
library(phangorn)

source("coalescentExpFuncB1.R")
source("coalescentNonExpFuncB1.R")
source("coalescentBothFuncNested1.R")
source("AnalysisFunct.R")
source("like0.3.R")


############# function usefulNodes ##################################
##########we will use only the nodes that have at least 4 (or maxChild) children ##############
########## if the expansion happens in nodes that have less than 4 children is very difficult to detect  
####### but the nodes that have fewer children are more numerous in some kinds of trees
###### (in the caterpillar tree it is not true)

usefulNodes<-function(tree, maxChilds){
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

############ function LikelihoodBoth  ##############
######## helps to build the MCMC ####################
LikelihoodBoth<-function(tree,theta,alpha, branch){
	ETree<-ExpTree(tree, branch)
	NTree<-NonExpTree(tree,branch)
	Lnon<-LikeliHood(NTree, theta)
	Lexp<-LikeliHoodExp(ETree,theta,alpha)
	LBoth<-log(Lnon)+log(Lexp)
	return(list(likelihood=LBoth,rate=alpha,branch=branch,Lnon=log(Lnon)))
	}


########## function checkNonexp #############
########## checks that the tree is not just a normal costant pop coalescet tree ###########

checkNonExp<-function (tree,tot){
	first<-c()
	for (i in 1:tot){
	tmp<-LikeliHood(tree, theta)
	first<-c(first, tmp)
	}
	med<-median(first)
	return(list(median=med,all=first))
}

########### function extractNode  ############
extractNodes<-function(nodes,branch){
	l<-length(nodes)
	me<-which(nodes==branch)
	use<-4
	if (((me-use)>0)&((me+use)<=l)){
	subs<-nodes[(me-use):(me+use)]}
	else{
		if(((me-use)<=0)){
		missing=me-use
		subs<-c(nodes[1:(me+use)],nodes[(l+missing):l])
		}
		if((me+use)>l){
		missing=me+use-l
		subs<-c(nodes[(me-use):l],nodes[1:missing])
		}
	}
	one<- sample(1:((2*use)+1),1)  
	node<-subs[one]
	return(node)
}

########### function extractNode  ############
extractNode1<-function(nodes,branch){
	l<-length(nodes)
	one<- sample(1:l,1)  
	node<-nodes[one]
	return(node)
}


##############  function accept ###############

accept<-function(this){
	return(list(likelihood=this$likelihood,rate=this$rate, branch=this$branch, Lnon=this$Lnon))
}

############ function RandomAccept #################
RandomAccept<-function(one,two){
	if (one$likelihood<0.000000000001){
		one$likelihood<-(-40)}
	if  (two$likelihood<0.000000000001){
		two$likelihood<-(-40.1)}
	if (exp(one$likelihood)>exp(two$likelihood)){
		p<-exp(two$likelihood)/exp(one$likelihood)
		}
	else{
		p<-exp(one$likelihood)/exp(two$likelihood)
		}
	q<-runif(1, min=0,max=1)
	if(p>q){
		accepted<-accept(one)
		}
	else{
		accepted<-accept(two)
		}
	return(accepted)
}

###############  MCMC  #######################
######### input: a tree (phylo object) ###############
######### output: a table where we have two columns
#######      "    a branch where we think the expansion started
########     "    and a value for the expected expansion rate

Posterior<-function (tree) {
	tot<-100
	first<-checkNonExp(tree,tot)
	if (first$median<0.01){ 
		table<-first$all
	}
	else{
		table<-MCMC(tree,expRates,maxChilds)
	}
	return(table)
}
	
core<-function(ETree, NTree,theta, expRates, branch){
	allLexp<-c()
	for (j in 1:length(expRates)){
		alpha<-expRates[j]			
		Lexp<-LikeliHoodExp(ETree,theta,alpha)
		allLexp<-c(allLexp, Lexp)
	}
	like<-max(allLexp)
	Erate<-expRates[which(allLexp==like)]
	Lnon<-LikeliHood(NTree,theta)
	li<-log(like)+log(Lnon)
	return(list(likelihood=li,rate=Erate,branch=branch, Lnon=log(Lnon)))	
}


MCMC<-function(tree, expRates, maxChilds){
	######## if expRates is not given, set a default
	####### if maxChilds is not given, set a default
	nRec<-10000
	alpha<-20
	nodes<-usefulNodes(tree,maxChilds)
	branch<-extractNode1(nodes)
	ETree<-ExpTree(tree, branch)
	
	write.tree(ETree, file="tmp2.tre")
	ETree<-read.tree("tmp2.tre")

	NTree<-NonExpTree(tree,branch)
	
	write.tree(NTree, file="tmp1.tre")
	NTree<-read.tree("tmp1.tre")

	Lnon<-LikeliHood(NTree, theta)
	
	best<-LikelihoodBoth(ETree,NTree, theta, alpha,branch)
	table<-matrix(ncol=7, nrow=nRec)
	i<-1
	
	while (i<nRec){
		oldbranch<-branch
		branch<-extractNodes(nodes,oldbranch)
		ETree<-ExpTree(tree, branch)
		write.tree(ETree, file="tmp2.tre")
		ETree<-read.tree("tmp2.tre")

		NTree<-NonExpTree(tree,branch)
		write.tree(NTree, file="tmp1.tre")
		NTree<-read.tree("tmp1.tre")

		
		Lnon<-LikeliHood(NTree, theta)
		table[i,5]<-best$rate
		table[i,6]<-best$branch
		table[i,7]<-best$likelihood
		tocompare<-LikelihoodBoth(ETree,NTree, theta, alpha, branch)
		#if ((tocompare$likelihood<=0)&(tocompare$likelihood>0.000000001)){
		if (tocompare$likelihood>best$likelihood){
				best<-accept(tocompare)
			}else{
			best<-RandomAccept(best,tocompare)
		}
		if((length(best$rate)==1)&&(is.numeric(best$rate))){
		table[i,1]<-best$likelihood
		table[i,2]<-best$Lnon ########## no devo inserire Lnon nella list
		table[i,3]<-best$rate
		table[i,4]<-best$branch ########### no devo inserirlo nel list che va in giro
		i<-i+1
		}
		#}
		#else{i<-i-1}
	}
	
	return(table)
}
	

	
	
	
