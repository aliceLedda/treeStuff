library(ape)

#Ttree<-tree

#to find the branch from the tree image use command trex
############## You need to specify a tree and a branch where you think the expansion started

########################  FUNCTIONS  #####################################################

##ExpTree extract the expanding tree given the total tree and the branch
#parameters: the total tree and the branch where we want to cut the subtree (from that branch to leaves)
#returns: the subtree from that branch on

ExpTree<-function(tree,branch){

	Etree<-extract.clade(tree,branch)
	return (Etree)
	}

########## NonExpTree exract the non expanding tree given the total tree and the branch 
#parameters: the total tree and the branch where we want to cut the subtree (excluding the tree from that branch to leaves)
#returns: the subtree excluding all the descendants of the branch

NonExpTree<-function(tree, branch){
	Etree<-extract.clade(tree,branch)

	Ntree<-drop.tip(tree, Etree$tip.label,subtree=TRUE,trim.internal=TRUE) 
	
	return(Ntree)
	}
	
######### PrepareTable extracts from the tree in phylo format the informations that we need to compute the likelihood
#preparatory function to get a table to compute the likelihood
#parameters: the tree (phylo object)
#returns: a table to be used to compute the likelihood.
#	the table contain 4 columns and as many rows as events (a node or a leaf).
#	the columns are: 1. the number of the event, 2. the time at which the event happens, time goes backwards, with 0 at present
#	3.the number of lineages present between the event and the next event (from root to leaves), 
#	4. the event type: r (root), l (leaf apperance), c (coalescente event)

PrepareTable<-function(tree){
	
	attr(tree,"order")<-"postorder"

	tab<-RootDist(tree)
	SortedTab<-tab[order(tab[,2]),]
	SortedTab[1,3]<-2
	LineagesTab<-NofLin(SortedTab)
	LineagesTab[,2]<-(-(as.numeric(LineagesTab[,2])-as.numeric(LineagesTab[length(LineagesTab[,2]),2])))
	return(LineagesTab)
}	
	
########## RootDist transforms the $edge matrix in a matrix displaying each node and its distance from the root
######### and its status as leaf (l) o coalescente event (c)
#parameter: tree
#returns: a 4 columns matrix. The columns are the same as the previous but the 3rd column is empty (will be filled by the next function)

RootDist<-function(tree){
	tab<-matrix(ncol=4, nrow=(nrow(tree$edge)+1))
	for(i in 1:length(tree$edge[,1])){
		tab[i,1]<-tree$edge[i,2]
		tab[i,2]<-node.depth.edgelength(tree)[tree$edge[i,2]]
		if (tab[i,1] %in% tree$tip){tab[i,4]<-"l"
		}else{tab[i,4]<-"c"}
	}
	root<-setdiff(unique(tree$edge[,1]),unique(tree$edge[,2]))
	i<-nrow(tree$edge)+1
	tab[i,1]<-root

	tab[i,2]<-node.depth.edgelength(tree)[root]
	tab[i,4]<-"r"
	return(tab)
}

######### NofLin adds the third column of the matrix in which we specify how many lineages are present in the interval i
###### the intervals are defined in the reverse direction with respect to the direction of time
######## ( i.e. the first interval is between the root and the first event beckwards in time and contains 2 lineages)
####### thus every coalescent event ADDS one lineage and every leaf event counts -1 lineages
#parameter: sorted tab from the PrepareTable function
#returns: fills the 3rd column of the matrix

NofLin<-function (SortedTab){
	#if ((SortedTab[1,4]=="r")&(SortedTab[2:length(SortedTab[,1]),4]==("c"|"l")))
	for (i in 2:length(SortedTab[,1])){
		if (SortedTab[i,4]=="c"){
			SortedTab[i,3]=as.numeric(SortedTab[(i-1),3])+1
		}
		else if(SortedTab[i,4]=="l"){
			SortedTab[i,3]=as.numeric(SortedTab[(i-1),3])-1
		}
	}
	return(SortedTab)
}

########### Likelihood Non Expanding Population  #############3
####### implements the product in formula 1 Drummond et al 2005 ##############3
###### can have as an input both the tree or the Matrix that comes from PrepareTable
# theta is the factor converting "coalescent time" into calendar time. 
#In our simulated trees is 1 because time is already in "coalescent time"
#### returns the likelihood of the tree under the no expansion assumption

LikeliHood<-function (input, theta){
	if (class(input)=="phylo"){
	table<-PrepareTable(input)
	}else{
	table<-input}
	Likeli<-1
	for (j in 1:(length(table[,1])-1)){
		i<-length(table[,1])-j
		ToM<-ToMultiply(table,i,theta)
		Likeli<-Likeli*ToM
		}
	Ninternal<-length(which(table[,4]!="l"))
	Likelihood<-(theta^(Ninternal))*Likeli
	return(Likelihood)
}
	
############  ToMultiply  ######################
########## does the calculation for the previous function ##################
#########  needs the table and the line that we are analising  ###############
# theta is the factor converting "coalescent time" into calendar time. 
#In our simulated trees is 1 because time is already in "coalescent time"
###### returns the contribution of that event to the overall likelihood 

ToMultiply<-function(table,i,theta){
	ki<-as.numeric(table[(i+1),3])
	ti<-as.numeric(table[i,2])
	timinus1<-as.numeric(table[(i+1),2])
	Expon<-(-ki*(ki-1)/(2*theta))*(ti-timinus1)
	ToM<-exp(Expon)
	return(ToM)
	}

########### Likelihood Expanding Population  #############3
####### implements the product in formula 8 Drummond et al 2005 ##############3
###### can have as an input both the tree or the Matrix that comes from PrepareTable 
##### and the rate alpha at which the population expands 
# theta is the factor converting "coalescent time" into calendar time. 
#In our simulated trees is 1 because time is already in "coalescent time"
#### returns the likelihood of the tree under the  assumption that the population is expanding at a rate alpha

LikeliHoodExp<-function (input, theta, alpha){
	if (class(input)=="phylo"){
	table<-PrepareTable(input)
	}else{
	table<-input}
	Likeli<-1
	for (j in 1:(length(table[,1])-1)){
		#j<-j+1
		i<-length(table[,1])-j
		ToM<-ToMultiplyExp(table, i, theta, alpha)
		Likeli<-Likeli*ToM
		ToM
		Likeli
		j
		}
	Ninternal<-length(which(table[,4]!="l"))
	Likelihood<-(theta^(Ninternal))*Likeli
	return(Likelihood)
}
	
############  ToMultiply Expansion  ######################
########## does the calculation for the previous function ##################
#########  needs the table and the line that we are analising  ###############
##### and the rate alpha at which the population expands 
# theta is the factor converting "coalescent time" into calendar time. 
#In our simulated trees is 1 because time is already in "coalescent time"
###### returns the contribution of that event to the overall likelihood 


ToMultiplyExp<-function(table,i,theta, alpha){
	ki<-as.numeric(table[(i+1),3])
	ti<-as.numeric(table[i,2])
	timinus1<-as.numeric(table[(i+1),2])
	eti<-exp(alpha*ti)
	etiminus1<-exp(alpha*timinus1)
	Expon<-(-ki*(ki-1)/(2*theta*alpha))*(eti-etiminus1)
	ToM<-eti*exp(Expon)
	return(ToM)
	}

################## Likelihood Both #############
######### computes the likelihood of a tree with an expanding subtree ######
######## needs the tree (in phylo format) 
######### and the number of the branch where we expect the expansion to have started 
######## and the expansion rate (expected?)
# theta is the factor converting "coalescent time" into calendar time. 
#In our simulated trees is 1 because time is already in "coalescent time"
#### returns the likelihood of the tree under the  assumption that 
### a specific part of the tree is expanding at a rate alpha


LikelihoodBoth<-function(tree, branch, alpha,theta){
	ETree<-ExpTree(tree, branch)
	NTree<-NonExpTree(tree,branch)
	Lnon<-LikeliHood(NTree, theta)
	Lexp<-LikeliHoodExp(ETree,theta,alpha)
	LBoth<-Lnon*Lexp
	return(LBoth)
	}
