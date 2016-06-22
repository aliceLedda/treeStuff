library(ape)
#########FUNCTIONS TO GET PROGRAM WORKING

####  function Evolution -------------------------------------------------------
#starts the evolution matrix
#parameters: vector of extant leaves, vector of leaves timesExp

evolutionExpBoth<-function(lineagesExp, tipsExp, timesExp,treeHeightExp,alphaExp,PopR){
	#lineagesExp<-as.matrix(lineagesExp)
	Nl<-length(lineagesExp[,1])
	coalTimeExp<-extrTimeExpBoth(l=Nl,treeHeightExp=treeHeightExp,alphaExp=alphaExp,PopR=PopR)  
	coalTimeExp
	#while(((treeHeightExp+coalTimeExp)<timesExp[1])&(length(lineagesExp[,1])==1)){
	#	coalTimeExp<-extrTimeExp(Nl,treeHeightExp) 
	#}
	
	if(((treeHeightExp+coalTimeExp)<timesExp[1])|(identical(timesExp,numeric(0)))){
		coal<-chooseLeaves(length(lineagesExp[,1]))
		ma<-matrix(nrow=2,ncol=3,data=0)
		ma[,2]<-c(lineagesExp[coal[1],1],lineagesExp[coal[2],1])
		if(max(lineagesExp[,1])<(length(timesExp)+(max(tipsExp))+1)){
			ma[,1]<-rep((length(timesExp)+(max(tipsExp))+1),2)
			lineagesExp<-rbind(lineagesExp,c(0,0))
			lineagesExp[length(lineagesExp[,1]),1]<-(length(timesExp)+(max(tipsExp))+1)
			lineagesExp[length(lineagesExp[,1]),2]<-treeHeightExp+coalTimeExp
		}
		
		else{
			ma[,1]<-rep((max(lineagesExp[,1])+1),2)
			lineagesExp<-rbind(lineagesExp,c(0,0))
			lineagesExp[length(lineagesExp[,1]),1]<-(max(lineagesExp[,1])+1)
			lineagesExp[length(lineagesExp[,1]),2]<-treeHeightExp+coalTimeExp
		}
		ma[(ma[,2]==lineagesExp[coal[1],1]),3]<-(treeHeightExp+coalTimeExp-(lineagesExp[coal[1],2])) #this need to be fixed
		ma[(ma[,2]==lineagesExp[coal[2],1]),3]<-(treeHeightExp+coalTimeExp-(lineagesExp[coal[2],2]))
		lineagesExp<-lineagesExp[-(coal),,drop=FALSE]
		treeHeightExp<-treeHeightExp+coalTimeExp
		
		return(list(lineages=lineagesExp,tips=tipsExp,evo=ma,height=treeHeightExp,times=timesExp))
	}else{
		lineagesExp<-rbind(lineagesExp,c(0,0))
		lineagesExp[length(lineagesExp[,1]),1]<-(max(tipsExp)+1)
		lineagesExp[length(lineagesExp[,1]),2]<-timesExp[1]
		tipsExp<-c(tipsExp,(max(tipsExp)+1))
		treeHeightExp<-timesExp[1]
		timesExp<-timesExp[-1]
		return(list(lineages=lineagesExp,tips=tipsExp,evo=0,height=treeHeightExp,times=timesExp))
	}
	
}


#### function Choose Leaves --------------------------------------------------
#chooses 2 leaves out of the extant k to coalesce
#parameters: number of extant lineagesExp
chooseLeaves<-function(l){
	choosen<-c()
	cho<-sample(1:l,2)
	#choosen[1]<-leaves[cho[1]]
	#choosen[2]<-leaves[cho[2]]
	return(cho)
}

#### function Extract Time ---------------------------------------------------
#extracts "time" from an exponential distribution with rate k(k-1)/2
#parameters: the number of extant lineagesExp
extrTimeExpBoth<-function(l, treeHeightExp, alphaExp,PopR){
	if(l>1){
		par<-l*(l-1)/(2*PopR)
		timeO<-rexp(1,rate=par)
		if(alphaExp>0){
		timeNew<-((1/alphaExp)*log(1+(alphaExp*timeO)*exp(-(alphaExp*treeHeightExp))))
		}else{
			timeNew<-timeO
		}

		#timeNew<-abs(((1/alphaExp)*log(1+(alphaExp*timeO)*exp(-(alphaExp*treeHeightExp))))-treeHeightExp)
		#timeNew<-((1/alphaExp)*log(1+(alphaExp*timeO)*exp(-(alphaExp*treeHeightExp))))-treeHeightExp
	}else{timeNew<-100000}
	#if(timeNew>0)
	#{
		return(timeNew)
	#}
	#else{
	#	extrTimeExp(l,treeHeightExp,alphaExp)
	#}
}





#############function makephylo

# makephylo<-function(evoR){	
# 	#now that the matrix is ready I transform it in a phylo object. 
# 	#If you get seg fault, try plot.phylo(alb), if it gives a mistake the problem is that the matrix is not properly declared.
# 	tipsExp<-setdiff(unique(evoR[,2]),unique(evoR[,1]))
# 	#if(identical(evoR[1,],c(0,0,0),num.eq=TRUE)){
# 	evoR<-evoR[-1,]
# 	#}

# 	#tipsExp<-setdiff(unique(evo[,2]),unique(evo[,1]))

# 	evoR[,1]<-evoR[,1]+1

# 	evoR[!(evoR[,2] %in% tipsExp),2]<-evoR[!(evoR[,2] %in% tipsExp),2]+1

# 	evoR[evoR[,1]==max(evoR[,1]),1]<-max(tipsExp)+1

	
# 	tmp<-evoR[,1:2]
# 	tree<-list(edge=evoR[,1:2])
# 	class(tree)<-"phylo"
	
# 	attr(tree,"order")<-"postorder"
	
# 	tree$edge.length<-evoR[,3]
# 	tree$tip.label<-seq(1,length(tipsExp),1)
# 	tree$Nnode<-length(unique(evoR[,1]))
# 	#back<-list(tree=tree,cBR=cBR)
# 	return(tree)
# }
 
##########################################################
########### THE PROGRAM ##########################
##################################################

###idenfy the tipsExp....
######### the tipsExp are all the branches from which there is no other branch
#tipsExp<-setdiff(unique(evo[,2]),unique(evo[,1]))

#hei<-c()
#for(i in 1:1000){

################## PARAMETERS ###############################
# to have the main function working we have to specify two sets of parameters:
# 1. alphaExp that is the rate at which the population expands
# 2. timesExp the vector of the times at which the leaves appear. It also specifies the NUMBER of leaves present.

#####times at which leaves appear
####the number of leaves that appear is total the number of leaves that you will have in the tree
####their appearance time "determines" the length of the tree
#times<-extrLeaves2(total)
#times<-extrLeaves3(total)
#timesExp<-c(seq(0.5,1,0.1),2)
##########################################################################
#timesExp<-rep(0,50)+0.0001*c(1:50)

#alphaExp<-20

################# MAIN FUNCTION 

# MakeTreeExp<-function(alphaExp,timesExp){
# timesExp<-as.numeric(timesExp)
# ######SETTING THE STARTING PARAMETERS
# #lineagesExp is the matrix containing the extant lineagesExp and the time they appear in the genealogy either by coalescent event or by addition of a lineage
# lineag<-matrix(nrow=1,ncol=2,data=0)
# timesExp<-timesExp-min(timesExp) #first leaf at 0
# lineag[1,]<-c(1,timesExp[1]) #it is initialised with one lineage at time 0
# timesExp<-timesExp[-1]
# Nl<-1 # is the numbe of lineagesExp
# tipsExp<-1
# #total number of lineagesExp in the simulation
# #total<-10

# #alphaExp<-20

# treeHeightExp<-0
# evo<-c(0,0,0)


# ######preparing evolutionExp
# evolvedExp<-list(lineages=as.matrix(lineag),times=timesExp,height=0,tips=tipsExp,evo=evo)


# ######## evolving!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# 	#while((length(evolvedExp$lineagesExp[,1])+length(evolvedExp$timesExp)>1){	
#   while((length(evolvedExp$lineages[,1])>1)|(length(evolvedExp$times)>0)){
		
# 	evolvedExp<-evolutionExp(lineagesExp=evolvedExp$lineages,
# 				tipsExp=evolvedExp$tips,
# 				timesExp=evolvedExp$times,
# 				treeHeightExp=evolvedExp$height,
# 				alphaExp=alphaExp)
	
# 	if(is.matrix(evolvedExp$evo)){
# 		evo<-rbind(evo,evolvedExp$evo)
# 	}
# }
# tree<-makephylo(evo)
# return(tree)
# }

##################### EXAMPLE #################

#timesExp<-c(rep(0,20)+0.00001*c(1:20),2)

#alphaExp<-20
#tree<-MakeTreeExp(alphaExp=alphaExp,timesExp=timesExp)
#plot.phylo(tree)

###################### END PROGRAM-----------------------------------------------
####### here you should get an evolutionExp matrix altready transformed 
####### into a phylo object that we can then plot and analyse

###just to debug, remove afterwardslv
#lineagesExp<-evolvedExp$lineagesExp

#treeHeightExp<-evolvedExp$height
#tipsExp<-evolvedExp$tipsExp
#evo<-evolvedExp$evo

#timesExp<-evolvedExp$timesExp


#evolvedExp<-list(lineagesExp=lineagesExp,timesExp=timesExp,height=treeHeightExp,tipsExp=tipsExp,evo=evo)
