library(ape)
#########FUNCTIONS TO GET PROGRAM WORKING

####  function Evolution -------------------------------------------------------
#starts the evolution matrix
#parameters: vector of extant leaves, vector of leaves timesExp

evolution<-function(lineages, tips, times,treeHeight,alpha,ratio, expns){
	if(missing(ratio)){ratio<-0}
	#lineagesExp<-as.matrix(lineagesExp)
	Nl<-length(lineages[,1])

	coalTime<-extrTime(l=Nl,treeHeight=treeHeight,alpha=alpha, ratio=ratio, expns=expns)  
	
	coalTime
	#while(((treeHeightExp+coalTimeExp)<timesExp[1])&(length(lineagesExp[,1])==1)){
	#	coalTimeExp<-extrTimeExp(Nl,treeHeightExp) 
	#}
	
	if(((treeHeight+coalTime)<times[1])|(identical(times,numeric(0)))){
		coal<-chooseLeaves(length(lineages[,1]))
		ma<-matrix(nrow=2,ncol=3,data=0)
		ma[,2]<-c(lineages[coal[1],1],lineages[coal[2],1])
		if(max(lineages[,1])<(length(times)+(max(tips))+1)){
			ma[,1]<-rep((length(times)+(max(tips))+1),2)
			lineages<-rbind(lineages,c(0,0))
			lineages[length(lineages[,1]),1]<-(length(times)+(max(tips))+1)
			lineages[length(lineages[,1]),2]<-treeHeight+coalTime
		}
		
		else{
			ma[,1]<-rep((max(lineages[,1])+1),2)
			lineages<-rbind(lineages,c(0,0))
			lineages[length(lineages[,1]),1]<-(max(lineages[,1])+1)
			lineages[length(lineages[,1]),2]<-treeHeight+coalTime
		}
		ma[(ma[,2]==lineages[coal[1],1]),3]<-(treeHeight+coalTime-(lineages[coal[1],2])) #this need to be fixed
		ma[(ma[,2]==lineages[coal[2],1]),3]<-(treeHeight+coalTime-(lineages[coal[2],2]))
		lineages<-lineages[-(coal),,drop=FALSE]
		treeHeight<-treeHeight+coalTime
		
		return(list(lineages=lineages,tips=tips,evo=ma,height=treeHeight,times=times))
	}else{
		lineages<-rbind(lineages,c(0,0))
		lineages[length(lineages[,1]),1]<-(max(tips)+1)
		lineages[length(lineages[,1]),2]<-times[1]
		tips<-c(tips,(max(tips)+1))
		treeHeight<-times[1]
		times<-times[-1]
		return(list(lineages=lineages,tips=tips,evo=0,height=treeHeight,times=times))
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
extrTime<-function(l, treeHeight, alpha,ratio,expns){
	if ((missing(ratio))|(ratio==0)) {ratio=1}
	if(l>1){
		par<-l*(l-1)/(2*ratio)
		timeO<-rexp(1,rate=par)
		if((alpha==0)|(expns=="none")){timeNew<-timeO
		}else{
			if(expns=="exp"){
				timeNew<-((1/alpha)*log((alpha*timeO)+exp(alpha*treeHeight)))-treeHeight}
	}
	}else{timeNew<-100000}
	
	return(timeNew)
}





#############function makephylo

makephylo<-function(evoR){	
	#now that the matrix is ready I transform it in a phylo object. 
	#If you get seg fault, try plot.phylo(alb), if it gives a mistake the problem is that the matrix is not properly declared.
	tips<-setdiff(unique(evoR[,2]),unique(evoR[,1]))
	#if(identical(evoR[1,],c(0,0,0),num.eq=TRUE)){
	evoR<-evoR[-1,]
	#}

	#tipsExp<-setdiff(unique(evo[,2]),unique(evo[,1]))

	evoR[,1]<-evoR[,1]+1

	evoR[!(evoR[,2] %in% tips),2]<-evoR[!(evoR[,2] %in% tips),2]+1

	evoR[evoR[,1]==max(evoR[,1]),1]<-max(tips)+1

	
	tmp<-evoR[,1:2]
	tree<-list(edge=evoR[,1:2])
	class(tree)<-"phylo"
	
	attr(tree,"order")<-"postorder"
	
	tree$edge.length<-evoR[,3]
	tree$tip.label<-seq(1,length(tips),1)
	tree$Nnode<-length(unique(evoR[,1]))
	#back<-list(tree=tree,cBR=cBR)
	return(tree)
}
 
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
###########3timesExp<-rep(0,50)+0.0001*c(1:50)

###########alphaExp<-20

################# MAIN FUNCTION 

TreeMake<-function(times,alpha,ratio, expns){
	if(missing(times)){ratio<-0
		}
######SETTING THE STARTING PARAMETERS
#lineagesExp is the matrix containing the extant lineagesExp and the time they appear in the genealogy either by coalescent event or by addition of a lineage
lineag<-matrix(nrow=1,ncol=2,data=0)
#timesExp<-timesExp-min(timesExp) #first leaf at 0
lineag[1,]<-c(1,times[1]) #it is initialised with one lineage at time 0
times<-times[-1]
Nl<-1 # is the numbe of lineagesExp
tips<-1
#total number of lineagesExp in the simulation
#total<-10

#alphaExp<-20

treeHeight<-0
evo<-c(0,0,0)


######preparing evolutionExp
evolved<-list(lineages=as.matrix(lineag),times=times,height=0,tips=tips,evo=evo,ratio=ratio,alpha=alpha,expns=expns)


######## evolving!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	#while((length(evolvedExp$lineagesExp[,1])+length(evolvedExp$timesExp)>1){	
  while((length(evolved$lineages[,1])>1)|(length(evolved$times)>0)){
		
	evolved<-evolution(lineages=evolved$lineages,
				tips=evolved$tips,
				times=evolved$times,
				treeHeight=evolved$height,
				alpha=alpha,
				ratio=ratio,
				expns=expns)
	
	if(is.matrix(evolved$evo)){
		evo<-rbind(evo,evolved$evo)
	}
}
return(list(evo=evo, height=evolved$height,evolved=evolved))
}



MakeTree<-function(times,alpha,ratio,expns){

if(missing(times)) {stop("you need to specify at least times")}	
if((missing(ratio))&(missing(alpha))&(missing(expns))) {
	expns<-"none"
	ratio<-1
	alpha<-0
}	
if((missing(ratio))| (ratio==0)) {ratio<-1}

if(((alpha==0)&(expns!="none"))|((alpha!=0)&(expns=="none"))) {stop("your request is incoherent")}
if((missing(alpha))&(expns!="none")) {
	alpha<-0
}
expns1<-tolower(expns)
expns<-expns1

times<-times-min(times)  ##### Set the first time to zero

done<-TreeMake(times, alpha,ratio, expns)

tree<-makephylo(done$evo)
return(tree)
}

##################### EXAMPLE #################

timesExp<-c(rep(0,20)+0.00001*c(1:20),2)

alphaExp<-20
tree<-MakeTree(times=rep(0,20),alpha=0,ratio=1,expns="none")
plot.phylo(tree)

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
