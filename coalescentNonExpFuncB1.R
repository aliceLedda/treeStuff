library(ape)
#########FUNCTIONS TO GET PROGRAM WORKING

####  function Evolution -------------------------------------------------------
#starts the evolution matrix
#parameters: vector of extant leaves, vector of leaves timesNon

evolutionNon<-function(lineagesNon, tipsNon, timesNon,heightNon){
	#lineagesNon<-as.matrix(lineagesNon)
	Nl<-length(lineagesNon[,1])
	coalTime<-extrTime(Nl)  
	coalTime
	if(((heightNon+coalTime)<as.numeric(timesNon[1]))|(identical(as.numeric(timesNon),numeric(0)))){
		coal<-chooseLeaves(length(lineagesNon[,1]))
		ma<-matrix(nrow=2,ncol=3,data=0)
		ma[,2]<-c(lineagesNon[coal[1],1],lineagesNon[coal[2],1])
		if(max(lineagesNon[,1])<(length(timesNon)+(max(tipsNon))+1)){
			ma[,1]<-rep((length(timesNon)+(max(tipsNon))+1),2)
			lineagesNon<-rbind(lineagesNon,c(0,0))
			lineagesNon[length(lineagesNon[,1]),1]<-(length(timesNon)+(max(tipsNon))+1)
			lineagesNon[length(lineagesNon[,1]),2]<-heightNon+coalTime
		}else{
			ma[,1]<-rep((max(lineagesNon[,1])+1),2)
			lineagesNon<-rbind(lineagesNon,c(0,0))
			lineagesNon[length(lineagesNon[,1]),1]<-(max(lineagesNon[,1])+1)
			lineagesNon[length(lineagesNon[,1]),2]<-heightNon+coalTime
		}
		ma[(ma[,2]==lineagesNon[coal[1],1]),3]<-(heightNon+coalTime-(lineagesNon[coal[1],2])) #this need to be fixed
		ma[(ma[,2]==lineagesNon[coal[2],1]),3]<-(heightNon+coalTime-(lineagesNon[coal[2],2]))
		lineagesNon<-lineagesNon[-(coal),,drop=FALSE]
		heightNon<-heightNon+coalTime
		
		return(list(lineages=lineagesNon,tips=tipsNon,evo=ma,height=heightNon,times=timesNon))
	}else{
		lineagesNon<-rbind(lineagesNon,c(0,0))
		lineagesNon[length(lineagesNon[,1]),1]<-(max(tipsNon)+1)
		lineagesNon[length(lineagesNon[,1]),2]<-timesNon[1]
		tipsNon<-c(tipsNon,(max(tipsNon)+1))
		heightNon<-timesNon[1]
		timesNon<-timesNon[-1]
		return(list(lineages=lineagesNon,tips=tipsNon,evo=0,height=heightNon,times=timesNon))
	}
	
}


#### function Choose Leaves --------------------------------------------------
#chooses 2 leaves out of the extant k to coalesce
#parameters: number of extant lineages
chooseLeaves<-function(l){
	choosen<-c()
	cho<-sample(1:l,2)
	#choosen[1]<-leaves[cho[1]]
	#choosen[2]<-leaves[cho[2]]
	return(cho)
}

#### function Extract Time ---------------------------------------------------
#extracts "time" from an exponential distribution with rate k(k-1)/2
#parameters: the number of extant lineages
extrTime<-function(l){
	if(l>1){
		par<-l*(l-1)/2
		time<-rexp(1,rate=par)
	}else{time<-rexp(1,rate=0.0001)}
	return(time)
}





#############function makephylo

makephylo<-function(evoR){	
	#now that the matrix is ready I transform it in a phylo object. 
	#If you get seg fault, try plot.phylo(alb), if it gives a mistake the problem is that the matrix is not properly declared.
	tips<-setdiff(unique(evoR[,2]),unique(evoR[,1]))
	#if(identical(evoR[1,],c(0,0,0),num.eq=TRUE)){
	evoR<-evoR[-1,]
	#}

	#tips<-setdiff(unique(evo[,2]),unique(evo[,1]))

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

###idenfy the tips....
######### the tips are all the branches from which there is no other branch
#tips<-setdiff(unique(evo[,2]),unique(evo[,1]))

#hei<-c()
#for(i in 1:1000){

########## ####### PARAMETERS ##############################
# for the main function to work we have to specify the times at which the leaves are sampled
#####timesNon at which leaves appear
####the number of leaves that appear is total the number of leaves that you will have in the tree
####their appearance time "determines" the length of the tree
#timesNon<-extrLeaves2(total)
#timesNon<-extrLeaves3(total)
#timesNon<-c(seq(0.5,1,0.02),20)
#timesNon<-c(seq(0.5,1,0.1),2)
#timesNon<-rep(0,10)+0.3*c(1:10)


#################### MAIN FUNCTION #########################

MakeNonTree<-function(timesNon){

######SETTING THE STARTING PARAMETERS
#lineages is the matrix containing the extant lineages and the time they appear in the genealogy either by coalescent event or by addition of a lineage
lineag<-matrix(nrow=1,ncol=2,data=0)
timesNon<-timesNon-min(timesNon)  #first leaf at 0
lineag[1,]<-c(1,timesNon[1]) #it is initialised with one lineage at time 0
timesNon<-timesNon[-1]
Nl<-1 # is the numbe of lineages
tipsNon<-1
#total number of lineages in the simulation
#total<-10


heightNon<-0
evo<-c(0,0,0)


######preparing evolution
evolvedNon<-list(lineages=as.matrix(lineag),times=as.numeric(timesNon),height=0,tips=tipsNon,evo=evo)


######## evolving!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	#while((length(evolvedNon$lineagesNon[,1])+length(evolvedNon$timesNon)>1){	
  while((length(evolvedNon$lineages[,1])>1)|(length(evolvedNon$times)>0)){
		
	evolvedNon<-evolutionNon(lineagesNon=evolvedNon$lineages,
				tipsNon=evolvedNon$tips,
				timesNon=as.numeric(evolvedNon$times),
				heightNon=evolvedNon$height)
	
	if(is.matrix(evolvedNon$evo)){
		evo<-rbind(evo,evolvedNon$evo)
	}
}
return(list(evo=evo, height=evolvedNon$height ))
}


MakeTreeNon<-function(timesNon){

done<-MakeNonTree(timesNon)

tree<-makephylo(done$evo)
return(tree)
}

#################### EXAMPLE ################


#tree<-MakeTreeNon(timesNon=timesNon)
#plot.phylo(tree)

###################### END PROGRAM-----------------------------------------------
####### here you should get an evolution matrix altready transformed 
####### into a phylo object that we can then plot and analyse
