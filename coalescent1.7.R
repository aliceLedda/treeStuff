library(ape)
#########FUNCTIONS TO GET PROGRAM WORKING

####  function Evolution -------------------------------------------------------
#starts the evolution matrix
#parameters: vector of extant leaves, vector of leaves times

evolution<-function(lineages, tips, times,treeHeight){
	#lineages<-as.matrix(lineages)
	Nl<-length(lineages[,1])
	coalTime<-extrTime(Nl)  
	coalTime
	if(((treeHeight+coalTime)<times[1])|(identical(evolved$times,numeric(0)))){
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
	evoR<-evoR[-1,]

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

hei<-c()
#for(i in 1:1000){

#####times at which leaves appear
####the number of leaves that appear is total the number of leaves that you will have in the tree
####their appearance time "determines" the length of the tree
#times<-extrLeaves2(total)
#times<-extrLeaves3(total)
times<-c(seq(0.5,1,0.02),20)
#times<-rep(1,10)


######SETTING THE STARTING PARAMETERS
#lineages is the matrix containing the extant lineages and the time they appear in the genealogy either by coalescent event or by addition of a lineage
lineag<-matrix(nrow=1,ncol=2,data=0)
lineag[1,]<-c(1,times[1]) #it is initialised with one lineage at time 0
times<-times[-1]
Nl<-1 # is the numbe of lineages
tips<-1
#total number of lineages in the simulation
#total<-10


treeHeight<-0
evo<-c(0,0,0)


######preparing evolution
evolved<-list(lineages=as.matrix(lineag),times=times,height=0,tips=tips,evo=evo)


######## evolving!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	#while((length(evolved$lineages[,1])+length(evolved$times)>1){	
  while((length(evolved$lineages[,1])>1)|(length(evolved$times)>0)){
		
	evolved<-evolution(lineages=evolved$lineages,tips=evolved$tips,times=evolved$times,treeHeight=evolved$height)
	
	if(is.matrix(evolved$evo)){
		evo<-rbind(evo,evolved$evo)
	}
}
tree<-makephylo(evo)

plot.phylo(tree)

###################### END PROGRAM-----------------------------------------------
####### here you should get an evolution matrix altready transformed 
####### into a phylo object that we can then plot and analyse
