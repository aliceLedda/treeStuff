library(ape)
#########FUNCTIONS TO GET PROGRAM WORKING

####  function evolutionExp -------------------------------------------------------
#gives the evolution matrix of the exponentially growing population
#parameters: vector of extant leaves, vector of leaves times

evolutionExp<-function(lineages, tips, timesExp,treeHeight){
	#lineages<-as.matrix(lineages)
	Nl<-length(lineages[,1])
	coalTime<-extrTimeExp(Nl,treeHeight)  
	coalTime
	#while(((treeHeight+coalTime)<timesExp[1])&(length(lineages[,1])==1)){
	#	coalTime<-extrTime(Nl,treeHeight) 
	#}
	
	if(((treeHeight+coalTime)<timesExp[1])|(identical(timesExp,numeric(0)))){
		coal<-chooseLeaves(length(lineages[,1]))
		ma<-matrix(nrow=2,ncol=3,data=0)
		ma[,2]<-c(lineages[coal[1],1],lineages[coal[2],1])
		if(max(lineages[,1])<(length(timesExp)+(max(tips))+1)){
			ma[,1]<-rep((length(timesExp)+(max(tips))+1),2)
			lineages<-rbind(lineages,c(0,0))
			lineages[length(lineages[,1]),1]<-(length(timesExp)+(max(tips))+1)
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
		
		return(list(lineages=lineages,tips=tips,evo=ma,height=treeHeight,timesExp=timesExp))
	}else{
		lineages<-rbind(lineages,c(0,0))
		lineages[length(lineages[,1]),1]<-(max(tips)+1)
		lineages[length(lineages[,1]),2]<-timesExp[1]
		tips<-c(tips,(max(tips)+1))
		treeHeight<-timesExp[1]
		timesExp<-timesExp[-1]
		return(list(lineages=lineages,tips=tips,evo=0,height=treeHeight,timesExp=timesExp))
	}
	
}

####  function evolutionNon -------------------------------------------------------
#gives the evolution matrix of the non-exponentially expanding population
#parameters: vector of extant leaves, vector of leaves timesNon

evolutionNon<-function(lineages, tips, timesNon,treeHeight){
	#lineages<-as.matrix(lineages)
	Nl<-length(lineages[,1])
	coalTime<-extrTime(Nl)  
	coalTime
	if(((treeHeight+coalTime)<timesNon[1])|(identical(timesNon,numeric(0)))){
		coal<-chooseLeaves(length(lineages[,1]))
		ma<-matrix(nrow=2,ncol=3,data=0)
		ma[,2]<-c(lineages[coal[1],1],lineages[coal[2],1])
		if(max(lineages[,1])<(length(timesNon)+(max(tips))+1)){
			ma[,1]<-rep((length(timesNon)+(max(tips))+1),2)
			lineages<-rbind(lineages,c(0,0))
			lineages[length(lineages[,1]),1]<-(length(timesNon)+(max(tips))+1)
			lineages[length(lineages[,1]),2]<-treeHeight+coalTime
		}else{
			ma[,1]<-rep((max(lineages[,1])+1),2)
			lineages<-rbind(lineages,c(0,0))
			lineages[length(lineages[,1]),1]<-(max(lineages[,1])+1)
			lineages[length(lineages[,1]),2]<-treeHeight+coalTime
		}
		ma[(ma[,2]==lineages[coal[1],1]),3]<-(treeHeight+coalTime-(lineages[coal[1],2])) #this need to be fixed
		ma[(ma[,2]==lineages[coal[2],1]),3]<-(treeHeight+coalTime-(lineages[coal[2],2]))
		lineages<-lineages[-(coal),,drop=FALSE]
		treeHeight<-treeHeight+coalTime
		
		return(list(lineages=lineages,tips=tips,evo=ma,height=treeHeight,timesNon=timesNon))
	}else{
		lineages<-rbind(lineages,c(0,0))
		lineages[length(lineages[,1]),1]<-(max(tips)+1)
		lineages[length(lineages[,1]),2]<-timesNon[1]
		tips<-c(tips,(max(tips)+1))
		treeHeight<-timesNon[1]
		timesNon<-timesNon[-1]
		return(list(lineages=lineages,tips=tips,evo=0,height=treeHeight,timesNon=timesNon))
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

#### function Extract Time Exp ---------------------------------------------------
#extracts "time" from an exponential distribution with rate k(k-1)/2
#and rescales it for an exponentially expanding population with rate alpha
#parameters: the number of extant lineages
extrTimeExp<-function(l, treeHeight){
	if(l>1){
		par<-l*(l-1)/2
		timeO<-rexp(1,rate=par)
		timeNew<-(1/alpha)*log((alpha*timeO)+exp(alpha*treeHeight))-treeHeight
	}else{timeNew<-100000}
	
	return(timeNew)
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

#########function joinevolved

joinevolved<-function(popExp,popNon){
	if (length(popExp$tips)>=length(popNon$tips)){
	
		newlin<-popExp$lineages
		newlin[,1]<-newlin[,1]+max(popNon$lineages) 
		ancestralLineages<-rbind(popNon$lineages,newlin)
		ancestralTimes<-c(popNon$timesNon,popExp$timesExp)  #maybe not so easy
		ancestralHeight<-max(popNon$height,popExp$height)   #maybe not so easy too
		NewExpTips<-popExp$tips+max(popNon$lineages)     #maybe not so easy either
		ancestralTips<-c(popNon$tips, NewExpTips)
		ancestralEvo<-0
	}
	else{
		newlin<-popNon$lineages
		newlin[,1]<-newlin[,1]+max(popExp$lineages) 
		ancestralLineages<-rbind(popExp$lineages,newlin)
		ancestralTimes<-c(popExp$timesExp,popNon$timesNon)  #maybe not so easy
		ancestralHeight<-max(popExp$height,popNon$height)   #maybe not so easy too
		NewNonTips<-popNon$tips+max(popExp$lineages)     #maybe not so easy either
		ancestralTips<-c(popExp$tips, NewNonTips)
		ancestralEvo<-0
	}
	
	return (list(lineagesA=ancestralLineages,timesA=ancestralTimes,heightA=ancestralHeight,tipsA=ancestralTips,evoA=ancestralEvo))
	}


#############function joinevo

joinevo<-function(evoExp,evoNon){	
	#now that the matrix is ready I transform it in a phylo object. 
	#If you get seg fault, try plot.phylo(alb), if it gives a mistake the problem is that the matrix is not properly declared.
	tipsExp<-setdiff(unique(evoExp[,2]),unique(evoExp[,1]))
	tipsNon<-setdiff(unique(evoNon[,2]),unique(evoNon[,1]))
	evoExp<-evoExp[-1,]
	evoNon<-evoNon[-1,]
	evoNew<-matrix(ncol=ncol(evoExp),nrow=sum(nrow(evoExp),nrow(evoNon)))
	#tips<-setdiff(unique(evo[,2]),unique(evo[,1]))
	if (max(evoExp)>=max(evoNon)){
		evoExp[,1]<-evoExp[,1]+max(evoNon)
		evoExp[,2]<-evoExp[,2]+max(evoNon)
	
		evoNew<-rbind(evoNon,evoExp)
	}
	else{
		evoNon[,1]<-evoNon[,1]+max(evoExp)
		evoNon[,2]<-evoNon[,2]+max(evoExp)
		
		evoNew<-rbind(evoExp,evoNon)
	}
	return(evoNew)
}


#############function makephylo

makephylo<-function(evoR){	
	#now that the matrix is ready I transform it in a phylo object. 
	#If you get seg fault, try plot.phylo(alb), if it gives a mistake the problem is that the matrix is not properly declared.
	tips<-setdiff(unique(evoR[,2]),unique(evoR[,1]))
	evoNew<-matrix(nrow=nrow(evoR),ncol=ncol(evoR),data=0)
	#evoR<-evoR[-1,]
	evoNew[(evoR[,2] %in% tips),2]<-c(1:length(tips))
	
	#tips<-setdiff(unique(evo[,2]),unique(evo[,1]))
	minimum<-length(tips)+2-min(evoR[evoR[,1]<length(tips),1])
	#minimum<-min(evoR[evoR[,1]<length(tips),1])
	#evoR[,1]<-evoR[,1]+1
	evoNew[evoR[,1]>=length(tips),1]<-evoR[evoR[,1]>=length(tips),1]+1
	evoNew[evoR[,1]<length(tips),1]<-evoR[evoR[,1]<length(tips),1]+minimum
	
	evoNew[((evoR[,2]>=length(tips))&(!(evoR[,2] %in% tips))),2]<-evoR[((evoR[,2]>=length(tips))&(!(evoR[,2] %in% tips))),2]+1
	
	evoNew[((evoR[,2]<length(tips))&(!(evoR[,2] %in% tips))),2]<-evoR[((evoR[,2]<length(tips))&(!(evoR[,2] %in% tips))),2]+minimum
	
	#evoR[!(evoR[,2] %in% tips),2]<-evoR[!(evoR[,2] %in% tips),2]+1

	evoNew[evoR[,1]==max(evoR[,1]),1]<-length(tips)+1
	
	#tips[tips>length(tips)]<-tips[tips>length(tips)]+1
	
	
	tmp<-evoNew[,1:2]
	tree<-list(edge=evoNew[,1:2])
	class(tree)<-"phylo"
	
	#attr(tree,"order")<-"postorder"
	#attr(tree,"order")<-"cladewise"
	tree$edge.length<-evoR[,3]
	tree$tip.label<-seq(1,length(tips),1)
	tree$tip.label<-tips
	tree$Nnode<-length(unique(evoNew[,1]))
	#back<-list(tree=tree,cBR=cBR)
	return(tree)
}
 	
	
 
##########################################################
########### THE PROGRAM ##########################
##################################################

########## GLOBAL PARAMETERS
####### The parameters of the problem

######### times at which we sample the leaves in the Exponentially growing pop
#
timesExp<-c(seq(0.5,1,0.02))
#timesExp<-rep(0,10)

####### times at which we sample the leaves in the NON-Exponentially growing pop
timesNon<-c(seq(0.5,1,0.01))

####### the rate at which the Exponentially expanding pop is expanding
alpha=10

######## The times for the tree:
########## Ts>=Texp       
#check this!!
############# 1. the time at which the two populations split
Ts<-5
############ 2. the time at which the Exp populations start expanding
Texp<-3


######SETTING THE STARTING PARAMETERS
#lineages is the matrix containing the extant lineages 
#and the time they appear in the genealogy either by coalescent event or by addition of a lineage
lineagNon<-matrix(nrow=1,ncol=2,data=0)
lineagExp<-matrix(nrow=1,ncol=2,data=0)

lineagNon[1,]<-c(1,timesNon[1]) #it is initialised with one lineage at time 0
lineagExp[1,]<-c(1,timesExp[1]) #it is initialised with one lineage at time 0

timesNon<-timesNon[-1]
timesExp<-timesExp[-1]
NlNon<-1 # is the numbe of lineages
NlExp<-1
tipsNon<-1
tipsExp<-1
#total number of lineages in the simulation
#total<-10


treeHeightNon<-0
treeHeightExp<-0
evoNon<-c(0,0,0)
evoExp<-c(0,0,0)

######preparing evolution
evolvedNon<-list(lineages=as.matrix(lineagNon),timesNon=timesNon,height=0,tips=tipsNon,evo=evoNon)

evolvedExp<-list(lineagesExp=as.matrix(lineagExp),timesExp=timesExp,heightExp=0,tipsExp=tipsExp,evoExp=evoExp)

evolvedA<-list(lineagesA=as.matrix(lineagExp),timesA=0,heightA=0,tipsA=0,evoA=0)

######## evolving!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
control=0
while((evolvedExp$height<Ts)&(control==0)){
	while((evolvedExp$height<Texp)&(control==0)){
		if((length(evolvedExp$lineages[,1])>1)|(length(evolvedExp$timesExp)>0)){
		 evolvedExp<-evolutionExp(lineages=evolvedExp$lineages,tips=evolvedExp$tips,timesExp=evolvedExp$timesExp,treeHeight=evolvedExp$height)
			if(is.matrix(evolvedExp$evo)){
				evoExp<-rbind(evoExp,evolvedExp$evo)
			}
		}else{control=1}
	}
	if((length(evolvedExp$lineages[,1])>1)|(length(evolvedExp$timesExp)>0)){
		evolvedExp<-evolutionNon(lineages=evolvedExp$lineages,tips=evolvedExp$tips,timesExp=evolvedExp$timesExp,treeHeight=evolvedExp$height)
		if(is.matrix(evolvedExp$evo)){
			evoExp<-rbind(evoExp,evolvedExp$evo)
			}
		}
	else{control=1}
}

control=0
while((evolvedNon$height<Ts)&(control==0)){
	if((length(evolvedNon$lineages[,1])>1)|(length(evolvedNon$timesNon)>0)){
	evolvedNon<-evolutionNon(lineages=evolvedNon$lineages,tips=evolvedNon$tips,timesNon=evolvedNon$timesNon,treeHeight=evolvedNon$height)
		if(is.matrix(evolvedNon$evo)){
			evoNon<-rbind(evoNon,evolvedNon$evo)
		}
	}
	else{control=1}
}

evolvedA<-joinevolved(evolvedExp,evolvedNon)

evoA<-joinevo(evoExp,evoNon)

#while((length(evolved$lineages[,1])+length(evolved$times)>1){	
  while((length(evolvedA$lineages[,1])>1)|(length(evolvedA$times)>0)){
		
	evolvedA<-evolutionNon(lineages=evolvedA$lineages,tips=evolvedA$tips,times=evolvedA$times,treeHeight=evolvedA$height)
	
	if(is.matrix(evolvedA$evo)){
		evoA<-rbind(evoA,evolvedA$evo)
	}
}

tree<-makephylo(evoA)

plot.phylo(tree)


#treeExp<-makephylo(evoExp)
#plot.phylo(treeExp)

###################### END PROGRAM-----------------------------------------------
####### here you should get an evolution matrix altready transformed 
####### into a phylo object that we can then plot and analyse

###just to debug, remove afterwardslv
#lineages<-evolvedExp$lineages

#treeHeight<-evolvedExp$height
#tips<-evolvedExp$tips
#evo<-evolvedExp$evo

#timesExp<-evolvedExp$timesExp


#evolved<-list(lineages=lineages,times=times,height=treeHeight,tips=tips,evo=evo)
