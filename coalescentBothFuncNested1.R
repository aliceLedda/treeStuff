library(ape)
#########FUNCTIONS TO GET PROGRAM WORKING
source("coalescentExpFuncB1.R")
source("coalescentNonExpFuncB1.R")

#########function joinevolved

joinevolved<-function(popExp,popNon,Ts){
	if (length(popExp$tips)>=length(popNon$tips)){
	
		newlin<-popExp$lineages
		newlin[,1]<-newlin[,1]+max(popNon$lineages) 
		ancestralLineages<-rbind(popNon$lineages,newlin)
		if((identical(popNon$times,numeric(0)))&(identical(popExp$times,numeric(0)))){
		ancestralTimes<-numeric(0)
		}else{
		ancestralTimes<-c(popNon$times,popExp$times)}  #maybe not so easy
		ancestralHeight<-max(popNon$height,popExp$height,Ts)   #maybe not so easy too
		NewExpTips<-popExp$tips+max(popNon$lineages[,1])     #maybe not so easy either
		ancestralTips<-c(popNon$tips, NewExpTips)
		ancestralEvo<-0
	}
	else{
		newlin<-popNon$lineages
		newlin[,1]<-newlin[,1]+max(popExp$lineages) 
		ancestralLineages<-rbind(popExp$lineages,newlin)
		if((identical(popNon$times,numeric(0)))&(identical(popExp$times,numeric(0)))){
		ancestralTimes<-numeric(0)
		}else{
		ancestralTimes<-c(popExp$times,popNon$times)} #maybe not so easy
		ancestralHeight<-max(popExp$height,popNon$height,Ts)   #maybe not so easy too
		NewNonTips<-popNon$tips+max(popExp$lineages[,1])     #maybe not so easy either
		ancestralTips<-c(popExp$tips, NewNonTips)
		ancestralEvo<-0
	}
	
	return (list(lineages=ancestralLineages,times=ancestralTimes,height=ancestralHeight,tips=ancestralTips,evo=ancestralEvo))
	}

# joinevolved<-function(popExp,popNon){
# 	if (length(popExp$tips)>=length(popNon$tips)){
	
# 		newlin<-popExp$lineages
# 		newlin[,1]<-newlin[,1]+max(popNon$lineages) 
# 		ancestralLineages<-rbind(popNon$lineages,newlin)
# 		if((identical(popNon$times,numeric(0)))&(identical(popExp$times,numeric(0)))){
# 		ancestralTimes<-numeric(0)
# 		}else{
# 		ancestralTimes<-c(popNon$times,popExp$times)}  #maybe not so easy
# 		ancestralHeight<-max(popNon$height,popExp$height)   #maybe not so easy too
# 		NewExpTips<-popExp$tips+max(popNon$lineages[,1])     #maybe not so easy either
# 		ancestralTips<-c(popNon$tips, NewExpTips)
# 		ancestralEvo<-0
# 	}
# 	else{
# 		newlin<-popNon$lineages
# 		newlin[,1]<-newlin[,1]+max(popExp$lineages) 
# 		ancestralLineages<-rbind(popExp$lineages,newlin)
# 		if((identical(popNon$times,numeric(0)))&(identical(popExp$times,numeric(0)))){
# 		ancestralTimes<-numeric(0)
# 		}else{
# 		ancestralTimes<-c(popExp$times,popNon$times)} #maybe not so easy
# 		ancestralHeight<-max(popExp$height,popNon$height)   #maybe not so easy too
# 		NewNonTips<-popNon$tips+max(popExp$lineages[,1])     #maybe not so easy either
# 		ancestralTips<-c(popExp$tips, NewNonTips)
# 		ancestralEvo<-0
# 	}
	
# 	return (list(lineages=ancestralLineages,times=ancestralTimes,height=ancestralHeight,tips=ancestralTips,evo=ancestralEvo))
# 	}


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

makephyloAll<-function(evoR){	
	#now that the matrix is ready I transform it in a phylo object. 
	#If you get seg fault, try plot.phylo(alb), if it gives a mistake the problem is that the matrix is not properly declared.
	tips<-setdiff(unique(evoR[,2]),unique(evoR[,1]))
	evoNew<-matrix(nrow=nrow(evoR),ncol=ncol(evoR),data=0)
	
	evoNew[(evoR[,2] %in% tips),2]<-c(1:length(tips))
	
	
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
	tree$tip.label<-tips   #this is the line causing the funny tip names (not consecutive). It can be deleted with no further consequences than the tips getting consecutive numbers names.
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
timesExp<-c(seq(0.5,0.7,0.01))
#timesExp<-rep(0,10)

####### times at which we sample the leaves in the NON-Exponentially growing pop
timesNon<-c(seq(0.6,1,0.02))

####### the rate at which the Exponentially expanding pop is expanding
alpha=10

######## The times for the tree:
########## Ts>=Texp       
#check this!!
############# 1. the time at which the two populations split
Ts<-10
############ 2. the time at which the Exp populations start expanding
Texp<-3

# ############# MAIN FUNCTION ###############

# MakeTreeBoth<-function(timesExp,timesNon,alpha,Ts,Texp){

# ######SETTING THE STARTING PARAMETERS
# #lineages is the matrix containing the extant lineages 
# #and the time they appear in the genealogy either by coalescent event or by addition of a lineage
# lineagNon<-matrix(nrow=1,ncol=2,data=0)
# lineagExp<-matrix(nrow=1,ncol=2,data=0)

# timesNon<-timesNon-min(timesNon)   #first Leave at 0
# timesExp<-timesExp-min(timesExp)   # first leave at 0

# lineagNon[1,]<-c(1,timesNon[1]) #it is initialised with one lineage at time 0
# lineagExp[1,]<-c(1,timesExp[1]) #it is initialised with one lineage at time 0

# timesNon<-timesNon[-1]
# timesExp<-timesExp[-1]
# NlNon<-1 # is the numbe of lineages
# NlExp<-1
# tipsNon<-1
# tipsExp<-1
# #total number of lineages in the simulation
# #total<-10


# heightNon<-0
# height<-0
# evoNon<-c(0,0,0)
# evoExp<-c(0,0,0)

# ######preparing evolution
# evolvedNon<-list(lineages=as.matrix(lineagNon),
# 		times=timesNon,
# 		height=0,
# 		tips=tipsNon,
# 		evo=evoNon)

# evolvedExp<-list(lineages=as.matrix(lineagExp),
# 		times=timesExp,
# 		height=0,
# 		tips=tipsExp,
# 		evo=evoExp)

# evolvedA<-list(lineages=as.matrix(lineagExp),
# 		times=0,
# 		height=0,
# 		tips=0,
# 		evo=0)

# ######## evolving!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# control=0
# while((evolvedExp$height<Ts)&(control==0)){
# 	while((evolvedExp$height<Texp)&(control==0)){
# 		if((length(evolvedExp$lineages[,1])>1)|(length(evolvedExp$times)>0)){
# 		 evolvedExp<-evolutionExp(lineagesExp=evolvedExp$lineages,
# 					tipsExp=evolvedExp$tips,
# 					timesExp=evolvedExp$times,
# 					treeHeightExp=evolvedExp$height,
# 					alphaExp=alpha)
# 			if(is.matrix(evolvedExp$evo)){
# 				evoExp<-rbind(evoExp,evolvedExp$evo)
# 			}
# 		}else{control=1}
# 	}
# 	if((length(evolvedExp$lineages[,1])>1)|(length(evolvedExp$times)>0)){
# 		evolvedExp<-evolutionNon(lineagesNon=evolvedExp$lineages,
# 					tipsNon=evolvedExp$tips,
# 					timesNon=evolvedExp$times,
# 					heightNon=evolvedExp$height)
# 		if(is.matrix(evolvedExp$evo)){
# 			evoExp<-rbind(evoExp,evolvedExp$evo)
# 			}
# 		}
# 	else{control=1}
# }

# #treeExp<-makephylo(evoExp)

# control=0
# while((evolvedNon$height<Ts)&(control==0)){
# 	if((length(evolvedNon$lineages[,1])>1)|(length(evolvedNon$times)>0)){
# 	evolvedNon<-evolutionNon(lineagesNon=evolvedNon$lineages,
# 				tipsNon=evolvedNon$tips,
# 				timesNon=evolvedNon$times,
# 				heightNon=evolvedNon$height)
# 		if(is.matrix(evolvedNon$evo)){
# 			evoNon<-rbind(evoNon,evolvedNon$evo)
# 		}
# 	}
# 	else{control=1}
# }

# #treeNon<-makephylo(evoNon)

# evolvedA<-joinevolved(evolvedExp,evolvedNon,Ts)

# evoA<-joinevo(evoExp,evoNon)

# #while((length(evolved$lineages[,1])+length(evolved$times)>1){	
#   while((length(evolvedA$lineages[,1])>1)|(length(evolvedA$times)>0)){
		
# 	evolvedA<-evolutionNon(lineagesNon=evolvedA$lineages,
# 				tipsNon=evolvedA$tips,
# 				timesNon=evolvedA$times,
# 				heightNon=evolvedA$height)
	
# 	if(is.matrix(evolvedA$evo)){
# 		evoA<-rbind(evoA,evolvedA$evo)
# 	}
# }

# tree<-makephyloAll(evoA)
# return(tree)
# }

######################   example ################

# timesExp<-c(seq(0.5,4,0.1))
# ####### times at which we sample the leaves in the NON-Exponentially growing pop
# timesNon<-c(seq(0.1,3,0.1))
# ####### the rate at which the Exponentially expanding pop is expanding
# alpha=10
# ############# 1. the time at which the two populations split
# Ts<-5
# ############ 2. the time at which the Exp populations start expanding
# Texp<-3


# tree<-MakeTreeBoth(timesExp,timesNon,alpha,Ts,Texp)
# plot.phylo(tree)


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

MakeTreeBoth<-function(timesExp,timesNon,alpha,Ts,Texp){

######SETTING THE STARTING PARAMETERS
#lineages is the matrix containing the extant lineages 
#and the time they appear in the genealogy either by coalescent event or by addition of a lineage
lineagNon<-matrix(nrow=1,ncol=2,data=0)
lineagExp<-matrix(nrow=1,ncol=2,data=0)

timesNon<-timesNon-min(timesNon)   #first Leave at 0
timesExp<-timesExp-min(timesExp)   # first leave at 0

lineagNon[1,]<-c(1,timesNon[1]) #it is initialised with one lineage at time 0
lineagExp[1,]<-c(1,timesExp[1]) #it is initialised with one lineage at time 0

timesNon<-timesNon[-1]
#timesExp<-timesExp[-1]
NlNon<-1 # is the numbe of lineages
#NlExp<-1
tipsNon<-1
#tipsExp<-1
#total number of lineages in the simulation
#total<-10


heightNon<-0
height<-0
evoNon<-c(0,0,0)
#evoExp<-c(0,0,0)

######preparing evolution
evolvedNon<-list(lineages=as.matrix(lineagNon),
		times=timesNon,
		height=0,
		tips=tipsNon,
		evo=evoNon)

#evolvedExpB<-list(lineages=as.matrix(lineagExp),
		#times=timesExp,
		#height=0,
		#tips=tipsExp,
		#evo=evoExp)

 evolvedA<-list(lineages=as.matrix(lineagExp),
 		times=0,
 		height=0,
 		tips=0,
 		evo=0)

######## evolving!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
done<-MakeExpTree(alpha,timesExp,timesNon)
while(done$height>Ts){
	done<-MakeExpTree(alphaExp,timesExp)
}

#treeExp<-makephylo(done$evo)
# plot.phylo(treeExp)

control=0
while((evolvedNon$height<Ts)&(control==0)){
	if((length(evolvedNon$lineages[,1])>1)|(length(evolvedNon$times)>0)){
	evolvedNon<-evolutionNon(lineagesNon=evolvedNon$lineages,
				tipsNon=evolvedNon$tips,
				timesNon=evolvedNon$times,
				heightNon=evolvedNon$height)
		if(is.matrix(evolvedNon$evo)){
			evoNon<-rbind(evoNon,evolvedNon$evo)
		}
	}
	else{control=1}
}

#treeNon<-makephylo(evoNon)

evolvedA<-joinevolved(done$evolvedExp,evolvedNon,Ts)

evoA<-joinevo(done$evo,evoNon)


#while((length(evolved$lineages[,1])+length(evolved$times)>1){	
  while((length(evolvedA$lineages[,1])>1)|(length(evolvedA$times)>0)){
		
	evolvedA<-evolutionNon(lineagesNon=evolvedA$lineages,
				tipsNon=evolvedA$tips,
				timesNon=evolvedA$times,
				heightNon=evolvedA$height)
	
	if(is.matrix(evolvedA$evo)){
		evoA<-rbind(evoA,evolvedA$evo)
	}
}

tree<-makephyloAll(evoA)
#plot.phylo(tree)
#axisPhylo()


return(tree)
}



timesExp<-c(rep(0,10))
####### times at which we sample the leaves in the NON-Exponentially growing pop
timesNon<-c(rep(0,20))
####### the rate at which the Exponentially expanding pop is expanding
alpha=0
############# 1. the time at which the two populations split
Ts<-1
############ 2. the time at which the Exp populations start expanding
Texp<-1

tree<-MakeTreeBoth(timesExp,timesNon,alpha,Ts,Texp)
plot.phylo(tree)
axisPhylo()

