library(ape)
#########FUNCTIONS TO GET PROGRAM WORKING

source("coalescentEvoLL.R")

#########function joinevolved

joinevolved<-function(popExp,popNon,Ts){
	if ((identical(popExp$times,numeric(0)))&(identical(popNon$times,numeric(0)))) {
		NtipsExp<-length(popExp$tips)
		NtipsNon<-length(popNon$tips)
	}
	else{
		if(identical(popExp$times,numeric(0))) {
			NtipsExp<-length(popExp$tips)
			NtipsNon<-length(popNon$tips)+length(popNon$times)
		}
		else{
			NtipsExp<-length(popExp$tips)+length(popExp$times)
			NtipsNon<-length(popNon$tips)
		}
	}
	if (NtipsExp>=NtipsNon){
	
		newlin<-popExp$lineages
		newlin[,1]<-newlin[,1]+max(popNon$lineages) 
		ancestralLineages<-rbind(popNon$lineages,newlin)
		if((identical(popNon$times,numeric(0)))&(identical(popExp$times,numeric(0)))){
		ancestralTimes<-numeric(0)
		}else{
		ancestralTimes<-c(popNon$times,popExp$times)}  #maybe not so easy
		ancestralHeight<-max(popNon$treeHeight,popExp$treeHeight,Ts)   #maybe not so easy too
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
		ancestralHeight<-max(popExp$treeHeight,popNon$treeHeight,Ts)   #maybe not so easy too
		NewNonTips<-popNon$tips+max(popExp$lineages[,1])     #maybe not so easy either
		ancestralTips<-c(popExp$tips, NewNonTips)
		ancestralEvo<-0
	}
	
	return (list(lineages=ancestralLineages,times=ancestralTimes,treeHeight=ancestralHeight,tips=ancestralTips,evo=ancestralEvo))
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

	NtipsExp<-min(evoExp[,1])
	NtipsNon<-min(evoNon[,1])
	#tips<-setdiff(unique(evo[,2]),unique(evo[,1]))
	if (NtipsExp>=NtipsNon){
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
	if(evoR[1,1]<1) {
		evoR<-evoR[-1,]
	}
	evoNew<-matrix(nrow=nrow(evoR),ncol=ncol(evoR),data=0)
	
	if(all(evoR[,1]>=length(tips))){
		evoR[,1]<-evoR[,1]+1

		evoR[!(evoR[,2] %in% tips),2]<-evoR[!(evoR[,2] %in% tips),2]+1

		evoR[evoR[,1]==max(evoR[,1]),1]<-max(tips)+1
		evoNew<-evoR
	}
	else{
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
	}
	
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
 	

TreeMakeBoth<-function(timesExp,timesNon,alpha,Ts,Texp,RatioNow,RatioTexp,expns, bounded){

#if(missing(ratio)){ratio<-((length(timesExp))/(length(timesNon)))}

######SETTING THE STARTING PARAMETERS
#lineages is the matrix containing the extant lineages 
#and the time they appear in the genealogy either by coalescent event or by addition of a lineage
lineagNon<-matrix(nrow=1,ncol=2,data=0)
lineagExp<-matrix(nrow=1,ncol=2,data=0)

# timesNon<-timesNon-min(c(timesNon,timesExp))   #first Leave at 0
# timesExp<-timesExp-min(c(timesExp,timesNon))   # 

lineagNon[1,]<-c(1,timesNon[1]) #it is initialised with one lineage at time 0
lineagExp[1,]<-c(1,timesExp[1]) #it is initialised with one lineage at time 0

timesNon<-timesNon[-1]
#timesExp<-timesExp[-1]
NlNon<-1 # is the numbe of lineages
#NlExp<-1
tipsNon<-1
tipsExp<-1
#total number of lineages in the simulation
#total<-10


heightNon<-0
height<-0
evoNon<-c(0,0,0)
#evoExp<-c(0,0,0)

if (RatioNow>=1){
	ratioE<-RatioNow
	ratioN<-1
	ratioA<-1+RatioTexp
}
else{
	ratioE<-1
	ratioN<-1/RatioNow
	if(RatioTexp>0){
		ratioA<-1+(1/RatioTexp)
	}
	else{ratioA<-ratioN}
}

######preparing evolution
evolvedNon<-list(lineages=as.matrix(lineagNon),
		times=timesNon,
		treeHeight=0,
		tips=tipsNon,
		evo=evoNon,
		alpha=1,
		ratio=ratioN,
		Texp=Texp,
		expns="none")

#evolvedExpB<-list(lineages=as.matrix(lineagExp),
		#times=timesExp,
		#height=0,
		#tips=tipsExp,
		#evo=evoExp)

 evolvedA<-list(lineages=as.matrix(lineagExp),
 		times=0,
 		treeHeight=0,
 		tips=0,
 		evo=0,
 		alpha=1,
		ratio=ratioA,
		Texp=Texp,
		expns="none")

######## evolving!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#done<-TreeMake(timesExp,alpha,ratio,expns)

if (bounded=="TRUE"){
	done<-TreeMake(timesExp,alpha,ratioE,expns)
	while(done$height>Ts){

		done<-TreeMake(timesExp,alpha,ratioE,expns)
		evolvedExp<-done$evolved
		evoExp<-done$evo
	}
}else{
	evoExp<-c(0,0,0)
	evolvedExp<-list(lineages=as.matrix(lineagExp),
		times=timesExp,
		treeHeight=0,
		tips=tipsExp,
		evo=evoExp,
		alpha=alpha,
		ratio=ratioE,
		Texp=Texp,
		expns="exp")


	control=0
	while((evolvedExp$treeHeight<Ts)&(control==0)){
	if((length(evolvedExp$lineages[,1])>1)|(length(evolvedExp$times)>0)){
	evolvedExp<-evolution(lineages=evolvedExp$lineages,
				tips=evolvedExp$tips,
				times=evolvedExp$times,
				treeHeight=evolvedExp$treeHeight,
				alpha=alpha,
				ratio=ratioE, 
				Texp=Texp,  
				expns="exp")

		if(is.matrix(evolvedExp$evo)){
			evoExp<-rbind(evoExp,evolvedExp$evo)
		}
	}
	else{control=1}
}

}



#treeExp<-makephylo(done$evo)
# plot.phylo(treeExp)

control=0
while((evolvedNon$treeHeight<Ts)&(control==0)){
	if((length(evolvedNon$lineages[,1])>1)|(length(evolvedNon$times)>0)){
	evolvedNon<-evolution(lineages=evolvedNon$lineages,
				tips=evolvedNon$tips,
				times=evolvedNon$times,
				treeHeight=evolvedNon$treeHeight,
				alpha=0,
				ratio=ratioN,
				Texp=Texp,
				expns="none")

		if(is.matrix(evolvedNon$evo)){
			evoNon<-rbind(evoNon,evolvedNon$evo)
		}
	}
	else{control=1}
}

#treeNon<-makephylo(evoNon)

evolvedA<-joinevolved(evolvedExp,evolvedNon,Ts)

evoA<-joinevo(evoExp,evoNon)


#while((length(evolved$lineages[,1])+length(evolved$times)>1){	
  while((length(evolvedA$lineages[,1])>1)|(length(evolvedA$times)>0)){
		
	evolvedA<-evolution(lineages=evolvedA$lineages,
				tips=evolvedA$tips,
				times=evolvedA$times,
				treeHeight=evolvedA$treeHeight,
				alpha=0,
				ratio=ratioA,
				Texp=Texp,
				expns="none")
	
	if(is.matrix(evolvedA$evo)){
		evoA<-rbind(evoA,evolvedA$evo)
	}
}

#tree<-makephyloAll(evoA)
#plot.phylo(tree)
#axisPhylo()

return(list(evo=evoA, height=evolvedA$treeHeight,evolved=evolvedA))
#return(tree)
}



MakeTreeAny<-function(timesExp,timesNon,alpha,Ts,Texp,ratio,expns, bounded,ratioPos){

if(missing(timesExp)&&missing(timesNon)) {stop("you need to specify at least times")}	
if(missing(timesExp)||missing(timesNon)) {
	if (missing(timesExp)){times<-timesNon}
	else{times<-timesExp}

if(missing(ratioPos)){#stop("there's just one population")
ratioPos<-"now"}
RatioPos<-tolower(ratioPos)

if((missing(ratio))||(ratio==0)) {ratio<-1}


if((missing(ratio))&(missing(alpha))&(missing(expns))) {
	expns<-"none"
	ratio<-1
	alpha<-1
}	

if (missing(Ts)){bounded<-"FALSE"
Ts<-1000000
}
if (missing(Texp)){Texp<-Ts}

if(((alpha==0)&(expns!="none"))|((alpha!=0)&(expns=="none"))) {stop("your request is incoherent")}
if((missing(alpha))&(expns!="none")) {
	alpha<-0    
}

#if (alpha==0){stop("alpha=0 is nonsense!")}
if (missing(bounded)){bounded<-"FALSE"}
if (bounded=="F"){bounded<-"FALSE"}

if (bounded=="T"){bounded<-"TRUE"}

if ((bounded=="TRUE")&&(missing (Ts))){stop("when is the tree bounded to coalesce?")}


expns1<-tolower(expns)
expns<-expns1

times<-sort(times-min(times))  ##### Set the first time to zero

if (bounded=="TRUE"){
done<-TreeMake(times,alpha,ratio,expns)
	while(done$height>Ts){

		done<-TreeMake(times,alpha,ratio,expns)
		# evolvedExp<-done$evolved
		# evoExp<-done$evo
	}
}	
else{
	done<-TreeMake(times,alpha,ratio,expns)

}

}	

else{


if((missing(ratio))&(missing(alpha))&(missing(expns))&(missing(Texp))) {
	expns<-"none"
	ratio<-1
	alpha<-0
	Texp<-Ts
}	
if((missing(ratio))| (ratio==0)) {ratio<-1}

if(((alpha==0)&(expns!="none"))|((alpha!=0)&(expns=="none"))) {stop("your request is incoherent")}
if((missing(alpha))&(expns!="none")) {
	alpha<-0
}

#if (alpha==0){stop("alpha=0 is nonsense!")}

if (missing(Texp)){Texp<-Ts}

if (missing(bounded)){bounded<-"FALSE"}
if (bounded=="F"){bounded<-"FALSE"}

if (bounded=="T"){bounded<-"TRUE"}

if (missing(ratioPos)){ratioPos<-"now"}
RatioPos<-tolower(ratioPos)
if (expns=="exp"){
	if(RatioPos=="texp"){
		RatioNow=ratio*exp(alpha*Texp)
		RatioTexp=ratio
	}else{
		RatioNow=ratio
		RatioTexp=ratio*exp(-alpha*Texp)
	}
}else{
	if (expns=="none"){
		RatioNow=ratio
		RatioTexp=ratio
	}else{
		if(expns=="lin"){
			if (missing(ratio)){ratio<-1}
			if(RatioPos=="texp"){
				RatioNow=ratio+alpha*Texp
				RatioTexp=ratio
				stop("Unsupported at the moment")
			}else{
				RatioNow=alpha*Texp
				RatioTexp=0
			}
			if (RatioTexp<0){stop("Incoherent: Negative population at Texp")}
		}
		else{
			if(expns=="Xav"){
				RatioNow=ratio
				RatioTexp=0

			}
		}
	}
}

print(paste("Your population is expanding ",expns," Nexp at present is ", RatioNow," Nexp at Texp is ", RatioTexp, "Na is ",(1+RatioTexp),  sep=""))

expns1<-tolower(expns)
expns<-expns1


timesNon<-sort(timesNon-min(c(timesNon,timesExp)))   #first Leave at 0
 timesExp<-sort(timesExp-min(c(timesExp,timesNon)))   # 


done<-TreeMakeBoth(timesExp,timesNon,alpha,Ts,Texp,RatioNow,RatioTexp,expns, bounded)
}

tree<-makephyloAll(done$evo)
return(tree)
}




# timesExp<-c(rep(0,10))
# ####### times at which we sample the leaves in the NON-Exponentially growing pop
# timesNon<-c(rep(0,20))
# ####### the rate at which the Exponentially expanding pop is expanding
# alpha=0
# ############# 1. the time at which the two populations split
# Ts<-1
# ############ 2. the time at which the Exp populations start expanding
# Texp<-1

# tree<-MakeTreeAny(timesExp,timesNon,alpha,Ts,Texp)
# plot.phylo(tree)
# axisPhylo()

# treeB<-MakeTreeAny(timesExp=rep(0,10),timesNon=rep(0,5),alpha=50,Ts=1,Texp=1,ratio=10,expns="exp")
#  plot.phylo(treeB)
# axisPhylo()

timesExp<-rep(0,10)
 timesNon<-rep(0,5)
 alpha<-50
 Ts<-2
 Texp<-2
 ratio<-10
 expns<-"exp"
