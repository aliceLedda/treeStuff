library(ape)
#########FUNCTIONS TO GET PROGRAM WORKING

source("progs/coalescentEvoLL1.R")
source("progs/mergEvo.R")



TreeMakeBoth<-function(timesExp,timesNon,alpha,Ts,Texp,ratio,expns, bounded){

#lineages is the matrix containing the extant lineages 
#and the time they appear in the genealogy either by coalescent event or by addition of a lineage
lineagNon<-matrix(nrow=1,ncol=2,data=0)
lineagExp<-matrix(nrow=1,ncol=2,data=0)

# timesNon<-timesNon-min(c(timesNon,timesExp))   #first Leave at 0
# timesExp<-timesExp-min(c(timesExp,timesNon))   # 

lineagNon[1,]<-c(1,timesNon[1]) #it is initialised with one lineage at time 0
lineagExp[1,]<-c(1,timesExp[1]) #it is initialised with one lineage at time 0

timesNon<-timesNon[-1]
timesExp<-timesExp[-1]
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

	ratioE<-ratio
	ratioN<-1
	ratioA<-1

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
		expns=expns)


	control=0
	while((evolvedExp$treeHeight<Ts)&(control==0)){
	if((length(evolvedExp$lineages[,1])>1)|(length(evolvedExp$times)>0)){
#print(evolvedExp$lineages)
	evolvedExp<-evolution(lineages=evolvedExp$lineages,
				tips=evolvedExp$tips,
				times=evolvedExp$times,
				treeHeight=evolvedExp$treeHeight,
				alpha=alpha,
				ratio=ratioE, 
				Texp=Texp,  
				expns=expns)

		if(is.matrix(evolvedExp$evo)){
			evoExp<-rbind(evoExp,evolvedExp$evo)
		}
	}
	else{control=1}
}

}

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

evolvedA<-mergEvolved(evolvedExp,evolvedNon,Ts,evoExp,evoNon)
evoA<-evolvedA$evoM
evolvedA$evoM<-NULL


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

return(list(evo=evoA, height=evolvedA$treeHeight,evolved=evolvedA))

}



MakeTreeAny<-function(timesExp,timesNon,alpha,Ts,Texp,ratio,expns, bounded,ratioPos){

if(missing(timesExp)&&missing(timesNon)) {stop("you need to specify at least times")}	
if(missing(timesExp)||missing(timesNon)) {
	if (missing(timesExp)){times<-timesNon}
	else{times<-timesExp}
RatioPos<-tolower(ratioPos)
if((RatioPos=="texp")|(RatioPos=="now")){stop("there's just one population")}

if((missing(ratio))&(missing(alpha))&(missing(expns))) {
	expns<-"none"
	ratio<-0
	alpha<-1
}	
if((missing(ratio))| (ratio==0)) {ratio<-1}

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

expns1<-tolower(expns)
expns<-expns1


timesNon<-sort(timesNon-min(c(timesNon,timesExp)))   #first Leave at 0
 timesExp<-sort(timesExp-min(c(timesExp,timesNon)))   # 


done<-TreeMakeBoth(timesExp,timesNon,alpha,Ts,Texp,ratio,expns, bounded)
}

tree<-makePhyloMerged(done$evo)
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

# timesExp<-rep(0,10)
#  timesNon<-rep(0,5)
#  alpha<-50
#  Ts<-2
#  Texp<-2
#  ratio<-10
#  expns<-"exp"

###expA tree
 #treeB<-MakeTreeAny(timesExp=rep(0,50),timesNon=rep(0,30),alpha=50,Ts=0.5,Texp=0.5,ratio=0.1,expns="exp")


# plot.phylo(treeB,show.tip.label=FALSE,edge.width=3)
#  axisPhylo(cex=4,lwd=2)

# write.tree(treeB,file="tmp1.tre")
#  treeB<-read.tree("tmp1.tre")
#  plot.phylo(treeB,show.tip.label=FALSE,edge.width=3)
#  axisPhylo(cex=4,lwd=2)
#  nodelabels()
#  sub<-extract.clade(treeB,110)
#  plot.phylo(sub)
#  hist(sub$edge.length)
#  max(sub$edge.length)

#  hist(sub$edge.length,breaks=18)

Pop.plot<-function(tree,alpha, ratio ,Texp){
	pop<-0
	times<-node.depth.edgelength(tree)
	revTimes<-max(signif(times,digits=6))-signif(times,digits=6)
	for (i in 1:length(revTimes)){
		if (revTimes[i]<Texp){
			popT<-ratio*(tanh(alpha*(Texp-revTimes[i])/2))/(tanh(alpha*Texp/2))
		}else{
			popT<-0
		}
		pop<-c(pop,popT)
	}
	pop<-pop[-1]
	par(mfcol=c(2,1))
	plot.phylo(tree)
	axisPhylo()
	plot(sort(times),sort(pop),type="l",lwd=3, xaxt="n",xlab="times", ylab="pop size")
	 axisPhylo()

}

