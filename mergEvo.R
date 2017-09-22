library(ape)

mergEvolved<-function(popExp,popNon, Ts, evoExp,evoNon){
	#if ((identical(popExp$times,numeric(0)))&(identical(popNon$times,numeric(0)))) {
		## might not be needed because length(numeric(0))=0    Check if it works always. If troubles this might be reason
	#popExp<-evolvedExp
	#popNon<-evolvedNon
	totaltips<-sum(length(popExp$tips),length(popNon$tips),length(popNon$times),length(popExp$times))
	root<-totaltips+1
	
	##prepare the expanding part for merging

	newExp<-evoExp
	if(newExp[1,1]<1) {
		newExp<-newExp[-1,]
	}
	#identify tips and internal nodes
	tipsExp<-popExp$tips
	internalExp<-unique(newExp[,1])
	
	newTipsExp<-tipsExp
	newInternalExp<-internalExp+root-min(internalExp)+1

	#create new Evo matrix
	newExp[(newExp[,2] %in% internalExp),2]<-newExp[(newExp[,2] %in% internalExp),2]+root-min(internalExp)+1
	newExp[(newExp[,1] %in% internalExp),1]<-newExp[(newExp[,1] %in% internalExp),1]+root-min(internalExp)+1
	
	
	#update the expanding evolution list
	popExp$lineages[popExp$lineages[,1] %in% internalExp,1]<-popExp$lineages[popExp$lineages[,1] %in% internalExp,1]+root-min(internalExp)+1

	if(is.matrix(popExp$evo)){
		popExp$evo[popExp$evo[,1] %in% internalExp,1]<-popExp$evo[popExp$evo[,1] %in% internalExp,1]+root-min(internalExp)+1
		popExp$evo[popExp$evo[,2] %in% internalExp,2]<-popExp$evo[popExp$evo[,2] %in% internalExp,2]+root-min(internalExp)+1
	}
	
	#information needed to prepare the non expanding prt for merging
	maxTipsExp<-max(newTipsExp)
	maxIntExp<-max(newInternalExp)

	#prepare the non expanding part for merging
	newNon<-evoNon
	if(newNon[1,1]<1) {
		newNon<-newNon[-1,]
	}

	#identify tips and internal nodes
	tipsNon<-popNon$tips
	internalNon<-unique(newNon[,1])

	newTipsNon<-tipsNon+maxTipsExp
	newInternalNon<-internalNon+maxIntExp-min(internalNon)+1

	#create new Evo matrix
	newNon[(newNon[,2] %in% internalNon),2]<-newNon[(newNon[,2] %in% internalNon),2]+maxIntExp-min(internalNon)+1
	newNon[(newNon[,1] %in% internalNon),1]<-newNon[(newNon[,1] %in% internalNon),1]+maxIntExp-min(internalNon)+1
	newNon[(newNon[,2] %in% tipsNon),2]<-newNon[(newNon[,2] %in% tipsNon),2]+maxTipsExp
	
	
	#update the non-expanding evolution list
	popNon$lineages[popNon$lineages[,1] %in% internalNon,1]<-popNon$lineages[popNon$lineages[,1] %in% internalNon,1]+maxIntExp-min(internalNon)+1
	popNon$lineages[popNon$lineages[,1] %in% tipsNon,1]<-popNon$lineages[popNon$lineages[,1] %in% tipsNon,1]+maxTipsExp
	
	if(is.matrix(popNon$evo)){
		popNon$evo[popNon$evo[,1] %in% internalNon,1]<-popNon$evo[popNon$evo[,1] %in% internalNon,1]+maxIntExp-min(internalNon)+1
		popNon$evo[popNon$evo[,2] %in% internalNon,2]<-popNon$evo[popNon$evo[,2] %in% internalNon,2]+maxIntExp-min(internalNon)+1
		popNon$evo[popNon$evo[,1] %in% tipsNon,1]<-popNon$evo[popNon$evo[,1] %in% tipsNon,1]+maxTipsExp
		popNon$evo[popNon$evo[,2] %in% tipsNon,2]<-popNon$evo[popNon$evo[,2] %in% tipsNon,2]+maxTipsExp
	}

	#create the evoA matrix
	evoA<-rbind(newExp,newNon)

	#create the evolutionA list
	ancestralTimes<-sort(c(popExp$times,popNon$times))
	ancestralHeight<-max(popExp$treeHeight,popNon$treeHeight,Ts)   #maybe not so easy too
	ancestralTips<-c(newTipsExp,newTipsNon)
	ancestralLineages<-rbind(popNon$lineages,popExp$lineages)
	ancestralEvo<-rbind(popNon$evo,popExp$evo)


	return (list(lineages=ancestralLineages,times=ancestralTimes,treeHeight=ancestralHeight,tips=ancestralTips,evo=ancestralEvo, evoM=evoA))
}


makePhyloMerged<-function(evoM){
	tips<-setdiff(unique(evoM[,2]),unique(evoM[,1]))
	root<-setdiff(seq((max(tips)+1),max(evoM[,1]),1),unique(evoM[,1]))
	evoT<-evoM
	evoT[evoT[,1]==max(evoT[,1]),1]<-root
	tmp<-evoT[,1:2]
	tree<-list(edge=evoT[,1:2])
	class(tree)<-"phylo"
	
	#attr(tree,"order")<-"postorder"
	#attr(tree,"order")<-"cladewise"
	
	tree$edge.length<-evoT[,3]
	
	#tree$tip.label<-seq(1,length(tips),1)
	tree$tip.label<-tips   #this is the line causing the funny tip names (not consecutive). It can be deleted with no further consequences than the tips getting consecutive numbers names.
	tree$Nnode<-length(unique(evoT[,1]))
	#back<-list(tree=tree,cBR=cBR)
	return(tree)
}
