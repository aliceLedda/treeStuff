library(ape)

times<-1000   #times is the number of simulations

logli<-matrix(nrow=times, ncol=18)    #logli is the results matrix

total<-100    #total is the number of leaves in the end of the simulation
NEx<-0
changeT<-5
pippo<-1
while(pippo<(times)){
#while((100-changeT)>=NEx){
tip<-c()

tipE<-c()
tipN<-c()

T<-0
t<-0
#t<-rexp(1,rate=1)
evo<-matrix(nrow=2,ncol=3,data=0)
evo[1,]<-c(1,2,t)
evo[2,]<-c(1,3,t)
tip<-c(2,3)
last<-3
nbirths<-total-2
#changeT<-sample(90:(nbirths-3),1)
#changeT<-sample(3:(nbirths-3),1)
#changeT<-25
endT<-runif(1,min=20, max=500) #this part need to be redone. The min and max parameters are just choosen by me, no reason.

changeT<-runif(1,min=20, max=endT)

Nrate<-1
#Erate<-1
Erate<-sample(1:100,1)   #the expansion rate
#total<-100  #total number of individuals

while (T<changeT){

t<-rexp(1,rate=1)
l<-sample(1:length(tip),1)
evo[evo[,2] %in% tip,3]<-evo[evo[,2] %in% tip,3]+t
new<-matrix(nrow=2,ncol=3)
new[,1]<-rep(tip[l], length(new[,1]))
new[,3]<-rep(0, length(new[,1]))
new[1,2]<-last+1
new[2,2]<-last+2
tip<-tip[-l]
tip<-c(tip,last+1,last+2)
#evo[evo[,2] %in% tip,3]<-evo[evo[,2] %in% tip,3]+t
last<-last+2
evo<-rbind(evo,new)
T<-T+t
}

#length(tip)
#I want the mutation to happen at most at 3 times from the end.
if (length(tip)<total-3){

#I want the switch to happen really at changeT, so if the branch length I extracted is larger, I get rid of what's more
DeltaT<-T-changeT
evo[evo[,2] %in% tip,3]<-evo[evo[,2] %in% tip,3]-DeltaT
change<-sample(1:length(tip),1)
cB<-tip[change]
tipE<-tip[change]
tipN<-tip[-change]


#for(j in (changeT+1):nbirths){
for(j in (length(tip)-1):(total-2)){

mut<-(length(tipE)*Erate)+(Nrate*length(tipN))
p<-(length(tipE)*Erate)/((length(tipE)*Erate)+(Nrate*length(tipN)))
t<-rexp(1,rate=mut)
m<-runif(1, min=0,max=1)

#I assign randomly where the split is going to happen, whether in the expanding branch or in the other
if(m<p){
le<-sample(1:length(tipE),1)
l<-which(tip==tipE[le])
tipE<-tipE[-le]
tipE<-c(tipE,last+1,last+2)
}else{lN<-sample(1:length(tipN),1)
l<-which(tip==tipN[lN])
tipN<-tipN[-lN]
tipN<-c(tipN,last+1,last+2)
}

evo[evo[,2] %in% tip,3]<-evo[evo[,2] %in% tip,3]+t
new<-matrix(nrow=2,ncol=3)
new[,1]<-rep(tip[l], length(new[,1]))
new[,3]<-rep(0, length(new[,1]))
new[1,2]<-last+1
new[2,2]<-last+2
tip<-tip[-l]
tip<-c(tip,last+1,last+2)
#evo[evo[,2] %in% tip,3]<-evo[evo[,2] %in% tip,3]+t
last<-last+2
evo<-rbind(evo,new)

}

#now I have to extract the las branch length
tMax<-rexp(1,rate=mut)   # I extract the time to next event
t<-runif(1,min=0,max=tMax)  #I extract a time within this time interval where I "sample"
evo[evo[,2] %in% tip,3]<-evo[evo[,2] %in% tip,3]+t      #update branch length

#now I transform the tree into a phylo object. The tips must have numbers between 1 and length(tip), 
#the root is length(tip)+1 and then I go on numbering the splits in time order.
evoR<-matrix(ncol=length(evo[1,]), nrow=length(evo[,1]))

evoR[,3]<-evo[,3]
evoR[,1]<-evo[,1]+total
evoR[,2]<-evo[,2]+total
# I do not want the exponential growing branch to consist of only one tip
if (cB %in% tip){
}else{
pippo<-pippo+1
cBR<-cB+total
#cB1<-cB+total+1
#cB2<-cB+total-1


tipR<-tip+total
#cB<-change+total

keep<-evoR

leaves<-c()
for (k in 1:length(tipR)){
j<-length(tipR)+1-k
old<-tipR[j]
new<-j
leaves<-c(leaves,j)
evoR[evoR[,1]==old,1]<-j
evoR[evoR[,2]==old,2]<-j
evoR[evoR[,1]>old,1]<-evoR[evoR[,1]>old,1]-1
evoR[evoR[,2]>old,2]<-evoR[evoR[,2]>old,2]-1
if (cBR>old){cBR<-cBR-1}
if(cB2>old){cB2<-cB2-1}
if(cB1>old){cB1<-cB1-1}
}

#now that the matrix is ready I transform it in a phylo object. 
#If you get seg fault, try plot.phylo(alb), if it gives a mistake the problem is that the matrix is not properly declared.
tmp<-evoR[,1:2]
alb<-list(edge=evoR[,1:2])
class(alb)<-"phylo"

attr(alb,"order")<-"postorder"

alb$edge.length<-evoR[,3]
alb$tip.label<-seq(1,total,1)
alb$Nnode<-nbirths+1

#plot.phylo(alb)

# a<-extract.clade(alb,cB+1)
# b<-extract.clade(alb,cB)
# Na<-length(a$tip)
# Nb<-length(b$tip)
# if(Na>Nb){CB<-cB+1
# }else{CB<-cB}
# 
# albEx<-extract.clade(alb,CB)
# albNon<-drop.tip(alb,CB)


#now I divide the exponential growing tree (E) from the tree growing at rate 1 (N)
albEx<-extract.clade(alb,cBR)
albNon<-drop.tip(alb,albEx$tip.label)

#you need the next two commands to be sure it does not make mistakes in node.depth.edgelength
prova<-reorder(alb, order="cladewise")
alb2<-reorder(prova, order="postorder")

u<-evoR[evoR[,2]==cBR,1]
tu<-node.depth.edgelength(alb2)[u]
ts<-node.depth.edgelength(alb2)[cBR]

LNon<-sum(albNon$edge.length)
Nnon<-length(albNon$tip)

#LEx<-node.depth.edgelength(albEx)[1]
LEx<-sum(albEx$edge.length)
NEx<-length(albEx$tip)

L0<-sum(alb$edge.length)
N0<-length(alb$tip)

#plot.phylo(albEx)

#in the test you use the number of splits not the number of tips, so -1
N0<-N0-1
NEx<-NEx-1
Nnon<-Nnon-1



#compute lambdas
lambda0<-N0/L0
#lambdaS<-NEx/(LEx+tu-ts)
lambdaS<-NEx/LEx
lambdaNon<-Nnon/(LNon+tu-ts)
#lambdaNon<-Nnon/(LNon)

#lambdaNon<-Nnon/(LNon)

#compute log likelyhoods
LogLi0<-((-lambda0)*L0)+(N0*log(lambda0))

LogLi<-((-lambdaNon)*LNon)+(Nnon*log(lambdaNon))-(lambdaS*LEx)+(NEx*log(lambdaS))-(lambdaNon*(tu-ts))


#LogLi<-((-lambdaNon)*LNon)+(Nnon*log(lambdaNon))-(lambdaS*LEx)+(NEx*log(lambdaS))-(lambdaNon*(ts-tu))

#compute the difference
DeltaLogLi<-LogLi-LogLi0

#compute the p value
pchi<-1-(pchisq(DeltaLogLi,df=2))

#print and save
logli[pippo,]<-c(Nrate,Erate,changeT,L0,N0,LEx,NEx,LNon,Nnon,tu,ts,lambda0,lambdaS,lambda, LogLi0, LogLi, DeltaLogLi,pchi)
write.table(t(c(Nrate,Erate,changeT,L0,N0,LEx,NEx,LNon,Nnon,tu,ts,lambda0,lambdaS,lambda, LogLi0, LogLi, DeltaLogLi,pchi)),file="lnonProvaNo.txt",append = TRUE, quote = FALSE, sep = "\t" ,row.names=FALSE, col.names=FALSE)
}
}
}


