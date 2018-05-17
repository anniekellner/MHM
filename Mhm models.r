#####################
#Hemoplasma analysis#
#####################
#cadn = California (all four sites) Domestic cat Number tested
#cadi = California Domestic cat Infected number
#cads = California Domestic cat susceptible number
#cadd = California Domestic cat Density
#cadm = California Domestic cat Male proportion
#b = bobcat
#p = puma

###############################################
##Estimation of unknown parameters at CA-SDRC##
###############################################
rm(list=ls())
require(stats4)
setwd("C:/Documents and Settings/scarver/My Documents/CSU/Hemoplasma/Analysis 3-10-11/Analysis 4-6-11/Analysis 24-3-16/Analyses 24-10-17") #Desktop
setwd("C:/Users/Scott Carver/Dcauments/CSU/Hemoplasma/Analysis 3-10-11/Analysis 4-6-11/Analysis 24-3-16/Analyses 24-10-17") #Laptop
setwd("D:/CSU/Hemoplasma/Analysis 3-10-11/Analysis 4-6-11/Analysis 24-3-16/Analyses 24-10-17")
AllSites=read.csv("All.Sites new.csv") ####
attach(AllSites)
AllSites

#######################
#######################
##Transmission models##
#######################
#######################

##NB: The numbering of the models is one off, due to adding Aggressive 
##encounters model on 5-13-11

##################
##Social contact##
##################
##Model 1a
##Direct trans is density dependent, no cross species trans
caTv1a=function(Td,Tb,Tp){    
cadinf=exp(Td*cadd/cadn)/(1+exp(Td*cadd/cadn))
cabinf=exp(Tb*cabd/cabn)/(1+exp(Tb*cabd/cabn))
capinf=exp(Tp*capd/capn)/(1+exp(Tp*capd/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv1a=list(Td=1,Tb=1,Tp=1)
modelcaTv1a=mle(caTv1a,start=guesscaTv1a,method="BFGS")
summary(modelcaTv1a) 
#confint(modelcaTv1a)

coef(modelcaTv1a)
Td=coef(modelcaTv1a)[1]
Tb=coef(modelcaTv1a)[2]
Tp=coef(modelcaTv1a)[3]
D=exp(Td*cadd/cadn)/(1+exp(Td*cadd/cadn))
B=exp(Tb*cabd/cabn)/(1+exp(Tb*cabd/cabn))
P=exp(Tp*capd/capn)/(1+exp(Tp*capd/capn))
m1a=c(D,B,P)
m1a


#########################
##Aggressive encounters##
#########################
##Model 2a1
##Direct trans by male aggression is density dependent,no cross species trans
caTv2a1=function(Tmd,Tmb,Tmp){    
cadinf=exp((Tmd*cadm*cadd)/cadn)/(1+exp((Tmd*cadm*cadd)/cadn))
cabinf=exp((Tmb*cabm*cabd)/cabn)/(1+exp((Tmb*cabm*cabd)/cabn))
capinf=exp((Tmp*capm*capd)/capn)/(1+exp((Tmp*capm*capd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv2a1=list(Tmd=-5,Tmb=60,Tmp=0)
modelcaTv2a1=mle(caTv2a1,start=guesscaTv2a1,method="BFGS")
summary(modelcaTv2a1) 
#confint(modelcaTv2a1)

coef(modelcaTv2a1)
Tmd=coef(modelcaTv2a1)[1]
Tmb=coef(modelcaTv2a1)[2]
Tmp=coef(modelcaTv2a1)[3]
D=exp((Tmd*cadm*cadd)/cadn)/(1+exp((Tmd*cadm*cadd)/cadn))
B=exp((Tmb*cabm*cabd)/cabn)/(1+exp((Tmb*cabm*cabd)/cabn))
P=exp((Tmp*capm*capd)/capn)/(1+exp((Tmp*capm*capd)/capn))
m2a1=c(D,B,P)
m2a1


##########################################
##Social contact + Aggressive encounters##
##########################################
##Model 2a
##Direct trans is density dependent, male aggression,no cross species trans
caTv2a=function(Td,Tb,Tp,Tmd,Tmb,Tmp){    
cadinf=exp((Td*cadd+Tmd*cadm*cadd)/cadn)/(1+exp((Td*cadd+Tmd*cadm*cadd)/cadn))
cabinf=exp((Tb*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tb*cabd+Tmb*cabm*cabd)/cabn))
capinf=exp((Tp*capd+Tmp*capm*capd)/capn)/(1+exp((Tp*capd+Tmp*capm*capd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv2a=list(Td=-10,Tb=-10,Tp=1,Tmd=10,Tmb=200,Tmp=-2.5)
modelcaTv2a=mle(caTv2a,start=guesscaTv2a,method="BFGS")
summary(modelcaTv2a) 
#confint(modelcaTv2a)

coef(modelcaTv2a)
Td=coef(modelcaTv2a)[1]
Tb=coef(modelcaTv2a)[2]
Tp=coef(modelcaTv2a)[3]
Tmd=coef(modelcaTv2a)[4]
Tmb=coef(modelcaTv2a)[5]
Tmp=coef(modelcaTv2a)[6]
D=exp((Td*cadd+Tmd*cadm*cadd)/cadn)/(1+exp((Td*cadd+Tmd*cadm*cadd)/cadn))
B=exp((Tb*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tb*cabd+Tmb*cabm*cabd)/cabn))
P=exp((Tp*capd+Tmp*capm*capd)/capn)/(1+exp((Tp*capd+Tmp*capm*capd)/capn))
m2a=c(D,B,P)
m2a


##############################
##Social contact + Predation##
##############################
##Model 3a
##Direct trans is density dependent, cross species trans from d to p by predation 
caTv3a=function(Td,Tb,Tp,Tprpd){    
cadinf=exp((Td*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tb*cabd)/cabn)/(1+exp((Tb*cabd)/cabn))
capinf=exp((Tp*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tp*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv3a=list(Td=1,Tb=9,Tp=1500,Tprpd=0)
modelcaTv3a=mle(caTv3a,start=guesscaTv3a,method="BFGS")
summary(modelcaTv3a) 
#confint(modelcaTv3a) 

coef(modelcaTv3a)
Td=coef(modelcaTv3a)[1]
Tb=coef(modelcaTv3a)[2]
Tp=coef(modelcaTv3a)[3]
Tprpd=coef(modelcaTv3a)[4]
D=exp((Td*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tb*cabd)/cabn)/(1+exp((Tb*cabd)/cabn))
P=exp((Tp*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tp*capd+Tprpd*cadi/cadn*cadd)/capn))
m3a=c(D,B,P)
m3a


##Model 3b
##Direct trans is density dependent, cross species trans from d to b and p by predation
caTv3b=function(Td,Tb,Tp,Tprbd,Tprpd){    
cadinf=exp((Td*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tb*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tb*cabd+Tprbd*cadi/cadn*cadd)/cabn))
capinf=exp((Tp*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tp*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv3b=list(Td=1,Tb=450,Tp=1200,Tprbd=1,Tprpd=1)
modelcaTv3b=mle(caTv3b,start=guesscaTv3b,method="BFGS")
summary(modelcaTv3b) 
#confint(modelcaTv3b) 

coef(modelcaTv3b)
Td=coef(modelcaTv3b)[1]
Tb=coef(modelcaTv3b)[2]
Tp=coef(modelcaTv3b)[3]
Tprbd=coef(modelcaTv3b)[4]
Tprpd=coef(modelcaTv3b)[5]
D=exp((Td*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tb*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tb*cabd+Tprbd*cadi/cadn*cadd)/cabn))
P=exp((Tp*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tp*capd+Tprpd*cadi/cadn*cadd)/capn))
m3b=c(D,B,P)
m3b


##Model 3c
##Direct trans is density dependent, cross species trans from d and b to p by predation
caTv3c=function(Td,Tb,Tp,Tprpd,Tprpb){    
cadinf=exp((Td*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tb*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tb*cabd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv3c=list(Td=1,Tb=0,Tp=1200,Tprpd=1,Tprpb=-150)
modelcaTv3c=mle(caTv3c,start=guesscaTv3c,method="BFGS")
summary(modelcaTv3c) 
#confint(modelcaTv3c) 

coef(modelcaTv3c)
Td=coef(modelcaTv3c)[1]
Tb=coef(modelcaTv3c)[2]
Tp=coef(modelcaTv3c)[3]
Tprpd=coef(modelcaTv3c)[4]
Tprpb=coef(modelcaTv3c)[5]
D=exp((Td*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tb*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tb*cabd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m3c=c(D,B,P)
m3c


##Model 3d
##Direct trans is density dependent, cross species trans from d to b and p and b to p by predation
caTv3d=function(Td,Tb,Tp,Tprbd,Tprpd,Tprpb){    
cadinf=exp((Td*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tb*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tb*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv3d=list(Td=1,Tb=1,Tp=1000,Tprbd=1,Tprpd=1,Tprpb=1)
modelcaTv3d=mle(caTv3d,start=guesscaTv3d,method="BFGS")
summary(modelcaTv3d) 
#confint(modelcaTv3d) 

coef(modelcaTv3d)
Td=coef(modelcaTv3d)[1]
Tb=coef(modelcaTv3d)[2]
Tp=coef(modelcaTv3d)[3]
Tprbd=coef(modelcaTv3d)[4]
Tprpd=coef(modelcaTv3d)[5]
Tprpb=coef(modelcaTv3d)[6]
D=exp((Td*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tb*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tb*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m3d=c(D,B,P)
m3d


######################################################
##Social contact + Aggressive encounters + Predation##
######################################################
##Model 4a
##Direct trans is density dependent, male aggression, cross species trans from d to p by predation
caTv4a=function(Td,Tb,Tp,Tprpd,Tmd,Tmb,Tmp){    
cadinf=exp((Td*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tb*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tb*cabd+Tmb*cabm*cabd)/cabn))
capinf=exp((Tp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv4a=list(Td=-10,Tb=-50,Tp=1300,Tprpd=1,Tmd=10,Tmb=1300,Tmp=1)
modelcaTv4a=mle(caTv4a,start=guesscaTv4a,method="BFGS")
summary(modelcaTv4a) 
#confint(modelcaTv4a) 

coef(modelcaTv4a)
Td=coef(modelcaTv4a)[1]
Tb=coef(modelcaTv4a)[2]
Tp=coef(modelcaTv4a)[3]
Tprpd=coef(modelcaTv4a)[4]
Tmd=coef(modelcaTv4a)[5]
Tmb=coef(modelcaTv4a)[6]
Tmp=coef(modelcaTv4a)[7]
D=exp((Td*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tb*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tb*cabd+Tmb*cabm*cabd)/cabn))
P=exp((Tp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
m4a=c(D,B,P)
m4a


##Model 4b
##Direct trans is density dependent, male aggression, cross species trans from d to b and p by predation
caTv4b=function(Td,Tb,Tp,Tprbd,Tprpd,Tmd,Tmb,Tmp){    
cadinf=exp((Td*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn))
capinf=exp((Tp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv4b=list(Td=-10,Tb=-230,Tp=740,Tprbd=-10,Tprpd=1,Tmd=10,Tmb=670,Tmp=60)
modelcaTv4b=mle(caTv4b,start=guesscaTv4b,method="BFGS")
summary(modelcaTv4b) 
#confint(modelcaTv4b) 

coef(modelcaTv4b)
Td=coef(modelcaTv4b)[1]
Tb=coef(modelcaTv4b)[2]
Tp=coef(modelcaTv4b)[3]
Tprbd=coef(modelcaTv4b)[4]
Tprpd=coef(modelcaTv4b)[5]
Tmd=coef(modelcaTv4b)[6]
Tmb=coef(modelcaTv4b)[7]
Tmp=coef(modelcaTv4b)[8]
D=exp((Td*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn))
P=exp((Tp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
m4b=c(D,B,P)
m4b


##Model 4c
##Direct trans is density dependent, male aggression, cross species trans from d and b to p by predation
caTv4c=function(Td,Tb,Tp,Tmd,Tmb,Tmp,Tprpd,Tprpb){    
cadinf=exp((Td*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tb*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tb*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv4c=list(Td=-10,Tb=-160,Tp=60,Tmd=15,Tmb=-45,Tmp=2100,Tprpd=8,Tprpb=-350)
modelcaTv4c=mle(caTv4c,start=guesscaTv4c,method="BFGS")
summary(modelcaTv4c) 
#confint(modelcaTv4c) 

coef(modelcaTv4c)
Td=coef(modelcaTv4c)[1]
Tb=coef(modelcaTv4c)[2]
Tp=coef(modelcaTv4c)[3]
Tmd=coef(modelcaTv4c)[4]
Tmb=coef(modelcaTv4c)[5]
Tmp=coef(modelcaTv4c)[6]
Tprpd=coef(modelcaTv4c)[7]
Tprpb=coef(modelcaTv4c)[8]
D=exp((Td*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tb*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tb*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m4c=c(D,B,P)
m4c


##Model 4d
##Direct trans is density dependent, male aggression, cross species trans from d to b and p and b to p by predation
caTv4d=function(Td,Tb,Tp,Tprbd,Tprpbd,Tmd,Tmb,Tmp){    
cadinf=exp((Td*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpbd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpbd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpbd*cabi/cabn*cabd)/cabn)/(1+exp((Tb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpbd*cabi/cabn*cabd)/cabn))
capinf=exp((Tp*capd+Tmp*capm*capd+Tprpbd*cadi/cadn*cadd+Tprpbd*cabi/cabn*cabd)/capn)/(1+exp((Tp*capd+Tmp*capm*capd+Tprpbd*cadi/cadn*cadd+Tprpbd*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv4d=list(Td=-10,Tb=-250,Tp=840,Tprbd=-10,Tprpbd=1,Tmd=10,Tmb=690,Tmp=60)
modelcaTv4d=mle(caTv4d,start=guesscaTv4d,method="BFGS")
summary(modelcaTv4d) 
#confint(modelcaTv4d) 

coef(modelcaTv4d)
Td=coef(modelcaTv4d)[1]
Tb=coef(modelcaTv4d)[2]
Tp=coef(modelcaTv4d)[3]
Tprbd=coef(modelcaTv4d)[4]
Tprpbd=coef(modelcaTv4d)[5]
Tmd=coef(modelcaTv4d)[6]
Tmb=coef(modelcaTv4d)[7]
Tmp=coef(modelcaTv4d)[8]
D=exp((Td*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpbd*cadi/cadn*cadd)/cadn)/(1+exp((Td*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpbd*cadi/cadn*cadd)/cadn))
B=exp((Tb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpbd*cabi/cabn*cabd)/cabn)/(1+exp((Tb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpbd*cabi/cabn*cabd)/cabn))
P=exp((Tp*capd+Tmp*capm*capd+Tprpbd*cadi/cadn*cadd+Tprpbd*cabi/cabn*cabd)/capn)/(1+exp((Tp*capd+Tmp*capm*capd+Tprpbd*cadi/cadn*cadd+Tprpbd*cabi/cabn*cabd)/capn))
m4d=c(D,B,P)
m4d


################
##Vector-borne##
################
##Model 5a
##vector trans, population density is important, Tv is uniform
caTv5a=function(Tv){    
cadinf=exp((Tv*cadd)/cadn)/(1+exp((Tv*cadd)/cadn))
cabinf=exp((Tv*cabd)/cabn)/(1+exp((Tv*cabd)/cabn))
capinf=exp((Tv*capd)/capn)/(1+exp((Tv*capd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv5a=list(Tv=1)
modelcaTv5a=mle(caTv5a,start=guesscaTv5a,method="BFGS")
summary(modelcaTv5a) 
#confint(modelcaTv5a) 

coef(modelcaTv5a)
Tv=coef(modelcaTv5a)[1]
D=exp((Tv*cadd)/cadn)/(1+exp((Tv*cadd)/cadn))
B=exp((Tv*cabd)/cabn)/(1+exp((Tv*cabd)/cabn))
P=exp((Tv*capd)/capn)/(1+exp((Tv*capd)/capn))
m5a=c(D,B,P)
m5a


##Model 5b 
##vector trans, population density is important, b and p aquire vectors from d
caTv5b=function(Tvd,Tvb,Tvp){    
cadinf=exp((Tvd*cadd)/cadn)/(1+exp((Tvd*cadd)/cadn))
cabinf=exp(((Tvb+Tvd)*cabd)/cabn)/(1+exp(((Tvb+Tvd)*cabd)/cabn))
capinf=exp(((Tvp+Tvd)*capd)/capn)/(1+exp(((Tvp+Tvd)*capd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv5b=list(Tvd=1,Tvb=1,Tvp=1)
modelcaTv5b=mle(caTv5b,start=guesscaTv5b,method="BFGS")
summary(modelcaTv5b) 
#confint(modelcaTv5b)

coef(modelcaTv5b)
Tvd=coef(modelcaTv5b)[1]
Tvb=coef(modelcaTv5b)[2]
Tvp=coef(modelcaTv5b)[3]
D=exp((Tvd*cadd)/cadn)/(1+exp((Tvd*cadd)/cadn))
B=exp(((Tvb+Tvd)*cabd)/cabn)/(1+exp(((Tvb+Tvd)*cabd)/cabn))
P=exp(((Tvp+Tvd)*capd)/capn)/(1+exp(((Tvp+Tvd)*capd)/capn))
m5b=c(D,B,P)
m5b


##Model 5c 
##vector trans, population density is important, b and p share vectors and aquire vectors from d also
caTv5c=function(Tvd,Tvbp){    #not able to split Tvbp on this model
cadinf=exp((Tvd*cadd)*cadn/cadn)/(1+exp((Tvd*cadd)*cadn/cadn))
cabinf=exp(((Tvd+Tvbp)*cabd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd)/cabn))
capinf=exp(((Tvd+Tvbp)*capd)/capn)/(1+exp(((Tvd+Tvbp)*capd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv5c=list(Tvd=1,Tvbp=1)
modelcaTv5c=mle(caTv5c,start=guesscaTv5c,method="BFGS")
summary(modelcaTv5c) 
#confint(modelcaTv5c)

coef(modelcaTv5c)
Tvd=coef(modelcaTv5c)[1]
Tvbp=coef(modelcaTv5c)[2]
D=exp((Tvd*cadd)*cadn/cadn)/(1+exp((Tvd*cadd)*cadn/cadn))
B=exp(((Tvd+Tvbp)*cabd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd)/cabn))
P=exp(((Tvd+Tvbp)*capd)/capn)/(1+exp(((Tvd+Tvbp)*capd)/capn))
m5c=c(D,B,P)
m5c


##Model 5d
##vector trans, population density is important, b and p share vectors, d aquires vectors from b and p
caTv5d=function(Tvd,Tvb,Tvp){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvb+Tvp)*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd)/cadn))
cabinf=exp(((Tvb)*cabd)/cabn)/(1+exp(((Tvb)*cabd)/cabn))
capinf=exp(((Tvp)*capd)/capn)/(1+exp(((Tvp)*capd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv5d=list(Tvd=-250,Tvb=1,Tvp=0)
modelcaTv5d=mle(caTv5d,start=guesscaTv5d,method="BFGS")
summary(modelcaTv5d) 
#confint(modelcaTv5d) 

coef(modelcaTv5d)
Tvd=coef(modelcaTv5d)[1]
Tvb=coef(modelcaTv5d)[2]
Tvp=coef(modelcaTv5d)[3]
D=exp(((Tvd+Tvb+Tvp)*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd)/cadn))
B=exp(((Tvb)*cabd)/cabn)/(1+exp(((Tvb)*cabd)/cabn))
P=exp(((Tvp)*capd)/capn)/(1+exp(((Tvp)*capd)/capn))
m5d=c(D,B,P)
m5d


##Model 5e
##vector trans, population density is important, b and p share vectors, d aquires vectors from b and p
caTv5e=function(Tvd,Tvbp){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvbp)*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd)/cadn))
cabinf=exp(((Tvbp)*cabd)/cabn)/(1+exp(((Tvbp)*cabd)/cabn))
capinf=exp(((Tvbp)*capd)/capn)/(1+exp(((Tvbp)*capd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv5e=list(Tvd=1,Tvbp=1)
modelcaTv5e=mle(caTv5e,start=guesscaTv5e,method="BFGS")
summary(modelcaTv5e) 
#confint(modelcaTv5e) 

coef(modelcaTv5e)
Tvd=coef(modelcaTv5e)[1]
Tvbp=coef(modelcaTv5e)[2]
D=exp(((Tvd+Tvbp)*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd)/cadn))
B=exp(((Tvbp)*cabd)/cabn)/(1+exp(((Tvbp)*cabd)/cabn))
P=exp(((Tvbp)*capd)/capn)/(1+exp(((Tvbp)*capd)/capn))
m5e=c(D,B,P)
m5e


########################################
##Vector-borne + Aggressive encounters##
########################################
##Model 6a
##vector trans, population density is important, male aggression, Tv is uniform
caTv6a=function(Tv,Tmd,Tmb,Tmp){    
cadinf=exp((Tv*cadd+Tmd*cadm*cadd)/cadn)/(1+exp((Tv*cadd+Tmd*cadm*cadd)/cadn))
cabinf=exp((Tv*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tv*cabd+Tmb*cabm*cabd)/cabn))
capinf=exp((Tv*capd+Tmp*capm*capd)/capn)/(1+exp((Tv*capd+Tmp*capm*capd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv6a=list(Tv=1,Tmd=1,Tmb=1,Tmp=1)
modelcaTv6a=mle(caTv6a,start=guesscaTv6a,method="BFGS")
summary(modelcaTv6a) 
#confint(modelcaTv6a) 

coef(modelcaTv6a)
Tv=coef(modelcaTv6a)[1]
Tmd=coef(modelcaTv6a)[2]
Tmb=coef(modelcaTv6a)[3]
Tmp=coef(modelcaTv6a)[4]
D=exp((Tv*cadd+Tmd*cadm*cadd)/cadn)/(1+exp((Tv*cadd+Tmd*cadm*cadd)/cadn))
B=exp((Tv*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tv*cabd+Tmb*cabm*cabd)/cabn))
P=exp((Tv*capd+Tmp*capm*capd)/capn)/(1+exp((Tv*capd+Tmp*capm*capd)/capn))
m6a=c(D,B,P)
m6a


##Model 6b
##vector trans, population density is important, male aggression, b and p aquire vectors from d
caTv6b=function(Tvd,Tvb,Tvp,Tmd,Tmb,Tmp){    
cadinf=exp((Tvd*cadd+Tmd*cadm*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd)/cadn))
cabinf=exp(((Tvb+Tvd)*cabd+Tmb*cabm*cabd)/cabn)/(1+exp(((Tvb+Tvd)*cabd+Tmb*cabm*cabd)/cabn))
capinf=exp(((Tvp+Tvd)*capd+Tmp*capm*capd)/capn)/(1+exp(((Tvp+Tvd)*capd+Tmp*capm*capd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv6b=list(Tvd=-1,Tvb=1,Tvp=162,Tmd=1,Tmb=1110,Tmp=-100)
modelcaTv6b=mle(caTv6b,start=guesscaTv6b,method="BFGS")
summary(modelcaTv6b) 
#confint(modelcaTv6b) 

coef(modelcaTv6b)
Tvd=coef(modelcaTv6b)[1]
Tvb=coef(modelcaTv6b)[2]
Tvp=coef(modelcaTv6b)[3]
Tmd=coef(modelcaTv6b)[4]
Tmb=coef(modelcaTv6b)[5]
Tmp=coef(modelcaTv6b)[6]
D=exp((Tvd*cadd+Tmd*cadm*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd)/cadn))
B=exp(((Tvb+Tvd)*cabd+Tmb*cabm*cabd)/cabn)/(1+exp(((Tvb+Tvd)*cabd+Tmb*cabm*cabd)/cabn))
P=exp(((Tvp+Tvd)*capd+Tmp*capm*capd)/capn)/(1+exp(((Tvp+Tvd)*capd+Tmp*capm*capd)/capn))
m6b=c(D,B,P)
m6b


##Model 6c
##vector trans, population density is important, male aggression, b and p share vectors and aquire vectors from d also
caTv6c=function(Tvd,Tvbp,Tmd,Tmb,Tmp){    #not able to split Tvbp on this model
cadiv=exp((Tvd*cadd+Tmd*cadm*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd)/cadn))
cabinf=exp(((Tvbp+Tvd)*cabd+Tmb*cabm*cabd)/cabn)/(1+exp(((Tvbp+Tvd)*cabd+Tmb*cabm*cabd)/cabn))
capinf=exp(((Tvbp+Tvd)*capd+Tmp*capm*capd)/capn)/(1+exp(((Tvbp+Tvd)*capd+Tmp*capm*capd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv6c=list(Tvd=1,Tvbp=160,Tmd=1,Tmb=1110,Tmp=-100)
modelcaTv6c=mle(caTv6c,start=guesscaTv6c,method="BFGS") #unable to run successfully !!!!!!!!!!!!!! Make same as 6a
summary(modelcaTv6c) 
#confint(modelcaTv6c) 

#coef(modelcaTv6c)
#Tvd=coef(modelcaTv6c)[1]
#  Tvbp=coef(modelcaTv6c)[2]
#  Tmd=coef(modelcaTv6c)[3]
#  Tmb=coef(modelcaTv6c)[4]
#  Tmp=coef(modelcaTv6c)[5]
#  D=exp((Tvd*cadd+Tmd*cadm*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd)/cadn))
#B=exp(((Tvbp+Tvd)*cabd+Tmb*cabm*cabd)/cabn)/(1+exp(((Tvbp+Tvd)*cabd+Tmb*cabm*cabd)/cabn))
#P=exp(((Tvbp+Tvd)*capd+Tmp*capm*capd)/capn)/(1+exp(((Tvbp+Tvd)*capd+Tmp*capm*capd)/capn))
#m6c=c(D,B,P)
m6c=m6a
m6c


##Model 6d
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv6d=function(Tvd,Tvb,Tvp,Tmd,Tmb,Tmp){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd)/cadn))
cabinf=exp((Tvb*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tvb*cabd+Tmb*cabm*cabd)/cabn))
capinf=exp((Tvp*capd+Tmp*capm*capd)/capn)/(1+exp((Tvp*capd+Tmp*capm*capd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv6d=list(Tvd=40,Tvb=-250,Tvp=200,Tmd=10,Tmb=1100,Tmp=-410)
modelcaTv6d=mle(caTv6d,start=guesscaTv6d,method="BFGS")
summary(modelcaTv6d) 
#confint(modelcaTv6d) 

coef(modelcaTv6d)
Tvd=coef(modelcaTv6d)[1]
Tvb=coef(modelcaTv6d)[2]
Tvp=coef(modelcaTv6d)[3]
Tmd=coef(modelcaTv6d)[4]
Tmb=coef(modelcaTv6d)[5]
Tmp=coef(modelcaTv6d)[6]
D=exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd)/cadn))
B=exp((Tvb*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tvb*cabd+Tmb*cabm*cabd)/cabn))
P=exp((Tvp*capd+Tmp*capm*capd)/capn)/(1+exp((Tvp*capd+Tmp*capm*capd)/capn))
m6d=c(D,B,P)
m6d


##Model 6e
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv6e=function(Tvd,Tvbp,Tmd,Tmb,Tmp){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd)/cadn))
cabinf=exp((Tvbp*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tvbp*cabd+Tmb*cabm*cabd)/cabn))
capinf=exp((Tvbp*capd+Tmp*capm*capd)/capn)/(1+exp((Tvbp*capd+Tmp*capm*capd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv6e=list(Tvd=50,Tvbp=-60,Tmd=11,Tmb=150,Tmp=1)
modelcaTv6e=mle(caTv6e,start=guesscaTv6e,method="BFGS")
summary(modelcaTv6e) 
#confint(modelcaTv6e) 

coef(modelcaTv6e)
Tvd=coef(modelcaTv6e)[1]
Tvbp=coef(modelcaTv6e)[2]
Tmd=coef(modelcaTv6e)[3]
Tmb=coef(modelcaTv6e)[4]
Tmp=coef(modelcaTv6e)[5]
D=exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd)/cadn))
B=exp((Tvbp*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tvbp*cabd+Tmb*cabm*cabd)/cabn))
P=exp((Tvbp*capd+Tmp*capm*capd)/capn)/(1+exp((Tvbp*capd+Tmp*capm*capd)/capn))
m6e=c(D,B,P)
m6e


############################
##Vector-borne + Predation##
############################
##Model 7a1 
##vector trans, population density is important, male aggression, cross species trans from d to p by predation, Tv is uniform
caTv7a1=function(Tv,Tprpd){    
cadinf=exp((Tv*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tv*cabd)/cabn)/(1+exp((Tv*cabd)/cabn))
capinf=exp((Tv*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tv*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7a1=list(Tv=0.1,Tprpd=1)
modelcaTv7a1=mle(caTv7a1,start=guesscaTv7a1,method="BFGS")
summary(modelcaTv7a1) 
#confint(modelcaTv7a1)

coef(modelcaTv7a1)
Tv=coef(modelcaTv7a1)[1]
Tprpd=coef(modelcaTv7a1)[2]
D=exp((Tv*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tv*cabd)/cabn)/(1+exp((Tv*cabd)/cabn))
P=exp((Tv*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tv*capd+Tprpd*cadi/cadn*cadd)/capn))
m7a1=c(D,B,P)
m7a1


##Model 7b1
##vector trans, population density is important, male aggression, cross species trans from d to b and p by predation, Tv is uniform
caTv7b1=function(Tv,Tprbd,Tprpd){    
cadinf=exp((Tv*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tv*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tv*cabd+Tprbd*cadi/cadn*cadd)/cabn))
capinf=exp((Tv*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tv*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7b1=list(Tv=1,Tprbd=1,Tprpd=1)
modelcaTv7b1=mle(caTv7b1,start=guesscaTv7b1,method="BFGS")
summary(modelcaTv7b1) 
#confint(modelcaTv7b1) 

coef(modelcaTv7b1)
Tv=coef(modelcaTv7b1)[1]
Tprbd=coef(modelcaTv7b1)[2]
Tprpd=coef(modelcaTv7b1)[3]
D=exp((Tv*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tv*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tv*cabd+Tprbd*cadi/cadn*cadd)/cabn))
P=exp((Tv*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tv*capd+Tprpd*cadi/cadn*cadd)/capn))
m7b1=c(D,B,P)
m7b1


##Model 7c1
##vector trans, population density is important, male aggression, cross species trans from d and b to p by predation, Tv is uniform
caTv7c1=function(Tv,Tprpd,Tprpb){    
cadinf=exp((Tv*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tv*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tv*cabd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tv*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tv*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7c1=list(Tv=1,Tprpd=1,Tprpb=1)
modelcaTv7c1=mle(caTv7c1,start=guesscaTv7c1,method="BFGS")
summary(modelcaTv7c1) 
#confint(modelcaTv7c1) 

coef(modelcaTv7c1)
Tv=coef(modelcaTv7c1)[1]
Tprpd=coef(modelcaTv7c1)[2]
Tprpb=coef(modelcaTv7c1)[3]
D=exp((Tv*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tv*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tv*cabd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tv*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tv*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m7c1=c(D,B,P)
m7c1


##Model 7d1
##vector trans, population density is important, male aggression, cross species trans from d to b and p and b to p by predation, Tv is uniform
caTv7d1=function(Tv,Tprbd,Tprpd,Tprpb){    
cadinf=exp((Tv*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tv*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tv*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tv*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tv*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7d1=list(Tv=1,Tprbd=1,Tprpd=1,Tprpb=1)
modelcaTv7d1=mle(caTv7d1,start=guesscaTv7d1,method="BFGS")
summary(modelcaTv7d1) 
#confint(modelcaTv7d1) 

coef(modelcaTv7d1)
Tv=coef(modelcaTv7d1)[1]
Tprbd=coef(modelcaTv7d1)[2]
Tprpd=coef(modelcaTv7d1)[3]
Tprpb=coef(modelcaTv7d1)[4]
D=exp((Tv*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tv*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tv*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tv*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tv*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m7d1=c(D,B,P)
m7d1


##Model 7a2 
##vector trans, population density is important, male aggression, cross species trans from d to p by predation, b and p aquire vectors from d
caTv7a2=function(Tvd,Tvb,Tvp,Tprpd){    
cadinf=exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvb)*cabd)/cabn)/(1+exp(((Tvd+Tvb)*cabd)/cabn))
capinf=exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7a2=list(Tvd=1,Tvb=1,Tvp=1725,Tprpd=1)
modelcaTv7a2=mle(caTv7a2,start=guesscaTv7a2,method="BFGS")
summary(modelcaTv7a2) 
#confint(modelcaTv7a2)

coef(modelcaTv7a2)
Tvd=coef(modelcaTv7a2)[1]
Tvb=coef(modelcaTv7a2)[2]
Tvp=coef(modelcaTv7a2)[3]
Tprpd=coef(modelcaTv7a2)[4]
D=exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvb)*cabd)/cabn)/(1+exp(((Tvd+Tvb)*cabd)/cabn))
P=exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd)/capn))
m7a2=c(D,B,P)
m7a2


##Model 7b2
##vector trans, population density is important, male aggression, cross species trans from d to b and p by predation, b and p aquire vectors from d
caTv7b2=function(Tvd,Tvb,Tvp,Tprbd,Tprpd){    
cadinf=exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvb)*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp(((Tvd+Tvb)*cabd+Tprbd*cadi/cadn*cadd)/cabn))
capinf=exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7b2=list(Tvd=1,Tvb=1,Tvp=1300,Tprbd=1,Tprpd=1)
modelcaTv7b2=mle(caTv7b2,start=guesscaTv7b2,method="BFGS")
summary(modelcaTv7b2) 
#confint(modelcaTv7b2) 

coef(modelcaTv7b2)
Tvd=coef(modelcaTv7b2)[1]
Tvb=coef(modelcaTv7b2)[2]
Tvp=coef(modelcaTv7b2)[3]
Tprbd=coef(modelcaTv7b2)[4]
Tprpd=coef(modelcaTv7b2)[5]
D=exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvb)*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp(((Tvd+Tvb)*cabd+Tprbd*cadi/cadn*cadd)/cabn))
P=exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd)/capn))
m7b2=c(D,B,P)
m7b2


##Model 7c2
##vector trans, population density is important, male aggression, cross species trans from d and b to p by predation, b and p aquire vectors from d
caTv7c2=function(Tvd,Tvb,Tvp,Tprpd,Tprpb){    
cadinf=exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvb)*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvb)*cabd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7c2=list(Tvd=-1,Tvb=-85,Tvp=1500,Tprpd=3,Tprpb=-195)
modelcaTv7c2=mle(caTv7c2,start=guesscaTv7c2,method="BFGS")
summary(modelcaTv7c2) 
#confint(modelcaTv7c2) 

coef(modelcaTv7c2)
Tvd=coef(modelcaTv7c2)[1]
Tvb=coef(modelcaTv7c2)[2]
Tvp=coef(modelcaTv7c2)[3]
Tprpd=coef(modelcaTv7c2)[4]
Tprpb=coef(modelcaTv7c2)[5]
D=exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvb)*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvb)*cabd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m7c2=c(D,B,P)
m7c2


##Model 7d2
##vector trans, population density is important, male aggression, cross species trans from d to b and p and b to p by predation, b and p aquire vectors from d
caTv7d2=function(Tvd,Tvb,Tvp,Tprbd,Tprpd,Tprpb){    
cadinf=exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvb)*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvb)*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7d2=list(Tvd=1,Tvb=300,Tvp=1000,Tprbd=1,Tprpd=1,Tprpb=-250)
modelcaTv7d2=mle(caTv7d2,start=guesscaTv7d2,method="BFGS")
summary(modelcaTv7d2) 
#confint(modelcaTv7d2) 

coef(modelcaTv7d2)
Tvd=coef(modelcaTv7d2)[1]
Tvb=coef(modelcaTv7d2)[2]
Tvp=coef(modelcaTv7d2)[3]
Tprbd=coef(modelcaTv7d2)[4]
Tprpd=coef(modelcaTv7d2)[5]
Tprpb=coef(modelcaTv7d2)[6]
D=exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvb)*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvb)*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m7d2=c(D,B,P)
m7d2


##Model 7a3 
##vector trans, population density is important, male aggression, cross species trans from d to p by predation, b and p aquire vectors from d
caTv7a3=function(Tvd,Tvbp,Tprpd){    
cadinf=exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvbp)*cabd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd)/cabn))
capinf=exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7a3=list(Tvd=1,Tvbp=1,Tprpd=1)
modelcaTv7a3=mle(caTv7a3,start=guesscaTv7a3,method="BFGS")
summary(modelcaTv7a3) 
#confint(modelcaTv7a3)

coef(modelcaTv7a3)
Tvd=coef(modelcaTv7a3)[1]
Tvbp=coef(modelcaTv7a3)[2]
Tprpd=coef(modelcaTv7a3)[3]
D=exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvbp)*cabd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd)/cabn))
P=exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd)/capn))
m7a3=c(D,B,P)
m7a3


##Model 7b3
##vector trans, population density is important, male aggression, cross species trans from d to b and p by predation, b and p aquire vectors from d
caTv7b3=function(Tvd,Tvbp,Tprbd,Tprpd){    
cadinf=exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvbp)*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd+Tprbd*cadi/cadn*cadd)/cabn))
capinf=exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7b3=list(Tvd=1,Tvbp=1,Tprbd=1,Tprpd=1)
modelcaTv7b3=mle(caTv7b3,start=guesscaTv7b3,method="BFGS")
summary(modelcaTv7b3) 
#confint(modelcaTv7b3) 

coef(modelcaTv7b3)
Tvd=coef(modelcaTv7b3)[1]
Tvbp=coef(modelcaTv7b3)[2]
Tprbd=coef(modelcaTv7b3)[3]
Tprpd=coef(modelcaTv7b3)[4]
D=exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvbp)*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd+Tprbd*cadi/cadn*cadd)/cabn))
P=exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd)/capn))
m7b3=c(D,B,P)
m7b3


##Model 7c3
##vector trans, population density is important, male aggression, cross species trans from d and b to p by predation, b and p aquire vectors from d
caTv7c3=function(Tvd,Tvbp,Tprpd,Tprpb){    
cadinf=exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvbp)*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7c3=list(Tvd=1,Tvbp=1,Tprpd=1,Tprpb=-250)
modelcaTv7c3=mle(caTv7c3,start=guesscaTv7c3,method="BFGS")
summary(modelcaTv7c3) 
#confint(modelcaTv7c3) 

coef(modelcaTv7c3)
Tvd=coef(modelcaTv7c3)[1]
Tvbp=coef(modelcaTv7c3)[2]
Tprpd=coef(modelcaTv7c3)[3]
Tprpb=coef(modelcaTv7c3)[4]
D=exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvbp)*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m7c3=c(D,B,P)
m7c3


##Model 7d3
##vector trans, population density is important, male aggression, cross species trans from d to b and p and b to p by predation, b and p aquire vectors from d
caTv7d3=function(Tvd,Tvbp,Tprbd,Tprpd,Tprpb){    
cadinf=exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvbp)*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7d3=list(Tvd=-5,Tvbp=295,Tprbd=-15,Tprpd=5,Tprpb=-250)
modelcaTv7d3=mle(caTv7d3,start=guesscaTv7d3,method="BFGS")
summary(modelcaTv7d3) 
#confint(modelcaTv7d3) 

coef(modelcaTv7d3)
Tvd=coef(modelcaTv7d3)[1]
Tvbp=coef(modelcaTv7d3)[2]
Tprbd=coef(modelcaTv7d3)[3]
Tprpd=coef(modelcaTv7d3)[4]
Tprpb=coef(modelcaTv7d3)[5]
D=exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvbp)*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m7d3=c(D,B,P)
m7d3


##Model 7a4
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv7a4=function(Tvd,Tvb,Tvp,Tprpd){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvb+Tvp)*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvb*cabd)/cabn)/(1+exp((Tvb*cabd)/cabn))
capinf=exp((Tvp*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvp*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7a4=list(Tvd=-500,Tvb=1,Tvp=500,Tprpd=1)
modelcaTv7a4=mle(caTv7a4,start=guesscaTv7a4,method="BFGS")
summary(modelcaTv7a4) 
#confint(modelcaTv7a4) 

coef(modelcaTv7a4)
Tvd=coef(modelcaTv7a4)[1]
Tvb=coef(modelcaTv7a4)[2]
Tvp=coef(modelcaTv7a4)[3]
Tprpd=coef(modelcaTv7a4)[4]
D=exp(((Tvd+Tvb+Tvp)*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvb*cabd)/cabn)/(1+exp((Tvb*cabd)/cabn))
P=exp((Tvp*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvp*capd+Tprpd*cadi/cadn*cadd)/capn))
m7a4=c(D,B,P)
m7a4


##Model 7b4
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv7b4=function(Tvd,Tvb,Tvp,Tprbd,Tprpd){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvb+Tvp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvb*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tvb*cabd+Tprbd*cadi/cadn*cadd)/cabn))
capinf=exp((Tvp*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvp*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7b4=list(Tvd=-990,Tvb=400,Tvp=600,Tprbd=1,Tprpd=1)
modelcaTv7b4=mle(caTv7b4,start=guesscaTv7b4,method="BFGS")
summary(modelcaTv7b4) 
#confint(modelcaTv7b4) 

coef(modelcaTv7b4)
Tvd=coef(modelcaTv7b4)[1]
Tvb=coef(modelcaTv7b4)[2]
Tvp=coef(modelcaTv7b4)[3]
Tprbd=coef(modelcaTv7b4)[4]
Tprpd=coef(modelcaTv7b4)[5]
D=exp(((Tvd+Tvb+Tvp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvb*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tvb*cabd+Tprbd*cadi/cadn*cadd)/cabn))
P=exp((Tvp*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvp*capd+Tprpd*cadi/cadn*cadd)/capn))
m7b4=c(D,B,P)
m7b4


##Model 7c4
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv7c4=function(Tvd,Tvb,Tvp,Tprpd,Tprpb){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvb+Tvp)*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvb*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvb*cabd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tvp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7c4=list(Tvd=-180,Tvb=-106,Tvp=320,Tprpd=1,Tprpb=-220)
modelcaTv7c4=mle(caTv7c4,start=guesscaTv7c4,method="BFGS")
summary(modelcaTv7c4) 
# confint(modelcaTv7c4) 

coef(modelcaTv7c4)
Tvd=coef(modelcaTv7c4)[1]
Tvb=coef(modelcaTv7c4)[2]
Tvp=coef(modelcaTv7c4)[3]
Tprpd=coef(modelcaTv7c4)[4]
Tprpb=coef(modelcaTv7c4)[5]
D=exp(((Tvd+Tvb+Tvp)*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvb*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvb*cabd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tvp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m7c4=c(D,B,P)
m7c4


##Model 7d4
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv7d4=function(Tvd,Tvb,Tvp,Tprbd,Tprpd,Tprpb){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvb+Tvp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvb*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvb*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tvp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7d4=list(Tvd=-330,Tvb=300,Tvp=-75,Tprbd=1,Tprpd=1,Tprpb=250)
modelcaTv7d4=mle(caTv7d4,start=guesscaTv7d4,method="BFGS")
summary(modelcaTv7d4) 
#confint(modelcaTv7d4) 

coef(modelcaTv7d4)
Tvd=coef(modelcaTv7d4)[1]
Tvb=coef(modelcaTv7d4)[2]
Tvp=coef(modelcaTv7d4)[3]
Tprbd=coef(modelcaTv7d4)[4]
Tprpd=coef(modelcaTv7d4)[5]
Tprpb=coef(modelcaTv7d4)[6]
D=exp(((Tvd+Tvb+Tvp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvb*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvb*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tvp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m7d4=c(D,B,P)
m7d4


##Model 7a5
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv7a5=function(Tvd,Tvbp,Tprpd){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvbp)*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvbp*cabd)/cabn)/(1+exp((Tvbp*cabd)/cabn))
capinf=exp((Tvbp*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvbp*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7a5=list(Tvd=1,Tvbp=1,Tprpd=1)
modelcaTv7a5=mle(caTv7a5,start=guesscaTv7a5,method="BFGS")
summary(modelcaTv7a5) 
#confint(modelcaTv7a5) 

coef(modelcaTv7a5)
Tvd=coef(modelcaTv7a5)[1]
Tvbp=coef(modelcaTv7a5)[2]
Tprpd=coef(modelcaTv7a5)[3]
D=exp(((Tvd+Tvbp)*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvbp*cabd)/cabn)/(1+exp((Tvbp*cabd)/cabn))
P=exp((Tvbp*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvbp*capd+Tprpd*cadi/cadn*cadd)/capn))
m7a5=c(D,B,P)
m7a5


##Model 7b5
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv7b5=function(Tvd,Tvbp,Tprbd,Tprpd){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvbp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvbp*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tvbp*cabd+Tprbd*cadi/cadn*cadd)/cabn))
capinf=exp((Tvbp*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvbp*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7b5=list(Tvd=-250,Tvbp=250,Tprbd=1,Tprpd=1)
modelcaTv7b5=mle(caTv7b5,start=guesscaTv7b5,method="BFGS")
summary(modelcaTv7b5) 
#confint(modelcaTv7b5) 

coef(modelcaTv7b5)
Tvd=coef(modelcaTv7b5)[1]
Tvbp=coef(modelcaTv7b5)[2]
Tprbd=coef(modelcaTv7b5)[3]
Tprpd=coef(modelcaTv7b5)[4]
D=exp(((Tvd+Tvbp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvbp*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tvbp*cabd+Tprbd*cadi/cadn*cadd)/cabn))
P=exp((Tvbp*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvbp*capd+Tprpd*cadi/cadn*cadd)/capn))
m7b5=c(D,B,P)
m7b5


##Model 7c5
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv7c5=function(Tvd,Tvbp,Tprpd,Tprpb){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvbp)*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvbp*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvbp*cabd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tvbp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvbp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7c5=list(Tvd=100,Tvbp=-100,Tprpd=0,Tprpb=-185)
modelcaTv7c5=mle(caTv7c5,start=guesscaTv7c5,method="BFGS")
summary(modelcaTv7c5) 
#confint(modelcaTv7c5) 

coef(modelcaTv7c5)
Tvd=coef(modelcaTv7c5)[1]
Tvbp=coef(modelcaTv7c5)[2]
Tprpd=coef(modelcaTv7c5)[3]
Tprpb=coef(modelcaTv7c5)[4]
D=exp(((Tvd+Tvbp)*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvbp*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvbp*cabd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tvbp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvbp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m7c5=c(D,B,P)
m7c5


##Model 7d5
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv7d5=function(Tvd,Tvbp,Tprbd,Tprpd,Tprpb){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvbp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvbp*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvbp*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tvbp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvbp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv7d5=list(Tvd=-300,Tvbp=300,Tprbd=-15,Tprpd=5,Tprpb=-245)
modelcaTv7d5=mle(caTv7d5,start=guesscaTv7d5,method="BFGS")
summary(modelcaTv7d5) 
#confint(modelcaTv7d5) 

coef(modelcaTv7d5)
Tvd=coef(modelcaTv7d5)[1]
Tvbp=coef(modelcaTv7d5)[2]
Tprbd=coef(modelcaTv7d5)[3]
Tprpd=coef(modelcaTv7d5)[4]
Tprpb=coef(modelcaTv7d5)[5]
D=exp(((Tvd+Tvbp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvbp*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvbp*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tvbp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvbp*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m7d5=c(D,B,P)
m7d5


####################################################
##Vector-borne + Aggressive encounters + Predation##
####################################################
##Model 8a1 
##vector trans, population density is important, male aggression, cross species trans from d to p by predation, Tv is uniform
caTv8a1=function(Tv,Tmd,Tmb,Tmp,Tprpd){    
cadinf=exp((Tv*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tv*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tv*cabd+Tmb*cabm*cabd)/cabn))
capinf=exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8a1=list(Tv=1,Tmd=1,Tmb=1,Tmp=1500,Tprpd=1)
modelcaTv8a1=mle(caTv8a1,start=guesscaTv8a1,method="BFGS")
summary(modelcaTv8a1) 
#confint(modelcaTv8a1)

coef(modelcaTv8a1)
Tv=coef(modelcaTv8a1)[1]
Tmd=coef(modelcaTv8a1)[2]
Tmb=coef(modelcaTv8a1)[3]
Tmp=coef(modelcaTv8a1)[4]
Tprpd=coef(modelcaTv8a1)[5]
D=exp((Tv*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tv*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tv*cabd+Tmb*cabm*cabd)/cabn))
P=exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
m8a1=c(D,B,P)
m8a1


##Model 8b1
##vector trans, population density is important, male aggression, cross species trans from d to b and p by predation, Tv is uniform
caTv8b1=function(Tv,Tmd,Tmb,Tmp,Tprbd,Tprpd){    
cadinf=exp((Tv*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tv*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tv*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn))
capinf=exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8b1=list(Tv=1,Tmd=1,Tmb=420,Tmp=1100,Tprbd=1,Tprpd=1)
modelcaTv8b1=mle(caTv8b1,start=guesscaTv8b1,method="BFGS")
summary(modelcaTv8b1) 
#confint(modelcaTv8b1) 

coef(modelcaTv8b1)
Tv=coef(modelcaTv8b1)[1]
Tmd=coef(modelcaTv8b1)[2]
Tmb=coef(modelcaTv8b1)[3]
Tmp=coef(modelcaTv8b1)[4]
Tprbd=coef(modelcaTv8b1)[5]
Tprpd=coef(modelcaTv8b1)[6]
D=exp((Tv*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tv*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tv*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn))
P=exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
m8b1=c(D,B,P)
m8b1


##Model 8c1
##vector trans, population density is important, male aggression, cross species trans from d and b to p by predation, Tv is uniform
caTv8c1=function(Tv,Tmd,Tmb,Tmp,Tprpd,Tprpb){    
cadinf=exp((Tv*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tv*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tv*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8c1=list(Tv=-10,Tmd=15,Tmb=-285,Tmp=2500,Tprpd=5,Tprpb=-375)
modelcaTv8c1=mle(caTv8c1,start=guesscaTv8c1,method="BFGS")
summary(modelcaTv8c1) 
#confint(modelcaTv8c1) 

coef(modelcaTv8c1)
Tv=coef(modelcaTv8c1)[1]
Tmd=coef(modelcaTv8c1)[2]
Tmb=coef(modelcaTv8c1)[3]
Tmp=coef(modelcaTv8c1)[4]
Tprpd=coef(modelcaTv8c1)[5]
Tprpb=coef(modelcaTv8c1)[6]
D=exp((Tv*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tv*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tv*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m8c1=c(D,B,P)
m8c1

##NB: larger positive values for derived parameters mean that they
##contribute more to transmission...see inverse logit transformations
##of parameter values below
#exp(Tv)/(1+exp(Tv))
#exp(Tmd)/(1+exp(Tmd))
#exp(Tmb)/(1+exp(Tmb))
#exp(Tmp)/(1+exp(Tmp))
#exp(Tprpd)/(1+exp(Tprpd))
#exp(Tprpb)/(1+exp(Tprpb))


##Model 8d1
##vector trans, population density is important, male aggression, cross species trans from d to b and p and b to p by predation, Tv is uniform
caTv8d1=function(Tv,Tmd,Tmb,Tmp,Tprbd,Tprpd,Tprpb){    
cadinf=exp((Tv*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tv*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tv*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8d1=list(Tv=-10,Tmd=10,Tmb=55,Tmp=2000,Tprbd=-10,Tprpd=5,Tprpb=-330)
modelcaTv8d1=mle(caTv8d1,start=guesscaTv8d1,method="BFGS")
summary(modelcaTv8d1) 
#confint(modelcaTv8d1) 

coef(modelcaTv8d1)
Tv=coef(modelcaTv8d1)[1]
Tmd=coef(modelcaTv8d1)[2]
Tmb=coef(modelcaTv8d1)[3]
Tmp=coef(modelcaTv8d1)[4]
Tprbd=coef(modelcaTv8d1)[5]
Tprpd=coef(modelcaTv8d1)[6]
Tprpb=coef(modelcaTv8d1)[7]
D=exp((Tv*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tv*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tv*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tv*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tv*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m8d1=c(D,B,P)
m8d1


##Model 8a2 
##vector trans, population density is important, male aggression, cross species trans from d to p by predation, b and p aquire vectors from d
caTv8a2=function(Tvd,Tvb,Tvp,Tmd,Tmb,Tmp,Tprpd){    
cadinf=exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd)/cabn)/(1+exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd)/cabn))
capinf=exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8a2=list(Tvd=-10,Tvb=-550,Tvp=1700,Tmd=10,Tmb=900,Tmp=-900,Tprpd=1)
modelcaTv8a2=mle(caTv8a2,start=guesscaTv8a2,method="BFGS")
summary(modelcaTv8a2) 
#confint(modelcaTv8a2)

coef(modelcaTv8a2)
Tvd=coef(modelcaTv8a2)[1]
Tvb=coef(modelcaTv8a2)[2]
Tvp=coef(modelcaTv8a2)[3]
Tmd=coef(modelcaTv8a2)[4]
Tmb=coef(modelcaTv8a2)[5]
Tmp=coef(modelcaTv8a2)[6]
Tprpd=coef(modelcaTv8a2)[7]
D=exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd)/cabn)/(1+exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd)/cabn))
P=exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
m8a2=c(D,B,P)
m8a2


##Model 8b2
##vector trans, population density is important, male aggression, cross species trans from d to b and p by predation, b and p aquire vectors from d
caTv8b2=function(Tvd,Tvb,Tvp,Tmd,Tmb,Tmp,Tprbd,Tprpd){    
cadinf=exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn))
capinf=exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8b2=list(Tvd=1,Tvb=-140,Tvp=900,Tmd=1,Tmb=585,Tmp=100,Tprbd=1,Tprpd=1)
modelcaTv8b2=mle(caTv8b2,start=guesscaTv8b2,method="BFGS")
summary(modelcaTv8b2) 
#confint(modelcaTv8b2) 

coef(modelcaTv8b2)
Tvd=coef(modelcaTv8b2)[1]
Tvb=coef(modelcaTv8b2)[2]
Tvp=coef(modelcaTv8b2)[3]
Tmd=coef(modelcaTv8b2)[4]
Tmb=coef(modelcaTv8b2)[5]
Tmp=coef(modelcaTv8b2)[6]
Tprbd=coef(modelcaTv8b2)[7]
Tprpd=coef(modelcaTv8b2)[8]
D=exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn))
P=exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
m8b2=c(D,B,P)
m8b2


##Model 8c2
##vector trans, population density is important, male aggression, cross species trans from d and b to p by predation, b and p aquire vectors from d
caTv8c2=function(Tvd,Tvb,Tvp,Tmd,Tmb,Tmp,Tprpd,Tprpb){    
cadinf=exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8c2=list(Tvd=-10,Tvb=-210,Tvp=-5,Tmd=15,Tmb=15,Tmp=2600,Tprpd=5,Tprpb=-380)
modelcaTv8c2=mle(caTv8c2,start=guesscaTv8c2,method="BFGS")
summary(modelcaTv8c2) 
#confint(modelcaTv8c2) 

coef(modelcaTv8c2)
Tvd=coef(modelcaTv8c2)[1]
Tvb=coef(modelcaTv8c2)[2]
Tvp=coef(modelcaTv8c2)[3]
Tmd=coef(modelcaTv8c2)[4]
Tmb=coef(modelcaTv8c2)[5]
Tmp=coef(modelcaTv8c2)[6]
Tprpd=coef(modelcaTv8c2)[7]
Tprpb=coef(modelcaTv8c2)[8]
D=exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m8c2=c(D,B,P)
m8c2


##Model 8d2
##vector trans, population density is important, male aggression, cross species trans from d to b and p and b to p by predation, b and p aquire vectors from d
caTv8d2=function(Tvd,Tvb,Tvp,Tmd,Tmb,Tmp,Tprbd,Tprpd,Tprpb){    
cadinf=exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8d2=list(Tvd=-10,Tvb=-50,Tvp=185,Tmd=10,Tmb=95,Tmp=150,Tprbd=-5,Tprpd=10,Tprpb=-320)
modelcaTv8d2=mle(caTv8d2,start=guesscaTv8d2,method="BFGS")
summary(modelcaTv8d2) 
#confint(modelcaTv8d2) 

coef(modelcaTv8d2)
Tvd=coef(modelcaTv8d2)[1]
Tvb=coef(modelcaTv8d2)[2]
Tvp=coef(modelcaTv8d2)[3]
Tmd=coef(modelcaTv8d2)[4]
Tmb=coef(modelcaTv8d2)[5]
Tmp=coef(modelcaTv8d2)[6]
Tprbd=coef(modelcaTv8d2)[7]
Tprpd=coef(modelcaTv8d2)[8]
Tprpb=coef(modelcaTv8d2)[9]
D=exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvb)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m8d2=c(D,B,P)
m8d2


##Model 8a3 
##vector trans, population density is important, male aggression, cross species trans from d to p by predation, b and p aquire vectors from d
caTv8a3=function(Tvd,Tvbp,Tmd,Tmb,Tmp,Tprpd){    
cadinf=exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd)/cabn))
capinf=exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8a3=list(Tvd=1,Tvbp=1,Tmd=1,Tmb=1000,Tmp=1500,Tprpd=1)
modelcaTv8a3=mle(caTv8a3,start=guesscaTv8a3,method="BFGS")
summary(modelcaTv8a3) 
#confint(modelcaTv8a3)

coef(modelcaTv8a3)
Tvd=coef(modelcaTv8a3)[1]
Tvbp=coef(modelcaTv8a3)[2]
Tmd=coef(modelcaTv8a3)[3]
Tmb=coef(modelcaTv8a3)[4]
Tmp=coef(modelcaTv8a3)[5]
Tprpd=coef(modelcaTv8a3)[6]
D=exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd)/cabn))
P=exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
m8a3=c(D,B,P)
m8a3


##Model 8b3
##vector trans, population density is important, male aggression, cross species trans from d to b and p by predation, b and p aquire vectors from d
caTv8b3=function(Tvd,Tvbp,Tmd,Tmb,Tmp,Tprbd,Tprpd){    
cadinf=exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn))
capinf=exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8b3=list(Tvd=1,Tvbp=100,Tmd=1,Tmb=10,Tmp=1,Tprbd=1,Tprpd=1)
modelcaTv8b3=mle(caTv8b3,start=guesscaTv8b3,method="BFGS")
summary(modelcaTv8b3) 
#confint(modelcaTv8b3) 

coef(modelcaTv8b3)
Tvd=coef(modelcaTv8b3)[1]
Tvbp=coef(modelcaTv8b3)[2]
Tmd=coef(modelcaTv8b3)[3]
Tmb=coef(modelcaTv8b3)[4]
Tmp=coef(modelcaTv8b3)[5]
Tprbd=coef(modelcaTv8b3)[6]
Tprpd=coef(modelcaTv8b3)[7]
D=exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn))
P=exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
m8b3=c(D,B,P)
m8b3


##Model 8c3
##vector trans, population density is important, male aggression, cross species trans from d and b to p by predation, b and p aquire vectors from d
caTv8c3=function(Tvd,Tvbp,Tmd,Tmb,Tmp,Tprpd,Tprpb){    
cadinf=exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8c3=list(Tvd=-10,Tvbp=-215,Tmd=15,Tmb=35,Tmp=3000,Tprpd=5,Tprpb=-370)
modelcaTv8c3=mle(caTv8c3,start=guesscaTv8c3,method="BFGS")
summary(modelcaTv8c3) 
#confint(modelcaTv8c3) 

coef(modelcaTv8c3)
Tvd=coef(modelcaTv8c3)[1]
Tvbp=coef(modelcaTv8c3)[2]
Tmd=coef(modelcaTv8c3)[3]
Tmb=coef(modelcaTv8c3)[4]
Tmp=coef(modelcaTv8c3)[5]
Tprpd=coef(modelcaTv8c3)[6]
Tprpb=coef(modelcaTv8c3)[7]
D=exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m8c3=c(D,B,P)
m8c3


##Model 8d3
##vector trans, population density is important, male aggression, cross species trans from d to b and p and b to p by predation, b and p aquire vectors from d
caTv8d3=function(Tvd,Tvbp,Tmd,Tmb,Tmp,Tprbd,Tprpd,Tprpb){    
cadinf=exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8d3=list(Tvd=-10,Tvbp=-50,Tmd=10,Tmb=20,Tmp=2500,Tprbd=-7,Tprpd=7,Tprpb=-360)
modelcaTv8d3=mle(caTv8d3,start=guesscaTv8d3,method="BFGS")
summary(modelcaTv8d3) 
#confint(modelcaTv8d3) 

coef(modelcaTv8d3)
Tvd=coef(modelcaTv8d3)[1]
Tvbp=coef(modelcaTv8d3)[2]
Tmd=coef(modelcaTv8d3)[3]
Tmb=coef(modelcaTv8d3)[4]
Tmp=coef(modelcaTv8d3)[5]
Tprbd=coef(modelcaTv8d3)[6]
Tprpd=coef(modelcaTv8d3)[7]
Tprpb=coef(modelcaTv8d3)[8]
D=exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp((Tvd*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp(((Tvd+Tvbp)*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp(((Tvd+Tvbp)*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m8d3=c(D,B,P)
m8d3


##Model 8a4
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv8a4=function(Tvd,Tvb,Tvp,Tmd,Tmb,Tmp,Tprpd){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvb*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tvb*cabd+Tmb*cabm*cabd)/cabn))
capinf=exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8a4=list(Tvd=-800,Tvb=-650,Tvp=1200,Tmd=10,Tmb=1010,Tmp=-30,Tprpd=1)
modelcaTv8a4=mle(caTv8a4,start=guesscaTv8a4,method="BFGS")
summary(modelcaTv8a4) 
#confint(modelcaTv8a4) 

coef(modelcaTv8a4)
Tvd=coef(modelcaTv8a4)[1]
Tvb=coef(modelcaTv8a4)[2]
Tvp=coef(modelcaTv8a4)[3]
Tmd=coef(modelcaTv8a4)[4]
Tmb=coef(modelcaTv8a4)[5]
Tmp=coef(modelcaTv8a4)[6]
Tprpd=coef(modelcaTv8a4)[7]
D=exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvb*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tvb*cabd+Tmb*cabm*cabd)/cabn))
P=exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
m8a4=c(D,B,P)
m8a4


##Model 8b4
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv8b4=function(Tvd,Tvb,Tvp,Tmd,Tmb,Tmp,Tprbd,Tprpd){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tvb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn))
capinf=exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8b4=list(Tvd=-1000,Tvb=-290,Tvp=770,Tmd=10,Tmb=760,Tmp=-5,Tprbd=-10,Tprpd=-1)
modelcaTv8b4=mle(caTv8b4,start=guesscaTv8b4,method="BFGS")
summary(modelcaTv8b4) 
#confint(modelcaTv8b4) 

coef(modelcaTv8b4)
Tvd=coef(modelcaTv8b4)[1]
Tvb=coef(modelcaTv8b4)[2]
Tvp=coef(modelcaTv8b4)[3]
Tmd=coef(modelcaTv8b4)[4]
Tmb=coef(modelcaTv8b4)[5]
Tmp=coef(modelcaTv8b4)[6]
Tprbd=coef(modelcaTv8b4)[7]
Tprpd=coef(modelcaTv8b4)[8]
D=exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tvb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn))
P=exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
m8b4=c(D,B,P)
m8b4


##Model 8c4
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv8c4=function(Tvd,Tvb,Tvp,Tmd,Tmb,Tmp,Tprpd,Tprpb){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvb*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvb*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8c4=list(Tvd=150,Tvb=-205,Tvp=50,Tmd=15,Tmb=-1,Tmp=2100,Tprpd=5,Tprpb=-370)
modelcaTv8c4=mle(caTv8c4,start=guesscaTv8c4,method="BFGS")
summary(modelcaTv8c4) 
#confint(modelcaTv8c4) 

coef(modelcaTv8c4)
Tvd=coef(modelcaTv8c4)[1]
Tvb=coef(modelcaTv8c4)[2]
Tvp=coef(modelcaTv8c4)[3]
Tmd=coef(modelcaTv8c4)[4]
Tmb=coef(modelcaTv8c4)[5]
Tmp=coef(modelcaTv8c4)[6]
Tprpd=coef(modelcaTv8c4)[7]
Tprpb=coef(modelcaTv8c4)[8]
D=exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvb*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvb*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m8c4=c(D,B,P)
m8c4


##Model 8d4
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv8d4=function(Tvd,Tvb,Tvp,Tmd,Tmb,Tmp,Tprbd,Tprpd,Tprpb){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8d4=list(Tvd=-105,Tvb=-120,Tvp=210,Tmd=10,Tmb=170,Tmp=280,Tprbd=-5,Tprpd=5,Tprpb=-325)
modelcaTv8d4=mle(caTv8d4,start=guesscaTv8d4,method="BFGS")
summary(modelcaTv8d4) 
#confint(modelcaTv8d4) 

coef(modelcaTv8d4)
Tvd=coef(modelcaTv8d4)[1]
Tvb=coef(modelcaTv8d4)[2]
Tvp=coef(modelcaTv8d4)[3]
Tmd=coef(modelcaTv8d4)[4]
Tmb=coef(modelcaTv8d4)[5]
Tmp=coef(modelcaTv8d4)[6]
Tprbd=coef(modelcaTv8d4)[7]
Tprpd=coef(modelcaTv8d4)[8]
Tprpb=coef(modelcaTv8d4)[9]
D=exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvb+Tvp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvb*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m8d4=c(D,B,P)
m8d4


##Model 8a5
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv8a5=function(Tvd,Tvbp,Tmd,Tmb,Tmp,Tprpd){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvbp*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tvbp*cabd+Tmb*cabm*cabd)/cabn))
capinf=exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8a5=list(Tvd=80,Tvbp=-800,Tmd=3,Tmb=570,Tmp=2000,Tprpd=1)
modelcaTv8a5=mle(caTv8a5,start=guesscaTv8a5,method="BFGS")
summary(modelcaTv8a5) 
#confint(modelcaTv8a5) 

coef(modelcaTv8a5)
Tvd=coef(modelcaTv8a5)[1]
Tvbp=coef(modelcaTv8a5)[2]
Tmd=coef(modelcaTv8a5)[3]
Tmb=coef(modelcaTv8a5)[4]
Tmp=coef(modelcaTv8a5)[5]
Tprpd=coef(modelcaTv8a5)[6]
D=exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvbp*cabd+Tmb*cabm*cabd)/cabn)/(1+exp((Tvbp*cabd+Tmb*cabm*cabd)/cabn))
P=exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
m8a5=c(D,B,P)
m8a5


##Model 8b5
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv8b5=function(Tvd,Tvbp,Tmd,Tmb,Tmp,Tprbd,Tprpd){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvbp*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tvbp*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn))
capinf=exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8b5=list(Tvd=1,Tvbp=1,Tmd=1,Tmb=1,Tmp=800,Tprbd=1,Tprpd=1)
modelcaTv8b5=mle(caTv8b5,start=guesscaTv8b5,method="BFGS")
summary(modelcaTv8b5) 
#confint(modelcaTv8b5) 

coef(modelcaTv8b5)
Tvd=coef(modelcaTv8b5)[1]
Tvbp=coef(modelcaTv8b5)[2]
Tmd=coef(modelcaTv8b5)[3]
Tmb=coef(modelcaTv8b5)[4]
Tmp=coef(modelcaTv8b5)[5]
Tprbd=coef(modelcaTv8b5)[6]
Tprpd=coef(modelcaTv8b5)[7]
D=exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvbp*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn)/(1+exp((Tvbp*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd)/cabn))
P=exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn)/(1+exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd)/capn))
m8b5=c(D,B,P)
m8b5


##Model 8c5
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv8c5=function(Tvd,Tvbp,Tmd,Tmb,Tmp,Tprpd,Tprpb){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvbp*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvbp*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8c5=list(Tvd=195,Tvbp=-205,Tmd=14,Tmb=2,Tmp=2200,Tprpd=5,Tprpb=-370)
modelcaTv8c5=mle(caTv8c5,start=guesscaTv8c5,method="BFGS")
summary(modelcaTv8c5) 
#confint(modelcaTv8c5) 

coef(modelcaTv8c5)
Tvd=coef(modelcaTv8c5)[1]
Tvbp=coef(modelcaTv8c5)[2]
Tmd=coef(modelcaTv8c5)[3]
Tmb=coef(modelcaTv8c5)[4]
Tmp=coef(modelcaTv8c5)[5]
Tprpd=coef(modelcaTv8c5)[6]
Tprpb=coef(modelcaTv8c5)[7]
D=exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvbp*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvbp*cabd+Tmb*cabm*cabd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m8c5=c(D,B,P)
m8c5


##Model 8d5
##vector trans, population density is important, male aggression, b and p share vectors, d aquires vectors from b and p
caTv8d5=function(Tvd,Tvbp,Tmd,Tmb,Tmp,Tprbd,Tprpd,Tprpb){    #not able to split Tvbp on this model
cadinf=exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
cabinf=exp((Tvbp*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvbp*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
capinf=exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv8d5=list(Tvd=30,Tvbp=-40,Tmd=10,Tmb=25,Tmp=2165,Tprbd=-5,Tprpd=5,Tprpb=-350)
modelcaTv8d5=mle(caTv8d5,start=guesscaTv8d5,method="BFGS")
summary(modelcaTv8d5) 
#confint(modelcaTv8d5) 

coef(modelcaTv8d5)
Tvd=coef(modelcaTv8d5)[1]
Tvbp=coef(modelcaTv8d5)[2]
Tmd=coef(modelcaTv8d5)[3]
Tmb=coef(modelcaTv8d5)[4]
Tmp=coef(modelcaTv8d5)[5]
Tprbd=coef(modelcaTv8d5)[6]
Tprpd=coef(modelcaTv8d5)[7]
Tprpb=coef(modelcaTv8d5)[8]
D=exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn)/(1+exp(((Tvd+Tvbp)*cadd+Tmd*cadm*cadd-Tprbd*cadi/cadn*cadd-Tprpd*cadi/cadn*cadd)/cadn))
B=exp((Tvbp*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn)/(1+exp((Tvbp*cabd+Tmb*cabm*cabd+Tprbd*cadi/cadn*cadd-Tprpb*cabi/cabn*cabd)/cabn))
P=exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn)/(1+exp((Tvbp*capd+Tmp*capm*capd+Tprpd*cadi/cadn*cadd+Tprpb*cabi/cabn*cabd)/capn))
m8d5=c(D,B,P)
m8d5


#################
##Environmental##
#################
##Model 9a
##Environmental trans is density dependent
caTv9a=function(Ted,Teb,Tep){    
cadinf=exp((Ted*(cadd+cabd+capd))/cadn)/(1+exp((Ted*(cadd+cabd+capd))/cadn))
cabinf=exp((Teb*(cadd+cabd+capd))/cabn)/(1+exp((Teb*(cadd+cabd+capd))/cabn))
capinf=exp((Tep*(cadd+cabd+capd))/capn)/(1+exp((Tep*(cadd+cabd+capd))/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv9a=list(Ted=1,Teb=1,Tep=1)
modelcaTv9a=mle(caTv9a,start=guesscaTv9a,method="BFGS")
summary(modelcaTv9a) 
#confint(modelcaTv9a) 

coef(modelcaTv9a)
Ted=coef(modelcaTv9a)[1]
Teb=coef(modelcaTv9a)[2]
Tep=coef(modelcaTv9a)[3]
D=exp((Ted*(cadd+cabd+capd))/cadn)/(1+exp((Ted*(cadd+cabd+capd))/cadn))
B=exp((Teb*(cadd+cabd+capd))/cabn)/(1+exp((Teb*(cadd+cabd+capd))/cabn))
P=exp((Tep*(cadd+cabd+capd))/capn)/(1+exp((Tep*(cadd+cabd+capd))/capn))
m9a=c(D,B,P)
m9a


#########################################
##Environmental + Aggressive encounters##
#########################################
##Model 10a
##Environmental trans is density dependent
caTv10a=function(Ted,Teb,Tep,Tmd,Tmb,Tmp){    
cadinf=exp((Ted*(cadd+cabd+capd)+Tmd*cadm*cadd)/cadn)/(1+exp((Ted*(cadd+cabd+capd)+Tmd*cadm*cadd)/cadn))
cabinf=exp((Teb*(cadd+cabd+capd)+Tmb*cabm*cabd)/cabn)/(1+exp((Teb*(cadd+cabd+capd)+Tmb*cabm*cabd)/cabn))
capinf=exp((Tep*(cadd+cabd+capd)+Tmp*capm*capd)/capn)/(1+exp((Tep*(cadd+cabd+capd)+Tmp*capm*capd)/capn))
L1=sum(dbinom(cadi,size=cadn,prob=cadinf,log=TRUE))  
L2=sum(dbinom(cabi,size=cabn,prob=cabinf,log=TRUE))
L3=sum(dbinom(capi,size=capn,prob=capinf,log=TRUE))
L=L1+L2+L3                       #add the likelihoods
return(-L)                       #return negative likelihood
}

guesscaTv10a=list(Ted=1,Teb=1,Tep=1,Tmd=1,Tmb=1,Tmp=1)
modelcaTv10a=mle(caTv10a,start=guesscaTv10a,method="BFGS")
summary(modelcaTv10a) 
#confint(modelcaTv10a) 

coef(modelcaTv10a)
Ted=coef(modelcaTv10a)[1]
Teb=coef(modelcaTv10a)[2]
Tep=coef(modelcaTv10a)[3]
Tmd=coef(modelcaTv10a)[4]
Tmb=coef(modelcaTv10a)[5]
Tmp=coef(modelcaTv10a)[6]
D=exp((Ted*(cadd+cabd+capd)+Tmd*cadm*cadd)/cadn)/(1+exp((Ted*(cadd+cabd+capd)+Tmd*cadm*cadd)/cadn))
B=exp((Teb*(cadd+cabd+capd)+Tmb*cabm*cabd)/cabn)/(1+exp((Teb*(cadd+cabd+capd)+Tmb*cabm*cabd)/cabn))
P=exp((Tep*(cadd+cabd+capd)+Tmp*capm*capd)/capn)/(1+exp((Tep*(cadd+cabd+capd)+Tmp*capm*capd)/capn))
m10a=c(D,B,P) 
m10a




########################################
##Colating the results from the models##
########################################

modnames=c("1a","2a","3a","4a","4b","4c","4d","5a","5b","5c","5d","6a","6b","6c","6d","6e","7a","7b","7c","7d","7e","8a1","8b1",
           "8c1","8d1","8a2","8b2","8c2","8d2","8a3","8b3","8c3","8d3","8a4","8b4","8c4","8d4","8a5","8b5","8c5","8d5","9a1",
           "9b1","9c1","9d1","9a2","9b2","9c2","9d2","9a3","9b3","9c3","9d3","9a4","9b4","9c4","9d4","9a5","9b5","9c5","9d5","10a","11a")

mods=list(modelcaTv1a, modelcaTv2a1, modelcaTv2a, modelcaTv3a, modelcaTv3b, modelcaTv3c, modelcaTv3d, modelcaTv4a, modelcaTv4b, 
          modelcaTv4c, modelcaTv4d, modelcaTv5a, modelcaTv5b, modelcaTv5c, modelcaTv5d, modelcaTv5e, modelcaTv6a, modelcaTv6b,modelcaTv6a, 
          modelcaTv6d, modelcaTv6e, modelcaTv7a1, modelcaTv7b1, modelcaTv7c1, modelcaTv7d1, modelcaTv7a2,modelcaTv7b2, modelcaTv7c2, modelcaTv7d2, 
          modelcaTv7a3, modelcaTv7b3, modelcaTv7c3, modelcaTv7d3, modelcaTv7a4, modelcaTv7b4, modelcaTv7c4, modelcaTv7d4, modelcaTv7a5, modelcaTv7b5, 
          modelcaTv7c5, modelcaTv7d5, modelcaTv8a1, modelcaTv8b1, modelcaTv8c1, modelcaTv8d1, modelcaTv8a2, modelcaTv8b2, modelcaTv8c2, modelcaTv8d2, 
          modelcaTv8a3, modelcaTv8b3, modelcaTv8c3, modelcaTv8d3, modelcaTv8a4, modelcaTv8b4, modelcaTv8c4, modelcaTv8d4, modelcaTv8a5, modelcaTv8b5, 
          modelcaTv8c5, modelcaTv8d5, modelcaTv9a, modelcaTv10a)

K=c(NA)
neg2logL=c(NA)
aic=c(NA)
for(i in 1:length(mods)){
  K[i]=length(coef(mods[[i]]))#K
  neg2logL[i]=-2*logLik(mods[[i]])#-2Log(L)
  aic[i]=AIC(mods[[i]])#AIC
}
K
neg2logL
aic
deltaAIC=aic-min(aic)
weight=exp(-1/2*deltaAIC)/sum(exp(-1/2*deltaAIC))
weightpct=weight*100
table1=cbind(1:63,K,neg2logL,aic,deltaAIC,weight,weightpct)
colnames(table1)=c("Model","K","-2LOG(L)","AIC","deltaAIC","Weight","Weight (%)")
rownames(table1)=modnames
table1=table1[order(table1[,4]),]
table1=data.frame(table1)
table1
write.csv(table1,file="AIC table.csv")


allmodelresults=cbind(m1a,m2a1,m2a,m3a,m3b,m3c,m3d,m4a,m4b,m4c,m4d,m5a,m5b,m5c,m5d,m5e,m6a,m6b,m6c,m6d,m6e
                 ,m7a1,m7b1,m7c1,m7d1,m7a2,m7b2,m7c2,m7d2,m7a3,m7b3,m7c3,m7d3,m7a4,m7b4,m7c4,m7d4,m7a5,m7b5,m7c5,m7d5
                 ,m8a1,m8b1,m8c1,m8d1,m8a2,m8b2,m8c2,m8d2,m8a3,m8b3,m8c3,m8d3,m8a4,m8b4,m8c4,m8d4,m8a5,m8b5,m8c5,m8d5
                 ,m9a,m10a)

allmodelres=list(m1a,m2a1,m2a,m3a,m3b,m3c,m3d,m4a,m4b,m4c,m4d,m5a,m5b,m5c,m5d,m5e,m6a,m6b,m6c,m6d,m6e
                      ,m7a1,m7b1,m7c1,m7d1,m7a2,m7b2,m7c2,m7d2,m7a3,m7b3,m7c3,m7d3,m7a4,m7b4,m7c4,m7d4,m7a5,m7b5,m7c5,m7d5
                      ,m8a1,m8b1,m8c1,m8d1,m8a2,m8b2,m8c2,m8d2,m8a3,m8b3,m8c3,m8d3,m8a4,m8b4,m8c4,m8d4,m8a5,m8b5,m8c5,m8d5
                      ,m9a,m10a)

allmodelaverage=list(NA)
for(i in 1:length(mods)){
  allmodelaverage[[i]]=allmodelres[[i]]*weight[i]
}
allmodelaverage
allmodavg=rowSums(mapply(c,allmodelaverage))
table2=cbind(allmodelresults,allmodavg)
colnames(table2)=c(append(modnames,"modavg"))
rownames(table2)=c("NLA domestic","SLA domestic","WS domestic","FR domestic",
                   "NLA bobcat","SLA bobcat","WS bobcat","FR bobcat",
                   "NLA puma","SLA puma","WS puma","FR puma")
table2
write.csv(table2,file="Model prevalence table.csv")



####################################################
##Evaluate how well the models fit prevalence data##
####################################################

dat=c(AllSites[,5],AllSites[,11],AllSites[,17])
dat

jpeg(file="Model Fit Figure.jpg",width=2000,height=2000,res=300)  #specifies the name of the image file, file type, dimensions and resolution
par(mar=c(5,5,2,1))#b,l,t,r

plot(dat,allmodavg,type="n",cex.lab=1.25,cex=2,xlim=c(0,0.7),ylim=c(0,0.7),ylab="Predicted infection prevalence (all models averaged)",xlab="Observed infection prevalence")
points(dat[1:4],allmodavg[1:4],pch=21,cex=2,col="black",bg="black")
points(dat[5:8],allmodavg[5:8],pch=21,cex=2,col="black",bg="grey")
points(dat[9:12],allmodavg[9:12],pch=21,cex=2,col="black",bg="white")
fit2=lm(allmodavg~dat)
abline(fit2)
summary(fit2)
cor.test(dat,allmodavg,method=c("spearman"))
legend("topleft",expression(paste(rho," = 0.874, ",italic('P')," <0.001")),bty='n')

dev.off() #prints the image file



#############################
##Variable importance plots## 
#############################

importance=read.csv("Variable importance.csv")
importance

jpeg(file="Proportional Variable Importance Figure.jpg",width=2000,height=2000,res=300)  #specifies the name of the image file, file type, dimensions and resolution

par(mar=c(5,5,2,1),oma=c(0,2,0,0))#b,l,t,r
plot(importance$deltaAIC,importance$S,type="l",lty=1,ylim=c(0,0.8),xlim=c(0,4),cex.lab=1.25
,xlab=expression(paste(Delta, "AIC (number of models)")),ylab=''
,xaxt='n') #removes x-axis values
mtext("Proporation of models occupied \n by each transmission route",side=2,line=3,cex=1.25)
axis(1,at=0:4,lab=c("0 (1)","1 (3)","2 (6)","3 (10)","4 (16)"))#,"5 (21)"))
points(importance$deltaAIC,importance$A,type="l",lty=2)
points(importance$deltaAIC,importance$Pa,type="l",lty=3)
points(importance$deltaAIC,importance$Pb,type="l",lty=4)
points(importance$deltaAIC,importance$Pc,type="l",lty=5)
points(importance$deltaAIC,importance$Pd,type="l",lty=6)
points(importance$deltaAIC,importance$Va,type="l",lty=1,lwd=2)
points(importance$deltaAIC,importance$Vb,type="l",lty=2,lwd=2)
points(importance$deltaAIC,importance$Vc,type="l",lty=3,lwd=2)
points(importance$deltaAIC,importance$Vd,type="l",lty=4,lwd=2)
points(importance$deltaAIC,importance$Ve,type="l",lty=5,lwd=2)
points(importance$deltaAIC,importance$E,type="l",lty=6,lwd=2)
legend("topleft",bty="n",title="Transmission route (number of models occuring in)"
,c("1a (10)","1b (32)*","2a (12)*","2b (12)","2c (12)*","2d (12)","3a (10)*","3b (10)*","3c (10)*","3d (10)*","3e (10)*","3f (2)")
,lty=1:12,lwd=c(1,1,1,1,1,1,2,2,2,2,2,2))

dev.off() #prints the image file



##Parameter importance using regression trees##
##i.e., which parameter predict the best models deltaAIC best##

rm(list=ls())
require(rpart)
require(randomForest)
require(rpart.plot)
?rpart.plot
?prp
ParamImp=read.csv("Regression tree.csv")
ParamImp

## randomForest Regression:
set.seed(131)
deltaAIC.rf <- randomForest(ParamImp$deltaAIC ~ ParamImp$Transmission_mode_1a+
                              ParamImp$Transmission_mode_1b+
                              ParamImp$Transmission_mode_2a+
                              ParamImp$Transmission_mode_2b+
                              ParamImp$Transmission_mode_2c+
                              ParamImp$Transmission_mode_2d+
                              ParamImp$Transmission_mode_3a+
                              ParamImp$Transmission_mode_3b+
                              ParamImp$Transmission_mode_3c+
                              ParamImp$Transmission_mode_3d+
                              ParamImp$Transmission_mode_3e+
                              ParamImp$Transmission_mode_3f, 
				mtry=8,ntree=100000,importance=TRUE, na.action=na.omit)
print(deltaAIC.rf)
## Show "importance" of variables: higher value mean more important:
round(importance(deltaAIC.rf), 2)
write.csv(round(importance(deltaAIC.rf), 2)[,2],file="Variable importance values.csv")
## Plot variable importance
varImpPlot(deltaAIC.rf,sort=FALSE,type=1)
Model=c("1a","1b","2a","2b","2c","2d","3a","3b","3c","3d","3e","3f")
IncNodePurity=round(importance(deltaAIC.rf), 2)[,2]

jpeg(file="Variable Importance Figure.jpg",width=2000,height=2000,res=300)  #specifies the name of the image file, file type, dimensions and resolution
par(mar=c(5,5,1,1),oma=c(0,0,0,0))#b,l,t,r

barplot(rev(IncNodePurity),horiz=TRUE,names.arg=rev(Model),xlim=c(0,max(IncNodePurity)),cex.lab=1.1,las=1,xlab="Variable importance (mean decrease in Gini index)",ylab="Transmission mode")

dev.off() #prints the image file

##Notes##
#Mean Decrease Accuracy (%IncMSE): It is constructed by
#permuting the values of each variable of the test set,
#recording the prediction and comparing it with the unpermuted
#test set prediction of the variable (normalised by
#the standard error). For classification, it is the increase in
#the percentage of times a test set tuple is misclassified
#when the variable is permuted. For regression, it is the
#average increase in squared residuals of the test set when
#the variable is permuted. A higher %IncMSE value represents
#a higher variable importance.

#Mean Decrease Gini (IncNodePurity): Measures the
#quality (NodePurity) of a split for every variable (node) of
#a tree by means of the Gini Index. Every time a split of a
#node is made on a variable the gini impurity criterion for
#the two descendent nodes is less than the parent node.
#Adding up the gini decreases for each individual variable
#over all trees in the forest gives a fast variable importance
#that is often very consistent with the permutation importance
#measure. A higher IncNodePurity value represents a
#higher variable importance, i.e. nodes are much 'purer'.





########
##Junk##
########

#See varible importance values at the bottom of "Variable importance.xls"
Transmode=c("S","A","Ptotal","Pa","Pb","Pc","Pd","Vtotal","Va","Vb","Vc","Vd","Ve","E")
VarImp=c(0.069702883,0.998307035,0.74979853,0.249412872,0.044928436,0.333002326,0.122454896,0.913730448,0.10880619,0.146435368,0.310533726,0.130332903,0.217622261,0.016564724)
WgtdVarImp=c(0.01106395,0.049519198,0.024794925,0.032991121,0.005942915,0.044047927,0.016197738,0.029007316,0.017270824,0.023243709,0.049291068,0.020687762,0.034543216,0.013146606)
VI=cbind(VarImp,WgtdVarImp)
rownames(VI)=Transmode
VI

barplot(VI[c(1:2,4:7,9:14),1],beside=T,horiz=T,legend.text=c("1a","1b","2a","2b","2c","2d","3a","3b","3c","3d","3e","3f"),xlab="Variable importance")



########################
##Parameter importance##
########################


## rpart tree
deltaAIC.rp <- rpart(deltaAIC ~ X1a+X1b+X2a+X2b+X2c+X2d+X3a+X3b+X3c+X3d+X3e+X3f) 
deltaAIC.rp
summary(deltaAIC.rp)
prp(deltaAIC.rp,type=1,extra=1,branch=0.1,uniform=FALSE,left=FALSE,round=1,varlen=0,faclen=0,xflip=TRUE,cex=1.1,tweak=1,nn.cex=1.1,box.col="lightgrey")

##Notes##
#The splitting process continues until further subdivision no longer
#reduces the Gini index. Such classification tree is said to be fully
#grown, and the final regaions are called terminal nodes
