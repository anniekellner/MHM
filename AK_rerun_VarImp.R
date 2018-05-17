rm(list=ls())
require(rpart)
require(randomForest)
require(rpart.plot)
ParamImp=read.csv("Regression tree_AK.csv")
ParamImp

## randomForest Regression:
set.seed(131)
deltaAICc.rf <- randomForest(ParamImp$deltaAICc ~ ParamImp$Transmission_mode_1a+
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
print(deltaAICc.rf)

## Show "importance" of variables: higher value mean more important:
round(importance(deltaAICc.rf), 2)
write.csv(round(importance(deltaAICc.rf), 2)[,2],file="Variable importance values.csv")

## Plot variable importance
varImpPlot(deltaAICc.rf,sort=FALSE,type=1)
Model=c("1a","1b","2a","2b","2c","2d","3a","3b","3c","3d","3e","3f")
IncNodePurity=round(importance(deltaAICc.rf), 2)[,2]

jpeg(file="Variable Importance Figure.jpg",width=2000,height=2000,res=300)  #specifies the name of the image file, file type, dimensions and resolution
par(mar=c(5,5,1,1),oma=c(0,0,0,0))#b,l,t,r

barplot(rev(IncNodePurity),horiz=TRUE,names.arg=rev(Model),xlim=c(0,max(IncNodePurity)),cex.lab=1.1,las=1,xlab="Variable importance (mean decrease in Gini index)",ylab="Transmission mode")

#dev.off() #prints the image file