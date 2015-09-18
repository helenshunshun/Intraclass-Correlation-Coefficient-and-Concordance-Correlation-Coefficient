#### This is for making a big R function that can create table, when imput ralated ICC and CCC parameters, it can 
# make simulaitons and output result.
require(psych)
require(MASS)
require(epiR)

#####################################################
#### this function is for making a summary of ICC result based on once simulaiton
iccsum<- function(avg.icc,avg.rankicc,cov.icc,cov.rankicc,low.icc,up.icc,low.rankicc,up.rankicc,trueicc){
  
  mean.icc<-mean(avg.icc) #### this is point estimate
  mean.rankicc<-mean(avg.rankicc)
  ## now to get coverage prob
  prob.icc<-length(cov.icc[cov.icc!=0])/100
  prob.rankicc<-length(cov.rankicc[cov.rankicc!=0])/100
  icc.matrix<-cbind(low.icc,up.icc,low.rankicc,up.rankicc)
  CI.icc<-apply(icc.matrix,2,mean)
  ## get bias
  bias.icc<-mean.icc-trueicc
  bias.rankicc<-mean.rankicc-trueicc
  sum.icc<-c(trueicc, round(c(mean.icc, bias.icc, mean.rankicc, bias.rankicc),4), c(prob.icc, prob.rankicc)*100, round(CI.icc,4))
  names(sum.icc)<-c("trueicc", "esti.icc","bias.icc", "esti.rankicc","bias.rankicc", "cov.icc","cov.rankicc",
                    "low.icc","up.icc","low.rankicc","up.rankicc")
  return(sum.icc)
}

#### this function is for making a summary of CCC result based on once simulaiton
### similar to fucntion cccsum
########################
# make a fucntion to get different types of outlier

######## 3 types of outliers, need 3 ways to choose
genoutlier<-function(type,outlier,percent,ybi11,ybi12,ybi21,ybi22){
  if(outlier=="reading"){
    out<-c(ybi11,ybi12,ybi21,ybi22)
    index<-sample(1:400,40)
    for (n in 1:400){
      for(j in 1:40){
        if(n==index[j]){
          out[n]=out[index[j]]*percent
        }
      } 
    }
    y<-as.data.frame(cbind(out[1:100],out[101:200],out[201:300],out[301:400])) ##for unify,here don't use y1, only use y
  }
  
  if(outlier=="subject"){
    y<-as.data.frame(cbind(ybi11,ybi12,ybi21,ybi22))
    index<-sample(1:100,10)###choose randomly 1% of data
    for(m in 1:100){
      for(j in 1:10){
        if(m==index[j]){
          y[m,]<-y[index[j],]*percent   ###inflate by 10 and get new y
        }
      }
    }
  }
  
  
  if(outlier=="observer"){
    y<-as.data.frame(cbind(ybi11,ybi12,ybi21,ybi22))
    index<-sample(1:100,10)
    col<-sample(1:4,10,replace=T)
    for(m in 1:100){
      for(c in 1:4){
        for(j in 1:10){
          if(m==index[j] & c==col[j]){
            y[m,c]<-y[index[j],col[j]]*percent
          }
        }
      }
    }
  }
  return(y)
}

##### this is used to get true icc & true ccc
sig.a<-0.7 ### this all need to be changed due to true value each time 
sig.e<-2.8
rho<-0.8
#### the way to calculate true ICC and CCC is according to reference and dataset type
true1<-function(sig.a,sig.e){
  icc1t<-sig.a/(sig.a+(sig.e/(2*2)))
  ccc.inter<-sig.a/(sig.a+sig.e)
  print(icc1t)
  print(ccc.inter)
}
true1(sig.a,sig.e)

true2<-function(sig.a,sig.e,rho){
  sigma<-1-rho^2*1
  icc2t<-sig.a/(sig.a+((1-rho^2)+sig.e)/(2*2))
  ccc.total<-(sig.a+sig.a)/(sig.a*2+sig.e+sigma)## here is special, because I add replication of a as a reading 
  print(icc2t)
  print(ccc.total)
}
true2(sig.a,sig.e,rho)


##################################################################  be written before formal function
#### it can output a row of table for need of this research

row<-function(sig.a,sig.e,rho,trueicc,trueccc,type,exp,percent,outlier){
  icc<-NULL
  rankicc<-NULL
  ccc<-NULL
  rankccc<-NULL
  avg.icc<-rep(0,100)
  avg.rankicc<-rep(0,100)
  low.icc<-rep(0,100)
  up.icc<-rep(0,100)
  low.rankicc<-rep(0,100)
  up.rankicc<-rep(0,100)
  cov.icc<-rep(0,100)
  cov.rankicc<-rep(0,100)
  
  ### similar to ccc related vectors
  #..........
  
  for(i in 1:100){
    a<-rnorm(100,0,sig.a)### subject effect
    e<-rnorm(400,0,sig.e)### reading effect
    
    #### change Y simualation dataset
    if (type=="ICC1"){
      if (exp==F){        
        ybi11<-a+e[1:100]
        ybi12<-a+e[101:200]
        ybi21<-a+e[201:300]
        ybi22<-a+e[301:400]
      }
      if(exp==T){
        ybi11<-exp(a+e[1:100])
        ybi12<-exp(a+e[101:200])
        ybi21<-exp(a+e[201:300])
        ybi22<-exp(a+e[301:400])
      }
      
      ##### exponential should be done before outlier.  if both exist
      if(!is.null(outlier)){
        y1<-genoutlier(type,outlier,percent,ybi11,ybi12,ybi21,ybi22)
        ybi11<-y1[,1]
        ybi12<-y1[,2]
        ybi21<-y1[,3]
        ybi22<-y1[,4]
      }
      
      icc<-ICC(as.data.frame(cbind(ybi11,ybi12,ybi21,ybi22))) # get icc1
      rrr<-rank(c(ybi11,ybi12,ybi21,ybi22))
      ryb<-as.data.frame(cbind(rrr[1:100],rrr[101:200],rrr[201:300],rrr[301:400]))
      rankicc<-ICC(ryb)   # get rank icc1
      
      ob1<-c(ybi11,ybi12)
      ob2<-c(ybi21,ybi22)
      ccc<-epi.ccc(x=ob1,y=ob2) # get ccc.intra
      rankccc<-epi.ccc(x=rrr[1:200],y=rrr[201:400]) # get rank ccc.intrta
    }
    
    
    if (type=="ICC2"){
      var=c(1,1)
      mu=c(0,0)
      Sigma<-matrix(c(var[1],rho,rho,var[2]),ncol=2)
      bij1<-mvrnorm(100,mu,Sigma) 
      bij2<-mvrnorm(100,mu,Sigma) 
      
      
      ##  condition of if exponential and if outlier is similar to those of ICC1
      #............
      ### same way to get ranked ICC and CCC as well as original ICC and CCC
      #..............
     
    }
    
    ############### this is to get some icc vector
    avg.icc[i]<-icc[[1]][2][[1]][4] 
    avg.rankicc[i]<-rankicc[[1]][2][[1]][4]
    low.icc[i]<-icc[[1]][7][[1]][4]  
    up.icc[i]<-icc[[1]][8][[1]][4]
    low.rankicc[i]<-rankicc[[1]][7][[1]][4]
    up.rankicc[i]<-rankicc[[1]][8][[1]][4]
    
    if(trueicc>low.icc[i] & trueicc<up.icc[i]){   ### get probability
    }
    if(trueicc>low.rankicc[i] & trueicc<up.rankicc[i]){
      cov.rankicc[i]<-trueicc
    }
    
    
    ####### this is to get some CCC vectors
    avg.ccc[i]<-ccc[[1]][1][[1]][1]
    avg.rankccc[i]<-rankccc[[1]][1][[1]][1]
    low.ccc[i]<-ccc[[1]][2][[1]][1]
    up.ccc[i]<-ccc[[1]][3][[1]][1]
    low.rankccc[i]<-rankccc[[1]][2][[1]][1]
    up.rankccc[i]<-rankccc[[1]][3][[1]][1]
    
   ### if statement to get covarage probability of original and ranked CCC
   #............
   
  } ################################################################ ### for loop should end here
  
  ######## this can get result of a row of table
  sicc<- iccsum(avg.icc,avg.rankicc,cov.icc,cov.rankicc,low.icc,up.icc,low.rankicc,up.rankicc,trueicc)
  sccc<- cccsum(avg.ccc,avg.rankccc,cov.ccc,cov.rankccc,low.ccc,up.ccc,low.rankccc,up.rankccc,trueccc)
  
  return(c(sicc,sccc))
  
}#### this is the end of function row

##############################################################      apply
###  (1)trueicc=0.8
## 
set.seed(0916)
row1<-row(sig.a,sig.e,rho,0.8,0.5,type="ICC1",exp=F,percent=NULL,outlier=NULL)
set.seed(0916)
row2<-row(sig.a,sig.e,rho,0.8,0.5,type="ICC1",exp=F,percent=10,outlier="reading")
set.seed(0916)
row3<-row(sig.a,sig.e,rho,0.8,0.5,type="ICC1",exp=F,percent=5,outlier="reading")
set.seed(0916)
row5<-row(sig.a,sig.e,rho,0.8,0.5,type="ICC1",exp=F,percent=10,outlier="subject")
set.seed(0916)
#........................................................
