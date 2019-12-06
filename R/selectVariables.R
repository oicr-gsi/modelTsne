doANOVA<-function(group,data.row){
  if(max(data.row, na.rm=TRUE) - min(data.row, na.rm=TRUE) != 0){
    names(group)<-names(data.row)
    data.sub<-data.frame(variable=data.row,group=group)
    colnames(data.sub)<-c("variable","group")
    data.sub<-data.sub[is.finite(data.sub$variable),]
    p.value<-summary(aov(variable~group,data=data.sub))[[1]]['group','Pr(>F)']
  }else{
    p.value<-999
  }
  return(p.value)
}

doLogistic<-function(group,data.row){
  if(max(data.row, na.rm=TRUE) - min(data.row, na.rm=TRUE) != 0){
    names(group)<-names(data.row)
    data.sub<-data.frame(variable=data.row,group=group)
    colnames(data.sub)<-c("variable","group")
    data.sub<-data.sub[is.finite(data.sub$variable),]
    fit.binomial<-glm(formula=as.formula("group~variable"),data=data.sub,family=binomial)
    p.value<-summary(fit.binomial)$coefficients['variable','Pr(>|z|)']
  }else{
    p.value<-999
  }
  return(p.value)
}

chop<-function(df,block.size=10000){
  if(!is.data.frame(df)){df<-as.data.frame(df)}
  if(block.size>=nrow(df)){
    block.size<-nrow(df) %/% 2 #integer division
  }
  no.blocks<-nrow(df) %/% block.size #integer division
  block.map<-NULL
  c<-1
  for(i in 1:no.blocks){
    block.map<-c(block.map,rep(c,block.size))
    c<-c+1
  }
  block.map<-c(block.map,rep(c,nrow(df)-(block.size*no.blocks)))
  chunks<-split(df, block.map)
  return(chunks)
}

selectVariables<-function(data=NULL,no.variables=10000,method="IQR",group=NULL,filter="none",threads=1,file.stats=NULL){
  if(no.variables>=nrow(data)){
    message("\tno.variables >= variables in data. All the variables supplied in data were returned.")
    return(data)
  }
  if(!method %in% c("IQR","ANOVA","SD","logistic")){
    stop("\tAccepted methods are ANOVA, IQR, logistic and SD")
  }
  if(method=="ANOVA" | method=="logistic"){
    if(is.null(group) | length(group)!=ncol(data)){
      msg<-paste0("\tThe method ",method," requires a valid \"group\" parameter")
      stop(msg)
    }
  }
  message("modelTsne::selectVariables")
  message("\tdata: ",ncol(data)," samples, ",nrow(data)," variables")
  message("\tmethod: ",method)
  if(method=="ANOVA" | method=="logistic"){
    message("\tgroup: ",paste(unique(group),collapse=","))
  }
  message("\tfilter: ",filter)
  if(filter=="meth_array"){
    data<-data[row.names(data) %in% probes.intersect,]
    message("\tSelected variables present in both 450k and EPIC arrays")
    message("\tRemoved variables mapped to chromosomes X, Y, and frequent SNPs positions")
    pass<-apply(data,1,function(x){sum(is.finite(x))/ncol(data)>0.8})
    data<-data[pass,]
    if(nrow(data)<=no.variables){
      message("\tWARNING: Filtering reduced the number of variables below the requested no.variables")
      return(data)
    }
  }
  data<-data[order(row.names(data)),]
  data<-apply(data,2,function(x) {x[!is.finite(x)]<-NA;return (x)}) #formats different types of missing values (-Inf and other)
  data.sel<-NULL
  message("\tRunning with ",threads," cores")
  if(method %in% c("IQR","SD")){
    myFunc<-IQR
    if(method=="SD"){myFunc<-sd}
    if(threads==1){
      result<-foreach(chunks=chop(data),.combine='c') %do% { apply(chunks,1,function(x) myFunc(x,na.rm = TRUE))}
    }else{
      cl <- makeCluster(threads)
      registerDoParallel(cl)
      result<-foreach(chunks=chop(data),.combine='c') %dopar% { apply(chunks,1,function(x) myFunc(x,na.rm = TRUE))}
      stopCluster(cl)
    }
    threshold.result<-sort(result,decreasing=TRUE)[no.variables]
    data.sel<-data[result>=threshold.result,]
  }
  
  if(method %in% c("ANOVA","logistic")){
    myFunc<-doANOVA
    if(method=="logistic"){myFunc<-doLogistic}
    if(threads==1){
      p.value<-foreach(chunks=chop(data),.combine='c') %do% { apply(chunks,1,function(x) myFunc(group,x))}
    }else{
      cl <- makeCluster(threads)
      registerDoParallel(cl)
      p.value<-foreach(chunks=chop(data),.combine='c') %dopar% { apply(chunks,1,function(x) myFunc(group,x))}
      stopCluster(cl)
    }
    cutoff<-sort(p.value)[no.variables]
    data.sel<-data[p.value<=cutoff,]
    
    if(!is.null(file.stats)){
      df<-data.frame(sample=row.names(data),p.value)
      groups<-unique(group)
      for(group.curr in groups){
        df<-cbind(df,rowMeans(data[,group==group.curr],na.rm=TRUE))
      }
      colnames(df)[3:(length(groups)+2)]<-paste0("avg_",groups)
      write.csv(df,file=file.stats,row.names=F)
    }
  }
  
  message("\tOutput: ",ncol(data.sel)," samples, ",nrow(data.sel)," variables")
  return(data.frame(data.sel))
}
