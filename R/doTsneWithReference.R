doTsneWithReference<-function(data,model,group=NULL){
  if(!is.null(group)){
    if(ncol(data)!=length(group)){
      stop("methTsne::doTsneWithReference error: ncol(data)!=length(group)")
    }
  }else{
    group<-rep("unknown_sample",nrow(data))
  }
  list.doTsneWithReference<-list()
  for(pos.sample in 1:ncol(data)){
    data.model<-model$data
    data.sample<-data[row.names(data) %in% row.names(data.model),pos.sample,drop=FALSE]
    data.model<-data.model[row.names(data.model) %in% row.names(data.sample),]
    message("Probes found in sample: ",nrow(data.model),"/",nrow(model$data))
    data.merge<-cbind(data.model,data.sample)
    for(c in 1:ncol(data.merge)){
      col<-data.merge[,c]
      col.mean<-mean(as.numeric(col),na.rm=TRUE)
      col[!is.finite(as.numeric(col))]<-col.mean
      data.merge[,c]<-col
    }
    message("running t-SNE...",pos.sample,": ",colnames(data.sample)[1])
    row.new<-matrix(rep(0,model$dims),nrow=1,ncol=model$dims)
    colnames(row.new)<-paste0("tsne.dim",1:model$dims)
    r<-as.matrix(rbind(model$df.tsne[,c(2:(1+model$dims)),drop=FALSE],row.new))
    r2<-as.numeric(r[,1])
    if(model$dims>1){
      for(i in 2:model$dims){
        r2<-cbind(r2,as.numeric(r[,i]))
      }
    }
    r<-as.matrix(r2)
    tsne_out <- Rtsne(t(data.merge),perplexity=model$perplexity,check_duplicates = FALSE,Y_init=r,dims=model$dims,eta=model$eta,theta=model$theta)
    group.merge<-c(as.character(model$df.tsne$group),group[pos.sample])
    df.tsne<-data.frame(sample=colnames(data.merge),tsne_out$Y,group=group.merge)
    colnames(df.tsne)[c(2:(model$dims+1))]<-paste0("tsne.dim",c(1:model$dims))
    doTsneWithReference<-list(df.tsne=df.tsne)
    class(doTsneWithReference)<-"doTsneWithReference"
    list.doTsneWithReference[[length(list.doTsneWithReference)+1]]<-doTsneWithReference
  }
  names(list.doTsneWithReference)<-colnames(data)
  df.tsneFinal<-model$df.tsne
  for(pos.sample in 1:length(list.doTsneWithReference)){
    df.tsne.curr<-list.doTsneWithReference[[pos.sample]]$df.tsne
    df.tsneFinal<-rbind(df.tsneFinal,df.tsne.curr[nrow(df.tsne.curr),])
  }
  modelFinal<-model
  modelFinal$df.tsne<-df.tsneFinal
  modelFinal[[length(modelFinal)+1]]<-list.doTsneWithReference
  names(modelFinal)[length(modelFinal)]<-"list.doTsneWithReference"
  return(modelFinal)
}
