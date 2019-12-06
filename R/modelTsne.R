modelTsne<-function(data=NULL,group=NULL,perplexity=NULL,dims=2,eta=200,theta=0,persistent=FALSE){
  if(is.null(perplexity)){
    stop("A valid perplexity value must be provided")
  }
  if(is.null(group)){
    group<-rep("unknown",ncol(data))
  }
  mean.row<-rowMeans(data,na.rm=T)
  for(c in 1:ncol(data)){
    do.fix<-!is.finite(data[,c])
    data[do.fix,c]<-mean.row[do.fix]
  }
  tsne_out <- Rtsne(t(data),perplexity=perplexity,check_duplicates = FALSE,dims=dims,eta=eta,theta=theta)
  df.tsne<-data.frame(sample=colnames(data),tsne_out$Y,group)
  colnames(df.tsne)<-c("sample",paste0("tsne.dim",1:dims),"group")
  modelTsne<-list(df.tsne=df.tsne, data=NULL,perplexity=perplexity,dims=dims,eta=eta,theta=theta)
  if(persistent){modelTsne$data<-data}
  class(modelTsne)<-"modelTsne"
  return(modelTsne)
}
