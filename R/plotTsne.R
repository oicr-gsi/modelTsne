plotTsne<-function(model=NULL,method="plotly",dims1.bw="nrd0", title=NULL, group=NULL){
  if(class(model)!="modelTsne"){
    stop("plotTsne: a valid modelTsne object must be provided")
  }
  df<-model$df.tsne
  if(!is.null(group)){
    if(length(group)!=nrow(df)){
      stop("plotTsne: group contains an incorrect number of elements")
    }else{
      df$group<-group
    }
    
  }
  dims<-model$dims
  p<-NULL
  if(method=="plotly"){
    #Color palette from the 'scales' package matching ggplot's. The default palette in Plotly produces warnings if no. groups <3 or >8
    col.palette=hue_pal()(length(unique(df$group)))
    
    if(dims==1){
      gg<-ggplot(df, aes(x = tsne.dim1)) +
        stat_density(bw=dims1.bw,aes(group = group, color = group),position="identity",geom="line") +
        theme_light() +
        ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
      p<-ggplotly(gg)
    }
    if(dims==2){
      p<-plot_ly(x=df$tsne.dim1, y=df$tsne.dim2, type="scatter", mode="markers", color=df$group,text = df$sample, colors = col.palette) %>%
        layout(
          title = title,
          xaxis = list(title="tsne.dim1",zeroline = FALSE),
          yaxis = list(title="tsne.dim2",zeroline = FALSE)
        )
    }
    if(dims==3){
      p<-plot_ly(x=df$tsne.dim1, y=df$tsne.dim2, z=df$tsne.dim3, type="scatter3d", mode="markers", color=df$group,text = df$sample, colors = col.palette) %>%
        layout(
          title = title,
          scene = list(
            xaxis = list(title = "tsne.dim1"),
            yaxis = list(title = "tsne.dim2"),
            zaxis = list(title = "tsne.dim3")
          ))
    }
  }
  if(method=="ggplot"){
    if(dims==1){
      p<-ggplot(df, aes(x = tsne.dim1)) +
        stat_density(bw=dims1.bw,aes(group = group, color = group),position="identity",geom="line") +
        theme_light() +
        ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
    }
    if(dims==2){
      p<-ggplot(df, aes(x=tsne.dim1, y=tsne.dim2, color=group)) +geom_point() + theme_light() +
        ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
    }
    if(dims==3){
      stop("modelTsne::plotTsne: error, ggplot does not support 3d scatterplots")
    }
  }
  return(p)
}
