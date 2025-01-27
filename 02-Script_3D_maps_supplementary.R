
library(sf)
depth_layers<-read.csv('depth_layers.csv',sep=";", row.names=1)
depth_layers[,1]<-factor(est_prof[,1],levels=c('SRF','DCM','250M', '900M','1200M','2300M'),labels=c('SRF','DCM','250M', '900M','1200M','2300M'))

## Loading  arguments for the 3D Maps
args3<-readRDS('args_map3D.rds')
coords<-depth_layers[,c("coordN", 'coordE')]
cruise<-depth_layers$cruise
map_letters<-  rep(LETTERS[1:5],5)
depth<-depth_layers['sample_depth']



for( l in seq_along(indicator_taxa)) {


  var<-indicator_taxa[l]
  print(paste0(l,"/",length(indicator_taxa),"-",var))
  data<-family_sqrt4[,var,drop=F]
  data$hc5<-hc5
  data$size<-scales::rescale(data[,var],c(0.5,5))
  dl<-split(data,cruise)
  file=paste0(l,"-",map_letters[[l]],"-",var,".png")
  png(file,height=13,width=30,res=500,units="cm")
  par(mfrow=c(1,2))
  par(mar=c(0,0,0,0),xpd=T)
  hc<-hc5
  for(j in seq_along(dl)){
    if(j==1){
      par(mar=c(0,3,2,0),xpd=T)
      args3$legend<-F
    } else{
      par(mar=c(0,0,2,3),xpd=T)
      args3$legend<-T
    }
    dd<-dl[[j]]
    ids<-rownames(dd)
    attr(dd,"coords")<-coords[rownames(dd),]

    attr(dd,"breaks")<-round(seq(round(min(family_sqrt4[,var]),1),max(family_sqrt4[,var]),len=5),2)

    if(j==1){
      args3$main<-paste(var,"-","2019")
    } else{
      args3$main<-paste(var,"-","2021")
    }

    ids=rownames(dd)
    co<-coords[ids,]
    colnames(co)<-c("x",'y')
    z_value<-depth[ids,]
    res<-data.frame(data=dd[,var],co,z_value,size=dd$size)
    colzones<-viridis::turbo(nlevels(hc5))
    res$color<-colzones[hc5[rownames(dd)]]

    res$size_value<-  res$data
    resdf<-res[c("x","y",'z_value',"color","size_value","size")]
    resdf<-resdf[order(resdf$z_value),]
    args3$addpoints<-resdf
    args3$br<- attr(dd,"breaks")
    pmat<-do.call(get_4D,args3)

  }
  if(length(attr(pmat,"legend"))>0){
    do.call(legend,attr(pmat,"legend"))
  }
  graphics.off()
}



