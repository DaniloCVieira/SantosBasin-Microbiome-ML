
## Dendrogram functions
.midDend<-function (x){
  if (is.null(mp <- attr(x, "midpoint")))
    0  else mp
}
.memberDend<-function (x){
  r <- attr(x, "x.member")
  if (is.null(r)) {
    r <- attr(x, "members")
    if (is.null(r)) {
      r <- 1L
    }
  }
  r
}
add_ggtheme<-function(p,theme,base_size){
  p<-switch(theme,
            'theme_grey'={p+theme_grey(base_size)},
            'theme_bw'={p+theme_bw(base_size)},
            'theme_linedraw'={p+theme_linedraw(base_size)},
            'theme_light'={p+theme_light(base_size)},
            'theme_dark'={p+theme_dark(base_size)},
            'theme_minimal'={p+theme_minimal(base_size)},
            'theme_classic'={p+theme_classic(base_size)},
            'theme_void'={p+theme_void(base_size)})
  p
}
gg_dendrogram<-function(hc.clusters,hc.object, palette=viridis::turbo, labels=NULL, lwd=2, main="", xlab="Observations", ylab="Height", base_size=12, theme='theme_grey',offset_labels=-.1,xlab_adj=20, legend=c("outside","inside"),base_color="black",angle_label=0,log=F){
  {
    tree_data <- imesc_dendrogram_data(as.dendrogram(hc.object))
    cols<-palette(nlevels(hc.clusters))
    segtree<-segment_dd(tree_data)

    tree_data$labels$group<-hc.clusters[tree_data$labels$label]

    ranges<-data.frame(do.call(rbind,lapply(split(tree_data$labels,tree_data$labels$group),function(x) range(x$x)))
    )

    levs<-levels(hc.clusters)

    group<-lapply(1:nrow(segtree),function(i){
      pic<-as.numeric(which(apply(ranges,1,function(x) between(segtree$xend[i],x[1],x[2]))))
      if(length(pic)==0){pic<-NA}
      pic
    })

    segtree$group<-do.call(c,group)
    segtree$group<-factor(segtree$group,levels=levs)
    #segtree[tree_data$labels$x,"group"]<-hc.clusters[tree_data$labels$label]


    num_clusters<-nlevels(hc.clusters)
    heights <- sort(hc.object$height, decreasing = TRUE)
    cut_height <- heights[num_clusters - 1]


    #  e<-sapply(1:nrow(segtree) ,function(i) segtree$y[i]!=segtree$yend[i])
    #segtree$group[segtree$y>=cut_height&e]<-NA
    segtree2<-segtree
    segtree2$group[segtree$yend>=cut_height]<-NA
  }

  {
    labels_dend<-do.call(rbind,lapply(split(segtree2,segtree2$group),function(x){
      data.frame(x=x$x[which.max(x$y)],y=max(x$y),group=x$group[1])
    }))

    labels_dend$y<-min(labels_dend$y)
    if(!is.null(labels)){tree_data$labels$label<-labels}
    segtree_head<-segtree_a<-segtree[segtree$y>=cut_height,]

    minv<-segtree_a[segtree_a$yend<cut_height,]
    rownames(labels_dend)<-labels_dend$group
    rownames(minv)<-minv$group
    minv[rownames(labels_dend),"y"]<-labels_dend$y
    maxv<-segtree_a[sapply(1:nrow(segtree_a),function(i){
      segtree_a$x[i]!=segtree_a$xend[i]
    }),]
    maxv2<-segtree_a[segtree_a$y>=cut_height,]
    maxv<-maxv[maxv$y<=cut_height,]
    segtree_a<-minv
    segtree_b<-segtree[!is.na(segtree$group),]
    segtree_b<-segtree_b
    segtree_a<-rbind(segtree_a,maxv)
    maxv2$yend[which(maxv2$yend<=cut_height)]<-labels_dend$y[1]
    maxv2<-maxv2[maxv2$y!=cut_height,]

    segtree_head<-maxv2
    p<-ggplot()
    # if(FALSE)
    p<-p+geom_segment(
      data=segtree_b,
      lineend="square",
      aes(x = x,
          y = y,
          xend = xend,
          yend = yend,
          color=group),
      linewidth=lwd
    )
    #if(FALSE)
    p<-p+
      geom_segment(
        data=segtree_head,
        lineend="square",
        aes(x = x,
            y = y,
            xend = xend,
            yend = yend),
        color=base_color,
        linewidth=lwd


      )
    #if(FALSE)
    p<-p+geom_segment(
      data=segtree_a,
      lineend="square",
      aes(x = x,
          y = y,
          xend = xend,
          yend = yend,
          color=group),

      linewidth=lwd)

    p<-p+scale_color_manual(values=cols,name="")


  }

  p<-p+geom_text(aes(x=x,y=y,label=label),data=tree_data$labels,size=base_size*.2,nudge_y=offset_labels,angle=angle_label)



  p<-p+geom_label(data=labels_dend,aes(x,y,label=group),color=cols,nudge_y=-0.1)

  p<-add_ggtheme(p,theme,base_size)+xlab(xlab)+ylab(ylab)+ggtitle(main)
  if(isTRUE(log)){
    p<-p+scale_y_continuous(transform="log")
  }
  p
}
imesc_dendrogram_data<-function (x, type = c("rectangle", "triangle"), ...) {
  leaflab <- "perpendicular"
  center <- FALSE
  xlab <- ""
  ylab <- ""
  horiz <- FALSE
  xaxt <- "n"
  yaxt <- "s"
  nodePar <- NULL
  edgePar <- list()
  dLeaf <- NULL
  edge.root <- is.leaf(x) || !is.null(attr(x, "edgetext"))
  type <- match.arg(type)
  hgt <- attr(x, "height")
  if (edge.root && is.logical(edge.root)) {
    edge.root <- 0.0625 * if (is.leaf(x)) {
      1}    else {hgt}
  }
  mem.x <- .memberDend(x)
  yTop <- hgt + edge.root
  if (center) {
    x1 <- 0.5
    x2 <- mem.x + 0.5
  }  else {
    x1 <- 1
    x2 <- mem.x
  }
  xl. <- c(x1 - 1/2, x2 + 1/2)
  yl. <- c(0, yTop)
  if (edge.root) {
    if (!is.null(et <- attr(x, "edgetext"))) {
      my <- mean(hgt, yTop)
    }
  }
  ret <- plotNode(x1, x2, x, type = type, center = center,
                  leaflab = leaflab, dLeaf = dLeaf, nodePar = nodePar,                   edgePar = edgePar, horiz = FALSE)
  ret$segments <- as.data.frame(matrix(ret$segments, ncol = 4,
                                       byrow = TRUE, dimnames = list(NULL, c("x", "y", "xend",
                                                                             "yend"))))
  ret$labels <- cbind(as.data.frame(matrix(ret$labels$xy, ncol = 2,byrow = TRUE, dimnames = list(NULL, c("x", "y")))), data.frame(label = ret$labels$text))
  ret
}
plotNode<-function (x1, x2, subtree, type, center, leaflab, dLeaf, nodePar,edgePar, horiz = FALSE){
  ddsegments <- NULL
  ddlabels <- list()
  wholetree <- subtree
  depth <- 0L
  llimit <- list()
  KK <- integer()
  kk <- integer()
  repeat {
    inner <- !is.leaf(subtree) && x1 != x2
    yTop <- attr(subtree, "height")
    bx <- plotNodeLimit(x1, x2, subtree, center)
    xTop <- bx$x
    depth <- depth + 1L
    llimit[[depth]] <- bx$limit
    hasP <- !is.null(nPar <- attr(subtree, "nodePar"))
    if (!hasP) {
      nPar <- nodePar}

    Xtract <- function(nam, L, default, indx) rep(if (nam %in%
                                                      names(L)) {L[[nam]] }else {default}, length.out = indx)[indx]
    asTxt <- function(x) if (is.character(x) || is.expression(x) ||is.null(x)) {x}  else {as.character(x)}
    i <- if (inner || hasP) {
      1} else {2}
    if (!is.null(nPar)) {
      pch <- Xtract("pch", nPar, default = 1L:2, i)
      cex <- Xtract("cex", nPar, default = c(1, 1), i)
      col <- Xtract("col", nPar, default = par("col"),
                    i)
      bg <- Xtract("bg", nPar, default = par("bg"), i)
      points(if (horiz)
        cbind(yTop, xTop)
        else cbind(xTop, yTop), pch = pch, bg = bg, col = col,
        cex = cex)
    }
    if (leaflab == "textlike")
      p.col <- Xtract("p.col", nPar, default = "white",
                      i)
    lab.col <- Xtract("lab.col", nPar, default = par("col"),
                      i)
    lab.cex <- Xtract("lab.cex", nPar, default = c(1, 1),
                      i)
    lab.font <- Xtract("lab.font", nPar, default = par("font"),
                       i)
    lab.xpd <- Xtract("xpd", nPar, default = c(TRUE, TRUE),
                      i)
    if (is.leaf(subtree)) {
      if (leaflab == "perpendicular") {
        if (horiz) {
          X <- yTop + dLeaf * lab.cex
          Y <- xTop
          srt <- 0
          adj <- c(0, 0.5)
        }        else {
          Y <- yTop - dLeaf * lab.cex
          X <- xTop
          srt <- 90
          adj <- 1
        }
        nodeText <- asTxt(attr(subtree, "label"))
        ddlabels$xy <- c(ddlabels$xy, X, 0)
        ddlabels$text <- c(ddlabels$text, nodeText)
      }
    }    else if (inner) {
      for (k in seq_along(subtree)) {
        child <- subtree[[k]]
        yBot <- attr(child, "height")
        if (getOption("verbose"))
          cat("ch.", k, "@ h=", yBot, "; ")
        if (is.null(yBot))
          yBot <- 0
        xBot <- if (center) {
          mean(bx$limit[k:(k + 1)])}else {bx$limit[k] + .midDend(child)}
        hasE <- !is.null(ePar <- attr(child, "edgePar"))
        if (!hasE)
          ePar <- edgePar
        i <- if (!is.leaf(child) || hasE) {
          1}  else {2}
        col <- Xtract("col", ePar, default = par("col"),
                      i)
        lty <- Xtract("lty", ePar, default = par("lty"),
                      i)
        lwd <- Xtract("lwd", ePar, default = par("lwd"),
                      i)
        if (type == "triangle") {
          ddsegments <- c(ddsegments, xTop, yTop, xBot,
                          yBot)
        }        else {
          ddsegments <- c(ddsegments, xTop, yTop, xBot,
                          yTop)
          ddsegments <- c(ddsegments, xBot, yTop, xBot,
                          yBot)
        }
        vln <- NULL
      }
    }
    if (inner && length(subtree)) {
      KK[depth] <- length(subtree)
      if (storage.mode(kk) != storage.mode(KK))
        storage.mode(kk) <- storage.mode(KK)
      kk[depth] <- 1L
      x1 <- bx$limit[1L]
      x2 <- bx$limit[2L]
      subtree <- subtree[[1L]]
    }    else {
      repeat {
        depth <- depth - 1L
        if (!depth || kk[depth] < KK[depth])
          break
      }
      if (!depth)
        break
      length(kk) <- depth
      kk[depth] <- k <- kk[depth] + 1L
      x1 <- llimit[[depth]][k]
      x2 <- llimit[[depth]][k + 1L]
      subtree <- wholetree[[kk]]
    }
  }
  list(segments = ddsegments, labels = ddlabels)
}
plotNodeLimit<-function (x1, x2, subtree, center){
  inner <- !is.leaf(subtree) && x1 != x2
  if (inner) {
    K <- length(subtree)
    mTop <- .memberDend(subtree)
    limit <- integer(K)
    xx1 <- x1
    for (k in 1L:K) {
      m <- .memberDend(subtree[[k]])
      xx1 <- xx1 + (if (center) {
        (x2 - x1) * m/mTop
      }      else {
        m
      })
      limit[k] <- xx1
    }
    limit <- c(x1, limit)
  }  else {
    limit <- c(x1, x2)
  }
  mid <- attr(subtree, "midpoint")
  center <- center || (inner && !is.numeric(mid))
  x <- if (center) {
    mean(c(x1, x2))
  }  else {
    x1 + (if (inner) {
      mid
    }    else {
      0
    })
  }
  list(x = x, limit = limit)
}
segment_dd<-function (x) {
  x$segments
}


# Function to plot feature importance based on minimal depth analysis
# It highlights significant features and generates a plot for minimal depth distribution
# Parameters allow filtering by significance and customization of plot aesthetics
plot_feature_importance<-function(res_multi, sigs=T, sig.value=0.05,size_plot=10,
              min_no_of_trees = 0,
              mean_sample = "top_trees",
              mean_scale = FALSE,
              newcolhabs=newcolhabs,
              palette=palette,
              size_mmd=5){
  multi_imps=res_multi[[2]]
  multi_imps$col<-"gray"
  min_depth_frame=res_multi[[1]]


  multi_imps<-multi_imps[order(multi_imps$mean_min_depth),]
  multi_imps[which(multi_imps$p_value<sig.value),"col"]<-"red"


  BoldImportance<-multi_imps[which(multi_imps$p_value<sig.value),]


  sig_vars<-which(min_depth_frame$variable%in%BoldImportance$variable)

  if(isTRUE(sigs)){min_depth_frame_sigs<-min_depth_frame[sig_vars,]
  } else{    min_depth_frame_sigs<-min_depth_frame
  }

  if(isTRUE(sigs)){pick<-nrow(BoldImportance)} else{pick=sigs}




  MDD.plot<-plot_mdd(min_depth_frame_sigs, k=pick,
                     min_no_of_trees = min_no_of_trees,
                     mean_sample = mean_sample,
                     mean_scale = mean_scale,
                     newcolhabs=newcolhabs,palette=palette,size_mmd=size_mmd,size_plot=size_plot)


  return(MDD.plot)
}

# Function to create a plot for the distribution of minimal depth values
# Computes the mean minimal depth and overlays it on the histogram of feature occurrences in trees
# Allows filtering of variables based on a minimum number of occurrences in trees
# Parameters support custom color palettes and plot size adjustments

plot_mdd<-function (min_depth_frame, k = 10, min_no_of_trees = 0, mean_sample = "top_trees",mean_scale = FALSE, mean_round = 2, main = "Feature Importance",newcolhabs,palette="viridis",
                    size_mmd=5,size_plot=10) {
  minimal_depth <- NULL
  mean_minimal_depth_label <- NULL
  mean_minimal_depth <- NULL
  if (any(c("randomForest", "ranger") %in% class(min_depth_frame))) {
    min_depth_frame <- randomForestExplainer::min_depth_distribution(min_depth_frame)
  }

  min_depth_count_list <- randomForestExplainer:::min_depth_count(min_depth_frame)
  min_depth_means <- randomForestExplainer:::get_min_depth_means(min_depth_frame, min_depth_count_list,
                                                                 mean_sample)
  frame_with_means <- merge(min_depth_count_list[[1]], min_depth_means)
  frame_with_means[is.na(frame_with_means$minimal_depth), "count"] <- frame_with_means[is.na(frame_with_means$minimal_depth),
                                                                                       "count"] - min(frame_with_means[is.na(frame_with_means$minimal_depth),
                                                                                                                       "count"])
  if (mean_scale == TRUE) {
    frame_with_means$mean_minimal_depth <- (frame_with_means$mean_minimal_depth -
                                              min(frame_with_means$mean_minimal_depth))/(max(frame_with_means$mean_minimal_depth) -
                                                                                           min(frame_with_means$mean_minimal_depth))
  }
  frame_with_means$mean_minimal_depth_label <- (frame_with_means$mean_minimal_depth -
                                                  min(frame_with_means$mean_minimal_depth))/(max(frame_with_means$mean_minimal_depth) -
                                                                                               min(frame_with_means$mean_minimal_depth)) * max(min_depth_count_list[[2]]$no_of_occurrences)
  variables <- min_depth_count_list[[2]][min_depth_count_list[[2]]$no_of_occurrences >=
                                           min_no_of_trees, "variable"]
  frame_with_means <- frame_with_means[frame_with_means$variable %in%
                                         variables, ]
  frame_with_means <- within(frame_with_means, variable <- factor(variable,
                                                                  levels = unique(frame_with_means[order(frame_with_means$mean_minimal_depth),
                                                                                                   "variable"])))
  data <- frame_with_means[frame_with_means$variable %in% levels(frame_with_means$variable)[1:min(k,
                                                                                                  length(unique(frame_with_means$variable)))], ]
  data$variable <- droplevels(data$variable)
  data_for_labels <- unique(data[, c("variable", "mean_minimal_depth",
                                     "mean_minimal_depth_label")])
  data_for_labels$mean_minimal_depth <- round(data_for_labels$mean_minimal_depth, digits = mean_round)


  #print(aggregate(data["count"],data["variable"],sum))
  # print(unique(data[order(data$count),"variable"]))


  #colors<-colors[mdepth]
  plot<- ggplot(data, aes(x = variable, y = count)) +
    geom_col(position = position_stack(reverse = TRUE),aes(fill = as.factor(minimal_depth)))+
    scale_fill_manual(values=viridis::viridis(nlevels(as.factor(data$minimal_depth))))+
    coord_flip() +
    scale_x_discrete(limits = rev(levels(data$variable))) +
    geom_errorbar(aes(ymin = mean_minimal_depth_label, ymax = mean_minimal_depth_label),linewidth = 1.5) +
    xlab("Variable") +
    ylab("Number of trees") +
    guides(fill = guide_legend(title = "Minimal depth")) +
    theme_bw(base_size = size_plot) +
    geom_label(data = data_for_labels, aes(y = mean_minimal_depth_label,
                                           label = mean_minimal_depth), size=size_mmd)
  if (!is.null(main)) {
    plot <- plot + ggtitle(main)
  }
  attr(plot,"sigs")<-  levels(data$variable)
  return(plot)
}
# Function to calculate feature importance and minimal depth distribution for a random forest model
# Uses the randomForestExplainer package to extract feature importance measures
# Returns a list containing the minimal depth distribution and the calculated feature importance

multipimp<-function(rf, mean_sample = "top_trees", measures = NULL){


  if(class(rf)[1]=="randomForest"){
    forest=rf } else {forest=rf$finalModel}

  multi_imps = randomForestExplainer::measure_importance(forest,mean_sample=mean_sample,measures=measures)
  indicadores<-as.character(multi_imps[order(multi_imps$mean_min_depth),][,"variable"])
  min_depth_frame <- randomForestExplainer::min_depth_distribution(forest)
  res_multi<-list(min_depth_frame,multi_imps)
  return(res_multi)
}



# Define a function to perform multivariate indicator species analysis
indicator_multipart<-function(data,hc,npic=5,func="IndVal.g",nperm=199){
    indicator_r.g = indicspecies::multipatt(data, hc, func=func, duleg=T,control=permute::how(nperm=nperm))
  attr(indicator_r.g,"y")<-hc
  indicator_r.g
}
# Define a summary function for multivariate indicator species analysis
imesc_summary_multipatt<-function (object, alpha = 0.05, minstat = NULL, At = NULL, Bt = NULL,indvalcomp = FALSE, ...){
  x <- object
  ncomb = ncol(x$str)  # Number of group combinations
  ncolsign = ncol(x$sign)  # Number of columns in the significance table
  nsps = nrow(x$sign)  # Total number of species
  # Print analysis details
  c("\n Multilevel pattern analysis")
  c("\n ---------------------------\n")
  c("\n Association function:", x$func)
  c("\n Significance level (alpha):", alpha)
  if (!is.null(minstat))
    c("\n Minimum statistic value (minstat):", minstat)
  if (x$func == "IndVal" || x$func == "IndVal.g") {
    if (!is.null(At))
      c("\n Minimum positive predictive value (At):",
        At)
    if (!is.null(Bt))
      c("\n Minimum sensitivity (Bt):", Bt)
  }
  c("\n\n Total number of species:", nsps)
  # Apply selection criteria for significant associations
  sel = !is.na(x$sign$p.value) & x$sign$p.value <= alpha
  if (!is.null(minstat))
    sel = sel & (x$sign$stat >= minstat)
  if (!is.null(Bt) && !is.null(x$B)) {
    for (i in 1:nrow(x$sign)) sel[i] = sel[i] && (x$B[i,
                                                      x$sign$index[i]] >= Bt)
  }
  if (!is.null(At) && !is.null(x$A)) {
    for (i in 1:nrow(x$sign)) sel[i] = sel[i] && (x$A[i,
                                                      x$sign$index[i]] >= At)
  }
  # Extract significant species and print summary
  a = x$sign[sel, ]
  c("\n Selected number of species:", nrow(a), "\n")
  cols = (ncolsign - 1):ncolsign
  if (indvalcomp && !is.null(x$B) && !is.null(x$A)) {
    As = numeric(nrow(x$sign))
    Bs = numeric(nrow(x$sign))
    for (i in 1:nrow(x$sign)) {
      As[i] = x$A[i, x$sign$index[i]]
      Bs[i] = x$B[i, x$sign$index[i]]
    }
    y = cbind(x$sign, As, Bs)
    cols = c(ncol(y) - 1, ncol(y), cols)
    names(y) = c(names(x$sign), "A", "B")
  }  else y = x$sign
  for (k in 1:(ncolsign - 4)) {
    c(" Number of species associated to", k, if (k == 1)
      "group:"
      else "groups:", sum(rowSums(a[, 1:(ncolsign - 3)]) ==
                            k), "\n")
  }
  c("\n List of species associated to each combination: \n")
  # Identify species associated with each group combination
  com_result<-list()
  for (i in 1:ncomb) {
    sel = x$sign$index == i & !is.na(x$sign$p.value) & x$sign$p.value <=
      alpha
    if (!is.null(minstat))
      sel = sel & (x$sign$stat >= minstat)
    if (!is.null(Bt) && !is.null(x$B)) {
      for (j in 1:nrow(x$sign)) sel[j] = sel[j] && (x$B[j,
                                                        x$sign$index[j]] >= Bt)
    }
    if (!is.null(At) && !is.null(x$A)) {
      for (j in 1:nrow(x$sign)) sel[j] = sel[j] && (x$A[j,
                                                        x$sign$index[j]] >= At)
    }
    m = y[sel, ]
    if (nrow(m) > 0) {
      c("\n Group", colnames(x$comb)[i], " #sps. ", nrow(m),
        "\n")
      m = m[order(m$stat, decreasing = TRUE), cols]
      com_result[[colnames(x$comb)[i]]]<-m
      #  printCoefmat(m, signif.stars = TRUE, signif.legend = FALSE,                   digits = 4, P.values = TRUE, has.Pvalue = TRUE)
    }
  }
  Signif <- symnum(x$sign$p.value, corr = FALSE, na = FALSE,
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***",
                                                                            "**", "*", ".", " "))
  c("---\nSignif. codes: ", attr(Signif, "legend"), "\n")
  # names(com_result)<-colnames(x$comb)
  return(com_result)
}
# Summarize results of multivariate indicator species analysis
summary_multipatt<-function(indicator_r.g){
  hcs<-attr(indicator_r.g,"y")
  y0<-hcs
  # sss<-summary(indicator_r.g)
  sss<-imesc_summary_multipatt(indicator_r.g)
  sssdf<-do.call(rbind,sss)
  hcs<-do.call(c, lapply(names(sss),function(x) rep(x,nrow(sss[[x]]))))
  indresult<-data.frame(group=hcs,var=do.call(c,lapply(sss,function(x) rownames(x))),sssdf)

  indresult$group<-factor(indresult$group,levels=levels(y0))
  indresult
}
# Filter top results from multivariate indicator species analysis
filter_multipart_result<-function(indicator_r.g,npic){
  sss<-imesc_summary_multipatt(indicator_r.g)

  top_result<-do.call(c,lapply(  seq_along(sss),function(i){

    x<-sss[[i]]
    vec<-rownames(na.omit(x[1:npic,]))
    vec1<-rep(i,length(vec))
    names(vec1)<-vec
    vec1
  }))
  top_result
}

# Crop and align two raster objects to the same extent and resolution
crop_raters<-function(my_rst1,my_rst2){
  library(raster)
  # Crop rasters to their overlapping extent
  r2 <- raster::crop(my_rst1,extent(my_rst2))
  r3<-raster::crop(my_rst2,extent(my_rst1))
  extent(r2)<-extent(r3)
  r3e<-r3@ncols*r3@nrows
  r2e<-r2@ncols*r2@nrows
  if(r3e!=r2e){

    r23<-list(r2=r2,r3=r3)

    min_res<-which.min(c(r2e,r3e))
    max_res<-which.max(c(r2e,r3e))
    r_change<-r23[[max_res]]
    r_base<-r23[[min_res]]
    coarser=min_res
    e <- extent(r_base)
    r_new <- raster(e, ncol=r_base@ncols, nrow=r_base@nrows)
    r_tochange<-rasterToPoints(r_change)

    r_new2 <- rasterize( r_tochange[,1:2], r_new, r_tochange[,3], fun=mean)
    r23[[max_res]]<-r_new2
    r2<-r23$r2
    r3<-r23$r3
  }
  s<-raster::stack(r2,r3)

  s  # Return a stack of aligned rasters

}
# Add 3D axes to a plot
add_3dAxes <- function(x, y, z, pmat, xlab, ylab, zlab,xadj=c(0.5,0.5,0),yadj=c(0.5,0.5,0),zadj=c(0.5,0.5,90),arrow_adj=1.001,arrow_len=1,cex.lab=1,font.axes.lab=2) {
  # Define axis ranges and tick marks
  min.x <- min(x, na.rm = TRUE)
  max.x <- max(x, na.rm = TRUE)
  min.y <- min(y, na.rm = TRUE)
  max.y <- max(y, na.rm = TRUE)
  min.z <- round(min(z, na.rm = TRUE), 3)
  max.z <- round(max(z, na.rm = TRUE), 3)
  x.axis <- seq(min.x, max.x, len = 5)
  y.axis <- seq(min.y, max.y, len = 5)
  z.axis <- seq(min.z, max.z, len = 5)

  # Draw X-axis and label
  x_labels_pos_start <- trans3d(min.x, min.y*arrow_adj, min.z, pmat)
  x_labels_pos_end <- trans3d(max.x, min.y*arrow_adj, min.z, pmat)
  arrows(x_labels_pos_start$x, x_labels_pos_start$y, x_labels_pos_end$x, x_labels_pos_end$y, lwd = 1, length = 0.1)

  label_X <- trans3d(mean(x, na.rm=T), min.y, min.z, pmat)
  text(label_X$x, label_X$y, labels = ylab, cex = cex.lab,srt=xadj[3], adj=c(xadj[1:2]),font=font.axes.lab)

  # Draw Y-axis and label
  y_labels_pos_start <- trans3d(min.x, min.y, min.z, pmat)
  y_labels_pos_end <- trans3d(min.x, max.y, min.z, pmat)
  arrows(y_labels_pos_start$x, y_labels_pos_start$y, y_labels_pos_end$x, y_labels_pos_end$y, lwd = 1, length = 0.1)
  label_Y <- trans3d(min.x, mean(y, na.rm=T), min.z, pmat)

  text(label_Y$x, label_Y$y, labels = xlab, cex = cex.lab,srt=yadj[3],adj=yadj[1:2],font=font.axes.lab)


  # Draw Z-axis and label
  z_labels_pos_start <- trans3d(min.x, min.y*arrow_adj, min.z*arrow_len, pmat)
  z_labels_pos_end <- trans3d(min.x, min.y*arrow_adj, max.z*(1-arrow_len), pmat)
  arrows( z_labels_pos_end$x, z_labels_pos_end$y,z_labels_pos_start$x, z_labels_pos_start$y,  lwd = 1, length = 0.1)
  label_Z <- trans3d(min.x, min.y, mean(z,na.rm=T), pmat)

  text(label_Z$x, label_Z$y, labels = zlab, cex = cex.lab,srt=zadj[3], adj=c(zadj[1:2]),font=font.axes.lab)
}

# Function to add polygons to a 3D plot
add_3dpolygon<-function(my_rst1,args_shape, pmat,which_shape="base_shape"){
  # Check if the shape argument is set to TRUE
  if(isTRUE(args_shape$shape)) {
    # Extract the shape attribute from the raster
    shape<-attr(my_rst1,which_shape)
    fill<-args_shape$fill
    if(isTRUE(fill)){
      fill_color<-adjustcolor(args_shape$color,args_shape$fillOpacity)
    } else{
      fill_color<-NA
    }

    if(!is.null(shape)){
      # Get polygon coordinates from the shape
      poly_list<-get_polygon_coordinates(shape)
      poly_list<-poly_list[ sapply(sapply(poly_list,dim),length)>0]
      poly_list<-poly_list[sapply(poly_list,length)>0]
      # Loop through each polygon and add it to the 3D plot
      lapply(poly_list, function(poly){

        polygon(trans3d(poly[,1] , y=poly[,2], z=args_shape$z,  pmat=pmat),col=fill_color, border=args_shape$border_col, lwd=args_shape$weight)
      })

    }


  }
}
# Function to extract polygon coordinates from shape data
get_polygon_coordinates<-function(shape){
  # Extract geometry and retrieve coordinates
  result<-lapply(shape$geometry,function(x) {
    lapply(x,function(xx){
      xx[[1]]
    })
  })
  # Flatten the nested list and return the result
  result<-unlist(result,recursive=F)
  result
}
# Function to generate labels for breaks in a range
break_labels<-function(br){
  start<-0
  res<-c()
  for(i in seq_along(br)[-1]){
    start<-start+1
    res[i-1]<-  paste(br[start],"-",br[i])
  }
  res
}

# Function to generate a 3D plot with multiple layers and optional legend
get_4D<-function(my_rst1,my_rst2=NULL,colors,exp=0.2,wlegend=20,hlegend=20,tictype='detailed',shape_args=NULL,addpoints=NULL,addtext=NULL,legend=T,slices=F,bg_points=T,theta=0,phi=40,r=1.73,d=1,xlab="",ylab="",
                 zlab="",main=NULL,plot.title_size=1,
                 axis.title_size=1,
                 axis.text_size=1,
                 y.intersp=1,
                 inset_x=-0.25,
                 inset_y=0,
                 x.intersp=1,
                 color_legend_point=NULL,
                 lines_list=NULL,
                 col_lines_labels="gray",
                 line_label_cex=.8,
                 line_adj_x=1,
                 line_srt =45,
                 line_label_adj=c(1,1),
                 colkey=F,
                 bty="b",
                 line_color="white",
                 line_lwd=0.1,
                 br=NULL,
                 leg_position="topright",
                 return_legend=F,
                 alpha_points=1,
                 light_points=0,
                 xadj=c(0.5,0.5,0),
                 yadj=c(0.5,0.5,0),
                 zadj=c(0.5,0.5,0),
                 arrow_adj,arrow_len=1,
                 cex.lab=1,
                 cex.leg.label=1.2,
                 font.axes.lab=2,
                 ...)
{



  col0<-colors
  my_rst2_0<-my_rst2
  # Crop and align the two input raster layers
  stack_rasters<-crop_raters(my_rst1,my_rst2)
  my_rst1<-stack_rasters[[1]]
  my_rst2=stack_rasters[[2]]
  # Convert rasters to matrices for plotting
  ex<-extent(my_rst1)
  r1<-raster::as.matrix(my_rst1,transpose=F)
  z1<-t(apply(r1,2,rev))
  r2<-raster::as.matrix(my_rst2,transpose=F)
  z2<-t(apply(r2,2,rev))
  x<-seq(ex[1],ex[2],len=nrow(z1))
  y<-seq(ex[3],ex[4],len=ncol(z1))
  colfunc = colorRampPalette(colors)
  color_apl<-colfunc(100)
  # Calculate center values for each facet of the raster for coloring
  z.facet.center <- (z2[-1, -1] + z2[-1, -ncol(z2)] + z2[-nrow(z2), -1] + z2[-nrow(z2), -ncol(z2)])/4
  # Range of the facet center on a 100-scale (number of colors)
  z.facet.range<-cut(z.facet.center, 100, include.lowest=T)

  # Create the 3D base plot using perspbox
  pmat<-plot3D::perspbox(x=x,y=y,z=z1, exp=exp,ticktype=tictype,theta=theta,phi=phi,r=r,d=d,xlab=xlab,ylab=ylab,zlim=range(z1[!is.na(z1)]),axes=F,
                         zlab=zlab, main=main,cex.main=plot.title_size/12,cex.axis=axis.text_size/12,cex.lab=axis.title_size/12)

  # Add 3D axes to the plot
  add_3dAxes(x=x,y=y,z=z1, pmat,  xlab=xlab, ylab=ylab, zlab=zlab,xadj=xadj,yadj=yadj,zadj=zadj,arrow_adj=arrow_adj,arrow_len=arrow_len,cex.lab=cex.lab,font=font.axes.lab)

  # Extract and add polygons for base and layer shapes
  base_args<-shape_args$base_shape_args
  layer_args<-shape_args$layer_shape_args
  add_3dpolygon(my_rst1,base_args, pmat,which_shape="base_shape")
  add_3dpolygon(my_rst1,layer_args, pmat,which_shape="layer_shape")


  # Add the second raster layer as a 3D surface with color mapping
  plot3D::persp3D(x=x,y=y,z=z1, colvar=z2,
                  col=colfunc(length(z2)),
                  border=NA,exp=exp,ticktype=tictype,theta=theta,phi=phi,r=r,d=d,xlab=xlab,ylab=ylab,add=T,clab=attr(my_rst2_0,'z_name'),colkey=F)

  # Additional points, text, and legend can be added based on user inputs
  if(!is.null(lines_list)){
    for(d in 1:length(lines_list)){
      line<-trans3d(lines_list[[d]]$x,lines_list[[d]]$y,lines_list[[d]]$z,pmat)
      lines(line,lwd=line_lwd,col=line_color)

    }
    ldeph<-do.call(rbind,lines_list)
    linexy<-do.call(rbind,lapply(split(ldeph,ldeph$z),function(xxx){
      xxx[which.min(rowSums(xxx[,1:2])),]
    }))
    text(trans3d(linexy$x*line_adj_x,linexy$y,linexy$z,pmat), labels=rownames(linexy),cex=line_label_cex,col=col_lines_labels,srt=line_srt,adj=line_label_adj)
  }

  if(!is.null(addpoints)){
    addpoints$color<-adjustcolor(colorspace::lighten(addpoints$color,light_points),alpha_points)

    if(isFALSE(bg_points)){
      points(trans3d(addpoints[,1] , y=addpoints[,2], z=addpoints[,3],  pmat=pmat), col=addpoints[,4], pch=16,cex=addpoints[,"size"])
    } else{
      pal<-colfunc(length(unique(addpoints[,3])))[cut(addpoints[,3],length(unique(addpoints[,3])))]
      points(trans3d(addpoints[,1] , y=addpoints[,2], z=addpoints[,3],  pmat=pmat),col=pal ,bg=addpoints[,4], pch=21,cex=addpoints[,"size"])
    }

  }

  if(!is.null(addtext)){
    text(trans3d(addtext$x , y=addtext$y, z=addtext$z,  pmat=pmat),col=addtext$color , cex=addtext$size/5, labels=addtext$label)

  }
  if(!is.null(addpoints)){

    par(xpd=T)
    if(isTRUE(legend)){

      if(length(unique(addpoints$color))>0){
        usr<-par("usr")

        zsize<-addpoints$size
        lev_legend<-levels(cut(zsize,br))
        lev_legend<-gsub("\\(","",lev_legend)
        lev_legend<-gsub("\\]","",lev_legend)
        x<-lev_legend[1]
        lev_legend<-sapply(lev_legend, function(x) {
          paste0(rev(strsplit(x,",")[[1]]),collapse=" - ")
        })
        y_legend<-seq(usr[2],usr[3],len=10)[1:length(lev_legend)]
        x_legend<-rep(usr[2],length(lev_legend))

        zval<-addpoints$size_value

        if(is.null(br)){
          pretty_z<-pretty(zval)
          r_br<-max(sapply(pretty_z,decimal_places))
          br<-round(seq(min(zval),max(zval),len=5),1)
        }


        cut_var<-cut(zval,breaks=br,include.lowest =T)

        libr<-split(zval,cut_var)


        libr<-libr[sapply(libr,length)>0]
        zval_legend<-break_labels(br)

        pt.cex=unique(seq(min(zsize,na.rm=T),max(zsize,na.rm=T),len=length(libr)))
        col=rev(unique(addpoints$color))
        if(!is.null(color_legend_point)){
          col<-color_legend_point
        }
        legend=zval_legend

        args_legend<-list(
          x=leg_position,legend=legend,pch=1,pt.cex=pt.cex,col=col,y.intersp=y.intersp,inset=c(inset_x,inset_y),bty="n",x.intersp=x.intersp,cex=cex.leg.label
        )
        if(!return_legend){
          legend(leg_position,legend=legend,pch=1,pt.cex=pt.cex,col=col,y.intersp=y.intersp,inset=c(inset_x,inset_y),bty="n",x.intersp=x.intersp,cex=cex.leg.label)
        }

      }
      attr(pmat,"legend")<-args_legend
    }



  }



  return(pmat)

}
