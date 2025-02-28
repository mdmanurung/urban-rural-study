require(vegan)

`sites.long` <-
  function(x, env.data=NULL){
    if (is.null(env.data)) {
      result <- data.frame(axis1=x$sites[, 1], axis2=x$sites[, 2], labels=rownames(x$sites))
    }else{
      result <- data.frame(cbind(env.data, axis1=x$sites[, 1], axis2=x$sites[, 2], labels=rownames(x$sites)))
    }
    return(result)
  } 



`ordiellipse.long` <-
  function(z, grouping.name="Grouping") {
    # copied from vegan 2.5-6 (a function that is not exported)
    `veganCovEllipse2` <-
      function(cov, center = c(0,0), scale = 1, npoints = 100)
      {
        ## Basically taken from the 'car' package: The Cirlce
        theta <- (0:npoints) * 2 * pi/npoints
        Circle <- cbind(cos(theta), sin(theta))
        ## scale, center and cov must be calculated separately
        Q <- chol(cov, pivot = TRUE)
        ## pivot takes care of cases when points are on a line
        o <- attr(Q, "pivot")
        t(center + scale * t(Circle %*% Q[,o]))
      }
    grouping <- names(z)
    for (g in 1:length(grouping)) {
      g.arg <- z[[grouping[g]]]
      result.g <- data.frame(veganCovEllipse2(g.arg$cov, 
                                              center=g.arg$center, 
                                              scale=g.arg$scale, 
                                              npoints=100))
      result.g <- data.frame(cbind(rep(grouping[g], nrow(result.g)), result.g))
      names(result.g) <- c(grouping.name, "axis1", "axis2")
      result.g[, 1] <- as.factor(result.g[, 1])
      if (g == 1) {
        result <- result.g
      }else{
        result <- rbind(result, result.g)
      }
    }
    return(result)
  }    


`axis.long` <-
  function(w, choices=c(1, 2), cmdscale.model=FALSE, CAPdiscrim.model=FALSE){
    if (cmdscale.model==FALSE && CAPdiscrim.model==FALSE) {
      eigs <- NULL
      if ("cca" %in% class(w)) {
        if (is.null(w$CCA$eig)) {
          eigs.all <- w$CA$eig
        }else{
          eigs.all <- c(w$CCA$eig, w$CA$eig)          
        }
        eigs <- round(100 * eigs.all / w$tot.chi, digits=1)[choices]  
      }
      if ("monoMDS" %in% class(w)) {
        labels <- paste0("NMS", choices)
        return(data.frame(axis=c(1:2), ggplot=c("xlab.label", "ylab.label"), label=labels))
      }
      if ("wcmdscale" %in% class(w)) {
        eigs <- round(100 * w$eig / sum(w$eig), digits=1)[choices]
        names(eigs) <- paste0("WMDS", choices)        
      }
      if ("decorana" %in% class(w)) {eigs <- round(100 * w$evals / sum(w$evals), digits=1)[choices]}
      if (is.null(eigs)) {
        labels <- paste0("DIM", choices)
        return(data.frame(axis=c(1:2), ggplot=c("xlab.label", "ylab.label"), label=labels))       
      }
    }
    if (cmdscale.model==TRUE) {
      eigs <- round(100 * w$eig / sum(w$eig), digits=1)[choices]
      names(eigs) <- paste0("MDS", choices)        
    }
    if (CAPdiscrim.model==TRUE) {    
      eigs <- round(100 * w$F / sum(w$F), digits=1)[choices]
      names(eigs) <- names(data.frame(w$x))[choices]       
    }
    xlab.label <- paste0(names(eigs)[1], " (", eigs[1], "%)")
    ylab.label <- paste0(names(eigs)[2], " (", eigs[2], "%)")
    result <- data.frame(axis=c(1:2), ggplot=c("xlab.label", "ylab.label"), label=c(xlab.label, ylab.label))
    return(result)
  }


`centroids.long` <-
  function(y, grouping, FUN=mean, centroids.only=FALSE){
    gr.name <- rlang::enexpr(grouping)
    cent.means <- stats::aggregate(cbind(axis1, axis2) ~ grouping, data = y, FUN = FUN)
    names(cent.means) <- c(gr.name, "axis1c", "axis2c")
    cent.means$Centroid <- cent.means[, 1]
    if (centroids.only == TRUE) {return(cent.means)}
    by.name <- names(y)[names(y) == gr.name]
    result <- dplyr::full_join(y, cent.means, by=by.name)
    return(result)
  } 