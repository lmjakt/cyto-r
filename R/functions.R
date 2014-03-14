## library MASS has a function kde2d which is a two-dimensional kernel density estimator..
## that seems a very useful thing to have..  using a bivariate normal kernel.
## cool!

## library sm also has a smoothing function (this time nonparametric?) sm.density
## which we can try in the future.

library(MASS)

vline <- function(x, col='black', lwd=1, lty=1){
  lines(c(x, x), c(0, 1e6), col=col, lwd=lwd)
}

hline <- function(y, col='black', lwd=1, lty=1){
  lines(c(0, 1e6), c(y, y), col=col, lwd=lwd)
}

compensate <- function(d, c1, c2, cf, auto.adjust=TRUE){
  for(i in 1:length(data)){
    c1.min <- min(d[[i]][,c1])
    d[[i]][,c1] <- d[[i]][,c1] - (cf * d[[i]][,c2])
    if(auto.adjust){
      c1.c.min <- min(d[[i]][,c1])
      d[[i]][,c1] <- d[[i]][,c1] + (c1.min - c1.c.min)
    }
  }
  d
}


## filter.simple takes a matrix containing minimum and maximum values
## for named columns in d. These are then used iteratively to remove
## values
filter.simple <- function(d, f){
  d.f <- d
  f.n <- rownames(f)
  for(i in 1:nrow(f)){
    d.f <- d.f[ as.vector(d.f[,f.n[i]]) > f[i,1], ]
    d.f <- d.f[ as.vector(d.f[,f.n[i]]) < f[i,2], ]
  }
  d.f
}

## this is for single gate defined for 2 columns.
## the gate is a list containing: points, xcol, ycol and a boolean indicating right handedness
## the gate must be ciruclar!! 
filter.2D <- function(d, gate){
  b <- rep(TRUE, nrow(d))
  if(gate$xlog) d[,gate$xcol] <- log10(d[,gate$xcol])
  if(gate$ylog) d[,gate$ycol] <- log10(d[,gate$ycol])
  for(i in 2:nrow(gate$points))
    b <- b & (is.right.of.m( gate$points[i-1,], gate$points[i,], d[,c(gate$xcol, gate$ycol)], gate$right.handed ))
  b
}

## combines a set of filters specified by a list of
## gates which are combined in an OR manner.
## (i.e. a set of convex polygons combined to allow concavity)
## returns the filtered data not the boolean vector
filter.2D.complex <- function(d, gates){
  b <- filter.2D(d, gates[[1]])
  if(length(gates) > 1){
    for(i in 2:length(gates))
      b <- (b | filter.2D(d, gates[[i]]))
  }
  d[b,]
}

filter.2D.complex.list <- function(d.l, gates){
  for(i in 1:length(d.l))
    d.l[[i]] <- filter.2D.complex(d.l[[i]], gates)
  d.l
}

## the gate defines its own logness:
## (this is to avoid rounding errors when converting between linear
## and log space)
plot.gates <- function(gates, xlog=FALSE, ylog=FALSE, col="gold", lwd=2, type='b'){
  if(!length(gates)) return(invisible(1))
  for(i in 1:length(gates)){
    if(xlog) gates[[i]]$points[,1] <- log10(gates[[i]]$points[,1])
    if(ylog) gates[[i]]$points[,2] <- log10(gates[[i]]$points[,2])
    points(gates[[i]]$points, type=type, col=col, lwd=lwd)
  }
}
filter.simple.list <- function(d.l, f){
  for(i in 1:length(d.l))
    d.l[[i]] <- filter.simple(d.l[[i]], f)
  d.l
}


## make a list of values from a list of matrices
## and a column name
make.v.list <- function(v, c, use.log=TRUE, filt=NA){
  l <- list()
  samples <- names(v)
  for(i in 1:length(v)){
    if(is.matrix(filt)){
      v[[i]] <- filter.simple( v[[i]], filt )
    }
    if(use.log){
      l[[i]] <-  log10(v[[i]][,c])
    }else{
      l[[i]] <-  v[[i]][,c]
    }
  }
  names(l) <- samples
  l
}

## take a named list of variables to calculate histograms forn
plot.histograms <- function(v.l, breaks=b, cex=0.5){
  al <- as.numeric(unlist(v.l))
  al.hist <- hist(al, breaks=breaks, plot=FALSE)
  h.density <- matrix(ncol=length(v.l), nrow=length(al.hist$mids))
  for(i in 1:length(v.l)){
    t.h <- hist(v.l[[i]], breaks=al.hist$breaks, plot=FALSE)
    h.density[,i] <- t.h$density
  }
  plot(al.hist$mids, h.density[,1], type='n', ylim=range(h.density))
  for(i in 1:ncol(h.density))
    points(al.hist$mids, h.density[,i], col=i, type='b', pch=i, cex=cex)
  legend('topright', legend=names(v.l), col=1:length(v.l), pch=1:length(v.l))
}

## take a list of dataframes, plot x vs y for each of them, using
## default colors (i.e. col=numeric value.. )
plot.xy <- function(v.l, xcol, ycol, filt=NA, use.log=TRUE, pch=20, cex=0.1){
  if(is.matrix(filt)){
    v.l <- filter.simple.list(v.l, filt)
  }
  xval <- c()
  yval <- c()
  for(i in 1:length(v.l)){
    if(use.log){
      v.l[[i]][,xcol] <- log10(v.l[[i]][,xcol])
      v.l[[i]][,ycol] <- log10(v.l[[i]][,ycol])
    }
    xval <- c(xval, v.l[[i]][,xcol])
    yval <- c(yval, v.l[[i]][,ycol])
  }
  
  plot(1,1, type='n', xlim=range(xval), ylim=range(yval), xlab=xcol, ylab=ycol)
  for(i in 1:length(v.l))
    points(v.l[[i]][,xcol], v.l[[i]][,ycol], pch=pch, cex=cex, col=i)
}

## take a list of data frames, and two column headers, return
## list of kde2d results
## uses global limits
make.kde2d.list <- function(d.l, xcol, ycol, n=50, use.log=TRUE){
  if(use.log){
    for(i in 1:length(d.l)){
      d.l[[i]][,xcol] <- log10(d.l[[i]][,xcol])
      d.l[[i]][,ycol] <- log10(d.l[[i]][,ycol])
    }
  }
  xlim <- range(d.l[[1]][,xcol])
  ylim <- range(d.l[[1]][,ycol])
  for(i in 2:length(d.l)){
    xlim <- range(c(xlim, range(d.l[[i]][,xcol])))
    ylim <- range(c(ylim, range(d.l[[i]][,ycol])))
  }
  de.l <- list()
  for(i in 1:length(d.l))
      de.l[[i]] <- kde2d(d.l[[i]][,xcol], d.l[[i]][,ycol], n=n, lims=c(xlim, ylim))

  names(de.l) <- names(d.l)
  de.l
}

                       
## kde2dplot
## Perspective plot and contour plot
## code originally from Romain Francois copied from:
## http://addictedtor.free.fr/graphiques/graphcode.php?graph=1
kde2dplot <- function(d,                 # a 2d density computed by kde2D
                      ncol=50,                   # the number of colors to use
                      zlim=c(0,max(z)),      # limits in z coordinates
                      nlevels=20,              # see option nlevels in contour
                      theta=30,                # see option theta in persp
                      phi=30){                   # see option phi in persp
  
  z   <- d$z
  nrz <- nrow(z)
  ncz <- ncol(z)
  
  couleurs  <- tail(topo.colors(trunc(1.4 * ncol)),ncol)
  fcol      <- couleurs[trunc(z/zlim[2]*(ncol-1))+1]
  dim(fcol) <- c(nrz,ncz)
  fcol      <- fcol[-nrz,-ncz]
  
  par(mfrow=c(1,2),mar=c(0.5,0.5,0.5,0.5))
  persp(d,col=fcol,zlim=zlim,theta=theta,phi=phi,zlab="density")
  
  par(mar=c(2,2,2,2))
  image(d,col=couleurs)
  contour(d,add=T,nlevels=nlevels)
  box()
}

## plot several density / contour plots next to each othern
## de.l contains a list of density estimates (i.e. lists containing x, y, z
## components
## specify:
##  1. number of columns (column.no)
##  2. number of colors (50)
##  3. number of contous (nlevels)
x.y.density <- function(de.l, column.no=3, ncol=50, nlevels=20, xlines=NULL, ylines=NULL){
  column.no <- min(length(de.l), column.no)
  row.no <- ceiling(length(de.l)/column.no)
  par(mfrow=c(row.no,column.no), mar=c(2,2,1,1))

  colors <- topo.colors(ncol)
#  colors <-  tail(topo.colors(trunc(1.4 * ncol)),ncol)
  zlim <- c(0, 0)
  ## make sure to use the same color scheme for each
  ## plot
  for(i in 1:length(de.l))
    zlim <- range(zlim, de.l[[i]]$z)
  for(i in 1:length(de.l)){
    xlim <- range(de.l[[i]]$x)
    ylim <- range(de.l[[i]]$y)
    image(de.l[[i]], col=colors)
    contour(de.l[[i]], add=TRUE, nlevels=nlevels, zlim=zlim)
    if(!(is.null(xlines))){
      for(j in 1:length(xlines))
        lines(c(xlines[j], xlines[j]), ylim)
    }
    if(!(is.null(ylines))){
      for(j in 1:length(ylines))
        lines(xlim, c(ylines[j], ylines[j]))
    }
  }
}
                        
read.files <- function(file.names){
  d.l <- list();
  d.l.names <- c()
  for(i in 1:length(file.names)){
    d.l[[i]] <- read.table(file.names[i], header=TRUE, sep="\t")
    d.l.names <- c(d.l.names, file.names[i])
  }
  names(d.l) <- d.l.names
  d.l
}

## log transform columns in all the dataframes
list.log <- function(d.l, columns){
  for(i in 1:length(d.l)){
    for(j in 1:length(columns))
      d.l[[i]][,columns[j]] <- log10(d.l[[i]][,columns[j]])
  }
  d.l
}

## do a principal components analysis the specified colums and sorts
pca.facs <- function(d.l, columns, use.log=TRUE){
  if(use.log)
    d.l <- list.log(d.l, columns)
  ## make a matrix containing all rows
  ## for the specified columns
  d.sub <- d.l[[1]][,columns]
  d.ranges <- matrix(nrow=length(d.l), ncol=2)
  d.ranges[1,] <- c(1, nrow(d.l[[1]]))
  for(i in 2:length(d.l)){
    d.ranges[i, 1] <- 1+nrow(d.sub)
    d.sub <- rbind(d.sub, d.l[[i]][,columns])
    d.ranges[i, 2] <- nrow(d.sub)
  }
  ## and then simply
  pc <- prcomp(d.sub)
  pc.res <- list(pca=pc, sample.rows=d.ranges, data=d.sub)
}


## a test for which side of a line a point lies

## p1 and p2 define a line running from p1 to p2
## q is a point which is tested for lying to the right
## of the vector defined by p1 and p2
is.right.of <- function(p1, p2, q){
  p.delta <- (p2 - p1)
  h <- sqrt(p.delta[1]^2 + p.delta[2]^2)
  ## the end result of a 2d rotation matrix is simply
  ## x' = xcos(a) - ysin(a)
  ## y' = xsin(a) + ycos(a)
  ## where a is the angle in a clockwise direction

  ## not sure quite why terms work when
  ## outside of the positive x and y values, but it does
  ## since for an angle to the y-axis:
  ## y / h = cos(a)
  ## x / h = sin(a)
  ## a rotation using those values should give
  ## a line lying on the y-axis.

  ## all we need to determine is the x' of q
  p.q <- (q - p1)
  ## x' = xcos(a) - ysin(a)
  ## and cos(a) = y / h  : sin(a) = x/h
  #cat(p.delta)
  return( ((p.q[1] * p.delta[2]/h) - (p.q[2] * p.delta[1]/h)) > 0)
}  

## a function similar to above, but takes a matrix of points to test
## q is a matrix of points (col1 = x, col2 = y)
## for details of how it works, see is.right.of above
is.right.of.m <-  function(p1, p2, q, check.right=TRUE){
  p.delta <- (p2 - p1)
  h <- sqrt(p.delta[1]^2 + p.delta[2]^2)
  p.q <- t(apply(q, 1, '-', p1))  ## ?? so complex
  if(check.right)
    return( ((p.q[,1] * (p.delta[2]/h)) - (p.q[,2] * (p.delta[1]/h))) > 0)
  ## we're actually looking for is.left.of...
  return( ((p.q[,1] * (p.delta[2]/h)) - (p.q[,2] * (p.delta[1]/h))) < 0)
}

## a function for getting mouse input..
## taken from : http://stat.ethz.ch/R-manual/R-devel/library/grDevices/html/getGraphicsEvent.html
## written by Duncan Murdoch
## (though I may modify it for playing arount with). ## to use this
## 1. savepar <- par()
## 2. dragplot( ) # with something that can be plotted
## press q to quit
## 3. getGraphicsEvent()  ## reset the handlers to null..  
## 4. par(savepar)
dragplot <- function(..., xlim=NULL, ylim=NULL, xaxs="r", yaxs="r") {
    plot(..., xlim=xlim, ylim=ylim, xaxs=xaxs, yaxs=yaxs)
    startx <- NULL
    starty <- NULL
    usr <- NULL

    devset <- function()
        if (dev.cur() != eventEnv$which) dev.set(eventEnv$which)
        
    dragmousedown <- function(buttons, x, y) {
        startx <<- x
        starty <<- y
        devset()
        usr <<- par("usr")
        eventEnv$onMouseMove <- dragmousemove
        NULL
    }
    
    dragmousemove <- function(buttons, x, y) {
        devset()
        deltax <- diff(grconvertX(c(startx,x), "ndc", "user"))
        deltay <- diff(grconvertY(c(starty,y), "ndc", "user"))
        plot(..., xlim=usr[1:2]-deltax, xaxs="i",
                  ylim=usr[3:4]-deltay, yaxs="i")
        NULL
    }
    
    mouseup <- function(buttons, x, y) {    
    	eventEnv$onMouseMove <- NULL
    }	
        
    keydown <- function(key) {
        if (key == "q") return(invisible(1))
        eventEnv$onMouseMove <- NULL
        NULL
    }
    
    setGraphicsEventHandlers(prompt="Click and drag, hit q to quit",
                     onMouseDown = dragmousedown,
                     onMouseUp = mouseup,
                     onKeybd = keydown)
    eventEnv <- getGraphicsEventEnv()
}

## a function to split a multiline gate such into gates that do not
## have any concavities, and which cover the original area.
## hmm, my thinking doesn't work; so at the moment this
## function DOES NOT WORK!!! DON'T USE !! 
split.2D.gate <- function(gate){
  check.hands <- function(gt, p, s=(p-1)){
    forward.hand <- is.right.of(gt$points[(s-1),], gt$points[(s),], gt$points[p,])
    closing.hand <- is.right.of(gt$points[s,], gt$points[p,], gt$points[1,])
    closing2.hand <- is.right.of(gt$points[p,], gt$points[1,], gt$points[2,])
    #cat(paste("check.hands", length(gt$points), " p: ", p, "\n"))
    #cat(paste("check.hands", forward.hand, closing.hand, closing2.hand, "\n"))
    c(forward.hand, closing.hand, closing2.hand)
  }
  equal.hands <- function(h1, h2){  ## because identical seems broken (not supposed to be, but)
    id <- TRUE
    for(k in 1:length(h2)){
      if(h1 != h2[k])
        id <- FALSE
    }
    id
  }
  right.handed <- gate$right.handed
  daughters <- list(gate);
  if(identical(gate$points[1,], gate$points[ nrow(gate$points), ]))  ## we need to avoid circularity
    daughters[[1]]$points <- gate$points[1:(nrow(gate$points)-1), ]

  i <- 0
  while(i < length(daughters)){
    i <- i+1
    d.points <- daughters[[i]]$points[1:2,]
    j <- 2
    while(j < nrow(daughters[[i]]$points)){
      j <- j + 1
      ## if we encounter an incorrect handedness.
      hands <- check.hands(daughters[[i]], j)
      cat(paste(j, ": "))
      cat(hands)
      cat("\n")
      cat(paste(j, "d.points:", nrow(d.points), "\n"))
      if(equal.hands(right.handed, hands)){
        d.points <- rbind(d.points, daughters[[i]]$points[j,])
        next;
      }
      s <- (j-1)  ## the first point of the new gate
      dd.points <- daughters[[i]]$points[s:j,] ## then increase j, until we find one thats ok
      while(j < nrow(daughters[[i]]$points)){
        j <- j + 1
        dd.points <- rbind(dd.points, daughters[[i]]$points[j,])
        hands <- check.hands(daughters[[i]], j, s)
        cat(paste(j, "d.points:", nrow(d.points), "\n"))
        if(equal.hands(right.handed, hands)){
          #dd.points <- rbind(dd.points, dd.points[1,]) ## make circular
          daughters[[length(daughters)+1]] <- list(points=dd.points,  right.handed=gate$right.handed, # a bit funny, but necessary
                                 xcol=gate$xcol, ycol=gate$ycol)
          d.points <- rbind(d.points, daughters[[i]]$points[j,])
          ## ??
          cat(paste(j, ": d.points ", nrow(d.points), " daughter", i, ":", nrow(daughters[[i]]$points), "\n"))
          if(j < nrow(daughters[[i]]$points))
            daughters[[i]]$points <- rbind(d.points, daughters[[i]]$points[ (j+1):nrow(daughters[[i]]$points), ] )
          cat(paste("and daughter now:", nrow(daughters[[i]]$points), "\n"))
          j <- nrow(d.points)
          break
        }
      }
    }
    cat(paste("And surely we didn't get here all of a sudden?", i, " and j is:", j,"\n"))
    daughters[[i]]$points <- d.points
  }
  dg <- list(parent=gate, daughters=daughters)
  dg
}

## a function to make a multiline gate
## arguments are
## df : a dataframe containing the values
## xcol, ycol: the names of the x and y columns
## xlog, ylog: whether to use log values for the data
## normally used with FSC and SSC, hence xlog and ylog
## defaulting to false
make.2D.gate <- function(df, xcol, ycol, xlog=FALSE, ylog=FALSE,
                         pch=1, cex=0.5){
  ## make a plot for xcol and ycol
  if(xlog)
    df[,xcol] <- log10(df[,xcol])
  if(ylog)
    df[,ycol] <- log10(df[,ycol])
  
  plot.graph <- function()  plot(df[,xcol], df[,ycol], xlab=xcol, ylab=ycol, pch=pch, cex=cex)
  plot.points <- function() {
    points(eventEnv$m.points, col="red", type='b', lwd=1)
  }
  ## this only works with x11(type="Xlib")
  ## check for a current device supporting events
  if(dev.cur() == 1 || length(dev.capabilities()$events) == 0)
    x11(type="Xlib")
  plot.graph()
                                        #  plot(df[,xcol], df[,ycol], xlab=xcol, ylab=ycol, pch=pch, cex=cex)
  ## then we need to set up a number of functions that
  ## respond to mouse_input
  ## eventEnv is set by calling getGraphicsEvent further down. (i.e. not defined in the functions) 
  ## we store the points in eventEnv. Not sure if this is a good idea, but it seems to work.
  
  devset <- function()
    if(dev.cur() != eventEnv$which) dev.set(eventEnv$which)
  
  cg.mouse.down <- function(buttons, x, y){
    new.point <- c(grconvertX(x, "ndc", "user"), grconvertY(y, "ndc", "user"))
    nr <- nrow(eventEnv$m.points)
    ## if nr == 2, then we can determine the handedness of the gate
    if(nr == 2)
      eventEnv$is.right <- is.right.of(eventEnv$m.points[1,], eventEnv$m.points[2,], new.point)
    
    point.is.right <- eventEnv$is.right  ## defaults to true (needed when nr <= 2)
    gate.is.closable <- TRUE                ## can the new point be connected to close the gate?
    closing.is.ok <- TRUE
    if(nr > 2){
      point.is.right <- is.right.of(eventEnv$m.points[(nr-1),], eventEnv$m.points[nr,], new.point)
      gate.is.closable <- (point.is.right == is.right.of(eventEnv$m.points[nr,], new.point, eventEnv$m.points[1,]))
      closing.is.ok <- (point.is.right == is.right.of(new.point, eventEnv$m.points[1,], eventEnv$m.points[2,]))
    }
    if((point.is.right == eventEnv$is.right) && gate.is.closable && closing.is.ok){
      eventEnv$m.points <- rbind(eventEnv$m.points, new.point)
      nr <- nrow(eventEnv$m.points)
      points(eventEnv$m.points[nr, 1], eventEnv$m.points[nr, 2], col="red")
      if(nr > 1)
        lines(eventEnv$m.points[(nr-1):nr,1], eventEnv$m.points[(nr-1):nr,2], col="red")
    }else{
      cat(paste("Illegal Point. Try again.\n"))
      lines(c(eventEnv$m.points[(nr-1),1], new.point[1]),
            c(eventEnv$m.points[(nr-1),2], new.point[2]), col="blue")
    }
    devset()  ## not sure what is the point of this..
    NULL  ## ???
  }
  
  keydown <- function(key){
    reset.points <- function(){
      eventEnv$m.points <- matrix(ncol=2, nrow=0)
      colnames(eventEnv$m.points) <- c(xcol, ycol)
    }
    add.gate <- function(){
      if(nrow(eventEnv$m.points) > 2){
        eventEnv$m.points <- rbind(eventEnv$m.points, eventEnv$m.points[1,])
        eventEnv$gates[[length(eventEnv$gates)+1]] <-
          list(points=eventEnv$m.points, right.handed=eventEnv$is.right, xcol=xcol, ycol=ycol, xlog=xlog, ylog=ylog)
      }
      plot.graph()
      plot.gates(eventEnv$gates)
    }
    if (key == "c"){
      reset.points()
      plot.graph()
      plot.gates(eventEnv$gates)
    }
    if(key == "a"){
      add.gate()
      reset.points()
    }
    if(key == "f"){
      add.gate()
      eventEnv$gate.is.accepted <- TRUE
      return(invisible(1))
    }
    if(key == "Up")
      cat("Up arrow\n")
    if(key == "Down")
      cat("Down Arrow\n")
    if(key == "Left"){
      cat("left Arrow\n")
      if(nrow(eventEnv$m.points) > 1){
        eventEnv$m.points <- eventEnv$m.points[1:(nrow(eventEnv$m.points)-1),]
        plot.graph()
        plot.gates(eventEnv$gates)
        plot.points()
      }
    }
    if(key == "Right")
      cat("right arrow\n")
    if (key == "q") return(invisible(1))
    eventEnv$onMouseMove <- NULL
    NULL
  }

  setGraphicsEventHandlers(prompt="click to add gate nodes",
                           onMouseDown = cg.mouse.down,
                           onKeybd = keydown)
  eventEnv <- getGraphicsEventEnv()
  eventEnv$m.points <- matrix(ncol=2, nrow=0)
  colnames(eventEnv$m.points) <- c(xcol, ycol)
  eventEnv$gates <- list()
  eventEnv$is.right <- TRUE
  eventEnv$gate.is.accepted <- FALSE
  ## and then rather dangerously:
  prompt.text <- paste(
                       "\tClick to create the points of the gates\n",
                       "\n\tThe gate is represented by one or more convex polygons.\n",
                       "\t(i.e. all turn must be either right handed or left handed)\n",
                       "\n",
                       "Press:\n\t\tc to clear the curret gate points\n\t\ta to add the gate\n",
                       "\t\tf to finally accept the gate selection\n",
                       "\t\tq to quit witout accepting the gate\n")
  ## setting prompt in getGraphicsEvent() seems to cause some sort of problem, so:
  cat(paste(prompt.text, "\n"))
  getGraphicsEvent()  ## this will hopefully return on pressing q
  ## return the points here
  gates <- list()
  if(eventEnv$gate.is.accepted) gates <- eventEnv$gates
#  for(i in 1:length(gates)){
#    if(xlog) gates[[i]]$points[,1] <- 10^gates[[i]]$points[,1]
#    if(ylog) gates[[i]]$points[,2] <- 10^gates[[i]]$points[,2]
#  }
  gates
}
