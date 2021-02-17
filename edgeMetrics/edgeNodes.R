
# Function for detecting nodes along a given ispleth

#  Input:
#    val2d     input array (2D)
#    edgeVal   edge detection (isopleth) value

#  Edge nodes are nodes  with  val2d >= edgeVal  and at least
#   one x or y neighbor  with  val2d <  edgeVal

#  Output:
#   $nEdge  - no. of edge nodes detected
#   $xpEdge - x grid position of edge nodes (length $nEdge)
#   $ypEdge - y grid position of edge nodes (length $nEdge)
#  Note that lower left corner is [1,1]

setEdgeNodes <- function(val2d,edgeVal) {

    setEdgeNodesOut <- list()

    nxy <- dim(val2d)
    nx  <- nxy[1]
    ny  <- nxy[2]

    xpEdge  <- vector(mode="numeric")
    ypEdge  <- vector(mode="numeric")
    nEdge   <- 0

    # Search along constant y lines
    #  Note that 'ym', 'yp', 'xm', 'xp'
    #  are defined so that boundary issues are avoided
    for (y in 1:ny) {
        ym   <- max(c( 1,y-1))
        yp   <- min(c(ny,y+1))
        indx <- which(val2d[1:nx,y] >= edgeVal)
        # ...nodes [indx,y] are potential edge nodes;
        #    check if neighbor(s) have sub-isopleth value:
        if (length(indx) > 0)
            for (x in indx) {
                xm  <- max(c( 1,x-1))
                xp  <- min(c(nx,x+1))
                # ...set 'val' to up-down-left-right neighbors
                val <- c(val2d[xm,y],val2d[x,ym],
                         val2d[xp,y],val2d[x,yp])
                if (any(!is.na(val)))
                    # ... check for sub-isopleth value(s):
                    if (min(val, na.rm=TRUE) < edgeVal) {
                        nEdge         <- nEdge + 1
                        xpEdge[nEdge] <- x
                        ypEdge[nEdge] <- y
                    }
            }
    }

    setEdgeNodesOut$nEdge  <- nEdge
    setEdgeNodesOut$xpEdge <- xpEdge
    setEdgeNodesOut$ypEdge <- ypEdge
    return(setEdgeNodesOut)

}
