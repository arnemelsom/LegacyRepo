

# Function  edgeDistro  returns edge distances from a set of grid nodes

# Elsewhere,
#  grid nodes in the primary   edge is referred to as 'from' nodes
#  grid nodes in the secondary edge is referred to as  'to'  nodes

#  Note: requires vectors with grid positions (xp1[] etc.) to be non-empty!
#         (not checked here)

# Input:
#  dx, dy             : grid resolution
#  (xp1[],yp1[])      : edge node grid cell coordinates, primary   edge [int]
#  (xp2[],yp2[])      : edge node grid cell coordinates, secondary edge [int]
#  SIC1[1:nx,1:ny]    : concentration (0-1), primary field
#  edgeVal            : concentration isopleth value that defines the edge
#  gDistro[]          : selected, ordered grid distance
#                        (e.g. =2: return the second largest distance from all edge
#                                  nodes in primary product to secondary edge;
#                                  examination restricted to positive values, i.e.,
#                                  when secondary product is on open ocean side)
#                       NOTE! must be monotonic (lowest to highest)

# Output:
#  For each ordered distance entry in  gDistro :
#  $x,  $y     :  x,y position of secondary edge grid node
#  $x0, $y0    :  x,y position of primary edge grid node
#  $d          :  distance between node $x,$y (on secondary edge)
#                  and primary edge (at $x0,$y0)
# If there are no positive distances, all values will be set to -9


#  NOTE: This version requires that boundary edge nodes have been removed!



edgeDistro <- function(dx,dy,xp1,yp1,SIC1,xp2,yp2,edgeVal,gDistro) {

    edgeDistroOut <- list()

    n1  <- length(xp1)
    n2  <- length(xp2)
    if (n2 == 0 | n1 == 0) {
        print("Error in function edgeDistro. One or both edges are empty.")
        print(paste("Edge 1 length:",n1,", edge 2 length",n2))
        return(NA)
    }
    nxy <- dim(SIC1)
    nx  <- nxy[1]
    ny  <- nxy[2]
    ng  <- length(gDistro)
    if (ng > 1)
        # ...check monotonicity
        if (any(gDistro[2:ng] - gDistro[1:(ng-1)] <= 0)) return(NA)

    # Adjustment factors if dx != dy
    dxSq     <- dx*dx
    dySq     <- dy*dy

    # For each edge node in secondary edge, determine
    #  dn      - distance to primary edge
    #  xn, yn  - grid node x,y position in primary edge
    #              to which  dn  is measured
    dn  <- vector(mode="numeric", length=n2)
    xn  <- vector(mode="numeric", length=n2)
    yn  <- vector(mode="numeric", length=n2)
    for (n in 1:n2) {
        xp    <- xp2[n]
        yp    <- yp2[n]
        dSq   <- (xp1 - xp)*(xp1 - xp)*dxSq +
                 (yp1 - yp)*(yp1 - yp)*dySq
        dd    <- sqrt(dSq)
        dm    <- min(dd)
        nn    <- which(dd == dm)
        xn[n] <- xp1[nn[1]]
        yn[n] <- yp1[nn[1]]
        # If (primary) concentration  SIC1  is below isopleth value
        #  at secondary edge node, the sea ice has expanded,
        #  store distance as positive:
        if (SIC1[xp,yp] <= edgeVal) dn[n] <-  dm else
                                    dn[n] <- -dm
        if (any(is.na(dd))) quit("no")  # This should never happen.
    }

    # Sort results for distances:
    di    <- sort.int(dn,index.return=TRUE)
    # ...consider non-negative displacements only
    ipos  <- which(di$x >= 0)
    ii    <- length(ipos)

    # Define vectors to hold output values for each of the
    #  gDistro[1:ng] elements
    x0pos <- vector(mode="numeric", length=ng)
    y0pos <- vector(mode="numeric", length=ng)
    xpos  <- vector(mode="numeric", length=ng)
    ypos  <- vector(mode="numeric", length=ng)
    dist  <- vector(mode="numeric", length=ng)
    x0pos[] <- -9
    y0pos[] <- -9
    xpos [] <- -9
    ypos [] <- -9
    dist [] <- -9

    # For each  gDistro  element, overwrite default value  -9 
    #  with properly selected results:
    for (i in 1:ng)
        if (ii >= gDistro[i]) {
            ig       <- n2 - gDistro[i] + 1
            ii0      <- di$ix[ig]
            x0pos[i] <- xn [ii0]
            y0pos[i] <- yn [ii0]
            xpos [i] <- xp2[ii0]
            ypos [i] <- yp2[ii0]
            dist [i] <- dn [ii0]
        }

    edgeDistroOut$x0 <- x0pos
    edgeDistroOut$y0 <- y0pos
    edgeDistroOut$x  <- xpos
    edgeDistroOut$y  <- ypos
    edgeDistroOut$d  <- dist
    return(edgeDistroOut)

}



# edgeDistroCoast  is like  edgeDistro but includes coastline,
#                   boundary in search for closest node

edgeDistroCoast <- function(dx,dy,xp1,yp1,SIC1,xp2,yp2,edgeVal,gDistro) {

    edgeDistroOut <- list()

    n1  <- length(xp1)
    n2  <- length(xp2)
    if (n2 == 0 | n1 == 0) quit("no")
    nxy <- dim(SIC1)
    nx  <- nxy[1]
    ny  <- nxy[2]
    ng  <- length(gDistro)
    if (ng > 1)
        # ...check monotonicity
        if (any(gDistro[2:ng]-gDistro[1:(ng-1)] <= 0)) return(NA)

    # Find edge of dry nodes (xpNA[],ypNA[])
    xpNA <- vector(mode="numeric")
    ypNA <- vector(mode="numeric")
    nNA  <- 0
    indx <- which(is.na(SIC1))
    for (i in indx) {
        y  <- 1 + floor((i-1)/nx)
        x  <- i - (y - 1)*nx
        xa <- max(x-1, 1)
        xe <- min(x+1,nx)
        ya <- max(y-1, 1)
        ye <- min(y+1,ny)
        if (any(!is.na(c(SIC1[xa,y],SIC1[xe,y],
                         SIC1[x,ya],SIC1[x,ye])))) {
            nNA       <- nNA + 1
            xpNA[nNA] <- x
            ypNA[nNA] <- y
        }
    }

    # Consider distance between edges:
    dxSq     <- dx*dx
    dySq     <- dy*dy

    # ...distances to primary edge:
    xn  <- vector(mode="numeric", length=n2)
    yn  <- vector(mode="numeric", length=n2)
    dn  <- vector(mode="numeric", length=n2)
    if (nNA > 0) {
        xp1 <- c(xp1,xpNA)
        yp1 <- c(yp1,ypNA)
    }
    for (n in 1:n2) {
        xp    <- xp2[n]
        yp    <- yp2[n]
        dSq   <- (xp1 - xp)*(xp1 - xp)*dxSq +
                 (yp1 - yp)*(yp1 - yp)*dySq
        dd    <- sqrt(dSq)
        dm    <- min(dd)
        nn    <- which(dd == dm)
        xn[n] <- xp1[nn[1]]
        yn[n] <- yp1[nn[1]]
        if (SIC1[xp,yp] <= edgeVal) dn[n] <-  dm else
                                    dn[n] <- -dm
        if (any(is.na(dd))) quit("no")
    }

    # Sort results:
    di    <- sort.int(dn,index.return=TRUE)
    # ...consider non-negative displacements only
    ipos  <- which(di$x >= 0)
    ii    <- length(ipos)
    x0pos <- vector(mode="numeric", length=ng)
    y0pos <- vector(mode="numeric", length=ng)
    xpos  <- vector(mode="numeric", length=ng)
    ypos  <- vector(mode="numeric", length=ng)
    dist  <- vector(mode="numeric", length=ng)
    x0pos[] <- -9
    y0pos[] <- -9
    xpos [] <- -9
    ypos [] <- -9
    dist [] <- -9
    for (i in 1:ng)
        if (ii >= gDistro[i]) {
            ig       <- n2 - gDistro[i] + 1
            ii0      <- di$ix[ig]
            x0pos[i] <- xn [ii0]
            y0pos[i] <- yn [ii0]
            xpos [i] <- xp2[ii0]
            ypos [i] <- yp2[ii0]
            dist [i] <- dn [ii0]
        }

    edgeDistroOut$x0 <- x0pos
    edgeDistroOut$y0 <- y0pos
    edgeDistroOut$x  <- xpos
    edgeDistroOut$y  <- ypos
    edgeDistroOut$d  <- dist
    return(edgeDistroOut)

}

