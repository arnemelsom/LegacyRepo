
# Function  allLinesPlain  returns a lists of index vectors
#  for continuous edge segments

# Input:
#  xp[1:np],yp[1:np]  : x, y grid coordinates of edge cells

# Output:
#  list of continous ice edge grid cell indices  lineList[[]] ,
#   first edge line is
#    (xp[lineList[[1]][1]],yp[lineList[[1]][1]]),
#    (xp[lineList[[1]][2]],yp[lineList[[1]][2]]),
#                ...,
#    (xp[lineList[[1]][L]],yp[lineList[[1]][L]])
#  where L is the length of the returned index vector  lineList[[1]]


# NOTE!!
#  This is  ** NOT **  a perfect algorithm, as it may split branching
#   segments in a sub-optimal way.


allLinesPlain <- function(xp,yp) {

    np <- length(xp)
    if (np == 0) return(NA)

    # lineList  will be the results returned from this function,
    #            with each list entry containing a vector of indices
    #            indices for input grid nodes xp[], yp[] that belongs
    #            to a continuous line segment
    lineList <- list()
    listNo   <- 0

    preList  <- list()  # preList will be used for initial storeage;
    preNo    <- 0       #  checked & re-shuffled before returned as  lineList
    postList <- list()  # for discarded segment appendices
    postNo   <- 0


    # Set (xp0,yp0) to the next (here: first) position for
    #  a continuous edge search:
    xp0 <- xp[1]
    yp0 <- yp[1]


    #  pList  will hold indices of all edge node no.s that have
    #          been assigned to a line segment
    pList    <- vector(mode="numeric")
    pList[1] <- 1

    # l   will be the length of the present edge segment
    # l0  will be the length of the present edge segment
    #      in the alternative direction
    l  <- 1  # initialize (need to re-initialize for each segment)
    l0 <- 0  # initialize (need to re-initialize for each segment)

    empty <- vector(mode="numeric")
    pDone <- vector(mode="numeric")
    n        <- length(pDone)
    pAll     <- 1:np
    pRest    <- setdiff(pAll,pDone)  # list of non-assigned
                                     #  edge grid node indices

    #  n  is no of edge nodes that have been assigned to a line
    #      segment. Keep adding new segments as long as there are
    #      remaining (non-assigned) edge nodes
    while (n < np) {

        # Set  pRest  to the set of cell indices yet to be examined
        pRest <- setdiff(pRest,pList[1])  # remove pList[1] from pRest
        if (length(pRest) > 0) {

            # Find neighbour node(s),
            #  defined as at a distance of 1 or sqrt(2) (diagonal)
            ip <- which(abs(xp[pRest] - xp0) < 2 &
                        abs(yp[pRest] - yp0) < 2)
            while (length(ip) > 0) {
                # ...found a neighbour
                l        <- l + 1
                pList[l] <- pRest[ip[1]]  # include pRest[ip[1]] in line pList
                pRest    <- setdiff(pRest,pList[l]) # remove pList[l] from pRest
                if (length(pRest) == 0) ip <- pRest else {
                    # Look for another neighbour not in  pList
                    xpr      <- xp[pList[l]]
                    ypr      <- yp[pList[l]]
                    ip       <- which(abs(xp[pRest] - xpr) < 2 &
                                      abs(yp[pRest] - ypr) < 2)
                }
            }

            # ... pList  is now a continuos ice edge segment, but will be
            #             lacking closure at (xp0,yp0) if this is not the
            #             end node in the line, but in the middle
            # ...so check for another neighbour at (xp0,yp0)
            if (length(pRest) > 0) {
                ip <- which(abs(xp[pRest] - xp0) < 2 &
                            abs(yp[pRest] - yp0) < 2)
                while (length(ip) > 0) {
                    l0          <- l0 + 1   # as  l  but other direction
                    pList[l+l0] <- pRest[ip[1]]
                    pRest       <- setdiff(pRest,pList[l+l0])
                    if (length(pRest) == 0) ip <- pRest else {
                        xpr         <- xp[pList[l+l0]]
                        ypr         <- yp[pList[l+l0]]
                        ip          <- which(abs(xp[pRest] - xpr) < 2 &
                                             abs(yp[pRest] - ypr) < 2)
                    }
                }
            }      #  <-  if (length(pRest) > 0)
        }          #  <-  if (length(pRest) > 0)

        # Update index list for continuous ice edge:
        preNo <- preNo + 1
        preList[[preNo]] <- 0
        # ...if  l0 > 0  shuffle latest  l0  nodes
        #                 to make  preList  continous
        if (l0 > 0)
            # ...swap by mapping  (l+l0):(l+1)  to  1:l0
            for (i in 1:l0) preList[[preNo]][i] <- pList[l+l0-i+1]
        preList[[preNo]][(l0+1):(l0+l)] <- pList[1:l]

        # Update list of nodes that have been checked:
        pDone <- c(pDone,pList[1:(l+l0)])
        # Initialize next search, if relevant:
        n <- length(pDone)
        if (n < np) {
            pRest    <- setdiff(pAll,pDone)
            ip       <- pRest[1]
            xp0      <- xp[ip]
            yp0      <- yp[ip]
            pList[1] <- ip
            l        <- 1   # reinitialized for next segment
            l0       <- 0   # reinitialized for next segment
        }
    }

    # What we wish to return is now  preList  but for
    #  convenience we shuffle the line segments so that
    #  the longest segment is returned as list no.  [[1]]
    #  the second longest as list no.  [[2]]  etc.

    # Get list elements length (=no. grid nodes in segment)
    segLength <- vector(mode="numeric", length=preNo)
    for (p in 1:preNo) segLength[p] <- length(preList[[p]])

    # Re-arrange so that line segments are returned longest to shortest
    segSort <- sort(segLength,index.return=TRUE)
    n       <- preNo + 1
    for (i in segSort$ix) {
        n <- n - 1
        lineList[[n]] <- preList[[i]]
    }

    return(lineList)

}
