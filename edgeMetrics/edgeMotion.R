

# PURPOSE
#  To output results for distances between edges, including
#   information relevant for 'to' edge cells that are separated by
#   the decorrelation length scale.


# Sample sketch:

#  +:  'to'  edge (T is a sample considered edge node,
#                  t are other nodes considered for
#                   'Dfile'/'xFile'/'Yfile' if decorr. scale is 4 nodes)
#  *: 'from' edge (tilted line corresponds to shortest distance)

#      t+
#        ++T+++t+   +++t  'to' line
#           \    ++t
#            \
#             \
#              ***
#             *   **
#          ***      ******   'from' line
#
#  Here, the distance will be sqrt(3**2 + 4**2) = 5
#   (= -5  if sea ice has retreated /
#          high concentration values on the upper side)


# Output line contents on 'outFile':
#  col  1:  date (floating point year) -NOT USED HERE- (output is set to 1)
#  col  2:  no edge grids in "from" ice edge
#  col  3:  no edge grids in  "to"  ice edge
#  col  4:  "from node", x pos. / largest offset
#  col  5:  "from node", y pos.
#  col  6:  "to node", x pos.
#  col  7:  "to node", y pos.
#  col  8:  separation distance (unit as dx/dy, e.g. non-dim or km)
#  col  9:  spatial decorrelation length of distances
# ..units of distance in col 8 will correpond to spec. for dx, dy below


# Output contents on 'Dfile':
#  vector of displacements from 'to' grid nodes to 'from' edge line,
#   with 'to' nodes separated by the decorrelation length
# Output contents on 'xFile', 'yFile':
#  grid node positions for 'to' nodes for which displacements
#   are stored in 'Dfile'


# --- User's specifications start here ---

# Set input file name
fileName <- "SICdemo.nc4"
# ...edge information must have been added to  fileName  w/  addEdgePos.R
#     prior to executing this script (---not checked---)

# Set dimension names:
xd <- "x"
yd <- "y"

# Set partial variable names (must add '1'/'2', done in loop below):
npVarN <- "npEdge"
xpVarN <- "xpEdge"
ypVarN <- "ypEdge"
varN   <- "sic"

# Name of global attribute w/ value for edge defining isopleth
isoPlethName <- "sea_ice_concentration_at_edge"

# Autocorrelation (AC) value used to define decorrelation distance
#  (the distance at which AC drops below this value)
AClim      <- exp(-1)

# Autocorrelation will be computed for offsets from 1 to  maxGridOff
maxGridOff <- 100

# Line segments that are composed of fewer than  minLen  grid nodes
#  will be discarded from the analysis below:
minLen <- 20

# --- User's specifications end here ---



library(RNetCDF)

source("IEdistance.R")
source("allLinesPlain.R")

noList   <- c(1,2)
typeList <- c("ref","mod")


# This version uses the non-dimensional units of grid cell no.
dx <- 1.
dy <- 1.
dA <- dx*dy
ds <- sqrt(dA)                  # normalized, retaining dA


# Open file with concentrations and specifications of edge lines:
id <- open.nc(fileName)

# Read domain size:
nx  <- dim.inq.nc(id, xd)$length
ny  <- dim.inq.nc(id, yd)$length

# Strings defining format of output vectors:
fmtRef  <- "%3d %3d %3d %3d %5.1f %3d"
fmtStr  <- "%10.4f %4d %4d"


# Loop products as given by  typeList  entries above

for (ni in 1:length(noList)) {
    no <- noList[ni]

    # Set names of output ascii files
    type    <- typeList[ni]
    outFile <- paste("edgeMotionStats_",minLen,"Lgrd_",
                                        type,".dat", sep="")
    xFile   <- paste(           "xpos_",type,".dat", sep="")
    yFile   <- paste(           "ypos_",type,".dat", sep="")
    Dfile   <- paste(       "distance_",type,".dat", sep="")

    # Set full names of variables on netCDF file  'fileName'
    vNm <- paste(npVarN,no, sep="")
    xNm <- paste(xpVarN,no, sep="")
    yNm <- paste(ypVarN,no, sep="")
    sNm <- paste(  varN,no, sep="")

    # Read edge grid nodes
    # ...no. edge grid nodes
    np   <- var.get.nc(id,vNm)
    # ...x,y coordinates of edge grid nodes, 'from' edge
    fnp  <- np[1]
    fxp  <- var.get.nc(id,xNm,start=c(1,1),count=c(fnp,1))
    fyp  <- var.get.nc(id,yNm,start=c(1,1),count=c(fnp,1))
    # ...x,y coordinates of edge grid nodes,  'to'  edge
    tnp  <- np[2]
    txpa <- var.get.nc(id,xNm,start=c(1,2),count=c(tnp,1))
    typa <- var.get.nc(id,yNm,start=c(1,2),count=c(tnp,1))

    # ...read concentration values corresponding to 'from' edge:        
    start   <- c( 1, 1,1)
    count   <- c(nx,ny,1)
    edgeVal <- att.get.nc(id,"NC_GLOBAL",isoPlethName)
    sic     <- var.get.nc(id,sNm, start=start,count=count)


    # Prepare output:
    outLine <- vector(mode="character")
    lNo     <- 0
    
    
    # Re-distribute 'to' grid nodes in separated line segments:
    lineSeg <- allLinesPlain(txpa,typa)
    nLines  <- length(lineSeg)

    #  txp/typ  will exclude 'to' edge grid nodes from txpa/typa
    #             from line segments shorter than  minLen  nodes
    txp <- vector(mode="numeric")
    typ <- vector(mode="numeric")


    # *** PART 1  starts ***

    #  Determine the spatial decorrelation length scale in
    #   edge grid node length, to be stored as  cn

    dpAll <- vector(mode="numeric")
    tnp   <- 0   # dynamic length of  txp/typ
    wSum  <- 0   # weight sum for decorrelation length
    dLsum <- 0   # length-weighted sum of decorrelation length
    for (n in 1:nLines) {
        dp  <- vector(mode="numeric")
        tnn <- length(lineSeg[[n]])
        # ...only keep line segments longer than  minLen :
        if (tnn >= minLen) {
            tna          <- tnp + 1
            tnp          <- tnp + tnn
            indx         <- lineSeg[[n]]
            txp[tna:tnp] <- txpa[indx]
            typ[tna:tnp] <- typa[indx]
            ii           <- 0
            # Compute, store shortest distance from 'to' grid node
            #  to 'from' edge line:
            for (i in tna:tnp) {
                ii       <- ii + 1
                dd       <- sqrt((txp[i] - fxp)*(txp[i] - fxp) +
                                 (typ[i] - fyp)*(typ[i] - fyp))
                dp  [ii] <- min(dd)
                dpAll[i] <- dp[ii]
            }

            # Compute spatial auto correlation for this line segment:
            corrList   <- vector(mode="numeric", length=maxGridOff)
            for (offset in 1:min(maxGridOff,(tnn-2))) {
                le <- tnn - offset
                v1 <- dp[1:le]
                v2 <- dp[(1+offset):tnn]
                corrList[offset] <- cor(v1,v2) # de-corr. values
            }
            if (any(corrList <= AClim)) {
                c0 <- 1
                while (corrList[c0] > AClim) c0 <- c0 + 1
                # ...now, c0 is the auto decorrelation grid length
                wSum  <-  wSum + tnn
                dLsum <- dLsum + tnn*c0
            }

        }  #  <-  if (tnn >= minLen)
    }      #  <-  for (n in 1:nLines)

    # Compute weighted average spatial auto decorrelation  cn :
    if (wSum == 0) {
        print("Decorrelation length scale is too large.")
        print("Try increasing value of  maxGridOff  and re-run.")
        close.nc(id)
        quit("no")
    } else cn <- round(dLsum/wSum)

    # *** PART 1  ends ***


    # Compute stat.s:
    if (fnp >= 1 & tnp >= 1) {
        distro <- edgeDistro(dx,dy,fxp,fyp,sic,txp,typ,edgeVal,1)
        if (is.na(distro[1])) {
            close.nc(id)
            quit("no")
        }
        # Store results for output; see top of script for details
        #  (and also specification for list of output for function
        #    distro  on a separate file)
        lNo <- lNo + 1
        outLine[lNo] <- sprintf(fmtStr,1.,fnp,tnp)
        outLine[lNo] <- paste(outLine[lNo],
                              sprintf(fmtRef,distro$x0[1],distro$y0[1],
                                             distro$x[1], distro$y[1],
                                             distro$d[1],cn))
    }

    # Create strings for distances and grid cell coordinates,
    #  to be stored for the selection of 'to' grid cells that are
    #  separated by the decorrelation distance  cn
    Dline <- vector(mode="character")
    xLine <- vector(mode="character")
    yLine <- vector(mode="character")
    dn    <- 0      # dynamic index of Dline/xLine/yLine

    # Make sure that the 'to' grid node with the largest
    #  displacement is included in the output:
    d1 <- which(txp == distro$x[1] & typ == distro$y[1])
    if (length(d1) < 1) quit("no")  # This should never happen!
    d0 <- d1[1]
    # ...move back  cn  nodes at a time to the "beginning"
    #     (we wish to include  ..., d1[1]-cn, d1[1], d1[1]+cn, ...)
    while (d0-cn > 0) d0 <- d0 - cn
    dn        <- dn + 1
    xLine[dn] <- ""
    yLine[dn] <- ""
    Dline[dn] <- ""
    while (d0 <= tnp) {
        xLine[dn] <- paste(xLine[dn],txp  [d0])
        yLine[dn] <- paste(yLine[dn],typ  [d0])
        Dline[dn] <- paste(Dline[dn],dpAll[d0])
        d0        <- d0 + cn
    }


    # Write results to files
    cat(xLine,   file=xFile,   sep='\n')
    cat(yLine,   file=yFile,   sep='\n')
    cat(Dline,   file=Dfile,   sep='\n')
    cat(outLine, file=outFile, sep='\n')

}   #  <-  for (ni in 1:length(noList)


close.nc(id)

warnings()

quit("no")
