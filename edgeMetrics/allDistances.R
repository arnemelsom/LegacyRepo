
# This script computes shortest distances for two sets of two edges,
#  results are considered separtely for each of the sets, and written
#  to separate output files.

# Time step one has the 'from' edge line/nodes,
#      step two has the   'to' edge line/nodes
# We seek to find the shortest distance between the 'from' edge line
#  and each of the 'to' edge nodes. This is a quantification of the
#  displacement of the ice from one time to another in order to reach
#  the 'to' edge node.


# Sample sketch:

#  +:   'to' edge (T is a sample considered edge node)
#  *: 'from' edge (tilted line corresponds to shortest distance)


#      ++
#        ++T+++++   ++++  'to' line
#           \    +++
#            \
#             \
#              ***
#             *   **
#          ***      ******  'from' line
#
#  Here, the distance will be sqrt(3**2 + 4**2) = 5
#   (= -5  if sea ice has retreated /
#          high concentration values on the upper side)


# --- User's specifications start here ---

# Set input file name
fileName <- "SICdemo.nc4"
# ...edge information must have been added to  fileName  w/  addEdgePos.R
#     prior to executing this script (---not checked---)


# Set partial variable names (must add '1'/'2'):
npVarN <- "npEdge"
xpVarN <- "xpEdge"
ypVarN <- "ypEdge"
varN   <- "sic"

# Set dimension names:
xd <- "x"
yd <- "y"
td <- "t"

# Name of global attribute w/ value for edge defining isopleth
isoPlethName <- "sea_ice_concentration_at_edge"

# --- User's specifications end here ---


library(RNetCDF)


# This version uses the non-dimensional units of grid cell no.
dx   <- 1
dy   <- 1
dA   <- dx*dy
ds   <- sqrt(dA)                  # normalized, retaining dA
dxSq <- dx*dx
dySq <- dy*dy


# Open input file, set dimensions of 2d field:
id  <- open.nc(fileName)
nx  <- dim.inq.nc(id, xd)$length
ny  <- dim.inq.nc(id, yd)$length

for (edgeNo in 1:2) {

    # Name of output file with results for distances between
    DallFile <- paste("distance_all_edge",edgeNo,".dat", sep="")

    # Set names of variables for this edge no.:
    npVar  <- paste(npVarN, edgeNo, sep="")
    xpVar  <- paste(xpVarN, edgeNo, sep="")
    ypVar  <- paste(ypVarN, edgeNo, sep="")
    sicVar <- paste(varN,   edgeNo, sep="")

    # Read field values for the 'from' time step (no. 1)
    #  and the ispoleth value that defines the edge
    sic <- var.get.nc(id,sicVar, start=c( 1, 1,1),
                                 count=c(nx,ny,1))
    edgeVal <- att.get.nc(id,"NC_GLOBAL",isoPlethName)

    # Read no. edge grid cells (two time steps):
    np  <- var.get.nc(id,npVar)
    fnp <- np[1]
    tnp <- np[2]

    # Define vector that will hold all distances,
    #  i.e., one distance for each 'to' edge node
    dAll <- vector(mode="numeric", length=tnp)

    if (tnp > 0 & fnp > 0) {

        # Read grid nodes for 'from', 'to' line:
        fxp <- var.get.nc(id,xpVar,start=c(1,1),count=c(fnp,1))
        fyp <- var.get.nc(id,ypVar,start=c(1,1),count=c(fnp,1))
        txp <- var.get.nc(id,xpVar,start=c(1,2),count=c(tnp,1))
        typ <- var.get.nc(id,ypVar,start=c(1,2),count=c(tnp,1))

        # Find distance to 'from' line for each of the  tnp  'to' grids:
        for (n in 1:tnp) {
            dp <- (fxp - txp[n])*(fxp - txp[n])*dxSq +
                  (fyp - typ[n])*(fyp - typ[n])*dySq
            dd <- sqrt(dp)
            ip <- which(dd == min(dd))
            # Store values, positive for expanding cover,
            #               negative for retreating
            # (recall that  sic  is for the 'from' time step)
            if (sic[txp[n],typ[n]] <= edgeVal) dAll[n] <-  dd[ip[1]] else
                                               dAll[n] <- -dd[ip[1]]
        }

    }      #  <-  if (tnp > 0 & fnp > 0)

    # Write results to file:
    cat(dAll, file=DallFile)

}          #  <-  for (edgeNo in 1:2)

# To un-scramble irregular data on x[ft]File, y[ft]File:
#  xData <- scan(x[ft]File,what="character", sep="\n")
#  xNode <- as.numeric(unlist(strsplit(xData[1], " "))) etc.


close.nc(id)

warnings()

quit("no")
