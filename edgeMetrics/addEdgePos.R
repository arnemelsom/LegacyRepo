
# PURPOSE
#  To add grid information that specifies an edge
#   along a specific isopleth ('SICedgeVal')
#  This is done for two separate variables ('varName1', 'varName2')
#  Note that dimension names are also hard coded ('x','y', 't')

# Variables that will be appended to the netCDF file are
#  npEdge1[t], npEdge2[t]
#      - no. edge nodes at each time for var.s 'varName1', 'varName2'
#  xpEdge1[t,npEdge1], ypEdge1[t,npEdge1]
#      - x,y grid positions of edge nodes for 'varName1'
#        (lower left corner is [1,1])
#  xpEdge2[t,npEdge2], ypEdge2[t,npEdge2]
#      - x,y grid positions of edge nodes for 'varName2'
# Note that trailing non-used entries of '[xy]pEdge[12]'
#  will have NA values


# --- User's specifications start here ---

# Set file name:
fileName <- "SICdemo.nc4"

# Set variable names:
varName1 <- "sic1"
varName2 <- "sic2"

# Set sea ice concentration value used to define ice edge:
isoPlethValue <- 0.15

# Set dimension names:
xd <- "x"
yd <- "y"
td <- "t"

# Name of global attribute w/ value for edge defining isopleth
isoPlethName <- "sea_ice_concentration_at_edge"

# --- User's specifications end here ---


library(RNetCDF)

source("edgeNodes.R")


id <- open.nc(fileName, write=TRUE)

nx <- dim.inq.nc(id,xd)$length
ny <- dim.inq.nc(id,yd)$length
nt <- dim.inq.nc(id,td)$length

start <- c( 1, 1,1)
count <- c(nx,ny,1)

edgeNodes1 <- list()
edgeNodes2 <- list()
np1        <- vector(mode="numeric", length=nt)
np2        <- vector(mode="numeric", length=nt)


# Detect ice edge nodes:
for (t in 1:nt) {
    start[3]   <- t

    # Read, handle 'varName1'
    SIC1            <- var.get.nc(id,varName1,start=start,count=count)
    edgeNodes1[[t]] <- setEdgeNodes(SIC1,isoPlethValue)
    np1        [t]  <- edgeNodes1[[t]]$nEdge

    # Read, handle 'varName2'
    SIC2            <- var.get.nc(id,varName2,start=start,count=count)
    edgeNodes2[[t]] <- setEdgeNodes(SIC2,isoPlethValue)
    np2        [t]  <- edgeNodes2[[t]]$nEdge

}

# Set max no. of edge grids for each of the two variables:
nMax1 <- max(np1)
nMax2 <- max(np2)


# Create variables to store results, add some attributes:

att.put.nc(id,"NC_GLOBAL",isoPlethName,"NC_FLOAT",isoPlethValue)

dim.def.nc(id,"npEdge1",max(nMax1,1))
var.def.nc(id,"npEdge1","NC_INT", td)
att.put.nc(id,"npEdge1","long_name","NC_CHAR","no. edge grids")
var.def.nc(id,"xpEdge1","NC_INT", c("npEdge1",td))
att.put.nc(id,"xpEdge1","long_name", "NC_CHAR","edge grid x-position")
att.put.nc(id,"xpEdge1","_FillValue","NC_INT", -1e6)
var.def.nc(id,"ypEdge1","NC_INT", c("npEdge1",td))
att.put.nc(id,"ypEdge1","long_name","NC_CHAR","edge grid y-position")
att.put.nc(id,"ypEdge1","_FillValue","NC_INT", -1e6)

dim.def.nc(id,"npEdge2",max(nMax2,1))
var.def.nc(id,"npEdge2","NC_INT", td)
att.put.nc(id,"npEdge2","long_name","NC_CHAR","no. edge grids")
var.def.nc(id,"xpEdge2","NC_INT", c("npEdge2",td))
att.put.nc(id,"xpEdge2","long_name", "NC_CHAR","edge grid x-position")
att.put.nc(id,"xpEdge2","_FillValue","NC_INT", -1e6)
var.def.nc(id,"ypEdge2","NC_INT", c("npEdge2",td))
att.put.nc(id,"ypEdge2","long_name","NC_CHAR","edge grid y-position")
att.put.nc(id,"ypEdge2","_FillValue","NC_INT", -1e6)


# Store edge nodes:

var.put.nc(id,"npEdge1",np1)
var.put.nc(id,"npEdge2",np2)

count1 <- c(nMax1,1)
count2 <- c(nMax2,1)
for (t in 1:2) {
    start <- c(1,t)

    # First, store edge node grid coordinates in 'xp', 'yp',
    #  fill non-used (trailing) entries with  NA ,
    #  and write to file

    # edge no. 1:
    n <- np1[t]
    if (n > 0) {
        xp <- edgeNodes1[[t]]$xp[1:n]
        yp <- edgeNodes1[[t]]$yp[1:n]
    }
    if (n < nMax1) {
        xp[(n+1):nMax1] <- NA
        yp[(n+1):nMax1] <- NA
    }
    var.put.nc(id,"xpEdge1",xp,start=start,count=count1)
    var.put.nc(id,"ypEdge1",yp,start=start,count=count1)
    rm(xp)
    rm(yp)

    # edge no. 2:
    n <- np2[t]
    if (n > 0) {
        xp <- edgeNodes2[[t]]$xp[1:n]
        yp <- edgeNodes2[[t]]$yp[1:n]
    }
    if (n < nMax2) {
        xp[(n+1):nMax2] <- NA
        yp[(n+1):nMax2] <- NA
    }
    var.put.nc(id,"xpEdge2",xp,start=start,count=count2)
    var.put.nc(id,"ypEdge2",yp,start=start,count=count2)
    rm(xp)
    rm(yp)

}

close.nc(id)

warnings()

quit("no")
