
source("R/node.R")

BeamNode <- setClass("BeamNode", representation(width="numeric", depth="numeric"), contains = "Node")

setGeneric("beamSearch", function(object) { standardGeneric("beamSearch") } )

setMethod("beamSearch", signature("BeamNode"), function(beamNode) {

  depth <- beamNode@depth
  width <- beamNode@width

  # Generate child nodes

  # Find top n child nodes

  # Recursively search

  sapply(top.child.nodes, beamSearch)

})

# PrimNode <- setClass("PrimNode", representation(quality="numeric"), prototype(width = 1), contains = "BeamNode")

