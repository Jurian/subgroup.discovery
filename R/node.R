
Node <- setClass("Node", representation(parent="Node"), prototype(parent=NULL))

setGeneric("isRoot", function(object) { standardGeneric("isRoot")} )

setMethod("isRoot", signature("Node"), function(object) {is.null(object@parent)})
