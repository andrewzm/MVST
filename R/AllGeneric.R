#### GENERIC FUNCTIONS ######
#' @title Subset an MVST object
#' @description This function can be used to subset observations, meshes, and latent fields in the usual manner of \code{subset}.
#' @param x the MVST object
#' @param ... relations pertaining to the data frame within the MVST object
#' @return object of the same class as \code{x} unless \code{x} is a mesh, in which case a data frame is returned.
#' @export
#' @examples 
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' icesat_sub <- subset(icesat_obs,t==3)
#' print(class(icesat_sub))
#' 
#' data(surf_fe)
#' Mesh <- initFEbasis(p=surf_fe$p, t = surf_fe$t, M = surf_fe$M, K = surf_fe$K)
#' Mesh_sub <- subset(Mesh,x < 0) 
#' print(class(Mesh_sub))
setGeneric("subset")


#' @title Plot meshes and observations
#' @description This function is the general one used for plotting/visualising objects with MVST.
#' @param x an MVST object of class \code{Obs} or \code{Basis_GMRF}
#' @param y a character indicated which column to plot
#' @param ... other parameter used to configure the plot. These include \code{max} (upperbound on colour scale) and \code{min} (lowerbound on colour scale)
#' @return a ggplot2 object
#' @export
#' @examples 
#' \dontrun{
#' data(surf_fe)
#' Mesh <- initFEbasis(p=surf_fe$p, t = surf_fe$t, M = surf_fe$M, K = surf_fe$K)
#' Mesh["z"] <- sin(Mesh["x"]/1000)*cos(Mesh["y"]/1000)
#' g <- plot(Mesh,"z")
#' 
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' plot(subset(icesat_obs,t==2),"z",min=-0.5,max=0.5)}
setGeneric("plot")

#' @title Plots an interpolated field from a mesh
#' @docType methods
#' @description This function takes a Mesh and a column name in the mesh, to generate an interpolated
#' field of the indicated column
#' @param x a Mesh object
#' @param y a character indicated which column to plot
#' @param ds the resolution of the plot (defaulted to 40 x 40)
#' @param max upperbound on colour scale
#' @param min lowerbound on colour scale
#' @return a ggplot2 object
#' @export
#' @examples 
#' \dontrun{
#' data(surf_fe)
#' Mesh <- initFEbasis(p=surf_fe$p, t = surf_fe$t, M = surf_fe$M, K = surf_fe$K)
#' Mesh["z"] <- sin(Mesh["x"]/1000)*cos(Mesh["y"]/1000)
#' g <- plot_interp(Mesh,"z",ds=200)}
setGeneric("plot_interp", function(x,y,ds,...) standardGeneric("plot_interp"))

#' @title Get data frame
#' @description Returns data frame from MVST block object. Can be used on objects of class \code{FEBasis,GMRF} and \code{Obs}. If used on an object of class \code{Graph}, the observation data frame is returned.
#' @param .Object object from which data frame will be extracted
#' @export
#' @examples 
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' df <- getDf(icesat_obs)
setGeneric("getDf", function(.Object) standardGeneric("getDf"))


#' @title Get precision matrix
#' @description Returns precition matrix from MVST object. 
#' @param .Object object from which precision matrix will be extracted
#' @export
#' @examples 
#' my_RW <- GMRF_RW(n=10, order=1, precinc =2, name="my_first_RW")
#' Q <- getPrecision(my_RW)
setGeneric("getPrecision", function(.Object) standardGeneric("getPrecision"))

#' @title Get mean vector
#' @description Returns the mean from MVST object. 
#' @param .Object object from which mean vector will be extracted
#' @export
#' @examples 
#' my_RW <- GMRF_RW(n=10, order=1, precinc =2, name="my_first_RW")
#' mu <- getMean(my_RW)
setGeneric("getMean", function(.Object) standardGeneric("getMean"))

#' @title Get mass matrix
#' @description returns the mass matrix of a finite element basis
#' @param B an object of class \code{FEBasis}
#' @export
#' @examples
#' data(surf_fe)
#' Mesh <- initFEbasis(p=surf_fe$p, t = surf_fe$t, M = surf_fe$M, K = surf_fe$K)
#' M <- mass_matrix(Mesh)
setGeneric("mass_matrix", function(B) standardGeneric("mass_matrix"))

#' @title Get stiffness matrix
#' @description returns the stiffness matrix of a finite element basis
#' @param B an object of class \code{FEBasis}
#' @export
#' @examples
#' data(surf_fe)
#' Mesh <- initFEbasis(p=surf_fe$p, t = surf_fe$t, M = surf_fe$M, K = surf_fe$K)
#' M <- stiffness_matrix(Mesh)
setGeneric("stiffness_matrix", function(B) standardGeneric("stiffness_matrix"))


#' @title Return size of object Get stiffness matrix
#' @description Returns the size of the object, where the 'size' is determined by the class of the object being passed. For a \code{GMRF} and \code{GMRF_basis}, this is the number of variables of the associated MRF, for an \code{FEBasis} this is the number of vertices in the mesh, and for an \code{Obs} this is the number of observations.
#' @param x an MVST object
#' @export
#' @examples
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' nrow(icesat_obs)
setGeneric("nrow", function(x) standardGeneric("nrow"))

#' @title Return head object data frame
#' @description Short for \code{head(getDf())}. Returns the top of the object data frame which contains most information on the object, see \code{getDf} for details.
#' @param x an MVST object
#' @param ... other parameters passed to \code{base::head}
#' @export
#' @examples
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' head(icesat_obs,n=10)
setGeneric("head", function(x,...) standardGeneric("head"))

#' @title Return tail object data frame
#' @description Short for \code{tail(getDf())}. Returns the tail of the object data frame which contains most information on the object, see \code{getDf} for details.
#' @param x an MVST object
#' @param ... other parameters passed to \code{base::tail}
#' @export
#' @examples
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' tail(icesat_obs,n=10)
setGeneric("tail", function(x,...) standardGeneric("tail"))

#' @title Concatenation
#' @description Concatenates MVST objects of the same class together. This is primarily used to join up \code{GMRF_basis} blocks and \code{Obs} blocks together.
#' @param ... a series of \code{MVST} objects
#' @export
#' @examples
#' data(icesat)
#' icesat_obs <- Obs(df=icesat)
#' icesat_2x <- concat(icesat_obs,icesat_obs)
setGeneric("concat", function(...) standardGeneric("concat"))




setGeneric(".exist", function(L,to,from) standardGeneric(".exist"))
setGeneric("compress", function(Graph) standardGeneric("compress"))
setGeneric("getData", function(.Object) standardGeneric("getData"))
setGeneric("getC",function(.Object) standardGeneric("getC"))
setGeneric("setGMRF", function(Graph,obj) standardGeneric("setGMRF"))
setGeneric("setalpha", function(.Object,alpha,av_dist) standardGeneric("setalpha"))
setGeneric("Infer", function(Graph,...) standardGeneric("Infer"))

setGeneric("validate", function(Results,G,sim_obs=F,...) standardGeneric("validate"))
setGeneric("pred_variance_large", function(Results,G) standardGeneric("pred_variance_large"))
setGeneric("split_validation", function(.Object,samples,common,...) standardGeneric("split_validation"))
setGeneric("sample_GMRF", function(G,L=NULL,reps=1,P=NULL) standardGeneric("sample_GMRF"))


setGeneric("basisinterp", function(G,s,weights) standardGeneric("basisinterp"))
setGeneric(".extractClass", function(L,Cl) standardGeneric(".extractClass"))
setGeneric(".find_inc_matrix", function(basis,obs,...) standardGeneric(".find_inc_matrix"))

