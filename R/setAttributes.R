#' Polygon-attribution
#' 
#' Takes a data frame \code{df} with fields \code{x} and \code{y} and a shape in data frame format, \code{shape_table}, and determines whether, and in which, polygon each point in \code{df} lies.
#' @param df data frame to which the polygon id is attributed. Must contain fields \code{x} and \code{y}.
#' @param shape_df the polygon-containing data frame. Must contain fields \code{x}, \code{y} and \code{id}. Multiple polygons/land-masses can be treated simultaneously through \code{id}.
#' @return a vector of length \code{nrow(df)} indicating the id in which each point in \code{df} lies.
#' @export
#' @examples
#' df <- data.frame(x=c(0.2,1),y=c(0.2,1))
#' shape_df <- data.frame(x=c(0,0,0.5,0.5),
#' y = c(0,0.5,0.5,0),
#' id = 1)
#' df$id <- attribute_polygon(df,shape_df)
attribute_polygon <- function(df,shape_df) {

  df$n <- 1:nrow(df)
  vals <- rep(0,nrow(df))
  
  if (length(unique(shape_df$id)) == 1) {
      df2 <- subset(df,x > min(shape_df$x) & x < max(shape_df$x) & y > min(shape_df$y) & y < max(shape_df$y)) # Find enclosing box
      myind <- df2[which(pnt.in.poly(cbind(df2$x,df2$y),shape_df[c("x","y")])$pip == 1),]$n
      vals[myind] <- T
  } else {
    for(i in unique(shape_df$id)){
      my_sub <- subset(shape_df,id==i)
      df2 <- subset(df,x > min(my_sub$x) & x < max(my_sub$x) & y > min(my_sub$y) & y < max(my_sub$y)) # Find enclosing box
      if(nrow(df2) > 0) {
        myind <- df2[which(pnt.in.poly(cbind(df2$x,df2$y),my_sub[c("x","y")])$pip == 1),]$n
        vals[myind] <- i
      }
    }
    
  }
  return(vals)
}

#' Data-attribution
#'
#' Attempts to attribute data in a data frame describing data on a grid to a another data frame (with irregular points in \code{x} and \code{y}. The data frame is assumed to have coordinate values as integer values. Hence, If the data-containing data frame describes a spacing of 1, a simple merge is carried out. If not a nearest-grid interpolation is carried out.
#' @param df data frame to which data is attributed. Must contain fields \code{x} and \code{y}.
#' @param info_df the data-containing data frame. Must contain fields \code{x}, \code{y} and \code{z} and points must lie on an ordered grid.
#' @param miss_value the value attributed to data points which cannot be mapped.
#' @param averaging_buffer a distance radius around the values which cannot be matched to the data. The missing value will then be replaced with the average of data found in this radius. The number should be in the same lenght units as the dataset.
#' @return a vector of the data mapped onto the coordinates in \code{df}
#' @export
#' @examples
#' df <- data.frame(x=runif(100),y=runif(100))
#' info_df <- cbind(expand.grid(x=seq(0,1,0.1),y=seq(0,1,0.1)),z=runif(121))
#' df$z <- attribute_data(df=df,info_df = info_df)
attribute_data <- function(df,info_df,miss_value = 0,averaging_buffer=NA) {
  require("raster")

  xgrid <- sort(unique(info_df$x))
  ygrid <- sort(unique(info_df$y))
  xdiff <- stat_mode(diff(xgrid))
  ydiff <- stat_mode(diff(ygrid))
  
  if(!(xdiff == ydiff)) stop("Currently method only implemented for symmetric grids")
  
  # The points are a regular grid, so go an perfom a merge
  # Firstly do a simple merge (the mesh node takes the value of the cell it
  # is situated in)

  info_df_ras <- rasterFromXYZ(info_df)
  df$cont_vals <- extract(info_df_ras, df[c("x", "y")], method="simple")


  # Now check if there is an averaging buffer when dealing with missing values
  if(!(is.na(averaging_buffer))){
    stopifnot(is.numeric(averaging_buffer) & length(averaging_buffer == 1))
    
    # Perform the merge again with just missing values and an averging radius buffer
    # the mean of the values within the point radius (excluding NA's) are returned
    df[is.na(df$cont_vals), "cont_vals"] <- extract(info_df_ras, df[is.na(df$cont_vals),c("x", "y")],
     method="simple", buffer=averaging_buffer, fun=mean, na.rm=TRUE)
    
    df$cont_vals[is.na(df$cont_vals)] = miss_value
  
  } else {
    df$cont_vals[is.na(df$cont_vals)] = miss_value
  }
  
  return(df$cont_vals)
  
}


#' @rdname setPrecision
#' @aliases setPrecision,GMRF-method
setMethod("setPrecision",signature(.Object="GMRF"),function(.Object,Q) {
  .Object@Q <- Q
  .Object
})

#' @rdname setPrecision
#' @aliases setPrecision,GMRF_basis-method
setMethod("setPrecision",signature(.Object="GMRF_basis"),function(.Object,Q) {
  .Object@G <- setPrecision(.Object@G,Q)
  .Object
})

