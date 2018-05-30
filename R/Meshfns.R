
#' @title Export INLA Mesh
#' @description Exports the 'inla.mesh' class object to CSV files that can be read by 
#' the 'Load_meshes' function of the MVST package. The MVST R package requires the FEM mesh 
#' components to be stored in a specific
#' csv file format in order to be properly imported and constructed. 
#'
#' @param mesh The INLA mesh to export, of class 'inla.mesh'
#' @param directory directory in which to save output files, must end with a '/'
#' @return creates four csv files ‘p.csv, t.csv,
#' M.csv’ and ‘K.csv’ which are csv files containing the vertices,
#' triangulations, mass matrix and stiffness matrix details
#' respectively. The matrices are stored in the three-column format ‘i,j,k’.
#' @export
write_INLA_mesh_MVST <- function(mesh, directory){
  # First perform check to make sure the right object is passed
  if (class(mesh) != "inla.mesh"){stop("Mesh object must be of class 'inla.mesh'")}
  # Now perform checks to see if the directoryectory exists and is passed correctly
  
  if (file.exists(directory) == TRUE & substring(directory, nchar(directory)) == "/"){
     # Passed arguments are correct, now go and write the output files
     # Get P matrix (vertex locations)
     p <- mesh$loc[,1:2]
     write.table(p, file = paste0(directory, 'p.csv'), quote = F, sep=',', row.names = F,
                 col.names=F)
     cat(paste0("vertices location file written to ", directory, 'p.csv \n'))
     # Get t Matrix
     t <- mesh$graph$tv
     write.table(t, file = paste0(directory, 't.csv'), quote = F, sep=',', row.names=F, 
                 col.names=F)
     cat(paste0("triangulation matrix written to ", directory, 't.csv \n'))
     # Now calculate FEM matrices (to get stiffness and mass matrix)
     fem <- inla.mesh.fem(mesh)
     # Get M matrix
     # Get a summary of the mass matrix (as a 3xn i, j, k matrix)
     M <- summary(fem$c1)
     write.table(M, file = paste0(directory, 'M.csv'), quote = F, sep=',', row.names=F,
                 col.names = F)
     cat(paste0("mass matrix written to ", directory, 'M.csv \n'))
     # Get k matrix
     # Get a summary of the stiffness matrix (as a 3xn i, j, k matrix)
     k <- summary(fem$g1)
     write.table(k, file = paste0(directory, 'K.csv'), quote = F, sep=',', row.names=F,
                 col.names=F)
     cat(paste0("stiffness matrix written to ", directory, 'K.csv \n'))
     
  } else if (file.exists(directory) == FALSE){
    stop("Mesh output directory doesn't exist")
  } else if (substring(directory, nchar(directory)) != "/") {
    stop("directory name must end with a '/'")
  }
}

#' @title Convert INLA mesh object MVST FEM class
#'
#' @description The MVST R package requires the FEM mesh components to be in a different object class 
#' to that of 'inla.mesh'. This function converts the 'inla.mesh' class objects to required 
#' MVST class, without export of the mesh to disk. 
#' @param mesh The INLA mesh to convert (inla.mesh))
#' @return MVST finite element mesh object
#' @export
convert_INLA_mesh_MVST <- function(mesh){
  require("MVST")
  
  # First perform check to make sure the right object is passed
  if (class(mesh) != "inla.mesh"){stop("Mesh object must be of class 'inla.mesh'")}
  # Passed arguments are correct, now go and write the output files
  # Get P matrix (vertex locations)
  p <- mesh$loc[,1:2]
  # Get t Matrix
  t <- mesh$graph$tv
  # Now calculate FEM matrices (to get stiffness and mass matrix)
  fem <- inla.mesh.fem(mesh)
  # Get M matrix
  # Get a summary of the mass matrix (as a 3xn i, j, k matrix)
  M <- summary(fem$c1)
  # Get k matrix
  # Get a summary of the stiffness matrix (as a 3xn i, j, k matrix)
  K <- summary(fem$g1)
  
  # Now create MVST mesh object 
  p <-  round(as.matrix(p))
  tri <- as.matrix(t)
  M <-  as.matrix(M)
  K <-  as.matrix(K)
  n <- nrow(p)
  M <- sparseMatrix(i=M[,1],j=M[,2],x=M[,3],dims=c(n,n))
  K <- sparseMatrix(i=K[,1],j=K[,2],x=K[,3],dims=c(n,n))
  
  return(initFEbasis(p=p, t = tri, M = M, K = K))
}