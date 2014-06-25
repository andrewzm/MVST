## ----lib-load------------------------------------------------------------
library(Matrix)
library(ggplot2)
library(MVST)

## ----data-load-----------------------------------------------------------
data(icesat)
data(surf_fe)
data(shapefiles)

## ----taxis-set-----------------------------------------------------------
t_axis=0:6

## ----md5-----------------------------------------------------------------
md5_wrapper <- md5_cache("~/cache/")

## ----mesh-init-----------------------------------------------------------
Mesh <- md5_wrapper(initFEbasis,p=surf_fe$p,
                                t=surf_fe$t,
                                M=surf_fe$M,
                                K=surf_fe$K)

## ----mesh_attribute------------------------------------------------------
Mesh["island"] <- md5_wrapper(attribute_polygon,Mesh[c("x","y")],shapefiles$Islands)
Mesh["in_land"] <- md5_wrapper(attribute_polygon,Mesh[c("x","y")],shapefiles$grounding_sub) | 
                                                 Mesh["island"]
Mesh["in_coast"] <-md5_wrapper(attribute_polygon,Mesh[c("x","y")],shapefiles$coast_sub) | 
                                                 Mesh["island"]

## ----mesh_mass-----------------------------------------------------------
Mesh["mass_GT_per_year"] <- Mesh["area_tess_km2"] * 0.5/1000

## ----PlotAnt-------------------------------------------------------------
 PlotAntarctica <- function(g,shapefiles) {
   AIS_x <- c(-2800,2800)
   AIS_y <- c(-2500,2500)
   g <- g +
     geom_path(data = shapefiles$grounding_sub,aes(x,y),colour="black") +
     geom_path(data= shapefiles$coast_sub,aes(x,y),linetype=2,size=0.5) +
     geom_path(data=shapefiles$Islands,aes(x,y,group=id))  + xlab("x (km)") + ylab("y (km)") +
     coord_fixed(xlim=AIS_x,ylim=AIS_y)
 }

## ----plot_coast,fig.height=4---------------------------------------------
g <- plot(Mesh,"in_coast",size=1)
print(PlotAntarctica(g,shapefiles))

## ----icesat_obs----------------------------------------------------------
icesat_obs <- Obs(df=icesat,
                  abs_lim = 5,
                  avr_method = "median",
                  box_size=100,
                  name="icesat")

## ----plot_obs,fig.height=4-----------------------------------------------
g <- plot(subset(icesat_obs,t==3),"z",pt_size=1,max=0.5,min=-0.5)
print(PlotAntarctica(g,shapefiles))

## ----pseudo1-------------------------------------------------------------
pseudo_single_frame <- subset(Mesh,in_coast==F,select=c(x,y))
pseudo_single_frame$z <- 0
pseudo_single_frame$std <- 1e-6

## ----pseudo2-------------------------------------------------------------
pseudo_df <- plyr::adply(t_axis,1,function(x) {
                            return(cbind(pseudo_single_frame,data.frame(t=x)))
                            })

## ----pseudo3-------------------------------------------------------------
pseudo <- Obs(df=pseudo_df,name="pseudo")

## ----mu------------------------------------------------------------------
 mu <- function(k) { # time-varying matrix for mu
   return(matrix(0,nrow(Mesh),1)) 
 }

## ----A-------------------------------------------------------------------
A <- function(k) {  # time-varying matrix for A
   0.2*Imat(nrow(Mesh))
 }

## ----Q-------------------------------------------------------------------
 Q <- function(k) {  # time-varying matrix for Q
    Prec_from_SPDE_wrapper(M = mass_matrix(Mesh),
                           K = stiffness_matrix(Mesh),
                           nu = 1,
                           desired_prec = 1, 
                           l = 1200) 
}

## ----SURF_VAR------------------------------------------------------------
SURF_VAR <- VAR_Gauss( mu = mu,A=A, Qw = Q,t_axis = t_axis,name="SURF")

## ----SURF----------------------------------------------------------------
SURF <- GMRF_basis(G = SURF_VAR, Basis = Mesh)

## ----Cmat----------------------------------------------------------------
L1 <- link(SURF,icesat_obs)
L2 <- link(SURF,pseudo)

## ----Graph---------------------------------------------------------------
e <- link_list(list(L1,L2))
v <- block_list(list(G1=SURF,O1=icesat_obs,O2 = pseudo))
G <- Graph(e=e,v=v)

## ----Graph_reduced-------------------------------------------------------
G_reduced <- compress(G)

## ----Infer,cache=TRUE----------------------------------------------------
 Results <- Infer(G_reduced)

## ----Plot----------------------------------------------------------------
Mesh["to_plot"] <- c(subset(Results$Post_GMRF,t==2,select=x_mean))

## ----Plot2---------------------------------------------------------------
 Mesh['marg_std'] <- sqrt(subset(Results$Post_GMRF,t==2,select=x_margvar))

## ----Plot3---------------------------------------------------------------
g <- plot_interp(Mesh,"to_plot",ds=500,min=-0.5,max=0.5,leg_title="m/yr")

## ----Plot4---------------------------------------------------------------
g <- plot(Mesh,g=PlotAntarctica(g,shapefiles),plot_dots=F)

## ----Plot5,fig.height=6,dev='png',dpi=800--------------------------------
g <- g +  geom_point(data = subset(Mesh,abs(to_plot) > (marg_std)),
                     aes(x,y),size=1,colour=scales::muted("green")) +
    theme(text = element_text(size=15))
print(g)

## ----getDf---------------------------------------------------------------
graph_df <- getDf(G@v[[1]])

## ----Comb----------------------------------------------------------------
Comb <- rbind(((graph_df$t== 2) & (graph_df$in_land == T))*graph_df$mass_GT_per_year,
                       ((graph_df$t== 4) & (graph_df$in_land == T))*graph_df$mass_GT_per_year)

## ----Results_comb,cache=TRUE---------------------------------------------
Results_linear_comb <- Infer(G_reduced,Comb = Comb)

## ----Results_comb2-------------------------------------------------------
print(Results_linear_comb$Comb_results)

