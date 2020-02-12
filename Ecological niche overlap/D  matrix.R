# 
# # The script was modified from the supplementary scripts of
# # 'Broennimann, O., Fitzpatrick, M.C., Pearman, P.B., Petitpierre, B., Pellissier, L., Yoccoz, N.G.,
# # Thuiller, W., Fortin, M.J., Randin, C., Zimmermann, N.E. and Graham, C.H., 2012. 
# # Measuring ecological niche overlap from occurrence and spatial environmental data.
# # Global ecology and biogeography, 21(4), pp.481-497.'

## functions were written by Olivier Broennimann. Departement of Ecology and Evolution (DEE). 
## 2nd of November 2010. University of Lausanne. Department of ecology and evolution, Switzerland

rm(list = ls())
source("niche.overlap.functions.R")
source("occ.prep.functions.R")

library(BIOMOD)
library(ade4)
library(adehabitat)
library(sp)
library(gam)
library(MASS)
library(mvtnorm)
library(gbm)
library(dismo)
library(ggplot2)  
library(reshape2)

#################################################################################################
############################## preparation of datasets ##########################################
#################################################################################################

load('climatic.Rdata')

climate<-pder[,c(19,20,1:16)]
occ<-read.csv('Occurrence_probability_matrix.csv')
sp.name<-unique(occ$species)
n<-length(sp.name)
D=matrix(NA,nrow = n,ncol = n)
zz<-list()

for (i in 1:n) {
  for (j in 1:n) {
    
    # load climate variable for all site of the study area  (column names should be x,y,X1,X2,...,Xn)
    clim<-na.exclude(climate)
    sp1<-occ[which(occ$species==sp.name[i] & occ$Presence==1),c(3,4)]
    sp2<-occ[which(occ$species==sp.name[j] & occ$Presence==1),c(3,4)]    
    
    # loading occurrence sites for the species (column names should be x,y)
    occ.sp.aggr<-sp1
    occ.sp.aggr1<-sp2
    
    # remove occurrences closer than a minimum distance to each other (remove aggregation). Setting min.dist=0 will remove no occurrence.
    occ.sp1<-occ.desaggragation(df=occ.sp.aggr,colxy=1:2,min.dist=0.00833,plot=F)
    occ.sp2<-occ.desaggragation(df=occ.sp.aggr1,colxy=1:2,min.dist=0.00833,plot=F)
    
    # create sp occurrence dataset by adding climate variables from the global climate datasets
    # resolution should be the resolution of the climate data grid
    occ.sp1<-na.exclude(sample.sp.globvar(dfsp=occ.sp1,colspxy=1:2,colspkept=NULL,dfvar=clim,colvarxy=1:2,colvar="all",resolution=0.00833))
    occ.sp2<-na.exclude(sample.sp.globvar(dfsp=occ.sp2,colspxy=1:2,colspkept=NULL,dfvar=clim,colvarxy=1:2,colvar="all",resolution=0.00833))
    
    #################################################################################################
    ############################## ANALYSIS - selection of parameters ###############################
    #################################################################################################
    
    # selection of the type of analysis.
    # If PROJ =F, the models are calibrated on both ranges.
    # If PROJ =T, the models are calibrated on species 1 range only and projected to range 2. 
    # Analyses where both ranges are needed (ex: LDA) are not done
    PROJ = F
    
    # selection of variables to include in the analyses
    Xvar<-c(3:18)
    nvar<-length(Xvar)
    
    #number of interation for the tests of equivalency and similarity
    iterations<-100
    
    #resolution of the gridding of the climate space
    R=100
    
    #################################################################################################
    ################### row weigthing and grouping factors for ade4 functions  ######################
    #################################################################################################
    clim1<-clim
    clim2<-clim
    clim12<-rbind(clim,clim)
    # if PROJ = F
    row.w.1.occ<-1-(nrow(occ.sp1)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ1
    row.w.2.occ<-1-(nrow(occ.sp2)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ2
    row.w.occ<-c(rep(0, nrow(clim1)),rep(0, nrow(clim2)),rep(row.w.1.occ, nrow(occ.sp1)),rep(row.w.2.occ, nrow(occ.sp2)))
    
    row.w.1.env<-1-(nrow(clim1)/nrow(clim12))  # prevalence of clim1
    row.w.2.env<-1-(nrow(clim2)/nrow(clim12))  # prevalence of clim2
    row.w.env<-c(rep(row.w.1.env, nrow(clim1)),rep(row.w.2.env, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
    
    fac<-as.factor(c(rep(1, nrow(clim1)),rep(2, nrow(clim2)),rep(1, nrow(occ.sp1)),rep(2, nrow(occ.sp2))))
    
    # if PROJ = T
    
    row.w.occ.PROJT<-c(rep(0, nrow(clim1)),rep(0, nrow(clim2)),rep(1, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
    row.w.env.PROJT<-c(rep(1, nrow(clim1)),rep(0, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
    
    # global dataset for the analysis and rows for each sub dataset
    data.env.occ<-rbind(clim1,clim2,occ.sp1,occ.sp2)[Xvar]
    row.clim1<-1:nrow(clim1)
    row.clim2<-(nrow(clim1)+1):(nrow(clim1)+nrow(clim2))
    row.clim12<-1:(nrow(clim1)+nrow(clim2))
    row.sp1<-(nrow(clim1)+nrow(clim2)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1))
    row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+nrow(occ.sp2))
    
    #################################################################################################
    #################################### PCA-ENV ####################################################
    #################################################################################################
    
    # measures niche overlap along the two first axes of a PCA calibrated on all the pixels of the study areas
    
    if(PROJ == F){	#fit of the analyse using occurences from both ranges		
      pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2)
    }
    if(PROJ == T){	#fit of the analyse using occurences from range 1
      pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env.PROJT, center = T, scale = T, scannf = F, nf = 2)
    }
    # predict the scores on the axes
    scores.clim12 <- pca.cal$li[row.clim12,]
    scores.clim1 <- pca.cal$li[row.clim1,]
    scores.clim2 <- pca.cal$li[row.clim2,]
    scores.sp1 <- pca.cal$li[row.sp1,]
    scores.sp2 <- pca.cal$li[row.sp2,]
    
    # calculation of occurence density and test of niche equivalency and similarity 
    z1 <- grid.clim(scores.clim12,scores.clim1,scores.sp1,R)
    zz[[i]] <- z1
    z2 <- grid.clim(scores.clim12,scores.clim2,scores.sp2,R)
    D[i,j] <- round(as.numeric(niche.overlap(z1,z2,cor=T)[1]),4)
    }
}
colnames(D) <- sp.name
rownames(D) <- sp.name
D.ij <- round(D,3)
write.csv(D.ij,file = 'D matrix.csv')
# D<-read.csv(file = 'D matrix.csv')
save(zz,file = 'z.rdata')
load('z.rdata')


theme_Tian <- function() {
  theme_bw() +
    theme(panel.grid = element_blank()) + # no grid lines
    theme(strip.background = element_rect(fill = "grey95", colour = "grey30"))+
    theme(axis.title = element_text(colour = 'black',size=8,face = 'plain'))+
    theme(axis.text = element_text(colour = 'black',size=8,face = 'plain'))+
    theme(legend.text = element_text(colour = 'black',size = 8,face = 'plain'))+
    theme(legend.title = element_text(colour = 'black',size=8,face = 'plain'))+
    theme(strip.text = element_text(colour = 'black',size=8,face = 'plain'))
}
theme_set(theme_Tian())

rownames(zz[[1]]$z.uncor) <- zz[[1]]$x
colnames(zz[[1]]$z.uncor) <- zz[[1]]$y
pl <- melt(zz[[1]]$z.uncor)
pl$sp <- as.character(sp.name[1])

for (i in 2:7) {
  rownames(zz[[i]]$z.uncor) <- zz[[i]]$x
  colnames(zz[[i]]$z.uncor) <- zz[[i]]$y
  pl.i <- melt(zz[[i]]$z.uncor)
  pl.i$sp <- as.character(sp.name[i])
  pl <- rbind(pl, pl.i)
}

plot <- ggplot(pl, aes(Var1, Var2)) + geom_tile(aes(fill = value)) +
  scale_fill_gradient(low =  "#d11141", high = "darkblue") +
  facet_wrap(. ~ sp, nrow = 2) +
  ylim(min(pl$Var2), max(pl$Var2)) +
  xlim(min(pl$Var1), max(pl$Var1)) +
  theme(legend.justification = c(1, 0),
        legend.position = c(1, 0)) +
  xlab('PC1') + ylab('PC2') +
  labs(fill = "Density") +
  theme(strip.text = element_text(
    colour = 'black',
    size = 8,
    face = 'italic'
  ))

ggsave(
  plot,
  file = 'density.png',
  width = 20,
  height = 10,
  unit = 'cm'
)

    