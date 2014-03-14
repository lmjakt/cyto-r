source("../../R/functions.R")

data <- read.files(c("PI", "PE", "APC", "PE-Cy7", "549_1", "549_1_d", "549_5", "549_5_d", "734_1", "734_1_d", "763_2", "763_2_d",
                     "763_2_t", "763_2_dt"))

## mismatched colnames..
par.names <- colnames(data[[1]])
for(i in 1:length(data))
  colnames(data[[i]]) <- par.names

data.f <- filter.simple.list(data, filter.1) ## mainly to remove 0 values.. which make trouble for log transform

flk1.v <- make.v.list(data.f, 'FLK1')
pdgfra.v <- make.v.list(data.f, 'PDGFRA')
cdh5.v <- make.v.list(data.f, 'CDH5')
gfp.v <-  make.v.list(data.f, 'GFP')

plot.histograms(flk1.v[5:14], 30)
plot.histograms(pdgfra.v[5:14], 30)
plot.histograms(cdh5.v[5:14], 30)
plot.histograms(gfp.v[5:14], 30)

de.flk1.pdgfra <- make.kde2d.list(data.f, 'FLK1', 'PDGFRA', n=75)
de.flk1.cdh5 <- make.kde2d.list(data.f, 'FLK1', 'CDH5', n=75)
de.gfp.flk1 <- make.kde2d.list(data.f, 'GFP', 'FLK1', n=75)
de.gfp.pdgfra <- make.kde2d.list(data.f, 'GFP', 'PDGFRA', n=75)

x.y.density(de.flk1.pdgfra, column.no=4, xlines=c(2.8), ylines=c(2.4))
x.y.density(de.gfp.pdgfra, column.no=4, xlines=c(2.5), ylines=c(2.4))
x.y.density(de.gfp.flk1, column.no=4, xlines=c(2.5), ylines=c(2.4))

## on the facs machine I set a compensation for GFP / PE overlap
## that was specified as simply PE = PE - (0.15 * GFP)
## so we can try that here.
data.fc <- compensate(data.f, 'PDGFRA', 'GFP', 0.15)

de.flk1.pdgfra.c <- make.kde2d.list(data.fc, 'FLK1', 'PDGFRA', n=75)
de.flk1.cdh5.c <- make.kde2d.list(data.fc, 'FLK1', 'CDH5', n=75)
de.gfp.flk1.c <- make.kde2d.list(data.fc, 'GFP', 'FLK1', n=75)
de.gfp.pdgfra.c <- make.kde2d.list(data.fc, 'GFP', 'PDGFRA', n=75)

x.y.density(de.flk1.pdgfra.c, column.no=4, xlines=c(2.8), ylines=c(2.4))
x.y.density(de.gfp.pdgfra.c, column.no=4, xlines=c(2.5), ylines=c(2.4))
x.y.density(de.gfp.flk1.c, column.no=4, xlines=c(2.5), ylines=c(2.4))



# d = list of data tables,
# c1 and c2 = colnameas
# cf = compensation factor
# c1 <- c1 - (c2 * cf)
## since c1 can now be negative we can choose to auto.ajust it

## compensate function implemented in functions.R .. where it belongs
compensate <- function(d, c1, c2, cf, auto.adjust=TRUE){
  for(i in 1:length(data)){
    c1.min <- min(d[[i]][,c1])
    d[[i]][,c1] <- d[[i]][,c1] - (cf * d[[i]][,c2])
    if(auto.adjust){
      c1.c.min <- min(d[[i]][,c1])
      d[[i]][,c1] <- d[[i]][,c1] + (c1.min - c1.c.min)
    }
  }
  d
}
