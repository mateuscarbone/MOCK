file.choose()
source("codigos_MOCK.R")
library(plotly)
library(mclust)
# dados ficticios ---------------------------------------------------------


am1<-c(rnorm(100,2,.2),0,  rnorm(100,1,.2),rnorm(100,0,.2),2.5,1)
am2<-c(rnorm(100,2,.2),2.5,rnorm(100,1,.2),rnorm(100,0,.2),0,1.2)
am3<-c(rnorm(100,2,.2),2.5,rnorm(100,1,.2),rnorm(100,0,.2),0,1.2)
dados<-data.frame(cbind(am1,am2,am3))

plot_ly(dados,x=~am1,y=~am2,z=~am3,type = "scatter3d")

grupos<-c(rep(1,100),1,rep(2,100),rep(3,100),2,2)
plot_ly(dados,x=~am1,y=~am2,z=~am3,color=~factor(grupos),type = "scatter3d")

# MOCK --------------------------------------------------------------------

start.time <- Sys.time()
res.mock<-MOCK(dados,5)
end.time <- Sys.time()
end.time-start.time

plot_ly(dados,x=~am1,y=~am2,z=~am3,color=~factor(res.mock$group),type = "scatter3d")
adjustedRandIndex(res.mock$group,grupos)
