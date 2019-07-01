source("codigos\\MOCK.R")

# dados ficticios ---------------------------------------------------------


am1<-c(rnorm(600,2,.2),0,  rnorm(600,.5,.2),rnorm(600,0,.2),2.5,1)
am2<-c(rnorm(600,2,.2),2.5,rnorm(600,.5,.2),rnorm(600,0,.2),0,1.2)
am3<-c(rnorm(600,2,.2),2.5,rnorm(600,.5,.2),rnorm(600,0,.2),0,1.2)
dados<-data.frame(cbind(am1,am2,am3))

plot_ly(dados,x=~am1,y=~am2,z=~am3,type = "scatter3d")

grupos<-c(rep(1,600),1,rep(2,600),rep(3,600),2,2)
plot_ly(dados,x=~am1,y=~am2,z=~am3,color=~factor(grupos),type = "scatter3d")

# MOCK --------------------------------------------------------------------

start.time <- Sys.time()
res.mock<-MOCK(dados,5)
end.time <- Sys.time()
end.time-start.time

plot_ly(dados,x=~am1,y=~am2,z=~am3,color=~factor(res.mock$group),type = "scatter3d")
adjustedRandIndex(res.mock$group,grupos)