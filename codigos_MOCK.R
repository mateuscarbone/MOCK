library(parallel)
library(spatstat)
library(Rfast)

# Cruzamento --------------------------------------------------------------



crox <- function (li,cruz) {
  flh <- lapply(cruz,ifelse,li[[1]],li[[2]])
  return(flh)
}



# Algoritmo de Prim -------------------------------------------------------


prim <- function(diss){
  conec <- nrow(diss)
  vrep <- rep(NA,conec)
  part <- 1
  cheg <- which.min(diss[part,])
  vrep[part] <- cheg
  seq <- c(part,cheg)
  mtn <- diss[seq,-seq]
  while(!is.null(nrow(mtn))){
    cval <- colnames(mtn)
    lval <- rownames(mtn)
    vv<-which.min(apply(mtn,1,min))
    part<-as.numeric(lval[vv])
    cheg<-as.numeric(cval[which.min(mtn[vv,])])
    if(is.na(vrep[part])){
      vrep[part] <- cheg
    }else{
      vrep[cheg] <- part
    }
    seq <- c(seq,cheg)
    mtn <- diss[seq,-seq]
  }
  nfa <- which(is.na(vrep))
  nfalt <- nfa[!nfa%in%vrep]
  vrep[nfalt] <- as.numeric(names(which.min(diss[,nfalt])))
  vrep[nfa[nfa%in%vrep]] <- nfa[nfa%in%vrep]
  return(vrep)
}


# Link interessante -------------------------------------------------------


interes <- function(obj,diss){
  diag(diss)=NA
  dist <- vector("numeric")
  int <- vector("numeric")
  for (i in 1:length(obj)){
    a <- i
    b <- obj[i]
    c1<-order(diss[,a])[1:10]
    c2<-order(diss[,b])[1:10]
    if(!any(c1%in%c2)){
      int <- c(int,i)
      dist <- c(dist,diss[a,b])
    }
  }
  ord<-order(dist,decreasing = T)
  return(list(p.int=int[ord],dist=dist[ord]))
}


# Mutação -----------------------------------------------------------------
mut.f <- function (vrep,mdiss.ord) {
  n <- length(vrep)
  pos.d <- vector("integer",length = n)
  pos.d <- sapply(1:n,function(i)which(vrep[i]==mdiss.ord[,i]))
  pcm<-pm(pos.d,n)
  prob <- runif(n)
  v.mut<-which(ifelse(prob<pcm,T,F))
  if (length(v.mut)!=0) {
    vrep[v.mut] <- apply(as.matrix(mdiss.ord[1:10,v.mut]),2,sample,size=1)
  }
  return(vrep)
}

# mut.f <- function (vrep,diss.o) {
#   pm<-function (X) {1/n+(X/n)^2}
#   n <- length(vrep)
#   pos.d <- vector("integer",length = n)
#   for (i in 1:n){
#     pos.d[i] <- which(vrep[i]==diss.o[,i])
#   }
#   pcm<-pm(pos.d)
#   v.mut<-which(ifelse(runif(n)<pcm,T,F))
#   if (length(v.mut)==0) {
#     return(vrep)
#   }
#   else{
#     vrep[v.mut] <- apply(diss.o[1:10,v.mut],2,sample,size=1)
#   }
#   return(vrep)
# }



# Soluções ligação ---------------------------------------------------------


gsoll <- function(p.int,vrep){
  nsol <- (length(p.int)+1)
  sol <- vector("list",length=nsol)
  sol[[1]] <- vrep
  for (i in 2:nsol){
    vv <- sol[[1]]
    vv[p.int[1:(i-1)]] <- p.int[1:(i-1)]
    sol[[i]] <- vv
  }
  return(sol)
}


# Classificação -----------------------------------------------------------


cls <- function(g){
  N <- length(g)
  ac <- rep(-1,N)
  CA <- 1
  ant <-vector("numeric",length = N)
  for(i in 1:N){
    ctr <- 1
    if(ac[i]==-1){
      ac[i] <- CA
      viz <- g[i]
      ant[ctr] <- i
      ctr <- ctr+1
    }
    while(ac[viz]==-1){
      ant[ctr] <- viz
      ac[viz] <- CA
      viz <- g[viz]
      ctr=ctr+1
    }
    if(ac[viz]!=CA){
      ctr <- ctr-1
      while(ctr>=1){
        ac[ant[ctr]] <- ac[viz]
        ctr <- ctr-1
      }
    }
    else{
      CA <- CA+1
    }
  }
  return(ac)
}


# cls <- function(agr){
#   grupo <- rep(NA,length=length(agr))
#   pp <- 1
#   while(anyNA(grupo)){
#     pos <- min(which(is.na(grupo)))
#     gg <- unique(c(pos,agr[pos]))
#     repeat{
#       n <- length(gg)  
#       pos <- which(agr%in%gg)
#       gg <- unique(c(pos,agr[gg]))
#       if(n==length(gg))break
#     }
#     grupo[gg] <- pp
#     pp <- pp+1
#   }
#   return(as.factor(grupo))
# }


# Conectividade -----------------------------------------------------------

s.conect <- function(gp,mdiss.o) { 
  ns <- length(gp)
  ncol <- ncol(mdiss.o)
  conn <- vector("numeric",length=ns)
  for (i in 1:ns) {
    soma <- vector("numeric",length=ncol)
    g.grup <- gp[[i]]
    for(j in 1:ncol) {
      gpi <- g.grup[j]
      pos <- mdiss.o[,j]
      ss <- which(gpi!=g.grup[pos])
      soma[j] <- sum(1/ss)
    }
    conn[i] <- sum(soma)
  }
  return(conn)
}


# s.conect2 <- function(i,gp,mdiss.o) { 
#   pos <- which(gp[mdiss.o[,i]]!=gp[i])
#   value <- sum(1/pos)
#   return(value)
# }


# conect <- function(d.ord,vg){
#   return(sum(apply(d.ord[1:10,],2,function(X)sum(1/which(X%in%vg)))))
# }


# Compactação -------------------------------------------------------------



# compac <- function(data,grup){
#   mat<-apply(data,2,function(X)ave(X,grup,FUN = mean))
#   dt.eu<-rowSums((data-mat)^2)
#   dist<-sum(dt.eu)
#   return(dist)
# }


# compac <- function(data,grup){
#   mat<-apply(data,2,function(X)ave(X,grup,FUN = mean))
#   dt.eu<-sqrt(rowSums((data-mat)^2))
#   dist<-sum(dt.eu)
#   return(dist)
# }

# 
# compac <- function(data,grup){
#   sep <- split(data,grup)
#   cent <- lapply(sep,colMeans)
#   for (i in 1:length(cent)){
#     sep[[i]] <- (t(t(sep[[i]])-cent[[i]]))^2
#   }
#   sep<- plyr::ldply(sep,.id=NULL)
#   dt.eu<-sqrt(rowSums(sep))
#   dist <- sum(dt.eu)
#   return(dist)
# }


compac <- function(data,grup){
  sep <- split(data,grup)
  cent <- lapply(sep,colMeans)
  soma <- vector("numeric",length=length(cent))
  for (i in 1:length(cent)){
    soma[i] <-sum(Rfast::dista(t(cent[[i]]),sep[[i]]))
  }
  dist <- sum(soma)
  return(dist)
}






# Codificação K-Means -----------------------------------------------------



kme.s <- function(gg,snp,diss.o){
  nn <- length(gg)
  vv<-vector("logical",length = nn)
  for (i in 1:nn){
    vv[i] <- gg[i]!=gg[snp[i]]
  }
  kme.s <- snp
  pos <- which(vv)
  kme.s[pos] <- apply(as.matrix(diss.o[1:10,pos]),2,sample,size=1)
  return(kme.s)
}



# Soluções iniciais K-Means -----------------------------------------------



s.km <- function(data,fsize,n.si,diss.o,agm){
  rkm <- vector("list",length=(fsize-n.si-1))
  for ( i in 2:(fsize-n.si)){
    res<-kmeans(data,i)
    rkm[[i-1]] <- res$cluster
  }
  slk<-lapply(rkm,kme.s,snp=agm,diss.o=diss.o)
  return(slk)
}


# Soluções não dominadas --------------------------------------------------



# s.nd <- function(con,comp){
#   nsol <- length(con)
#   n.dom <- rep(F,nsol)
#   for( i in 1:nsol){
#     nd<-!any((con[i]<=con)+(comp[i]<=comp)==0)
#     if(nd){
#       n.dom[i] <- T
#     }
#   }
#   return(n.dom)
# }



# s.nd <- function(con,comp){
#   nsol <- length(con)
#   n.dom <- vector("logical",length=nsol)
#   for( i in seq_along(n.dom)){
#     n.dom[i] <- !any((con[i]<=con[-i])+(comp[i]<=comp[-i])==0)
#   }
#   dup <- which(duplicated(con)+duplicated(comp)==2)
#   n.dom[dup] <- F
#   return(n.dom)
# }


s.nd <- function(vet.obj){
  v2 <- vet.obj
  n <- length(vet.obj)
  snd<-vector("logical",length = n)
  v.snd <- vector("double")
  snd[1] <- T
  v.snd[1] <- v2[1]
  for( i in 2:n) {
    if (any(v2[i] >= v.snd)) {
      snd[i] <- F
    }else{
      snd[i] <- T
      v.snd <- append(v.snd,v2[i])
    }
  }
  return(snd)
}


# s.nd2 <- function(vet.obj){
#   v2 <- as.matrix(vet.obj)[,2]
#   n <- nrow(vet.obj)
#   snd<-vector("logical",length = n)
#   v.snd <- vector("double")
#   snd[1] <- T
#   v.snd[1] <- v2[1]
#   for( i in 2:n) {
#     if (any(v2[i] >= v.snd)) {
#       snd[i] <- F
#     }else{
#       snd[i] <- T
#       v.snd <- append(v.snd,v2[i])
#     }
#   }
#   return(snd)
# }

# cn <- rnorm(1000)
# cp <- rnorm(1000)
# plot(cn,cp)
# 
# 
# start.time <- Sys.time()
# s.nd(cn,cp)
# end.time <- Sys.time()
# print(end.time-start.time)



# squeze factor -----------------------------------------------------------



sqze <- function(x,y){
  int <- seq(0,1,.1)
  a<-findInterval(x, int)
  b<-findInterval(y, int)
  sqf <- ave(a,a,b,FUN = table) ### contagem em cada intervalo
  areac<-paste0(a,b) ### toda as areas
  res <-  !duplicated(areac)
  area <- areac[res]
  sqf <- sqf[res]
  return(list(sqf=sqf,area=area,areac=areac))
}


# dados simulados ---------------------------------------------------------

dsimu2 <- function(d){ 
  minimo <- sapply(d,min)
  maximo <- sapply(d,max)
  n <- nrow(d)
  mat <- matrix(0,ncol=ncol(d),nrow=nrow(d))
  for( i in 1:ncol(d)){
    mat[,i] <- runif(n,minimo[i],maximo[i])    
  }
  return(mat)
}



dsimu <- function(d) {
  d <- as.matrix(d)
  nr <- nrow(d)
  nc <- ncol(d)
  pca <- prcomp(d)
  min <- apply(pca$x,2,min)
  max <- apply(pca$x,2,max)
  ls<-vector("list",nc)
  for(i in 1:nc){
    ls[[i]] <- runif(nr,min[i],max[i])
  }
  simd<-matrix(unlist(ls),ncol=nc)
  simd<- simd %*% t(pca$rotation) + pca$center # tranformando para o espaÃ§o original
  return(simd)
}

# function(d) {
#   d <- as.matrix(d)
#   nr <- nrow(d)
#   nc <- ncol(d)
#   pca<- prcomp(d)
#   min <- apply(pca$scores,2,min)
#   max <- apply(pca$scores,2,max)
#   ls<-vector("list",nc)
#   for(i in 1:nc){
#     ls[[i]] <- runif(nr,min[i],max[i])
#   }
#   simd<-matrix(unlist(ls),ncol=nc)
#   simd<- simd %*% t(pca$loadings) + pca$center # tranformando para o espaÃ§o original
#   return(simd)
# }

# probabilidade mutação ---------------------------------------------------

pm<-function (X,n) 1/n+(X/n)^2



# nomalização -------------------------------------------------------------

normalize <- function(x) {
  return((x- min(x)) /(max(x)-min(x)))
}

############################################################
############################################################
############## inicialia??o ################################



sol.front <- function (data.cod,kcls=25) {
  data<-data.cod$data
  #data=dados
  n <- nrow(data)
  fsize <- 2*kcls
  ### diss = matriz de dissimilaridade
  diss <- as.matrix(dist(data))
  diag(diss)=NA
  ### diss.o = matriz de dissimilaridade ordenado
  diss.o<-apply(diss,2,order)
  ### vrep = vetor resposta algo. de prim
  vrep <- prim(diss)
  ### l.int = links interessantes
  l.int<-interes(vrep,diss)
  #n.si = numero de solu??es iniciais
  n.si <- min((fsize*.5-1),length(l.int$p.int))
  #cjl = geraÃ§Ã£o dos possiveis vetores removendo os links interessantes
  cjl <- gsoll(l.int$p.int,vrep)
  #solu??es iniciais k-meas
  sikm <-s.km(data,fsize,n.si,diss.o,cjl[[1]])
  ###solu??es iniciais geral
  sli <- c(sikm,cjl)
  #grup = classifica??o para cada sli()
  grup <- lapply(sli,cls)
  #s.con = resultados conectividade
  s.con <- s.conect(grup,diss.o[1:10,])
  ord <- order(s.con)
  s.con <- s.con[ord]
  #s.comp = resultado da compacta??o
  s.comp <- sapply(grup,compac,data=data)
  s.comp <- s.comp[ord]
  #snd = solu??es n?o dominads
  snd<-s.nd(s.comp)
  #b <- proc.time()
  #tdec <- b-a
  ######### popula??o externa
  PE <- sli[ord][snd]
  grup <- grup[ord][snd]
  s.con <- s.con[snd]
  s.comp <- s.comp[snd]
  s.conP<-normalize(s.con)
  s.compP <- normalize(s.comp)
  for( j in 1:1000){
    ##### esqueze factor
    sfnd <- sqze(s.conP,s.compP)
    ##problemas se area for de tamanho 1
    ################################# loop principlai
    PI <- vector("list",length=10)
    PIt <- vector("list",length=10)
    #####loop PIt
    for (i in 1:10){
      pos<-which(sfnd$area%in%sample(sfnd$area,2)) ### seleciona duas areas aleaorias
      msqz <- which.min(sfnd$sqf[pos]) ## menor squeeze factor
      pri <- sfnd$area[pos[msqz]] ## seleciona a Ã¡rea
      sas <- which(sfnd$areac%in%pri) ## seleciona as soluÃ§Ãµes
      if(length(sas)<=1){ ## se for unica Ã© ela mesma, senÃ£o sorteia-se uma
        PIt[[i]] <- PE[[sas]]
      }else{
        PIt[[i]] <- PE[[sample(sas,size=1)]]
      }
    }
    #####loop cruzamento
    for(i in 1:10){#ok
      u<-runif(1)
      if(u<=.7){
        ind <- sample(10,2)
        pcr<-list(sample(c(T,F),n,replace = T))
        flh <- unlist(crox(PIt[ind],pcr))
        PI[[i]]<-mut.f(flh,diss.o)
      }else{
        ind <- sample(10,1)  
        PI[[i]] <- mut.f(PIt[[ind]],diss.o)
      }
    }
    #### avaliuaÃ§Ã£o PI
    gpi <- lapply(PI,cls)
    con.pi <- s.conect(gpi,diss.o[1:10,])
    comp.pi <- sapply(gpi,compac,data=data)
    s.con <- append(s.con,con.pi)
    s.comp <- append(s.comp,comp.pi) # quase
    ord <- order(s.con)
    s.con <- s.con[ord]
    s.comp <- s.comp[ord]
    grup <- append(grup,gpi) 
    grup <- grup[ord]
    PE <- append(PE,PI)
    PE <- PE[ord]
    snd <- s.nd(s.comp) #ok
    PE <- PE[snd] #ok
    grup <- grup[snd]
    s.con <- s.con[snd]
    s.comp <- s.comp[snd]
    s.conP<-normalize(s.con)
    s.compP <- normalize(s.comp)
  }
  if (data.cod$cod.f) {
    return(list(df_med = data.frame(conec=sqrt(s.conP),compa=sqrt(s.compP)),
                SND=PE,group=grup))
  }else{
    return(df_med = data.frame(conec=sqrt(s.conP),compa=sqrt(s.compP)))    
  }
}







MOCK <- function(data,n.sim=5){
  tsim <- ncol(data)>nrow(data)
  func<-c("cls","compac","crox","gsoll","interes","kme.s",
          "mut.f","pm","prim","s.conect","s.km","s.nd","sqze","normalize")
  cont.list <- vector("list",length=6)
  cont.list[[1]] <- list(data=data,cod.f=T)
  print(cat("Alta dimensionalidade:",tsim))
  if(tsim){
    for (i in 2:(n.sim+1)){
      cont.list[[i]]  <- list(data=data.frame(dsimu2(data)),cod.f=F)
    }
  }else{
    for (i in 2:(n.sim+1)){
      cont.list[[i]]  <- list(data=data.frame(dsimu(data)),cod.f=F)
    }
  }
  clust <- makeCluster(5)
  clusterExport(cl=clust, varlist=func,envir=environment())
  res <- parLapply(clust,cont.list,fun=sol.front,kcls=25)
  stopCluster(clust)
  cont.front <- do.call(rbind,res[2:(n.sim+1)])
  sol <- res[[1]]
  euclid<-crossdist(sol$df_med$conec,sol$df_med$compa,cont.front$conec,cont.front$compa)
  vet2<-apply(euclid,1,min)
  min.dist<-which(euclid == max(vet2), arr.ind = TRUE)
  pos<- lapply(vet2,function(X)which(euclid == X, arr.ind = TRUE))
  group <- sol$group[[min.dist[1]]]
  gc()
  return(list(group=group,
              sol=min.dist,
              sol.data=sol$df_med,
              dol.pop=sol$group,
              cont.front=cont.front,
              all.dis=vet2,
              all.pos=pos))
}









# plot mock ---------------------------------------------------------------

plot.mock <- function(obj.mock,n.solu=1,nome=NULL){
  r.value <- which(rank(-obj.mock$all.dis)%in%n.solu)
  p.chose <- obj.mock$all.pos[seq_along(obj.mock$all.pos)%in%r.value]
  selec <- plyr::ldply(p.chose)
  front <- selec[,1]
  cont <- selec[,2]
  best <- rbind(obj.mock$sol.data[front,],obj.mock$cont.front[cont,])
  best$grupo <- rep(1:length(n.solu),2)
  obj <- ggplot(obj.mock$cont.front,aes(y=compa,x=conec))+
    geom_point(alpha=.6,col="black")+
    geom_point(data=obj.mock$sol.data,aes(y=compa,x=conec),alpha=.8,col="red")+
    geom_line(data=best,aes(y=compa,x=conec,group=grupo))+
    theme_bw()+
    xlab("Conectividade")+ylab("Compactação")
  if(!is.null(nome)){ 
    obj<- obj+ggtitle(nome)}
  print(obj)
}


