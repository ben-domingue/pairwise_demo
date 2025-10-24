####################################################
##generate data

##make item response data for recovering theta
gen.irt<-function(b,nn) {
    th<-rnorm(nn)
    k<-outer(th,b,'-')
    p<-1/(1+exp(-k))
    resp<-p
    for (i in 1:ncol(p)) resp[,i]<-rbinom(nrow(resp),1,p[,i])
    resp<-data.frame(resp)
    resp
}

##make pairwise comparisons for recovering theta
gen.pair<-function(th,np,nu= -1) {
    pr.pair<-function(th1,th2,nu=nu) { ##based on davidson model, https://link.springer.com/article/10.3758/s13428-021-01714-2
        K<-(exp(th1)+exp(th2)+exp(nu+(th1+th2)/2))
        ##pr tie, ##see eqn 8
        pr.tie<-exp(nu+(th1+th2)/2)/K
        ##see eqn 7
        pr.12<-exp(th1)/K
        pr.21<-exp(th2)/K
        x<-rmultinom(1,1,c(pr.tie,pr.12,pr.21))[,1]
        ii<-which(x==1)
        ii-1 ##0-draw, 1 1>2, 2 2>1
    }
    pairs<-list()
    for (i in 1:np) {
        ii<-sample(1:length(th),2)
        pairs[[i]]<-c(ii,pr.pair(th[ii[1]],th[ii[2]],nu=nu),th[ii[1]]-th[ii[2]])
    }
    pairs<-do.call("rbind",pairs)
    pairs<-data.frame(pairs)
    names(pairs)<-c("agent_a","agent_b","winner","delta")
    pairs
}

####################################################
##estimate

## ##irt
## est.irt<-function(resp) {
##     m<-mirt::mirt(resp,1,'Rasch')
##     th.irt<-mirt::fscores(m)[,1]
##     th.irt
## }
est.irt<-function(resp) {
    m<-mirt::mirt(resp,1,'Rasch')
    ##th.irt<-mirt::fscores(m)[,1]
    ##th.irt
    mirt::coef(m,IRTpars=TRUE,simplify=TRUE)$item[,2]
}

##davidson model
est.pair<-function(pairs) {
    pairs$winner<-ifelse(pairs$winner==2,-1,pairs$winner)
    pairs$agent_a<-as.factor(pairs$agent_a)
    pairs$agent_b<-as.factor(pairs$agent_b)
    pairs$delta<-NULL
    library(gnm)
    library(BradleyTerry2)
    pairs.tri <- expandCategorical(pairs, "winner", idvar = "pair")
    dav <- gnm(count ~ GenDavidson(
                   winner==1,winner==0,winner==-1,
                   player1=agent_a,
                   player2=agent_b) - 1,
               eliminate = pair, family = poisson, data = pairs.tri)
    th.dav<-dav$coef[-1]
    th.dav
}
    
####################################################
## sim.irt<-function(ni,n=100) {
##     th<-rnorm(n)
##     th.irt<-est.irt(gen.irt(th,ni))
##     e.irt<-cor(th,th.irt)
##     c(n=length(th),ni=ni,e.irt=e.irt)
## }
sim.irt<-function(nn,n=100) {
    b<-rnorm(n)
    est<-est.irt(gen.irt(b,nn))
    e.irt<-cor(b,est)
    c(n=length(b),nn=nn,n.irt=e.irt)
}
sim.pair<-function(np,n=100,nu=0) {
    th<-rnorm(n)
    th.dav<-est.pair(gen.pair(th,np,nu=nu))
    e.pair<-cor(th,th.dav)
    c(n=length(th),np=np,e.pair=e.pair)
}

n<-100
out<-list()
nn<-sort(runif(50,3,log10(10000)))
nn<-round(10^nn) #number of comparisons
for (n in c(100)) {
    library(parallel)
    for (np in c(100,250,500)) irt[[as.character(np)]]<-mclapply(rep(np,10),sim.irt,n=n,mc.cores=10) #np is number of people in irt
    ##
    pairs<-list()
    for (nu in c(-1,0,1)) {
        z<-mclapply(nn,sim.pair,n=n,mc.cores=10,nu=nu)
        pairs[[as.character(nu)]]<-do.call("rbind",z)
    }
    ##
    out[[as.character(n)]]<-list(irt=irt,pairs=pairs)
}




par(mgp=c(2,1,0),mar=c(3,3,1.5,2),oma=rep(.5,4))
pf<-function(x,y,...) {
    m<-loess(y~x)
    yy<-predict(m,se=TRUE)
    lines(x,yy$fit,col='red',lwd=3,...)
    cc<-col2rgb("red")
    ## polygon(c(x,rev(x)),
    ##         c(yy$fit+1.96*yy$se.fit,rev(yy$fit-1.96*yy$se.fit)),
    ##         col=rgb(cc[1],cc[2],cc[3],max=255,alpha=155),
    ##         border=NA)
}
for (ii in 1:length(out)) {
    plot(NULL,xlim=range(nn),ylim=c(.8,1),xlab="# comparisons",ylab="cor(est,true)")
    n<-as.character(names(out)[ii])
    mtext(side=3,line=0,n)
    ##
    irt<-out[[ii]]$irt
    hl<-numeric()
    for (i in 1:length(irt)) {
        z<-irt[[i]]
        ni<-unique(sapply(z,function (x) x[2]))
        z<-sapply(z,function(x) x[3])
        hl[i]<-mean(z)
        abline(h=hl[i])
        mtext(side=4,line=0,las=2,ni,at=hl[i])
    }
    ##
    pairs<-out[[ii]]$pairs
    for (i in 1:length(pairs)) pf(pairs[[i]][,2],pairs[[i]][,3],lty=i)
    legend("bottomright",bty='n',names(pairs),lwd=2,lty=1:length(names(pairs)),col='red',title=expression(nu))
    ni<-as.numeric(names(irt))
    points(ni*as.numeric(n),hl,cex=1.5,pch=19)
}
