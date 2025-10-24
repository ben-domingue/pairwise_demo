
##make pairwise comparisons for recovering theta
gen.pair<-function(th,np,nu= -1,pairsd) {
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
        ##ii<-sample(1:length(th),2) #naive approach
        ii<-sample(1:length(th),1)
        del<-th-th[ii]
        del<-ifelse(del==0,NA,del)
        p<-dnorm(del,mean=0,sd=pairsd)
        p<-ifelse(is.na(p),0,p)
        p<-rmultinom(1,1,p/sum(p))
        ii<-c(ii,which(p[,1]>0))
        pairs[[i]]<-c(ii,pr.pair(th[ii[1]],th[ii[2]],nu=nu),th[ii[1]]-th[ii[2]])
    }
    pairs<-do.call("rbind",pairs)
    pairs<-data.frame(pairs)
    names(pairs)<-c("agent_a","agent_b","winner","delta")
    pairs
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

sim.pair<-function(np,n=100,nu=0,pairsd) {
    th<-rnorm(n)
    pairs<-gen.pair(th,np,nu=nu,pairsd=pairsd)
    th.dav<-est.pair(pairs)
    e.pair<-cor(th,th.dav)
    c(n=length(th),np=np,e.pair=e.pair)
}

out<-list()
nn<-sort(runif(100,3,log10(20000)))
nn<-round(10^nn) #number of comparisons
nu<-0
for (n in c(100,250)) {
    library(parallel)
    ##
    pairs<-list()
    for (pairsd in c(.1,.5,1,2)) {
        z<-mclapply(nn,sim.pair,n=n,mc.cores=10,nu=nu,pairsd=pairsd)
        pairs[[as.character(pairsd)]]<-do.call("rbind",z)
    }
    ##
    out[[as.character(n)]]<-list(pairs=pairs)
}


cols<-colorRampPalette(c("blue", "red"))( length(out[[1]]$pairs) )
par(mfrow=c(2,1),mgp=c(2,1,0),mar=c(3,3,1.5,2),oma=rep(.5,4))
pf<-function(x,y,...) {
    m<-loess(y~x)
    yy<-predict(m,se=TRUE)
    lines(x,yy$fit,lwd=3,...)
}
for (ii in 1:length(out)) {
    plot(NULL,xlim=range(nn),ylim=c(0,1),xlab="# comparisons",ylab="cor(est,true)")
    n<-as.character(names(out)[ii])
    mtext(side=3,line=0,n)
    ##
    ##
    pairs<-out[[ii]]$pairs
    for (i in 1:length(pairs)) pf(pairs[[i]][,2],pairs[[i]][,3],col=cols[i])
}
legend("bottomright",bty='n',names(pairs),lwd=2,col=cols,title=expression(nu))
