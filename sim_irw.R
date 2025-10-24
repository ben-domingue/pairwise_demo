
##make pairwise comparisons for recovering theta
gen.pair<-function(N=100,N_pairs=10000,nu= 0) {
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
    th<-rnorm(N)
    pairs<-list()
    for (i in 1:N_pairs) {
        ii<-sample(1:length(th),2)
        pairs[[i]]<-c(ii,pr.pair(th[ii[1]],th[ii[2]],nu=nu))
    }
    pairs<-do.call("rbind",pairs)
    pairs<-data.frame(pairs)
    names(pairs)<-c("agent_a","agent_b","winner")
    tmp<-pairs$winner
    pairs$winner<-ifelse(tmp==0,'draw','agent_a')
    pairs$winner<-ifelse(tmp==2,'agent_b',pairs$winner)
    pairs
}

gen.pair()
