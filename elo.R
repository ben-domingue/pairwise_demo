##meant as an example for the webpage

df <- irw::irw_fetch("nfl_2010-2019", comp = TRUE)
df<-df[order(df$date),]

library(elo)
er<-elo.run(score(score_a, score_b) ~agent_a+agent_b , data = df, k = 20)

z<-as.data.frame(er)
L<-list()
for (i in 1:nrow(z)) L[[i]]<-data.frame(team=c(z$team.A[i],z$team.B[i]),
                                        elo=c(z$elo.A[i],z$elo.B[i]),
                                        date=df$date[i])
z<-data.frame(do.call("rbind",L))
L<-split(z,z$team)

par(mfrow=c(1,2),mgp=c(2,1,0),mar=c(3,3,1,1))
plot(NULL,xlim=c(min(z$date),1.05*max(z$date)),ylim=range(z$elo))
for (i in 1:length(L)) {
    lines(L[[i]]$date,L[[i]]$elo)
    n<-nrow(L[[i]])
    text(L[[i]]$date[n],L[[i]]$elo[n],L[[i]]$team[n],pos=4)
}

##maybe also split into seasons and then do bt analysis in each season? 
#km<-kmeans(df$date,10,iter.max=250)
#plot(df$date)
#plot(km$cluster,df$date)
del<-diff(df$date)
df$season<-c(0,cumsum(del>1e7))+2010

df<-df[df$winner!='draw',]
df$win<-ifelse(df$winner=='agent_a',1,0)
df$agent_a<-factor(df$agent_a)
df$agent_b<-factor(df$agent_b)

library(BradleyTerry2)
L<-split(df,df$season)
out<-list()
for (i in 1:length(L)) {
    mod <- BTm(win, agent_a, agent_b,
               data = L[[i]], id = "team")
    co<-coef(mod)
    names(co)<-gsub("team","",names(co))
    out[[i]]<-data.frame(season=names(L)[i],est=co,team=names(co))
}
z<-data.frame(do.call("rbind",out))
z$season<-as.numeric(z$season)

L<-split(z,z$team)
plot(NULL,xlim=c(2010,2020),ylim=c(-4,4))
for (i in 1:length(L)) {
    lines(L[[i]]$season,L[[i]]$est)
    n<-nrow(L[[i]])
    text(L[[i]]$season[n],L[[i]]$est[n],L[[i]]$team[n],pos=4)
}
