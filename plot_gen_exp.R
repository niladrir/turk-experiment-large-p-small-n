#Setting the seed of the random numbers

set.seed(1234)

library(nullabor)
library(tourr)
library(ggplot2)
library(fpc)


####=============================================================
##Function to generate 1D lineups
####=============================================================

generate_plot_1d<-function(n=30,p, noise=1, m=20){
	x<-matrix(rnorm(p*n),ncol=p)
	if(noise==0){
x[1:15,p]<-x[1:15,p]-3
x[16:30,p]<-x[16:30,p]+3
}
colnames(x)<-paste("X",1:(p),sep="")
x<-scale(x)
x<-data.frame(x, cl=factor(c(rep(1,n/2),rep(2,n/2))))
d=1

optima <- save_history(x[,-(p+1)], guided_tour(index_f=pda_pp(cl=x[,(p+1)], lambda=0.2), d=1, max.tries=500), max_bases=1000, rescale=F)
nbases<-dim(optima)[3]

optima.global<-unclass(optima)[,,nbases]
projdata<-data.frame(x=as.matrix(x[,-(p+1)])%*%optima.global, cl=x[,(p+1)],nbases=rep(nbases,n))
projdata.true<-projdata[with(projdata, order(cl)), ] 
if(mean(projdata.true$x[1:n/2])>0)
projdata.true$cl<-c(rep(2,n/2),rep(1,n/2))

projdata.samples<-NULL
for (i in 1:(m-1)) {
 x[,(p+1)]<-sample(x[,(p+1)])
 optima <- save_history(x[,-(p+1)], guided_tour(index_f=pda_pp(cl=x[,(p+1)], lambda=0.2),d=1,max.tries=500), max_bases=1000, rescale=F)
 nbases<-dim(optima)[3]
 optima.global<-unclass(optima)[,,nbases]
 projdata<-data.frame(x=as.matrix(x[,-(p+1)])%*%optima.global, cl=x[,(p+1)], nbases=rep(nbases,n), .n=i)
 projdata.sort<-projdata[with(projdata, order(cl)), ] 
 if(mean(projdata.sort$x[1:n/2])>0)
 projdata.sort$cl<-c(rep(2,n/2),rep(1,n/2))
 projdata.samples<-rbind(projdata.samples, projdata.sort)
# cat(i,"\n")
}
y<-rep(1.5,n)
pos<-sample(m,1)
lineup.data<-lineup(true=projdata.true, samples=projdata.samples,pos=pos)
#write.table(lineup.data,"1d_data_example3.txt")
fig<-qplot(x,y, data=lineup.data, colour=factor(cl),geom="jitter") + facet_wrap(~ .sample) + scale_colour_discrete(name="Group") + scale_y_continuous("",breaks=c(0,1.5,3),limits=c(0,3))  + scale_x_continuous("", breaks=c(-1.5,0.0,1.5)) + opts(legend.position="none")
result<-c(n,"p_noise_d",paste(p,"_",noise,"_",d,sep=""),pval=NA,pos)
#ggsave("1d_example3.png",width=7, height=2.3, dpi=75) 
return(list(lineup.data=lineup.data,result=result,fig=fig))
}

#res1<-generate_plot_1d(n=30,p=40,noise=0,m=20)

####=============================================================
##Function to generate 2D lineups
####=============================================================


generate_plot_2d<-function(n=30,p, noise=1, m=20){
	x<-matrix(rnorm(p*n),ncol=p)
	if(noise==0){
x[1:10,(p-1)]<-x[1:10,(p-1)]+3
x[11:20,(p-1)]<-x[11:20,(p-1)]-3
x[21:30,p]<-x[21:30,p]+sqrt(27)
}
colnames(x)<-paste("X",1:(p),sep="")
x<-scale(x)
x<-data.frame(x, cl=factor(c(rep(1,n/3),rep(2,n/3),rep(3,n/3))))
d=2

optima <- save_history(x[,-(p+1)], tour_path=guided_tour(index_f=pda_pp(cl=x[,(p+1)], lambda=0.2), max.tries=500), max_bases=1000, rescale=F)
nbases<-dim(optima)[3]
optima.global<-unclass(optima)[,,nbases]

projdata.true<-data.frame(as.matrix(x[,-(p+1)])%*%optima.global, cl=x[,(p+1)], nbases=rep(nbases,n))

lamb.true<-summary(manova(cbind(X1, X2)~cl, data=projdata.true), test="Wilks")[[4]][3]
projdata.true<-cbind(projdata.true,wilks.lambda=lamb.true)


projdata.samples<-NULL
lamb.vals<-NULL
for (i in 1:30) {
 x[,(p+1)]<-sample(x[,(p+1)])
 optima <- save_history(x[,-(p+1)], guided_tour(index_f=pda_pp(cl=x[,(p+1)], lambda=0.2),max.tries=500), max_bases=1000, rescale=F)
 nbases<-dim(optima)[3]
 optima.global<-unclass(optima)[,,nbases]
 projdata<-data.frame(as.matrix(x[,-(p+1)])%*%optima.global, cl=x[,(p+1)], nbases=rep(nbases,30), .n=i)
lamb<-summary(manova(cbind(X1, X2)~cl, data=projdata), test="Wilks")[[4]][3]
lamb.vals<-c(lamb.vals,lamb)
 projdata.samples<-rbind(projdata.samples, projdata)

}
 wilks.lambda<-rep(lamb.vals,each=n)
projdata.samples<-cbind(projdata.samples,wilks.lambda)
projdata.new<-head(projdata.samples[order(projdata.samples$wilks.lambda),],(m-1)*30)
projdata.new$.n<-rep(1:(m-1),each=n)

pos<-sample(m,1)
lineup.data<-lineup(true=projdata.true, samples=projdata.new,pos=pos)
fig<-qplot(X1, X2, data=lineup.data, colour=cl) + facet_wrap(~ .sample) + scale_colour_discrete(name="Group") + opts(legend.position="none") + scale_y_continuous("") + scale_x_continuous("")
result<-c(n,"p_noise_d",paste(p,"_",noise,"_",d,sep=""),pval=NA,pos)
return(list(lineup.data=lineup.data,result=result,fig=fig))
}

#res2<-generate_plot_2d(n=30,p=40,noise=0,m=20)

###===============================================================
##generate lineups for large p, small n 
###===============================================================

n <- 30
dimen <- c(40,80)
pres_noise<-c(0,1)

results1 <- NULL
results2 <- NULL
difficulty <- 0
for( p in dimen){
	for( noise in pres_noise){
		if(noise==0){
	difficulty <- difficulty + 1
	}
		if(noise==1){
	difficulty <- 5
	}
	for (rep in 1){
		res1 <- generate_plot_1d(n=30,p=p,noise=noise,m=20)
		fname <- paste("large_p_small_n_",n, "_", p, "_",noise,"_1_", rep, sep="")
		dat.name <- paste("dat_", fname, ".txt", sep="")
		file.name <- paste("plot_", fname, ".png", sep="")
		results1 <- rbind(results1,c(res1$result, file.name,"large_p_small_n",difficulty))
		ggsave(paste(file.name,sep=""), plot=res1$fig, width=7.17, height=7.17, dpi=75)
		write.table(res1$lineup.data, file=paste(dat.name,sep=""), row.names =F)
		res2 <- generate_plot_2d(n=30,p=p,noise=noise,m=20)
		fname <- paste("large_p_small_n_",n, "_", p, "_",noise,"_2_", rep, sep="")
		dat.name <- paste("dat_", fname, ".txt", sep="")
		file.name <- paste("plot_", fname, ".png", sep="")
		results2 <- rbind(results2,c(res2$result, file.name,"large_p_small_n",difficulty))
		ggsave(paste(file.name,sep=""), plot=res2$fig, width=7.17, height=7.17, dpi=75)
		write.table(res2$lineup.data, file=paste(dat.name,sep=""), row.names =F)
		}
	}
}
results<-rbind(results1,results2)
pic_id<-1:dim(results)[1]
results<-data.frame(pic_id,results)

colnames(results) <- c("pic_id","sample_size","test_param","param_value","p_value","location","pic_name","experiment","difficulty")
write.table(results, file="result_table.txt",row.names =FALSE)
