rm(list=ls())
# load and clean data ---------------------------------------------------------------
library(rio)
d0=import("data/STEMI_data.xlsx")
d1=d0[d0$Outcome=="cardiovascmortality",]
d2=d1[, c("Study", "year", "hr", "lci", "uci")]
d2=d2[complete.cases(d2),]
d2=d2[(d2$Study %in% c("IPD_META-ANALYSIS", "SENIOR-RITA")),]
d2$hr=as.numeric(d2$hr)
d2$lci=as.numeric(d2$lci)
d2$uci=as.numeric(d2$uci)
d2$TE=log(d2$hr)
d2$seTE=(log(d2$uci)-log(d2$lci))/3.92


# frequentist MA
library(meta)
m1=metagen(TE=d2$TE, seTE=d2$seTE, studlab = d2$Study, sm="HR")
forest(m1, leftcols = c("studlab"), prediction=T, rightcols = 	c("effect", "ci"), addrows.below.overall=4, 
       label.right = "Favours control", col.label.right = "black",
       label.left = "Favours active", col.label.left = "black", random=F)

# Bayesian 
library(R2jags)
modelFE="
model {
ipd.prec<-1/ipdseTE^2
mu~dnorm(ipdTE, ipd.prec)

seniorita.prec<-1/senioritaseTE^2
senioritaTE~dnorm(mu, seniorita.prec)
HR<-exp(mu)
}
"
model1.specFE<-textConnection(modelFE) 
data <- list(ipdseTE=d2$seTE[1], ipdTE=d2$TE[1], senioritaseTE=d2$seTE[2], senioritaTE=d2$TE[2])
jags.mFE=0
n.chains=6
jags.mFE <- jags.model(model1.specFE, data = data, n.chains =n.chains, n.adapt = 500)

params <- c("HR","mu")
closeAllConnections()
sampsFE<- coda.samples(jags.mFE, params, n.iter =10000)

library(MCMCvis)
MCMCsummary(sampsFE)
MCMCtrace(sampsFE,pdf = FALSE, params = "HR") 

all.sampsFE=sampsFE[[1]]
for(i in 2:n.chains){ all.sampsFE=rbind(all.sampsFE, sampsFE[[i]])}
mean(all.sampsFE[,1]>1)
all.sampsHRFE=data.frame(HR=all.sampsFE[,1])


### create a plot for the posterior
library(dplyr)
library(gglot2)
cutoff=1
hist.y <- density(all.sampsHRFE$HR, 
                  from = 0.6, to = 1.6) %$% 
  data.frame(x = x, y = y) %>% 
  mutate(area = x >= cutoff)

mean(all.sampsHRFE$HR>1)


hist.y$area2="Control better"
hist.y$area2[hist.y$area==F]="Active better"

CompoFE <- ggplot(data = hist.y, aes(x = x, ymin = 0, ymax = y, fill = area2)) +
  geom_ribbon() +
  geom_line(aes(y = y)) +
  geom_vline(xintercept = cutoff, color = 'red') +
  annotate(geom = 'text', x = 1.03, y = 0.2, color = 'Black', label = '67%', hjust = -0.1)+
  annotate(geom = 'text', x = 0.85, y = 0.2, color = 'Black', label = '33%', hjust = -0.1)+
  labs(fill = "")+
  labs(y= "Posterior probability density", x = "Hazard ratio")+
  labs(title = "Cardiovascular mortality")+
  theme(
    plot.title = element_text(color="black", size=18, face="bold",hjust = 0.5),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    plot.subtitle = element_text(color="black", size=16,hjust = 0.5)
  )                                             

print(CompoFE)

# create plot showing distributions from all sources and combined
Nsamps=20000
IPD_MA=exp(rnorm(Nsamps, mean = d2$TE[1],sd=d2$seTE[1] ))
SENIOR_RITA=exp(rnorm(Nsamps, mean = d2$TE[2],sd=d2$seTE[2]))
combined=exp(rnorm(Nsamps, mean = MCMCsummary(sampsFE)[2,1],sd=MCMCsummary(sampsFE)[2,2]))

library(ggplot2);library(reshape2)
combine<- data.frame(IPD_MA,SENIOR_RITA,combined)
data<- melt(combine)
ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)+ 
  labs(y= "Posterior probability density", x = "Hazard ratio")+ labs(fill = "")+
  labs(title = "Cardiovascular mortality")+
  geom_vline(xintercept = cutoff, color = 'red')+ 
  theme(
    plot.title = element_text(color="black", size=18, face="bold",hjust = 0.5),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    plot.subtitle = element_text(color="black", size=16,hjust = 0.5)
  )  


