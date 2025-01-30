# load and clean data ---------------------------------------------------------------
library(rio)
d0=import("data/ACS_data.xlsx")
d1=d0[d0$Outcome=="comp_allcausemortality_mi",]
d2=d1[, c("Study", "year", "hr", "lci", "uci")]
d2=d2[complete.cases(d2),]
d2[,1:5]
d2=d2[!(d2$Study %in% c("IPD_META-ANALYSIS")),]
d2$hr=as.numeric(d2$hr)
d2$lci=as.numeric(d2$lci)
d2$uci=as.numeric(d2$uci)
d2$TE=log(d2$hr)
d2$seTE=(log(d2$uci)-log(d2$lci))/3.92
d2$studyear=paste(d2$Study, " (", d2$year, ")",sep="")

### frequentist meta-analysis ----------------------------------
library(meta)
m1=metagen(TE=d2$TE, seTE=d2$seTE, studlab = d2$studyear, sm="HR")
forest(m1, leftcols = c("studlab"), prediction=T, rightcols = 	c("effect", "ci"), addrows.below.overall=4, 
       label.right = "Favours control", col.label.right = "black",
       label.left = "Favours active", col.label.left = "black",)

# Bayesian Fixed effect meta-analysis ----------------------------------
library(R2jags)
modelFE="
model {
for (i in 1: NS){
prec[i]<-1/seTE[i]^2
TE[i]~dnorm(mu, prec[i])
}
mu~dnorm(0,0.01)
HR<-exp(mu)
}
"
model1.specFE<-textConnection(modelFE) 
data <- list(NS=length(d2$Study), TE=d2$TE, seTE=d2$seTE)
jags.mFE=0
n.chains=6
jags.mFE <- jags.model(model1.specFE, data = data, n.chains =n.chains, n.adapt = 500)

params <- "HR"
closeAllConnections()
sampsFE<- coda.samples(jags.mFE, params, n.iter =10000)

library(MCMCvis)
MCMCsummary(sampsFE)
MCMCtrace(sampsFE,pdf = FALSE, params = "HR") 

all.sampsFE=sampsFE[[1]]
for(i in 2:n.chains){ all.sampsFE=rbind(all.sampsFE, sampsFE[[i]])}
mean(all.sampsFE[,1]<1)
all.sampsHRFE=data.frame(HR=all.sampsFE[,1])


# create a plot - fixed effects
library(dplyr)
cutoff=1
hist.y <- density(all.sampsHRFE$HR, from = 0.6, to = 1.4) %$% 
  data.frame(x = x, y = y) %>% 
  mutate(area = x >= cutoff)

mean(all.sampsHRFE$HR>1)


hist.y$area2="Control better"
hist.y$area2[hist.y$area==F]="Active better"

library(ggplot2)
CompoFE <- ggplot(data = hist.y, aes(x = x, ymin = 0, ymax = y, fill = area2)) +
  geom_ribbon() +
  geom_line(aes(y = y)) +
  geom_vline(xintercept = cutoff, color = 'red') +
  annotate(geom = 'text', x = 1.01, y = 0.2, color = 'Black', label = '7%', hjust = -0.1)+
  annotate(geom = 'text', x = 0.8, y = 0.2, color = 'Black', label = '93%', hjust = -0.1)+
  labs(fill = "")+labs(y= "Posterior probability density", x = "Hazard ratio")+
  labs(title = "Composite outcome all-cause mortality & MI ",
       subtitle = "fixed effects model")+
  theme(
    plot.title = element_text(color="black", size=18, face="bold",hjust = 0.5),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    plot.subtitle = element_text(color="black", size=16,hjust = 0.5)
  )                                             

print(CompoFE)

# sensitivity analysis - changing the prior 
library(R2jags)
modelFE.sens="
model {
for (i in 1: NS){
prec[i]<-1/seTE[i]^2
TE[i]~dnorm(mu, prec[i])
}
mu~dnorm(0,1)
HR<-exp(mu)
}
"
model1.specFE.sens<-textConnection(modelFE.sens) 
data <- list(NS=length(d2$Study), TE=d2$TE, seTE=d2$seTE)
jags.mFE.sens=0
n.chains=6
jags.mFE.sens <- jags.model(model1.specFE.sens, data = data, n.chains =n.chains, n.adapt = 500)

params <- "HR"
closeAllConnections()
sampsFE.sens<- coda.samples(jags.mFE.sens, params, n.iter =10000)

library(MCMCvis)
MCMCsummary(sampsFE.sens)



#### Bayesian random effects meta-analysis ----------------------------------
library(R2jags)
model="
model {
for (i in 1: NS){
prec[i]<-1/seTE[i]^2
TE[i]~dnorm(theta[i], prec[i])
theta[i]~dnorm(mu, tau.prec)
}

mu~dnorm(0,0.01)
tau.prec<-1/tau^2
tau~dnorm(0,1) T(0,)
tau.s<-tau*tau

HR<-exp(mu)
}
"
model1.spec<-textConnection(model) 
data <- list(NS=length(d2$Study), TE=d2$TE, seTE=d2$seTE)
jags.m=0
n.chains =6
jags.m <- jags.model(model1.spec, data = data, n.chains =n.chains, n.adapt = 500)

params <- c("tau", "HR", "tau.s", "mu") 
closeAllConnections()
sampsRE<- coda.samples(jags.m, params, n.iter =20000)

library(MCMCvis)
MCMCsummary(sampsRE)
MCMCtrace(sampsRE,pdf = FALSE, params = "HR") 


all.sampsRE=sampsRE[[1]]
for(i in 2:n.chains){ all.sampsRE=rbind(all.sampsRE, sampsRE[[i]])}
mean(all.sampsRE[,1]<1)

library(ggplot2)
all.sampsHRRE=data.frame(HR=all.sampsRE[,1])


### create a plot random effects
library(magrittr)
library(dplyr)
 
cutoff=1

hist.y <- density(all.sampsHRRE$HR, from = 0.6, to = 1.4) %$% 
  data.frame(x = x, y = y) %>% 
  mutate(area = x >= cutoff)

mean(all.sampsHRRE$HR<1)
mean(all.sampsHRRE$HR<0.975)
mean(all.sampsHRRE$HR<0.95)
mean(all.sampsHRRE$HR<0.90)

mean(all.sampsHRRE$HR>1)
mean(all.sampsHRRE$HR>1.025)
mean(all.sampsHRRE$HR>1.05)
mean(all.sampsHRRE$HR>1.10)


hist.y$area2="Control better"
hist.y$area2[hist.y$area==F]="Active better"

CompoRE <- ggplot(data = hist.y, aes(x = x, ymin = 0, ymax = y, fill = area2)) +
  geom_ribbon() +
  geom_line(aes(y = y)) +
  geom_vline(xintercept = cutoff, color = 'red') +
  annotate(geom = 'text', x = 1.01, y = 0.2, color = 'Black', label = '14%', hjust = -0.1)+
  annotate(geom = 'text', x = 0.8, y = 0.2, color = 'Black', label = '86%', hjust = -0.1)+
  labs(fill = "")+labs(y= "Posterior probability density", x = "Hazard ratio")+
  labs(title = "Composite outcome all-cause mortality & MI ",
       subtitle = "random effects model")+
  theme(
    plot.title = element_text(color="black", size=18, face="bold",hjust = 0.5),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    plot.subtitle = element_text(color="black", size=16,hjust = 0.5)
  )                                           


print(CompoRE)

library(gridExtra)
grid.arrange(CompoFE, CompoRE, nrow = 1)


# create plot showing distributions from all sources and combined
N.samps=20000
IPD_MA=exp(rnorm(N.samps, mean = log(as.numeric(d1$hr[1])),
                 sd=(log(as.numeric(d1$uci[1]))-log(as.numeric(d1$lci[1])))/3.92
                 ))
SENIOR_RITA=exp(rnorm(N.samps, mean = d2$TE[7],sd=d2$seTE[7]))
updated_IPD_MA=exp(rnorm(N.samps, mean = MCMCsummary(sampsRE)[2,1],sd=MCMCsummary(sampsRE)[2,2]))

library(ggplot2);library(reshape2)
combine<- data.frame(IPD_MA,SENIOR_RITA,updated_IPD_MA)
data<- melt(combine)
ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)+ 
  labs(y= "Posterior probability density", x = "Hazard ratio")+ labs(fill = "")+
  labs(title = "Composite outcome all-cause mortality & MI")+
  geom_vline(xintercept = cutoff, color = 'red')+ 
  theme(
    plot.title = element_text(color="black", size=18, face="bold",hjust = 0.5),
    axis.title.x = element_text(color="black", size=12, face="bold"),
    axis.title.y = element_text(color="black", size=12, face="bold"),
    plot.subtitle = element_text(color="black", size=16,hjust = 0.5)
  )  


# sensitivity analysis
library(R2jags)
model.RE.sens="
model {
for (i in 1: NS){
prec[i]<-1/seTE[i]^2
TE[i]~dnorm(theta[i], prec[i])
theta[i]~dnorm(mu, tau.prec)
}

mu~dnorm(0,1) # dnorm(0,0.01)
tau.prec<-1/tau^2
tau~ dt(0, pow(2.5,-2), 1)I(0,) 
tau.s<-tau*tau

HR<-exp(mu)
}
"
model1.spec.RE.sens<-textConnection(model.RE.sens) 
data <- list(NS=length(d2$Study), TE=d2$TE, seTE=d2$seTE)
jags.m=0
n.chains =6
jags.m.RE.sens <- jags.model(model1.spec.RE.sens, data = data, n.chains =n.chains, n.adapt = 500)

params <- c("tau", "HR", "tau.s", "mu") 
closeAllConnections()
sampsRE.sens<- coda.samples(jags.m.RE.sens, params, n.iter =20000)

library(MCMCvis)
MCMCsummary(sampsRE)
MCMCsummary(sampsRE.sens)
