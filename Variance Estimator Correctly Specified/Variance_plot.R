###
cov_prob<-function(lo,up){
  x=mean((lo[,3]<0.5) & (up[,3]>0.5))
  return(x)
}
##
ind1<-readRDS("ind1.rds")
sr1=cov_prob(ind1[[1]],ind1[[2]])
pl1=cov_prob(ind1[[3]],ind1[[4]])

ind2<-readRDS("ind2.rds")
sr2=cov_prob(ind2[[1]],ind2[[2]])
pl2=cov_prob(ind2[[3]],ind2[[4]])

ind3<-readRDS("ind3.rds")
sr3=cov_prob(ind3[[1]],ind3[[2]])
pl3=cov_prob(ind3[[3]],ind3[[4]])

ind4<-readRDS("ind4.rds")
sr4=cov_prob(ind4[[1]],ind4[[2]])
pl4=cov_prob(ind4[[3]],ind4[[4]])

summ1<-readRDS("summ1.rds")
ps1=cov_prob(summ1[[1]],summ1[[2]])
cal1=cov_prob(summ1[[3]],summ1[[4]])

summ2<-readRDS("summ2.rds")
ps2=cov_prob(summ2[[1]],summ2[[2]])
cal2=cov_prob(summ2[[3]],summ2[[4]])

summ3<-readRDS("summ3.rds")
ps3=cov_prob(summ3[[1]],summ3[[2]])
cal3=cov_prob(summ3[[3]],summ3[[4]])

summ4<-readRDS("summ4.rds")
ps4=cov_prob(summ4[[1]],summ4[[2]])
cal4=cov_prob(summ4[[3]],summ4[[4]])

####
library(ggplot2)
library(tidyverse)
library(reshape2)

data=NULL
data$values=c(pl1,pl2,pl3,pl4,
              sr1,sr2,sr3,sr4,
              ps1,ps2,ps3,ps4,
              cal1,cal2,cal3,cal4)

data$method=c(rep("PL",times=4),rep("SR",times=4),rep("PS",times=4),rep("CL",times=4))
data$method=factor(c(rep("PL",times=4),rep("SR",times=4),rep("PS",times=4),rep("CL",times=4)),
                   levels=c("PL","SR","PS","CL"))
data$setup=rep(c(rep(c("DAG 1","DAG 2","DAG 3","DAG 4"),each=1)),times=4)
data=data.frame(data)
         

g <- ggplot(data = data, aes(x=method,y = values, fill = method)) +
  geom_bar(width=0.5,position = "dodge", color = "black", stat = "summary", fun = "mean") +
  scale_fill_manual(values = c("gray65", "deepskyblue1","slateblue2","aquamarine2")) +
  facet_wrap(~setup)+
  theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1),
        legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 10, face = "bold"),
        legend.background = element_rect(size=.5, linetype="solid", color = "black"))

g<-g + theme(legend.position = "none")+labs(y="Coverage Probability",x="Method")
g <-g + theme(axis.text = element_text(size = 16))+ 
  theme(title = element_text(size = 16))+
  theme(axis.title = element_text(size = 16)) +  theme(legend.title = element_text(size = 16)) + 
  theme(legend.text = element_text(size = 16)) + theme(strip.text.x = element_text(size = 16))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
# +
#   geom_point(data = filter(df, replicate == "B"), aes(shape = "B"), position = position_dodge(width = 0.9)) +
#   guides(fill = guide_legend(override.aes = list(shape = c(NA, NA, NA, NA))))
ggsave("cp.pdf",g,width = 10, height = 7)
ggsave()

