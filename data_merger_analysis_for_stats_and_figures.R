###one more time, we're gonna publish on citri pop-gen, I swear!
#this is the ancillary script to end em all, generate a one-stop genome mat that holds all relevant info
#we want to start with the 7mat, swap in the SNV curation data
library(data.table)
library(rcompanion)
library(FSA)
setwd("C:/Users/andre/Documents/Documents/Documents/CitriPopgen/mac ported data and code/")

#we need a modded version with the variants clipped out:
pc12<-as.data.table(read.csv("p_citri_matrix_v7_April_2024_varless.csv",stringsAsFactors=F,header=T))

#let's deal with ericerus orthologs

ortho<-as.data.frame(fread("C:/Users/andre/Documents/Documents/Documents/CitriPopgen/mac ported data and code/Ericerus/PlciErpa.proteinortho_mod.txt"),stringsAsFactors=F)
consortho<-ortho[which(ortho$Genes==2),4:5]
#5,179 1-to-1's
colnames(consortho)<-c("Ericerus_ortholog","gene")




pc13<-merge(pc12,consortho,by="gene",all.x=T)

Ericerus_ortho<-rep("None",nrow(pc13))
pc14<-cbind(pc13,Ericerus_ortho)

pc14[which(!is.na(pc14$Ericerus_ortholog)),18]<-"Ericerus"

table(pc14$Ericerus_ortho,pc14$ortho)



#and we need the curated variants:
cur<-as.data.table(read.csv("C:/Users/andre/Documents/Documents/Documents/CitriPopgen/mac ported data and code/SNP_curation/curatedSNVs_and_freqs.csv",stringsAsFactors = F))
colnames(cur)[1]<-"gene"
pc15<-merge(pc14,cur,by="gene",all.X=T)


#write.csv(pc15,"p_citri_matrix_v8_april_2024.csv",quote=F,row.names=F)
#that's where we'll start for now
pc15<-as.data.table(read.csv("p_citri_matrix_v8_April_2024_fresh.csv",stringsAsFactors=F,header=T))

#while we're here, we can precalculate some stuff
dN.dS<-(pc15$Dn/pc15$Non.sites)/(pc15$Ds/pc15$Syn.sites)
pN.pS<-(pc15$Pn/pc15$Non.sites)/(pc15$Ps/pc15$Syn.sites)
pN2.pS<-(pc15$Pn.2/pc15$Non.sites)/(pc15$Ps/pc15$Syn.sites)
pN3.pS<-(pc15$Pn.3/pc15$Non.sites)/(pc15$Ps/pc15$Syn.sites)
pN4.pS<-(pc15$Pn.4/pc15$Non.sites)/(pc15$Ps/pc15$Syn.sites)
pN5.pS<-(pc15$Pn.5/pc15$Non.sites)/(pc15$Ps/pc15$Syn.sites)
pN6.pS<-(pc15$Pn.6/pc15$Non.sites)/(pc15$Ps/pc15$Syn.sites)
pN7.pS<-(pc15$Pn.7/pc15$Non.sites)/(pc15$Ps/pc15$Syn.sites)
pN8.pS<-(pc15$Pn.8/pc15$Non.sites)/(pc15$Ps/pc15$Syn.sites)
pN9.pS<-(pc15$Pn.9/pc15$Non.sites)/(pc15$Ps/pc15$Syn.sites)
alpha.simple<-1-(pN.pS/dN.dS)
alpha.impmkt2<-1-(pN2.pS/dN.dS)
alpha.impmkt3<-1-(pN3.pS/dN.dS)
alpha.impmkt4<-1-(pN4.pS/dN.dS)
alpha.impmkt5<-1-(pN5.pS/dN.dS)
alpha.impmkt6<-1-(pN6.pS/dN.dS)
alpha.impmkt7<-1-(pN7.pS/dN.dS)
alpha.impmkt8<-1-(pN8.pS/dN.dS)
alpha.impmkt9<-1-(pN9.pS/dN.dS)

pc16<-cbind(pc15,dN.dS,pN.pS,pN2.pS,pN3.pS,pN4.pS,pN5.pS,pN6.pS,pN7.pS,pN8.pS,pN9.pS,alpha.simple,alpha.impmkt2,alpha.impmkt3,
            alpha.impmkt4,alpha.impmkt5,alpha.impmkt9)

#table(pc16$sex.bias.adult,pc16$sex.bias.nymph)

#              female_biased male_biased  unbiased
#female_biased          1164         189      2050
#male_biased             143        1003      2097
#unbiased                453         343      5155




kruskal.test((pc16$Ps/pc16$Syn.sites)~pc16$consistency)
#Kruskal-Wallis chi-squared = 234.57, df = 2, p-value < 2.2e-16
dt<-dunnTest((pc17$Pn/pc17$Non.sites)~pc17$consistency)
dtr<-dt$res


##let's simplify things for the main manuscript


pc17<-pc16[which(pc16$consistency=="female_biased"|pc16$consistency=="unbiased"|pc16$consistency=="male_biased"),]
pc17$consistency<-factor(pc17$consistency,levels=c("female_biased","unbiased","male_biased"))

table(pc17$consistency)
#female_biased      unbiased   male_biased 
#         1164          5155          1003 


####what do we ultimately want to show?
#1. Pn, Ps, pN/pS, (and how pN drops out at higher freqs)
#1.1 pN
kruskal.test((pc17$Pn/pc17$Non.sites)~pc17$consistency)
#Kruskal-Wallis chi-squared = 234.57, df = 2, p-value < 2.2e-16
dt<-dunnTest((pc17$Pn/pc17$Non.sites)~pc17$consistency)
dtr<-dt$res
dt
#                  Comparison         Z      P.unadj        P.adj
# female_biased - male_biased  9.247833 2.290904e-20 4.581807e-20
#    female_biased - unbiased 15.314270 6.139756e-53 1.841927e-52
#      male_biased - unbiased  2.855575 4.295902e-03 4.295902e-03

cldList(comparison = dtr$Comparison, p.value    = dtr$P.adj, threshold  = 0.05)
          #Group Letter MonoLetter
# female_biased      a        a  
#   male_biased      b         b 
#      unbiased      c          c
#1.2 pS
kruskal.test((pc17$Ps/pc17$Syn.sites)~pc17$consistency)
#Kruskal-Wallis chi-squared = 44.735, df = 2, p-value = 1.932e-10
dt<-dunnTest((pc17$Ps/pc17$Syn.sites)~pc17$consistency)
dtr<-dt$res
dt
#                  Comparison        Z      P.unadj        P.adj
# female_biased - male_biased 3.117291 1.825215e-03 3.650430e-03
#    female_biased - unbiased 6.584042 4.578282e-11 1.373485e-10
#      male_biased - unbiased 2.299575 2.147233e-02 2.147233e-02

cldList(comparison = dtr$Comparison, p.value    = dtr$P.adj, threshold  = 0.05)
#         Group Letter MonoLetter
# female_biased      a        a  
#   male_biased      b         b 
#      unbiased      c          c

#why does pS differ?
nrow(pc17[which(pc17$consistency=="female_biased"&pc17$Ps==0),])/nrow(pc17[which(pc17$consistency=="female_biased"),])
#49.1%
nrow(pc17[which(pc17$consistency=="unbiased"&pc17$Ps==0),])/nrow(pc17[which(pc17$consistency=="unbiased"),])
#58.3%
nrow(pc17[which(pc17$consistency=="male_biased"&pc17$Ps==0),])/nrow(pc17[which(pc17$consistency=="male_biased"),])
#54.3%

 
#1.3 pN/pS
kruskal.test((pc17$pN.pS)~pc17$consistency)
#Kruskal-Wallis chi-squared = 148.43, df = 2, p-value < 2.2e-16
dt<-dunnTest((pc17$pN.pS)~pc17$consistency)
dtr<-dt$res
dt
#                  Comparison         Z      P.unadj        P.adj
# female_biased - male_biased  7.662233 1.827278e-14 3.654556e-14
#    female_biased - unbiased 12.157702 5.220669e-34 1.566201e-33
#      male_biased - unbiased  1.585706 1.128060e-01 1.128060e-01
cldList(comparison = dtr$Comparison, p.value    = dtr$P.adj, threshold  = 0.05)
#Group Letter MonoLetter
# female_biased      a         a 
#   male_biased      b          b
#      unbiased      b          b



#1.4 the figure
layout( matrix(c(1,2,3,3), nrow=2, byrow=TRUE) )
#par(mai=c(bottom,left,top.right))
par(mai=c(0.65,0.9,0.7,0.06))
boxplot((pc17$Pn/pc17$Non.sites)~pc17$consistency,outline=F,notch=T,ylab="",xlab="Sex-bias in expression",las=1,
        col=c("firebrick3","grey33","dodgerblue3"),cex.lab=1.2,cex.axis=1.2,
        cex.main=1.6,main="",names=c("F","UB","M"),ylim=c(0,0.014))
mtext("Variants per site", side = 2, line = 4)
mtext("pN", side = 3, line = 1, cex = 2)
text(c("a","c","b"),x=c(1,2,3),y=c(0.0014,0.0006,0.0008),cex=1.5,col="white")
par(mai=c(0.65,0.9,0.7,0.06))
boxplot((pc17$Ps/pc17$Syn.sites)~pc17$consistency,outline=F,notch=T,ylab="",xlab="Sex-bias in expression",las=1,
        col=c("firebrick3","grey33","dodgerblue3"),cex.lab=1.2,cex.axis=1.2,
        cex.main=1.6,main="",names=c("F","UB","M"),ylim=c(0,0.09))
mtext("pS", side = 3, line = 1, cex = 2)
text(c("a","c","b"),x=c(1,2,3),y=c(0.012,0.005,0.008),cex=1.7,col="white")
par(mai=c(0.7,0.9,0.6,0.06))
boxplot((pc17$pN.pS)~pc17$consistency,outline=F,notch=T,ylab="",xlab="Sex-bias in expression",las=1,
        col=c("firebrick3","grey33","dodgerblue3"),cex.lab=1.2,cex.axis=1.2,
        cex.main=1.6,main="",names=c("F","UB","M"),ylim=c(0,1.38))
text(c("a","b","b"),x=c(1,2,3),y=c(0.29,0.15,0.15),cex=1.7,col="white")
mtext("pN/pS", side = 3, line = 1, cex = 2)
mtext("Scaled polymorphism rate", side = 2, line = 4)

#double checking that factor levels are correct
#abline(h=median(pc17f$pN.pS,na.rm=T),col="firebrick3")
#abline(h=median(pc17u$pN.pS,na.rm=T),col="grey33")
#abline(h=median(pc17m$pN.pS,na.rm=T),col="dodgerblue3")
#yup they are now

#1.5 Evidence of how pN.x/pS decreases by gene class..
pc17f<-pc17[which(pc17$consistency=="female_biased"),]
#1164 genes
pc17u<-pc17[which(pc17$consistency=="unbiased"),]
#5155 genes
pc17m<-pc17[which(pc17$consistency=="male_biased"),]
#1003 genes
par(mfrow=c(1,1))
plot(x=1,y=1)
par(mfrow=c(2,2))
#female-biased
boxplot(pc17u$Pn.1/pc17u$Non.sites,notch=T, col="firebrick3",outline=F,xlim=c(1,10),ylim=c(0,0.0015),las=1)
boxplot(pc17u$Pn.2/pc17u$Non.sites,notch=T, col="firebrick3",outline=F,add=T,at=2,las=1)
boxplot(pc17u$Pn.3/pc17u$Non.sites,notch=T, col="firebrick3",outline=F,add=T,at=3,las=1)
boxplot(pc17u$Pn.4/pc17u$Non.sites,notch=T, col="firebrick3",outline=F,add=T,at=4,las=1)
boxplot(pc17u$Pn.5/pc17u$Non.sites,notch=T, col="firebrick3",outline=F,add=T,at=5,las=1)
boxplot(pc17u$Pn.6/pc17u$Non.sites,notch=T, col="firebrick3",outline=F,add=T,at=6,las=1)
boxplot(pc17u$Pn.7/pc17u$Non.sites,notch=T, col="firebrick3",outline=F,add=T,at=7,las=1)
boxplot(pc17u$Pn.8/pc17u$Non.sites,notch=T, col="firebrick3",outline=F,add=T,at=8,las=1)
boxplot(pc17u$Pn.9/pc17u$Non.sites,notch=T, col="firebrick3",outline=F,add=T,at=9,las=1)
mtext("Variants per site",side=2, line=3.5, cex=1.5)
mtext("Female-biased",side=3,line=1.2, cex=1.8)
mtext(text="0.1",side=1,adj=0.0,las=0,cex=1.5,line=1)
mtext(text="0.2",side=1,adj=0.11,las=0,cex=1.5,line=1)
mtext(text="0.3",side=1,adj=0.22,las=0,cex=1.5,line=1)
mtext(text="0.4",side=1,adj=0.33,las=0,cex=1.5,line=1)
mtext(text="0.5",side=1,adj=0.44,las=0,cex=1.5,line=1)
mtext(text="0.6",side=1,adj=0.55,las=0,cex=1.5,line=1)
mtext(text="0.7",side=1,adj=0.66,las=0,cex=1.5,line=1)
mtext(text="0.8",side=1,adj=0.77,las=0,cex=1.5,line=1)
mtext(text="0.9",side=1,adj=0.88,las=0,cex=1.5,line=1)
mtext("Pn frequency >X",side=1,adj=0.5,line=3.1,cex=1.6)

#unbiased
boxplot(pc17u$Pn.1/pc17u$Non.sites,notch=T, col="grey33",outline=F,xlim=c(1,10),las=1,ylim=c(0,0.0015))
boxplot(pc17u$Pn.2/pc17u$Non.sites,notch=T, col="grey33",outline=F,add=T,at=2,las=1)
boxplot(pc17u$Pn.3/pc17u$Non.sites,notch=T, col="grey33",outline=F,add=T,at=3,las=1)
boxplot(pc17u$Pn.4/pc17u$Non.sites,notch=T, col="grey33",outline=F,add=T,at=4,las=1)
boxplot(pc17u$Pn.5/pc17u$Non.sites,notch=T, col="grey33",outline=F,add=T,at=5,las=1)
boxplot(pc17u$Pn.6/pc17u$Non.sites,notch=T, col="grey33",outline=F,add=T,at=6,las=1)
boxplot(pc17u$Pn.7/pc17u$Non.sites,notch=T, col="grey33",outline=F,add=T,at=7,las=1)
boxplot(pc17u$Pn.8/pc17u$Non.sites,notch=T, col="grey33",outline=F,add=T,at=8,las=1)
boxplot(pc17u$Pn.9/pc17u$Non.sites,notch=T, col="grey33",outline=F,add=T,at=9,las=1)
mtext("Variants per site",side=2, line=3.5, cex=1.5)
mtext("Unbiased",side=3,line=1.2, cex=1.8)
mtext(text="0.1",side=1,adj=0.0,las=0,cex=1.5,line=1)
mtext(text="0.2",side=1,adj=0.11,las=0,cex=1.5,line=1)
mtext(text="0.3",side=1,adj=0.22,las=0,cex=1.5,line=1)
mtext(text="0.4",side=1,adj=0.33,las=0,cex=1.5,line=1)
mtext(text="0.5",side=1,adj=0.44,las=0,cex=1.5,line=1)
mtext(text="0.6",side=1,adj=0.55,las=0,cex=1.5,line=1)
mtext(text="0.7",side=1,adj=0.66,las=0,cex=1.5,line=1)
mtext(text="0.8",side=1,adj=0.77,las=0,cex=1.5,line=1)
mtext(text="0.9",side=1,adj=0.88,las=0,cex=1.5,line=1)
mtext("Pn frequency >X",side=1,adj=0.5,line=3.1,cex=1.6)


#male-biased
boxplot(pc17m$Pn.1/pc17m$Non.sites,notch=T, col="dodgerblue3",outline=F,xlim=c(1,10),las=1,ylim=c(0,0.0015))
boxplot(pc17m$Pn.2/pc17m$Non.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=2,las=1)
boxplot(pc17m$Pn.3/pc17m$Non.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=3,las=1)
boxplot(pc17m$Pn.4/pc17m$Non.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=4,las=1)
boxplot(pc17m$Pn.5/pc17m$Non.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=5,las=1)
boxplot(pc17m$Pn.6/pc17m$Non.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=6,las=1)
boxplot(pc17m$Pn.7/pc17m$Non.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=7,las=1)
boxplot(pc17m$Pn.8/pc17m$Non.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=8,las=1)
boxplot(pc17m$Pn.9/pc17m$Non.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=9,las=1)
mtext("Variants per site",side=2, line=3.5, cex=1.5)
mtext("Male-biased",side=3,line=1.2, cex=1.8)
mtext(text="0.1",side=1,adj=0.0,las=0,cex=1.5,line=1)
mtext(text="0.2",side=1,adj=0.11,las=0,cex=1.5,line=1)
mtext(text="0.3",side=1,adj=0.22,las=0,cex=1.5,line=1)
mtext(text="0.4",side=1,adj=0.33,las=0,cex=1.5,line=1)
mtext(text="0.5",side=1,adj=0.44,las=0,cex=1.5,line=1)
mtext(text="0.6",side=1,adj=0.55,las=0,cex=1.5,line=1)
mtext(text="0.7",side=1,adj=0.66,las=0,cex=1.5,line=1)
mtext(text="0.8",side=1,adj=0.77,las=0,cex=1.5,line=1)
mtext(text="0.9",side=1,adj=.88,las=0,cex=1.5,line=1)
mtext("Pn frequency >X",side=1,adj=0.5,line=3.1,cex=1.6)

#
#let's check Ps as well
#

par(mfrow=c(1,1))
plot(x=1,y=1)
par(mfrow=c(2,2))
#female-biased
boxplot(pc17u$Ps.1/pc17u$Syn.sites,notch=T, col="firebrick3",outline=F,xlim=c(1,10),ylim=c(0,0.015),las=1)
boxplot(pc17u$Ps.2/pc17u$Syn.sites,notch=T, col="firebrick3",outline=F,add=T,at=2,las=1)
boxplot(pc17u$Ps.3/pc17u$Syn.sites,notch=T, col="firebrick3",outline=F,add=T,at=3,las=1)
boxplot(pc17u$Ps.4/pc17u$Syn.sites,notch=T, col="firebrick3",outline=F,add=T,at=4,las=1)
boxplot(pc17u$Ps.5/pc17u$Syn.sites,notch=T, col="firebrick3",outline=F,add=T,at=5,las=1)
boxplot(pc17u$Ps.6/pc17u$Syn.sites,notch=T, col="firebrick3",outline=F,add=T,at=6,las=1)
boxplot(pc17u$Ps.7/pc17u$Syn.sites,notch=T, col="firebrick3",outline=F,add=T,at=7,las=1)
boxplot(pc17u$Ps.8/pc17u$Syn.sites,notch=T, col="firebrick3",outline=F,add=T,at=8,las=1)
boxplot(pc17u$Ps.9/pc17u$Syn.sites,notch=T, col="firebrick3",outline=F,add=T,at=9,las=1)
mtext("Variants per site",side=2, line=3.5, cex=1.5)
mtext("Female-biased",side=3,line=1.2, cex=1.8)
mtext(text="0.1",side=1,adj=0.0,las=0,cex=1.5,line=1)
mtext(text="0.2",side=1,adj=0.11,las=0,cex=1.5,line=1)
mtext(text="0.3",side=1,adj=0.22,las=0,cex=1.5,line=1)
mtext(text="0.4",side=1,adj=0.33,las=0,cex=1.5,line=1)
mtext(text="0.5",side=1,adj=0.44,las=0,cex=1.5,line=1)
mtext(text="0.6",side=1,adj=0.55,las=0,cex=1.5,line=1)
mtext(text="0.7",side=1,adj=0.66,las=0,cex=1.5,line=1)
mtext(text="0.8",side=1,adj=0.77,las=0,cex=1.5,line=1)
mtext(text="0.9",side=1,adj=0.88,las=0,cex=1.5,line=1)
mtext("Ps frequency >X",side=1,adj=0.5,line=3.1,cex=1.6)

#unbiased
boxplot(pc17u$Ps.1/pc17u$Syn.sites,notch=T, col="grey33",outline=F,xlim=c(1,10),las=1,ylim=c(0,0.015))
boxplot(pc17u$Ps.2/pc17u$Syn.sites,notch=T, col="grey33",outline=F,add=T,at=2,las=1)
boxplot(pc17u$Ps.3/pc17u$Syn.sites,notch=T, col="grey33",outline=F,add=T,at=3,las=1)
boxplot(pc17u$Ps.4/pc17u$Syn.sites,notch=T, col="grey33",outline=F,add=T,at=4,las=1)
boxplot(pc17u$Ps.5/pc17u$Syn.sites,notch=T, col="grey33",outline=F,add=T,at=5,las=1)
boxplot(pc17u$Ps.6/pc17u$Syn.sites,notch=T, col="grey33",outline=F,add=T,at=6,las=1)
boxplot(pc17u$Ps.7/pc17u$Syn.sites,notch=T, col="grey33",outline=F,add=T,at=7,las=1)
boxplot(pc17u$Ps.8/pc17u$Syn.sites,notch=T, col="grey33",outline=F,add=T,at=8,las=1)
boxplot(pc17u$Ps.9/pc17u$Syn.sites,notch=T, col="grey33",outline=F,add=T,at=9,las=1)
mtext("Variants per site",side=2, line=3.5, cex=1.5)
mtext("Unbiased",side=3,line=1.2, cex=1.8)
mtext(text="0.1",side=1,adj=0.0,las=0,cex=1.5,line=1)
mtext(text="0.2",side=1,adj=0.11,las=0,cex=1.5,line=1)
mtext(text="0.3",side=1,adj=0.22,las=0,cex=1.5,line=1)
mtext(text="0.4",side=1,adj=0.33,las=0,cex=1.5,line=1)
mtext(text="0.5",side=1,adj=0.44,las=0,cex=1.5,line=1)
mtext(text="0.6",side=1,adj=0.55,las=0,cex=1.5,line=1)
mtext(text="0.7",side=1,adj=0.66,las=0,cex=1.5,line=1)
mtext(text="0.8",side=1,adj=0.77,las=0,cex=1.5,line=1)
mtext(text="0.9",side=1,adj=0.88,las=0,cex=1.5,line=1)
mtext("Ps frequency >X",side=1,adj=0.5,line=3.1,cex=1.6)


#male-biased
boxplot(pc17m$Ps.1/pc17m$Syn.sites,notch=T, col="dodgerblue3",outline=F,xlim=c(1,10),las=1,ylim=c(0,0.015))
boxplot(pc17m$Ps.2/pc17m$Syn.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=2,las=1)
boxplot(pc17m$Ps.3/pc17m$Syn.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=3,las=1)
boxplot(pc17m$Ps.4/pc17m$Syn.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=4,las=1)
boxplot(pc17m$Ps.5/pc17m$Syn.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=5,las=1)
boxplot(pc17m$Ps.6/pc17m$Syn.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=6,las=1)
boxplot(pc17m$Ps.7/pc17m$Syn.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=7,las=1)
boxplot(pc17m$Ps.8/pc17m$Syn.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=8,las=1)
boxplot(pc17m$Ps.9/pc17m$Syn.sites,notch=T, col="dodgerblue3",outline=F,add=T,at=9,las=1)
mtext("Variants per site",side=2, line=3.5, cex=1.5)
mtext("Male-biased",side=3,line=1.2, cex=1.8)
mtext(text="0.1",side=1,adj=0.0,las=0,cex=1.5,line=1)
mtext(text="0.2",side=1,adj=0.11,las=0,cex=1.5,line=1)
mtext(text="0.3",side=1,adj=0.22,las=0,cex=1.5,line=1)
mtext(text="0.4",side=1,adj=0.33,las=0,cex=1.5,line=1)
mtext(text="0.5",side=1,adj=0.44,las=0,cex=1.5,line=1)
mtext(text="0.6",side=1,adj=0.55,las=0,cex=1.5,line=1)
mtext(text="0.7",side=1,adj=0.66,las=0,cex=1.5,line=1)
mtext(text="0.8",side=1,adj=0.77,las=0,cex=1.5,line=1)
mtext(text="0.9",side=1,adj=.88,las=0,cex=1.5,line=1)
mtext("Ps frequency >X",side=1,adj=0.5,line=3.1,cex=1.6)




#2. dN, dS, dN/dS
#2.1 dN
kruskal.test((pc17$Dn/pc17$Non.sites)~pc17$consistency)
#Kruskal-Wallis chi-squared = 59.751, df = 2, p-value = 1.06e-13
dt<-dunnTest((pc17$Dn/pc17$Non.sites)~pc17$consistency)
dtr<-dt$res
dt
#                  Comparison         Z      P.unadj        P.adj
# female_biased - male_biased  7.040967 1.909108e-12 5.727323e-12
#    female_biased - unbiased  6.843241 7.742126e-12 1.548425e-11
#      male_biased - unbiased -2.354940 1.852571e-02 1.852571e-02

cldList(comparison = dtr$Comparison, p.value    = dtr$P.adj, threshold  = 0.05)
#Group Letter MonoLetter
# female_biased      a        a  
#   male_biased      b         b 
#      unbiased      c          c

#2.2 dS
kruskal.test((pc17$Ds/pc17$Syn.sites)~pc17$consistency)
#Kruskal-Wallis chi-squared = 153.85, df = 2, p-value < 2.2e-16
dt<-dunnTest((pc17$Ds/pc17$Syn.sites)~pc17$consistency)
dtr<-dt$res
dt
#                 Comparison          Z      P.unadj        P.adj
# female_biased - male_biased  -1.722877 8.491082e-02 8.491082e-02
#    female_biased - unbiased -10.749506 5.958003e-27 1.787401e-26
#      male_biased - unbiased  -7.957229 1.759356e-15 3.518712e-15

cldList(comparison = dtr$Comparison, p.value    = dtr$P.adj, threshold  = 0.05)
# female_biased      a         a 
#   male_biased      a         a 
#      unbiased      b          b


#is this due to zeroes?
nrow(pc17[which(pc17$consistency=="female_biased"&pc17$Ds==0),])/nrow(pc17[which(pc17$consistency=="female_biased"),])
#17.7%
nrow(pc17[which(pc17$consistency=="unbiased"&pc17$Ds==0),])/nrow(pc17[which(pc17$consistency=="unbiased"),])
#4.1%
nrow(pc17[which(pc17$consistency=="male_biased"&pc17$Ds==0),])/nrow(pc17[which(pc17$consistency=="male_biased"),])
#8.1%

#2.3 dN/dS
kruskal.test((pc17$dN.dS)~pc17$consistency)
#Kruskal-Wallis chi-squared = 509.11, df = 2, p-value < 2.2e-16
dt<-dunnTest((pc17$dN.dS)~pc17$consistency)
dtr<-dt$res
dt
#                 Comparison         Z       P.unadj         P.adj
# female_biased - male_biased 13.583396  5.024180e-42  1.004836e-41
#    female_biased - unbiased 22.547158 1.431899e-112 4.295696e-112
#      male_biased - unbiased  4.552627  5.298027e-06  5.298027e-06
> 
cldList(comparison = dtr$Comparison, p.value    = dtr$P.adj, threshold  = 0.05)
#          Group Letter MonoLetter
# female_biased      a        a  
#   male_biased      b         b 
#      unbiased      c          c

#2.4 the Figure 
layout( matrix(c(1,2,3,3), nrow=2, byrow=TRUE) )
#par(mai=c(bottom,left,top.right))
par(mai=c(0.65,0.9,0.7,0.06))
boxplot((pc17$Dn/pc17$Non.sites)~pc17$consistency,outline=F,notch=T,ylab="",xlab="Sex-bias in expression",las=1,
        col=c("firebrick3","grey33","dodgerblue3"),cex.lab=1.2,cex.axis=1.2,
        cex.main=1.6,main="",names=c("F","UB","M"),ylim=c(0,0.024))
mtext("Variants per site", side = 2, line = 4)
mtext("dN", side = 3, line = 1, cex = 2)
text(c("a","b","c"),x=c(1,2,3),y=c(0.008,0.007,0.006),cex=1.7,col="white")
par(mai=c(0.65,0.9,0.7,0.06))
boxplot((pc17$Ds/pc17$Syn.sites)~pc17$consistency,outline=F,notch=T,ylab="",xlab="Sex-bias in expression",las=1,
        col=c("firebrick3","grey33","dodgerblue3"),cex.lab=1.2,cex.axis=1.2,
        cex.main=1.6,main="",names=c("F","UB","M"),ylim=c(0,0.19))
mtext("dS", side = 3, line = 1, cex = 2)
text(c("b","a","b"),x=c(1,2,3),y=c(0.042,0.062,0.042),cex=1.7,col="white")
par(mai=c(0.7,0.9,0.6,0.06))
boxplot((pc17$dN.dS)~pc17$consistency,outline=F,notch=T,ylab="",xlab="Sex-bias in expression",las=1,
        col=c("firebrick3","grey33","dodgerblue3"),cex.lab=1.2,cex.axis=1.2,
        cex.main=1.6,main="",names=c("F","UB","M"),ylim=c(0,0.8))
text(c("a","c","b"),x=c(1,2,3),y=c(0.25,0.12,0.15),cex=1.7,col="white")
mtext("dN/dS", side = 3, line = 1, cex = 2)
mtext("Scaled divergence rate", side = 2, line = 4)

#double checking that factor levels are correct
#abline(h=median(pc17f$dN.dS,na.rm=T),col="firebrick3")
#abline(h=median(pc17u$dN.dS,na.rm=T),col="grey33")
#abline(h=median(pc17m$dN.dS,na.rm=T),col="dodgerblue3")
#yup they are now

# alpha, imp
#3.1 simple alpha
kruskal.test((pc17$alpha.simple)~pc17$consistency)
#Kruskal-Wallis chi-squared = 7.9234, df = 2, p-value = 0.01903
dt<-dunnTest((pc17$alpha.simple)~pc17$consistency)
dtr<-dt$res
dt
#                 Comparison           Z     P.unadj      P.adj
# female_biased - male_biased -1.99002118 0.046588602 0.09317720
#    female_biased - unbiased -2.78269092 0.005391014 0.01617304
#      male_biased - unbiased -0.08755435 0.930230884 0.93023088

cldList(comparison = dtr$Comparison, p.value    = dtr$P.adj, threshold  = 0.05)
#         Group Letter MonoLetter
# female_biased      a         a 
#   male_biased     ab         ab
#      unbiased      b          b
median(pc17f$alpha.simple,na.rm=T)
#0.1264535
> median(pc17u$alpha.simple,na.rm=T)
#0.3033991
> median(pc17m$alpha.simple,na.rm=T)
#0.2142857
#here it's F < M < U but not strongly significant

#3.2 what if we use the impMKT Pn>0.2?

kruskal.test((pc17$alpha.impmkt2)~pc17$consistency)
#Kruskal-Wallis chi-squared = 21.54, df = 2, p-value = 2.102e-05
dt<-dunnTest((pc17$alpha.impmkt2)~pc17$consistency)
dtr<-dt$res
dt
#                 Comparison          Z      P.unadj        P.adj
# female_biased - male_biased -3.1231775 1.789098e-03 3.578196e-03
#    female_biased - unbiased -4.6156485 3.918697e-06 1.175609e-05
#      male_biased - unbiased -0.3776528 7.056885e-01 7.056885e-01

cldList(comparison = dtr$Comparison, p.value    = dtr$P.adj, threshold  = 0.05)
#              Group Letter MonoLetter
# female_biased      a         a 
#   male_biased      b          b
#      unbiased      b          b
median(pc17f$alpha.impmkt2,na.rm=T)
#0.4222222
median(pc17u$alpha.impmkt2,na.rm=T)
#0.6781609
median(pc17m$alpha.impmkt2,na.rm=T)
#0.7083333

#3.3 but male-biased have a lot of polys up to 0.3, let's do a stricter cutoff >0.4
kruskal.test((pc17$alpha.impmkt4)~pc17$consistency)
#Kruskal-Wallis chi-squared = 26.814, df = 2, p-value = 1.505e-06
dt<-dunnTest((pc17$alpha.impmkt4)~pc17$consistency)
dtr<-dt$res
dt
#                  Comparison          Z      P.unadj        P.adj
# female_biased - male_biased -3.4054880 6.604592e-04 1.320918e-03
#    female_biased - unbiased -5.1587884 2.485531e-07 7.456592e-07
#      male_biased - unbiased -0.5185853 6.040500e-01 6.040500e-01

cldList(comparison = dtr$Comparison, p.value    = dtr$P.adj, threshold  = 0.05)
#              Group Letter MonoLetter
# female_biased      a         a 
#   male_biased      b          b
#      unbiased      b          b
median(pc17f$alpha.impmkt4,na.rm=T)
#0.6666667
median(pc17u$alpha.impmkt4,na.rm=T)
#0.8815236
median(pc17m$alpha.impmkt4,na.rm=T)
#0.8194254

#3.3 the plot
par(mfrow=c(1,2))
boxplot(pc17$alpha.simple~pc17$consistency,outline=F,notch=T,ylab="",xlab="",las=1,
        col=c("firebrick3","grey33","dodgerblue3"),cex.lab=1.2,cex.axis=1.2,
        cex.main=1.6,main="",names=c("F","UB","M"),ylim=c(-4.2,1))
abline(h=0)
mtext(expression(alpha), side = 2, line = 3, cex = 3)
mtext("Considering all polymorphisms", side = 3, line = 1.5, cex = 1.8)
mtext("Sex bias in expression",side=1, line=2,adj=0.5,cex=1.5)
text(c("b","a","ab"),x=c(1,2,3),y=c(0.25,0.45,0.35),cex=1.7,col="white")
boxplot(pc17$alpha.impmkt2~pc17$consistency,outline=F,notch=T,ylab="",xlab="",las=1,
        col=c("firebrick3","grey33","dodgerblue3"),cex.lab=1.2,cex.axis=1.2,
        cex.main=1.6,main="",names=c("F","UB","M"),ylim=c(-1,1))
abline(h=0)
mtext("Sex bias in expression",side=1, line=2,adj=0.5,cex=1.5)
text(c("b","a","a"),x=c(1,2,3),y=c(0.5,0.75,0.75),cex=1.7,col="white")
mtext("Considering Pn > 0.2", side = 3, line = 1.5, cex = 1.8)

#for the supplement
par(mfrow=c(1,1))
boxplot(pc17$alpha.impmkt4~pc17$consistency,outline=F,notch=T,ylab="",xlab="",las=1,
        col=c("firebrick3","grey33","dodgerblue3"),cex.lab=1.2,cex.axis=1.2,
        cex.main=1.6,main="",names=c("F","UB","M"),ylim=c(-1,1))
abline(h=0)
mtext("Sex bias in expression",side=1, line=2,adj=0.5,cex=1.5)
text(c("b","a","a"),x=c(1,2,3),y=c(0.75,0.95,0.95),cex=1.7,col="white")
mtext("Considering Pn > 0.4", side = 3, line = 1.5, cex = 1.8)


#orthology
erico<-table(pc17$Ericerus_ortho,pc17$consistency)
#         female_biased unbiased male_biased
#Ericerus            67     2542         259
#None              1097     2613         744
chisq.test(erico)
#X-squared = 842.94, df = 2, p-value < 2.2e-16


#4 Hollie's biparental
###Hollie's updated biparentals
hol<-as.data.frame(fread("andres_imp_genes_output_update_march2024.txt"),stringsAsFactors=F)
colnames(hol)[1]<-"gene"
colnames(pc17)[11]<-"dN.dS.dupe"
pc18<-merge(pc17,hol,by="gene",all.x=T)
plot(pc18$maternal_proportion,pc18$alpha.impmkt2,ylim=c(-1,1))
is.na(pc18$maternal_proportion)<-1

pc18um<-pc18[which(pc18$consistency!="female_biased"),]
pc18u<-pc18[which(pc18$consistency=="unbiased"),]
pc18m<-pc18[which(pc18$consistency=="male_biased"),]
pc18m<-pc18m[is.finite(pc18m$alpha.impmkt2),]
pc18u<-pc18u[is.finite(pc18u$alpha.impmkt2),]
pc18um<-pc18um[is.finite(pc18um$alpha.impmkt2),]
plot(pc18m$maternal_proportion,pc18m$alpha.impmkt2,ylim=c(-1,1))
test<-lm(pc18m$alpha.impmkt2~pc18m$maternal_proportion,na.action=na.omit)
summary(test)
#                          Estimate Std. Error t value Pr(>|t|)  
#(Intercept)                  3.401      1.427   2.384   0.0410 *
 # pc18m$maternal_proportion   -3.290      1.662  -1.979   0.0791 .
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Residual standard error: 0.6094 on 9 degrees of freedom
#(338 observations deleted due to missingness)
#Multiple R-squared:  0.3033,	Adjusted R-squared:  0.2259 
#F-statistic: 3.918 on 1 and 9 DF,  p-value: 0.07915

test<-lm(pc18u$alpha.impmkt2~pc18u$maternal_proportion,na.action=na.omit)
summary(test)
#Estimate Std. Error t value Pr(>|t|)
#(Intercept)                 -1.495      1.182  -1.265    0.211
#pc18u$maternal_proportion    2.374      1.513   1.570    0.122

#Residual standard error: 1.467 on 53 degrees of freedom
#(1663 observations deleted due to missingness)
#Multiple R-squared:  0.04442,	Adjusted R-squared:  0.02639 
#F-statistic: 2.463 on 1 and 53 DF,  p-value: 0.1225


pc18m<-pc18m[is.finite(pc18m$dN.dS),]
test<-lm(pc18m$dN.dS~pc18m$maternal_proportion,na.action=na.omit)
summary(test)

#                           Estimate Std. Error t value Pr(>|t|)
#(Intercept)               0.0005415  0.3019627   0.002    0.999
#pc18m$maternal_proportion 0.1488614  0.3517114   0.423    0.682

#Residual standard error: 0.129 on 9 degrees of freedom
#(329 observations deleted due to missingness)
#Multiple R-squared:  0.01952,	Adjusted R-squared:  -0.08943 
#F-statistic: 0.1791 on 1 and 9 DF,  p-value: 0.682


pc18u<-pc18u[is.finite(pc18u$dN.dS),]
test<-lm(pc18u$dN.dS~pc18u$maternal_proportion,na.action=na.omit)
summary(test)
#                           Estimate Std. Error t value Pr(>|t|)
#(Intercept)                0.141779   0.146340   0.969    0.337
#pc18u$maternal_proportion -0.004177   0.187060  -0.022    0.982

#Residual standard error: 0.1812 on 52 degrees of freedom
#(1643 observations deleted due to missingness)
#Multiple R-squared:  9.589e-06,	Adjusted R-squared:  -0.01922 
#F-statistic: 0.0004987 on 1 and 52 DF,  p-value: 0.9823

#Figure S3
hist(pc18$maternal_proportion,main="Distribution of maternal expression bias", xlab="Matern expression proportion",ylab="Number of genes")








