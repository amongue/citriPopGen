#citri popgen analyses for publication. This is an abridged file with the most relevant analyses and results only.
library(data.table)
library(FSA)
setwd("C:/Users/andre/Documents/Documents/Documents/CitriPopgen/mac ported data and code/")
pcit6<-fread("Pcit.DTOL.matrix.v2.1.csv",stringsAsFactors = F)

#adding Hollie's ASE data
ase<-fread("hollie_ase_base.txt",header=T,stringsAsFactors = F)
table(ase$lineage,ase$actually_sig)
#            no yes
#biparental   9   0
#cp1          0  36
#maternal     0 483
#unbiased   187   0
#wye          0  14

#for now I think we want to subset out the biparentals and the lineage biased genes?
#biparental are the obvious case of diploid expression in a male
#lineage specific are weirder..based on Hollie's plot they could be maternally or paternally biased
#it depends on the direction of the cross...so let's remove those 50 genes?
abr<-as.data.table(cbind(ase$Gene,ase$lineage),stringsAsFactors=F)
colnames(abr)<-c("Gene","Lineage")
pcit7<-merge(pcit6,abr,by="Gene",all.x=T)
#pcit7$Lineage[is.na(pcit7$Lineage)] <- "pass"
table(pcit7$Lineage,pcit7$DE_bias)
table(pcit7$Lineage,pcit7$SPM.Sex.Class)


kruskal.test(pN.pS.P ~ SPM.Sex.Class, data = pcit7)
#Kruskal-Wallis chi-squared = 72.787, df = 2, p-value < 2.2e-16
dunnTest(pN.pS.P ~ SPM.Sex.Class, data = pcit7)

#Comparison        Z      P.unadj        P.adj
#     Female - Male 5.703002 1.177156e-08 2.354311e-08
# Female - Unbiased 7.993741 1.309050e-15 3.927150e-15
#   Male - Unbiased 2.167733 3.017899e-02 0.03017899
boxplot(pN.pS.P ~ SPM.Sex.Class, data = pcit7,notch=T,outline=F)
#within pop variation follows UB < M < F
#would be consistent with 
#ub = constitutatively expressed + some haploid selection
#male haploid selected but only in half the population
#female = diploid selection and only in half the pop

#this pattern should be most apparent with pN and not affect pS?
kruskal.test((pN.P/Non.sites) ~ SPM.Sex.Class, data = pcit7)
#Kruskal-Wallis chi-squared = 38.923, df = 2, p-value = 3.531e-09
dunnTest((pN.P/Non.sites) ~ SPM.Sex.Class, data = pcit7)
#Comparison        Z      P.unadj        P.adj
#     Female - Male 4.378181 1.196740e-05 2.393481e-05
# Female - Unbiased 5.735826 9.703814e-09 2.911144e-08
#   Male - Unbiased 1.314064 1.888246e-01 0.1888246
#so female higher than male, unbiased. M = UB

boxplot((pN.P/Non.sites) ~ SPM.Sex.Class, data = pcit7,notch=T,outline=F,ylim=c(0,0.001))


kruskal.test((pS.P/Syn.sites) ~ SPM.Sex.Class, data = pcit7)
#Kruskal-Wallis chi-squared = 0.021544, df = 2, p-value = 0.9893

boxplot((pS.P/Syn.sites) ~ SPM.Sex.Class, data = pcit7,notch=T,outline=F)

#sanity check, is pN.pS.P calculated correctly? Yes
#head((pcit7$pN.S/pcit7$Non.sites)/(pcit7$pS.S/pcit7$Syn.sites))

#now divergence
kruskal.test(dN.dS ~ SPM.Sex.Class, data = pcit7)
#Kruskal-Wallis chi-squared = 253.93, df = 2, p-value < 2.2e-16
dunnTest(dN.dS ~ SPM.Sex.Class, data = pcit7)
#Comparison         Z      P.unadj        P.adj
#1     Female - Male  7.057846 1.691033e-12 1.691033e-12
#2 Female - Unbiased 15.899891 6.347542e-57 1.904262e-56
#3   Male - Unbiased  7.884405 3.160371e-15 6.320742e-15

boxplot(dN.dS ~ SPM.Sex.Class, data = pcit7,notch=T,outline=F)
#U < M < F
#same as above

#this pattern should be most apparent with pN and not affect pS?
kruskal.test((dN/Non.sites) ~ SPM.Sex.Class, data = pcit7)
#Kruskal-Wallis chi-squared = 59.028, df = 2, p-value = 1.522e-13
dunnTest((dN/Non.sites) ~ SPM.Sex.Class, data = pcit7)
#         Comparison         Z      P.unadj        P.adj
#     Female - Male  6.325048 2.531535e-10 5.063069e-10
# Female - Unbiased -1.756214 7.905194e-02 0.07905194
#   Male - Unbiased -7.119480 1.083348e-12 3.250043e-12
#so M < F = U

boxplot((dN/Non.sites) ~ SPM.Sex.Class, data = pcit7,notch=T,outline=F)


kruskal.test((dS/Syn.sites) ~ SPM.Sex.Class, data = pcit7)
#Kruskal-Wallis chi-squared = 224.34, df = 2, p-value < 2.2e-16
dunnTest((dS/Syn.sites) ~ SPM.Sex.Class, data = pcit7)
#Comparison          Z      P.unadj        P.adj
#     Female - Male   2.977294 2.908046e-03 2.908046e-03
# Female - Unbiased -12.584378 2.573632e-36 5.147264e-36
#   Male - Unbiased -13.900320 6.305897e-44 1.891769e-43
boxplot((dS/Syn.sites) ~ SPM.Sex.Class, data = pcit7,notch=T,outline=F)
#M < F < U

#####now let's try with Tamsin's data
kruskal.test(pN.pS.P ~ DE_bias, data = pcit7)
#Kruskal-Wallis chi-squared = 66.865, df = 2, p-value = 3.023e-15
dunnTest(pN.pS.P ~ DE_bias, data = pcit7)
#Comparison         Z      P.unadj        P.adj
#     Female - Male  7.552523 4.269066e-14 1.280720e-13
# Female - Unbiased  6.239162 4.399218e-10 8.798435e-10
#   Male - Unbiased -1.140709 2.539912e-01 0.02539912
#M < U < F
boxplot(pN.pS.P ~ DE_bias, data = pcit7,notch=T,outline=F)

kruskal.test((pN.P/Non.sites) ~ DE_bias, data = pcit7)
#Kruskal-Wallis chi-squared = 81.989, df = 2, p-value < 2.2e-16
dunnTest((pN.P/Non.sites) ~ DE_bias, data = pcit7)
#Comparison          Z      P.unadj        P.adj
#     Female - Male  7.8688328 3.579652e-15 1.073896e-14
# Female - Unbiased  7.6589677 1.874336e-14 3.748672e-14
#   Male - Unbiased -0.1151231 9.083475e-01 0.09083475
#so female higher than male, unbiased. M = UB

boxplot((pN.P/Non.sites) ~ DE_bias, data = pcit7,notch=T,outline=F,ylim=c(0,0.01))


kruskal.test((pS.P/Syn.sites) ~ DE_bias, data = pcit7)
#Kruskal-Wallis chi-squared = 12.784, df = 2, p-value = 0.001675

dunnTest((pS.P/Syn.sites) ~ DE_bias, data = pcit7)
#Comparison         Z      P.unadj       P.adj
#     Female - Male 2.6177908 0.0088501033 0.017700207
# Female - Unbiased 3.3814650 0.0007210043 0.002163013
#   Male - Unbiased 0.7768625 0.4372399090 0.437239909
#M<F, U<F, M=U
boxplot((pS.P/Syn.sites) ~ DE_bias, data = pcit7,notch=T,outline=F)


kruskal.test(dN.dS ~ DE_bias, data = pcit7)
#Kruskal-Wallis chi-squared = 200.33, df = 2, p-value < 2.2e-16
dunnTest(dN.dS ~ DE_bias, data = pcit7)
#Comparison         Z      P.unadj        P.adj
#     Female - Male 13.464364 2.534915e-41 7.604746e-41
# Female - Unbiased 10.255518 1.117751e-24 2.235502e-24
#   Male - Unbiased -2.994195 2.751701e-03 2.751701e-03
#M < U < F
boxplot(dN.dS ~ DE_bias, data = pcit7,notch=T,outline=F)

kruskal.test((dN/Non.sites) ~ DE_bias, data = pcit7)
#Kruskal-Wallis chi-squared = 161.8, df = 2, p-value < 2.2e-16
dunnTest((dN/Non.sites) ~ DE_bias, data = pcit7)
#Comparison         Z      P.unadj        P.adj
#     Female - Male 12.269720 1.317079e-34 3.951237e-34
# Female - Unbiased  8.746955 2.191872e-18 4.383744e-18
#   Male - Unbiased -3.344057 8.256279e-04 8.256279e-04
boxplot((dN/Non.sites) ~ DE_bias, data = pcit7,notch=T,outline=F)

kruskal.test((dS/Syn.sites) ~ DE_bias, data = pcit7)
#Kruskal-Wallis chi-squared = 26.153, df = 2, p-value = 2.094e-06
dunnTest((dS/Syn.sites) ~ DE_bias, data = pcit7)
#Comparison         Z      P.unadj        P.adj
#     Female - Male  5.012918 5.361062e-07 1.608319e-06
# Female - Unbiased  3.257205 1.125152e-03 2.250304e-03
#   Male - Unbiased -1.675402 9.385546e-02 0.09385546
boxplot((dS/Syn.sites) ~ DE_bias, data = pcit7,notch=T,outline=F)
#M=U<F

###now we try alpha
NItgCalc<-function(dn,ds,pn,ps)
{
  unbiasedNI<-sum((ds*pn)/(ps+ds))/sum((ps*dn)/(ps+ds))
  return(unbiasedNI)
}

mcit<-pcit7[which(pcit7$SPM.F < 0.3),]
ucit<-pcit7[which(pcit7$SPM.F > 0.3 & pcit7$SPM.F < 0.7),]
fcit<-pcit7[which(pcit7$SPM.F > 0.7),]

#or
mcit<-pcit7[which(pcit7$DE_bias=="Male"),]
ucit<-pcit7[which(pcit7$DE_bias=="Unbiased"),]
fcit<-pcit7[which(pcit7$DE_bias=="Female"),]

mcit<-mcit[which((mcit$dS + mcit$pS.P)!=0),]
ucit<-ucit[which((ucit$dS + ucit$pS.P)!=0),]
fcit<-fcit[which((fcit$dS + fcit$pS.P)!=0),]
ma<-1-NItgCalc(mcit$dN,mcit$dS,mcit$pN.P,mcit$pS.P)
ua<-1-NItgCalc(ucit$dN,ucit$dS,ucit$pN.P,ucit$pS.P)
fa<-1-NItgCalc(fcit$dN,fcit$dS,fcit$pN.P,fcit$pS.P)

#let's test permute male vs female:
testor<-abs((ma-fa))
size<-2000
bsstata<-rep(0,size)
mwhole<-rbind(mcit,fcit)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$pN.P
mwhole$pS<-mwhole$pS.P
pa<-nrow(mcit)
pb<-nrow(mwhole)
for(i in 1:size)
{
  #next we'll shuffle our combined geneset
  bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#p = 0.2955, male = female

testor<-abs((ma-ua))
size<-1000
bsstata<-rep(0,size)
mwhole<-rbind(mcit,ucit)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$pN.P
mwhole$pS<-mwhole$pS.P
pa<-nrow(mcit)
pb<-nrow(mwhole)
for(i in 1:size)
{
  #next we'll shuffle our combined geneset
  bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.030 M < U

#Female vs unbiased
testor<-abs((fa-ua))
size<-1000
bsstata<-rep(0,size)
mwhole<-rbind(fcit,ucit)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$pN.P
mwhole$pS<-mwhole$pS.P
pa<-nrow(mcit)
pb<-nrow(mwhole)
for(i in 1:size)
{
  #next we'll shuffle our combined geneset
  bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.155

#let's test permute male vs female:
testor<-abs((ma-fa))
size<-1000
bsstata<-rep(0,size)
mwhole<-rbind(mcit,fcit)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$pN.P
mwhole$pS<-mwhole$pS.P
pa<-nrow(mcit)
pb<-nrow(mwhole)
for(i in 1:size)
{
  #next we'll shuffle our combined geneset
  bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.333



#what about DoS?
kruskal.test(DoS.P ~ SPM.Sex.Class, data = pcit7)

boxplot(DoS.P ~ SPM.Sex.Class, data = pcit7,notch=T,outline=F)
abline(h=0)
crosscheck<-rbind(mcit,ucit,fcit)
kruskal.test(DoS.P ~ DE_bias, data = crosscheck)

###Hollie's updated biparentals
hol<-as.data.frame(fread("andres_imp_genes_output_update_march2024.txt"),stringsAsFactors=F)
pcit8<-merge(pcit7,hol,by="Gene")
pcit8<-pcit8[which(is.finite(pcit8$dN.dS)),]
plot(pcit8$maternal_proportion,pcit8$dN.dS,ylim=c(0,1),pch=16)
check<-lm(pcit8$dN.dS~pcit8$maternal_proportion,na.action=na.omit)
summary(check)


plot(pcit8$maternal_proportion,pcit8$DoS.P,ylim=c(-1,1),pch=16)
check<-lm(pcit8$DoS.P~pcit8$maternal_proportion,na.action=na.omit)
summary(check)
abline(lm(pcit8$DoS.P~pcit8$maternal_proportion,na.action=na.omit))
#Estimate Std. Error t value Pr(>|t|)   
#(Intercept)               -0.21737    0.07663  -2.837  0.00477 **
#pcit8$maternal_proportion  0.24812    0.09532   2.603  0.00956 **


pcit8[which(pcit8$maternal_proportion<0.4&pcit8$SPM.F>0.8),]
#jg564
cor(pcit8$maternal_proportion,pcit8$dN.dS,use='pairwise.complete.obs')



pcit9<-pcit8[which(pcit8$maternal_proportion>0.45),]
plot(pcit9$maternal_proportion,pcit9$DoS.P,ylim=c(-1,1),pch=16)
check<-lm(pcit9$DoS.P~pcit9$maternal_proportion,na.action=na.omit)
summary(check)
#Estimate Std. Error t value Pr(>|t|)   
#(Intercept)               -0.20841    0.07907  -2.636  0.00869 **
#  pcit9$maternal_proportion  0.23731    0.09817   2.417  0.01605 * 
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.2586 on 435 degrees of freedom
#(157 observations deleted due to missingness)
#Multiple R-squared:  0.01325,	Adjusted R-squared:  0.01099 
#F-statistic: 5.843 on 1 and 435 DF,  p-value: 0.01605
#is significant even without the paternally-biased points for leverage