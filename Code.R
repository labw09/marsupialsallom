library(geomorph)
library(ggplot2)
library(ggfortify)
library(convevol)
library(phytools)
library(phangorn)
library(RRphylo)
library(ape)
library(RRPP)
library(mvMORPH)
library(smatr)
library(RColorBrewer)

#### Load data
data<-read.csv("Dataset.csv")
tree<-read.tree("Giannini_marsupials_tree.txt")
tree$tip.label[3]="rynrap"

#### Extract morpho data
morph_data<-data[,c(8:22)]

#### Extract GeoMean
GeoMean<-nrow(morph_data)
for (i in 1:GeoMean) {
      GeoMean[i]<-exp(mean(as.numeric(log(morph_data[i,]))))
}
morph_data$GM<-GeoMean

### Extract Mosimman's ratios
GMdata<-list()
for (i in 1:nrow(morph_data)) {
      GMdata[[i]]<-morph_data[i,1:15]/morph_data[i,16]
}
GMdata<-do.call(rbind, GMdata)
GMdata$GM<-GeoMean

## log-transform shape-ratio data
GMdata<-log(GMdata)
GMdata$Genus<-data$Genera
GMdata$Species<-paste(data$Genera,data$Species)
GMdata$Family<-data$Family
GMdata$Suborder<-data$Suborder
GMdata$Order<-data$Order
GMdata$Diet<-data$Diet

# Assign Superorder 
Superorder<-ifelse(GMdata$Order=="Dasyuromorphia"|GMdata$Order=="Diprotodontia"|GMdata$Order=="Microbiotheria"|GMdata$Order=="Notonycteromorphia"|GMdata$Order=="Peramelemorphia","Australiadelphia","Ameridelphia")
GMdata$Superorder<-Superorder

# Assign superorder-diet interaction
GMdata$SDi<-as.character(interaction(GMdata$Superorder,GMdata$Diet))

### Upload phylogeny and prune to match complete dataset
tree2<-keep.tip(tree,gsub(" ", "",paste(tolower(substr(sub(" .*", "", sort(levels(as.factor(GMdata$Species)))),1,3)),tolower(substr(sub(".* ", "", sort(levels(as.factor(GMdata$Species)))),1,3)))))
tree3<-keep.tip(tree,levels(as.factor(rownames(PCAdf))))

### Perform individual shape~size regressions
lmdata<-list()
for (i in 1:length(unique(GMdata$Species))) {
      lmdata[[i]]<-lm(GM ~ CIL + ZB + BB + OH + ORB + LN + LPAL + BPAL + Upos + LD + HD + Lpos + HM + HC + LC,data = subset(GMdata,Species==unique(GMdata$Species)[i]))
}

### Extract n values of individual regressions
n<-vector()
for (i in 1:length(unique(GMdata$Species))) {
   n[i]<-lmdata[[i]]$qr$rank
}

### Extract coefficients from individual regressions
coefficients<-list()
for (i in 1:length(unique(GMdata$Species))) {
      coefficients[[i]]<-lmdata[[i]]$coefficients
}
names(coefficients)<-unique(GMdata$Species)
### Extract number of coefficients per individual regression
ncoef<-vector()
for (i in 1:length(unique(GMdata$Species))) {
      ncoef[i]<-length(na.omit(coefficients[[i]]))
}
names(ncoef)<-unique(GMdata$Species)

## Create dataframe of coefficients
coef<-do.call(cbind, coefficients[c(names(ncoef[ncoef==15]))])
coef<-as.data.frame(coef[-c(1,16),])
coef<-as.data.frame(t(coef))

# Add Order variable
coeforder<-list()
for (i in 1:nrow(coef)) {
   coeforder[i]<-unique(GMdata$Order[GMdata$Species==rownames(coef)[i]])
}
coef$Order<-coeforder
coef$Order<-as.character(coef$Order)

# Add Superorder variable
coefsuperorder<-list()
for (i in 1:nrow(coef)) {
   coefsuperorder[i]<-unique(GMdata$Superorder[GMdata$Species==rownames(coef)[i]])
}
coef$Superorder<-coefsuperorder
coef$Superorder<-as.character(coef$Superorder)

# Add species variable
coef$Species<-rownames(coef)

# Add diet variable
coefdiet<-list()
for (i in 1:nrow(coef)) {
   coefdiet[i]<-unique(GMdata$Diet[GMdata$Species==rownames(coef)[i]])
}
coef$Diet<-coefdiet
coef$Diet<-as.character(coef$Diet)

# Add superorder*diet interaction term
coefSDi<-list()
for (i in 1:nrow(coef)) {
   coefSDi[i]<-unique(GMdata$SDi[GMdata$Species==rownames(coef)[i]])
}
coef$SDi<-coefSDi
coef$SDi<-as.character(coef$SDi)

# Reorder and rename rows to match phylo tips
coef<-coef[sort(rownames(coef)),]
rownames(coef)<-gsub(" ", "",paste(tolower(substr(sub(" .*", "", sort(levels(as.factor(rownames(coef))))),1,3)),tolower(substr(sub(".* ", "", sort(levels(as.factor(rownames(coef))))),1,3))))

# Reorder rows to match tip order
coef<-coef[match(tree3$tip.label,rownames(coef)),]

## Create allometric morphospace (allometric space)
PCA<-prcomp(coef[,1:14],scale. = TRUE)
PCAdf<-PCA$x
PCAdf<-as.data.frame(PCAdf)
PCAdf$Order<-coef$Order
PCAdf$Superorder<-coef$Superorder
PCAdf$Diet<-coef$Diet
PCAdf$SDi<-coef$SDi

### Prune phylogeny to match PCA dataset
tree3<-keep.tip(tree,levels(as.factor(rownames(PCAdf))))

# Reorder row to match tip order in tree
PCAdf<-PCAdf[match(tree3$tip.label,rownames(PCAdf)),]

# Plot allometric morphospace 
autoplot(PCA,data=coef,colour="Superorder",frame=TRUE)

# Test of morphological disparity
morphol.disparity(PCAdf[,1:4]~1,groups = PCAdf$Superorder,data = PCAdf,iter = 9999,partial = TRUE)

### SMA regressions
smaD1<-sma(CIL + ZB + BB + OH + ORB + LN + LPAL + BPAL + Upos + LD + HD + Lpos + HM + HC + LC ~GM + Diet,log = "xy",data = GMdata,multcomp = TRUE,multcompmethod = "adjusted")
smaD2<-sma(CIL + ZB + BB + OH + ORB + LN + LPAL + BPAL + Upos + LD + HD + Lpos + HM + HC + LC ~GM * Diet,log = "xy",data = GMdata,multcomp = TRUE,multcompmethod = "adjusted")
smaO1<-sma(CIL + ZB + BB + OH + ORB + LN + LPAL + BPAL + Upos + LD + HD + Lpos + HM + HC + LC ~GM * Order,log = "xy",data = GMdata,multcomp = TRUE,multcompmethod = "adjusted")
smaO2<-sma(CIL + ZB + BB + OH + ORB + LN + LPAL + BPAL + Upos + LD + HD + Lpos + HM + HC + LC ~GM + Order,log = "xy",data = GMdata,multcomp = TRUE,multcompmethod = "adjusted")
smaS1<-sma(CIL + ZB + BB + OH + ORB + LN + LPAL + BPAL + Upos + LD + HD + Lpos + HM + HC + LC ~GM + Superorder,log = "xy",data = GMdata,multcomp = TRUE,multcompmethod = "adjusted")
smaS2<-sma(CIL + ZB + BB + OH + ORB + LN + LPAL + BPAL + Upos + LD + HD + Lpos + HM + HC + LC ~GM * Superorder,log = "xy",data = GMdata,multcomp = TRUE,multcompmethod = "adjusted")
smaSD1<-sma(CIL + ZB + BB + OH + ORB + LN + LPAL + BPAL + Upos + LD + HD + Lpos + HM + HC + LC ~GM * SDi,log = "xy",data = GMdata,multcomp = TRUE,multcompmethod = "adjusted")
smaSD2<-sma(CIL + ZB + BB + OH + ORB + LN + LPAL + BPAL + Upos + LD + HD + Lpos + HM + HC + LC ~GM + SDi,log = "xy",data = GMdata,multcomp = TRUE,multcompmethod = "adjusted")

### plot allometric trajectories
lmunique<-lm.rrpp(CIL + ZB + BB + OH + ORB + LN + LPAL + BPAL + Upos + LD + HD + Lpos + HM + HC + LC ~GM * SDi,data = GMdata,iter = 9999)
plotAllometry(lmunique,size = GMdata$GM,logsz = TRUE,pch=as.numeric(as.factor(GMdata$Superorder)),col=as.numeric(as.factor(GMdata$Diet)),method = "PredLine")

lmuniqueS_Animal<-lm.rrpp(CIL + ZB + BB + OH + ORB + LN + LPAL + BPAL + Upos + LD + HD + Lpos + HM + HC + LC ~GM * Species, data = Animal,iter = 9999)
lmuniqueS_Herb<-lm.rrpp(CIL + ZB + BB + OH + ORB + LN + LPAL + BPAL + Upos + LD + HD + Lpos + HM + HC + LC ~GM * Species, data = Herb,iter = 9999)
lmuniqueS_Myco<-lm.rrpp(CIL + ZB + BB + OH + ORB + LN + LPAL + BPAL + Upos + LD + HD + Lpos + HM + HC + LC ~GM * Species, data = Myco,iter = 9999)
lmuniqueS_Omni<-lm.rrpp(CIL + ZB + BB + OH + ORB + LN + LPAL + BPAL + Upos + LD + HD + Lpos + HM + HC + LC ~GM * Species, data = Omni,iter = 9999)

saveAnimal<-plotAllometry(lmuniqueS_Animal,size = Animal$GM,logsz = TRUE,pch=19,
                          col=as.numeric(as.factor(Animal$Superorder)),method = "PredLine")
saveHerb<-plotAllometry(lmuniqueS_Herb,size = Herb$GM,logsz = TRUE,pch=19,
                        col=as.numeric(as.factor(Herb$Superorder)),method = "PredLine")
saveMyco<-plotAllometry(lmuniqueS_Myco,size = Myco$GM,logsz = TRUE,pch=19,
                        col=as.numeric(as.factor(Myco$Superorder)),method = "PredLine")
saveOmni<-plotAllometry(lmuniqueS_Omni,size = Omni$GM,logsz = TRUE,pch=19,
                        col=as.numeric(as.factor(Omni$Superorder)),method = "PredLine")
## Plot raw disparity values
disparity_diet<-c(Animalivory=15.364393, Herbivory=7.363011,Mycophagy=45.904377, Omnivory=4.729101)
barplot(disparity_diet, ylab="Procrustes Variance (PV)")

disparity_super<-c("Ameridelphia"=10.23439, Australidelphia=12.98418)
barplot(disparity_super, ylab="Procrustes Variance (PV)")

disparity_order<-c(Dasyuromorphia=20.753145, Diprotodontia=10.365109, Microbiotheria=3.699027, Peramelemorphia=6.799506, Didelphimorphia= 11.549129, Paucituberculata=3.003332)
barplot(disparity_order, ylab="Procrustes Variance (PV)",width=10)
text(srt=45)
legend.text=c("Dasyuromorphia", "Didelphimorphia", "Diprotodontia",
              "Microbiotheria","Paucituberculata", "Peramelemorphia")


### Plot phylomorphospace
col<-gg_color_hue(4)
names(col)<-levels(as.factor(PCAdf$Diet))
col<-col[match(PCAdf$Diet, names(col))]
names(col)<-rownames(PCAdf)
col<-c(col[tree3$tip.label],rep("black",tree3$Nnode))
names(col)<-1:(length(tree3$tip.label)+tree3$Nnode)
phylomorphospace(tree = tree3,PCAdf[,1:2],control = list(col.node=col),node.size=c(0.3,1.5),xlab="PC1 (44.67%)",ylab="PC2 (14.67%)")

## Pairwise comparison in residual variance - similar allometry-corrected disparity test
summary(pairwise(lmunique,groups = as.factor(GMdata$Superorder)),test.type="var")

## Create Tanglegram with allometric space
col.gp<-brewer.pal(length(levels(as.factor(PCAdf$SDi))),"Set1")
names(col.gp) <- levels(as.factor(PCAdf$SDi))
col.gp <- col.gp[match(PCAdf$SDi, names(col.gp))]
plot(cophylo(tree3,upgma(dist(PCAdf[,1:14]),method = "average")),link.lty="solid",link.type="curved",link.col=col.gp)

## Convergence C1-C4
convsig<-convratsig(tree3,as.matrix(PCAdf[,2:5]),rownames(PCAdf[PCAdf$Diet=="Animalivory",]),nsim=300)
convsigHerb<-convratsig(tree3,as.matrix(PCAdf[,2:5]),rownames(PCAdf[PCAdf$Diet=="Herbivory",]),nsim=300)
convsigOmni<-convratsig(tree3,as.matrix(PCAdf[,2:5]),rownames(PCAdf[PCAdf$Diet=="Omnivory",]),nsim=300)
convsigOmni<-convratsig(tree3,as.matrix(PCAdf[,2:5]),rownames(PCAdf[PCAdf$Diet=="Omnivory",]),nsim=1000)

## Create simmaps for competing evol model testing
Superorder<-as.factor(PCAdf$Superorder)
names(Superorder)<-rownames(PCAdf)
simmapS<-make.simmap(tree3,Superorder)

Diet<-as.factor(PCAdf$Diet)
names(Diet)<-rownames(PCAdf)
simmapD<-make.simmap(tree3,Diet)

SDi<-as.factor(PCAdf$SDi)
names(SDi)<-rownames(PCAdf)
simmapSDi<-make.simmap(tree3,SDi)
#BM1
mvBM(tree,PCAdf[,1:4],model = "BM1")
#BMMs
mvBM(simmapS,PCAdf[,1:4],model = "BMM")
#BMMd
mvBM(simmapD,PCAdf[,1:4],model = "BMM")
#BMMsd
mvBM(simmapSDi,PCAdf[,1:4],model = "BMM")
#BMMsm
mvBM(simmapS,PCAdf[,1:4],model = "BMM",param=list(smean=FALSE))
#BMMdm
mvBM(simmapD,PCAdf[,1:4],model = "BMM",param=list(smean=FALSE))
#BMMsdm
mvBM(simmapSDi,PCAdf[,1:4],model = "BMM",param=list(smean=FALSE))
#EB
mvEB(tree3,PCAdf[,1:4])
#OU1
mvOU(tree3,PCAdf[,1:4])
#OUs
mvOU(simmapS,PCAdf[,1:4],model = "OUM")
#OUd
mvOU(simmapD,PCAdf[,1:4],model = "OUM")
#OUsd
mvOU(simmapSDi,PCAdf[,1:4],model = "OUM")

#Broken stick model (per: https://en.proft.me/2016/11/15/principal-component-analysis-pca-r/)
ev<-PCA$sdev^2
function(ev) { 
  # Broken stick model (MacArthur 1957)
  n = length(ev)
  bsm = data.frame(j=seq(1:n), p=0)
  bsm$p[1] = 1/n
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p = 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}

