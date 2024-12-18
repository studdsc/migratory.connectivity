
# packages
library(ade4)
library(dendextend)
library(gdata)
library(fossil)

# d2h isotope data
bahamas <- read.csv("~bahamas.csv",header=TRUE,na.strings=NA)
belize <- read.csv("~belize.csv",header=TRUE,na.strings=NA)
cuba <- read.csv("~cuba.csv",header=TRUE,na.strings=NA)
dominica <- read.csv("~dominica.csv",header=TRUE,na.strings=NA)
dr <- read.csv("~dr.csv",header=TRUE,na.strings=NA)
jamaica <- read.csv("~jamaica.csv",header=TRUE,na.strings=NA)
mexico <- read.csv("~mexico.csv",header=TRUE,na.strings=NA)
nicaragua <- read.csv("~nicaragua.csv",header=TRUE,na.strings=NA)
puerto.rico <- read.csv("~puerto rico2.csv",header=TRUE,na.strings=NA)
st.martin <- read.csv("~st martin.csv",header=TRUE,na.strings=NA)
united.states <- read.csv("~united states2.csv",header=TRUE,na.strings=NA)

# subset data by site
BAWW.GB <- subset(bahamas,SITE=="Grand Bahama" & SPECIES=="BAWW",select=d2H)
BAWW.LI <- subset(bahamas,SITE=="Long Island" & SPECIES=="BAWW",select=d2H)
BAWW.BZ <- subset(belize,SPECIES=="BAWW",select=d2H)
BAWW.CC <- subset(cuba,SITE=="Cayo Coco" & SPECIES=="BAWW",select=d2H)
BAWW.UB <- subset(cuba,SITE=="Guantanamo Bay" & SPECIES=="BAWW",select=d2H)
BAWW.DO <- subset(dominica,SPECIES=="BAWW",select=d2H)
BAWW.SB <- subset(dr,SITE=="Sierra de Baharuco" & SPECIES=="BAWW",select=d2H)
BAWW.CA <- subset(dr,SITE=="La Canela" & SPECIES=="BAWW",select=d2H)
BAWW.PE <- subset(dr,SITE=="Parque del Este" & SPECIES=="BAWW",select=d2H)
BAWW.FH <- subset(jamaica,SITE=="Font Hill" & SPECIES=="BAWW",select=d2H)
BAWW.BM <- subset(jamaica,SITE=="Blue Mountains" & SPECIES=="BAWW",select=d2H)
BAWW.JM <- subset(jamaica,SPECIES=="BAWW",select=d2H)
BAWW.MX <- subset(mexico,SPECIES=="BAWW",select=d2H)
BAWW.NI <- subset(nicaragua,SPECIES=="BAWW",select=d2H)
BAWW.CR <- subset(puerto.rico,SITE=="Cabo Rojo NWR" & SPECIES=="BAWW",select=d2H)
BAWW.RR <- subset(puerto.rico,SITE=="Roosevelt Roads" & SPECIES=="BAWW",select=d2H)
BAWW.VI <- subset(puerto.rico,SITE=="Vieques" & SPECIES=="BAWW",select=d2H)
BAWW.SM <- subset(st.martin,SPECIES=="BAWW",select=d2H)
BAWW.CH <- subset(united.states,SITE=="Curry Hammock SP, FL" & SPECIES=="BAWW",select=d2H)
BAWW.EV <- subset(united.states,SITE=="Everglades NP, FL" & SPECIES=="BAWW",select=d2H)

# potential origin data
baww.origin <- read.csv("~baww.origin.csv",header=TRUE)
origin <- with(baww.origin,origin2) # spatial growing season H scaled to baww
length(origin) #3,790 cells of potential origin

# create empty matrix and name rows and columns
feather <- with(BAWW.GB, d2H) # separate assignments for each site
assign <- data.frame(matrix(ncol=68, nrow=3970)) #grand bahama has 68 feathers
colnames(assign) <- feather
tassign <- t(assign)
colnames(tassign) <- origin
assign <- t(tassign)
head(assign)

# fill matrix with likelihoods for d2H assignments
for(k in 1:length(origin)){
  for(i in 1:length(feather)){
    assign[k,i]=(1 / (sqrt(2*3.14*14.4))) * exp(-1*((feather[i]-origin[k])^2) / (2*(14.4^2)))
  } 
}

# divide elements in each column by the column sum to create proper probability for d2H assignments
baww.assign <- as.data.frame(rowSums(t(t(assign)/colSums(assign)))/max(rowSums(t(t(assign)/colSums(assign))))) 
names(baww.assign)[1] <- "probability"

## bayes rule assignments with weights from rushing et al. 2017,ecology & evolution DOI: 10.1002/ece3.2605

# abundance weight
10^-0.9

#isotope weight
unwweighted

head(baww.origin)
baww.assign<-within(baww.assign,abundance <- with(baww.origin, abundance/max(abundance)))
baww.assign<-within(baww.assign,posterior <- with(baww.assign, probability*(abundance^0.1258925)/
                                                    sum(probability*(abundance^0.1258925))))
baww.assign<-within(baww.assign,scaled <- with(baww.assign,posterior/max(posterior)))
baww.assign<-within(baww.assign,origin <- with(baww.origin, origin))
baww.assign<-within(baww.assign,x <- with(baww.origin, X))
baww.assign<-within(baww.assign,y <- with(baww.origin, Y))
head(baww.assign)


# posterior assingnment probabilities based on d2H values and abundance
grand.bahama<-with(baww.assign,posterior) 
long.island<-with(baww.assign,posterior)
belize<-with(baww.assign,posterior)
cayo.coco<-with(baww.assign,posterior)
guantanamo.bay<-with(baww.assign,posterior)
dominica<-with(baww.assign,posterior)
la.canela<-with(baww.assign,posterior)
parque.del.este<-with(baww.assign,posterior)
sierra.de.baharcuo<-with(baww.assign,posterior)
blue.mountains<-with(baww.assign,posterior)
font.hill<-with(baww.assign,posterior)
mexico<-with(baww.assign,posterior)
nicaragua<-with(baww.assign,posterior)
cabo.rojo<-with(baww.assign,posterior)
roosevelt.roads<-with(baww.assign,posterior)
vieques<-with(baww.assign,posterior)
st.martin<-with(baww.assign,posterior)
curry.hammock<-with(baww.assign,posterior)
everglades<-with(baww.assign,posterior)

### baww.post for cluster analysis
baww.post<-(cbind(font.hill,mexico,cayo.coco,guantanamo.bay,everglades,
                               curry.hammock,grand.bahama,long.island,sierra.de.baharcuo,la.canela,parque.del.este,
                               vieques,roosevelt.roads,st.martin,belize,nicaragua))

#baww cluster with 3 groups
d <- dist(as.matrix(t.baww.post))  
hc <- hclust(d,method="ward.D")                
plot(hc,axes=TRUE,xlab="",ylab="",sub="",main="")                      
rect.hclust(hc, k=3, border="black")
mtext(side=2,line=2.5,"Height")

# rotated baww solution with groups 
d <- dist(as.matrix(t.baww.post))  
hc <- hclust(d,method="ward.D")  
dend<-as.dendrogram(hc)
plot(dendextend::rotate(dend, c(1:3,1)))
mtext(side=2,line=2.5,"Height")

# repeat assignments for 3 groups and map (fig 1 in publication)
BAWW.MBN.d2H <- rbind(BAWW.MX,BAWW.BZ,BAWW.NI)
BAWW.UBCJD.d2H <- rbind(BAWW.CH,BAWW.EV,BAWW.GB,BAWW.LI,BAWW.CC,BAWW.UB,BAWW.FH,BAWW.SB,BAWW.CA,BAWW.PE)
BAWW.PS.d2H <- rbind(BAWW.RR,BAWW.VI,BAWW.SM)

# wBoot removed from cran on 2022-04-04

## migratory connectivity analysis (see Estimating-Migratory-Connectivity.Rmd for details)

# matrix of pairwise correlations among assignments
(baww.cor.mat<-round(cor(baww.post),2))

# geographic distance among sites data
baww.sites <- read.csv("~baww.sites.csv",header=TRUE,na.strings=NA)

# matrix of pairwise distances among sites
baww.long.lat<-cbind(with(baww.sites,lon),with(baww.sites,lat))
baww.sites.mat<-round(as.matrix(fossil::earth.dist(baww.long.lat,dist=T)),2)

# lower triangle of matrices
baww.cor.low <- gdata::lowerTriangle(baww.cor.mat)
baww.sites.low <- gdata::lowerTriangle(baww.sites.mat)

#  pearson correlation among matrices (high correlation = stronger connectivity)
wBoot::boot.cor.per(baww.cor.low, baww.sites.low, conf.level = 0.95, R = 1000)

# euclidean distances for original matrices
baww.cor.dist <- round(dist(baww.cor.mat, method = "euclidean"),2)
baww.sites.dist <- round(dist(baww.sites.mat, method = "euclidean"),2)

# p=value for migratory connectivity correlation (see Estimating-Migratory-Connectivity for details)
ade4::mantel.rtest(baww.cor.dist, baww.sites.dist, nrepet=1000)
