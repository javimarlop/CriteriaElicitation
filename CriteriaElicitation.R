# Author: Javier Martínez-López (javi.martinez.lopez # gmail.com)

library(googlesheets)
library(vegan)
library(ade4)
library(Prize)
library(ggplot2)
library(scales)

gs_ls()
data<-gs_title('PT AQ stakeholder assessment of ES in CS5 (Responses)')
df<-gs_read(data)

ndim<-11

# add stakeholders groups
g1<-'Política / Governança (ex: Governança Ambiental, da Pesca e Agricultura, Marinha, agências nacionais)'
g2<-'Administração Pública (ex: Administração Regional, Municípios, Freguesias)'
g3<-'Cidadãos (ex: moradores, proprietários de residências, grupos sub-representados e vulneráveis)'
g4<-'Ciência (ex: professores, investigadores ou técnicos em universidades / instituições locais ou independentes)'
g5<-'Grupos de interesse (ex: associações locais, organizações não governamentais (ONGs), organizações profissionais)'
g6<-'Negócios (ex: Indústria, Turismo, Agricultura, Pesca, pequenas empresas, empresas nacionais ou multinacionais com interesses locais)'
g7<-'Serviços (ex: Portos e transporte comercial, transporte público de passageiros)'

dfnum<-as.data.frame(df)

dfnum[df==g1]<-'PolGov'
dfnum[df==g2]<-'AdmPub'
dfnum[df==g3]<-'Cidad'
dfnum[df==g4]<-'Cienc'
dfnum[df==g5]<-'GrpInt'
dfnum[df==g6]<-'Neg'
dfnum[df==g7]<-'Serv'

mli<-'muito menos importante'#'much less important'
li<-'menos importante'#'less important'
ei<-'igualmente importante'#'equally important'
mi<-'mais importante'#'more important'
mmi<-'muito mais importante'#'much more important'

dfnum[df==mli]<-1/4
dfnum[df==li]<-1/2
dfnum[df==ei]<-1
dfnum[df==mi]<-2
dfnum[df==mmi]<-4

dff<-as.data.frame(matrix(NA,1,(ndim+1)))
names(dff)<-c('id','s1','s2','s3','s4','s5','s6','s7','s8','s9','s10','s11')

#CorrectScores<-function(){


# To be done for each individual in a loop
for(n in 1:dim(dfnum)[1]){
#n<-1
print(n)
dfnum2<-dfnum[n,-(1:2)]

a<-matrix(NA,ndim,ndim,byrow=T)
combs<-factorial(ndim)/(factorial(ndim-2)*factorial(2))
ind<-(ndim-1):1
l<-NULL
for(k in 1:length(ind)){l[k]<-sum(ind[1:k])}
init0<-l[1:(length(l)-1)]+1
init<-c(1,init0)
init[length(init)]<-combs
for(j in 1:(ndim-1)){
a[j,j]<-1
a[j,-(1:j)]<-as.numeric(dfnum2[1,init[j]:l[j]])
}
a[ndim,ndim]<-1

dfa<-as.data.frame(a)
names(dfa)<-names(dff)[2:12]
rownames(dfa)<-names(dff)[2:12]
write.table(dfa,paste('ind_',n,'.tsv',sep=''),sep='\t',quote=F)

A<-ahmatrix(a)

ea<-eigen(A@ahp_matrix)
dff[n,1]<-dfnum[n,2]
dff[n,2:(ndim+1)]<-rescale(as.numeric(ea$vectors[,1]),to=c(1,5)) # 2:6
}


     mat <- matrix(nrow = 17, ncol = 1, data = NA)
     mat[,1] <- c('ind_1.tsv',
                 'ind_2.tsv',
                 'ind_3.tsv',
                 'ind_4.tsv',
                 'ind_5.tsv',
                 'ind_6.tsv',
                 'ind_7.tsv',
                 'ind_8.tsv',
                 'ind_9.tsv',
                 'ind_10.tsv',
                 'ind_11.tsv',
                 'ind_12.tsv',
                 'ind_13.tsv',
                 'ind_14.tsv',
                 'ind_15.tsv',
                 'ind_16.tsv',
                 'ind_17.tsv'
			)
     rownames(mat) <- c('ind1','ind2','ind3', 'ind4', 'ind5', 'ind6', 'ind7', 'ind8', 'ind9', 'ind10', 'ind11', 'ind12', 'ind13', 'ind14', 'ind15', 'ind16', 'ind17')
     colnames(mat) <- c('individual_judgement')
     
     # non-weighted aggregation
     res <- gaggregate(srcfile = mat, method = 'geometric', simulation = 500)

# consistency ratio of the aggregated group judgement
gcr<-GCR(res)
incons<-ICR(res)
okind<-incons<=0.15
okind[2]<-FALSE

crplot(ICR(res), angle = 45)
ggsave('individual_consistency_ratio.png')

# Distance between individual opinions and the aggregated group judgement
dplot(IP(res))
ggsave('Distance2Group.png')

#dmh<-dist(dff[,-1])
dff2<-dff[okind,]
dmh<-dist(dff2[,-1])

mds<-metaMDS(dmh)

hclust(dmh,"ward.D2")->mds_hclust

coph<-cophenetic(mds_hclust)
cophval<-cor(coph,dmh)

q25<-quantile(mds_hclust$height)[2]
q50<-quantile(mds_hclust$height)[3]
q75<-quantile(mds_hclust$height)[4]

#cutree(mds_hclust,h=q75)->mds_hclust_mean # Dissimilarity threshold
cutree(mds_hclust,k=2)->mds_hclust_mean # Number of classes

#km<-kmeans(dmh,3)
#mds_hclust_mean<-km$cluster

ncl<-length(unique(mds_hclust_mean))

png('hclust_classes.png')
plot(mds_hclust,hang=-1)
rect.hclust(mds_hclust,k=ncl)
dev.off()

#mds_km<-kmeans(dmh,3)

mds_xy <- data.frame(mds$points)
#mds_xy$cluster<-as.vector(mds_km$cluster)
mds_xy$cluster<-as.vector(mds_hclust_mean)
mds_xy$groups<-dff2$id

png('nmds_classes.png')
ggplot(mds_xy, aes(MDS1, MDS2, color = as.factor(cluster),label = rownames(mds_xy))) + geom_point() + theme_bw() + geom_label()
dev.off()

png('nmds_initial_groups.png')
ggplot(mds_xy, aes(MDS1, MDS2, color = as.factor(groups),label = rownames(mds_xy))) + geom_point() + theme_bw() + geom_label()
dev.off()

dff2$classes<-mds_hclust_mean
lig<-length(unique(dff2$id))

h<-0
z<-matrix(NA,ndim,ncl,byrow=T)
z2<-matrix(NA,ndim,lig,byrow=T)
z4<-matrix(NA,ndim,1,byrow=T)
for(n in 2:(dim(dff2)[2]-1)){
h<-h+1
x<-aggregate(dff2[,n]~classes,data=dff2, FUN=function(x) mean(x))
x2<-aggregate(dff2[,n]~id,data=dff2, FUN=function(x) mean(x))
z[h,]<-x[,2]
z2[h,]<-x2[,2]
z4[h,]<-mean(dff2[,n])
}

z3<-as.data.frame(z2)
names(z3)<-x2[,1]

z5<-as.data.frame(z)
z6<-as.data.frame(z4)

write.table(z5[-2,],'weights_new_groups.csv',sep=',',row.names=F)
write.table(z3[-2,],'weights_initial_groups.csv',sep=',',row.names=F)
write.table(z6[-2,],'weights_nogroups.csv',sep=',',row.names=F)

#png('nmds.png')
#plot(mds)
#s.class(mds$points,fac=as.factor(1:dim(dff)[1]))#,col=1:length(unique(dff$id)))
#dev.off()


