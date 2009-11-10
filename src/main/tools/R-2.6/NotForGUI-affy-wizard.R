# R script for Affymetrix wizard served as a Web Service
# Created by JTT on 6.6.2008
# Modified by JTT on 18.6.2008
#
# Analysis steps:
# 1. Read in CEL-files
# 2. Normalize the data
# 3. Write out normalized data
# 4. Run a statistical test on the normalized data
# 5. Write out at least a hundred most significant genes
# 6. Perform hierarchical clustering
# 7. Write out two Newick-formatted trees
# 8. Write out a heatmap image
# 9. Annotate results
# 10. Write out a text file containing the annotations
# 11. Run the term enrichment analysis for KEGG pathways
# 12. Write the result to an HTML file
# 13. Run the term enrichment analysis for GO pathways
# 14. Write the result to an HTML file
#
# Output files are:
# normalized.tsv     normalized data
# two-sample.tsv     DEG
# hc-genes.txt       gene clustering
# hc-chips.txt       chip clustering
# heatmap.png        heatmap image
# annotations.txt    annotations
# hypergeo.html      term enrichment analysis


# Normalization
library(affy)
dat0<-ReadAffy()
dat1<-exprs(rma(dat0))
dat2<-as.data.frame(round(dat1, digits=2))
write.table(dat2, file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)

# Statistical testing
# Phenodata (groups) is needed here!
library(limma)
groups<-c(1,1,2,2)
design<-model.matrix(~groups)
fit<-lmFit(dat2, design)
fit<-eBayes(fit)
tab<-toptable(fit, coef=2, number=nrow(fit), adjust.method="BH")
rows<-as.numeric(row.names(tab))
rows<-rows[tab$adj.P.Val<=0.05]
p<-tab$adj.P.Val[tab$adj.P.Val<=0.05]
M<-tab$logFC[tab$adj.P.Val<=0.05]
dat3<-dat2[rows,]
if(nrow(dat3)==0) {
   dat3<-data.frame(dat2[order(p),], p.raw=round(p[order(p)], digits=4))
}
write.table(na.omit(data.frame(dat3, p.adjusted=round(p, digits=4), FoldChange=M)), file="two-sample.tsv", sep="\t", row.names=T, col.names=T, quote=F)

# Hierarchical clustering
library(amap)
library(ape)
clustg<-hcluster(x=dat3, method="pearson", link="average")
clustc<-hcluster(x=t(dat3), method="pearson", link="average")
heatcol<-colorRampPalette(c("Blue", "Yellow"))(32)
# Writing the Newick files is slow (6000 genes takes almost three minutes)
write.tree(as.phylo(clustg), "hc-genes.txt")
write.tree(as.phylo(clustc), "hc-chips.txt")
bitmap(file="heatmap.png", width=1500/72, height=1500/72)
heatmap(x=as.matrix(dat3), Rowv=as.dendrogram(clustg), Colv=as.dendrogram(clustc), col=heatcol)
dev.off()

# Annotation
# Annotation is very slow (6000 genes takes almost 4 minutes)
library(package=dat0@annotation, character.only=T)
library(annaffy)
annot.cols<-aaf.handler()
annot.table<-aafTableAnn(row.names(dat3), dat0@annotation, annot.cols)
expression<-data.frame(FoldChange=M)
rownames(expression)<-rownames(dat3)
expression<-aafTable(items=expression, signed=T)
expression@probeids<-rownames(dat3)
annot.table2<-merge(annot.table, expression)
pvalues<-data.frame(P.Adjusted=p)
rownames(pvalues)<-rownames(dat3)
pvalues<-aafTable(items=pvalues)
pvalues@probeids<-rownames(dat3)
annot.table3<-merge(annot.table2, pvalues)
# saveHTML(annot.table3, "annotations.html", title="Annotations for diffentially expressed genes")
saveText(annot.table3, "annotations.tsv")

# KEGG and GO enrichment analysis
library(GOstats)
library(KEGG)
library(GO)
pcut<-0.05
choisedirec<-"over"
ontology<-"KEGG"
env<-paste(dat0@annotation, "ENTREZID", sep="")
alleg<-get(env)
alleg<-as.list(alleg)
alleg<-unlist(alleg)
alleg<-as.data.frame(alleg)
myids<-unique(alleg[rownames(dat3),])
params<-new("KEGGHyperGParams", geneIds=myids, annotation=dat0@annotation, pvalueCutoff=pcut, testDirection=choisedirec)
result<-hyperGTest(params)
if(ontology=="KEGG" & sum(pvalues(result)<pcut)>=1) {
   ID<-names(pvalues(result)[pvalues(result)<=pcut])
   Pvalue<-pvalues(result)[pvalues(result)<=pcut]
   OddsRatio<-oddsRatios(result)[pvalues(result)<=pcut]
   ExpCount<-expectedCounts(result)[pvalues(result)<=pcut]
   Count<-geneCounts(result)[pvalues(result)<=pcut]
   Size<-universeCounts(result)[pvalues(result)<=pcut]
   Term<-as.vector(t(as.data.frame(mget(names(result@oddsRatios[result@pvalues<=pcut]), KEGGPATHID2NAME)))[,1])   
   write(x="<HTML>", file="hypergeo.html", append=T) 
   write(x="<BODY>", file="hypergeo.html", append=T)
   write(x="<TABLE border=1>", file="hypergeo.html", append=T)
   write(x="<CAPTION> KEGG - test for over-representation </CAPTION>", file="hypergeo.html", append=T)
   write(x="<TR> <TH>  </TH> <TH> ID </TH> <TH> Pvalue </TH> <TH> OddsRatio </TH> <TH> ExpCount </TH> <TH> Count </TH> <TH> Size </TH> <TH> Term </TH>  </TR>", file="hypergeo.html", append=T)
   for(i in 1:length(ID)) {
      write(x=paste("<TR> <TD> ", i, " </TD> <TD> ", ID[i], " </TD> <TD> ", round(Pvalue[i], digits=2), " </TD> <TD> ", round(OddsRatio[i], digits=2), " </TD> <TD> ", round(ExpCount[i], digits=2), " </TD> <TD> ", Count[i], " </TD> <TD> ", Size[i], " </TD> <TD> ", "<a href=", cat("\""), "http://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=gn&dbkey=pathway&keywords=", ID[i], "&page1", cat("\""), ">", Term[i], "</a>", "</TD> </TR>", sep=""), file="hypergeo.html", append=T)
   }
   write(x="</BODY>", file="hypergeo.html", append=T)
   write(x="</HTML>", file="hypergeo.html", append=T)      
} 
if(ontology=="KEGG" & sum(pvalues(result)<pcut)<1) {
   write(x="<HTML>", file="hypergeo.html", append=T) 
   write(x="<BODY>", file="hypergeo.html", append=T)
   write(x="No significant results found! <br>", file="hypergeo.html", append=T)
   write(x="</BODY>", file="hypergeo.html", append=T)
   write(x="</HTML>", file="hypergeo.html", append=T)      
}
ontology<-"GO"
params<-new("GOHyperGParams", geneIds=myids, annotation=dat0@annotation, ontology="BP", pvalueCutoff=pcut, conditional=FALSE,testDirection=choisedirec)
resultBP<-hyperGTest(params)
params<-new("GOHyperGParams", geneIds=myids, annotation=dat0@annotation, ontology="MF", pvalueCutoff=pcut, conditional=FALSE,testDirection=choisedirec)
resultMF<-hyperGTest(params)
params<-new("GOHyperGParams", geneIds=myids, annotation=dat0@annotation, ontology="CC", pvalueCutoff=pcut, conditional=FALSE,testDirection=choisedirec)
resultCC<-hyperGTest(params)
if( sum(pvalues(resultCC)<pcut) + sum(pvalues(resultBP)<pcut) + sum(pvalues(resultMF)<pcut) >=1) {
   htmlReport(resultBP, "hypergeo.html", append=T)
   htmlReport(resultMF, "hypergeo.html", append=T)
   htmlReport(resultCC, "hypergeo.html", append=T)
}
if( sum(pvalues(resultCC)<pcut) + sum(pvalues(resultBP)<pcut) + sum(pvalues(resultMF)<pcut) <1) {
   write(x="<HTML>", file="hypergeo.html", append=T) 
   write(x="<BODY>", file="hypergeo.html", append=T)
   write(x="No significant results found! <br>", file="hypergeo.html", append=T)
   write(x="</BODY>", file="hypergeo.html", append=T)
   write(x="</HTML>", file="hypergeo.html", append=T)   
}