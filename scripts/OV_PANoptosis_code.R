
if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
}

library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)

library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(text = element_text(family = 'Times'),panel.grid = element_blank())
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 
  return(geneList)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  # set.seed(2024)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Partial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
my_riskplot=function(cli_dat,cols=c("red","blue"),xlab='Samples',
                     a.ylab="Risk score",b.labs="Survival time(year)",cutoff=0,labs=c('A','B')){
  #cli_dat=tcga.risktype.cli
  cli.dat.order=cli_dat[order(cli_dat$Riskscore),c('OS.time','Status','Riskscore','Risktype')]
  fp_dat=data.frame(Samples=1:nrow(cli_dat),cli.dat.order)
  p1=ggplot(fp_dat,aes(x=Samples,y=Riskscore))+geom_point(aes(color=Risktype))+
    scale_colour_manual(values =cols)+
    theme_bw()+labs(x=xlab,y=a.ylab)+
    geom_hline(yintercept=cutoff,colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p1
  p2=ggplot(fp_dat,aes(x=Samples,y=OS.time))+geom_point(aes(col=Status))+theme_bw()+
    scale_colour_manual(values =cols)+
    labs(x=xlab,y=b.labs)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p2
  p.merge=mg_merge_plot(p1,p2,nrow=2,ncol=1,labels = labs)
  return(p.merge)
}
##########

table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
head(mrna_genecode)

##############
tcga.cli<-read.delim('origin_datas/TCGA/Merge_OV_clinical.txt',sep='\t',header = T)
head(tcga.cli)
tcga.cli=data.frame(Samples=paste0(tcga.cli$A0_Samples,'-01'),
                    Age=tcga.cli$A17_Age, Grade=tcga.cli$A7_Grade,Stage=tcga.cli$A6_Stage,
                    Status=tcga.cli$A2_Event,OS.time=tcga.cli$A1_OS)
rownames(tcga.cli)=tcga.cli$Samples
head(tcga.cli)
tcga.cli$OS.time
tcga.cli=tcga.cli %>% drop_na(OS.time)
tcga.cli=tcga.cli[tcga.cli$OS.time>0,]
tcga.cli$OS=ifelse(tcga.cli$Status=='Alive',0,1)

dim(tcga.cli)
table(tcga.cli$Stage)
tcga.cli$Stage[tcga.cli$Stage=='']=NA
tcga.cli$Stage=gsub('Stage ','',tcga.cli$Stage)
tcga.cli$Stage=gsub('[ABC]','',tcga.cli$Stage)

table(tcga.cli$Grade)
tcga.cli$Grade[tcga.cli$Grade%in%c('GB','GX','Not Available')]=NA
dim(tcga.cli)


tcga.data=read.delim('origin_datas/TCGA/TCGA_OV_TPM.txt',row.names = 1,check.names = F)
tcga.data[1:4,1:4]
table(substr(colnames(tcga.data),14,15))
dim(tcga.data)

sample_T=colnames(tcga.data)[which(as.numeric(substr(colnames(tcga.data),14,15))==1)]#
sample_T=intersect(sample_T,tcga.cli$Samples)

tcga.exp=tcga.data[rownames(tcga.data) %in% mrna_genecode$SYMBOL,sample_T]
range(tcga.exp)
tcga.exp=log2(tcga.exp+1)
dim(tcga.exp)

tcga.cli=tcga.cli[sample_T,]
head(tcga.cli)
fivenum(tcga.cli$Age)
tcga.cli$Age1=ifelse(tcga.cli$Age>59,'>59','<=59')

save(tcga.exp,file = 'results/tcga.exp.RData')



#01.PANoptosis ############
dir.create('results/01.PANoptosis_landscape')
PANoptosis.geneset=read.xlsx('origin_datas/PANoptosis.geneset.xlsx')
PANoptosis.geneset_list <- setNames(split(unlist(PANoptosis.geneset, use.names = FALSE), 
                                          rep(1:ncol(PANoptosis.geneset),
                                              each = nrow(PANoptosis.geneset))), 
                                    names(PANoptosis.geneset))
PANoptosis.genesets=Reduce(union,PANoptosis.geneset_list)
length(PANoptosis.genesets)
#486

####
PANoptosis.cox=cox_batch(dat = tcga.exp[PANoptosis.genesets,tcga.cli$Samples],
                   time = tcga.cli$OS.time,event = tcga.cli$OS)
PANoptosis.cox=na.omit(PANoptosis.cox)
head(PANoptosis.cox)
rownames(PANoptosis.cox)=gsub('-','_',rownames(PANoptosis.cox))
p_cutoff=0.05
table(PANoptosis.cox$p.value<p_cutoff)
PANoptosis.cox.fit=PANoptosis.cox[PANoptosis.cox$p.value<p_cutoff,]
dim(PANoptosis.cox.fit)
write.csv(PANoptosis.cox.fit,'results/01.PANoptosis_landscape/PANoptosis.cox.fit.csv')

bioForest=function(rt=null,col){
  #
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}

pdf('results/01.PANoptosis_landscape/PANoptosis_forest.pdf',height = 10,width = 8,onefile = F)
bioForest(rt = PANoptosis.cox.fit[order(PANoptosis.cox.fit$HR),],col=c('blue','orange'))
dev.off()
#####
tcga.maf=getTCGAMAFByCode('OV')
pdf('results/01.PANoptosis_landscape/PANoptosis_MAF.pdf',height = 10,width = 12,onefile = F)
oncoplot(maf = tcga.maf,genes = rownames(PANoptosis.cox.fit))
dev.off()


##PPI###########
ppi.res=read.csv('results/01.PANoptosis_landscape/ppi.res2.csv')
pre.genes=ppi.res$name
length(pre.genes)
#22

#02.############
dir.create('results/02.model')
tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[pre.genes, tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))

######
library(glmnet)
set.seed(2024)
tcga.lasso=get_riskscore.lasso(dat = tcga_model_data[,pre.genes],os = tcga_model_data$OS,os.time = tcga_model_data$OS.time )
length(tcga.lasso$lasso.gene)
tcga.lasso$plot



# tcga.lasso$lasso.gene

######
fmla <- as.formula(paste0("Surv(OS.time, OS) ~"
                          ,paste0(tcga.lasso$lasso.gene,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
cox=step(cox)

lan <- coef(cox)
lan
length(lan)
paste0(round(lan, 3), '*', names(lan),collapse = '+')


gene.coef=data.frame(gene=names(lan),coef=as.numeric(lan))
gene.coef
gene.coef$Type=ifelse(gene.coef$coef>0,'Risk','Protective')
gene.coef.fig=ggplot(gene.coef, aes(x = coef, y = reorder(gene,coef), fill =Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#C75DAA", "#009B9F")) +
  labs(x = 'coefficient', y = "") +
  geom_text(aes(label = round(coef,2),hjust =2), data = subset(gene.coef, coef > 0))+ 
  geom_text(aes(label = round(coef,2), hjust = -1), data = subset(gene.coef, coef < 0))+  
  theme_bw()+ theme(text = element_text(family = 'Times',size = 14,face = 'bold'),legend.position = 'top')
gene.coef.fig


gene.forest=ggforest(cox, data = tcga_model_data, 
                     main = "Hazardratio", fontsize =1.0, 
                     noDigits = 2)
gene.forest




risktype.col=c("orange", "#009B9F")
###########
risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.cli,Riskscore=risk.tcga)
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>median(risk.tcga),'High','Low')


tcga.roc=ggplotTimeROC(time = tcga.risktype.cli$OS.time,
                        status = tcga.risktype.cli$OS,
                      score = tcga.risktype.cli$Riskscore,mks = c(1:5))
tcga.roc

tcga.km.OS=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                   data = tcga.risktype.cli),
                      data=tcga.risktype.cli,
                      conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                      surv.median.line = 'hv',title='TCGA',
                      linetype = c("solid", "dashed","strata")[1],
                      palette = risktype.col,ggtheme = custom_theme(),
                      legend = c(0.8,0.85), # 
                      legend.title = "Risktype",legend.labs=c('High','Low'))
tcga.km.OS=mg_merge_plot(tcga.km.OS$plot,tcga.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
tcga.km.OS

tcga.model.fig=my_riskplot(cli_dat = tcga.risktype.cli,cols = risktype.col,xlab = 'sample',
            b.labs = 'Survival time(days)',cutoff = median(tcga.risktype.cli$Riskscore),labs = c('E',''))
tcga.model.fig


save.image(file = 'project.RData')

##GSE32062 【OK】#######
load('origin_datas/GEO/GSE32062.RData')
GSE32062.cli=GSE32062$Sample
head(GSE32062.cli)
table(GSE32062.cli$tissue)
GSE32062.cli=data.frame(Samples=GSE32062.cli$Acc,Grade=paste0('G',GSE32062.cli$grading),Stage=GSE32062.cli$Stage,
                        OS=GSE32062.cli$`death (1)`,OS.time=GSE32062.cli$`os (m)`)
GSE32062.cli=GSE32062.cli%>%drop_na(OS)
rownames(GSE32062.cli)=GSE32062.cli$Samples
head(GSE32062.cli)
dim(GSE32062.cli)
# 
GSE32062.exp=GSE32062$Exp$GPL6480_41093_Data_col1
GSE32062.exp[1:5,1:5]
dim(GSE32062.exp)
GSE32062.exp=exp_probe2symbol_v2(datExpr = GSE32062.exp,GPL = 'GPL6480')
range(GSE32062.exp)
GSE32062.exp.log=log2(2^GSE32062.exp+1)

GSE32062_model_data <- cbind(GSE32062.cli[, c("OS.time", "OS")],
                             t(GSE32062.exp.log[pre.genes, GSE32062.cli$Samples]))
colnames(GSE32062_model_data) <- gsub('-', '_', colnames(GSE32062_model_data))

risk.GSE32062=as.numeric(lan%*%as.matrix(t(GSE32062_model_data[GSE32062.cli$Samples,names(lan)])))
GSE32062.risktype.cli=data.frame(GSE32062.cli,Riskscore=risk.GSE32062)
GSE32062.risktype.cli$Risktype=ifelse(GSE32062.risktype.cli$Riskscore>median(risk.GSE32062),'High','Low')

GSE32062.roc=ggplotTimeROC(GSE32062.risktype.cli$OS.time,
                           GSE32062.risktype.cli$OS,
                           GSE32062.risktype.cli$Riskscore,mks = c(1:5))
GSE32062.roc



GSE32062.km.OS=ggsurvplot(fit=survfit( Surv(OS.time, OS) ~ Risktype,
                                       data = GSE32062.risktype.cli),
                          data=GSE32062.risktype.cli,
                          conf.int = F,pval = T,fun = "pct",risk.table =T, size = 0.7,
                          surv.median.line = 'hv',title='GSE32062',
                          linetype = c("solid", "dashed","strata")[1],
                          palette = risktype.col,ggtheme = custom_theme(),
                          legend = c(0.8,0.85), # 
                          legend.title = "Risktype",legend.labs=c('High','Low'))
GSE32062.km.OS=mg_merge_plot(GSE32062.km.OS$plot,GSE32062.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
GSE32062.km.OS

GSE32062.risktype.cli$Status=ifelse(GSE32062.risktype.cli$OS==0,'Alive','Death')
GSE32062.model.fig=my_riskplot(cli_dat = GSE32062.risktype.cli,cols = risktype.col,xlab = 'sample',
                           b.labs = 'Survival time(months)',cutoff = median(GSE32062.risktype.cli$Riskscore),labs = c('H',''))
GSE32062.model.fig


fig2=mg_merge_plot(mg_merge_plot(tcga.lasso$plot,gene.coef.fig,widths = c(2,1),labels = c('','C')),
                   mg_merge_plot(gene.forest,tcga.model.fig,widths = c(1.5,1),labels = c('D','')),
                   mg_merge_plot(tcga.km.OS,tcga.roc,labels = c('F','G')),
                   mg_merge_plot(GSE32062.model.fig,GSE32062.km.OS,GSE32062.roc,ncol=3,labels = c('','I','J')),
                   nrow=4)
ggsave('results/02.model/Fig2.pdf',fig2,height = 20,width = 16)




save.image(file = 'project.RData')

#03.#########
dir.create('results/03.nomogram')
########
head(tcga.risktype.cli)
tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)
table(tcga_cox_datas$Age1)

table(tcga_cox_datas$Stage)
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='I'|tcga_cox_datas$Stage=='II']<-'I+II'
tcga_cox_datas$Stage[tcga_cox_datas$Stage=='III'|tcga_cox_datas$Stage=='IV']<-'III+IV'

table(tcga_cox_datas$Grade)
tcga_cox_datas$Grade[tcga_cox_datas$Grade=='G1'|tcga_cox_datas$Grade=='G2']<-'G1+G2'
tcga_cox_datas$Grade[tcga_cox_datas$Grade=='G3'|tcga_cox_datas$Grade=='G4']<-'G3+G4'

#
#####################
head(tcga_cox_datas)
table(tcga_cox_datas$Age1)
cli.km1=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga_cox_datas[which(tcga_cox_datas$Age1=='<=59'),]),
                   data=tcga_cox_datas[which(tcga_cox_datas$Age1=='<=59'),],
                   conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                   title='Age<=59',ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "Risktype",
                   legend.labs = c("High","Low"))
cli.km1
cli.km1=mg_merge_plot(cli.km1$plot,cli.km1$table,nrow=2,heights = c(2.5,1),align = 'v')

cli.km2=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga_cox_datas[which(tcga_cox_datas$Age1=='>59'),]),
                   data=tcga_cox_datas[which(tcga_cox_datas$Age1=='>59'),],
                   conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                   title='Age>59',ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "Risktype",
                   legend.labs = c("High","Low"))
cli.km2
cli.km2=mg_merge_plot(cli.km2$plot,cli.km2$table,nrow=2,heights = c(2.5,1),align = 'v')

table(tcga_cox_datas$Stage)
cli.km3=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga_cox_datas[which(tcga_cox_datas$Stage=='I+II'),]),
                   data=tcga_cox_datas[which(tcga_cox_datas$Stage=='I+II'),],
                   conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                   title='Stage I+II',ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "Risktype",
                   legend.labs = c("High","Low"))
cli.km3
cli.km3=mg_merge_plot(cli.km3$plot,cli.km3$table,nrow=2,heights = c(2.5,1),align = 'v')

cli.km4=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga_cox_datas[which(tcga_cox_datas$Stage=='III+IV'),]),
                   data=tcga_cox_datas[which(tcga_cox_datas$Stage=='III+IV'),],
                   conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                   title='Stage III+IV',ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "Risktype",
                   legend.labs = c("High","Low"))
cli.km4
cli.km4=mg_merge_plot(cli.km4$plot,cli.km4$table,nrow=2,heights = c(2.5,1),align = 'v')


table(tcga_cox_datas$Grade)
cli.km5=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga_cox_datas[which(tcga_cox_datas$Grade=='G1+G2'),]),
                   data=tcga_cox_datas[which(tcga_cox_datas$Grade=='G1+G2'),],
                   conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                   title='G1+G2',ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "Risktype",
                   legend.labs = c("High","Low"))
cli.km5
cli.km5=mg_merge_plot(cli.km5$plot,cli.km5$table,nrow=2,heights = c(2.5,1),align = 'v')

cli.km6=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga_cox_datas[which(tcga_cox_datas$Grade=='G3+G4'),]),
                   data=tcga_cox_datas[which(tcga_cox_datas$Grade=='G3+G4'),],
                   conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,
                   title='G3+G4 ',ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "Risktype",
                   legend.labs = c("High","Low"))
cli.km6
cli.km6=mg_merge_plot(cli.km6$plot,cli.km6$table,nrow=2,heights = c(2.5,1),align = 'v')

pdf('results/03.nomogram/Fig3.pdf',height = 10,width = 15)
mg_merge_plot(cli.km1,cli.km2,cli.km3,cli.km4,cli.km5,cli.km6,nrow=2,ncol=3,labels =LETTERS[1:6])
dev.off()

#######
#Age
tcga_cox_datas=crbind2DataFrame(tcga_cox_datas)
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age1,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat



#Grade
Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat

#Grade
Grade_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Grade,
                               data=tcga_cox_datas))
Grade_sig_cox_dat <- data.frame(Names=rownames(Grade_sig_cox[[8]]),
                                HR = round(Grade_sig_cox[[7]][,2],3),
                                lower.95 = round(Grade_sig_cox[[8]][,3],3),
                                upper.95 = round(Grade_sig_cox[[8]][,4],3),
                                p.value=round(Grade_sig_cox[[7]][,5],3))
Grade_sig_cox_dat


#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=tcga_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Stage_sig_cox_dat,
                     Grade_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Features=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
rownames(data.sig) <- c("Age","Stage",'Grade',"RiskScore")
data.sig$Features=rownames(data.sig) 
data.sig
data.sig$p.value=ifelse(data.sig$p.value<0.001,'<0.001',data.sig$p.value)
pdf('results/03.nomogram/Univariate.pdf',height = 4,width = 6,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap =5,lineheight = 9,zero = 1,
                 boxsize = .2,lwd.zero=1,lwd.ci=1.5,lwd.xaxis=1,
                 box_col='skyblue',summary_col="black",lines_col='black',zero_col='grey',
                 xlab='Hazard Ratio',lty.ci = 6,graph.pos =4)
dev.off()

write.csv(data.sig,'results/03.nomogram/Univariate analysis.csv')

#########
#T.stage +  N.stage +M.stage
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age+Riskscore, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Features=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
data.muti
rownames(data.muti) <- c("Age","RiskScore")
data.muti$Features=rownames(data.muti)
data.muti
data.muti$p.value=ifelse(data.muti$p.value<0.001,'<0.001',data.muti$p.value)
pdf('results/03.nomogram/Multivariate.pdf',height = 4,width =6,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap =5,lineheight = 9,zero = 1,
                 boxsize = .2,lwd.zero=1,lwd.ci=1.5,lwd.xaxis=1,
                 box_col='skyblue',summary_col="black",lines_col='black',zero_col='grey',
                 xlab='Hazard Ratio',lty.ci = 6,graph.pos =4)
dev.off()
write.csv(data.muti,'results/03.nomogram/Multivariate analysis.csv')

# # ##
# 
# pdf('results/05.nomogram/nomogram.pdf', width = 12, height = 10)
# nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$Riskscore,
#                                 Age=tcga_cox_datas$Age,
#                                 pathologic_M=tcga_cox_datas$pathologic_M,
#                                 Grade=tcga_cox_datas$Grade),
#                      os = tcga_cox_datas$OS.time,
#                      status = tcga_cox_datas$OS,
#                      mks = c(1,3,5)
# )
# dev.off()
# mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,3*365,5*365))
# 

#04.###############
dir.create('results/04.TME')
#######

                    check.names = F,row.names = 1)
tcga.TIMER[1:5,1:5]
tcga.TIMER=tcga.TIMER[tcga.risktype.cli$Samples,]
dim(tcga.TIMER)
table(str_split_fixed(colnames(tcga.TIMER),'_',2)[,2])


tme.df2=tcga.TIMER[tcga.risktype.cli$Samples,str_split_fixed(colnames(tcga.TIMER),'_',2)[,2]=='CIBERSORT']
colnames(tme.df2)=gsub('_CIBERSORT','',colnames(tme.df2))
tme.df2$Risktype=tcga.risktype.cli$Risktype
tme.df2=melt(tme.df2)
head(tme.df2)
ggplot(tme.df2,aes(x=variable,y=value,fill=Risktype))+
  geom_boxplot()+stat_compare_means(aes(group=Risktype), label = "p.signif", method = 'wilcox.test')+
  scale_fill_manual(values =risktype.col)+
  xlab('')+ylab('Fraction')+ggtitle('CIBERSORT')+
  theme_bw()+theme(text = element_text(family = 'Times',size = 12),legend.position = 'top',
                   axis.text.x = element_text(color = "black", size = 12,angle = 30,hjust = 1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
ggsave('results/04.TME/Fig5a.pdf',height = 5,width = 15)

#####
load('results/04.TME/TCGA_gsva_score.RData')
tcga.immu.ssgsea[1:5,1:5]
tcga.immu.ssgsea=tcga.immu.ssgsea[tcga.risktype.cli$Samples,]
tcga.immu.ssgsea=crbind2DataFrame(tcga.immu.ssgsea)
tcga.immu.ssgsea$Risktype=tcga.risktype.cli$Risktype
tcga.immu.ssgsea.df=melt(tcga.immu.ssgsea)
head(tcga.immu.ssgsea.df)
ggplot(tcga.immu.ssgsea.df,aes(x=variable,y=value,fill=Risktype))+
  geom_boxplot()+stat_compare_means(aes(group=Risktype), label = "p.signif", method = 'wilcox.test')+
  scale_fill_manual(values =risktype.col)+
  xlab('')+ylab('Score')+
  theme_bw()+theme(text = element_text(family = 'Times',size = 12),legend.position = 'none',
                   axis.text.x = element_text(color = "black", size = 12,angle = 30,hjust = 1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
ggsave('results/04.TME/Fig5b.pdf',height = 5,width = 15)

#####
tme.df3=tcga.TIMER[tcga.risktype.cli$Samples,str_split_fixed(colnames(tcga.TIMER),'_',2)[,2]=='MCPCOUNTER']
colnames(tme.df3)=gsub('_MCPCOUNTER','',colnames(tme.df3))
library(ggcorrplot)
library(psych)
mcpcounter_RS_cor <- corr.test(x =tcga.risktype.cli$Riskscore,
                               y = tme.df3[tcga.risktype.cli$Samples,],
                               method = "spearman",adjust = "BH",ci = F)

mcpcounter_RS_cor_res=data.frame(Immune_cell=colnames(tme.df3))
mcpcounter_RS_cor_res$cor<-as.numeric(mcpcounter_RS_cor$r)
mcpcounter_RS_cor_res$p.adj<-as.numeric(mcpcounter_RS_cor$p.adj)
head(mcpcounter_RS_cor_res)
mcpcounter_RS_cor_res=mcpcounter_RS_cor_res[order(mcpcounter_RS_cor_res$cor),]
head(mcpcounter_RS_cor_res)

ggplot(data=mcpcounter_RS_cor_res,aes(x=cor,y=reorder(Immune_cell,cor), color = -log10(p.adj))) +
  geom_point(aes(size=abs(cor)),show.legend = F) +
  scale_color_gradient(low = "#009B9F", high = "#C75DAA")+
  geom_segment(aes(yend=Immune_cell,xend=0),size=.5) +
  labs(x='spearman Correlation',y='Immune cell')+theme_bw()+
  theme(text = element_text(family = 'Times',size = 14))
ggsave('results/04.TME/Fig5c.pdf',height = 6,width = 7)



#####
tme.df5=data.frame(t(tcga.exp[c('CD274','PDCD1','CTLA4','HAVCR2','LGALS9','TIGIT','LAG3',names(lan)),tcga.risktype.cli$Samples]))
cor_res <- Hmisc::rcorr(as.matrix(tme.df5),type = 'pearson')
cor_res$P[is.na(cor_res$P)] <- 0
library(corrplot)
cor_res$P[1:5,1:5]
pdf('results/04.TME/Fig5d.pdf',height = 6,width = 7)
corrplot(as.matrix(cor_res$r[1:7,-c(1:7)]),
         p.mat = as.matrix(cor_res$P[1:7,-c(1:7)]),
         mar = c(0,0,1,1),
         col=colorRampPalette(c('purple', 'white','orange'))(100),
         tl.srt = 90,tl.cex = 1,tl.col = 'black',tl.offset = 0.5,
         cl.pos = c("b","r","n")[1],cl.align.text = 'l',cl.length = 5,
         cl.ratio = 0.1,cl.cex = 0.8,
         addgrid.col = 'white',
         method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")[6],
         insig = 'label_sig',
         sig.level=c(0.001,0.01,0.05),
         pch.cex=1,is.corr=T,xpd=T)
dev.off()

#05.#############
dir.create('results/05.TIDE_GSEA')
##TIDE
tcga_tide_res<-read.csv('results/05.TIDE_GSEA/tcga.tide.res.csv',row.names = 1,stringsAsFactors = F)
head(tcga_tide_res)

tme.TIDE=data.frame(tcga_tide_res[tcga.risktype.cli$Samples,c('TIDE','Responder')],tcga.risktype.cli)
head(tme.TIDE)
tme.TIDE %>%
  ggplot(aes(x=Risktype,y=TIDE,fill=Risktype))+
  geom_boxplot(notch = T)+ stat_compare_means(aes(group=Risktype), label = "p.format", method = 'wilcox.test')+
  scale_fill_manual(values =risktype.col)+theme_bw()+
  theme(text = element_text(family = 'Times',size = 15),legend.position = 'none',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
ggsave('results/05.TIDE_GSEA/Fig6a.pdf',height = 4,width =3.5)


#####IPS
library(IOBR)
IPS.res=IPS_calculation(project = 'OV',eset = tcga.exp,plot = F)
head(IPS.res)
dim(IPS.res)

tme.IPS=data.frame(IPS.res[tcga.risktype.cli$Samples,'IPS',drop=F],tcga.risktype.cli)
head(tme.IPS)
tme.IPS %>%
  ggplot(aes(x=Risktype,y=IPS,fill=Risktype))+
  geom_violin()+  
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = NA,fill="white")+
  stat_compare_means(aes(group=Risktype), label = "p.format", method = 'wilcox.test')+
  scale_fill_manual(values =risktype.col)+theme_bw()+ylab('Immune Phenotype Score')+
  theme(text = element_text(family = 'Times',size = 15),legend.position = 'none',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
ggsave('results/05.TIDE_GSEA/Fig6b.pdf',height = 4,width =3.5)



#06.GSEA###############
tcga.geneList.risktype=getGeneFC(gene.exp=tcga.exp[,tcga.risktype.cli$Samples], 
                                 group=tcga.risktype.cli$Risktype,ulab='High',dlab='Low')

set.seed(777)
tcga.geneList.risktype.gsea<-GSEA(tcga.geneList.risktype,TERM2GENE = h.all.gmt,seed=T)
tcga.geneList.risktype.gsea.res=tcga.geneList.risktype.gsea@result
write.xlsx(tcga.geneList.risktype.gsea.res,'results/05.TIDE_GSEA/TCGA.risktypeGSEA.res.xlsx',overwrite = T)

library(dotplotGsea)
risktype.gsea.dotplot=dotplotGsea(data = tcga.geneList.risktype.gsea,order.by = 'NES')

risktype.gsea.dotplot$plot+theme(text = element_text(family = 'Times'))+
  ggtitle('High risk VS Low risk')
ggsave('results/05.TIDE_GSEA/Fig6c.pdf',height = 8,width =10)



##
library(ggcorrplot)
library(psych)
HALLMARK_RS_cor <- corr.test(x =tcga.risktype.cli$Riskscore,
                             y = t(tcga.hall.ssGSEA[,tcga.risktype.cli$Samples]),
                             method = "spearman",adjust = "BH",ci = F)


HALLMARK_RS_cor_res=data.frame(pathway=rownames(tcga.hall.ssGSEA))
HALLMARK_RS_cor_res$cor<-as.numeric(HALLMARK_RS_cor$r)
HALLMARK_RS_cor_res$p.adj<-as.numeric(HALLMARK_RS_cor$p.adj)
head(HALLMARK_RS_cor_res)
table(HALLMARK_RS_cor_res$p.adj<0.05)
HALLMARK_RS_cor_res=HALLMARK_RS_cor_res[HALLMARK_RS_cor_res$p.adj<0.05,]
HALLMARK_RS_cor_res=HALLMARK_RS_cor_res[order(HALLMARK_RS_cor_res$cor),]
head(HALLMARK_RS_cor_res)
pdf('results/04.Immune_treatment/RS_cor_HALLMARK.pdf',height = 5,width = 6)
library(rcartocolor)
ggplot(data=HALLMARK_RS_cor_res,aes(x=cor,y=reorder(pathway,cor), color = -log10(p.adj))) +
  geom_point(aes(size=abs(cor)),show.legend = F) +
  scale_color_continuous(type = "gradient")+
  scale_color_gradient(low = "skyblue", high = "pink")+
  geom_segment(aes(yend=pathway,xend=0),size=.5) +
  labs(x='spearman Correlation',y='')+theme_bw()+
  theme(text = element_text(family = 'Times',size = 14))+
  scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))
dev.off()
ggsave('results/05.TIDE_GSEA/Fig6d.pdf',height = 8,width =8.5)


#07.#############
dir.create('results/06.Mutation')
# #####
#######
tcga.maf= getTCGAMAFByCode('OV')
tcga.risktype.use=tcga.risktype.cli[,c('Samples','Risktype')]
table(tcga.risktype.use$Risktype)
head(tcga.risktype.use)
colnames(tcga.risktype.use)[1]='Tumor_Sample_Barcode'
tcga.risktype.use$Tumor_Sample_Barcode=substr(tcga.risktype.use$Tumor_Sample_Barcode,1,12)
tcga.risktype.high=tcga.risktype.use[which(tcga.risktype.use$Risktype=='High'),]
tcga.risktype.low=tcga.risktype.use[which(tcga.risktype.use$Risktype=='Low'),]

write.table(tcga.risktype.high,file='results/06.Mutation/tcga.risktype.high.txt')
write.table(tcga.risktype.low ,file='results/06.Mutation/tcga.risktype.low.txt')


tcga.maf1=subsetMaf(tcga.maf,tsb=tcga.maf@data$Tumor_Sample_Barcode[substr(tcga.maf@data$Tumor_Sample_Barcode,1,12)%in%
                                                                      tcga.risktype.high$Tumor_Sample_Barcode])
tcga.maf1<-read.maf(tcga.maf1@data,isTCGA=T,clinicalData = 'results/06.Mutation/tcga.risktype.high.txt')
tcga.maf1@clinical.data

tcga.maf2=subsetMaf(tcga.maf,tsb=tcga.maf@data$Tumor_Sample_Barcode[substr(tcga.maf@data$Tumor_Sample_Barcode,1,12)%in%
                                                                      tcga.risktype.low$Tumor_Sample_Barcode])
tcga.maf2<-read.maf(tcga.maf2@data,isTCGA=T,clinicalData = 'results/06.Mutation/tcga.risktype.low.txt')
tcga.maf2@clinical.data



######
mdf_cp=mafCompare(tcga.maf1,tcga.maf2,m1Name = 'High',m2Name = 'Low')
maf.compare.res=mdf_cp$results

head(maf.compare.res)
table(maf.compare.res$pval<0.05)
maf.compare.res.filter=maf.compare.res[maf.compare.res$pval<0.05,]
head(maf.compare.res.filter)
write.csv(maf.compare.res.filter,'results/06.Mutation/maf_compare.csv')



pdf('results/06.Mutation/diff_mutation.pdf',height = 7,width = 5,onefile = F)
forestPlot(mafCompareRes = mdf_cp,color = c('royalblue', 'maroon'),pVal = 0.05)
dev.off()

save.image(file = 'project.RData')
