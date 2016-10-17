library(gridExtra)
library(ggplot2)
library(reshape2)


# Importer une matrice de comptage
# lignes = sample / colonnes = noeuds taxonomiques
# doit posséder un entête de colonne ainsi que 
# les noms des sample dans la première colonne
ImportCountMatrix <- function(file) {
  Matrix=read.csv(file = file,sep=",",header=T)
  rownames(Matrix)=Matrix$SampleID
  Matrix=Matrix[,-1]
  if(sum(grepl("Archaea",colnames(Matrix)))) {
    Matrix=Matrix[,-grep("Archaea",colnames(Matrix))]
  }
  return(Matrix)
}

# Filtrage de la matrice de comptage
# se base sur les répartition des Standard Deviation
FilterCountMatrix <- function(m) {
  SDval=c()
  for(i in 1:ncol(m)) { v=sd(m[,i]);   SDval[i]=v;}
  SDvalBOXPLOT=boxplot(SDval)
  m=m[,which(SDval>SDvalBOXPLOT$stats[4])]
  return(m)
}

# Normalisation de la mtrice de comptage
# doit trouver une colonne qui contient "Bacteria"
# puis tri selon les rownames
NormaliseMatrix <- function(m) {
  BactIndice=grep("Bacteria",colnames(m))
  Normed=m[,]/m[,BactIndice]
  Normed=Normed[order(rownames(Normed)),]
  return(Normed)
}

# Importer les MetaData des sample
# le nom de sample doit être dans la première colonne
# puis tri selon les rownames
ImportMetada <- function(f) {
  Metadata=read.csv(file = f,sep=",",header=T)
  rownames(Metadata)=Metadata$Sample
  Metadata=Metadata[,-1]
  Metadata=Metadata[order(rownames(Metadata)),]
}


# Test statistique (t.test & Wilcox)
# pour chaque niveau taxonomique on filtre les colonnes de ce niveau
# puis correction des p-valeurs (BH) 
MakeTestByLevel <- function(matrice,liste1,liste2) {
  Global_results=  data.frame()
  LesLevel=c("P","C","O","F","G","S")
  for(l in 1:length(LesLevel)) {
    LesCol=grep(paste0("^",LesLevel[l],"."),colnames(matrice))
    LevelMatrice=matrice[,LesCol]
    All_results_T=data.frame()
    All_results_W=data.frame()
    for(i in 1:ncol(LevelMatrice)) {
      res=t.test(LevelMatrice[liste1,i],LevelMatrice[liste2,i])
      test_result = data.frame(pval=res$p.value,test_name="student","taxon"=colnames(LevelMatrice)[i],"level"=LesLevel[l])
      All_results_T=rbind(All_results_T,test_result)
      res=wilcox.test(LevelMatrice[liste1,i],LevelMatrice[liste2,i])
      test_result= data.frame(pval=res$p.value,test_name="wilcox","taxon"=colnames(LevelMatrice)[i],"level"=LesLevel[l])
      All_results_W=rbind(All_results_W,test_result)
    }
    All_results_T$adj.p.val = p.adjust(All_results_T$pval,method="BH")
    All_results_W$adj.p.val = p.adjust(All_results_W$pval,method="BH")
    All_results_TW=rbind(All_results_T,All_results_W)
    Global_results=rbind(Global_results,All_results_TW)
  }
  return(Global_results)
}

# Test statistique (t.test & Wilcox)
# pas de filtrage des colonnes selon le niveau taxonomique
# puis correction des p-valeurs (BH) 
MakeTestAllLevel <- function(matrice,liste1,liste2) {
  Global_results=  data.frame()
  All_results_T=data.frame()
  All_results_W=data.frame()
  LesLevel=c("P","C","O","F","G","S")
  
  for(l in 1:length(LesLevel)) {
    LesCol=grep(paste0("^",LesLevel[l],"."),colnames(matrice))
    LevelMatrice=matrice[,LesCol]
    for(i in 1:ncol(LevelMatrice)) {
      res=t.test(LevelMatrice[liste1,i],LevelMatrice[liste2,i])
      test_result = data.frame(pval=res$p.value,test_name="student","taxon"=colnames(LevelMatrice)[i],"level"=LesLevel[l])
      All_results_T=rbind(All_results_T,test_result)
      res=wilcox.test(LevelMatrice[liste1,i],LevelMatrice[liste2,i])
      test_result= data.frame(pval=res$p.value,test_name="wilcox","taxon"=colnames(LevelMatrice)[i],"level"=LesLevel[l])
      All_results_W=rbind(All_results_W,test_result)
    }
  }
  
  All_results_T$adj.p.val = p.adjust(All_results_T$pval,method="BH")
  All_results_W$adj.p.val = p.adjust(All_results_W$pval,method="BH")
  All_results_TW=rbind(All_results_T,All_results_W)
  return(All_results_TW)
}


# Analyse graphique globalede la matrice de comptage 
# pour un niveau taxonomique spécifié
# facetage selon une ou deux metadonnées (paramètres info et info2)
# uniquement graphique, pas de calcul de p-valeur
GlobalAnalysis <- function(matrice,level,info,info2) {
  matrice$Sample=rownames(matrice)
  pat=paste("^",level,".",sep="")
  if(!missing(info2)) {
    matrice=cbind(matrice,as.character(info))
    matrice=cbind(matrice,as.character(info2))
    matriceM=melt(matrice)
    colnames(matriceM)=c("Sample","Info","Info2","Node","Value")
    ggplot(matriceM[grep(pat,matriceM$Node),],aes(x=Node,y=Value,color=Info,shape=Info))+geom_point(size=3)+geom_boxplot()+coord_flip()+facet_wrap(~Info2)
  }
  else {
    if(!missing(info)) {
      matrice=cbind(matrice,as.character(info))
      matriceM=melt(matrice)
      colnames(matriceM)=c("Sample","Info","Node","Value")
      ggplot(matriceM[grep(pat,matriceM$Node),],aes(x=Node,y=Value,color=Info,shape=Info))+geom_point(size=3)+geom_boxplot()+coord_flip()
    }
    else {
      matriceM=melt(matrice)
      colnames(matriceM)=c("Sample","Node","Value")
      ggplot(matriceM[grep(pat,matriceM$Node),],aes(x=Node,y=Value))+geom_point(size=3)+geom_boxplot()+coord_flip()
    }
  }
}


# Boxplot 
# compare deux conditions pour un noeud donné (colonne)
# les conditions sont nommées & les liste d'échantillons sont fournies
BoxplotNode <- function(matrice,node,list1,list2,info1,info2) {
  m1=matrice[list1,c(node)]
  m2=matrice[list2,c(node)]
  df=rbind(data.frame("val"=m1,"info"=info1),data.frame("val"=m2,"info"=info2))
  g=ggplot(df,aes(x=info,y=val,color=info))+geom_boxplot()+scale_x_discrete(name="")+scale_y_continuous(name=node)+theme(legend.position="none")
  return(g)
}

