#good place to put your top-level package documentation

.onAttach <- function (lib, pkgname="midaR") {
    ## Put stuff here you want to run when your package is loaded
    invisible()
}

# Convert unlabled data into vector
make.corrMat<-function(unlabeled_vector){
  l=length(unlabeled_vector)
  #Convert vector into matrix
  unlabeled<-matrix(1:l,nrow=l,ncol=l,byrow=TRUE)
  unlabeled
  for(i in 1:dim(unlabeled)[1])
  {
    for (j in 1:dim(unlabeled)[2])
    {
      unlabeled[i,j]<-unlabeled_vector[abs(j-i)+1]
    }
  }
  #Convert into triangular matrix
  unlabeled[lower.tri(unlabeled)] <- 0
  return(unlabeled)
}

# Get isotope name (M+1, M+2 ...) from compound names
getMs<-function(inpt,compd_){
  Ms=unlist(lapply(as.character(inpt),function(n){
    m=strsplit(strsplit(n,compd_)[[1]][2],"Results")[[1]][1]
    r=NULL
    if (m=="."){r="M+0"}
    if (m!="."){r=gsub("M","M+",gsub("\\.","",m))}
    return(r)}))
  return(Ms)
}

# Get name for replicate set (biological replicates prefix)
getRepSetNames<-function(rnames){
  unique(unlist(lapply(rnames,function(n){
    v=strsplit(n,"_")[[1]]
    r=paste(v[1:(length(v)-1)],collapse="_")
    return(r)})))
}

# Takes a matrix and generate a barplot
ggBar<-function(mtx,colz,main_="", ylab_=""){
  df<-melt(mtx)
  compd_n<-strsplit(as.character(df$X2)[[1]],".Results")[[1]]
  df$Isotope=getMs(df$X2,compd_n)
  df$Isotope=factor(df$Isotope,levels=rev(unique(df$Isotope)))
  df$X1=factor(df$X1,levels=rownames(mtx))
  options(scipen=10000)
  g=ggplot(df, aes(fill=Isotope, y=value, x=X1)) +
    geom_bar(position="stack", stat="identity",colour="black",size=.3) +
    scale_fill_manual(name="Isotope",values=colz)+
    ggtitle(paste0(compd_n," \n(",main_,")")) +
    theme(plot.title=element_text(face="bold",hjust=0.5),axis.text.x = element_text(angle = 90, vjust = .5,colour = "black"),axis.text.y = element_text(colour = "black"),)+theme(plot.title=element_text(face="bold",hjust=0.5),panel.border = element_blank(),panel.background =element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line=element_line(colour="black"))+
    xlab("Samples") +
    ylab(ylab_)
  return(g)
}

ggBar.sde<-function(mtx,sderr,colz,main_="", ylab_=""){
  df<-melt(mtx)
  compd_n<-strsplit(as.character(df$X2)[[1]],".Results")[[1]]
  df$Isotope=getMs(df$X2,compd_n)
  df$Isotope=factor(df$Isotope,levels=rev(unique(df$Isotope)))
  df$X1=factor(df$X1,levels=rownames(mtx))
  options(scipen=10000)
  g=ggplot(df, aes(fill=Isotope, y=value, x=X1)) +
    geom_bar(position="stack", stat="identity",colour="black",size=.3) +
    scale_fill_manual(name="Isotope",values=colz)+
    ggtitle(paste0(compd_n," \n(",main_,")")) +
    theme(plot.title=element_text(face="bold",hjust=0.5),axis.text.x = element_text(angle = 90, vjust = .5,colour = "black"),axis.text.y = element_text(colour = "black"),)+theme(plot.title=element_text(face="bold",hjust=0.5),panel.border = element_blank(),panel.background =element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line=element_line(colour="black"))+
    xlab("Samples") +
    ylab(ylab_)

  for (i in 1:dim(mtx)[1]){
    nm=rownames(mtx)[i]
    yval=sum(df$value[which(df$X1==nm)])
    g <-g +geom_segment(aes_string(x = i, y =yval , xend = i, yend = yval+sderr[i]),size=.3)
    g <-g +geom_segment(aes_string(x = i-.1, y =yval+sderr[i] , xend = i+.1, yend = yval+sderr[i]),size=.3)
  }

  return(g)
}

# Computing the standard error from pool size
get.sde<-function(df,mt){
  rnames.st=rownames(mt);cnames=colnames(mt)
  sub.df=df[,cnames]
  sde<-numeric()
  for (rnm in rnames.st) {
    ix=which(startsWith(rownames(sub.df),rnm))
    v=rowSums(data.matrix(sub.df[ix,]))

    if (length(v)==1){  # if we have one biological replicate
      sde=c(sde,0)
    }else{  # if we have several biological replicates
      sde=c(sde, sd(v)/sqrt(length(v)))
    }


  }
  return(sde)
}

# QC remove the compound with Fractional Enrichment <95 or >105.
M0_FE_QC<-function(df,compoundz,L_){
  for (i in 1:length(compoundz)){
    cix=which(startsWith(tolower(colnames(df)),tolower(compoundz[i])))
    m=as.numeric(df[L_[["unlabeled"]],cix[1]])
    #m=rowSums(data.matrix(df[L_[["unlabeled"]],cix]))
    if (length(which(m < 95 | m > 105))/length(m)>=2/3){
      df[,cix]<-NA
    }
  }
  return(df)
}

# Imputation of missing values using the minimum detection values per compound isotop per group
# Exclude compounds below baseline detection (1000):
#     if ( 2/3 of unlabeled samples M+0 are below 1000) or all labeled sample entries are < 1000
imputMin<-function(df,compoundz,baseLineD){    ### This routine expects the order of isotopes for a particular compounds as M+0, M+1, M+2, ... M+n
  #df=InputData;compoundz=compounds;baseLineD=baseLineDet
  snames=unique(df$Labeling)

  for (i in 1:length(compoundz)){

    cix=which(startsWith(tolower(colnames(df)),tolower(compoundz[i])))
    m=data.matrix(df[,cix])

    #rownames(m)=df$Labeling
    uns=which(tolower(df$Labeling)=="unlabeled")

    #Check if compound is bellow baslize detection
    #if (sum(m[uns,1])==0){
    if (length(which(m[uns,1]<baseLineD))/length(uns) > 2/3  ){
      cat(paste0(compoundz[i]," : falls bellow detection baseline ",baseLineD," for M+0 unlabeled sample. It will be removed\n"))
      df[,cix]=NA
      next

    }else if ( length(which(m[-uns,] >= baseLineD))==0){
      cat(paste0(compoundz[i]," : falls bellow baseline detection ",baseLineD," in all labeled samples. It will be removed\n"))
      df[,cix]=NA
      next

    }
    else{
      for (j in 1:length(snames)){
        rix=which(startsWith(tolower(df$Labeling),tolower(snames[j])))
        m1=m[rix,]
        if(sum(m1)==0){
          cat(compoundz[i]," x ",snames[j]," : all entries are zero thus are skipped\n")
          next
        }else{

          #min accross sample for individual isotopologue

          if (is.matrix(m1)){
            m2=apply(m1,2,function(cl){
              ret=c()
              if (all(cl==0)){
                ret=cl
              }else
              {
                min_=min(cl[which(cl>0)])
                cl[which(cl==0)]<- min_
                ret=cl
              }
              return(ret)
            })
            df[rix,cix]<-m2
          }
        }
      }
    }
  }
  return(df)
}

# Removing compound with low tectection
removeCompNA<-function(df){
  compdNA=which(colSums(is.na(df)) > 0)
  if (length(compdNA)>0){
    df <- df[ ,-compdNA]
    message(paste0("\t\t\t",length(compdNA)," colunm(s) removed"))
  }else{
    message("\t\t\tNo colunm removed")
  }
  return(df)
}

# Adjasting compound List
adj.CompList<-function(df,comps){
  colnames_<-colnames(df)[3:dim(df)[2]]
  nLcomp<-c()
  for (comp in comps){
    if (length(which(startsWith(tolower(colnames_),tolower(comp))))>0){
      nLcomp=c(nLcomp,comp)
    }
  }
  d=setdiff(comps,nLcomp)
  message("\t\t\t", length(comps)-length(nLcomp)," compound(s) removed: ", paste0(d,collapse=", "))
  return(nLcomp)
}

# Adjasting compound Input Dataframe
adj.CompDF.initial<-function(df,compdsDF){
  comps=compdsDF$Compound[compdsDF$Included==1]
  idx=unlist(lapply(comps,function(cmp){
    r= which(startsWith(tolower(colnames(df)),tolower(cmp)))
    return(r)}))
  df=df[,c(c(1:4),idx)]
  return(df)
}

# Submitutes NAs with Zeros
sub_Na_Zero<-function(df_){
  for (i in 1: nrow(df_)){
    for (j in 2:ncol(df_)){
      if (is.na(df_[i,j])){df_[i,j]="0"}
      if (df_[i,j]==""){df_[i,j]="0"}
    }
  }
  return(df_)
}

# Submitutes negative values with Zeros
sub_Ng_Zero<-function(df){
  for (i in 1:nrow(df)){
    for (j in 2:ncol(df)){
      if (df[i,j]< 0){df[i,j]=0}
    }
  }
  return(df)
}

# Calculates coefficient of variation
cv=function(X){r=(sd(X,na.rm=TRUE)/mean(X,na.rm=TRUE))*100;return(r)}

# Takes input or ouput data table and generate the coeff of variation (sum all isotopologues and computing cv)
coev_df<-function(df,compoundz){
  sNames=getRepSetNames(rownames(df))
  M.cv=matrix(NA,nrow=length(sNames), ncol=length(compoundz), dimnames=list(sNames, compoundz))
  for (i in 1:length(sNames)){
    r.ix=which(startsWith(rownames(df),sNames[i]))

    for (j in 1:length(compoundz)){
      c.ix=which(startsWith( colnames(df), compoundz[j]))
      M.cv[i,j-2]=cv(rowSums(data.matrix(df[r.ix,c.ix])))
    }
  }
  return(M.cv)
}

# Takes input or ouput data table and generate the coeff of variation
coev_df.old<-function(df){
  sNames=unique(as.character(df$Labeling ) )
  M.cv=matrix(NA,nrow=length(sNames),ncol=dim(df)[2] - 2); rownames(M.cv)=sNames; colnames(M.cv)=colnames(df)[3:dim(df)[2]]
  for (i in 1:length(sNames)){
    r.ix=which(df$Labeling==sNames[i])
    for (j in 3:dim(df)[2]){
      M.cv[i,j-2]=cv(as.numeric(df[r.ix,j]))
    }
  }
  return(M.cv)
}


# Draws Arrows
arrEnds<-function(LisFrom,LisTo,sideFrom=0,sideTo=0){
  #LisFrom=LFig[["Suc"]]; LisTo=LFig[["AKG"]]; sideFrom=4;sideTo=1

  if (sideFrom==0 |sideTo==0){stop("Please specify the sides")}
  vFrom=c()
  if (sideFrom==1){ vFrom=c(mean(c(LisFrom[c(1,3)])),LisFrom[2])}
  if (sideFrom==2){ vFrom=c(LisFrom[1],mean(LisFrom[c(2,4)]))}
  if (sideFrom==3){ vFrom=c( mean(LisFrom[c(1,3)]),LisFrom[4])}
  if (sideFrom==4){ vFrom=c( LisFrom[3],  mean(LisFrom[c(2,4)]))}
  vTo=c()
  if (sideTo==1){ vTo=c(mean(c(LisTo[c(1,3)])),LisTo[2])}
  if (sideTo==2){ vTo=c(LisTo[1],mean(LisTo[c(2,4)]))}
  if (sideTo==3){ vTo=c( mean(LisTo[c(1,3)]),LisTo[4])}
  if (sideTo==4){ vTo=c( LisTo[3],  mean(LisTo[c(2,4)]))}
  v=c(vFrom,vTo)
  return(v)
}

# Makes QC plot based on CV
cvPlt<-function(df,main=""){
  #s.colz=c( '#ffe119', '#4363d8', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')
  d.m<-melt(df);names(d.m)=c("SampleName","Compound","cv")
  d.m$SampleName=factor(d.m$SampleName,levels=rownames(df))

  g<-ggplot(d.m, aes(x=Compound, y=cv, fill=SampleName)) +
    geom_bar(stat="identity", colour="black", position="dodge", size=0.25, width=0.8, alpha=0.8) +
    #scale_fill_manual(values=s.colz[1:dim(df)[1]])+
    scale_fill_manual(values=c("dodgerblue4","green","darkorchid4","magenta","brown","gray"))+
    labs(fill="Groups" )+
    ggtitle(main)+
    theme(plot.title=element_text(face="bold",hjust=0.5),axis.text.x = element_text(angle = 90, vjust = .5,colour = "black"),axis.text.y = element_text(colour = "black"),)+theme(plot.title=element_text(face="bold",hjust=0.5),panel.border = element_blank(),panel.background =element_blank(),panel.grid.major = element_line(color="gray",size=0.15),panel.grid.minor = element_blank(),axis.line=element_line(colour="black"))+
    coord_flip()
  return(g)
}

# Makes a figure for the TCA cycle using generated bargraph in pdf
tcaPlt<-function(inpath1){
  mrg=1
  arrowL=1.7
  spaceH=1
  spaceV=1
  fig.W=11
  fig.H=10

  height_=fig.H*7 + mrg*4 + arrowL*7  + spaceV*7;  height_
  width_=fig.W*5 + mrg*2 + arrowL*5  + spaceH*5;  width_

  #row1
  stY=0
  stX=(width_-fig.W)/2 ;                                 Suc=c(stX,stY,stX+fig.W,stY+fig.H);Suc
  stX=(width_-fig.W)/2 + fig.W*2 + spaceH*2 + arrowL*2 ; Gln=c(stX,stY,stX+fig.W,stY+fig.H);Gln
  #row2
  stY= fig.H + spaceH + arrowL
  stX=(width_-fig.W)/2 - (fig.W + spaceH + arrowL)     ; Fum=c(stX,stY,stX+fig.W,stY+fig.H);Fum
  stX=(width_-fig.W)/2 + fig.W + spaceH + arrowL       ; AKG=c(stX,stY,stX+fig.W,stY+fig.H);AKG
  stX=(width_-fig.W)/2 + fig.W*2 + spaceH*2 + arrowL*2 ; Glt=c(stX,stY,stX+fig.W,stY+fig.H);Glt
  #row3
  stY= fig.H*2 + spaceH*2 + arrowL*2
  stX=(width_-fig.W)/2 - (fig.W + spaceH*2 + arrowL)   ;  Mal=c(stX,stY,stX+fig.W,stY+fig.H);Mal
  stX=(width_-fig.W)/2 + fig.W + spaceH*2 + arrowL     ;  IsoC=c(stX,stY,stX+fig.W,stY+fig.H);IsoC
  #row4
  stY= fig.H*3 + spaceH*3 + arrowL*3 ;
  stX=(width_- (2*fig.W+ spaceH + arrowL ))/2        ;  OAA=c(stX,stY,stX+fig.W,stY+fig.H);OAA
  stX=(width_+ (spaceH + arrowL ))/2                 ;  Cit=c(stX,stY,stX+fig.W,stY+fig.H);Cit
  #row5
  stY= fig.H*4 + spaceH*2+  arrowL*3
  stX=(width_-fig.W)/2 - (fig.W*2 + spaceH + arrowL*2)    ; Asn=c(stX,stY,stX+fig.W,stY+fig.H);Asn
  stX=(width_-fig.W)/2 - (fig.W + spaceH + arrowL)    ; Asp=c(stX,stY,stX+fig.W,stY+fig.H);Asp
  stX=(width_-fig.W)/2 + fig.W*2 + spaceH + arrowL*2  ; Ita=c(stX,stY,stX+fig.W,stY+fig.H);Ita

  #row6
  stY= fig.H*5 + spaceH*5+  arrowL*4
  stX=(width_-fig.W)/2 - (fig.W + spaceH + arrowL)  ; Lac=c(stX,stY,stX+fig.W,stY+fig.H);Lac
  stX=(width_-fig.W)/2                              ; Pyr=c(stX,stY,stX+fig.W,stY+fig.H);Pyr
  stX=(width_-fig.W)/2 + fig.W + spaceH + arrowL    ; Ala=c(stX,stY,stX+fig.W,stY+fig.H);Ala
  #row7
  stY= fig.H*6 + spaceH*6+  arrowL*5
  stX=(width_-fig.W)/2                              ; Glu6P=c(stX,stY,stX+fig.W,stY+fig.H);Glu6P
  #row8
  stY= fig.H*7 + spaceH*7+  arrowL*6
  stX=(width_-fig.W)/2                              ; Glu=c(stX,stY,stX+fig.W,stY+fig.H);Glu


  #
  LFig=list(Glu,Glu6P,Pyr,Lac,Ala,Cit,Ita,IsoC,AKG,Glt,Gln,Suc,Fum,Mal,OAA,Asp,Asn)
  LTca<-list("Glucose","Glucose-6-P","Pyruvate","Lactate","L.Alanine","Citrate","Itaconate","Isocitrate","alpha.Ketoglutarate","L.Glutamate","L.Glutamine","Succinate","Fumarate","Malate",
             "Oxalacetate","Aspartic.Acid","L.Asparagine")
  names(LFig)<-names(LTca)<-c("Glu","Glu6P","Pyr","Lac","Ala","Cit","Ita","IsoC","AKG","Glt","Gln","Suc","Fum","Mal","OAA","Asp","Asn")

  par(mar=rep(mrg,4))
  plot(0,0, type="n",xaxt="n", yaxt="n", xlab="", ylab="",main="TCA cycle",
       xlim=c(0,width_), ylim=c(0,height_))

  for (i in 1:length(LFig)){
    AbrComp=names(LFig)[i]
    comp=LTca[[AbrComp]]

    if (file.exists(paste0(inpath1,"/",comp,"_TE.repset.pdf"))){
      im=image_read_pdf(paste0(inpath1,"/",comp,"_TE.repset.pdf"),density=300)
      rasterImage(im, LFig[[i]][1],LFig[[i]][2], LFig[[i]][3],LFig[[i]][4])
    }else{
      text(LFig[[i]][1]+fig.W/2,LFig[[i]][2]+fig.H/2,labels =comp,cex=.5)
    }
  }

  arr_col='blue'
  LArr<-list()
  LArr[[1]]<-c(arrEnds(LFig[["AKG"]],LFig[["Suc"]],1,4),.4)
  LArr[[2]]<-c(arrEnds(LFig[["Suc"]],LFig[["Fum"]],2,1),.4)
  LArr[[3]]<-c(arrEnds(LFig[["AKG"]],LFig[["Glt"]],4,2),0)
  LArr[[4]]<-c(arrEnds(LFig[["Glt"]],LFig[["Gln"]],1,3),0)
  LArr[[5]]<-c(arrEnds(LFig[["Fum"]],LFig[["Mal"]],3,1),.2)
  LArr[[6]]<-c(arrEnds(LFig[["Mal"]],LFig[["OAA"]],3,2),.4)
  LArr[[7]]<-c(arrEnds(LFig[["Cit"]],LFig[["IsoC"]],4,3),.4)
  LArr[[8]]<-c(arrEnds(LFig[["OAA"]],LFig[["Asp"]],2,1),.4)
  LArr[[9]]<-c(arrEnds(LFig[["OAA"]],LFig[["Cit"]],4,2),.2)
  LArr[[10]]<-c(arrEnds(LFig[["Cit"]],LFig[["Ita"]],4,2),.2)
  LArr[[11]]<-c(arrEnds(LFig[["Pyr"]],LFig[["Cit"]],1,3),-.4)
  LArr[[12]]<-c(arrEnds(LFig[["Pyr"]],LFig[["Ala"]],4,2),0)
  LArr[[13]]<-c(arrEnds(LFig[["Pyr"]],LFig[["Lac"]],2,4),0)
  LArr[[14]]<-c(arrEnds(LFig[["IsoC"]],LFig[["AKG"]],1,3),.1)
  LArr[[15]]<-c(arrEnds(LFig[["Glu6P"]],LFig[["Pyr"]],1,3),0)
  LArr[[16]]<-c(arrEnds(LFig[["Glu"]],LFig[["Glu6P"]],1,3),0)
  LArr[[17]]<-c(arrEnds(LFig[["Asp"]],LFig[["Asn"]],2,4),0)

  for (i in 1:length(LArr)){
    iArrows( LArr[[i]][1],LArr[[i]][2],LArr[[i]][3],LArr[[i]][4], h.lwd=1, sh.lwd=1, curve=LArr[[i]][5], sh.col=arr_col,size=.2,width=1.2)
  }
}

### Generating my own progress bar
get.prog<-function(prg_,LEN_){
  n=ceiling(prg_*100/LEN_)
  bar=paste0(rep("|",floor(n/10)),collapse = "")

  strng=paste0("|",bar,paste0(rep(".",10-nchar(bar)),collapse = "")," ",n,"%|")
  assign("prg", prg_+1,envir=parent.frame(n=1))
  return(strng)
}


##====================== MAIN ======================#

## MIDA's main Function.
mainFun<-function(inFiles__, inPath_, selectComp_, baseLineDet_=1000, overWrite_="Skip"){

  options(stringsAsFactors = F, scipen=999)
  #testing line
  #inFiles__= inFiles ; inPath_= inPath; selectComp_=selectComp; baseLineDet_=baseLineDet; overWrite_=overWrite

  inFiles_=c()
  if (tolower(overWrite_)=="skip"){
    for (inFileX in inFiles__){
      if (!file.exists(paste0(inPath_,"/",strsplit(inFileX,"_uncorr.csv")[[1]],"_out"))){
        inFiles_=c(inFiles_,inFileX)
      }
    }
    if (length(inFiles_)==0){
      return(message(paste0("\n\t\tAll files in the directory were already processed!      \n\n                *---      ", format(Sys.time(), "%a, %b %Y %d, %X"),"      ---*        \n\n"))  )
    }
  }else if ( tolower(overWrite_)=="overwrite") {
    inFiles_=inFiles__
  }

  message("#---------------------  Current settings  ----------------------------------------#")
  message("\tWorking directory: ", inPath_)
  message(paste0("\tNumber of uncorrected data files:  ",length(inFiles_)))
  message("\tCompound selection: ", selectComp_)
  message("\tDirectory overwrite: ", overWrite_)
  message("\tBaseline detection: ", baseLineDet_)
  message("#---------------------------------------------------------------------------------#\n")

  #*** run signal check-----
  #message(runSignal)
  #if (runSignal=="cancel"){ runSignal<<-"run"; message(runSignal); return(message(paste0("\n\t\t*---  MIDA run aborted !   ", format(Sys.time(), "%a, %b %Y %d, %X"),"    ---*       \n\n")))}
  #if (runSignal=="cancel"){ assign("runSignal","run",envir=parent.frame(n=2));message(runSignal);return(message(paste0("\n\t\t*---  MIDA run aborted !   ", format(Sys.time(), "%a, %b %Y %d, %X"),"    ---*       \n\n")))}

  LEN=1+(length(inFiles_)*9)
  prg=1

  for (inFile_ in inFiles_){

    message(paste0("\n\n*-------------  Correction for ", inFile_ ," ...  --------------"))

    prefix=strsplit(inFile_,"uncorr.csv")[[1]]
    outFile=paste0(prefix,"corr.csv")
    outHtmp="Heatmap.pdf"
    mtdFile=paste0( prefix,"mtd.csv")
    cmpdFile=paste0( prefix,"cmpd.csv")

    ## check files : metada and compound table
    if(!file.exists(paste0(inPath_,"/",mtdFile))) {
      message(paste0(" \t*--- warning :", inFile_, " has no metadata associated with it! -- it is skipped! ---*"))
      next
    }
    if(!file.exists(paste0(inPath_,"/",cmpdFile))) {
      message(paste0(" \t*--- warning :", inFile_, " has no compound associated with it! -- it is skipped! ---*"))
      next
    }

    ### Creating the output directories
    subDir1=paste0(inPath_,"/",strsplit(outFile,"_corr.csv")[[1]],"_out")
    subDir2=paste0(subDir1,"/extra")


    ### Overwrite current output or not
    if (file.exists(subDir1)){
      if (tolower(overWrite_)==("overwrite")){
        if (!file.exists(subDir2)){dir.create(subDir2)}
      }
    }else{
      dir.create(subDir1)
      dir.create(subDir2)
    }


    ################################
    # Pre-processing of input file
    #reading the tab- or comma-separatated values data and metadata

    message(paste0("\n",get.prog(prg,LEN)," *  1 --- Reading the compound abundance table"))
    InputData <- read.csv(paste0(inPath_,"/",inFile_),row.names=1,header=T)
    compoundsDF<-read.csv(paste0(inPath_,"/",cmpdFile),header=T)

    ##Initial ckeck - see whether the compounds in the compound list match those in the data table
    ##And select only usefull colunm { Sample Name 'Name', and all compounds isotopogue columns }

    ex=c(); ex.i=c()

    for (i in 1:length(compoundsDF$Compound)){
      cmpd=compoundsDF$Compound[i]
      L=which(startsWith(tolower(colnames(InputData)),tolower(cmpd)))
      ex.i<-c(ex.i,L)
      if (length(L)==0){ex=c(ex,0)}
      if (length(L)>0){ex=c(ex,1)}
    }

    if (sum(ex)<length(compoundsDF$Compound)){
      message("\tThe following compounds are not found in the uncorrected input table,\nPlease check your files:")
      message(compoundsDF$Compound[which(ex==0)])
      message(paste (inFile_," is not processed \n------------------"))
      next
    }else{   #keep relavent columns
      rmv.idx=which(!(1:length(colnames(InputData))%in%ex.i))
      InputData=InputData[,-rmv.idx]
    }

    # Check ordering of Isotopes M+0 M+1 M+3
    allCols=colnames(InputData)
    ixs=c()
    for (namCPD in compoundsDF$Compound){
      ix=which(startsWith(allCols,namCPD))
      curCols=allCols[ix]
      Ms=getMs(curCols,namCPD)
      MsN=as.numeric(unlist(lapply(Ms,function(m){(strsplit(m,"+",fixed=T)[[1]][2])})))
      if (!identical(sort(MsN),MsN)){
        ixs=c(ixs,ix)
        #cat(namCPD)
        #break
      }
    }
    if (length(ixs)>0){
      message("\tThe wrong order of Isotopologues detected!\nPlease check your files:")
      for (n in allCols[ixs]){cat(n,"\n")}
      message(paste0(inFile_," : is not processed \n-----------------------"))
      next
    }

    ##Make selection of compounds or use all compounds
    compounds=c()
    if (tolower(selectComp_)=="l" | tolower(selectComp_)=="list"){
      InputData <- adj.CompDF.initial(InputData,compoundsDF)
      compounds=compoundsDF$Compound[compoundsDF$Included==1]
      message(paste0("\t\tList of included compounds:  ",paste0(compounds,collapse = ", ")))
    }
    if (tolower(selectComp_)=="a" | tolower(selectComp_)=="all"){
      compounds=compoundsDF$Compound
      message("\t\tAll compounds are used.")
    }

    message(paste0("\n",get.prog(prg,LEN)," *  2 --- Reading the associated metadata"))

    metaData<- read.csv(paste0(inPath_,"/",mtdFile),row.names=1,header=T)
    if (!identical(rownames(metaData), rownames(InputData)[2:dim(InputData)[1]])){message(paste0(" * Sample names in data and metadata tables must match and be in the same order. Please check files and try again: ",inFile_));next}

    ## Add metadata to the input dataframe
    InputData=cbind.data.frame(c("Labeling",metaData$Labeling),InputData)
    colnames(InputData)[1]="Labeling"

    ## also select relevent rows for  unlabeled and labaled
    rix<-which(startsWith(tolower(InputData$Labeling),"unlabeled") | startsWith(tolower(InputData$Labeling),"labeled"))
    InputData=InputData[rix,]

    message(paste0("\n",get.prog(prg,LEN)," *  3 --- Substitute missing values with zeros"))
    InputData=sub_Na_Zero(InputData)

    message(paste0("\n",get.prog(prg,LEN)," *  4 --- Imputation of missing values using the minimum for isotope in each group"))
    InputData=imputMin(InputData,compounds,baseLineDet_) ### remove compound bellow baseline in 2/3 of unlabeled M+0 the input by min by isotope

    message(paste0("\n",get.prog(prg,LEN)," *  5 --- Removing all low detection compound(s)"))
    InputData=removeCompNA(InputData)
    compounds<-adj.CompList(InputData,compounds)

    message(paste0("\n",get.prog(prg,LEN)," *  6 --- Correction for isotope natural abundance"))


    ### get indeces of labaled and unlabaled
    #sNames= unique(metaData$Labeling[which(startsWith(metaData$Labeling,"unlabeled") | startsWith(metaData$Labeling,"labeled"))])
    sNames=unique(InputData$Labeling)
    L=list()
    counter=1
    for (i in 1:length(sNames)){
      L[[counter]]<-which(InputData$Labeling==sNames[i])
      counter=counter+1
    }
    names(L)<-sNames

    OutputData=InputData
    OutputData[1:dim(OutputData)[1],2:dim(OutputData)[2]]=0

    lowDet=c()
    for (compound in compounds){

      message(paste0(" Correcting for  ",compound))
      #making the correction matrix based on unlabaled samples
      colIx=which(startsWith(tolower(colnames(InputData)),tolower(compound)))

      if (length(colIx)==0){next}

      #unlabaledV= colMeans(data.matrix(InputData[L[[1]],colIx]))
      unlabaledV= colMeans(data.matrix(InputData[L[["unlabeled"]],colIx]))
      LowDetection=FALSE
      if (sum(unlabaledV)==0 ){      #(length(which(unlabaledV==0))>0 ){
        LowDetection=TRUE; message(paste0("  * Warning * : Low detection for ",compound))
      }else{
        corrMat=make.corrMat(unlabaledV)
        if  (rankMatrix(corrMat)[1]< dim(corrMat)[1]){
          LowDetection=TRUE; message(paste0("  * Warning * : The rank natural abundance correction matrix [CM] is < the no. of its rows or columns. CM cannot be inverted for ",compound,"\n"))
        }else{
          Inv<-solve(corrMat)
        }
      }
      ## make correction of labaled samples
      for (i in 1:length(L)){
        #Take inverse of unlabled sample matrix
        for (ri in L[[i]]){
          if (LowDetection==TRUE){
            OutputData[ri,colIx]<- 0
            lowDet=c(lowDet,compounds)
          }else{
            labaledV=as.numeric(InputData[ri,colIx])

            if (sum(labaledV)==0){
              OutputData[ri,colIx]<- 0
              message(paste0("  * Warning * : Low detection for ",compound));
              lowDet=c(lowDet,compound)
            }else{
              #Matrix multiplication
              outVec<-labaledV%*%Inv
              #Normalize
              outVecN<-(outVec*100)/sum(outVec)
              OutputData[ri,colIx]<- round(outVecN,10)
            }
          }
        }
      }
    }

    message(paste0("\n",get.prog(prg,LEN)," *  7 --- Generating tables for ", inFile_ ," ..."))


    #QC remove the compound with Fractional Enrichment <95 or >105 in unlabaled samples.
    message("\t\tQC: removing the compound(s) with fractional enrichment <95 or >105 in unlabaled samples")
    OutputData<-M0_FE_QC( OutputData,compounds,L)
    OutputData=removeCompNA(OutputData)
    compounds<-adj.CompList(OutputData,compounds) ### remove compounds that were removed in the OutputData table

    ### remove the unlabled samples
    OutputData=OutputData[-L[["unlabeled"]],]

    #table
    write.csv(OutputData,file=paste0(subDir1,"/",outFile))

    message(paste0("\n",get.prog(prg,LEN)," *  8 --- Generating figures for ", inFile_ ," ..."))
    #heatmap
    htmp=data.matrix(OutputData[,2:dim(OutputData)[2]])
    #htmp[which(htmp<0)]<- 0

    ## Summing up M+1, M+2, ...., M+n for compounds and include sums > 5/%
    htmp.TL=matrix(0,nrow=nrow(htmp),ncol=length(compounds),dimnames=list(rownames(htmp),compounds))
    for( i in 1:dim(htmp.TL)[2]){
      comp=colnames(htmp.TL)[i]
      cix=which(startsWith(colnames(htmp),comp))
      subM=htmp[,cix];
      ix=which(colnames(subM)==paste0(comp,".Results")); if (length(ix)==0){stop(" - please check column names of ",comp)}else{subM=subM[,-ix]}
      #subM[which(subM < 0)]<-0
      htmp.TL[,i]=rowSums(subM)
    }

    zscore<-function(v){(v-mean(v,na.rm=T)) / sqrt(var(v,na.rm=T))}#z-scroe: (x-μ)/σ
    cix=which(apply(htmp.TL,2,max)>5)
    ##dynamic height and width
    #htmp.TL.MC=apply(htmp.TL[,cix],2,function(cl){cl/mean(cl)}) #a
    htmp.TL.ZS=apply(htmp.TL[,cix],2,function(cl){zscore(cl)}) #a

    pdf(file=paste0(subDir1,"/",outHtmp),onefile = T,height=10,width=5)
    pheatmap(t(htmp.TL.ZS),main="Total Fractional Enrichment\nsum of M+1, M+2, M+3, ...",border_color = "NA",cluster_rows = FALSE,cluster_cols = FALSE)
    dev.off()

    ### Within replicates groups coefficient of variation
    OutDf= sub_Ng_Zero(OutputData)
    InputCV=coev_df(InputData,compounds); write.csv(InputCV,file=paste0(subDir2,"/", "CV_uncorr.csv"))  # Uncorrected data
    OutCV=coev_df(OutDf,compounds); write.csv(OutCV,file=paste0(subDir2,"/", "CV_corr.FE.csv"))         # Corrected data

    #make cv plot
    if (length(which(is.na(OutCV)))< (dim(OutCV)[1]*dim(OutCV)[2])){
      g1=cvPlt(InputCV,main="Input data - total pool size\nCoefficient of variation")
      #g2=cvPlt(OutCV,main="Franctional Enrichment")
      ggsave(paste0(subDir1,"/CV.pdf"),  g1 , height=11, width=6)
      #ggsave(paste0(subDir1,"/CV.pdf"), ggarrange(g1, g2 , ncol = 2, nrow = 1), height=11, width=9)
    }

    #cvPlt(InputCV,subDir1)

    ###------------------------------------------------------------------------###
    ###  Generating figure for Total /Percentage Enrichmnent table and figure  ###
    ###------------------------------------------------------------------------###
    #sNames=as.character(OutputData$DataFile )
    sNames=as.character(rownames(OutputData))
    OutputData.TE=OutputData
    for (compound in compounds){

      #*** run signal check-----
      #if (runSignal=="cancel"){ runSignal<<-"run";return(message(paste0("\n\t\t*---  MIDA run aborted !   ", format(Sys.time(), "%a, %b %Y %d, %X"),"    ---*       \n\n")))}

      for (i in 1:length(sNames)){
        rIx1=which(tolower(rownames(InputData))==tolower(sNames[i]))
        cIx1=which(startsWith(tolower(colnames(InputData)),tolower(compound)))
        rIx2=which(tolower(rownames(OutputData))==tolower(sNames[i]))
        cIx2=which(startsWith(tolower(colnames(OutputData)),tolower(compound)))
        #Total Enrichment(sum(Labaled Sample enrichment) * Fractional Enrichmet) * 100
        OutputData.TE[rIx2,cIx2]<-(sum(as.numeric(InputData[rIx1,cIx1]))*as.numeric(OutputData[rIx2,cIx2]))/100
      }

    }
    write.csv(OutputData.TE,file=paste0(subDir2,"/", paste0(strsplit(inFile_,"uncorr.csv")[[1]],"_TE.csv")))
    colorz=c("white",colorRamps::primary.colors(),"grey30","grey45")[c(1,15,4,5,12,10,13,6,7,8,9,11,14,16,17,18,19,20,21,22,23,24,25,26,27,28,2,29,3)] #[sample(1:30,30)]#plot(1:30,1:30,pch=19,cex=1,col=colz)
    #plot(1:29,col=c("white",colorRamps::primary.colors(),"grey30")[c(1,16,5,6,13,11,14,  4,7,8,9,10,12,15,17,18,19,20,21,22,23,24,25,26,27,28,2,29,3)],pch=19,cex=2.5)
    names(colorz)= unlist(lapply(1:length(colorz),function(i){paste0("M+",as.character(i-1))}))

    #Total enrichment figures
    gL1<-gL3<-gL4<-c()
    for (compound in compounds){
      #cat(compound,"\n")
      # Total Enrichment matrix for each compound
      cIx=which(startsWith(tolower(colnames(OutputData.TE)),tolower(compound)))
      mt<-data.matrix(OutputData.TE[,cIx])
      mt[which(mt< 0)]=0
      g1=ggBar(mt,colorz,main_="Total Enrichment",ylab_="Total Enrichment")
      gL1[[compound]]<-g1
      #if (nrow(mt)<=12){ H=4;W=4}; if (nrow(mt)>12 &nrow(mt)<15){ H=4;W=8}; if (nrow(mt)>=15 &nrow(mt)<20){ H=4;W=20}; if (nrow(mt)>=20 &nrow(mt)<30){ H=4;W=25};if (nrow(mt)>=30){ H=4;W=30}
      #ggsave(file=paste0(subDir2,"/",compound,"_TE.pdf"),height=4,width=4)

      # Percentage enrichment matrix for each compound
      #mt.prcnt= t(apply(mt,1,function(row){round(100*(row/sum(row)),2)}));
      #g2=ggBar(mt.prcnt,colorz,main_="Percentage Enrichment",ylab_="Percentage Enrichnmet")
      ##if (nrow(mt.prcnt)<=12){ H=4;W=4}; if (nrow(mt.prcnt)>12 &nrow(mt.prcnt)<15){ H=4;W=8}; if (nrow(mt.prcnt)>=15 &nrow(mt.prcnt)<20){ H=4;W=20}; if (nrow(mt.prcnt)>=20 &nrow(mt.prcnt)<30){ H=4;W=28};if (nrow(mt.prcnt)>=30){ H=4;W=30}
      #ggsave(file=paste0(subDir2,"/",compound,"_PRCNT.pdf"),height=4,width=4)

      # --   Averaging biological replicates   -- #

      #  Total Enrichment matrix for each compound (Compined replciates)
      sNames.rep=getRepSetNames(rownames(mt))
      mt.rep=matrix(0,nrow=length(sNames.rep),ncol=ncol(mt),dimnames = list(sNames.rep,colnames(mt)))
      for (i in 1:length(sNames.rep)){
        rix=which(startsWith(rownames(mt),sNames.rep[i]))
        if (length(rix)==1){mt.rep[i,]<-mt[rix,]} else { mt.rep[i,]=colMeans(mt[rix,])}
      }
      #Standard Error calculation
      sde=get.sde(InputData,mt.rep)
      g3=ggBar.sde(mt.rep,sde,colorz,main_="Total Enrichment",ylab_="Total Enrichment")
      gL3[[compound]]<-g3

      #if (nrow(mt.rep)<=12){ H=4;W=4}; if (nrow(mt.rep)>12 &nrow(mt.rep)<15){ H=4;W=8}; if (nrow(mt.rep)>=15 &nrow(mt.rep)<20){ H=4;W=20}; if (nrow(mt.rep)>=20 &nrow(mt.rep)<30){ H=4;W=25};if (nrow(mt.rep)>=30){ H=4;W=30}
      ggsave(file=paste0(subDir2,"/",compound,"_TE.repset.pdf"),height=3,width=3)

      # Percentage enrichment matrix for each compund
      mt.prcnt.rep= t(apply(mt.rep,1,function(row){round(100*(row/sum(row)),2)}))
      g4=ggBar(mt.prcnt.rep,colorz,main_="Percentage Enrichment",ylab_="Percentage Enrichnmet")
      gL4[[compound]]<-g4

      #if (nrow(mt.prcnt.rep)<=12){ H=4;W=4}; if (nrow(mt.prcnt.rep)>12 &nrow(mt.prcnt.rep)<15){ H=4;W=8}; if (nrow(mt.prcnt.rep)>=15 &nrow(mt.prcnt.rep)<20){ H=4;W=20}; if (nrow(mt.prcnt.rep)>=20 &nrow(mt.prcnt.rep)<30){ H=4;W=28};if (nrow(mt.prcnt.rep)>=30){ H=4;W=30}
      #ggsave(file=paste0(subDir2,"/",compound,"_PRCNT.repset.pdf"),height=3,width=3)
    }

    message(" Generating figures for total enrichment (all replicates)")
    grobs = lapply(gL1, "+", theme(plot.margin = unit(c(.8,.9,.8,.9), "cm")))
    grobs_ <- lapply(grobs, ggplotGrob)
    ggsave(paste0(subDir1,"/Total.Enrichment.replicates.pdf"), marrangeGrob(grobs = grobs_ , nrow=3, ncol=2), height=11, width=9)

    message(" Generating figures for total enrichment (combining replicates)")
    grobs = lapply(gL3, "+", theme(plot.margin = unit(c(.8,.9,.8,.9), "cm")))
    grobs_ <- lapply(grobs, ggplotGrob)
    ggsave(paste0(subDir1,"/Total.Enrichment.pdf"), marrangeGrob(grobs = grobs_ , nrow=3, ncol=2), height=11, width=9)

    message(" Generating figures for total percentage enrichment (combining replicates)")
    grobs = lapply(gL4, "+", theme(plot.margin = unit(c(.8,.9,.8,.9), "cm")))
    grobs_ <- lapply(grobs, ggplotGrob)
    ggsave(paste0(subDir1,"/Total.Enrichment_prcnt.pdf"), marrangeGrob(grobs = grobs_ , nrow=3, ncol=2), height=11, width=9)

    message(" Generating TCA cycle figure")
    ##TCA cycle figure
    pdf(file=paste0(subDir1,"/","TCA.pdf"),width=9,height=11)
    tcaPlt(subDir2)
    dev.off()


    message(paste0("\n",get.prog(prg,LEN)," * -------------- ", inFile_ ," processed succesfully!     --------------"))

  }
  message(paste0("\n",get.prog(prg,LEN),"\t\t         MIDA run completed successfully!                "))

  message(paste0("\n\t\t\t\t*---      ", format(Sys.time(), "%a, %b %Y %d, %X"),"      ---*        \n\n"))

}

####Function to run MIDA in the current working directory
midaR <- function(inPath=getwd(), selectComp="All", overWrite=TRUE, baseLineDet=1000) {

  library(midaR)
  library(RColorBrewer)
  library(colorRamps)
  library(reshape)
  library(readxl)
  library(dplyr)
  library(kableExtra)
  library(ggplot2)
  library(gridExtra)
  library(pheatmap)
  library(Matrix)
  library(crayon)
  library(igraph)
  library(magick)
  library(pdftools)

  iArrows <- igraph:::igraph.Arrows

  inFiles=list.files(inPath , pattern= "\\_uncorr.csv$")

  if (overWrite == TRUE) {ow <- "overwrite"
  }else{
    ow <- "skip"
  }

  if (length(inFiles)==0){
    print("No input files found.")

  }else{

      mainFun(inFiles__ = inFiles, inPath_ = inPath, selectComp_ = selectComp,
              baseLineDet_ =  baseLineDet, overWrite_ =  ow)
  }

}


###shiny app
runMIDA <- function(){
  library(midaR)
  library(RColorBrewer)
  library(colorRamps)
  library(reshape)
  library(readxl)
  library(dplyr)
  library(kableExtra)
  library(ggplot2)
  library(gridExtra)
  library(pheatmap)
  library(Matrix)
  library(crayon)
  library(igraph)
  library(magick)
  library(pdftools)

  iArrows <- igraph:::igraph.Arrows

  library(shiny)
  library(shinyWidgets)
  library(shinyalert)
  library(shinyFiles)
  library(fs)


  js <- "
$(document).ready(function(){
  var objDiv = document.getElementById('log');
  // create an observer instance
  var observer = new MutationObserver(function(mutations) {
    objDiv.scrollTop = objDiv.scrollHeight - objDiv.clientHeight;
  });
  // configuration of the observer
  var config = {childList: true};
  // observe objDiv
  observer.observe(objDiv, config);
})
"

  MIDA_ui <- fluidPage(
    #includeCSS('midaR/mycss.css'),
    useShinyalert(),
    shinyjs::useShinyjs(),
    tags$head(tags$script(HTML(js))),

    h5("VAI Metabolism - Mass Isotopologue Distribution Analysis (MIDA)" ,style='padding:5px; color: #0046E7',img(src="midaR/VAI_Logo.jpg" ,align = "right",height=30,width=105)),

    tabsetPanel(

      tabPanel("Console",
               tags$style(HTML("#log{height:700px; font-size:85%; background-color:#4265a5; color:#ffffff}")),
               verbatimTextOutput("log",placeholder=T),

               fluidRow(
                 column(12, align="center", offset = 0,
                        tags$style(HTML("#dirPath{font-size:75%; color:#A801DB}")),
                        verbatimTextOutput("dirPath", placeholder = TRUE)
                 ),

                 column(12,align="center",# offset = 3,
                        shinyDirButton("dir", "Choose folder", "Upload",style='padding:5px; font-size:80%; color: #0046E7;width: 90px'),
                        actionButton("btn_MIDA", "Run MIDA",style='padding:5px; font-size:80%; color: #0046E7; width: 90px'),
                        actionButton("btn_Cancel", "Cancel",style='padding:5px; font-size:80%; color: #0046E7; width: 90px')
                 )
               ),
               icon = icon("laptop")
      ),

      tabPanel("Plots",

               uiOutput("pdfview"),

               column(12,align="center",
                      shinyFilesButton("filePDF", "Select PDF", "Please select a file", multiple = FALSE,style='padding:5px; font-size:80%; color: #0046E7;width: 90px')
               ),
               icon = icon("chart-bar")
      ),

      tabPanel("Tables",

               tableOutput("kableview"),

               column(12,align="center",
                      shinyFilesButton("fileKable", "Select CSV", "Please select a file", multiple = FALSE,style='padding:5px; font-size:80%; color: #0046E7;width: 90px')
               ),
               icon = icon("fal fa-table")

      ),

      tabPanel("Settings",

               sidebarPanel(
                 tags$style(".well {background-color:#4265a5; color:#ffffff;}"),
                 br(),
                 h5("Descrition of all the settings used in this program:"),
                 br(),

                 h5("Natural abundance correction method"),
                 p("This allows to specify whether labelled samples or theoretical mass isotopologue distrubution values should be used"),
                 br(),

                 h5("Baseline detection level"),
                 p("The threshold detection level at which compounds are filtered. By default the method excludes compounds
                  if  2/3 of unlabeled samples M+0 are below 1000 or if all labeled sample entries are < 1000"),
                 br(),

                 h5("Specifying the list of compounds"),
                 p("A list of compounds are generated within the {samplename_cmpd.csv} file in which the
                  user can mark compounds of interest for which to correct for natural abundance"),
                 br(),br(),

                 h5("Running mode"),
                 p("Determines whether to skip or overwrite the existing output directories of processed input files.")

               ),

               mainPanel(

                 br(),br(),
                 #tags$style(HTML("#midaCMPDList {font-size:80%; color:#0046E7")),
                 radioButtons("midaAnalisis", "Natural abundance correction method:",
                              c("Using labeled sample" = "Labeled samples","Theoretical MIDs" = "Theoretical MIDs")
                 ),

                 br(),
                 #tags$style("#midaBaseLDet {font-size:80%;width:120px;color:#0046E7;}"),
                 numericInput("midaBaseLDet", label="Baseline detection:", value = 1000,  min = 500, max = 3000, step = 100),

                 br(),
                 #tags$style(HTML("#midaCMPDList {font-size:80%; color:#0046E7")),
                 radioButtons("midaCMPDList", "Specify the compounds to analyse:",
                              c("All compound" = "All","Selected List" = "List")
                 ),
                 br(),
                 #tags$style(HTML("#midaOverWrite {font-size:80%; color:#0046E7")),
                 radioButtons("midaOverWrite", "Run mode:",
                              c("Skip processed data" = "Skip","Overwrite processed data"="Overwrite")
                 )

               ),

               icon = icon("sliders-h")
      ),

      tabPanel("Documentation",
               #tags$iframe(style="height:700px;  width:100%; scrolling=yes",src="MIDA.documentation.pdf"),

               uiOutput("docview"),

               column(12,align="center",
                      actionButton("btn_Documentation", "Change",style='padding:5px; font-size:80%; color: #0046E7; width: 90px'),
               ),

               icon = icon("book")
      )
    )
  )


  MIDA_server <- function(input, output,session) {
    #libPath = "C:/Program Files/R/R-3.6.2/library"
    library(shinyjs)
    #source("utils.R")

    shiny::addResourcePath("midaR", system.file("www", package="midaR"))

    #Initial setting
    if (file.exists("midaR/Plot.pdf")){file.remove("midaR/Plot.pdf")}
    if (file.exists("midaR/Kable.csv")){file.remove("midaR/Kable.csv")}

    output$docview <- renderUI({tags$iframe(style="height:700px; width:100%; scrolling=yes", src="MIDA.README.pdf")})

    volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())


    shinyDirChoose(input, "dir", roots = volumes, session = session, restrictions = system.file(package = "base"),filetypes = c('', "csv","tsv", "xls","xlsx") )
    shinyFileChoose(input, "filePDF", roots = volumes, session = session,filetypes = c('', "pdf"))
    shinyFileChoose(input, "fileKable", roots = volumes, session = session,filetypes = c('', "csv","tsv", "xls","xlsx") )


    #Get the input directory path
    output$dirPath <- renderPrint({
      if (is.integer(input$dir)) {
        cat("Please choose the input directory.")
      } else {

        #inPath<<-parseDirPath(volumes, input$dir)
        parseDirPath(volumes, input$dir)
      }
    })

    # Action Button for running MIDA
    observeEvent(input$btn_MIDA, {

      inFiles=list.files(parseDirPath(volumes, input$dir) , pattern= "\\_uncorr.csv$")

      if (length(inFiles)==0){

        shinyalert("Warning", "No directory containing valid input files was selected!", type = "warning")
      }else{
        # Initial parmeters
        selectComp=input$midaCMPDList
        overWrite =input$midaOverWrite
        baseLineDet=input$midaBaseLDet
        inPath=parseDirPath(volumes, input$dir)


        shinyjs::disable("filePDF")                    #disable some controls when the program is running
        shinyjs::enable("fileKable")
        shinyjs::disable("btn_Documentation")

        withCallingHandlers(

          mainFun(inFiles, inPath, selectComp, baseLineDet, overWrite),

          #can use "warning" instead/on top of "message" to catch warnings too
          message = function(m) {
            shinyjs::html("log", m$message, add=TRUE)
          }
        )
        shinyjs::enable("filePDF")                    #dAnable some controls when the program has completed running
        shinyjs::enable("fileKable")
        shinyjs::enable("btn_Documentation")
      }
    })


    #Plot visualization
    observe({
      req(input$filePDF)
      output$pdfview <- renderUI({
        #tags$iframe(style="height:700px; width:100%; scrolling=yes", src="")
        inPathPlot=parseFilePaths(volumes, input$filePDF)

        if(file.exists("midaR/Plot.pdf")) {file.remove("midaR/Plot.pdf")}

        file.copy(inPathPlot[[length(inPathPlot)]],"midaR/Plot.pdf", overwrite = T)

        if(file.exists("midaR/Plot.pdf")){
          tags$iframe(style="height:700px; width:100%; scrolling=yes", src="Plot.pdf")
        }

      })

    })


    #Plot visualization
    observe({
      req(input$fileKable)
      output$kableview <- function() {

        inPathTable=parseFilePaths(volumes, input$fileKable)

        if(file.exists("midaR/Kable.csv")) {file.remove("midaR/Kable.csv")}
        file.copy(inPathTable[[length(inPathTable)]],"midaR/Kable.csv", overwrite = T)
        Kable.nm=basename(as.character(inPathTable[[length(inPathTable)]]))

        if(file.exists("midaR/Kable.csv")){
          ktable<-read.csv(file="midaR/Kable.csv",stringsAsFactors = F)
          knitr::kable( ktable  ,caption=Kable.nm) %>%
            kable_styling(bootstrap_options =c("striped", "hover","condensed"),font_size = 12) %>% scroll_box(width = "100%px", height = "450px")
        }
      }

    })


    output$mtcars_kable <- function() {
      req(input$mpg)
      # mtcars %>%
      #         mutate(car = rownames(.)) %>%
      #         select(car, everything()) %>%
      #         filter(mpg <= input$mpg) %>%
      #         knitr::kable("html") %>%
      #         kable_styling("striped", full_width = F) %>%
      #         scroll_box(width = "300px", height = "300px"))%>%
      #         add_header_above(c(" ", "Group 1" = 5, "Group 2" = 6))

    }

    ### Documentation visualization
    cntr<<-1
    observeEvent(input$btn_Documentation, {
      if (cntr%%2==0){
        output$docview <- renderUI({tags$iframe(style="height:700px; width:100%; scrolling=yes", src="MIDA.README.pdf")})
      }else{
        output$docview <- renderUI({tags$iframe(style="height:700px; width:100%; scrolling=yes", src="MIDA.documentation.pdf")})
      }
      cntr<<-cntr+1
    })


    #Cancel button to interrupt the analysis
    observeEvent(input$btn_Cancel, {
      shinyalert("This functionality still under implementation.")
    }
    )


    # Automatically stop a Shiny app when closing the browser tab
    session$onSessionEnded(stopApp)

  }



  library(shiny)
  library(shinyjs)



  # Global setting and variables
  options(stringsAsFactors = F, scipen=999, browser="C:/Program Files (x86)/Google/Chrome/Application/chrome.exe")
  shinyApp(ui=MIDA_ui, server=MIDA_server)
  #shinyApp( ui=functions(request) MIDA_ui, server=MIDA_server, ...)
}

#runMIDA()


































































