path<-"/Users/eric/Google Drive/!DIPnet_DB/Repository"
setwd(path)
filecount<-0
for(file in (list.files(paste(path,"/1-cleaned_QC2_mdfasta_files", sep=""),pattern="mdfasta",full.names=F))){ #file<-"name_checked_mdfasta_trimax_CO1_td.txt"
  filecount<-filecount+1
  print(file)
  
  sp<-read.table(file=paste("1-cleaned_QC2_mdfasta_files", file, sep="/"),header=T,stringsAsFactors=F,sep="\t",na.strings=c("","NA","#N/A"),strip.white=T,fill=T,comment.char="",quote="", colClasses=c("materialSampleID"="character"))
  
  
  #Remove ambiguity codes
  sp$sequence<-gsub(x=sp$sequence,pattern = "[MRWSYKVHDB?]",replacement = "N")

  #Change terminal gaps to N
  for(g in 50:1){
    gaps<-rep.int(x="-",g)
    gaps<-paste(gaps,collapse="")
    ns<-rep.int(x="N",g)
    ns<-paste(ns,collapse="")
    sp$sequence<-gsub(x=sp$sequence,pattern = paste("^",gaps,sep=""),replacement = ns)
    sp$sequence<-gsub(x=sp$sequence,pattern = paste(gaps,"$",sep=""),replacement = ns)
  }
 
  #FORMAT CONVERSION
  cat("converting data to fasta", "\n")
  #convert to seqinr alignment (need to specify that we are using the as.alignment() function from seqinr rather than ape)
  spseqs<-seqinr::as.alignment(nb=length(sp[,1]),nam=sp$materialSampleID, seq=sp$sequence)
  
  #write to fasta format for input to FaBox
  seqinr::write.fasta(as.list(spseqs$seq),spseqs$nam,file.out=paste("./1_cleaned_nonaligned_fasta_files/",gsub(x=file,pattern=".txt",replacement=""),".fasta",sep=""))

}
print(paste(filecount,"files processed"))