
#Required Libraries
require(xlsx) #read/write xlsx
require(plyr) #"count" function

#Function to loop through properly formatted ISAData spreadsheet
#Only to use for ISA Data generated from a single DNA Source
#Input: Path to ISA Data File
#Output: Path to New Output File (or existing file--default is to append)
ISAData.Condense.SR.Loop<- function(inPath, outPath)
{
  numSheets<-length(getSheets(loadWorkbook(inPath)))
  for(index in 1:numSheets) #Loops function isadata.condense.singlerun for each sheet
    #in input file parameter
  {
    tmp.df<-read.xlsx(inPath,index)
    tmp.df.Counts<-isadata.condense.singlerun(tmp.df)
    SheetName<-getSheets(loadWorkbook(inPath))[[index]]$getSheetName()
    write.xlsx(tmp.df.Counts,outPath,sheetName=paste(SheetName,nrow(tmp.df),"sites",sep="_"),col.names = TRUE, row.names = FALSE,append=TRUE)
  }
}

#Function to loop through properly formatted ISAData spreadsheet
#Only to use for ISA Data generated from multiple DNA Sources
#Input: Path to ISA Data File
#Output: Path to New Output File (or existing file--default is to append)
ISAData.Condense.MR.Loop<- function(inPath, outPath)
{
  numSheets<-length(getSheets(loadWorkbook(inPath)))
  for(index in 1:numSheets) #Loops multirun condense function
  {
    tmp.df<-read.xlsx(inPath,index)
    tmp.df.Counts<-isadata.condense.multirun(tmp.df)
    SheetName<-getSheets(loadWorkbook(inPath))[[index]]$getSheetName()
    write.xlsx(tmp.df.Counts,outPath,sheetName=paste(SheetName,nrow(tmp.df),"sites",sep="_"),col.names = TRUE, row.names = FALSE,append=TRUE)
  }
}

#Condensing function to generate frequency table 
#Of detected IS's from data generated from a single
#DNA source.
#input: Dataframe with column names:
#"chr", "LTR","bp_pseudo3LTR","provirus_orientation" ("GeneName" optional)
#Output: data.frame() with frequency of detected IS's
#If only 3LTR or 5LTR is detected for a given IS, that count is true
#If both are detected, the count is max(count(3LTR),count(5LTR))
isadata.condense.singlerun<- function(ISA.df)
{
  ISA.df<-data.frame(lapply(ISA.df, as.character), stringsAsFactors = FALSE) #Make data of type 'chr'
  #Create "raw" frequency table using identical rows in supplied data.frame() with column names given by colnames(data.frame())
  #Note, the same IS at the 3LTR and 5LTR will be counted separately
  tmp.counts.raw<-count(ISA.df, colnames(ISA.df))  
  tmp.counts.raw<-tmp.counts.raw[order(-tmp.counts.raw$freq),] #Order by descending number of counts
  #create "binary" data.frame() by removing the "LTR" (and "freq") column(s) from the raw frequency table
  tmp.LTRbinary<-tmp.counts.raw[,!names(tmp.counts.raw) %in% c("LTR","freq")]
  tmp.LTRbinary.counts<-count(tmp.LTRbinary,colnames(tmp.LTRbinary)) #determine which IS's are detected at both ends
  tmp.has2LTRdetect<-tmp.LTRbinary.counts[which(tmp.LTRbinary.counts$freq==2),] #Move those to new data.frame()
  
  #For each IS detected at the 5 and 3 LTR, find counterparts in raw count dataframe
  #remove the row in raw count with the smaller number of counts
  for(row in 1:nrow(tmp.has2LTRdetect))
  {
    tmp.row<-tmp.has2LTRdetect[row,]
    tmp.counts.raw.2LTR<-tmp.counts.raw[which(tmp.counts.raw$chr == tmp.row$chr &
                                                tmp.counts.raw$bp_pseudo3LTR == tmp.row$bp_pseudo3LTR &
                                                tmp.counts.raw$provirus_orientation == tmp.row$provirus_orientation),]
    row.delete<-row.names(tmp.counts.raw.2LTR[which.min(tmp.counts.raw.2LTR$freq),])[1]
    tmp.counts.raw<-tmp.counts.raw[!rownames(tmp.counts.raw) %in% row.delete, ]
  }
  
  #return frequency table without LTR column, corrected for LTR detection
  return (tmp.counts.raw[,!names(tmp.counts.raw) %in% "LTR"])
}

#Condensing function to generate frequency table 
#Of detected IS's from data generated from Multiple DNA Sources.
#input: Dataframe with column names:
#"Sample_run", "chr", "LTR","bp_pseudo3LTR","provirus_orientation" ("GeneName" optional)
#Output: data.frame() with frequency of detected IS's
#data.frame() will have multiple "freq" columns, dependent on
#the number of runs in the original data.
#Dependent Function: isadata.condense.singlerun(ISA.df)
#See comments on this function for more info
isadata.condense.multirun<-function(ISA.df)
{
  ISA.df<-data.frame(lapply(ISA.df, as.character), stringsAsFactors = FALSE) #Make data of type 'chr'
  ISA.df$Sample_run<-factor(ISA.df$Sample_run) #Factor Sample_run column
  ISA.List.raw<-split(ISA.df, ISA.df$Sample_run) #Split on the above factor into List of dataframes
  ISA.List.raw<-lapply(ISA.List.raw, function(x) x[!(names(x) %in% c("Sample_run"))]) #Remove "Sample_run" column from all dataframes
  ISA.List.SRCond<-list() #Create empty list
  
  #For each data frame in the raw list
  #Put condensed single-run data.frame into created empty list
  #Rename "freq" column to "freq_'sample_run'" where 'sample_run' is the name of the data.frame
  #in the raw list
  for(name in names(ISA.List.raw))
  {
    ISA.List.SRCond[[name]]<-isadata.condense.singlerun(ISA.List.raw[[name]])
    colnames(ISA.List.SRCond[[name]])[ncol(ISA.List.SRCond[[name]])]<-
      paste(colnames(ISA.List.SRCond[[name]])[ncol(ISA.List.SRCond[[name]])], name,sep="_")
  }
  
  #Create master column name vector
  colname.vec<-colnames(ISA.List.raw[[1]][,!names(ISA.List.raw[[1]]) %in% "LTR"])
  #Add unique "freq" columns to column name vector
  for(i in 1:length(ISA.List.SRCond))
  {
    colname.vec<-c(colname.vec, colnames(ISA.List.SRCond[[i]])[length(ISA.List.SRCond[[i]])])
  }
  #Add "TotalFreq" column
  colname.vec<-c(colname.vec,"TotalFreq")
  #Create Final dataframe to return
  ISA.Master.df<-data.frame(matrix(ncol=length(colname.vec)))
  ISA.Master.df<-setNames(ISA.Master.df,colname.vec)

  #Make Indicator list for number of runs IS is detected
  tmp.list<-list()
  for(name in names(ISA.List.SRCond))
  {
    tmp.list[[name]]<-ISA.List.SRCond[[name]][,-length(ISA.List.SRCond[[name]])]
  }
  tmp.df<-do.call("rbind",tmp.list)
  tmp.df.counts<-count(tmp.df,colnames(tmp.df))
  tmp.df.counts<-tmp.df.counts[order(-tmp.df.counts$freq),]
  ###
  
  #Collate all into final dataframe
  for(row in 1:nrow(tmp.df.counts))
  {
    tmp.row<-tmp.df.counts[row,]
    numDetect<-0
    for(name in names(ISA.List.SRCond))
    {
      if(nrow(ISA.List.SRCond[[name]][which(ISA.List.SRCond[[name]]$chr == tmp.row$chr &
                                            ISA.List.SRCond[[name]]$bp_pseudo3LTR == tmp.row$bp_pseudo3LTR &
                                            ISA.List.SRCond[[name]]$provirus_orientation == tmp.row$provirus_orientation),])
         ==0)
      {
        ISA.Master.df[row,paste("freq",name,sep="_")]<-0
      } else {
        ISA.Master.df[row,paste("freq",name,sep="_")]<-ISA.List.SRCond[[name]][which(ISA.List.SRCond[[name]]$chr == tmp.row$chr &
                                                                                       ISA.List.SRCond[[name]]$bp_pseudo3LTR == tmp.row$bp_pseudo3LTR &
                                                                                       ISA.List.SRCond[[name]]$provirus_orientation == tmp.row$provirus_orientation),
                                                                               paste("freq",name,sep="_")]
      }
    }
    ISA.Master.df[row,"chr"]<-tmp.row$chr
    ISA.Master.df[row,"provirus_orientation"]<-tmp.row$provirus_orientation
    ISA.Master.df[row,"bp_pseudo3LTR"]<-tmp.row$bp_pseudo3LTR
    if("GeneName" %in% colnames(tmp.row))
    {
      ISA.Master.df[row,"GeneName"]<-tmp.row$GeneName
    }
  }
  #Sum "freq" columns for "TotalFreq" column
  ISA.Master.df$TotalFreq<-apply(ISA.Master.df[,(length(colnames(ISA.List.raw[[1]][,!names(ISA.List.raw[[1]]) %in% "LTR"]))+1):(length(ISA.Master.df)-1)],1,sum)
  

  #Return dataframe ordered on descending "TotalFreq" 
  return(ISA.Master.df[order(-ISA.Master.df$TotalFreq),])
}