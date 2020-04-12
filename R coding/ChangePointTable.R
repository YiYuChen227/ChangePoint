library(changepoint) 
CPA.Flow <- function(input,data,YName,Time,Para.change.number,multiple,AfterSize_threshold,Pvalue_threshold){
  # data preprocess 
  data[which(is.na(data),arr.ind = T)] <- ""
  data <- Filter.Empty(data = data,ColumnName = YName)$data
  data <- Filter.Empty(data = data,ColumnName = Time)$data
  data <- Order.Y(data = data,Time = Time)$data
  data[,Time] <- format(as.POSIXct(as.character(data[,Time])),format = "%Y/%m/%d %H:%M:%S")
  # Take Y value
  var_Y <- data[,YName]
  var_Y <- as.numeric(as.character(var_Y))
  
  # CPA
  m.time <- CPA.model(var_Y = var_Y,Para.change.number = Para.change.number)$m.time
  output_table <- Region.table(data = data,row = m.time@cpts,YName = YName,Time = Time)$output_table

  output_table <- CPA.Filter(data = data ,output_table = output_table,YName = YName,Time = Time,
                             multiple = multiple,AfterSize_threshold = AfterSize_threshold,Pvalue_threshold = Pvalue_threshold)$output_table
  
  output_table <- CPA.Rank(output_table = output_table)$output_table
  
  #plot(1:length(var_Y),var_Y)
  #abline(v = c(as.numeric(as.character(output_table[,"Row"]))),col = 2)
  
  return(list(output_table = output_table))
}
Filter.Empty <- function(data,ColumnName){
  if(length(which(data[,ColumnName]==""))!=0L){
    data <- data[-which(data[,ColumnName]==""),]
  }else{
    data <- data
  }
  return(list(data = data))
}
Order.Y <- function(data,Time){
  if(length(which(data[,Time]==""))==0L){
    data[,Time] <- as.POSIXct(as.character(data[,Time]))
    data <- data[order(data[,Time],decreasing = F),]
  }else{
    data <- NULL
  }
  return(list(data = data))
}
CPA.model <- function(var_Y,Para.change.number){
  if(length(var_Y) < 7L){
    m.time <- cpt.meanvar(data = rep(0,10))
  }else if(length(var_Y) >= 7L & length(var_Y) <= Para.change.number*2){
    m.time <- cpt.meanvar(data = var_Y,method = "BinSeg",penalty = "Asymptotic",pen.value = 0.01,class = T,Q = trunc(Para.change.number/2) - 1,param.estimates = T)
  }else{
    m.time <- cpt.meanvar(data = var_Y,method = "BinSeg",penalty = "Asymptotic",pen.value = 0.01,class = T,Q = Para.change.number,param.estimates = T)
  }
  return(list(m.time = m.time))
}
Generate.CPA.data <- function(data,X.loc,count,YName){
  # count is change point time (line) for one line
  # X.loc contained Start_time(0) , End_time(final) and  all change point time
  XGroup <- rep(NA,nrow(data))
  XGroup[seq(X.loc[count] + 1,X.loc[count + 1],by = 1)] <- "TimeBefore"
  XGroup[seq(X.loc[count + 1] + 1 ,X.loc[count + 2],by = 1)] <- "TimeAfter"
  output_data <- cbind(data,XGroup)
  output_data <- output_data[complete.cases(output_data),]
  colnames(output_data)[grep(YName,colnames(output_data))] <- "YPARAM"
  output_data[,"YPARAM"] <- as.numeric(as.character(output_data[,"YPARAM"]))
  data_TimeBefore <- output_data[output_data[,"XGroup"]=="TimeBefore",] ; 
  data_TimeAfter <- output_data[output_data[,"XGroup"]=="TimeAfter",]
  return(list(data_TimeBefore = data_TimeBefore,data_TimeAfter = data_TimeAfter))
}
Two.Mean.Test <- function(group1,group2){
  Var.Test <- var.test(group1,group2,alternative = "two.side")
  if(is.na(Var.Test$p.value)){
    F_Stat <- F_PValue <- T_Stat <- T_PValue <- NaN
  }else{
    F_Stat <- as.numeric(Var.Test$statistic) ; F_PValue <- as.numeric(Var.Test$p.value)
    if(Var.Test$p.value <= 0.05){
      Ttest <- t.test(group1,group2,alternative = "two.side",var.equal = F)
      T_Stat <- as.numeric(Ttest$statistic) ; T_PValue <- as.numeric(Ttest$p.value)
    }else{
      Ttest <- t.test(group1,group2,alternative = "two.side",var.equal = T) ;
      T_Stat <- as.numeric(Ttest$statistic) ; T_PValue <- as.numeric(Ttest$p.value)
    }
  }
  return(list(F_Stat = F_Stat,F_PValue = F_PValue ,T_Stat = T_Stat , T_PValue = T_PValue))
}
Region.table <- function(data,row,YName,Time){
  X.loc <- c(0,row) #row included the last time point.
  T_PValue <- F_PValue <- T_Stat <- F_Stat <- Mean <- Median <- Sigma <- Size <- NULL
  if(length(X.loc) > 2){ # Make sure that existing the Cut_Point by CPA model
    for(count in seq(1,length(X.loc)-2,by = 1)){
      # get the data set of two region
      output <- Generate.CPA.data(data = data,YName = YName,X.loc = X.loc,count = count)
      data_TimeBefore <- output$data_TimeBefore ; data_TimeAfter <- output$data_TimeAfter
      # Get hypothesis result for two population mean T test
      Meantest <- Two.Mean.Test(group1 = data_TimeBefore[,"YPARAM"],group2 = data_TimeAfter[,"YPARAM"])
      F_Stat <- c(F_Stat, as.numeric(Meantest$F_Stat)) ; T_Stat <- c(T_Stat ,as.numeric(Meantest$T_Stat)) 
      F_PValue <- c(F_PValue, Meantest$F_PValue) ; T_PValue <- c(T_PValue ,Meantest$T_PValue)
      # get statistic 
      Mean <- c(Mean,mean(data_TimeBefore[,"YPARAM"])) ; Median <- c(Median,median(data_TimeBefore[,"YPARAM"])) 
      Sigma <- c(Sigma,sd(data_TimeBefore[,"YPARAM"])) ; Size <- c(Size,length(data_TimeBefore[,"YPARAM"])) 
      if(count == length(seq(1,length(X.loc)-2,by = 1))){
        Mean <- c(Mean,mean(data_TimeAfter[,"YPARAM"])) ; Median <- c(Median,median(data_TimeAfter[,"YPARAM"])) 
        Sigma <- c(Sigma,sd(data_TimeAfter[,"YPARAM"])) ; Size <- c(Size,length(data_TimeAfter[,"YPARAM"])) 
        F_Stat <- c(F_Stat,NaN) ; F_PValue <- c(F_PValue,NaN);T_Stat <- c(T_Stat ,NaN);T_PValue <- c(T_PValue ,NaN)
      }
    }
    
    if(length(row)==2){
      Start <- as.character(data[1,"TIME_ID"])
    }else{
      Start <- as.character(data[c(1,row[seq(1,length(row)-2,by = 1)]),"TIME_ID"])
    }
    CutPoint <- as.character(data[row[seq(1,length(row)-1,by = 1)],"TIME_ID"])
    
    End <- as.character(data[c(row[seq(2,length(row),by = 1)]),"TIME_ID"])
    output_table <- as.data.frame(cbind(seq(1,length(row)-1,by = 1), 
                                        Start,#Start_time
                                        CutPoint,#as.character(data[row[seq(1,length(row)-1,by = 1)],grep(Time,colnames(data))]), #Cut_Point time
                                        End,#as.character(data[c(row[seq(2,length(row),by = 1)]-1),grep(Time,colnames(data))]), #End_time
                                        Mean[seq(1,length(row)-1,by = 1)],Mean[seq(2,length(row),by = 1)],
                                        Median[seq(1,length(row)-1,by = 1)],Median[seq(2,length(row),by = 1)],
                                        Sigma[seq(1,length(row)-1,by = 1)],Sigma[seq(2,length(row),by = 1)],
                                        Size[seq(1,length(row)-1,by = 1)],Size[seq(2,length(row),by = 1)],
                                        T_PValue[seq(1,length(row)-1,by = 1)],
                                        row[seq(1,length(row)-1,by = 1)]))
    colnames(output_table) <- c("No","Start_time","Cut_Point","End_time","Before_Mean","After_Mean","Before_Median","After_Median",
                                "Before_Sigma","After_Sigma","Before_Size","After_Size","T_PValue","Row")
    
    output_table <- CPA.type.adjust(output_table = output_table)$output_table
    
    output_table[,"Mean_Diff"] <- (output_table[,"After_Mean"] - output_table[,"Before_Mean"])
    output_table[,"Before"] <- output_table[,"Mean_Diff"]/ifelse(output_table[,"Before_Sigma"]==0,1e-10,output_table[,"Before_Sigma"])
    output_table[,"After"] <- output_table[,"Mean_Diff"]/ifelse(output_table[,"After_Sigma"]==0,1e-10,output_table[,"After_Sigma"])
    output_table[,"Difference"] <- ifelse(abs(output_table[,"After"]) > abs(output_table[,"Before"]),
                                          sqrt(output_table[,"After_Size"])*output_table[,"After"] ,sqrt(output_table[,"Before_Size"])*output_table[,"Before"])
    output_table <- output_table[,c("No","Start_time","Cut_Point","End_time","Before_Mean","After_Mean","Mean_Diff",
                                    "Before_Median","After_Median","Before_Sigma","After_Sigma","Before_Size","After_Size",
                                    "Before","After","T_PValue","Difference","Row")]
  }else{
    ColName <- c("No","Start_time","Cut_Point","End_time","Before_Mean","After_Mean","Mean_Diff",
                 "Before_Median","After_Median","Before_Sigma","After_Sigma","Before_Size","After_Size","Before","After","T_PValue","Difference","Row")
    output_table <- data.frame(matrix(ncol = length(ColName),nrow = 0))
    colnames(output_table) <- ColName
  }
  
  return(list(output_table = output_table))
}
CPA.Filter <- function(data,output_table,YName,Time,multiple,AfterSize_threshold,Pvalue_threshold){
  output_table <- CPA.type.adjust(output_table = output_table)$output_table
  count <- 1
  while(count==1){
    output_table <- Mean.diff.Iterate(output_table = output_table,data = data,YName = YName,Time = Time,multiple = multiple)$output_table
    output_table <- After.Size.Iterate(output_table = output_table,AfterSize_threshold = AfterSize_threshold,data = data,YName = YName,Time = Time)$output_table
    output_table <- P.Value.Iterate(output_table = output_table,Pvalue_threshold = Pvalue_threshold,data = data,YName = YName,Time= Time)$output_table
    if(all(abs(output_table[,"After_Size"]) > AfterSize_threshold) & all(abs(output_table[,"Before_Size"]) > AfterSize_threshold) 
       &  all(abs(output_table[,"T_PValue"]) <= Pvalue_threshold)) count <- 0
    if(nrow(output_table)==0L) count <- 0
  }
  if(nrow(output_table)!=0L){
    output_table <- output_table[,c("No","Start_time","Cut_Point","End_time","Before_Mean","After_Mean","Mean_Diff",
                                    "Before_Median","After_Median","Before_Sigma","After_Sigma","Before_Size","After_Size",
                                    "Before","After","T_PValue","Difference","Row")]
  }
  return(list(output_table = output_table))
}
Mean.diff.Iterate <- function(output_table,data,YName,Time,multiple){
  count <- 1
  if(nrow(output_table)!=0L){
    while(count==1){
      nr <- nrow(output_table) ; 
      Before_Sigma <- multiple*output_table[,"Before_Sigma"] ; 
      After_Sigma <- multiple*output_table[,"After_Sigma"]
      if(any(abs(output_table[,"Difference"]) <= Before_Sigma) | any(abs(output_table[,"Difference"]) <= After_Sigma)){
        row <- abs(output_table[,"Difference"]) > Before_Sigma | abs(output_table[,"Difference"]) > After_Sigma
        output_table <- output_table[which(row),]
        output_table <- Region.table(data = data,row = c(output_table[,"Row"],nrow(data)),YName = YName,Time = Time)$output_table
        if(nrow(output_table)==0L | nrow(output_table)==nr) count <- 0
      }else{
        count <- 0
      }
    }
  }else{
    ColName <- c("No","Start_time","Cut_Point","End_time","Before_Mean","After_Mean","Mean_Diff",
                 "Before_Median","After_Median","Before_Sigma","After_Sigma","Before_Size","After_Size","Before","After","T_PValue","Difference","Row")
    output_table <- data.frame(matrix(ncol = length(ColName),nrow = 0))
    colnames(output_table) <- ColName
  }
  return(list(output_table = output_table))
}
After.Size.Iterate <- function(output_table,AfterSize_threshold,data,YName,Time){
  count <- 1
  if(nrow(output_table)!=0L){
    while(count==1){
      if(any(abs(output_table[,"After_Size"]) <= AfterSize_threshold) | abs(output_table[1,"Before_Size"]) <= AfterSize_threshold){
        #if(nrow(output_table)!=0L) if(abs(as.numeric(as.character(output_table[1,"Before_Size"]))) <= AfterSize_threshold) output_table <- output_table[-1,]
        output_table <- output_table[abs(output_table[,"After_Size"]) > AfterSize_threshold,]
        if(nrow(output_table)!=0L) if(abs(output_table[1,"Before_Size"]) <= AfterSize_threshold) output_table <- output_table[-1,]
        output_table <- Region.table(data = data,row = c(as.numeric(as.character(output_table[,"Row"])),nrow(data)),YName = YName,Time = Time)$output_table
        if(nrow(output_table)==0L) count <- 0
      }else{
        count <- 0
      }
    }
  }else{
    ColName <- c("No","Start_time","Cut_Point","End_time","Before_Mean","After_Mean","Mean_Diff","Before_Median","After_Median","Before_Sigma","After_Sigma","Before_Size","After_Size","Before","After","T_PValue","Difference","Row")
    output_table <- data.frame(matrix(ncol = length(ColName),nrow = 0))
    colnames(output_table) <- ColName
  }
  return(list(output_table = output_table))
}
P.Value.Iterate <- function(output_table,Pvalue_threshold,data,YName,Time){
  count <- 1
  if(any(is.na(output_table[,"T_PValue"]))) output_table <- output_table[-which(is.na(output_table[,"T_PValue"])),]
  if(nrow(output_table)!=0L){
    while(count==1){
      if(any(abs(output_table[,"T_PValue"]) > Pvalue_threshold)){
        output_table <- output_table[abs(output_table[,"T_PValue"]) <= Pvalue_threshold,]
        output_table <- Region.table(data = data,row = c(output_table[,"Row"],nrow(data)),YName = YName,Time = Time)$output_table
        if(nrow(output_table)==0L) count <- 0
      }else{
        output_table <- Region.table(data = data,row = c(output_table[,"Row"],nrow(data)),YName = YName,Time = Time)$output_table
        count <- 0
      }
    }
  }else{
    ColName <- c("No","Start_time","Cut_Point","End_time","Before_Mean","After_Mean","Mean_Diff","Before_Median","After_Median","Before_Sigma","After_Sigma","Before_Size","After_Size","Before","After","T_PValue","Difference","Row")
    output_table <- data.frame(matrix(ncol = length(ColName),nrow = 0))
    colnames(output_table) <- ColName
  }
  return(list(output_table = output_table))
}
CPA.type.adjust <- function(output_table){
  if(nrow(output_table)!=0L){
    output_table[,"No"] <- as.numeric(as.character(output_table[,"No"]))
    output_table[,"Start_time"] <- as.character(output_table[,"Start_time"])
    output_table[,"Cut_Point"] <- as.character(output_table[,"Cut_Point"])
    output_table[,"End_time"] <- as.character(output_table[,"End_time"])
    output_table[,"Before_Mean"] <- as.numeric(as.character(output_table[,"Before_Mean"]))
    output_table[,"After_Mean"] <- as.numeric(as.character(output_table[,"After_Mean"]))
    output_table[,"Before_Median"] <- as.numeric(as.character(output_table[,"Before_Median"]))
    output_table[,"After_Median"] <- as.numeric(as.character(output_table[,"After_Median"]))
    output_table[,"Before_Sigma"] <- as.numeric(as.character(output_table[,"Before_Sigma"]))
    output_table[,"After_Sigma"] <- as.numeric(as.character(output_table[,"After_Sigma"]))
    output_table[,"Before_Size"] <- as.numeric(as.character(output_table[,"Before_Size"]))
    output_table[,"After_Size"] <- as.numeric(as.character(output_table[,"After_Size"]))
    output_table[,"T_PValue"] <- as.numeric(as.character(output_table[,"T_PValue"]))
    output_table[,"Row"] <- as.numeric(as.character(output_table[,"Row"]))
  }else{
    output_table <- output_table
  }
  return(list(output_table = output_table))
}
Format.modify <- function(ChangePointTable){
  if(nrow(ChangePointTable)!=0L){
    ChangePointTable[,"Before Mean"] <- formatC(ChangePointTable[,"Before Mean"], format = "e", digits = 3)
    ChangePointTable[,"After Mean"] <- formatC(ChangePointTable[,"After Mean"], format = "e", digits = 3)
    ChangePointTable[,"Mean_Diff"] <- formatC(ChangePointTable[,"Mean_Diff"], format = "e", digits = 3)
    ChangePointTable[,"Before Sigma"] <- formatC(ChangePointTable[,"Before Sigma"], format = "e", digits = 3)
    ChangePointTable[,"After Sigma"] <- formatC(ChangePointTable[,"After Sigma"], format = "e", digits = 3)
    ChangePointTable[,"Before Sigma"] <- formatC(ChangePointTable[,"Before Sigma"], format = "e", digits = 3)
    ChangePointTable[,"After Sigma"] <- formatC(ChangePointTable[,"After Sigma"], format = "e", digits = 3)
    ChangePointTable[,"Before Median"] <- formatC(ChangePointTable[,"Before Median"], format = "e", digits = 3)
    ChangePointTable[,"After Median"] <- formatC(ChangePointTable[,"After Median"], format = "e", digits = 3)
    ChangePointTable[,"T_PValue"] <- formatC(ChangePointTable[,"T_PValue"], format = "e", digits = 3)
    ChangePointTable[,"Before"] <- formatC(ChangePointTable[,"Before"], format = "e", digits = 3)
    ChangePointTable[,"After"] <- formatC(ChangePointTable[,"After"], format = "e", digits = 3)
    ChangePointTable[,"Difference"] <- formatC(ChangePointTable[,"Difference"], format = "e", digits = 3)
    ChangePointTable <- ChangePointTable[,c("No","Start time","Cut Point","End time","Before Mean","After Mean","Mean_Diff",
                                            "Before Median","After Median","Before Sigma","After Sigma","Before Size","After Size",
                                            "Before","After","T_PValue","Difference","Row")]
  }else{
    ColName <- c("No","Start time","Cut Point","End time","Before Mean","After Mean","Mean_Diff",
                 "Before Median","After Median","Before Sigma","After Sigma","Before Size","After Size",
                 "Before","After","T_PValue","Difference","Row")
    ChangePointTable <- data.frame(matrix(ncol = length(ColName),nrow = 0))
    colnames(ChangePointTable) <- ColName
  }
  return(list(ChangePointTable = ChangePointTable))
}
CPA.Rank <- function(output_table){
  if(nrow(output_table)!=0L){
    output_table <- output_table[order(abs(output_table[,"Difference"]),decreasing = T),]
    output_table$No <- 1:nrow(output_table)
  }
  return(list(output_table = output_table))
}

input <- "C:\\Users\\cathy.chen\\Desktop\\CPA_New_Version\\File"
output <- "C:\\Users\\cathy.chen\\Desktop\\CPA_New_Version\\Csv for result"
YName = "Y" ;Time <- "TIME"; Para.change.number <- 10 ;method = "Fitler"; multiple = 1.5; AfterSize_threshold = 25; Pvalue_threshold = 0.05

data <- read.csv(file.path(input,"Case_1.csv"),header = T,check.names = F)
output_table <- CPA.Flow(input = input,data = data,YName = YName,Time = Time,Para.change.number = 10,multiple = 1.5,AfterSize_threshold = 25,Pvalue_threshold = 0.05)$output_table
output_table
#write.csv(output_table,file.path(output,"Case1_result.csv"),row.names = F)


