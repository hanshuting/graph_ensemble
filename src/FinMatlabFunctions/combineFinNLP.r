ptm = proc.time()

nlp_hrly <- read.csv("C:/Users/310203736/Dropbox/DarpaProject/fwMatch/data/combined_data/nlp_hrly_jul16.csv")
nlp_hrly$X = as.character(nlp_hrly$X) 
  
 FinEconBnryHrly <- read.csv("C:/Users/310203736/Dropbox/DarpaProject/fwMatch/data/combined_data/FinEconBnryHrly(Problematic).csv")
 FinEconBnryHrly$Time.Stamps = as.character(FinEconBnryHrly$Time.Stamps)
 
 View(nlp_hrly)
 View(FinEconBnryHrly)
 
 nlp_hrly$X = gsub("-","/",nlp_hrly$X)
 colnames(nlp_hrly)[1] = "Time_Stamps"
 colnames(FinEconBnryHrly)[1] = "Time_Stamps"
 
 
 
finTimesStr = as.POSIXct(FinEconBnryHrly[,"Time_Stamps"],"%Y/%m/%d %H:%M", tz = "EST")
nlpTimesStr = as.POSIXct(nlp_hrly[,"Time_Stamps"],"%Y/%m/%d %H:%M:%S", tz = "EST")

finTimes = as.numeric(finTimesStr)
nlpTimes = as.numeric(nlpTimesStr)

FinEconBnryHrly$Time_Stamps = finTimes
nlp_hrly$Time_Stamps = nlpTimes
 
###It's checked that the order of columns are preserved in this re-ordering process. 
finNLPCombo = merge(FinEconBnryHrly, nlp_hrly, by = "Time_Stamps", all = TRUE)

comboFinCols = colnames(finNLPCombo)[2:892]
comboNlpCols = colnames(finNLPCombo)[893:1112]

#Add state NaN to NLP variables.
firstNewNLPName = gsub("1","NaN",comboNlpCols[1])
newNLPDataFrame = cbind(data.frame(firstNewNLPName = rep(0,nrow(finNLPCombo))),finNLPCombo[,893:(893+9)])
colnames(newNLPDataFrame)[1] = firstNewNLPName

#Handle the rest
for (i in 2:22){
thisNewNLPName = gsub("1","NaN",comboNlpCols[(i-1)*10 + 1])
newNLPDataFrame = cbind(newNLPDataFrame,data.frame(thisNewNLPName = rep(0,nrow(finNLPCombo))))
colnames(newNLPDataFrame)[ncol(newNLPDataFrame)] = thisNewNLPName

newNLPDataFrame = cbind(newNLPDataFrame,finNLPCombo[,(893 + 10*(i-1)):(902 + 10*(i-1))])

}


newFinNLPCombo = cbind(finNLPCombo[1:892], newNLPDataFrame)



#####Missing data, switch the state to NaN########################


for (t in 1:nrow(newFinNLPCombo)){
print(t/nrow(newFinNLPCombo))
###Handle FIN first.
###FIN columns are 2~892
ptr = 2
while(ptr <= 892){

#sufficient to check one
if (is.na(newFinNLPCombo[t,ptr+1])){
newFinNLPCombo[t,ptr] = 1 #turn on NaN
newFinNLPCombo[t,(ptr+1):(ptr+10)] = 0     #turn off other states

}

ptr = ptr + 11
}

###Handle NLP
#ptr = 893

while(ptr <= ncol(newFinNLPCombo)){
#sufficient to check one
if (is.na(newFinNLPCombo[t,ptr + 1])){
newFinNLPCombo[t,ptr] = 1 #turn on NaN
newFinNLPCombo[t,(ptr+1):(ptr+10)] = 0     #turn off other states

}

ptr = ptr + 11

}

}

proc.time() - ptm