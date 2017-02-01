
#Check if every 11 columns of each row sums up to 1.

checkResults = c()

ptr = 2

while (ptr <= ncol(comboTest2)){


checkResults = c(checkResults, ifelse(all( apply(comboTest2[,ptr:(ptr+10)], 1,sum)==1),0,1))

ptr = ptr + 11


}