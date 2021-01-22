

###############

#investigation: i = 64
identicalLists(mm_original, mm_dev, 
               ignore = c("call", "FunList", "dropped"), chatty = TRUE)


mm_original$clusz
mm_dev$clusz

mm_original$offset
mm_dev$offset

length(mm_original$eta)
length(mm_dev$eta)

all(mm_original$eta[5:444] == mm_dev$eta)

mm_original$weights
mm_dev$weights


res[i, ]
#############
#Look at combos with error for dev but not original
res[res$error_new & !res$error_original,] #none! 

#Look at combos with error for original, but not dev
res[!res$error_new & res$error_original, c("corstr", "weights", "treat", "waves")] 

#Look at combos with error for original
View(res[res$error_original, c("corstr", "weights", "treat", "waves")]) 


#Look at combos without error for original where results are not identical
View(res[!res$error_original & !res$identical,])
