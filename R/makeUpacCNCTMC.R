setwd("/Users/sohrab/project/conifercp/src/main/resources/")

generate <- function(i1, i2, i3, someString) 
{  
  s <- paste0(someString, '"', i1, ',', i2, ',', i3, '"')
  if (i1 * 5 + i2 == 18) s <- paste0(s, "\n  ")
  if (i1 != 5 | i2 != 5) s <- paste0(s, ", ")
  s
}

# add JSON structure
result <- paste0('{\n  "caseSensitive" : false,\n  "orderedSymbols" : [ ')

# add the states in \Sigma_0 (a=0, b=0)
for (i in 0:5) result <<- generate(i, 0, 0, result)

# add the states in \Sigma_1 (b=1)
for (i in 0:5) { for (j in 0:5) { result <<- generate(i, j, 1, result)}}

# close the array
result <- paste0(result, " ],\n")

# add ambiguity
result <- paste0(result, '  "ambiguousSymbols": {\n')
result <- paste0(result, '"-,1,0": [ "1,1,0", "2,1,0", "3,1,0" ] \n}\n}')

cat(result)

writeLines(result, "conifer/io/cn-ctmc-iupac-encoding.txt")
