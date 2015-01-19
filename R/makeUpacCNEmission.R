setwd("/Users/sohrab/project/conifercp/src/main/resources/")

# add JSON structure
result <- paste0('{\n  "caseSensitive" : false,\n  "orderedSymbols" : [ ')

generate <- function(i1, i2, someString) 
{  
  s <- paste0(someString, '"', i1, ',', i2,  '"')
  if (i1 * 5 + i2 == 18) s <- paste0(s, "\n  ")
  if (i1 != 5 | i2 != 5) s <- paste0(s, ", ")
  s
}

# add the states
for (i in 0:5) { for (j in 0:5) { result <<- generate(i, j, result)}}
#lapply(0:5, function(i) { lapply(0:5, function(j) { result <<- generate(i, j, result)  })  } )
result <- paste0(result, " ],\n")

# add ambiguity
result <- paste0(result, '  "ambiguousSymbols": {\n')
result <- paste0(result, '"-,1": [ "1,1", "2,1", "3,1" ] \n}\n}')

cat(result)

writeLines(result, "conifer/io/cn-emission-iupac-encoding.txt")