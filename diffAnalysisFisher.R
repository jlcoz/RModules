# define custom functions
rowMatrix <- function(x1,x2,x3,x4) {
  return(matrix(c(x1,x2,x3,x4),nrow=2))
}
rowFisherP <- function(x, ...) {
  return(fisher.test(x)$p.value)
}

## Retrieve arguments
args=commandArgs(TRUE)

## Smallest number before p-values are set to 0
## this is necessary as other programs (e.g. awk) report WRONG results
MINNUM=2.2250738585072014e-308

## Help
help <- function(){
  cat("\ndiffAnalysisFisher.R : Retrieve differential peaks from a count matrix\n")
  cat("Usage: diffAnalysisFisher.R -i - -a N1 -b N2 -o -\n")
  cat("-i : Count matrix as a file or stdin (-) [Required]\n")
  cat("-a : Total number of reads for sample 1 : N1 [Required]\n")
  cat("-b : Total number of reads for sample 2 : N2 [Required]\n")
  cat("-o : Output as a file or stdout (-) [Default: stdout]\n")
  cat("\n")
  q()
}

## Save values of each argument
if(length(args)==0 || !is.na(charmatch("-help",args))){
  help()
} else {
  for(ii in 1:length(args)){
    if(grepl("^-",args[ii]) && args[ii] != "-"){
      if(ii+1<=length(args) && (!grepl("^-",args[ii+1]) || args[ii+1]=="-")){
        assign(gsub("-","",args[ii]),args[ii+1])
      } else {assign(gsub("-","",args[ii]),NA) }
    }
  }}

## Load the matrix into a dataframe
## Load the table into a dataframe
if(exists("i")){
  if (is.na(i)){
    cat("Input file does not exist\n"); q()
  } else if(i=="stdin" || i=="-"){
    count_table=read.csv(pipe('cat /dev/stdin'), sep="\t", skip=0, header = T, comment.char = "", check.names = F)
  } else if (file.exists(i)){
    count_table=read.csv(i, sep="\t", skip=0, header = T, comment.char = "", check.names = F)
  }
  ## Test the second-to-end column to see if they contain only integers
  if(!all(sapply(count_table, function(x) class(x) %in% c("integer"))[-1])){
    cat("The counts does not contain only integers\n");q()
  }
} else { cat("No input specified\n"); q() }

## Load the total number of reads for each of the 2 samples
if(exists("a") && exists("b")){
  if(!is.na(a) && !is.na(b)){
    N1=strtoi(a)
    N2=strtoi(b)
    if(N1<max(count_table[-1]) || N2<max(count_table[-1])){
      cat("Total number of reads is inferior to the highest count in the table\n"); q()
    }
  } else{ cat("Specify total number of reads for both samples\n"); q() }
} else { cat("Specify total number of reads for both samples\n"); q() }

## Set the ouput path : File or STDOUT
if(exists("o")){
  if(o=="stdout" || o=="-"){
    output=stdout()
  } else {output=o}
} else { output=stdout() }

## Start of the analysis

## Design contingency table as matrix
## [,1]           [,2]
## [1,] count1    count2
## [2,] N1-count1 N2-count2

contingency_matrix <- lapply(c(1:nrow(count_table)), function(x) {y <- count_table[x,]; matrix(unlist(c(y[,2],N1-y[,2],y[,3],N2-y[,3])),nrow=2)})

## Run fisher's exact on contingency table
pval <- sapply(contingency_matrix, rowFisherP)

## Adjust p-value using BH
adj_pval <- p.adjust(pval, method = "BH")

## Calculate fold change

logFC <- apply(count_table, 1, function(x) log(strtoi(x[[3]])/strtoi(x[[2]]))/log(2))

## Shaping the result dataframe
total=count_table
total$logFC=logFC
total$pval=pval
total$FDR=adj_pval
names(total)[[1]]="#peaks"

## Order the results by pvalues
total=total[order(total$pval),]

## Filter the pvalues under the minimum value
total$pval[total$pval<=MINNUM] = 0
total$FDR[total$FDR<=MINNUM] = 0

## Writing down the results : counts logFC pval adj_pval (FDR)
write.table(total,
            file=output,
            sep="\t",
            quote=FALSE,
            row.names = FALSE,
            col.names = TRUE)