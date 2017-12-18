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
  cat("Usage: diffAnalysisFisher.R -i - -n N1,N2 -o -\n")
  cat("-i : Count matrix as a file or stdin (-)\n")
  cat("-n : Total number of reads : N1,N2\n")
  cat("-o : Output as a file or stdout (-)\n")
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
      } else {assign(gsub("-","",args[ii]),1) }
    }
  }}

## Load the matrix into a dataframe
if(exists("i")){
  if(i=="stdin" || i=="-"){
    matrix_counts=read.csv(pipe('cat /dev/stdin'), sep="\t", skip=0, header = T, comment.char = "", check.names = F)
  } else if (file.exists(i)){
    matrix_counts=read.csv(i, sep="\t", skip=0, header = T, comment.char = "", check.names = F)
  } else { cat("Input file does not exist\n"); q() }
  
  ##Test if the matrix has only 3 columns
  if(ncol(matrix_counts)!=3){
    cat("The matrix is not 3 columns formatted\n"); q()
  }
  
  ## Test the second-to-end column to see if they contain only numbers
  if(!all(sapply(matrix_counts, function(x) class(x) %in% c("integer","numeric"))[-1])){
    cat("The matrix does not contain only numbers\n");q()
  }
} else { cat("No input specified\n"); q() }

## Get the total number of reads for the two samples
if(exists("n")){
  ## Args are stored in a list ordered like : 1 arg1 value1
  ##                                          3 arg2 value2
  vector=args[which(args=="-n")+1]
  vector_total_number_reads = as.numeric(strsplit(vector, ",")[[1]])
  if(length(vector_total_number_reads)!=2){
    cat("Please enter only 2 numbers\n"); q()
  }
} else {cat("Please give the total numbers of reads for the 2 samples\n"); q()}

## Set the ouput path : File or STDOUT
if(exists("o")){
  if(o=="stdout" || o=="-"){
    output=stdout()
  } else {output=o}
} else { output=stdout() }

## Load the total number of reads for each of the 2 samples
N1=vector_total_number_reads[[1]]
N2=vector_total_number_reads[[2]]

## Design contingency table as matrix
contingency_matrix <- lapply(c(1:nrow(matrix_counts)), function(x) {y <- matrix_counts[x,]; matrix(unlist(c(y[,2],N1-y[,2],y[,3],N2-y[,3])),nrow=2)})

## Run fisher's exact on contingency table
pval <- sapply(contingency_matrix, rowFisherP)

## Adjust p-value using BH
adj_pval <- p.adjust(pval, method = "BH")

## Calculate fold change
logFC <- apply(matrix_counts, 1, function(x) log(strtoi(x[[3]])/strtoi(x[[2]]))/log(2))

## Shaping the result dataframe
total=matrix_counts
total$logFC=logFC
total$pval=pval
total$FDR=adj_pval
names(total)[[1]]="#peaks"

## Order the results by pvalues
total=total[order(total$pval),]

## Filter the pvalues under the minimum value
total$pval[total$pval<=MINNUM] = 0
total$FDR[total$FDR<=MINNUM] = 0

## Remove the possible NA created by the p.adjust method
total = na.omit(total)

## Writing down the results : counts logFC pval adj_pval (FDR)
write.table(total,
            file=output,
            sep="\t",
            quote=FALSE,
            row.names = FALSE,
            col.names = TRUE)