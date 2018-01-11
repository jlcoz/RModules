#!/biosw/debian7-x86_64/R/3.2.2/bin/Rscript --vanilla

## Retrieve arguments
args=commandArgs(TRUE)

## Smallest number before p-values are set to 0
## this is necessary as other programs (e.g. awk) report WRONG results
MINNUM=2.2250738585072014e-308

## Help
help <- function(){
  cat("\nAdjust the p-values of a given column of a file, and output this file added the adjusted p-value column\n")
  cat("Usage: adjustPval.R -i - -f F -c N -m BH -o -\n")
  cat("-i : Input table as a file or stdin (-) [Required]\n")
  cat("-f : First line of the input table is a header : T/F [Default: T]\n")
  cat("-c : Column with raw p-values (int) [Required]\n")
  cat("-m : P-value adjust method used : BH or bonferroni [Default: BH]\n")
  cat("-o : Output as a file or stdout (-) [Default: stdout]\n")
  cat("\n")
  q()
}

## Save values of each argument
## In case of args written but not initialiazed, value is set to NA
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

## Set the head boolean to load the file accordingly
## Input is considered without header by default
if(exists("f")){
  if(is.na(f)){
    head=F
  } else if (f=="T"){
    head=T
  } else {head = F}
} else{ head = F}

## Load the matrix into a dataframe
if(exists("i")){
  if(is.na(i)){
    cat("Input file not specified\n"); q()  
  }
  if(i=="stdin" || i=="-"){
    count_table=read.csv(pipe('cat /dev/stdin'), sep="\t", skip=0, header = head, comment.char = "", check.names = F)
  } else if (file.exists(i)){
    count_table=read.csv(i, sep="\t", skip=0, header = head, comment.char = "", check.names = F)
  } else { cat("Input file does not exist\n"); q() }
} else { cat("No input specified\n"); q() }


## Check if column is specified
if(exists("c")){
  
  col=strtoi(c)

  if(is.na(col) || col<1){ cat("P-value column is not int>=1 (-c)\n"); q() }

  ##Test if the matrix has at least col columns
  if(ncol(count_table)<col){
    cat("The selected column is out of range\n"); q()
  }
  
  ## Test column c to see if they contain only numerics
  if(!all(sapply(count_table, function(x) class(x) %in% c("numeric"))[col])){
    cat("Column c does not contain only numerics\n");q()
  }

} else { cat("P-value column not specified (-c)\n"); q() }


## Check the padjust method selected
if(exists("m")){
  ## If the selected methode is not a choice of the p.adjust methods
  if(!m %in% c("BH","bonferroni")){
    cat("Please enter a valid method\n");q()
  }
}else{m="BH"}

## Set the ouput path : File or STDOUT
if(exists("o")){
  if(o=="stdout" || o=="-"){
    output=stdout()
  } else {output=o}
} else { output=stdout() }

## Adjust p-value accordingly to the method chosed by the user
adj_pval <- p.adjust(count_table[[col]], method = m)
adj_pval[adj_pval<=MINNUM]=0

count_table$adj_pval=adj_pval

## Writing the results
## Print the header accordingly to -f parameter
write.table(count_table,
            file=output,
            sep="\t",
            quote=FALSE, 
            row.names=FALSE, 
            col.names=head)