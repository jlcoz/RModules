#!/biosw/debian7-x86_64/R/3.2.2/bin/Rscript --vanilla

suppressWarnings(suppressMessages(library("DESeq2")))

## Smallest number before p-values are set to 0
## this is necessary as other programs (e.g. awk) report WRONG results
MINNUM=2.2250738585072014e-308

## Retrieve arguments
args=commandArgs(TRUE)

## Help
help <- function(){
  cat("\ndiffAnalysisDESeq2.R : Retrieve differential peaks from a count matrix\n")
  cat("Usage: diffAnalysisDESeq2.R -i - -a n1 -b n2 -n x1,x2... -o - -F T -1 F\n")
  cat("-i : Count matrix as a file or stdin (-) [Required]\n")
  cat("-a : Number of samples in the first condition [Required]\n")
  cat("-b : Numbre of samples in the second condition [Required]\n")
  cat("-n : Normalization factor as a vector : x1,x2,... [Default: DESeq2 computation\n")
  cat("-o : Output as a file or stdout (-) [Default: stdout]\n")
  cat("-f : Peaks with a mean normalized count < optimal threshold get NA as adjusted pvalue : T/F [Default: T]\n")
  cat("-1 : Replace NA with 1 (pvalue/adjusted pvalue) and 0 (logFC) : T/F [Default: F]\n")
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
  if (i==1){
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

## Construct the design vector
if(exists("a") & exists("b")){
  nb_samples_condition1=strtoi(a)
  nb_samples_condition2=strtoi(b)
  ## Test if the matrix as the same number of columns as number of samples
  ## The +1 represents the "peaks" column
  if(ncol(count_table)==nb_samples_condition1+nb_samples_condition2+1){
    condition = factor(c(rep(1,nb_samples_condition1),rep(2,nb_samples_condition2)))
  } else { cat("The experimentals design does not match the count matrix\n");q()}
  
} else { cat("No design specified\n"); q() }

## Construct the normalization factor matrix
if(exists("n")){
  ## Args are stored in a list ordered like : 1 arg1 value1
  ##                                          3 arg2 value2
  vector=args[which(args=="-n")+1]
  norm_factor_vector = as.numeric(strsplit(vector, ",")[[1]])
  ## The NF boolean notify that NF computation is not needed
  NF=T
  if(length(norm_factor_vector)!=nb_samples_condition1+nb_samples_condition2){
    cat("Length of vector different of number of samples\n"); q()
  }
} else {norm_factor_vector = c(rep(1,(nb_samples_condition1+nb_samples_condition2))) ; NF=F}

## Set the ouput path : File or STDOUT
if(exists("o")){
  if(o=="stdout" || o=="-"){
    output=stdout()
  } else {output=o}
} else {output=stdout()}

## Set the independant filtering boolean for the DESeq2 p.adjust method
if(exists("f")){
  if (F=="F") {
    independant_filtering=F
  }
} else {independant_filtering=T}

##TODO 
##IMPLEMENT THE MINUS 1 OPTION FOR REPLACING NA

## Start of the analysis
## DESeq2 requires only numbers within the matrix
## Set peaks as rownames...
rownames(count_table) <- count_table[[1]]
## ...and then remove the peak column
count_table[[1]] = NULL

##Convert the count table as a count matrix
matrix_counts=as.matrix(count_table)

##Create a dataframe linking each sample with its normalization factor
colData=data.frame(norm_factor_vector, condition)

## Create a DESeq2 object required for the analysis
dds <- DESeqDataSetFromMatrix(matrix_counts,
                              colData = colData,
                              design = ~ condition)

## Create the normalization factor matrix
normFactors <- matrix(colData$norm_factor_vector,
                      ncol=ncol(dds),nrow=nrow(dds),
                      dimnames=list(1:nrow(dds),1:ncol(dds)),
                      byrow = TRUE)

## If no normalization factors are given by parameters, use DESeq2 computation method
if(NF){
  normalizationFactors(dds) = normFactors
} else{dds=estimateSizeFactors(dds)}

## Estimation of the dispersion parameter
## Use the local fitTyp to fit a local regression of log dispersions over log base mean
dds <- estimateDispersions(dds, quiet=T, fitType="parametric")

## Negative Binomial GLM fitting and Wald statistics
dds <- nbinomWaldTest(dds, quiet=T)

## Multiple testing correction : Benjamimi Hochberg method
if(independant_filtering){
  resDESeq2 <- results(dds, pAdjustMethod = "BH", independentFiltering =T)
} else {
  resDESeq2 <- results(dds, pAdjustMethod = "BH", independentFiltering =F)
}

## Shape the final result table
total = data.frame(rownames(count_table),
                   count_table,
                   logFC = resDESeq2$log2FoldChange,
                   pval = resDESeq2$pvalue,
                   FDR = resDESeq2$padj)
  
## Rename the column Row.names to #
names(total)[[1]]="#peaks"

## Order the table by pvalue
total=total[order(total$pval),]

## Filter the pvalues under the minimum value
total$pval[total$pval<=MINNUM] = 0
total$FDR[total$FDR<=MINNUM] = 0

## Writing down the results : counts peaks logFC pval adj_pval (FDR) values
write.table(total,
            file=output,
            sep="\t",
            quote=FALSE,
            row.names = FALSE,
            col.names = TRUE)