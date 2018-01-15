#!/biosw/debian7-x86_64/R/3.2.2/bin/Rscript --vanilla

suppressWarnings(suppressMessages(library(edgeR)))
suppressWarnings(suppressMessages(library(limma)))
suppressWarnings(suppressMessages(library(tools)))

## Smallest number before p-values are set to 0
## this is necessary as other programs (e.g. awk) report WRONG results
MINNUM=2.2250738585072014e-308

## Retrieve arguments
args=commandArgs(TRUE)

## Help
help <- function(){
  cat("\ndiffAnalysisEdgeR_GLM.R : Retrieve differential peaks from a count matrix using general linear model methods\n")
  cat("Usage: diffAnalysisEdgeR.R -i - -f F -a n1 -b n2 -n x1,x2... -o -\n")
  cat("-i : Count table as a file or stdin (-) [Required]\n")
  cat("-f : First line of the input table is a header : T/F [Default: F]\n")
  cat("-a : Number of samples in the first condition [Required]\n")
  cat("-b : Numbre of samples in the second condition [Required]\n")
  cat("-n : Normalization factor as a vector : x1,x2,... [Default: EdgeR computation]\n")
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

## Set the head boolean to load the file accordingly
## Input is considered without header by default
if(exists("f")){
  if(is.na(f)){
    head=F
  } else if (f=="T"){
    head=T
  } else {head = F}
} else{ head = F}

## Load the table into a dataframe
if(exists("i")){
  if (is.na(i)){
    cat("Input file does not exist\n"); q()
  } else if(i=="stdin" || i=="-"){
    count_table=read.csv(pipe('cat /dev/stdin'), sep="\t", skip=0, header = head, comment.char = "", check.names = F)
  } else if (file.exists(i)){
    count_table=read.csv(i, sep="\t", skip=0, header = head, comment.char = "", check.names = F)
  }
  ## Test the second-to-end column to see if they contain only integers
  if(!all(sapply(count_table, function(x) class(x) %in% c("integer"))[-1])){
    cat("The counts does not contain only integers\n");q()
  }
} else { cat("No input specified\n"); q() }

## Construct the design vector
if(exists("a") & exists("b")){
  if(!is.na(a) && !is.na(b)){
    nb_samples_condition1=strtoi(a)
    nb_samples_condition2=strtoi(b)
  } else{ cat("Specify the number of samples for both condition\n"); q()}
  ## Test if the matrix as the same number of columns as number of samples
  ## The +1 represents the "peaks" column
  if(ncol(count_table)==nb_samples_condition1+nb_samples_condition2+1){
    condition = factor(c(rep(1,nb_samples_condition1),rep(2,nb_samples_condition2)))
  } else { cat("The experimentals design does not match the count matrix\n");q()}
  
} else { cat("No design specified\n"); q() }

## Construct the normalization factor vector
if(exists("n")){
  ## Args are stored in a list ordered like : 1 arg1 value1
  ##                                          3 arg2 value2
  vector=args[which(args=="-n")+1]
  norm_factor_vector = as.numeric(strsplit(vector, ",")[[1]])
  #The NF boolean notify that NF computation is not needed
  NF=T
  if(length(norm_factor_vector)!=nb_samples_condition1+nb_samples_condition2){
    cat("Length of vector different of number of samples\n"); q()
  }
} else {NF=F}

## Set the ouput path : File or STDOUT
if(exists("o")){
  if(o=="stdout" || o=="-"){
    output=stdout()
  } else {output=o}
} else {output=stdout()}

## Start of the analysis
## EdgeR requires only numbers within the matrix
## Set peaks as rownames...
rownames(count_table) <- count_table[[1]]
## ...and then remove the peak column
count_table[[1]] = NULL

## Create an edgeR object required for the analysis
dge <- DGEList(counts=count_table,group=condition)

## If no normalization factors are given by parameters, TMM computation
if(!NF){
  dge=calcNormFactors(dge)
} else{dge$samples$norm.factors = norm_factor_vector}

## Create a design matrix needed for the dispersion estimation
design <- model.matrix(~condition)
colnames(design)<-levels(design)

## Estimation of the dispersion parameter
dge <- estimateDisp(dge, design)

dge <- estimateGLMCommonDisp(dge,design)  #Estimates a common negative binomial dispersion
dge <- estimateGLMTrendedDisp(dge,design) #Estimates the abundance-dispersion trend by Cox-Reid approximate profile likelihood
dge <- estimateGLMTagwiseDisp(dge,design) #Compute an empirical Bayes estimate of the negative binomial dispersion parameter for each transcript, with expre
fit <- glmFit(dge,design)    #Fit a negative binomial generalized log-linear model to counts
lrt <- glmLRT(fit,coef=fit$design) #Conducts likelihood ratio tests for one or more coefficients in the linear model

## Multiple testing correction : Benjamimi Hochberg method
adj_pval=p.adjust(lrt$table$PValue,method = "BH")

edgeR_results = data.frame(peaks=rownames(lrt),
                           logFC=lrt$table$logFC.1,
                           pval=lrt$table$PValue,
                           FDR=adj_pval)

## ...and then merge both dataframes by the peaks ID
total=merge(count_table, edgeR_results, by.x="row.names", by.y="peaks")

## Rename the column Row.names
names(total)[[1]]="#peaks"

## Order the table by pvalue
total=total[order(total$pval),]

## Filter the pvalues under the minimum value for awk
total$pval[total$pval<=MINNUM] = 0
total$FDR[total$FDR<=MINNUM] = 0

## Writing down the results : counts peaks logFC adj_pval (FDR) values
write.table(total,
            file=output,
            sep="\t",
            quote=FALSE,
            row.names = FALSE,
            col.names = TRUE)
