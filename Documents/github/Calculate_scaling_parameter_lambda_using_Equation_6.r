##########################################################################################################################
### Note: R_packages: "Biostring" and "PWMErich" should be installed before running the code below.
### To install these two packages:
### install.packages("Biostring")
### library("Biostring")
### install.packages("PWMErich")
### library("PWMErich")
###
### Lambda_estimation_simple_equation is a function to estimate the scaling parameter lambda to a specific PWM matrix by
### simply using background genomic sequences. It will return a 2-column matrix with the first column contains TF names and 
### the second column contains the estimated lambda values. 
###
### There are four input parameters to this function as listed below:
### 1. PFM_folder should specify the directory of the folder where PFMs(Position-specifc Frequency Matrices)
### are stored, e.g.  PFM_folder="~/home/PFM". Individual files for PFMs need to be text files in the following format:
### 1) first line specify the TF name or PFM matrix name;
### 2) 2-5th line is the Tab deliminated frequency of nucleotide A,C,G,T respectively in each position of the PFM motif.
### 3) the motif should be at least 7 bp long in length by using the default threshold parameters.
### following is an example:
### 
### CTCF
### 306  313	457	676	257	1534	202	987	2	0	2	124	1	79	231
### 876	1147	383	784	714	1	0	0	4	0	0	1645	0	1514	773
### 403	219	826	350	87	192	1700	912	311	1902	1652	3	1807	8	144
### 317	223	236	92	844	175	0	3	1585	0	248	130	94	301	754
### 
### There is no requirement for individual file name. 
### 
### 2.pseudo_count is used for adjusting PWM matrInformationContentes to avoid the logrithm of zero. 
### By default, you can put pseudo_count=1
### 
### 3.background_seq is the background genome sequences used, note it should be a DNAString object.
### 
### 4. GC content is the average GC content for the given input DNA sequence
##########################################################################################################################

Lambda_estimation_simple_equation<-function(PFM_folder,pseudo_count,background_seq,GC_content) {
### lambda_table will store TF names and their estimated lambda values. 
lambda_table=c()

### the default PWM top percentage threshold is defined as top 0.1% to represent the weakest binding sites,
### if there is other alternative experimental evidences suggesting a different threshold of binding for
### a specifInformationContent TF, please change the value below. 
top_score_threshold=0.001

Reverse_complement=reverseComplement(background_seq)

folder=dir(PFM_folder)
for(file in folder){ 
  lambda=c()
  fullpath=paste(PFM_folder,"/",file,sep="")
  PFM=read.table(fullpath,header=FALSE,stringsAsFactors=FALSE,skip=1,sep="\t")
  
  ### motif length need to be longer than 7 base-pair for default threshold, if changing the default threshold
  ### to be more strInformationContentt, the length threshold need to be changed accordingly. (e.g. if using 0.01% threshold,
  ### motifs should be >= 9bp, and if using 0.001%, motifs should be >= 10bp or 11bp). 
  if(ncol(PFM)>=7){
    
    ### read-in TF name.
    TF=readLines(fullpath, n=1)
    
    ### construct PWM(Position Weight Matrix) from PFM.
    PWM_matrix=matrix(0,ncol=ncol(PFM),nrow=4)
    for (j in 1:4){
       for (i in 1:ncol(PWM_matrix)){
         sum=PFM[1,i]+PFM[2,i]+PFM[3,i]+PFM[4,i]
         if (j==1 |j==4){
           PWM_matrix[j,i]=log2((PFM[j,i]+(1-GC_content)/2*pseudo_count)/(sum+pseudo_count)/((1-GC_content)/2))
          }
         else{PWM_matrix[j,i]=log2((PFM[j,i]+GC_content/2*pseudo_count)/(sum+pseudo_count)/(GC_content/2))}
       }
    }
    rownames(PWM_matrix)=c("A","C","G","T")
    
    ### motif length l:
    l=length(PWM_matrix)/4
    
    ### calculate the background PWM score for the input background sequences for one strand 
    BackGroundPWMScore=PWMscoreStartingAt(PWM_matrix,background_seq,1:(length(background_seq)-l+1))
    
    ### calculate the background PWM score for the input background sequences for another strand 
    BackGroundPWMScore_rev=rev(PWMscoreStartingAt(PWM_matrix,Reverse_complement,1:(length(background_seq)-l+1)))
    
    ### assign PWM score to each position of the genome
    BackGroundPWMScore=rbind(BackGroundPWMScore,BackGroundPWMScore_rev)
    BackGroundPWMScore=apply(BackGroundPWMScore,2,max)
    
    ### find the maximum PWM score   
    Max_PWMScore=max(BackGroundPWMScore)
    
    ### calculate the information content for the PWM matrix
    InformationContent=0
    for (j in 1:4){
      for (i in 1:ncol(PWM_matrix)){
        sum=PFM[1,i]+PFM[2,i]+PFM[3,i]+PFM[4,i]
        InformationContent=InformationContent+PWM_matrix[j,i]*PFM[j,i]/sum
      }
    }
   
    ### estimate the mismatch bits tolerance by the TF, 6 and 13.2 are parameter values estimated from Maerkl et al., 
    ### 2007, if there are further evidence for mismatch energy estimation for specific TF or TF families, please
    ### change these values below.
    Mismatch_bits=6*InformationContent/13.2
    
    ### calculate lambda values based on Equation 6 of Ma et al.,2015
    lambda=c(lambda,(Max_PWMScore-quantile(pwm_back,1-top_score_threshold))/Mismatch_bits)
    lambda_table=rbind(lambda_table,c(TF,lambda))
  
  }
  else {
    print("Motif input problem, please check motif length.")
   }
 }
colnames(lambda_table)=c("TF","lambda")
Lambda_estimation_simple_equation=lambda_table
}
# PWM_scale
