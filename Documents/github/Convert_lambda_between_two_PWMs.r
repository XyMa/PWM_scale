##########################################################################################################################
### Note: R_packages: "Biostring", "PWMErich", "Metrics" and "plot3D"(If you want to see the heatmap) should be installed before running the code below.
### To install these two packages:
### install.packages("Biostring")
### library("Biostring")
### install.packages("PWMErich")
### library("PWMErich")
### install.packages("Metrics")
### library("Metrics")
### install.packages("plot3D")
### library("plot3D")
###
### Convert_lambda_between_two_PWMs is a function to calculate the situable value of lambda for a new PWM matrix of the 
### same transcription factor based on one known value of lambda but for another PWM matrix. 
### 
###
### There are 6 input parameters in this function as listed below:
### 1. PFMSet1 and PFMSet2 should be lists of PFM matrices (Position-specifc Frequency Matrices) 
###    with individual element names in the format of "TFname_PFMname".
### Note 
### 1) the order of PFMs in Set1 and Set2 should match; 
### 2) TFname (first half of each element name before "_") in Set1 and Set2 should match.
### 3) each PFM matrix should contain the frequency of nucleotide A,C,G,T in row 1 to row 4, respectively, and each column represents 
### each nuceltide position in the PFM motif;
### 4) motifs at least 9 bp long in length in both set of PFMs are preferred
### 5) there is no requirement for colnames or rownames of each PFMs.
### 
### 2.lambda_value_PFMSet1 is a 2-column matrix with the first column saving TF names(matching with TFname in PFMSet) and the second
### colunm contains the known values for lambda corresponding to PFMSet1. These lambda values can be the ones obtained from fitting 
### experimetal data.  
###
### 3.pseudo_count is used for adjusting PWM matrInformationContentes to avoid the logrithm of zero. 
### By default, you can put pseudo_count=1
### 
### 4.background_seq is the background genome sequences used, note it should be a DNAString object.
### 
### 5. GC content is the average GC content for the given input DNA sequence
###
### This function will return a matrix with first column containing TF names and the second column containing lambda values for each 
### TF in PFMSet2. It will also generate heatmaps of log(residence time) for each pair of PFMs for visualization.
##########################################################################################################################



Convert_lambda_between_two_PWMs=function(PFMSet1,PFMSet2,lambda_value_PFMSet1,pseudo_count,background_seq,GC_content){
  
  ### lambda_PWM2 will be the matrix that stores the estimated lambda for PWM2
  lambda_PWM2=c()
  
  ### the default potential lambda range is 0.1 to 3.0, so lambda_range_times_ten is 1 to 30.
  lambda_range_times_ten=1:30
 
  ### quantile_range means the value range of -log_10(top_quantile)*10. e.g. if the top score range we choose is 
  ### top 0.1% to top 0.001%, then the -log_10(top_quantile)*10 is 30 to 50. In order to save running time, we commented out
  ### the lines to compute the PWM score in the reversecomplement strand, because for most of asymmetric motifs, top 0.1% scores for 
  ### scanning both strands are the same as top 0.05% by scanning only one strand, so the -log_10(0.05%)*10=33. For symmetric
  ### motifs, lines 87,90,91,166,169 and 170 need to be included to be more accurate, and change the quantile_range to 30:50.
  
  quantile_range=33:50
  
  Reverse_complement=reverseComplement(background_seq)
  
  for(PFM_num in 1:length(PFMSet2)){  
    PFM2=PFMSet2[[PFM_num]]
    
    ### TF name should be the first half of the PFM name
    TF=unlist(strsplit(names(PFMSet2[PFM_num]),"[_]"))[1]
    
    ### construct PWM2 from PFM2
    PWM2=matrix(0,ncol=ncol(PFM2),nrow=4)
    for (j in 1:4){
      for (i in 1:ncol(PWM2)){
        sum=PFM2[1,i]+PFM2[2,i]+PFM2[3,i]+PFM2[4,i]
        if (j==1 |j==4){
          PWM2[j,i]=log2((PFM2[j,i]+(1-GC_content)/2*pseudo_count)/(sum+pseudo_count)/((1-GC_content)/2))
        }
        else{PWM2[j,i]=log2((PFM2[j,i]+GC_content/2*pseudo_count)/(sum+pseudo_count)/(GC_content/2))}
      }
    }
    rownames(PWM2)=c("A","C","G","T")
    l<-length(PWM2)/4
    
    ### calculate the background PWM score using PWM2 for the input background sequences for one strand 
    BackGroundPWMScore=PWMscoreStartingAt(PWM2,background_seq,1:(length(background_seq)-l+1))
    
    ### calculate the background PWM score using PWM2 for the input background sequences for another strand 
    ### BackGroundPWMScore_rev=rev(PWMscoreStartingAt(PWM2,Reverse_complement,1:(length(background_seq)-l+1)))
    
    ### assign PWM score to each position of the genome
    ### BackGroundPWMScore=rbind(BackGroundPWMScore,BackGroundPWMScore_rev)
    ### BackGroundPWMScore=apply(BackGroundPWMScore,2,max)
    
    ### avg_exp is an array storing the average value for exp(BackGroundPWMScore/lambda) used in calculating
    ### the parameter tau_0 for each PWM matrix. 
    avg_exp = rep(0,times=50) 
    
    ### tau_0 is a parameter in estimating residence time based on PWM score. See Zabet et al., 2013
    tau_0=rep(0,times=50)
    
    ### tau is the estimated residence time in each position of the genome. log_tau is a matrix storing 
    ### log(tau) values for every potential lambda and each binding site strength level (top quantiles)
    log_tau=matrix(0,nrow=50,ncol=50)
    
    for(lambda in lambda_range_times_ten)
    {
      avg_exp[lambda] = mean(exp((1/(lambda*0.1))*BackGroundPWMScore));
      
      ### tR=5ms is chosen to be the average sliding time measured by Elf et al, 2007
      tR=0.005
      
      ### averge sliding length including hopping is 90 bp(Elf et al, 2007)
      sl_obs = 90;
      
      ### See Zabet et al., 2013
      tau_0[lambda] = (2*tR)/(avg_exp[lambda]*sl_obs^2)
      
      for (j in quantile_range){
        quant=quantile(BackGroundPWMScore,(1-10^(-j*0.1))) 
        log_tau[lambda,j]=log(tau_0[lambda])+quant/(lambda*0.1)
      }
    }
    
    ### residence_time_range saves the residence time range for each potential lambda value for PWM2
    residence_time_range=log_tau[lambda_range_times_ten,quantile_range];
    rownames(residence_time_range)=lambda_range_times_ten/10
    
    ### colnames of residence_time_range matrix is -log(quantile) 
    colnames(residence_time_range)=quantile_range/10
    
    residence_time_range_save=residence_time_range
    
    ### get the reference lambda value for PWM1 from lambda_value_PFMSet1
    find_TF = sapply(1:nrow(lambda_value_PFMSet1), function(i){grep(TF, lambda_value_PFMSet1[i,1])})
    find_TF_row_num=which((sapply(find_TF, length) == 1) )
    #print(find_TF_row_num)
    if (length(find_TF_row_num)>0){
      lambda_reference=lambda_value_PFMSet1[find_TF_row_num,2] 
    } else if (length(find_TF_row_num)==0){
      print("TF names don't match")
    }
    ### Following lines describe the similar process to estimate the residence time for PWM1, but it is based on the 
    ### known lambda values.  
    PFM1=PFMSet1[[PFM_num]]
    
    ### TF name should be the first half of the PFM name
    TF=unlist(strsplit(names(PFMSet1[PFM_num]),"[_]"))[1]
    
    ### construct PWM1 from PFM1
    PWM1=matrix(0,ncol=ncol(PFM1),nrow=4)
    for (j in 1:4){
      for (i in 1:ncol(PWM1)){
        sum=PFM1[1,i]+PFM1[2,i]+PFM1[3,i]+PFM1[4,i]
        if (j==1 |j==4){
          PWM1[j,i]=log2((PFM2[j,i]+(1-GC_content)/2*pseudo_count)/(sum+pseudo_count)/((1-GC_content)/2))
        }
        else{PWM1[j,i]=log2((PFM2[j,i]+GC_content/2*pseudo_count)/(sum+pseudo_count)/(GC_content/2))}
      }
    }
    rownames(PWM1)=c("A","C","G","T")
    l_PWM1<-length(PWM1)/4
    
    ### calculate the background PWM score using PWM1 for the input background sequences for one strand 
    BackGroundPWMScore=PWMscoreStartingAt(PWM1,background_seq,1:(length(background_seq)-l_PWM1+1))
    
    ### calculate the background PWM score using PWM1 for the input background sequences for another strand 
    ### BackGroundPWMScore_rev=rev(PWMscoreStartingAt(PWM1,Reverse_complement,1:(length(background_seq)-l_PWM1+1)))
    
    ### assign PWM score to each position of the genome
    ### BackGroundPWMScore=rbind(BackGroundPWMScore,BackGroundPWMScore_rev)
    ### BackGroundPWMScore=apply(BackGroundPWMScore,2,max)
    
    ### estimate the residence time for PWM1
    log_tau_PWM1=rep(0,50)
   
    avg_exp_PWM1 = mean(exp((1/(as.numeric(lambda_reference)))*BackGroundPWMScore));
      
    tau_0 = (2*tR)/(avg_exp_PWM1*sl_obs^2)
    for (j in quantile_range){
      quant=quantile(BackGroundPWMScore,(1-10^(-j*0.1))) 
      log_tau_PWM1[j]=log(tau_0)+quant/as.numeric(lambda_reference)
    } 
   
    ### ref_residence_time saves the reference -log(residence time) in each binding site strength level
    ### calculated using PWM1 and known lambda
    ref_residence_time=log_tau_PWM1[quantile_range]
    
    ### If the -log(residence time) of some lambda choices for PWM2 deviate too much from the reference 
    ### obtained from PWM1 (we set a threshold of +-1.5), these lambda values are removed
    reference_range=rbind((ref_residence_time-1.5),(ref_residence_time+1.5))
    unsuitable=c()
    for( i in 1:nrow(residence_time_range)){
      for (j in 1:ncol(residence_time_range)){
        if(residence_time_range[i,j]<reference_range[1,j] | residence_time_range[i,j]>reference_range[2,j]){
          unsuitable=c(unsuitable,i)
          break
        }
      }
    }
    suitable=1:nrow(residence_time_range)
    suitable=suitable[-unsuitable]
    residence_time_range=residence_time_range[suitable,]
    
    ### calculate the mean square error between the reference residence time from PWM1 in each binding site strength
    ### level and the estimated residence time from PWM2 using each potential lambda value 
    Mse=c()
    for( i in 1:nrow(residence_time_range)){
    Mse=c(Mse,mse(exp(as.numeric(residence_time_range[i,])),exp(as.numeric(ref_residence_time))))
    }
    
    ### select the lambda value for PWM2 that minimize the mean square error
    lambda_chose=(which(Mse==min(Mse))[1]+suitable[1]-1+lambda_range_times_ten[1]-1)/10
    
    ### An alternative approach to find the optimized lambda is using the sum of absolute differences between 
    ### the reference and the estimated residence times for -log(residence time)
    ### Code is as follows:
    
    ############################################################################################################
    # sum_of_absolute_dif=c()
    # for( i in 1:nrow(residence_time_range)){
    # sum_of_absolute_dif=c(sum_of_absolute_dif,sum(abs(residence_time_range[i,]-as.numeric(ref_residence_time))))
    # }
    # lambda_chose=(which(sum_of_absolute_dif==min(sum_of_absolute_dif))[1]+suitable[1]-1+lambda_range_times_ten[1]-1)/10
    #############################################################################################################
    
    
    lambda_PWM2=rbind(lambda_PWM2,c(TF,lambda_chose))
    
    ### plot the heatmap of the -log(residence time) for PWM2. 
    ### If you don't want to see the heatmap, comment out the following 2 lines.
    reference_range=c(min(ref_residence_time)-0.2,max(ref_residence_time)+0.2)
    image2D(x=lambda_range_times_ten/10,y=quantile_range/10,z=residence_time_range_save,zlim=reference_range,col =  heat.colors(100),xlab = expression(lambda),ylab="-log(top quantile)",main=TF)
  }
  
  colnames(lambda_PWM2)=c("TF_name","lambda_for_PWM2")
  ### This function returns a matrix with TF name and suitable lambda values for PFMs listed in PFMSet2
  Convert_lambda_between_two_PWMs=lambda_PWM2
}
