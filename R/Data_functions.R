##################Data#########################
# The data is publicly available on https://github.com/m-a-robinson/sample-size
####Muscle Force####

#Gomes, A.A., Ackermann, M., Ferreira, J.P., Orselli, M.I.V. and Sacco, I.C.,
#2017. Muscle force distribution of the lower limbs during walking in diabetic 
#individuals with and without polyneuropathy. Journal of neuroengineering and 
#rehabilitation, 14, pp.1-13.

#They completed to be 101 points (0,1,...,100) by linear interpolation.

#The force time series for ankle extensor muscle (Soleus) for control and diabetic

MF_data <- function(type){
  if (is.null(type)) {
    return(NULL)  # Return NULL or an empty dataset
  }
  data <- read.delim("data/Gomes2017_MF.txt", header=FALSE)
  Domain <- c(0,100)
  
  if (type=="control") {
    control <- completed_data(x = data[,1], y = data[,2],
                                  defined_domain = Domain)
    data_out <- control
    
  }else if(type=="diabetic"){
    diabetic <- completed_data(x = data[,3], y = data[,4],
                                  defined_domain = Domain)
    data_out <- diabetic
    
  }else if(type=="both"){
    
    control <- completed_data(x = data[,1], y = data[,2],
                                  defined_domain = Domain)
    diabetic <- completed_data(x = data[,3], y = data[,4],
                                  defined_domain = Domain)
    data_out <- cbind(control, diabetic)
    colnames(data_out) <- c("Control","Diabetic")
  }else{
        print("Invalid type");     return(NULL)  # Handle invalid type
  }
  return(data_out)
}

