#' Download all files from osf directory 
#' This function is intended to run inside dplyr do() verb. Downloads data locally
#' Created by Filip Dechtěrenko
#'
#' @param df row with directory on osf repository
#' @param local_data_pth Path where the data should be saved
#' @param should_overwrite whether the data should be overwritten. Default is True
#'
#' @return
#' @export
#'
#' @examples
#' local_data_pth <- file.path("data","Source")
#' create_local_structure(local_data_pth)
#' data_files <- 
#'   weissman_project %>% 
#'   osf_ls_files() %>% 
#'   filter(name == "Source") %>%
#'   osf_ls_files() 
#'   
#'   data_files %>% 
#'   group_by(name) %>% # for each experiment type
#'     do(download_files(.,local_data_pth))
download_files <- function(df, local_data_pth) {
  # we need to set correct class as the current version of osfr does not works with dplyr properly
  class(df) <- c("osf_tbl_file","osf_tbl", class(df)) 
  stopifnot(nrow(df) == 1)
  df %>% 
    osf_ls_files() %>% 
    rowwise() %>% 
    do(osf_retrieve_file(.$id) %>% 
         osf_download(path = file.path(local_data_pth, df$name, .$name)))
}

#' Create local file structure
#' As we have four experiments, we create individual directory for each of them
#' @param pth path, in which the file structure will be created
#'
#' @return nothing, just the side effects
#' @export
#'
#' @examples
#' local_data_pth <- file.path("data","Source")
#' create_local_structure(local_data_pth) 
#'
create_local_structure <- function(pth) {
  # simple function for test-and-create behavior
  verbose_create <- function(f) {
    if(!dir.exists(f)) {
      message(sprintf("Directory %s does not exist, creating...", f))
      dir.create(f, recursive = T)  
    }
  }
  
  # Create experiment level subfolder
  Exp1_dir <- file.path(pth,"Exp1")
  Exp2_dir <- file.path(pth,"Exp2")
  Exp3_dir <- file.path(pth,"Exp3")
  
  verbose_create(Exp1_dir)
  verbose_create(Exp2_dir)
  verbose_create(Exp3_dir)

}

#' Removes local data 
#' Just a simple wrapper around unlink function
#'
#' @param local_data_pth path to downloaded data
#'
#' @return
#' @export
#'
#' @examples
#' local_data_pth <- file.path("data","Source")
#' create_local_structure(local_data_pth)
#' remove_local_data(local_data_pth)
remove_local_data <- function(local_data_pth) {
  unlink(local_data_pth, recursive = T)
}

#' Function to not include a vector in another vector
#' 
#' Retrieved from https://stackoverflow.com/questions/5831794/opposite-of-in
#' 
`%ni%` <- Negate(`%in%`)

#' Function to caluclate Bayes factors
#' 
## The function is retrieved from 
## https://link.springer.com/article/10.3758/s13423-017-1266-z

Bf <- function(sd, obtained, dfdata, meanoftheory, sdtheory, dftheory, tail = 2)
  
{
  
  area <- 0
  
  normarea <- 0
  
  theta <- meanoftheory - 10 * sdtheory
  
  incr <- sdtheory/200
  
  for (A in -2000:2000){
    
    theta <- theta + incr
    
    dist_theta <- dt((theta-meanoftheory)/sdtheory, df=dftheory)
    
    if(identical(tail, 1)){
      
      if (theta <= 0){
        
        dist_theta <- 0
        
      } else {
        
        dist_theta <- dist_theta * 2
        
      }
      
    }
    
    height <- dist_theta * dt((obtained-theta)/sd, df = dfdata)
    
    area <- area + height * incr
    
    normarea <- normarea + dist_theta*incr
    
  }
  
  LikelihoodTheory <- area/normarea
  
  Likelihoodnull <- dt(obtained/sd, df = dfdata)
  
  BayesFactor <- LikelihoodTheory/Likelihoodnull
  
  BayesFactor
  
}