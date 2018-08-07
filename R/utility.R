#' Load an RData file
#'
#' \code{load_rdata} loads an RData file into a specified variable
#'
#' @param  file_name a character string giving the name of the file to load
#' @return The objects in the RData file. If there is more than one object, a list will be returned.
#'
#' @examples
#' \dontrun{
#' x = pi # to ensure there is some data
#' x_list = list(pi = pi, pi_sq = pi^2)
#' save(x, file= "pi.RData")
#' save(x_list, file= "pi_list.RData")
#' rm(x)
#' rm(x_list)
#' x = load_rdata('pi.RData')
#' x_list = load_rdata('pi_list.RData')
#' }
#' @keywords internal
load_rdata = function(file_name){
     load(file_name)
     obj_names = ls()
     obj_names = obj_names[obj_names != 'file_name']
     if(length(obj_names) > 1){
          #return a list of objects
          dat = sapply(obj_names, function(x)get(x), simplify=FALSE, USE.NAMES=TRUE)
     }else{
          #return a single object
          dat = get(obj_names)
     }
     return(dat)
}

#' Compute Bayes factor from powder objects
#'
#' \code{bayes.factor} computes the Bayes factor from the powder objects
#'
#' @param m1 a Powder.Hierarchical or Powder.Individual object
#' @param m2 a Powder.Hierarchical or Powder.Individual object
#' @param log When TRUE, the natural log of the Bayes factor will be returned.
#' @return The natural log of the Bayes factor in terms of m1 / m2.
#' @export
bayes.factor = function(m1,m2,log=TRUE){
     bf = data.frame(Method = c('TI','Log Steppingstone','Harmonic Mean'),
                     Log.Bayes.factor= c(m1$Value[m1$Method == 'TI'] - m2$Value[m1$Method == 'TI'],
                                 m1$Value[m1$Method == 'Log Steppingstone'] - m2$Value[m1$Method == 'Log Steppingstone'],
                                 m1$Value[m1$Method == 'Harmonic Mean'] - m2$Value[m1$Method == 'Harmonic Mean']))
     if (log==FALSE) {
          bf$Log.Bayes.factor = exp(bf$Log.Bayes.factor)
          colnames(bf) = c('Method','Bayes factor')
     }
     return(bf)
}
