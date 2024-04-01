get_seq_from_granges <- function(genome, gr_obj){
  for (i in 1:length(gr_obj)){
    gr <- gr_obj[i]
    if (i == 1){
      res <- genome[gr]
    }else{
      res <- append(res, genome[gr])
    }
  }
  names(res) <- names(gr_obj)
  return(res)
}