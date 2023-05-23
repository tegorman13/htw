

gen_trials <- function(blocks=3,numinputs=3,shuffle=FALSE) {
  if (shuffle) {
    examples <- as.vector(apply(replicate(blocks,seq(1, numinputs)), 2,
                                sample, numinputs))
  } else{
    examples <- as.vector(replicate(blocks, seq(1, numinputs)))
  }
}
