

#pacman::p_load(parallel, future, furrr)

tg_m1 <- '10.147.17.202'
macpro <- '10.147.17.83'

make_spec <- function() {

tg_m1 <- '10.147.17.202'
m1l1 <- '10.147.17.245'
m1l2 <- '10.147.17.69'
m1l3 <- '10.147.17.122'
m1l4 <- '10.147.17.46'
m1l5 <- '10.147.17.192'
m1l6 <- '10.147.17.66'
m1l7 <- '10.147.17.139'
m1l8 <- '10.147.17.8'
m1r1 <- '10.147.17.45'
m1r2 <- '10.147.17.211'
tg15 <- '10.147.17.222'
macpro <- '10.147.17.83'

# alternative m1l8: 10.20.40.99
#ips <- c(m1l2,tg_m1)


# machineAddresses <- list(
#   #list(host=tg_m1,user='thomasgorman', ncore=7),
#   list(host=macpro,user='thomasgorman', ncore=1),
#   list(host=m1l2,user='thomasgorman', ncore=8),
#   list(host=m1l4,user='thomasgorman', ncore=8),
#   list(host=m1l5,user='thomasgorman', ncore=8),
#    list(host=m1l6,user='thomasgorman', ncore=8),
#    list(host=m1l7,user='thomasgorman', ncore=8),
#    list(host=m1r1,user='thomasgorman', ncore=8),
#   list(host=m1r2,user='thomasgorman', ncore=8)
# )

machineAddresses <- list(
  list(host=tg_m1,user='thomasgorman', ncore=8),
  list(host=m1l1,user='thomasgorman', ncore=8),
  list(host=m1l3,user='thomasgorman', ncore=8),
  list(host=m1l8,user='thomasgorman', ncore=8)
)


spec <- lapply(machineAddresses,
               function(machine) {
                 rep(list(list(host=machine$host,
                               user=machine$user)),
                     machine$ncore)
               })
spec <- unlist(spec,recursive=FALSE)
return(spec)

}

# pc <- parallel::makeCluster(type='PSOCK', master=tg_m1, spec)
# plan(cluster, workers = pc)


#makePSOCKcluster(ips, verbose=TRUE)






#fr = future_map(1:390, ~ rnorm(30000) * .x)








# killall R
# parallel::stopCluster(pc)



# if(!is.null(parallelCluster)) {
#   parallel::stopCluster(parallelCluster)
#   parallelCluster <- c()
# }