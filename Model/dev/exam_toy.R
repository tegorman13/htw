
source(here::here("Model/dev/funs_cm.r"))




## Toy scenario 
# input space ranges from X=1:7. True Function is Y=X. 
# Trained from X=3:5
# Assume perfect learning - initalize weight matrix diagonals to 1 for trained pairs. 

input_layer <- c(1,2,3,4,5,6,7)
output_layer <- c(1,2,3,4,5,6,7)

# generate weight matrix reflective of learning from X=3:5
# value of 1 for trained pairs, 0 else where
weight_mat <- matrix(0, nrow=7, ncol=7)
# only portion of diag that was trained set to 1, e.g. (3,3), (4,4), (5,5). 
weight_mat[3,3] <- 1
weight_mat[4,4] <- 1
weight_mat[5,5] <- 1

# generate training data
#train_vec <- c(3,4,5)
train_vec <- c(3)
train_outputs <- train_vec

n=5

td <- tibble(x=rep(train_vec, each=n), y=rep(train_outputs, each=n)) |> 
  mutate(trial=sample(seq(1,n*length(train_vec)), replace=FALSE)) |> 
  arrange(trial)

train_dat <- alm.sim(dat=td,
   c=.5, lr=.01, input_layer, output_layer, weight_mat)
train_dat$d

# test alt_exam on extrapolation items x=1 & x=2. 
map_dbl(c(1,2,4,5), ~ alt_exam(.x, c=1, input_layer, output_layer, train_dat$wm, train_vec, train_dat$d$almResp))





input_layer <- c(100, 350, 600, 800, 1000, 1200)
output_layer <- c(100, 350, 600, 800, 1000, 1200)

ds <- readRDS(here::here("data/e1_md_11-06-23.rds")) |> mutate(sbj=id) |> relocate(sbj,.after=id)


split_data <- split(ds, ds$id)
split_data <- ds %>% split(.$id) |> head(1)
c = .0008; lr=2.0; 
#split_data[[1]] |> filter(expMode2=="Train")

train_results <- alm.sim(split_data[[1]] |> filter(expMode2=="Train"), 
                         c, lr, input_layer, output_layer)

trainVec <- unique(train_results$d$x) |> sort()
trainOutputs <- unique(train_results$d$y) |> sort()
# generate predictions from exam_infer.response
map_dbl(c(100, 350, 600, 800, 1000, 1200), ~ exam_infer.response(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec, trainOutputs))

test_prediction <- map_dbl(input_layer, ~ exam_infer.response(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec, trainOutputs))
map_dbl(input_layer, ~ alm.responseOnly(.x, c, input_layer, output_layer, train_results$wm))
map_dbl(input_layer, ~ exam.response(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec))




split_data <- split(ds, ds$id)
split_data <- split_data[c(3)]
c = .0008; lr=2.0; 
train_results <- alm.sim(split_data[[1]] |> filter(expMode2=="Train"), 
                         c, lr, input_layer, output_layer)

trainVec <- c(0,unique(train_results$d$x)) |> sort()
trainOutputs <- unique(train_results$d$y) |> sort()
map_dbl(350, ~ exam_infer.response(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec, trainOutputs))
map_dbl(input_layer, ~ alm.responseOnly(.x, c, input_layer, output_layer, train_results$wm))
map_dbl(input_layer, ~ exam.response(.x, c, input_layer, output_layer, train_results$wm,  trainVec=trainVec))

