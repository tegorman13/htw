




generate.data <- function(x, type = "linear", noise = NA) {
  y =switch(type,
         linear = round(2.2*x + 30,0),
         exponential = round(200*(1-exp(-x/25)),0),
         sinusoidal = sin(2 * pi * x),
         quadratic = round(210 - ((x-50)^2)/12,0),
         stop("type must be linear, exponential, quadratic, or sinusoidal"))
  
  if(!is.na(noise)) {
    y <- y + round(rnorm(length(y), 0, noise),2)
  }
  data.frame(x, y, type)
}





deLosh_data <- list(
  human_data_linear = data.frame(
    x = c(0.7,2.6,4.5,6.6,8.5,11.3,12.7,14.3,
          17.6,18.3,20.9,22.1,24.7,27.7,28,32.9,34.6,37.6,40,
          43.3,44.9,47.5,52.2,55,59.8,65.6,67.3,
          72,72.9,76.2,76.9,81.2,81.6,84,85.9,87.5,89.9,92.7,
          93.7,94.8,99.1,99.5),
    y = c(11.7,17.7,25.5,27.3,29.8,38.3,41.1,
          47.1,56.7,61,67.3,69.5,77.6,83.6,92.5,100.6,106.3,
          112.7,118.7,125.1,128.6,135,144.2,150.9,160.5,
          170.3,176.4,180.6,185.9,189.4,194.8,196.9,199.7,
          203.6,207.5,212.8,219.2,219.1,220.5,228.7,230.1,234.4)),
  
  human_data_exp = data.frame(
    x = c(0.7,2.3,5.5,6.9,8.5,10.6,15.9,16.2,
          13.6,18.7,21.7,22.4,24.9,27,29.1,33.3,34.6,38.1,40.2,
          43.4,44.8,47.8,49.9,52.2,55,58,60.5,63,66.7,71.4,
          72.7,76.7,79.7,81.1,85.7,87.1,90.1,93.1,94.5,96.3,98.2,
          99.3),
    y = c(62.1,73.9,82.5,87,90.8,93.8,103.1,
          105.9,108.7,114.9,114.1,119,124.8,127.6,137.6,145.5,
          149.3,155.1,159.3,163.4,165.8,168.5,170.9,173.3,175.4,
          177.1,182.6,184.3,186.3,187,188.4,192.8,194.2,197.3,200,
          200.4,206.9,206.9,207.6,209.3,212.7,211.7)),
  
  human_data_quad = data.frame(
    x = c(1.1,2.5,5.1,7.2,10.7,13.2,15.1,18.6,
          17,21.1,22.7,25.3,27.4,28.8,32.8,35.3,38.6,43,
          45,48,50.8,53.3,55.4,57.7,60,65.3,68.5,71.5,
          73.5,75.4,77.2,81.1,79.7,83.1,85.4,88,89.3,90.9,93.2,
          95.5,97.8,99.6),
    y = c(70.5,82.3,99.9,101,108.3,116.6,122.8,
          136.3,137.7,140.1,143.9,154.7,162.3,166.1,181.6,191.7,
          196.2,203.1,205.2,206.6,208.7,206.6,206.6,203.2,
          198,193.2,184.2,177.6,176.3,171.1,166.9,160.4,
          155.5,148.3,147.3,142.8,133.8,132.7,122.4,119.3,118.6,
          105.5))
  
)




envTypes <- c("linear", "exponential", "quadratic")

trainingBlocks <- list(
   low = c(30.5, 36.0, 41.0, 46.5, 53.5, 59.0, 64.0, 69.5),
  
  med = c(
    30.0, 31.5, 33.0, 34.5, 36.5, 38.5, 41.0, 43.5, 46.0,
    48.5, 51.5, 54.0, 56.5, 59.0, 61.5, 63.5, 65.5, 67.0, 68.5, 70.0
  ),
  high = c(
    30.0, 30.5, 31.0, 32.0, 33.0, 33.5, 34.5, 35.5,
    36.5, 37.0, 38.0, 38.5, 39.5, 40.5, 41.5, 42.0, 43.0,
    43.5, 44.5, 45.5, 46.5, 47.0, 48.0, 48.5, 49.0, 51.0, 51.5, 52.0,
    53.0, 53.5, 54.5, 55.5, 56.5, 57.0, 58.0, 58.5, 59.5, 60.5, 61.5, 
    62.0, 63.0,63.5, 64.5, 65.5, 66.5, 67.0, 68.0, 69.0, 69.5, 70.0
  )
)

# trainingBlocks <- list(
#   low = c(30.5, 36.0, 41.0, 46.5, 53.5, 59.0, 64.0, 69.5),
#   med = seq(30.0, 70.0, length.out=20),
#   high = seq(30.0, 70.0, length.out=50)
# )


# 
# lowTrain <- map_dfr(envTypes, ~ generate.data(rep(lowDensityTrainBlock,25), type = .x)) %>% group_by(type) %>% mutate(block = rep(1:25, each = 8),trial=seq(1,200))
# medTrain <- map_dfr(envTypes, ~ generate.data(rep(medDensityTrainBlock,10), type = .x)) %>% group_by(type) %>% mutate(block = rep(1:10, each = 20),trial=seq(1,200))
# highTrain <- map_dfr(envTypes, ~ generate.data(rep(highDensityTrainBlock,4), type = .x)) %>% group_by(type) %>% mutate(block = rep(1:4, each = 50),trial=seq(1,200))


# short definition version

# lowDensityTrainBlock <- c(30.5, 36.0, 41.0, 46.5, 53.5, 59.0, 64.0, 69.5)
# medDensityTrainBlock <- seq(30.0, 70.0, length.out=20)
# highDensityTrainBlock <- seq(30.0, 70.0, length.out=50)



# generate.data <- function(x, type = "linear", noise = NA) {
#   if (type == "linear") {
#     y <- round(2.2*x + 30,0)
#   }
#   else if (type == "exponential") {
#     y <- round(200*(1-exp(-x/25)),0)
#   }
#   else if (type == "sinusoidal") {
#     y <- sin(2 * pi * x) 
#   }
#   else if (type == "quadratic") {
#     y <- round(210 - ((x-50)^2)/12,0)
#   }
#   else {
#     stop("type must be linear, exponential, quadratic, or sinusoidal")
#   }
#   # if noise is specified, add noise to the y values
#   if(!is.na(noise)) {
#     y <- y + round(rnorm(length(y), 0, noise),2)
#   }
#   data.frame(x, y,type)
#   
# }