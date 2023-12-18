

pacman::p_load(tidyverse,here,patchwork, conflicted)
walk(c("fun_alm","fun_model", "Display_Functions"), ~ source(here::here(paste0("Functions/", .x, ".R"))))

act1 <- compose(~activation_function(.x,input_layer,.y))
act1(1000,1)

act1 <- compose(\(x,y) activation_function(x,input_layer,y))
k=act1(1000,1)
k
length(k)
toString(k)
length(toString(k))


act2 <- compose(\(x) toString(x), \(x,y) activation_function(x,input_layer,y))
act2(1000,1)


act2 <- compose(\(x,y) activation_function(x,input_layer,y), \(x) toString(x), .dir = "forward")
act2(1000,1)

