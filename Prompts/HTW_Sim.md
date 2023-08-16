

Please only reply using R code. Please concisely implement the experiment described with 30 lines or less.  -Sorry, no line breaks. Please refrain from leaving comment-outs.


Experiment design: Two groups of subjects, constant and varied. Constant trains with 270 trials from position (800) , varied trains with 90 trials each from 3 Positions (800, 1000, 1200). Then, both groups are tested from 6 positions (200, 400, 600, 800, 1000, 1200). All participants improve during training, but the constant group achieves the best performance overall. Participants generate responses that are fairly noisy. Write R code to simulate this experiment. 


Here's the description of an experiment: "90 training trials split evenly divided between velocity bands. Varied training with 3 velocity bands and Constant training with 1 band. No-feedback testing from 3 novel extrapolation bands. 15 trials each. No-feedback testing from the 3 bands used during the training phase (2 of which were novel for the constant group). 9 trials each. Feedback testing for each of the 3 extrapolation bands. 10 trials each. Varied Training 800-1000 1000-1200 1200-1400 Testing - No Feedback 100-300 350-550 600-800 Constant Training 800-1000 Test From Train 800-1000 1000-1200 1200-1400 Testing -Feedback 100-300 350-550 600-800"



_________
Write R code to simulate this experiment. 20 varied participants, and 20 constant participants. In the learning phase, performance changes as an exponential decay function, with learning rate, asymptotic performance, and starting performance parameters. Performance for constant and varied participants are drawn from their own gaussian distributions.  For each participant, training and testing performance should be correlated, the strength of the correlation is set by a parameter. There is also an effect of training band, such that higher test bands are more difficult. Use the data.table, tidyverse, and purrr packages.

