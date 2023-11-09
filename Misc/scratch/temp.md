Background info on the models I'm working with. 
### ALM Model
###### Input Activation
$$
a_i(X)=\exp \left|-\gamma \cdot\left[X-X_i\right]^2\right|
$$
###### Output activation
$$
o_j(X)=\Sigma_{i=1, M} w_{j i} \cdot a_i(X) 
$$
###### Output Probability
$$
P\left[Y_j \mid X\right]=o_j(X) / \Sigma_{k=1, L} o_k(X) 
$$
###### Mean Response
$$
m(X)=\Sigma_{j=1, L} Y_j \cdot P\left[Y_j \mid X\right] 
$$

###### Feedback Signal
$$
f_j(Z)=e^{-c\cdot(Z-Y_j)^2}
$$
###### Weight Updates
$$
w_{ji}(t+1)=w_{ji}(t)+\alpha \cdot {f_i(Z(t))-O_j(X(t))} \cdot a_i(X(t))
$$
:::
### Exam Generalization
###### Input node actvation
$$
P[X_i|X] = \frac{a_i(X)}{\\sum_{k=1}^Ma_k(X)}
$$
###### Slope Computation
$$
E[Y|X_i]=m(X_i) + \bigg[\frac{m(X_{i+1})-m(X_{i-1})}{X_{i+1} - X_{i-1}} \bigg]\cdot[X-X_i]
$$


# Kwantes 2006

### Kwantes version of EXAM

Note that the calculation of $\boldsymbol{\Delta}$ requires the existence of a previous trial. That is, there is no representation for $\Delta$ for the first trial when $t=1$. Also, $\Delta$ is undefined when the values of **X** for trials $t$ and $t$-$1$ are the same. In either case, the field that corresponds to $\Delta$ is set to "null" and not considered during retrieval.

After training, memory contains as many traces as there were trials. The traces are lined up in such a way that the contents of memory can be viewed as a matrix. It is worth noting at this point, that we do not propose that trainees store an analogous representation of numbers in episodic memory. Instead, we favour a representation scheme like the one suggested by Hintzman (1984) in which items are represented as a vector of features. However, because we have complete control over the properties of any vectors representing number information in memory, we can bypass some of the representation issues and deal directly with expected values. Hence, much of the model we describe does not actually use vector representations even though, as a whole, we subscribe to the basic idea.

Once trained, the model is tested by letting a vector representing the predictor resonate with, or activate, the contents of memory. The extent to which a memory trace is activated by the probe is a function of the similarity of the two. When probe and trace contain more than one predictor, the similarity is measured as the average similarity of all the predictors. Finally, as mentioned above, we assume that proximal numbers share similar representations, and that the similarity between numbers decreases as the distance between them increases. The model has a parameter, $\varphi$, which reflects the pre-existing similarity between adjacent numbers (a parameter we set to $0.92$). The activation, $A$, of a single memory trace, $T_{i}$, by the probe, $P$ is calculated as the similarity between the two magnitudes.

$$A_{i}=\varphi^{|X_{P}-X_{T_{i}}|}$$ (3)

Once activated, the system selects an instance from memory as the best match to the probe. The probability of selecting one memory trace is a function of the trace's activation.

Specifically, the probability of selecting trace $i$ is equal to the activation of a trace divided by the sum of the activations of the $M$ traces in memory. More formally,

$$P(T_{i})=\frac{A_{i}}{\sum_{j=1}^{M}A_{j}}$$ (4) Instead of selecting one trace from memory, and operating upon its value of the **X'**, **Y'** and $\boldsymbol{\triangle}$', we opted to calculate expected values of the retrieved information using the following formulas,

$$E(X^{\prime})=\sum_{j=1}^{M}X_{j}\times P(T_{j})$$ (5a) 
$$E(Y^{\prime})=\sum_{j=1}^{M}Y_{j}\times P(T_{j})$$ (5b) 
$$E(\Delta^{\prime})=\sum_{j=1}^{M}\Delta_{j}\times P(T_{j})$$ (5c)

It is at this point that the model diverges from a strict memory-based account of function learning. If the predictor value falls outside of the domain of the training set, the best match from memory is most likely to be one of the items from boundary of the training set. Clearly, like EXAM, the system needs some way to adjust **Y'** when the predictor values fall outside the range shown during training. Again, like EXAM, we allow **Y'** to be adjusted to a degree that reflects the similarity/disparity between **X'** and **X** Instead of calculating slope estimates as part of the output stage, the adjustment on **Y'** (now denoted **Y'**(_new_) ) is done in such a way that its new value satisfies a constraint imposed by each of the $\boldsymbol{\triangle}$'s retrieved from memory.

How the system settles on a value of Y depends on the task that subjects are asked to perform. If asked to report an estimate of Y, we assume that the trainee searches for a value of Y starting at the closest matching value it retrieves from memory. To save time, we can calculate the new value of **Y'** directly by rearranging the terms of an equation almost identical to Equation 2 above. The equation is formally equivalent to Equation 1, taken from the Delosh et al (1997), used to adjust the retrieved value of Y.

$$Y^{\prime}(new)=Y^{\prime}(retrieved)-\left(\Delta^{\prime}{}_{i}\times \left[X^{\prime}{}_{i}-X_{i}\right]\right)$$ (6)

Delosh et al's (1997) participants did not report their estimates of Y. Instead, they filled a horizontal bar containing numerically labelled tick marks. They filled the bar by starting at zero and moving the fill up to their estimate of Y. We treat the way that trainees perform the estimation task as a constraint on the search process. In other words, we propose that, trainees move up on a mental number line to their estimate of Y(_new_). As they search the number line, they evaluate the goodness of the estimate. At each evaluation (arbitrarily set to be done at every increase of 1), the system calculates the difference between the estimated and retrieved value of Y relative to the difference between the cue and retrieved value of X (i.e., the slope of the line between (**X'**,**Y'**) and (**X**, **Y'**(_new_).** See Figure 2). The goodness of the estimate is the discrepancy between the calculation and retrieved value of $\boldsymbol{\triangle}$. The system stops changing **Y'**(new) when the difference reaches a minimum criterion (a parameter we set to 0.1). In formal terms, the discrepancy (or fit) is calculated as, The more lax the criterion is, the farther the slope value will be from the retrieved value of $\Delta$ when it stops searching. Whether **Y'(new)** overestimates or underestimate the correct value of **Y** depends on the direction in which trainees move their estimate. When trainees start their estimate at zero and move up, to the extent that the criterion is not set to zero, **Y'(new)** will underestimate the correct value of **Y**. The opposite will be true if trainees start at a maximum value of **Y** on their mental number line and adjust it down to their estimate of **Y'(new)**--** trainees should then overestimate the correct value of Y.

1) Given a probe of X, the best match in memory is (X', Y', $\Delta$)
2) The system starts Y at zero and searches up to find a value of Y(new) where the slope of (X,Y(new)) is similar to the slope associated with (X’,Y’) in memory


Our interpretation for why trainees tended to underestimate Y is very different from the explanation offered by DeLosh et al (1997). Recall that their explanation placed responsibility for the underestimate on a poor representation of Y in memory. Our interpretation does not place responsibility on the quality of the representations. Instead, we place it at the output stage and consider the underestimation an artefact of way in which trainees search for their estimate of the criterion.
