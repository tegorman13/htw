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




###### Input node actvation
$$
P[X_i|X] = \frac{a_i(X)}{\\sum_{k=1}^Ma_k(X)}
$$

###### Slope Computation
$$
E[Y|X_i]=m(X_i) + \bigg[\frac{m(X_{i+1})-m(X_{i-1})}{X_{i+1} - X_{i-1}} \bigg]\cdot[X-X_i]
$$