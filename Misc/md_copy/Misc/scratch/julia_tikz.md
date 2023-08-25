---
title: "rbf.qmd"
jupyter: julia-1.9
---

```{r}
pacman::p_load(tidyverse,RSNNS,tikzDevice,knitr,exams)
```


tinytex::reinstall_tinytex(repository = "illinois")



```{r}
#| echo: false
# Necessary for using dvisvgm on macOS
# See https://www.andrewheiss.com/blog/2021/08/27/tikz-knitr-html-svg-fun/
#Sys.setenv(LIBGS = "/usr/local/share/ghostscript/9.53.3/lib/libgs.dylib.9.53")
font_opts <- list(dvisvgm.opts = "--font-format=woff")

library("tinytex")

# options(
#   tinytex.engine = "xelatex",
#   tikzDefaultEngine = "xetex",
#   tikzXelatexPackages = c(
#     "\\usepackage[fontset=fandol]{ctex}",
#     "\\usepackage{amsfonts,mathrsfs,amssymb,pgfplots}\n"
#   )
# )
```



```{julia}

using DifferentialEquations

# parameter values from the paper
π = 2133
μ = 1/30
c = 1.7
βᵤ = 0.1
βₜ = 0.025
σ = 0.5
g = 0.05
νᵤ = 1/12
νₜ = 1/27
function blower₀₀(du,u,P,t)
    π, μ,c, βᵤ, βₜ, σ, g, νᵤ, νₜ, i = P
    X, Yᵤ, Yₜ = u
    λ = (βᵤ*Yᵤ + βₜ*Yₜ)/(X + Yᵤ + Yₜ)

    du[1] = π - c*(1 + i)λ*X - μ*X
    du[2] = c*(1 + i)λ*X + g*Yₜ - σ*Yᵤ - μ*Yᵤ - νᵤ*Yᵤ
    du[3] = σ*Yᵤ - g*Yₜ - μ*Yₜ - νₜ*Yₜ
end


using QuadGK

function deaths(sol,t, νᵤ = 1/12, νₜ = 1/27)
    AUC, err = quadgk(sol, 0, t)
    AUC[2]νᵤ +AUC[3]νₜ
end


function get_z(i, tmax = 10, tmin = 1, step = 1)
    P = [π, μ,c, βᵤ, βₜ, σ, g, νᵤ, νₜ, i]
    prob = ODEProblem(blower₀₀, 
        [0.7,0.3,0.0] .*40000,  
        (0.0,10),
        P)  
    sol = solve(prob)
    ARTdeaths = [deaths(sol,t) for t in tmin:step:tmax]
    P = [π, μ,c, βᵤ, βₜ, 0, g, νᵤ, νₜ, 0]
    prob = ODEProblem(blower₀₀, 
        [0.7,0.3,0.0] .*40000, 
        (0.0,10),
        P)  
    sol = solve(prob)
    NULLdeaths = [deaths(sol,t) for t in tmin:step:tmax]
    z =  ARTdeaths ./ NULLdeaths 
    vcat(0,(1 .- z)*100) # zero at the begining, for plotting purposes
end

using GLMakie
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "Years", 
          ylabel = "Percent of AIDS deaths averted",
          xticks = 2:2:10, xgridvisible = false, ygridvisible = false)
# Create a slider widget
i_slider = Slider(fig[2,1],range = 0:0.01:3.0)
# Create a vector 'z' which is the return value from calling
# get_z() on the slider value
z = lift(i_slider.value) do i
    get_z(i)
end
lines!(0:10,z, color = :red, linewidth = 5)
limits!(ax, 0, 10, 0, 50)
fig
```

## pure latex

```{latex}

\usepackage{skull}
\usepackage{tikz}
\usetikzlibrary{positioning}
\begin{tikzpicture}
  \node (n)[thick, circle, align = center, draw = black] at (0,0) { Number of \\ mice \\ $n(t)$};
  \node[left = of n] (neighbours){};
  \node[below right = of n] (d) {$\skull$};
  \draw[thick, ->] (neighbours) -- node[above]{$m$} ++ (n);
  \draw[thick, ->] (n) -- node[sloped,above]{$d \, n(t)$} ++ (d);
  \draw [thick, ->] (n) edge [out=60,in=30,looseness=4] node[below right,label={[xshift=0.4cm, yshift=-0.9cm]$b \, n'(t)$}] {} (n);
\end{tikzpicture}

```


```{julia}
#| eval: false 
using TikzPictures
TikzPicture(L"""
  \def\Radius{2cm}
  \draw [ultra thick] (0cm,0cm) circle[radius=\Radius];
""")
```

```{julia}
using WGLMakie
s = Array(sol)'                
fig = Figure()
ax1 = Axis(fig[1,1],xticks = 0:30:120, 
           xgridvisible = false, ygridvisible = false,
           xlabel = "days", ylabel = "lymphocytes")
ax2 = Axis(fig[1,1],yaxisposition = :right, 
           xgridvisible = false, ygridvisible = false,
           ylabel = "virions")
lines!(ax1,(1000(1 - τ) .+ s[:,1] .+ s[:,2] .+ s[:,3])[1:end], 
        color = :black)
lines!(ax2,log10.(s[:,4]), color = :black)                     
ax2.yticks = (log10.([0.1,1,10,100,1000,10000]),              
              string.([0.1,1,10,100,1000,10000]))            
ylims!(ax2,(-1,log10(10000)))
ylims!(ax1,0.0,1200)
hidexdecorations!(ax2)
text!(ax1, "CD4 lymphocytes", position = (45,825))
text!(ax2, "Cell-free virus", position = (40,1))
fig


```


```{tikz engine.opts=font_opts}
#| eval: false


\begin{tikzpicture}[scale=1] 
  \begin{axis}[ 
      axis on top = true, 
      axis x line = bottom, 
      axis y line = left, 
      xlabel = $$, 
      ylabel = $$, 
      ymin = , 
      ymax =  
  ] 
  \addplot[ 
      domain=-10:10, 
      samples=100 
  ] {(x>=0)*x}; 
  \end{axis} 
\end{tikzpicture} 

```



```{tikz engine.opts=font_opts}
#| echo: true
#| out-width: 60%
#| caption: "dag"
#| label: fig-line-plot
#| fig-cap: "A caual graph "

\usetikzlibrary{positioning}
\usetikzlibrary{shapes.geometric}
\usetikzlibrary{arrows}
\usetikzlibrary{decorations}
\tikzstyle{Arrow} = [->, thin, preaction = {decorate}]
\tikzset{>=latex}

\begin{tikzpicture}[{every node/.append style}=draw]
  \node [ellipse, draw=white] (Age) at (0, 0) {$A$};
  \node [rectangle, draw=white] (Marriage) at (2, 0) {$M$};
  \node [rectangle, draw=white] (Happiness) at (4, 0) {$H$};
  \draw [-latex, draw=black] (Age) to (Marriage);
  \draw [-latex, bend left] (Age) to (Happiness);
\end{tikzpicture}
```


## Diagrams

```{tikz}
%| label: tikz-ex
%| fig-cap: "A simple tikz example"
%| out-width: "40%"
\usetikzlibrary{arrows}
\begin{tikzpicture}[node distance=2cm, auto,>=latex', thick, scale = 0.5]
\node (P) {$P$};
\node (B) [right of=P] {$B$};
\node (A) [below of=P] {$A$};
\node (C) [below of=B] {$C$};
\node (P1) [node distance=1.4cm, left of=P, above of=P] {$\hat{P}$};
\draw[->] (P) to node {$f$} (B);
\draw[->] (P) to node [swap] {$g$} (A);
\draw[->] (A) to node [swap] {$f$} (C);
\draw[->] (B) to node {$g$} (C);
\draw[->, bend right] (P1) to node [swap] {$\hat{g}$} (A);
\draw[->, bend left] (P1) to node {$\hat{f}$} (B);
\draw[->, dashed] (P1) to node {$k$} (P);
\end{tikzpicture}
```







```{r fig.width=11, fig.height=9}
#| eval: false 
n_iter=500
n_nodes <- c(c(1,3,6,9),seq(10,300,10))
train <- tibble(input=as.matrix(seq(0,100,1)),output=as.matrix(round(2.2*input + 30,0)))

fits <- map(n_nodes,~rbf(train$input,train$output,size=.x,maxit=n_iter))



crossing(n_nodes,train) %>% 
  cbind(.,fit=unlist(map(fits,fitted))) %T>% 
  {print(ggplot(d,aes(x=input,y=output))+geom_line()+
  geom_line(aes(x=input,y=fit),col="green") + facet_wrap(~n_nodes,scales="free")) } %>% 
  {. ->> d}
  
tibble(n_nodes=as.factor(rep(n_nodes,each=n_iter))) %>% 
  cbind(.,error=unlist(purrr::map(fits, ~.x$IterativeFitError))) %>%
  mutate(it=seq(1,n()),.by=n_nodes) %T>% 
  {print(ggplot(.,aes(it,error,col=n_nodes))+geom_line()+facet_wrap(~n_nodes,scales="free")) } %>% 
  {. ->> errVec}


preds <- map(fits,~predict(.x,as.matrix(seq(0,150,1))))

tibble(n_nodes=n_nodes,fit=map(fits,~predict(.x,as.matrix(seq(0,150,1)))))

tibble(n_nodes=rep(n_nodes,each=150),input=as.matrix(seq(0,150,1)),output=as.matrix(round(2.2*input + 30,0)))

```





```{r}
#| echo: false
#| eval: false 

# d <- crossing(n_nodes,train) %>% 
#   cbind(.,fit=unlist(map(fits,fitted)))

# purrr::map(fits, ~.x$IterativeFitError)
# unlist(purrr::map(fits, ~.x$IterativeFitError))

inputs <- as.matrix(seq(0,100,1))
outputs <- as.matrix(round(2.2*inputs + 30,0))

model <- rbf(inputs, outputs, size=10, maxit=5000,
                     initFuncParams=c(0, 1, 0, 0.01, 0.01),
                     learnFuncParams=c(1e-8, 0, 1e-8, 0.1, 0.8), linOut=TRUE)

ggplot(d,aes(x=input,y=output))+geom_line()+
  geom_line(aes(x=input,y=fit),col="green") + facet_wrap(~n_nodes,scales="free")


inputs <- as.matrix(seq(0,10,0.1))
outputs <- as.matrix(sin(inputs) + runif(inputs*0.2))
outputs <- normalizeData(outputs, "0_1")

model <- rbf(inputs, outputs, size=40, maxit=1000,
                     initFuncParams=c(0, 1, 0, 0.01, 0.01),
                     learnFuncParams=c(1e-8, 0, 1e-8, 0.1, 0.8), linOut=TRUE)

par(mfrow=c(2,1))
plotIterativeError(model)
plot(inputs, outputs)
lines(inputs, fitted(model), col="green")
```

```{r}




```

