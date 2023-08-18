---
title: OJS ALM
date-modified: last-modified
categories: [Simulation,ALM,OJS,Interactive]
#page-layout: full
format:
  html:
    grid:
      sidebar-width: 350px
      body-width: 1300px
      margin-width: 100px
      gutter-width: 1.5rem
    self-contained: true
    toc: false 
execute: 
   warning: false
code-fold: true
code-tools: true
---

```{r}
pacman::p_load(tidyverse)
d <- tibble(x=1:20,y=x^2)
ojs_define(d = d)

inputNodes = seq(1,7,1)  # 
outputNodes = seq(50,1600,50)
#wm=matrix(rnorm(length(inputNodes)*length(outputNodes),5,2),nrow=length(outputNodes),ncol=length(inputNodes))
#ojs_define(iN = inputNodes, outputNodes = outputNodes, wm = wm)
```

```{ojs}
d3 = require("d3@7")
math = require('mathjs')
// let inputNodes = Array.from({ length: 7 }, (_, i) => i + 1);
// let outputNodes = Array.from({ length: 32 }, (_, i) => (i + 1) * 50);
// let wm = Array.from({ length: outputNodes.length }, () =>
//   Array.from({ length: inputNodes.length }, () => 0.0)
// );

function inputActivation(xTarget, c) {
  console.log(inputNodes)
  return inputNodes.map((inputNode) =>
    Math.exp(-1 * c * Math.pow(xTarget - inputNode, 2))
  );
}


function outputActivation(xTarget, weights, c) {
  const inputAct = inputActivation(xTarget, c);
  return math.multiply(weights, inputAct);
}

function meanPrediction(xTarget, weights, c) {
  const outputAct = outputActivation(xTarget, weights, c);
  const probability = math.divide(outputAct, math.sum(outputAct));
  return math.multiply(outputNodes, probability);
}


function updateWeights(xNew, yNew, weights, c, lr) {
  const yFeedbackActivation = outputNodes.map(
    (outputNode) => Math.exp(-1 * c * Math.pow(yNew - outputNode, 2))
  );
 //console.log(yFeedbackActivation)
  const xFeedbackActivation = outputActivation(xNew, weights, c);
  const inputAct = inputActivation(xNew, c);
  const inputActReshaped = math.reshape(inputAct, [inputAct.length, 1]);
  const yFeedbackActivationReshaped = math.reshape(yFeedbackActivation, [yFeedbackActivation.length, 1]);
  const xFeedbackActivationReshaped = math.reshape(xFeedbackActivation, [xFeedbackActivation.length, 1]);
  const error = math.reshape(math.subtract(yFeedbackActivationReshaped, xFeedbackActivationReshaped), [yFeedbackActivation.length, 1]);
  // console.log(math.size(math.transpose(inputActReshaped)))
  const weightUpdate = math.multiply(error, math.transpose(inputActReshaped));
  //console.log(weightUpdate)
  const raw_Weights = math.add(weights, math.multiply(lr, weightUpdate));

  const new_Weights = raw_Weights
  //return JSON.parse(result);
  return(new_Weights)
}

function randomNormal(mean, sd) {
  let u = 0,
    v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  const z = Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
  return mean + z * sd;
}

function examPrediction(xTarget, weights, c, trainVec) {
  const nearestTrain = trainVec[math.argmin(math.abs(trainVec - xTarget))];
  const aResp = meanPrediction(nearestTrain, weights, c);
  const xUnder = math.min(trainVec) === nearestTrain ? nearestTrain : trainVec[math.findIndex(trainVec, (d) => d === nearestTrain) - 1];
  const xOver = math.max(trainVec) === nearestTrain ? nearestTrain : trainVec[math.findIndex(trainVec, (d) => d === nearestTrain) + 1];
  const mUnder = meanPrediction(xUnder, weights, c);
  const mOver = meanPrediction(xOver, weights, c);
  const examOutput = math.round(aResp + ((mOver - mUnder) / (xOver - xUnder)) * (xTarget - nearestTrain), 3);
  return examOutput;
}


function trainALM(dat, c, lr, weights) {
    console.log('training')
  const almTrain = new Array(dat.input.length).fill(NaN);
  for (let i = 0; i < dat.input.length; i++) {
    console.log('i: ', i, ' dat.input[i]: ', dat.input[i], ' dat.vx[i]: ', dat.vx[i], ' c: ', c, ' lr: ', lr)
    weights = updateWeights(dat.input[i], dat.vx[i], weights, c, lr);
    const resp = math.round(meanPrediction(dat.input[i], weights, c),0);
    // round resp to 1 decimal place
    almTrain[i] = resp;
      weights = math.map(weights, (value) => {
        return value < 0 ? 0 : value;
      });
  }
  console.log('almTrain: ', almTrain)
  console.log(weights)
  return {almTrain, weights};
}

function trainTestALM(dat, c = 0.05, lr = 0.5, weights, testVec) {
  const almTrain = new Array(dat.length).fill(NaN);
  
  for (let i = 0; i < dat.length; i++) {
    weights = updateWeights(dat[i].input, dat[i].vx, weights, c, lr);
    const resp = meanPrediction(dat[i].input, weights, c);
    almTrain[i] = resp;
    weights = math.map(weights, (value) => {
      return value < 0 ? 0 : value;
    });
  }

  const almPred = testVec.map((value) => {
    return meanPrediction(value, weights, c);
  });

  const examPred = testVec.map((value) => {
    return examPrediction(value, weights, c, [1, ...math.sort(math.unique(dat.map((d) => d.input)))]);
  });
    
  return { almTrain, almPred, examPred };
}

// Modify the sim_data function to accept the dataset as an argument
// function sim_data(dat, c=0.5, lr=0.2, inNodes=7, outNodes=32, trainVec=[5,6,7]) {
//   inputNodes = math.range(1,7,inNodes).toArray();  
//   outputNodes = math.range(50,1600,outNodes).toArray(); 
//   wm = math.zeros(outputNodes.length, inputNodes.length)._data;
//   tt = trainTest_alm(dat, c, lr, wm, trainVec);
// }

function gen_train(trainVec, trainRep, noise) {
   let bandVec=[0,100,350,600,800,1000,1200];
   let ts = [];
   for (let i=0; i<trainRep; i++) {
       ts.push(...trainVec);
   }
    let mean = 0;
    let stdDev = 1;
   //let noiseVec = math.random([ts.length])._data;
   //noiseVec = math.multiply(noiseVec, noise)._data;
   //if(noise==0) {noiseVec=noiseVec*0}
   let inputArr = [];
   let vxArr = [];
   for (let i=0; i<ts.length; i++) {
       inputArr.push(ts[i]);
       vxArr.push(bandVec[ts[i]]);
       //vxArr.push(bandVec[ts[i]]+noiseVec[i]);
   }
   return {input: inputArr, vx: vxArr};
}




```

### Simulation

```{ojs}
//| panel: sidebar
//| code-fold: true
viewof c = Inputs.range([.0001, 2], {value: .00005, step: .05, label: "c value:"})
viewof lr = Inputs.range([.001, 2], {value: .05, step: .01, label: "lr value:"})

viewof n_inputNodes = Inputs.range([1, 50], {value: 7, step: 1, label: "N Input Nodes:"})
viewof n_outputNodes = Inputs.range([1, 200], {value: 32, step: 1, label: "N Output Nodes:"})

viewof weight_mean = Inputs.range([0, 1], {value: 0, step: .0005, label: "initial weight mean:"})
viewof weight_sd = Inputs.range([.00000001, 1], {value: .000001, step: .0001, label: "initial weight sd:"})

viewof trainRep = Inputs.range([4, 50], {value: 1, step: 1, label: "Train Reps:"})

 //inputNodes = Array.from({ length: n_inputNodes }, (_, i) => i + 1);
 inputNodes = Array.from({ length: n_inputNodes }, (_, i) => 1 + i * (7 - 1) / (n_inputNodes - 1));



start = 0;
end = 1800;
N_Steps = n_outputNodes; // replace with desired length
stepSize = (end - start) / (N_Steps - 1);
outputNodes = Array.from({ length: N_Steps }, (_, i) => start + i * stepSize);

console.log(inputNodes)
console.log(outputNodes)

 wm = outputNodes.map(() => {
  return inputNodes.map(() => randomNormal(weight_mean, weight_sd));
});

noise=0
viewof trainVec = Inputs.checkbox([1, 2, 3, 4, 5, 6], {value: [4,5,6], label: "Select training examples:"});
gd = gen_train(trainVec, trainRep, noise);

//trainVec = [1,2,4,5,6];
//gd = gen_train(trainVec, trainRep, noise)
// w2= updateWeights(4, 800, wm,c,lr)

talm = trainALM(gd, c, lr, wm);

//inputNodes = transpose(iN)
ia = inputActivation(inputX, c, inputNodes)
oa = outputActivation(inputX, wm, c)
mp = meanPrediction(inputX, wm, c)
// subtract constant inputX from inputNodes Array
// Math.exp((-1 * c) * Math.pow(inputX - inputNodes, 2));

```

```{ojs}
tdat = gd.vx.map((value, index) => {
  return { Trial: index, Vx: value, Response: talm.almTrain[index], Error: Math.abs(value -  talm.almTrain[index]) };
});

```

::: {layout-ncol="2"}
#### Vx Across Training

```{ojs}

Plot.plot({
  marks: [
    Plot.line(tdat, {
      x: "Trial",      // feature for the x channel
      y: "Response",     // feature for the y channel
      stroke: "Vx",     
    }),
  ],
  x: {label: "Trial Number"},
  y: {label: "Vx", domain: [0, 1800],grid: true},
  color: {legend: true, scheme: "Turbo",type: "categorical"},
  width: 400,
  height: 400
});
  //caption: html`Figure 1. This chart has a <i>fancy</i> caption.`

```

#### Training Error

```{ojs}
Plot.plot({
  marks: [
    Plot.line(tdat, {
      x: "Trial",      // feature for the x channel
      y: "Error",     // feature for the y channel
      stroke: "Vx",     // feature for the fill channel
    }),
  ],
  y: {label: "Error",grid: true},
  color: {legend: true, scheme: "Turbo",type: "categorical"},
  width: 400,
  height: 400
});
```
:::

### Weight Matrices

```{ojs}
//| code-fold: true
//| message: false
//| include: false
Plotly = require("https://cdn.plot.ly/plotly-latest.min.js")
//div = DOM.element('div');
P1=Plotly.newPlot("plot-canvas", [{
  z: wm,
  x: outputNodes,
  y: inputNodes,
  type: 'heatmap',
  colorscale: 'Viridis'
}],{width:500});

console.log(inputNodes)
console.log(outputNodes)

P2=Plotly.newPlot("plot-tw", [{
  z: talm.weights,
  x: inputNodes,
  y: outputNodes,
  type: 'heatmap',
  colorscale: 'Viridis'
}],{width:600});

```

::: {#Weight-Matrices layout="[[1,1],[1]]"}
#### Starting Weights
[]{#plot-canvas}

#### Final Weights
[]{#plot-tw}
:::

#### Input and Output layer activations

```{ojs}

 in_data = ia.map((value, index) => {
  return { Node: inputNodes[index], Activation: value };
});

 out_data = oa.map((value, index) => {
  return { Node: outputNodes[index], Activation: value };
});

viewof inputX = Inputs.range([1, 7], {value: 4, step: 1, label: "input value:"})

```

```{ojs}
//| label: Activations
//| fig-cap: Charts
//| fig-subcap: 
//|   - "**Input Activation**"
//|   - "Second"
//| layout-ncol: 2



Plot.plot({
  marks: [
    Plot.dot(in_data, {
      x: "Node",      // feature for the x channel
      y: "Activation",     // feature for the y channel
    }),
  ],
  width: 400,
  height: 200,
  title: "Input Activation Plot",
});
Plot.plot({
  marks: [
    Plot.dot(out_data, {
      x: "Node",      // feature for the x channel
      y: "Activation",     // feature for the y channel
    }),
  ],
  width: 400,
  height: 200,
  title: "Output Activation Plot",
});

```


```{ojs}
console.log(wm)
console.log(talm.weights)
```


### Testing

```{ojs}
//| panel: sidebar
viewof xV = Inputs.range(
  [1, 20], 
  {value: 1, step: 1, label: "x range:"}
)
viewof yV = Inputs.range(
  [1, 400], 
  {value: 1, step: 10, label: "y range:"}
)

dO = transpose(d)
filtered = dO.filter(function(dO) {
  return dO.x>=xV && dO.y >= yV;
})



```

::: panel-tabset
## Tab 1

```{ojs}
Plot.plot({
  marks: [
    Plot.dot(filtered, 
      { x: "x", y: "y"}, 
      { stroke: "black" }
    )
  ]
})
```

x: \${xV} y: \${yV}

## Tab 2

```{ojs}
//Plot = require("plot")
Plot.plot({
  marks: [
    Plot.line(transpose(d), 
      { x: "x", y: "y"}, 
      { stroke: "black" }
    )
  ]
})
```
:::

``` {{ojs}}
//| include: false
```
