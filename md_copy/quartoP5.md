---
title: quarto p5
date: last-modified
categories: [Simulation,ALM,OJS,Interactive]
code-fold: true
code-tools: true
execute: 
  warning: false
--- 





```{ojs}

P5 = require("p5")
function* createSketch(sketch) {
  const element = DOM.element('div');
  yield element;
  const instance = new P5(sketch, element, true);
  try {
    while (true) {
      yield element;
    }
  } finally {
    instance.remove();
  }
}
createSketch(s2 => {
  
    s2.setup = function() {
      s2.createCanvas(746, 300);
      s2.textFont('Courgette');
      s2.textStyle(s2.BOLD);
      s2.textAlign(s2.CENTER, s2.CENTER)

      s2.button = s2.createButton('clear');
      s2.button.mousePressed(s2.clearCanvas);
        s2.text('Click and drag to draw', s2.width/2, s2.height/10);

    };
    s2.draw = function() {
    if (s2.mouseIsPressed) {
    s2.fill(0);
    s2.ellipse(s2.mouseX, s2.mouseY, 10, 10);
    } else {
   //s.fill(255);
    }
  // add text input

    };

  // add button to clear canvas
  s2.clearCanvas = function() {
    s2.clear();
  };
  // add text
  // add slider
  }
)

```



### ALM 

```{ojs}



// Javascript version of Shiny app for ALM
// P5 = require("p5");

// function inputActivation(x, weights, bias) {
//   let activation = 0;
//   for (let i = 0; i < x.length; i++) {
//     activation += x[i] * weights[i];
//   }
//   return activation + bias;
// }

// function outputActivation(input, threshold) {
//   if (input >= threshold) {
//     return 1;
//   } else {
//     return 0;
//   }
// }

function meanPrediction(x, weights, bias) {
}

function examPrediction(x, weights, bias) {

}

function updateWeights(x, y, weights, bias, associationParameter, updateParameter) {
  let prediction = meanPrediction(x, weights, bias);
  let error = y - prediction;
  for (let i = 0; i < x.length; i++) {
    weights[i] += associationParameter * error * x[i];
  }
  bias += updateParameter * error;
  return {weights, bias};
}

function learnALM(x, y, weights, bias, associationParameter, updateParameter, nRep) {
  for (let i = 0; i < nRep; i++) {
    let updated = updateWeights(x, y, weights, bias, associationParameter, updateParameter);
    weights = updated.weights;
    bias = updated.bias;
  }
  return {weights, bias};
}

// function* createSketch(sketch) {
//   const element = DOM.element('div');
//   yield element;
//   const instance = new P5(sketch, element, true);
//   try {
//     while (true) {
//       yield element;
//     }
//   } finally {
//     instance.remove();
//   }
// }

createSketch(s => {
//   let assocSlider, updateSlider, trainRepSlider, noiseSlider;
//   let trainItem1, trainItem2, trainItem3;
//   let nRepInput;
//   let runButton;

  s.setup = function() {
    s.createCanvas(500, 500);
    s.textFont('Courgette');
    s.textStyle(s.BOLD);
    s.textAlign(s.CENTER, s.CENTER)

    s.button = s.createButton('clear');
    s.button.mousePressed(s.clearCanvas);
    s.text('Click and drag to draw', s.width/2, s.height/10);

    s.input = s.createInput();
    s.input.position(s.width/2 - 50, s.height/10 + 20);
    s.input.size(100, 25);

    s.slider = s.createSlider(0, 255, 100);
    s.slider.position(s.width/2 - 50, s.height/10 + 50);
    s.slider.style('width', '100px');

    s.assocSlider = s.createSlider(0, 1, 0.5, 0.01);
    s.updateSlider = s.createSlider(0, 1, 0.5, 0.01);
    s.trainRepSlider = s.createSlider(0, 50, 25);
    s.noiseSlider = s.createSlider(0, 1, 0.1, 0.01);

    s.trainItem1 = s.createCheckbox('Item 1', false);
    s.trainItem2 = s.createCheckbox('Item 2', false);
    s.trainItem3 = s.createCheckbox('Item 3', false);

    s.nRepInput = s.createInput();

    s.runButton = s.createButton('Run Simulation');
    s.runButton.mousePressed(s.runSimulation);

    };

    s.draw = function() {
    s.background(s.slider.value());
    s.text(s.input.value(), s.width/2, s.height/2);
    };

    s.clearCanvas = function() {
    s.clear();
    };

    s.inputActivation = function(inputs) {
    // Code for input activation function
   let activatedInputs = [];
    for (let i = 0; i < inputs.length; i++) {
    let input = inputs[i];
    // Using the sigmoid activation function
    activatedInputs[i] = 1 / (1 + Math.exp(-input));
    }
    return activatedInputs;




    };

    s.outputActivation = function(output) {
    // Code for output activation function
    


    };

    s.meanPrediction = function(inputs) {
    // Code for mean prediction function
      return s.outputActivation(s.inputActivation(x, weights, bias), 0);

    };

    s.examPrediction = function(inputs) {
    // Code for exam prediction function
      let prediction = s.meanPrediction(x, weights, bias);
  let noise = Math.random() < s.noiseSlider.value() / 100 ? 1 : 0;
  return prediction ^ noise;
    };

    s.updateWeights = function(inputs, prediction, target) {
    // Code for updating weights
     let prediction = s.meanPrediction(x, weights, bias);
  let error = y - prediction;
  for (let i = 0; i < x.length; i++) {
    weights[i] += associationParameter * error * x[i];
  }
  bias += updateParameter * error;
  return {weights, bias};
    };

    s.learnALM = function(inputs, target) {
    // Code for learning in the ALM model
    };

    s.runSimulation = function() {
    // Code for running the simulation
    };
    });







```







s

//| echo: false
//| eval: false
//| include: false


// include = false
// javascript version of Shiny app for ALM

//p5 js

//   P5 = require("p5")

// function* createSketch(sketch) {
//   const element = DOM.element('div');
//   yield element;
//   const instance = new P5(sketch, element, true);
//   try {
//     while (true) {
//       yield element;
//     }
//   } finally {
//     instance.remove();
//   }
// }

// const sketch = s => {
//   var assocSlider, updateSlider, trainRepSlider, noiseSlider;
//   var trainItem1, trainItem2, trainItem3;
//   var nRepInput;
//   var runButton;

//   s.setup = function() {
//     var canvas = s.createCanvas(500, 500);
//     canvas.parent("canvas");

//     assocSlider = s.select("#assoc");
//     updateSlider = s.select("#update");
//     trainRepSlider = s.select("#trainRep");
//     noiseSlider = s.select("#Noise");
//     trainItem1 = s.select("#trainItem1");
//     trainItem2 = s.select("#trainItem2");
//     trainItem3 = s.select("#trainItem3");
//     nRepInput = s.select("#nRep");
//     runButton = s.select("#run");
//     runButton.mousePressed(s.runSimulation);
//   };

//   s.runSimulation = function() {
//     var assoc = assocSlider.value();
//     var update = updateSlider.value();
//     var trainRep = trainRepSlider.value();
//     var noise = noiseSlider.value();
//     var nRep = nRepInput.value();

//     var x_plotting = [];
//     var y_plotting = [];
//     for (var i = 0; i < nRep; i++) {
//       // run ALM simulation
//     var assoc = assocSlider.value() / 100;
//     var update = updateSlider.value() / 100;
//     var trainRep = trainRepSlider.value();
//     var noise = noiseSlider.value() / 100;
//     var nRep = parseInt(nRepInput.value());

// // Initialize activation values to 0
// var a = [];
// for (var i = 0; i < 12; i++) {
//   a[i] = 0;
// }

// // Train items
// a[10] = trainRep / 100;
// a[11] = trainRep / 100;
// a[12] = trainRep / 100;

// // Run the ALM simulation
// for (var i = 0; i < nRep; i++) {
//   for (var j = 0; j < 12; j++) {
//     if (j == 10) {
//       a[j] = assoc * a[j] + update * (a[11] + a[12]);
//     } else if (j == 11) {
//       a[j] = assoc * a[j] + update * (a[10] + a[12]);
//     } else if (j == 12) {
//       a[j] = assoc * a[j] + update * (a[10] + a[11]);
//     } else {
//       a[j] = assoc * a[j];
//     }
//     // Add noise
//     a[j] += (Math.random() - 0.5) * noise;
//     // Ensure activation is within 0 and 1
//     a[j] = Math.max(0, Math.min(1, a[j]));
//   }
// }
     
//       // store the results for plotting
// x_plotting[i] = [];
// y_plotting[i] = [];
// for (var j = 10; j <= 12; j++) {
// x_plotting[i][j] = i;
// y_plotting[i][j] = a[j];
// }
// }
//       }
//     }
// // plot the results
// s.background(255);
// s.strokeWeight(2);
// for (var i = 0; i < nRep; i++) {
//   for (var j = 10; j <= 12; j++) {
//     if (j == 10 && trainItem1.checked()) {
//       s.stroke(0, 0, 255);
//       s.point(x_plotting[i][j], y_plotting[i][j]);
//     } else if (j == 11 && trainItem2.checked()) {
//       s.stroke(0, 255, 0);
//       s.point(x_plotting[i][j], y_plotting[i][j]);
//     } else if (j == 12 && trainItem3.checked()) {
//       s.stroke(255, 0, 0);
//       s.point(x_plotting[i][j], y_plotting[i][j]);
//     }
//   }
// }
// };

// s.draw = function() {};
// };"
// }
// }
// }
// };
// });"

                      
// function inputActivation(xTarget, associationParameter) {
//   return Math.exp(-1 * associationParameter * Math.pow(100 * xTarget - 100 * xPlotting, 2));
// }

// function outputActivation(xTarget, weights, associationParameter) {
//   let result = 0;
//   for (let i = 0; i < weights.length; i++) {
//     result += weights[i] * inputActivation(xTarget, associationParameter);
//   }
//   return result;
// }

// function meanPrediction(xTarget, weights, associationParameter) {
//   let probability = outputActivation(xTarget, weights, associationParameter) / 
//     sum(outputActivation(xTarget, weights, associationParameter));
//   let result = 0;
//   for (let i = 0; i < yPlotting.length; i++) {
//     result += yPlotting[i] * probability[i];
//   }
//   return result;
// }

// // function to generate exam predictions
// function examPrediction(xTarget, weights, associationParameter) {
//   let trainVec = Array.from(new Set(xLearning)).sort();
//   let nearestTrain = trainVec[findClosestIndex(trainVec, xTarget)];
//   let aresp = meanPrediction(nearestTrain, weights, associationParameter);
//   let xUnder = (trainVec[0] === nearestTrain) ? nearestTrain : 
//     trainVec[trainVec.indexOf(nearestTrain) - 1];
//   let xOver = (trainVec[trainVec.length - 1] === nearestTrain) ? nearestTrain : 
//     trainVec[trainVec.indexOf(nearestTrain) + 1];
//   let mUnder = meanPrediction(xUnder, weights, associationParameter);
//   let mOver = meanPrediction(xOver, weights, associationParameter);
//   let examOutput = round((aresp + (mOver - mUnder) / (xOver - xUnder) * 
//     (xTarget - nearestTrain)), 3);
//   return examOutput;
// }

// function updateWeights(xNew, yNew, weights, associationParameter, updateParameter) {
//   let yFeedbackActivation = Math.exp(-1 * associationParameter * Math.pow(yNew - yPlotting, 2));
//   let xFeedbackActivation = outputActivation(xNew, weights, associationParameter);
//   let updatedWeights = [];
//   for (let i = 0; i < weights.length; i++) {
//     updatedWeights[i] = weights[i] + updateParameter * (yFeedbackActivation - xFeedbackActivation) *
//       inputActivation(xNew, associationParameter);
//     if (updatedWeights[i] < 0) {
//       updatedWeights[i] = 0;
//     }
//   }
//   return updatedWeights;
// }

// function learnAlm(yLearning, associationParameter = 0.05, updateParameter = 0.5) {
//   let weights = Array(yPlotting.length * xPlotting.length).fill(0);
//   for (let i = 0; i < yLearning.length; i++) {
//     weights = updateWeights(xLearning[i], yLearning[i], weights, associationParameter, updateParameter)
// for (let j = 0; j < weights.length; j++) {
// if (weights[j] < 0) {
// weights[j] = 0;
// }
// }
// }
// let almPredictions = [];
// for (let i = 0; i < xPlotting.length; i++) {
// almPredictions.push(meanPrediction(xPlotting[i], weights, associationParameter));
// }
// let examPredictions = [];
// for (let i = 0; i < xPlotting.length; i++) {
// examPredictions.push(examPrediction(xPlotting[i], weights, associationParameter));
// }
// return {
// almPredictions: almPredictions,
// examPredictions: examPredictions
// };
// }










