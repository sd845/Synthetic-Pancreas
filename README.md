#### About
This repository contains code for the synthetic pancreas circuit models written in the Julia programming language. The circuit contains a glucose resposive switch that expresses RFP or GFP based on glucose levels in the blood sample. This code was generated using the JuGRN code generator in Julia from the Varnerlab GitHub Repository. 

#### Requirements
Julia v1.4 must be installed to run this code. 

In addition, the following packages are required: `DifferentialEquations.jl`, `DiffEqSensitivity.jl`, `NumericalIntegration.jl`, `DelimitedFiles.jl`, `JSON.jl`, `DataFrames.jl`, `LinearAlgebra.jl`, `PyPlot.jl`.

#### Running the model

To run the model on your local machine, you can download this repository as a zip file, clone or pull it by using the command (using command line)
```
$ git pull https://github.com/sd845/Synthetic-Pancreas.git
```
or
```
$ git clone https://github.com/sd845/Synthetic-Pancreas.git
```
Start Julia and navigate to your working directory in the Julia REPL. The folder `Glucose_Switch` contains code for the complete circuit and the folder `Sub-circuit` contains code for the GntR-GFP sub-circuit.

Once in the working directory, run the command (using Julia REPL) 
```
include("Driver.jl")
```
