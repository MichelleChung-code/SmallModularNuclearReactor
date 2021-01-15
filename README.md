# Capstone_2020

For final year Chemical Engineering Capstone Project

ENCH 511 and ENCH 531, University of Calgary  

PFD_Final_DBM.vsdx : 

Visio file containing process flow diagram.Latest revision: Dec.7th 2020

## SYNGAS FOLDER:
Compositions&Data final.xlsx: 

Excel file containing information on material streams and calculations for PFR. 

Latest Revision: Dec.5th 2020

BFD_Final.vsdx: 

Visio file containing block flow diagram.
Latest revision: Dec.6th 2020

Syngas_Combined_Final_DBM.vsym: 

Integrated file of simulation of SMR output and syngas process. Latest revision: Dec. 7th 2020
## ECONOMICS FOLDER:
DBM_Economics.xlsx

Excel file containing the equipment costs, fixed capital investment, total capital investment, and timeline.

###### Econ_Matlab Folder
CashFlowDiagram.m

MATLAB file for plotting the DBM cashflow diagrams, plots the cashflow data from: 

Cashflows.csv containing the net cashflows to plot in the cashflow diagram

Latest revision: Dec.8th 2020

###### Econ_Python Folder
run_economics_main.py

Python file which builds base case cashflows and runs and plots sensitivity analysis results (includes IRR, NPV, and discounted payback period calculations) using supporting files.  These files are:
* base_case.csv containing the base case revenues and expenses.  To run the code as is, this should be placed in its own mfs folder, external to the source code folder.
* ProfitabilityAnalysis.py containing class with functions to calculate NPV, IRR, and discounted payback period
* SensitivityAnalysis.py containing class with functions to run the sensitivity analysis and produce the resulting graphs, used in the DBM

Latest revision: Dec.8th 2020

## SMR FOLDER:
nuclear_SMR_core_pressure_drop.py:

Python file calculating pressure drop across the SMR bed

Latest revision: Dec.2nd 2020

###### Nodal_Approach Folder
NeutronKineticsRun.m:

MATLAB file which runs the supporting files for the nodal approach SMR model.  These files are:
* NeutronKinetics.m containing the kinetics and thermal hydraulics models 
* Initial Values.csv containing the initial guesses for NeutronKinetics.m
* NeutronKineticsPlotting.m which plots model results

Latest revision: Dec.8th 2020

###### Prelim_Model Folder 
CelsiusToFahrenheit.m 

MATLAB file containing function to convert from degrees celcius to fahrenheit

Latest revision: Nov.10th 2020

UnitConversionsFactor.m 

MATLAB file containing various constants for typical unit conversions

Latest revision: Nov.10th 2020

simplified_reactor_core_main.m 

MATLAB file which runs the supporting files for the preliminary SMR Simulink model.  These files are:
* Concentrationlambdasum.slx contains the Simulink calculating the neutron Concentrations
* Pthmodel.slx  contains the Simulink calculating the Power output of the model
* rhomodel.slx contains the Simulink calculating the reactivity of the system
* SMR_simplified_core.slx contains the Simulink of the simplified reactor interations
* T_outmodel.slx contains the Simulink of the Temperature outlet of the model
* Tfmodel.slx contains the Simulink of the Temperature of the fuel elements
* Theta1model.slx contains the Simulink of the first node of the temperature model
* Theta2Model.slx contains the Simulink of the second node of the temperature model

Latest revision: Nov.19th 2020
