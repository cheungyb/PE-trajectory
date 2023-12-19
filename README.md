# PE-trajectory
-----------------------------------------------------------------------------------------------------------------------------------------
Titile: Estimation of trajectory of protective efficacy in infectious disease prevention trials using recurrent event times
Authors: Yin Bun Cheung, Xiangmei Ma, K. F. Lam, Chee Fu Yung, Paul Milligan
-----------------------------------------------------------------------------------------------------------------------------------------
Application:
"Functions4ParPKPD.R": program includes all functions needed in the model estimation.
"Estimation4ParPKPD.R": codes for estimation of an extended AG model with protective efficacy captured by a 4-parameter PK/PD function.
"Example_data_2023.csv": a (hypothetical) recurrent events dataset for demonstration of data format and for testing only.
-----------------------------------------------------------------------------------------------------------------------------------------
Simulation evaluation of the method:
"Simu_Functions4ParPKPD.R": program includes all functions needed in the simulation study
"Simu_Program4ParPKPD.R": program includes recurrent event times genaration from an extended AG model with a 4-par PK/PD function,
			                    and estimation of the extended AG model with a 4-parameter PK/PD function.
"Simu_Program5ParPKPD.R": program includes recurrent event times genaration from an extended AG model with a 5-par PK/PD function, 
			                    and estimation of the extended AG model with a 4-parameter PK/PD function. 
"Simu_ProgramLinear.R": program includes recurrent event times genaration from an extended AG model with a piecewise linear function, 
			                  and estimation of the extended AG model with a 4-parameter PK/PD model function. 
"simu_snow_setting_4parPKPD.R": codes for simulation scenarios under the setting of data generation based on a 4-par PK/PD function.
"simu_snow_setting_5parPKPD.R": codes for simulation scenarios under the setting of data generation based on a 5-par PK/PD function.
"simu_snow_setting_linear.R": codes for simulation scenarios under the setting of data generation based on a piecewise linear function.
-----------------------------------------------------------------------------------------------------------------------------------------

