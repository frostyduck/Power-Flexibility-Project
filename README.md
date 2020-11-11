# Power-Flexibility-Project
Some program codes for the RSF-DFG project "Development of Innovative Technologies and Tools for Flexibility Assessment and Enhancement of Future Power Systems"

## Bilevel_PowerSystem 

This folder consists of Python and Matlab program codes, which presents a bi-level optimization framework applied to study operating flexibility and optimize two-level power system performance with (1) increasing presence of distributed energy resources (DER) at the low-voltage level, and (2) variable wind power generation at the high-voltage level.

* **BilevelTest_Matlab** - the folder, which includes MATLAB program files for classical bilevel optimization approach. It requires Mathpower, YALMIP and Gurobi. 
* **BilevelTest_Python** - the folder, which includes Python program files for classical bilevel optimization approach. It requires Pypower, Pyomo and Gurobi.
* **BilevelRL_Python** - the folder, which includes Python program files for bilevel reinforcement learning approach.


##  EPSTwoLevel
The MATLAB software solver has been enhanced for modeling the two-level power system. This algorithm in this solver developed by analogy with [TDNetGen](https://github.com/apetros/TDNetGen) allows modeling the two-level EPS that includes transmission and distribution parts of the grid.  As of today, a computational model of upper level (transmission part of the grid) has been developed for a 123-bus grid (flexibility_test_system_FTS-213). This scheme of the grid is intended for studying the flexibility and apart from  conventional scheme-and-operating parameters includes financial indices (prices, costs, investments) for different types of generation (hydropower plants gas-turbine units, PVs, wind farms, etc.), large energy storages, AC/DC lines, etc. Six 12 kV 33-bus distribution networks can be connected to the 123-bus grid. 
Bus 1 is the feeding bus for every subsystem of the distribution grid.  Operating conditions are computed in two stages: First, the steady-state variables of distribution network are computed using Backward/Forward method, then the state variables of a 110 kV transmission network are computed using the Newton-Raphson method.   

##  Flexibility_HeatPowerSystem. 

This folder consists of Python program file - Flexibility_MILP.py, to investigate the potential for the enhancement of integrated (heat+power) system flexibility by using a heating network 

We try to evaluate several optimization approaches based on MILP. We used the model to study the combined heat and power dispatch mode implemented on the integrated 6-bus electricity and 3-bus district heating network. This system comprises a conventional electricity generator G1, a wind producer W1, HP and an extraction CHP. The objective of the integrated heat and power dispatch (IHPD) is to mini-mize production costs of the overall energy system. In order to assess the impact of the district heating network on the power system, IHPD model is compared to a conventional economic dispatch (CED), which presents a co-optimization model for heat and electricity in which the district heating network is not modeled. This simulation shows that the IHPD model improves economic efficiency of the overall energy system and increases wind utilization by enhancing the flexibility of the district heating network.  
