# Power-Flexibility-Project
Some program codes for the RSF-DFG project "Development of Innovative Technologies and Tools for Flexibility Assessment and Enhancement of Future Power Systems"

## Bilevel_PowerSystem 

This folder consists of Python and Matlab program codes, which presents a bi-level optimization framework applied to study operating flexibility and optimize two-level power system performance with (1) increasing presence of distributed energy resources (DER) at the low-voltage level, and (2) variable wind power generation at the high-voltage level.

![Alt-текст](https://www.dropbox.com/s/lr84auynxpx14nh/Flex.png?dl=0)


* **BilevelTest_Matlab** - the folder, which includes MATLAB program files for classical bilevel optimization approach. It requires Mathpower, YALMIP and Gurobi. 
* **BilevelTest_Python** - the folder, which includes Python program files for classical bilevel optimization approach. It requires Pypower, Pyomo and Gurobi.
* **BilevelRL_Python** - the folder, which includes Python program files for bilevel reinforcement learning approach.


##  EPSTwoLevel
The MATLAB software has been enhanced for modeling the two-level power system. This algorithm in this solver developed by analogy with [TDNetGen](https://github.com/apetros/TDNetGen) allows modeling the two-level EPS that includes transmission and distribution parts of the grid.  As of today, a computational model of upper level (transmission part of the grid) has been developed for a 123-bus grid (flexibility_test_system_FTS-213). This scheme of the grid is intended for studying the flexibility and apart from  conventional scheme-and-operating parameters includes financial indices (prices, costs, investments) for different types of generation (hydropower plants gas-turbine units, PVs, wind farms, etc.), large energy storages, AC/DC lines, etc. Six 12 kV 33-bus distribution networks can be connected to the 123-bus grid. 
Bus 1 is the feeding bus for every subsystem of the distribution grid.  Operating conditions are computed in two stages: First, the steady-state variables of distribution network are computed using Backward/Forward method, then the state variables of a 110 kV transmission network are computed using the Newton-Raphson method.   
