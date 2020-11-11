"""
Код реализует модель частично-целочисленного линейного программирования (MILP) для реализации оперативной гибкости для интегрированных систем электро- и теплоснабжения.
В примере рассмотрена 6-ти узловая схема, которая содержит
Узел 1 - ветрогенератор и тепловой насос
Узел 2 - ТЭЦ
Узел 6 - обычный электрический генератор


Технические характеристики сети и блоков генерации получены из следующей статьи:

Z. Li, W. Wu, M. Shahidehpour, J. Wang, and B. Zhang, “Combined
heat and power dispatch considering pipeline energy storage of district
heating network,” IEEE Trans. Sustain. Energy, vol. 7, no. 1, pp. 12–22,
2016.

Данные по ветрогенерации взяты из https://sites.google.com/site/datasmopf/.

Пиковые часы потребности в тепле возникают в вечернее и утреннее время из-за более низких температур окружающей среды,
тогда как часы пиковой нагрузки на электроэнергию присутствуют в течение всего дня .

Более подробно математическая постановка этой задачи в статье
L. Mitridati and J. A. Taylor, "Power Systems Flexibility from District Heating Networks," 2018 Power Systems Computation Conference (PSCC), Dublin, 2018, pp. 1-7.
doi: 10.23919/PSCC.2018.8442617

"""

import os
import pandas as pd
import scipy.stats as sp
import matplotlib.pyplot as plt
import seaborn as sb
sb.set_style('ticks')

import gurobipy as gb
import itertools as it

import numpy as np
import datetime as dtt
import time as tt

# -------------- ЗАДАНИЕ ИСХОДНЫХ ПАРАМЕТРОВ ДЛЯ 6-ТИ УЗЛОВОЙ СХЕМЫ ---------------------

# %% ЗАДАНИЕ ИНДЕКСОВ


#%%


# %% indexes

T = 24
N = 6
G = 2  # elec generator
W = 1  # wind
H = 0  # chp
HO = 1  # heat only
HS = 0  # storage
ES = 0
HP = 1  # heat pump
HES=2
P = 2
L = 7
LL = 500

# %%

time = ['t{0:02d}'.format(t + 0) for t in range(T)]
time_day = ['t{0:02d}'.format(t) for t in range(24)]  # optimization periods indexes (24 hours)
time_list = np.arange(T)
node = ['N{0}'.format(n + 1) for n in range(N)]
line = ['L{0}'.format(l + 1) for l in range(L)]
pipe_supply = ['PS{0}'.format(p + 1) for p in range(P)]
pipe_return = ['PR{0}'.format(p + 1) for p in range(P)]
gen = ['G{0}'.format(g + 1) for g in range(G)]
heat_storage = ['HS{0}'.format(s + 1) for s in range(HS)]
elec_storage = ['ES{0}'.format(s + 1) for s in range(ES)]
heat_pump = ['HP{0}'.format(h + 1) for h in range(HP)]
wind = ['W{0}'.format(w + 1) for w in range(W)]
heat_only = ['HO{0}'.format(ho + 1) for ho in range(HO)]
heat_exchanger_station = ['HES{0}'.format(h+1) for h in range(HES)]
CHP_sorted = {'ex': ['CHP1'], 'bp': []}  # CHPs indexes sorted by type: extraction or backpressure
CHP = list(it.chain.from_iterable(CHP_sorted.values()))


# %%

pipe = ['P{0}'.format(p + 1) for p in range(P)]

pipe_connexion = {'P1': ['N1', 'N2'], 'P2': ['N2', 'N3']}
pipe_start = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
pipe_end = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
for n in node:
    for p in pipe:
        if pipe_connexion[p][0] == n:
            pipe_start[n].append(p)
        if pipe_connexion[p][1] == n:
            pipe_end[n].append(p)

# %% DH network parameters

pipe_maxflow = {p: 1000 for p in pipe}


def produit(list):  # takes a list as argument
    p = 1
    for x in list:
        p = p * x
    return (p)



# heat market data
heat_station = CHP + heat_only + heat_pump
elec_station = CHP + gen + wind + heat_pump
producers = CHP + heat_only + heat_storage + heat_pump + gen + wind + elec_storage

producers_node = {'N1': ['HP1', 'W1'], 'N2': ['CHP1', 'HO1', 'G2'], 'N3': [], 'N4': [], 'N5': [], 'N6': ['G1']}
heat_exchanger_station_node = {'N1':[],'N2':[],'N3':['HES1'],'N4':[],'N5':[],'N6':['HES2']}
elec_station_node = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
heat_station_node = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
CHP_node = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
CHP_sorted_node = {'N1': {'ex': [], 'bp': []}, 'N2': {'ex': [], 'bp': []}, 'N3': {'ex': [], 'bp': []},
                   'N4': {'ex': [], 'bp': []}, 'N5': {'ex': [], 'bp': []}, 'N6': {'ex': [], 'bp': []}}
heat_pump_node = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
heat_only_node = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
heat_storage_node = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
gen_node = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
wind_node = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
elec_storage_node = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}

for n in node:
    for h in CHP:
        if h in producers_node[n]:
            CHP_node[n].append(h)
            heat_station_node[n].append(h)
            elec_station_node[n].append(h)
            if h in CHP_sorted['ex']:
                CHP_sorted_node[n]['ex'].append(h)
            if h in CHP_sorted['bp']:
                CHP_sorted_node[n]['bp'].append(h)
    for h in heat_only:
        if h in producers_node[n]:
            heat_only_node[n].append(h)
            heat_station_node[n].append(h)
    for h in heat_storage:
        if h in producers_node[n]:
            heat_storage_node[n].append(h)
            # heat_station_node[n].append(h)
    for h in heat_pump:
        if h in producers_node[n]:
            heat_pump_node[n].append(h)
            heat_station_node[n].append(h)
            elec_station_node[n].append(h)
    for h in gen:
        if h in producers_node[n]:
            gen_node[n].append(h)
            elec_station_node[n].append(h)
    for h in wind:
        if h in producers_node[n]:
            wind_node[n].append(h)
            elec_station_node[n].append(h)
    for h in elec_storage:
        if h in producers_node[n]:
            elec_storage_node[n].append(h)
            # elec_station_node[n].append(h)

pipe_supply_connexion = {'PS1': ['N1', 'N2'], 'PS2': ['N2', 'N3']}
pipe_supply_start = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
pipe_supply_end = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
for n in node:
    for p in pipe_supply:
        if pipe_supply_connexion[p][0] == n:
            pipe_supply_start[n].append(p)
        if pipe_supply_connexion[p][1] == n:
            pipe_supply_end[n].append(p)

pipe_return_connexion = {
pipe_return[p]: [pipe_supply_connexion[pipe_supply[p]][1], pipe_supply_connexion[pipe_supply[p]][0]] for p in range(P)}

pipe_return_start = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
pipe_return_end = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}

for n in node:
    for p in pipe_return:
        if pipe_return_connexion[p][0] == n:
            pipe_return_start[n].append(p)
        if pipe_return_connexion[p][1] == n:
            pipe_return_end[n].append(p)

line_connexion = {'L1': ['N1', 'N2'], 'L2': ['N2', 'N3'], 'L3': ['N3', 'N6'], 'L4': ['N6', 'N5'], 'L5': ['N5', 'N4'],
                  'L6': ['N4', 'N1'], 'L7': ['N3', 'N5']}
line_start = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
line_end = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
for n in node:
    for l in line:
        if line_connexion[l][0] == n:
            line_start[n].append(l)
        if line_connexion[l][1] == n:
            line_end[n].append(l)

        # %% loads

#heat=pd.read_csv("heat_demand_Tingbjerg_2015.csv",sep=";",decimal=",")
#heat_load = {'N1':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':},'N2':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':},'N3':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':}}
heat_load = {}
for t in time:
    heat_load['N1', t] = 0
    heat_load['N2', t] = 0
    heat_load['N4', t] = 0
    heat_load['N5', t] = 0
    heat_load['N6', t] = 0

#heat_load_csv = pd.read_csv("heat_demand_Tingbjerg_2015.csv", sep=";", decimal=",")

#heat_load_0 = {time_day[t]: heat_load_csv.ix[t + 1 * 24, 1] for t in range(24)}

#for t in time:
   # heat_load['N3', t] = int(heat_load_0[t]) * 100 / max(int(heat_load_0[t]) for t in time)

# heat_load['N3', 't07']= 93.354037267080734
# heat_load['N3', 't08']= 92.11180124223602
# heat_load['N3', 't09']= 91.987577639751549
# heat_load['N3', 't10']= 61.440993788819867
# heat_load['N3', 't11']= 51.459627329192543
# heat_load['N3', 't12']= 100

heat_load['N3', 't23'] = 82.422360248447205
heat_load['N3', 't00'] = 78.012422360248436
heat_load['N3', 't01'] = 72.298136645962728
heat_load['N3', 't02'] = 72.298136645962728
heat_load['N3', 't03'] = 72.385093167701854
heat_load['N3', 't04'] = 82.608695652173907 + 10
heat_load['N3', 't05'] = 83.354037267080734 + 10
heat_load['N3', 't06'] = 85.590062111801231 + 10
heat_load['N3', 't07'] = 86.086956521739125
heat_load['N3', 't08'] = 92.11180124223602
heat_load['N3', 't09'] = 71.987577639751549
heat_load['N3', 't10'] = 69.440993788819867
heat_load['N3', 't11'] = 66.459627329192543
heat_load['N3', 't12'] = 64.720496894409933
heat_load['N3', 't13'] = 63.354037267080734
heat_load['N3', 't14'] = 62.11180124223602
heat_load['N3', 't15'] = 61.987577639751549
heat_load['N3', 't16'] = 64.968944099378874
heat_load['N3', 't17'] = 64.472049689440993
heat_load['N3', 't18'] = 66.583850931677006
heat_load['N3', 't19'] = 110.987577639751549
heat_load['N3', 't20'] = 97.440993788819867
heat_load['N3', 't21'] = 95.459627329192543
heat_load['N3', 't22'] = 85.590062111801231

# heat_load[ 'N3', 't23'] = 82.422360248447205
# heat_load[ 'N3', 't00'] = 78.012422360248436
# heat_load[ 'N3', 't01'] = 72.298136645962728
# heat_load[ 'N3', 't02'] = 72.298136645962728
# heat_load[ 'N3', 't03'] = 72.385093167701854
# heat_load[ 'N3', 't04'] = 82.608695652173907+10
# heat_load[ 'N3', 't05' ] =  83.354037267080734+10
# heat_load[ 'N3', 't06'] = 85.590062111801231+10
# heat_load['N3', 't07']=  86.086956521739125
# heat_load['N3', 't08']= 92.11180124223602
# heat_load['N3', 't09']= 71.987577639751549
# heat_load['N3', 't10']= 69.440993788819867
# heat_load['N3', 't11']= 66.459627329192543
# heat_load['N3', 't12']= 64.720496894409933
# heat_load['N3', 't13']= 63.354037267080734
# heat_load['N3', 't14']= 62.11180124223602
# heat_load['N3', 't15']= 61.987577639751549
# heat_load['N3', 't16']= 64.968944099378874
# heat_load['N3', 't17']= 64.472049689440993
# heat_load['N3', 't18']= 66.583850931677006
# heat_load['N3', 't19']= 110.987577639751549
# heat_load['N3', 't20']= 97.440993788819867
# heat_load['N3', 't21']= 95.459627329192543
# heat_load['N3', 't22']= 85.590062111801231

# heat_load['N3','t00']=0
# heat_load['N3','t01']=0
# heat_load['N3','t02']=0
# heat_load['N3','t03']=0
# heat_load['N3','t04']=0
# heat_load['N3','t05']=0
# heat_load['N3','t06']=0
# heat_load['N3','t07']=0
# heat_load['N3','t08']=0
# heat_load['N3','t09']=0
# heat_load['N3','t10']=0
# heat_load['N3','t11']=0
# heat_load['N3','t12']=0
# heat_load['N3','t13']=0
# heat_load['N3','t14']=0
# heat_load['N3','t15']=0
# heat_load['N3','t16']=0
# heat_load['N3','t17']=0
# heat_load['N3','t18']=0
# heat_load['N3','t19']=0
# heat_load['N3','t20']=0
# heat_load['N3','t21']=0
# heat_load['N3','t22']=0
# heat_load['N3','t23']=0
#
# heat_load['N6','t00']=0
# heat_load['N6','t01']=0
# heat_load['N6','t02']=0
# heat_load['N6','t03']=0
# heat_load['N6','t04']=0
# heat_load['N6','t05']=0
# heat_load['N6','t06']=0
# heat_load['N6','t07']=0
# heat_load['N6','t08']=0
# heat_load['N6','t09']=0
# heat_load['N6','t10']=0
# heat_load['N6','t11']=0
# heat_load['N6','t12']=0
# heat_load['N6','t13']=0
# heat_load['N6','t14']=0
# heat_load['N6','t15']=0
# heat_load['N6','t16']=0
# heat_load['N6','t17']=0
# heat_load['N6','t18']=0
# heat_load['N6','t19']=0
# heat_load['N6','t20']=0
# heat_load['N6','t21']=0
# heat_load['N6','t22']=0
# heat_load['N6','t23']=0

# elec_load_mean = {'t00':780,'t01':750,'t02':730,'t03':770,'t04':800,'t05':850,'t06':1000,'t07':1200,'t08':1400,'t09':1300,'t10':1280,'t11':1250,'t12':1230,'t13':1100,'t14':1050,'t15':1000,'t16':980,'t17':950,'t18':1010,'t19':1100,'t20':980,'t21':930,'t22':850,'t23':830}
# elec_load = {'N1':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':},'N2':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':},'N3':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':}}

elec_load_IEEE = {'t00': 820, 't01': 820, 't02': 815, 't03': 815, 't04': 810, 't05': 850, 't06': 1000,
                  't07': 1150 + 100 * 1400 / 350, 't08': 1250 + 100 * 1400 / 350, 't09': 1250 + 100 * 1400 / 350,
                  't10': 1100 + 100 * 1400 / 350, 't11': 1000 + 100 * 1400 / 350, 't12': 1000 + 50 * 1400 / 350,
                  't13': 955 + 50 * 1400 / 350, 't14': 950 + 50 * 1400 / 350, 't15': 950 + 36 * 1400 / 350,
                  't16': 900 + 25 * 1400 / 350, 't17': 950, 't18': 1010, 't19': 1100, 't20': 1125, 't21': 1025,
                  't22': 950, 't23': 850}

elec_load = {}
for t in time:
    elec_load['N1', t] = 0
    elec_load['N2', t] = 0
    elec_load['N3', t] = elec_load_IEEE[t] * 0.20 * 350 / 1400
    elec_load['N4', t] = elec_load_IEEE[t] * 0.40 * 350 / 1400
    elec_load['N5', t] = elec_load_IEEE[t] * 0.40 * 350 / 1400
    elec_load['N6', t] = 0

# elec_load['N2','t00']=0
# elec_load['N2','t01']=0
# elec_load['N2','t02']=0
# elec_load['N2','t03']=0
# elec_load['N2','t04']=0
# elec_load['N2','t05']=0
# elec_load['N2','t06']=0
# elec_load['N2','t07']=0
# elec_load['N2','t08']=0
# elec_load['N2','t09']=0
# elec_load['N2','t10']=0
# elec_load['N2','t11']=0
# elec_load['N2','t12']=0
# elec_load['N2','t13']=0
# elec_load['N2','t14']=0
# elec_load['N2','t15']=0
# elec_load['N2','t16']=0
# elec_load['N2','t17']=0
# elec_load['N2','t18']=0
# elec_load['N2','t19']=0
# elec_load['N2','t20']=0
# elec_load['N2','t21']=0
# elec_load['N2','t22']=0
# elec_load['N2','t23']=0
#
# elec_load['N3','t00']=0
# elec_load['N3','t01']=0
# elec_load['N3','t02']=0
# elec_load['N3','t03']=0
# elec_load['N3','t04']=0
# elec_load['N3','t05']=0
# elec_load['N3','t06']=0
# elec_load['N3','t07']=0
# elec_load['N3','t08']=0
# elec_load['N3','t09']=0
# elec_load['N3','t10']=0
# elec_load['N3','t11']=0
# elec_load['N3','t12']=0
# elec_load['N3','t13']=0
# elec_load['N3','t14']=0
# elec_load['N3','t15']=0
# elec_load['N3','t16']=0
# elec_load['N3','t17']=0
# elec_load['N3','t18']=0
# elec_load['N3','t19']=0
# elec_load['N3','t20']=0
# elec_load['N3','t21']=0
# elec_load['N3','t22']=0
# elec_load['N3','t23']=0
#
# elec_load['N4','t00']=0
# elec_load['N4','t01']=0
# elec_load['N4','t02']=0
# elec_load['N4','t03']=0
# elec_load['N4','t04']=0
# elec_load['N4','t05']=0
# elec_load['N4','t06']=0
# elec_load['N4','t07']=0
# elec_load['N4','t08']=0
# elec_load['N4','t09']=0
# elec_load['N4','t10']=0
# elec_load['N4','t11']=0
# elec_load['N4','t12']=0
# elec_load['N4','t13']=0
# elec_load['N4','t14']=0
# elec_load['N4','t15']=0
# elec_load['N4','t16']=0
# elec_load['N4','t17']=0
# elec_load['N4','t18']=0
# elec_load['N4','t19']=0
# elec_load['N4','t20']=0
# elec_load['N4','t21']=0
# elec_load['N4','t22']=0
# elec_load['N4','t23']=0

# %% DH network parameters

Dt = 60 * 60  # sec/hour
# pipelines parameters

# length of pipes in in meter
length_pipe = {p: 500 for p in pipe_supply + pipe_return}  # in m

fanning_coeff = {p: 0.05 / 4 for p in pipe_supply + pipe_return}  # NO UNIT (moody chart)
radius_pipe = {p: 0.8 for p in pipe_supply + pipe_return}  # in m
water_density = 1000  # in kg/m3
water_thermal_loss_coeff = 20  # in W/m2.K
water_heat_capacity = 4.18 * 0.28  # in Wh/kg.K
pressure_loss_coeff = {
p: 2 * 16 / (np.pi ** 2) * fanning_coeff[p] * length_pipe[p] / (((radius_pipe[p] * 2) ** 5) * water_density) for p in
pipe_supply + pipe_return}

# pressure bounds in Pa
pr_return_max = {n: 5000 for n in node}
pr_return_min = {n: 50 for n in node}
pr_supply_max = {n: 5000 for n in node}
pr_supply_min = {n: 50 for n in node}
pr_diff_min = {n: 1000 for n in node}

# OUTER APPROX SLICES

pr_return_slice = [pr_return_min[n] + l * (pr_return_max[n] - pr_return_min[n]) / LL for l in range(LL)]
pr_supply_slice = [pr_supply_min[n] + l * (pr_supply_max[n] - pr_supply_min[n]) / LL for l in range(LL)]

# mass flow bounds

# mf_pipe_max = {p:2*0.062*water_density*np.pi*(radius_pipe[p]**2) for p in pipe_supply+pipe_return}
mf_pipe_min = {p: 50 for p in pipe_supply + pipe_return}

mf_pipe_max = {}
mf_pipe_max['PS1'] = 300  # 300
mf_pipe_max['PS2'] = 300  # 300

mf_pipe_max['PR1'] = mf_pipe_max['PS1']
mf_pipe_max['PR2'] = mf_pipe_max['PS2']

mf_HS_max = {}
for h in heat_station + heat_storage:
    mf_HS_max[h] = 0  # 500
mf_HS_max['CHP1'] = 301  # 500
mf_HS_max['HP1'] = 300  # 500
mf_HS_min = {}  # carefull for HEAT STORAGE!!!!
for h in heat_station:
    mf_HS_min[h] = 0
for h in heat_storage:
    mf_HS_min[h] = 0
mf_HES_max = {}
mf_HES_max['N1'] = 0
mf_HES_max['N2'] = 0
mf_HES_max['N3'] = 300  # 300
mf_HES_max['N4'] = 0
mf_HES_max['N5'] = 0
mf_HES_max['N6'] = 0
mf_HES_min = {}
mf_HES_min['N1'] = 0
mf_HES_min['N2'] = 0
mf_HES_min['N3'] = 150
mf_HES_min['N4'] = 0
mf_HES_min['N5'] = 0
mf_HES_min['N6'] = 0

# time delay
time_delay_max = {p: int(water_density * length_pipe[p] * np.pi * (radius_pipe[p] ** 2) / (Dt * mf_pipe_min[p])) + 1 for
                  p in pipe_supply + pipe_return}  # in hour
time_delay_range = {p: [x for x in range(time_delay_max[p] + 1)] for p in pipe_supply + pipe_return}

time_init = {p: ['t00-{0}'.format(time_delay_max[p] - k) for k in time_delay_range[p][:-1]] for p in
             pipe_supply + pipe_return}
time_extended = {p: time_init[p] + time for p in pipe_supply + pipe_return}

tau_new_min = {
p: min(0, 1 - 2 * time_delay_max[p] * water_thermal_loss_coeff / (water_heat_capacity * water_density * radius_pipe[p]))
for p in pipe_supply + pipe_return}
tau_new_max = {p: 1 for p in pipe_supply + pipe_return}

# time delay
time_delay_max = {p: int(water_density * length_pipe[p] * np.pi * (radius_pipe[p] ** 2) / (Dt * mf_pipe_min[p])) + 1 for
                  p in pipe_supply + pipe_return}  # in hour
time_delay_range = {p: [x for x in range(time_delay_max[p] + 1)] for p in pipe_supply + pipe_return}

time_init = {p: ['t00-{0}'.format(time_delay_max[p] - k) for k in time_delay_range[p][:-1]] for p in
             pipe_supply + pipe_return}
time_extended = {p: time_init[p] + time for p in pipe_supply + pipe_return}

tau_new_min = {
p: min(0, 1 - 2 * time_delay_max[p] * water_thermal_loss_coeff / (water_heat_capacity * water_density * radius_pipe[p]))
for p in pipe_supply + pipe_return}
tau_new_max = {p: 1 for p in pipe_supply + pipe_return}

# Initial data: from XX last hours of the day before

mf_pipe_init = {}
for p in pipe_supply + pipe_return:
    for t in time_init[p]:
        mf_pipe_init[p, t] = mf_pipe_min[p]

T_in_init = {}
for p in pipe_supply:
    for t in time_init[p]:
        T_in_init[p, t] = 90 + 273.15
for p in pipe_return:
    for t in time_init[p]:
        T_in_init[p, t] = 60 + 273.15

# temperature bounds
T_return_max = {}
T_return_min = {}
T_supply_max = {}
T_supply_min = {}

for t in time_extended['PS1']:
    for n in node:
        T_supply_max[n, t] = 120 + 273.15
        T_supply_min[n, t] = 90 + 273.15
        T_return_max[n, t] = 60 + 273.15
        T_return_min[n, t] = 30 + 273.15

T_in_max = {}
T_in_min = {}
T_out_max = {}
T_out_min = {}

for p in pipe_supply:
    for t in time_extended[p]:
        T_in_max[p, t] = T_supply_max[pipe_supply_connexion[p][0], t]
        T_in_min[p, t] = T_supply_min[pipe_supply_connexion[p][0], t]
        T_out_max[p, t] = T_supply_max[pipe_supply_connexion[p][1], t]
        T_out_min[p, t] = T_supply_min[pipe_supply_connexion[p][1], t]
for p in pipe_return:
    for t in time_extended[p]:
        T_in_max[p, t] = T_return_max[pipe_return_connexion[p][0], t]
        T_in_min[p, t] = T_return_min[pipe_return_connexion[p][0], t]
        T_out_max[p, t] = T_return_max[pipe_return_connexion[p][1], t]
        T_out_min[p, t] = T_return_min[pipe_return_connexion[p][1], t]

    # %% cost

alpha = {}
for t in time:
    alpha['CHP1', t] = 12.5
    alpha['HO1', t] = 80
    alpha['G1', t] = 11
    alpha['G2', t] = 33
    alpha['W1', t] = 0.0001

# %% technical characteristics

heat_maxprod = {'CHP1': 100, 'HP1': 150, 'HO1': 0}  # HP 150
rho_elec = {'CHP1': 2.4}  # efficiency of the CHP for electricity production
rho_heat = {'CHP1': 0.25, 'HO1': 1}  # efficiency of the CHP for heat production
r_min = {'CHP1': 0.6}  # elec/heat ratio (flexible in the case of extraction units)
CHP_maxprod = {'CHP1': 250}

COP = {'HP1': 2.5}

storage_loss = {h: 10 for h in heat_storage + elec_storage}
storage_init = {h: 1000 for h in heat_storage + elec_storage}
storage_rho_plus = {h: 1.1 for h in heat_storage + elec_storage}  # >=1
storage_rho_moins = {h: 0.9 for h in heat_storage + elec_storage}  # <=1
storage_maxcapacity = {h: 2000 for h in heat_storage + elec_storage}
storage_maxprod = {h: 500 for h in heat_storage + elec_storage}

elec_maxprod = {'CHP1': 1000, 'G1': 180, 'G2': 100, 'W1': 500}  # known

#pp = pd.Panel({'W' + str(i + 1): pd.read_csv('scen_zone{0}.out'.format(i + 1), index_col=0) for i in range(W)})

# %%

wind_scenario = {}

# for w in wind:
#    for t in range(T):
#        wind_scenario[w,time[t]] = pp[w,t+1,'V2']

wind_scenario['W1', 't00'] = 0.4994795107848
wind_scenario['W1', 't01'] = 0.494795107848
wind_scenario['W1', 't02'] = 0.494795107848
wind_scenario['W1', 't03'] = 0.505243011484
wind_scenario['W1', 't04'] = 0.53537368424
wind_scenario['W1', 't05'] = 0.555562455471
wind_scenario['W1', 't06'] = 0.628348636916
wind_scenario['W1', 't07'] = 0.6461954549
wind_scenario['W1', 't08'] = 0.622400860956
wind_scenario['W1', 't09'] = 0.580111023006
wind_scenario['W1', 't10'] = 0.714935503018
wind_scenario['W1', 't11'] = 0.824880140759
wind_scenario['W1', 't12'] = 0.416551027874
wind_scenario['W1', 't13'] = 0.418463919582
wind_scenario['W1', 't14'] = 0.39525842857
wind_scenario['W1', 't15'] = 0.523097379857
wind_scenario['W1', 't16'] = 0.476699300008
wind_scenario['W1', 't17'] = 0.626077589123
wind_scenario['W1', 't18'] = 0.684294396661
wind_scenario['W1', 't19'] = 0.0598119722706
wind_scenario['W1', 't20'] = 0.0446453658917
wind_scenario['W1', 't21'] = 0.485237701755
wind_scenario['W1', 't22'] = 0.49466503395
wind_scenario['W1', 't23'] = 0.4993958131342

# supply_maxtwmp = {'N1':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':},'N2':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':},'N3':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':}}
# supply_mintemp = {'N1':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':},'N2':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':},'N3':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':}}
# return_maxtwmp = {'N1':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':},'N2':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':},'N3':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':}}
# return_mintemp = {'N1':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':},'N2':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':},'N3':{'t00':,'t01':,'t02':,'t03':,'t04':,'t05':,'t06':,'t07':,'t08':,'t09':,'t10':,'t11':,'t12':,'t13':,'t14':,'t15':,'t16':,'t17':,'t18':,'t19':,'t20':,'t21':,'t22':,'t23':}}

# susceptance in siemens
B = {}
line_maxflow = {}

B['L1'] = 1 / 0.170
B['L2'] = 1 / 0.037
B['L3'] = 1 / 0.258
B['L4'] = 1 / 0.197
B['L5'] = 1 / 0.037
B['L6'] = 1 / 0.140
B['L7'] = 1 / 0.018

line_maxflow['L1'] = 400
line_maxflow['L2'] = 200
line_maxflow['L3'] = 200
line_maxflow['L4'] = 200
line_maxflow['L5'] = 200
line_maxflow['L6'] = 200
line_maxflow['L7'] = 200

# %%
fontsize = 10
print('NET elec load')
for t in time:
    print(t, sum(elec_load[n, t] for n in node) - wind_scenario['W1', t] * elec_maxprod['W1'])
print('NET elec load')
for t in time:
    print(t, sum(heat_load[n, t] for n in node))



for t in time:
    plt.bar(t, sum(elec_load[n, t] for n in node) - wind_scenario['W1', t] * elec_maxprod['W1'], color='b')

plt.ylabel("Электрическая нагрузка (- ветровая генерация) (МВтч)", fontsize=fontsize)
plt.xlabel("Время (ч)", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
#plt.title("Состояние аккумуляторной батареи электромобиля", fontsize=fontsize)
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()



for t in time:
    plt.bar(t, sum(heat_load[n, t] for n in node),color='b')

plt.ylabel("Тепловая нагрузка (МВтч)", fontsize=fontsize)
plt.xlabel("Время (ч)", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)
#plt.title("Состояние аккумуляторной батареи электромобиля", fontsize=fontsize)
#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


# %% building the STOCHSTICA MPEC optimization problem (GUROBI) = heat market clearing

# os.chdir("C:/Users/lemitri/Documents/phd/Combined Heat and Power Flow/python/data")


# %% building the STOCHSTICA MPEC optimization problem (GUROBI) = heat market clearing

class expando(object):
    '''


        A small class which can have attributes set
    '''
    pass


class integrated_dispatch_MILP:
    def __init__(self):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self._load_data()
        self._build_model()

    def optimize(self):
        self.model.optimize()

    def computeIIS(self):
        self.model.computeIIS()

    def _load_data(self):

        # indexes
        self.data.time = time
        self.data.time_list = time_list
        self.data.node = node
        self.data.line = line
        self.data.pipe_supply = pipe_supply
        self.data.pipe_return = pipe_return
        self.data.gen = gen
        self.data.heat_storage = heat_storage
        self.data.elec_storage = elec_storage
        self.data.heat_pump = heat_pump
        self.data.wind = wind
        self.data.heat_only = heat_only
        self.data.CHP_sorted = CHP_sorted
        self.data.CHP = CHP
        self.data.producers = producers
        self.data.heat_station = heat_station
        self.data.elec_station = elec_station
        self.data.heat_exchanger_station = heat_exchanger_station

        # epsilon
        self.data.epsilon = 0.000001

        # producers sorted per node
        self.data.producers_node = producers_node
        self.data.heat_station_node = heat_station_node
        self.data.elec_station_node = elec_station_node
        self.data.CHP_node = CHP_node
        self.data.CHP_sorted_node = CHP_sorted_node
        self.data.heat_pump_node = heat_pump_node
        self.data.heat_only_node = heat_only_node
        self.data.heat_storage_node = heat_storage_node
        self.data.gen_node = gen_node
        self.data.wind_node = wind_node
        self.data.elec_storage_node = elec_storage_node
        self.data.heat_exchanger_station_node = heat_exchanger_station_node

        # connexions between nodes
        self.data.pipe_supply_connexion = pipe_supply_connexion
        self.data.pipe_supply_start = pipe_supply_start
        self.data.pipe_supply_end = pipe_supply_end
        self.data.pipe_return_connexion = pipe_return_connexion
        self.data.pipe_return_start = pipe_return_start
        self.data.pipe_return_end = pipe_return_end
        self.data.line_connexion = line_connexion
        self.data.line_start = line_start
        self.data.line_end = line_end

        # heat network parameters
        self.data.pressure_loss_coeff = pressure_loss_coeff
        self.data.radius_pipe = radius_pipe
        self.data.water_density = water_density
        self.data.water_thermal_loss_coeff = water_thermal_loss_coeff
        self.data.water_heat_capacity = water_heat_capacity
        self.data.length_pipe = length_pipe

        # pressure bounds
        self.data.pr_return_max = pr_return_max
        self.data.pr_return_min = pr_return_min
        self.data.pr_supply_max = pr_supply_max
        self.data.pr_supply_min = pr_supply_min
        self.data.pr_diff_min = pr_supply_min
        self.data.pr_return_slice = pr_return_slice
        self.data.pr_supply_slice = pr_supply_slice

        # temperature bounds
        self.data.T_in_max = T_in_max
        self.data.T_in_min = T_in_min
        self.data.T_out_max = T_out_max
        self.data.T_out_min = T_out_min
        self.data.T_return_max = T_return_max
        self.data.T_return_min = T_return_min
        self.data.T_supply_max = T_supply_max
        self.data.T_supply_min = T_supply_min

        # mass flow bounds
        self.data.mf_pipe_max = mf_pipe_max
        self.data.mf_pipe_min = mf_pipe_min
        self.data.mf_HS_max = mf_HS_max
        self.data.mf_HS_min = mf_HS_min
        self.data.mf_HES_max = mf_HES_max
        self.data.mf_HES_min = mf_HES_min

        # LOADS
        self.data.heat_load = heat_load
        self.data.elec_load = elec_load

        # Heat station parameters
        self.data.CHP_maxprod = CHP_maxprod
        self.data.heat_maxprod = heat_maxprod
        self.data.rho_elec = rho_elec
        self.data.rho_heat = rho_heat
        self.data.r_min = r_min
        self.data.storage_rho_plus = storage_rho_plus
        self.data.storage_rho_moins = storage_rho_moins
        self.data.storage_maxcapacity = storage_maxcapacity
        self.data.storage_loss = storage_loss
        self.data.storage_init = storage_init
        self.data.COP = COP

        # Elec station parameters
        self.data.elec_maxprod = elec_maxprod
        self.data.wind_scenario = wind_scenario

        # Cost parameters
        self.data.alpha = alpha

        # time delay
        self.data.Dt = Dt
        self.data.time_delay_max = time_delay_max
        self.data.time_delay_range = time_delay_range
        self.data.mf_pipe_init = mf_pipe_init
        self.data.T_in_init = T_in_init

        self.data.time_init = time_init
        self.data.time_extended = time_extended

        self.data.tau_new_min = tau_new_min
        self.data.tau_new_max = tau_new_max

        self.data.big_M = 1 + max([self.data.mf_pipe_max[p] * self.data.Dt * (self.data.time_delay_max[p] + 1) / (
                    np.pi * self.data.radius_pipe[p] * self.data.radius_pipe[p] * self.data.water_density) -
                                   self.data.length_pipe[p] for p in pipe_supply + pipe_return])

        # elec transmission system
        self.data.B = B
        self.data.line_maxflow = line_maxflow

    def _build_model(self):

        self.model = gb.Model()
        self._build_variables()
        self._build_objective()
        self._build_constraints()

    def _build_variables(self):

        # indexes shortcuts
        time = self.data.time
        time_init = self.data.time_init
        time_extended = self.data.time_extended
        node = self.data.node
        line = self.data.line
        pipe_supply = self.data.pipe_supply
        pipe_return = self.data.pipe_return
        gen = self.data.gen
        heat_storage = self.data.heat_storage
        elec_storage = self.data.elec_storage
        heat_pump = self.data.heat_pump
        wind = self.data.wind
        heat_only = self.data.heat_only
        CHP_sorted = self.data.CHP_sorted
        CHP = self.data.CHP
        producers = self.data.producers
        heat_station = self.data.heat_station
        elec_station = self.data.elec_station
        heat_exchanger_station=self.data.heat_exchanger_station
        m = self.model

        # heat market optimization variables

        self.variables.Q = {}  # heat production from CHPs and HO units (first satge)
        for t in time:
            for h in heat_station:
                self.variables.Q[h, t] = m.addVar(lb=0, ub=self.data.heat_maxprod[h], name='Q({0},{1})'.format(h, t))

        self.variables.storage_plus = {}  # heat storage: heat discharged (first stage)
        for t in time:
            for h in heat_storage + elec_storage:
                self.variables.storage_plus[h, t] = m.addVar(lb=0, ub=self.data.storage_maxprod[h],
                                                             name='storage plus({0},{1})'.format(h, t))

        self.variables.storage_moins = {}  # heat storage: heat charged (first stage)
        for t in time:
            for h in heat_storage + elec_storage:
                self.variables.storage_moins[h, t] = m.addVar(lb=0, ub=self.data.storage_maxprod[h],
                                                              name='storage moins({0},{1})'.format(h, t))

        self.variables.storage_energy = {}  # heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_storage + elec_storage:
                self.variables.storage_energy[h, t] = m.addVar(lb=0, ub=self.data.storage_maxcapacity[h],
                                                               name='storage energy({0},{1})'.format(h, t))

                # electricity market optimization variables : primal variables

        self.variables.P = {}  # electricity production from electricity generators, CHPs and wind producers
        for t in time:
            for g in CHP + gen + wind:
                self.variables.P[g, t] = m.addVar(lb=0, ub=self.data.elec_maxprod[g],
                                                  name='P({0},{1})'.format(g, t))  # dispatch of electricity generators

            for g in heat_pump:
                self.variables.P[g, t] = m.addVar(lb=-gb.GRB.INFINITY, ub=0,
                                                  name='P({0},{1})'.format(g, t))  # elec consumption of HP

        # masse flow rates!!!

        self.variables.mf_pipe = {}
        for p in pipe_supply + pipe_return:
            for t in self.data.time_extended[p]:
                self.variables.mf_pipe[p, t] = m.addVar(lb=self.data.mf_pipe_min[p], ub=self.data.mf_pipe_max[p],
                                                        name='mf pipe({0},{1})'.format(p, t))

        self.variables.mf_HES = {}
        for t in time:
            for p in node:
                self.variables.mf_HES[p, t] = m.addVar(lb=self.data.mf_HES_min[p], ub=self.data.mf_HES_max[p],
                                                       name='mf HES({0},{1})'.format(p, t))

        self.variables.mf_HS = {}
        for t in time:
            for p in heat_station + heat_storage:
                self.variables.mf_HS[p, t] = m.addVar(lb=self.data.mf_HS_min[p], ub=self.data.mf_HS_max[p],
                                                      name='mf HS({0},{1})'.format(p, t))

                # temperatures
        self.variables.T_in = {}
        for p in pipe_supply + pipe_return:
            for t in time_extended[p]:
                self.variables.T_in[p, t] = m.addVar(lb=self.data.T_in_min[p, t], ub=self.data.T_in_max[p, t],
                                                     name='temp in({0},{1})'.format(p, t))

        self.variables.T_out = {}
        for p in pipe_supply + pipe_return:
            for t in time:
                self.variables.T_out[p, t] = m.addVar(lb=self.data.T_out_min[p, t], ub=self.data.T_out_max[p, t],
                                                      name='temp out({0},{1})'.format(p, t))

        self.variables.T_supply = {}
        for t in time:
            for n in node:
                self.variables.T_supply[n, t] = m.addVar(lb=self.data.T_supply_min[n, t],
                                                         ub=self.data.T_supply_max[n, t],
                                                         name='temp supply({0},{1})'.format(p, t))

        self.variables.T_return = {}
        for t in time:
            for n in node:
                self.variables.T_return[n, t] = m.addVar(lb=self.data.T_return_min[n, t],
                                                         ub=self.data.T_return_max[n, t],
                                                         name='temp return({0},{1})'.format(p, t))

        # pressure variables

        self.variables.pr_supply = {}
        for t in time:
            for n in node:
                self.variables.pr_supply[n, t] = m.addVar(lb=self.data.pr_supply_min[n], ub=self.data.pr_supply_max[n],
                                                          name='pr supply({0},{1})'.format(n, t))

        self.variables.pr_return = {}
        for t in time:
            for n in node:
                self.variables.pr_return[n, t] = m.addVar(lb=self.data.pr_return_min[n], ub=self.data.pr_return_max[n],
                                                          name='pr return({0},{1})'.format(n, t))

        # time delays in pipes (discrete)

        self.variables.tau = {}
        for t in time:
            for p in pipe_supply + pipe_return:
                self.variables.tau[p, t] = m.addVar(lb=0, ub=self.data.time_delay_max[p],
                                                    name='tau({0},{1})'.format(p, t))

        self.variables.u = {}
        for t in time:
            for p in pipe_supply + pipe_return:
                for n in range(self.data.time_delay_max[p] + 1):
                    self.variables.u[p, n, t] = m.addVar(vtype=gb.GRB.BINARY, name='u({0},{1},{2})'.format(p, n, t))

        self.variables.T_in_new = {}
        for p in pipe_supply + pipe_return:
            for t in time:
                for n in range(self.data.time_delay_max[p] + 1):
                    self.variables.T_in_new[p, n, t] = m.addVar(lb=-gb.GRB.INFINITY,
                                                                name='temp in new({0},{1},{2})'.format(p, n, t))

        self.variables.v = {}
        for t in time:
            for p in pipe_supply + pipe_return:
                for n in range(self.data.time_delay_max[p] + 1):
                    self.variables.v[p, n, t] = m.addVar(vtype=gb.GRB.BINARY, name='v({0},{1},{2})'.format(p, n, t))

        # electricity transmission system variables

        self.variables.node_angle = {}
        for t in time:
            for n in node:
                self.variables.node_angle[n, t] = m.addVar(lb=-gb.GRB.INFINITY, name='node angle({0},{1})'.format(n, t))

        self.variables.flow_line = {}
        for t in time:
            for l in line:
                self.variables.flow_line[l, t] = m.addVar(lb=-self.data.line_maxflow[l], ub=self.data.line_maxflow[l],
                                                          name='flow line({0},{1})'.format(l, t))

            #         #temp mixing relaxation: auxiliary variables
        #
        #        self.variables.w = {}
        #        for t in time:
        #            for p in pipe_supply+pipe_return:
        #                self.variables.w[p,t] = m.addVar(lb=-gb.GRB.INFINITY,name='w({0},{1})'.format(p,t))

        # SLACK VARIABLES (FEASIBILITY)

        self.variables.mf_supply_slack_plus = {}
        for t in time:
            for n in node:
                self.variables.mf_supply_slack_plus[n, t] = m.addVar(lb=0, ub=2000,
                                                                     name='mf slack supply plus({0},{1})'.format(n, t))

        self.variables.mf_supply_slack_moins = {}
        for t in time:
            for n in node:
                self.variables.mf_supply_slack_moins[n, t] = m.addVar(lb=0, ub=2000,
                                                                      name='mf slack supply moins({0},{1})'.format(n,
                                                                                                                   t))

        self.variables.mf_return_slack_plus = {}
        for t in time:
            for n in node:
                self.variables.mf_return_slack_plus[n, t] = m.addVar(lb=0, ub=2000,
                                                                     name='mf slack return plus({0},{1})'.format(n, t))

        self.variables.mf_return_slack_moins = {}
        for t in time:
            for n in node:
                self.variables.mf_return_slack_moins[n, t] = m.addVar(lb=0, ub=2000,
                                                                      name='mf slack return moins({0},{1})'.format(n,
                                                                                                                   t))

        m.update()

    def _build_objective(self):  # building the objective function for the heat maret clearing

        # indexes shortcuts
        time = self.data.time
        time_init = self.data.time_init
        time_extended = self.data.time_extended
        node = self.data.node
        line = self.data.line
        pipe_supply = self.data.pipe_supply
        pipe_return = self.data.pipe_return
        gen = self.data.gen
        heat_storage = self.data.heat_storage
        elec_storage = self.data.elec_storage
        heat_pump = self.data.heat_pump
        wind = self.data.wind
        heat_only = self.data.heat_only
        CHP_sorted = self.data.CHP_sorted
        CHP = self.data.CHP
        producers = self.data.producers
        heat_station = self.data.heat_station
        elec_station = self.data.elec_station
        heat_exchanger_station=self.data.heat_exchanger_station
        m = self.model

        m.setObjective(0.00001 * gb.quicksum(
            self.variables.tau[p, t] for t in time for p in pipe_supply + pipe_return) + 1000 * gb.quicksum(
            self.variables.mf_supply_slack_moins[n, t] + self.variables.mf_supply_slack_plus[n, t] +
            self.variables.mf_return_slack_moins[n, t] + self.variables.mf_return_slack_plus[n, t] for n in node for t
            in time) + gb.quicksum(
            self.data.alpha[g, t] * self.variables.P[g, t] for t in time for g in gen + wind) + gb.quicksum(
            self.data.alpha[g, t] * self.variables.Q[g, t] for t in time for g in heat_only) + gb.quicksum(
            self.data.alpha[g, t] * (
                        self.data.rho_elec[g] * self.variables.P[g, t] + self.data.rho_heat[g] * self.variables.Q[g, t])
            for t in time for g in CHP),
                       gb.GRB.MINIMIZE)

    def _build_constraints(self):

        # indexes shortcuts
        time = self.data.time
        time_init = self.data.time_init
        time_extended = self.data.time_extended
        node = self.data.node
        line = self.data.line
        pipe_supply = self.data.pipe_supply
        pipe_return = self.data.pipe_return
        gen = self.data.gen
        heat_storage = self.data.heat_storage
        elec_storage = self.data.elec_storage
        heat_pump = self.data.heat_pump
        wind = self.data.wind
        heat_only = self.data.heat_only
        CHP_sorted = self.data.CHP_sorted
        CHP = self.data.CHP
        producers = self.data.producers
        heat_station = self.data.heat_station
        elec_station = self.data.elec_station
        heat_exchanger_station=self.data.heat_exchanger_station
        m = self.model

        # HES equations

        #        self.constraints.heat_load = {}
        #
        #        for t in time:
        #            for n in node:
        #                self.constraints.heat_load[n,t] = m.addConstr(
        #                        10**(-6)*self.data.water_heat_capacity*self.data.Dt*self.variables.mf_HES[n,t]*(self.variables.T_supply[n,t]-self.variables.T_return[n,t]),
        #                        gb.GRB.EQUAL,
        #                        self.data.heat_load[n,t])

        self.constraints.heat_load_1 = {}

        for t in time:
            for n in node:
                self.constraints.heat_load_1[n, t] = m.addConstr(
                    self.data.heat_load[n, t],
                    gb.GRB.GREATER_EQUAL,
                    10 ** (-6) * self.data.water_heat_capacity * self.data.Dt * (self.data.mf_HES_max[n] * (
                                self.variables.T_supply[n, t] - self.variables.T_return[n, t]) + self.variables.mf_HES[
                                                                                     n, t] * (self.data.T_supply_max[
                                                                                                  n, t] -
                                                                                              self.data.T_return_min[
                                                                                                  n, t]) -
                                                                                 self.data.mf_HES_max[n] * (
                                                                                             self.data.T_supply_max[
                                                                                                 n, t] -
                                                                                             self.data.T_return_min[
                                                                                                 n, t])),
                    name='heat load 1({0},{1})'.format(n, t))

        self.constraints.heat_load_2 = {}

        for t in time:
            for n in node:
                self.constraints.heat_load_2[n, t] = m.addConstr(
                    self.data.heat_load[n, t],
                    gb.GRB.GREATER_EQUAL,
                    10 ** (-6) * self.data.water_heat_capacity * self.data.Dt * (self.data.mf_HES_min[n] * (
                                self.variables.T_supply[n, t] - self.variables.T_return[n, t]) + self.variables.mf_HES[
                                                                                     n, t] * (self.data.T_supply_min[
                                                                                                  n, t] -
                                                                                              self.data.T_return_max[
                                                                                                  n, t]) -
                                                                                 self.data.mf_HES_min[n] * (
                                                                                             self.data.T_supply_min[
                                                                                                 n, t] -
                                                                                             self.data.T_return_max[
                                                                                                 n, t])),
                    name='heat load 2({0},{1})'.format(n, t))

        self.constraints.heat_load_3 = {}

        for t in time:
            for n in node:
                self.constraints.heat_load_3[n, t] = m.addConstr(
                    self.data.heat_load[n, t],
                    gb.GRB.LESS_EQUAL,
                    10 ** (-6) * self.data.water_heat_capacity * self.data.Dt * (self.data.mf_HES_max[n] * (
                                self.variables.T_supply[n, t] - self.variables.T_return[n, t]) + self.variables.mf_HES[
                                                                                     n, t] * (self.data.T_supply_min[
                                                                                                  n, t] -
                                                                                              self.data.T_return_max[
                                                                                                  n, t]) -
                                                                                 self.data.mf_HES_max[n] * (
                                                                                             self.data.T_supply_min[
                                                                                                 n, t] -
                                                                                             self.data.T_return_max[
                                                                                                 n, t])),
                    name='heat load 3({0},{1})'.format(n, t))

        self.constraints.heat_load_4 = {}

        for t in time:
            for n in node:
                self.constraints.heat_load_4[n, t] = m.addConstr(
                    self.data.heat_load[n, t],
                    gb.GRB.LESS_EQUAL,
                    10 ** (-6) * self.data.water_heat_capacity * self.data.Dt * (self.data.mf_HES_min[n] * (
                                self.variables.T_supply[n, t] - self.variables.T_return[n, t]) + self.variables.mf_HES[
                                                                                     n, t] * (self.data.T_supply_max[
                                                                                                  n, t] -
                                                                                              self.data.T_return_min[
                                                                                                  n, t]) -
                                                                                 self.data.mf_HES_min[n] * (
                                                                                             self.data.T_supply_max[
                                                                                                 n, t] -
                                                                                             self.data.T_return_min[
                                                                                                 n, t])),
                    name='heat load 4({0},{1})'.format(n, t))

        self.constraints.pressure_diff = {}
        for t in time:
            for n in node:
                self.constraints.pressure_diff[n, t] = m.addConstr(
                    self.variables.pr_supply[n, t] - self.variables.pr_return[n, t],
                    gb.GRB.GREATER_EQUAL,
                    self.data.pr_diff_min[n], name='pressure diff({0},{1})'.format(n, t))

        # HS equations

        #        self.constraints.heat_station = {}
        #
        #        for t in time:
        #            for n in node:
        #
        #                for h in self.data.heat_station_node[n]:
        #                    self.constraints.heat_station[h,t] = m.addConstr(
        #                            10**(-6)*self.data.water_heat_capacity*self.data.Dt*self.variables.mf_HS[h,t]*(self.variables.T_supply[n,t]-self.variables.T_return[n,t]),
        #                            gb.GRB.EQUAL,
        #                            self.variables.Q[h,t])
        #
        #                for h in self.data.heat_storage_node[n]:
        #                    self.constraints.heat_station[h,t] = m.addConstr(
        #                            10**(-6)*self.data.water_heat_capacity*self.data.Dt*self.variables.mf_HS[h,t]*(self.variables.T_supply[n,t]-self.variables.T_return[n,t]),
        #                            gb.GRB.EQUAL,
        #                            self.variables.storage_plus[h,t]-self.variables.storage_moins[h,t])

        self.constraints.heat_station_1 = {}

        for t in time:
            for n in node:
                for h in self.data.heat_station_node[n]:
                    self.constraints.heat_station_1[h, t] = m.addConstr(
                        self.variables.Q[h, t],
                        gb.GRB.GREATER_EQUAL,
                        10 ** (-6) * self.data.water_heat_capacity * self.data.Dt * (self.data.mf_HS_max[h] * (
                                    self.variables.T_supply[n, t] - self.variables.T_return[n, t]) +
                                                                                     self.variables.mf_HS[h, t] * (
                                                                                                 self.data.T_supply_max[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_min[
                                                                                                     n, t]) -
                                                                                     self.data.mf_HS_max[h] * (
                                                                                                 self.data.T_supply_max[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_min[
                                                                                                     n, t])),
                        name='heat station 1({0},{1})'.format(h, t))

                for h in self.data.heat_storage_node[n]:
                    self.constraints.heat_station_1[h, t] = m.addConstr(
                        self.variables.storage_plus[h, t] - self.variables.storage_moins[h, t],
                        gb.GRB.GREATER_EQUAL,
                        10 ** (-6) * self.data.water_heat_capacity * self.data.Dt * (self.data.mf_HS_max[h] * (
                                    self.variables.T_supply[n, t] - self.variables.T_return[n, t]) +
                                                                                     self.variables.mf_HS[h, t] * (
                                                                                                 self.data.T_supply_max[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_min[
                                                                                                     n, t]) -
                                                                                     self.data.mf_HS_max[h] * (
                                                                                                 self.data.T_supply_max[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_min[
                                                                                                     n, t])),
                        name='heat station 1({0},{1})'.format(h, t))

        self.constraints.heat_station_2 = {}

        for t in time:
            for n in node:
                for h in self.data.heat_station_node[n]:
                    self.constraints.heat_station_2[h, t] = m.addConstr(
                        self.variables.Q[h, t],
                        gb.GRB.GREATER_EQUAL,
                        10 ** (-6) * self.data.water_heat_capacity * self.data.Dt * (self.data.mf_HS_min[h] * (
                                    self.variables.T_supply[n, t] - self.variables.T_return[n, t]) +
                                                                                     self.variables.mf_HS[h, t] * (
                                                                                                 self.data.T_supply_min[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_max[
                                                                                                     n, t]) -
                                                                                     self.data.mf_HS_min[h] * (
                                                                                                 self.data.T_supply_min[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_max[
                                                                                                     n, t])),
                        name='heat station 2({0},{1})'.format(h, t))

                for h in self.data.heat_storage_node[n]:
                    self.constraints.heat_station_2[h, t] = m.addConstr(
                        self.variables.storage_plus[h, t] - self.variables.storage_moins[h, t],
                        gb.GRB.GREATER_EQUAL,
                        10 ** (-6) * self.data.water_heat_capacity * self.data.Dt * (self.data.mf_HS_min[h] * (
                                    self.variables.T_supply[n, t] - self.variables.T_return[n, t]) +
                                                                                     self.variables.mf_HS[h, t] * (
                                                                                                 self.data.T_supply_min[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_max[
                                                                                                     n, t]) -
                                                                                     self.data.mf_HS_min[h] * (
                                                                                                 self.data.T_supply_min[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_max[
                                                                                                     n, t])),
                        name='heat station 2({0},{1})'.format(h, t))

        self.constraints.heat_station_3 = {}

        for t in time:
            for n in node:
                for h in self.data.heat_station_node[n]:
                    self.constraints.heat_station_3[h, t] = m.addConstr(
                        self.variables.Q[h, t],
                        gb.GRB.LESS_EQUAL,
                        10 ** (-6) * self.data.water_heat_capacity * self.data.Dt * (self.data.mf_HS_max[h] * (
                                    self.variables.T_supply[n, t] - self.variables.T_return[n, t]) +
                                                                                     self.variables.mf_HS[h, t] * (
                                                                                                 self.data.T_supply_min[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_max[
                                                                                                     n, t]) -
                                                                                     self.data.mf_HS_max[h] * (
                                                                                                 self.data.T_supply_min[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_max[
                                                                                                     n, t])),
                        name='heat station 3({0},{1})'.format(h, t))

                for h in self.data.heat_storage_node[n]:
                    self.constraints.heat_station_3[h, t] = m.addConstr(
                        self.variables.storage_plus[h, t] - self.variables.storage_moins[h, t],
                        gb.GRB.LESS_EQUAL,
                        10 ** (-6) * self.data.water_heat_capacity * self.data.Dt * (self.data.mf_HS_max[h] * (
                                    self.variables.T_supply[n, t] - self.variables.T_return[n, t]) +
                                                                                     self.variables.mf_HS[h, t] * (
                                                                                                 self.data.T_supply_min[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_max[
                                                                                                     n, t]) -
                                                                                     self.data.mf_HS_max[h] * (
                                                                                                 self.data.T_supply_min[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_max[
                                                                                                     n, t])),
                        name='heat station 3({0},{1})'.format(h, t))

        self.constraints.heat_station_4 = {}

        for t in time:
            for n in node:
                for h in self.data.heat_station_node[n]:
                    self.constraints.heat_station_4[h, t] = m.addConstr(
                        self.variables.Q[h, t],
                        gb.GRB.LESS_EQUAL,
                        10 ** (-6) * self.data.water_heat_capacity * self.data.Dt * (self.data.mf_HS_min[h] * (
                                    self.variables.T_supply[n, t] - self.variables.T_return[n, t]) +
                                                                                     self.variables.mf_HS[h, t] * (
                                                                                                 self.data.T_supply_max[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_min[
                                                                                                     n, t]) -
                                                                                     self.data.mf_HS_min[h] * (
                                                                                                 self.data.T_supply_max[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_min[
                                                                                                     n, t])),
                        name='heat station 4({0},{1})'.format(h, t))

                for h in self.data.heat_storage_node[n]:
                    self.constraints.heat_station_4[h, t] = m.addConstr(
                        self.variables.storage_plus[h, t] - self.variables.storage_moins[h, t],
                        gb.GRB.LESS_EQUAL,
                        10 ** (-6) * self.data.water_heat_capacity * self.data.Dt * (self.data.mf_HS_min[h] * (
                                    self.variables.T_supply[n, t] - self.variables.T_return[n, t]) +
                                                                                     self.variables.mf_HS[h, t] * (
                                                                                                 self.data.T_supply_max[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_min[
                                                                                                     n, t]) -
                                                                                     self.data.mf_HS_min[h] * (
                                                                                                 self.data.T_supply_max[
                                                                                                     n, t] -
                                                                                                 self.data.T_return_min[
                                                                                                     n, t])),
                        name='heat station 4({0},{1})'.format(h, t))

                    ##1) CHPs

        self.constraints.CHP_maxprod = {}
        self.constraints.CHP_ratio = {}

        for t in time:
            for h in self.data.CHP_sorted['ex']:
                self.constraints.CHP_maxprod[h, t] = m.addConstr(
                    self.data.rho_heat[h] * self.variables.Q[h, t] + self.data.rho_elec[h] * self.variables.P[h, t],
                    gb.GRB.LESS_EQUAL,
                    self.data.CHP_maxprod[h], name='CHP maxprod({0},{1})'.format(h, t))

                self.constraints.CHP_ratio[h, t] = m.addConstr(
                    self.variables.P[h, t],
                    gb.GRB.GREATER_EQUAL,
                    self.data.r_min[h] * self.variables.Q[h, t], name='CHP ratio({0},{1})'.format(h, t))

            for h in self.data.CHP_sorted['bp']:
                self.constraints.CHP_ratio[h, t] = m.addConstr(
                    self.variables.P[h, t],
                    gb.GRB.EQUAL,
                    self.data.r_min[h] * self.variables.Q[h, t], name='CHP ratio({0},{1})'.format(h, t))

                ##2) storage

        self.constraints.storage_update = {}

        for (t1, t2) in zip(time[:-1], time[1:]):
            for h in heat_storage + elec_storage:
                self.constraints.storage_update[h, t2] = m.addConstr(
                    self.variables.storage_energy[h, t2],
                    gb.GRB.EQUAL,
                    self.variables.storage_energy[h, t1] - self.data.storage_rho_plus[h] * self.variables.storage_plus[
                        h, t2] + self.data.storage_rho_moins[h] * self.variables.storage_moins[h, t2] -
                    self.data.storage_loss[h], name='storage update({0},{1})'.format(h, t2))

        self.constraints.storage_init = {}

        for h in heat_storage + elec_storage:
            self.constraints.storage_init[h] = m.addConstr(
                self.variables.storage_energy[h, time[0]],
                gb.GRB.EQUAL,
                self.data.storage_init[h] - self.data.storage_rho_plus[h] * self.variables.storage_plus[h, time[0]] +
                self.data.storage_rho_moins[h] * self.variables.storage_moins[h, time[0]] - self.data.storage_loss[h],
                name='storage init({0})'.format(h))

        self.constraints.storage_final = {}

        for h in heat_storage + elec_storage:
            self.constraints.storage_final[h] = m.addConstr(
                self.variables.storage_energy[h, time[-1]],
                gb.GRB.GREATER_EQUAL,
                self.data.storage_init[h], name='storage final({0})'.format(h))

        # 3) heat pumps

        self.constraints.heat_pump_COP = {}

        for t in time:
            for h in heat_pump:
                self.constraints.heat_pump_COP[h, t] = m.addConstr(
                    -self.data.COP[h] * self.variables.P[h, t],
                    gb.GRB.EQUAL,
                    self.variables.Q[h, t], name='heat pump COP({0},{1})'.format(h, t))

                # continuity of mass flow at nodes

        self.constraints.mass_flow_continuity_supply = {}
        self.constraints.mass_flow_continuity_return = {}

        for t in time:
            for n in node:
                self.constraints.mass_flow_continuity_supply[n, t] = m.addConstr(
                    gb.quicksum(self.variables.mf_HS[h, t] for h in self.data.heat_station_node[n]) + gb.quicksum(
                        self.variables.mf_HS[h, t] for h in self.data.heat_storage_node[n]) + gb.quicksum(
                        self.variables.mf_pipe[p, t] for p in self.data.pipe_supply_end[n]) - gb.quicksum(
                        self.variables.mf_pipe[p, t] for p in self.data.pipe_supply_start[n]) + (
                                self.variables.mf_supply_slack_plus[n, t] - self.variables.mf_supply_slack_moins[n, t]),
                    gb.GRB.EQUAL,
                    self.variables.mf_HES[n, t], name='mass flow continuity supply({0},{1})'.format(n, t))

                self.constraints.mass_flow_continuity_return[n, t] = m.addConstr(
                    self.variables.mf_HES[n, t] + gb.quicksum(
                        self.variables.mf_pipe[p, t] for p in self.data.pipe_return_end[n]) - gb.quicksum(
                        self.variables.mf_pipe[p, t] for p in self.data.pipe_return_start[n]) + (
                                self.variables.mf_return_slack_plus[n, t] - self.variables.mf_return_slack_moins[n, t]),
                    gb.GRB.EQUAL,
                    gb.quicksum(self.variables.mf_HS[h, t] for h in self.data.heat_station_node[n]) + gb.quicksum(
                        self.variables.mf_HS[h, t] for h in self.data.heat_storage_node[n]),
                    name='mass flow continuity return({0},{1})'.format(n, t))

        # temperature mixing at nodes

        # SIMPLIFIED VERSION (ONLY 1 pipe arriving!!!!!!!!!!!!!!!)

        self.constraints.T_mixing_supply = {}
        self.constraints.T_mixing_return = {}

        for t in time:
            for n in node:
                for p in self.data.pipe_supply_end[n]:
                    self.constraints.T_mixing_supply[p, t] = m.addConstr(
                        self.variables.T_supply[n, t],
                        gb.GRB.EQUAL,
                        self.variables.T_out[p, t], name='Temp mixing supply simplified ({0},{1},{2})'.format(n, p, t))

                for p in self.data.pipe_return_end[n]:
                    self.constraints.T_mixing_return[p, t] = m.addConstr(
                        self.variables.T_return[n, t],
                        gb.GRB.EQUAL,
                        self.variables.T_out[p, t], name='Temp mixing return simplified ({0},{1},{2})'.format(n, p, t))

                ##        self.constraints.T_mixing_supply = {}
        ##        self.constraints.T_mixing_return = {}
        ##
        ##        for t in time:
        ##                for n in node:
        ##
        ##                    self.constraints.T_mixing_supply[n,t] = m.addConstr(
        ##                        self.variables.T_supply[n,t]*gb.quicksum(self.variables.mf_pipe[p,t] for p in self.data.pipe_supply_end[n]),
        ##                        gb.GRB.EQUAL,
        ##                        gb.quicksum(self.variables.T_out[p,t]*self.variables.mf_pipe[p,t] for p in self.data.pipe_supply_end[n]))
        ##
        ##                    self.constraints.T_mixing_return[n,t] = m.addConstr(
        ##                        self.variables.T_return[n,t]*gb.quicksum(self.variables.mf_pipe[p,t] for p in self.data.pipe_return_end[n]),
        ##                        gb.GRB.EQUAL,
        ##                        gb.quicksum(self.variables.T_out[p,t]*self.variables.mf_pipe[p,t] for p in self.data.pipe_return_end[n]))
        #
        #
        #        self.constraints.T_mixing_supply_0 = {}
        #        self.constraints.T_mixing_return_0 = {}
        #
        #        for t in time:
        #                for n in node:
        #
        #                    self.constraints.T_mixing_supply_0[n,t] = m.addConstr(
        #                        gb.quicksum(self.variables.w[p,t] for p in self.data.pipe_supply_end[n]),
        #                        gb.GRB.EQUAL,
        #                        0,name='T mixing supply 0({0},{1})'.format(n,t))
        #
        #                    self.constraints.T_mixing_return_0[n,t] = m.addConstr(
        #                        gb.quicksum(self.variables.w[p,t] for p in self.data.pipe_return_end[n]),
        #                        gb.GRB.EQUAL,
        #                        0,name='T mixing return 0({0},{1})'.format(n,t))
        #
        ##        self.constraints.T_mixing_supply_1 = {}
        ##        self.constraints.T_mixing_return_1 = {}
        ##
        ##        for t in time:
        ##                for n in node:
        ##                    for p in pipe_supply_end[n]:
        ##                        self.constraints.T_mixing_supply_1[p,t] = m.addConstr(
        ##                            self.variables.w[p,t],
        ##                            gb.GRB.EQUAL,
        ##                            self.variables.mf_pipe[p,t]*(self.variables.T_supply[n,t]-self.variables.T_out[p,t]))
        ##
        ##                    for p in pipe_return_end[n]:
        ##                        self.constraints.T_mixing_return_1[n,t] = m.addConstr(
        ##                            self.variables.w[p,t],
        ##                            gb.GRB.EQUAL,
        ##                            self.variables.mf_pipe[p,t]*(self.variables.T_return[n,t]-self.variables.T_out[p,t]))
        #
        #        self.constraints.T_mixing_supply_1 = {}
        #        self.constraints.T_mixing_return_1 = {}
        #
        #        for t in time:
        #                for n in node:
        #                    for p in self.data.pipe_supply_end[n]:
        #                        self.constraints.T_mixing_supply_1[p,t] = m.addConstr(
        #                            self.variables.w[p,t],
        #                            gb.GRB.GREATER_EQUAL,
        #                            self.data.mf_pipe_max[p]*(self.variables.T_supply[n,t]-self.variables.T_out[p,t])+self.variables.mf_pipe[p,t]*(self.data.T_supply_max[n,t]-self.data.T_out_min[p,t])-self.data.mf_pipe_max[p]*(self.data.T_supply_max[n,t]-self.data.T_out_min[p,t]),name='T mixing supply 1({0},{1})'.format(p,t))
        #
        #                    for p in self.data.pipe_return_end[n]:
        #                        self.constraints.T_mixing_return_1[p,t] = m.addConstr(
        #                            self.variables.w[p,t],
        #                            gb.GRB.GREATER_EQUAL,
        #                            self.data.mf_pipe_max[p]*(self.variables.T_return[n,t]-self.variables.T_out[p,t])+self.variables.mf_pipe[p,t]*(self.data.T_return_max[n,t]-self.data.T_out_min[p,t])-self.data.mf_pipe_max[p]*(self.data.T_return_max[n,t]-self.data.T_out_min[p,t]),name='T mixing return 1({0},{1})'.format(p,t))
        #
        #        self.constraints.T_mixing_supply_2 = {}
        #        self.constraints.T_mixing_return_2 = {}
        #
        #        for t in time:
        #                for n in node:
        #                    for p in self.data.pipe_supply_end[n]:
        #                        self.constraints.T_mixing_supply_2[p,t] = m.addConstr(
        #                            self.variables.w[p,t],
        #                            gb.GRB.GREATER_EQUAL,
        #                            self.data.mf_pipe_min[p]*(self.variables.T_supply[n,t]-self.variables.T_out[p,t])+self.variables.mf_pipe[p,t]*(self.data.T_supply_min[n,t]-self.data.T_out_max[p,t])-self.data.mf_pipe_min[p]*(self.data.T_supply_min[n,t]-self.data.T_out_max[p,t]),name='T mixing supply 2({0},{1})'.format(p,t))
        #
        #                    for p in self.data.pipe_return_end[n]:
        #                        self.constraints.T_mixing_return_2[p,t] = m.addConstr(
        #                            self.variables.w[p,t],
        #                            gb.GRB.GREATER_EQUAL,
        #                            self.data.mf_pipe_min[p]*(self.variables.T_return[n,t]-self.variables.T_out[p,t])+self.variables.mf_pipe[p,t]*(self.data.T_return_min[n,t]-self.data.T_out_max[p,t])-self.data.mf_pipe_min[p]*(self.data.T_return_min[n,t]-self.data.T_out_max[p,t]),name='T mixing return 2({0},{1})'.format(p,t))
        #
        #        self.constraints.T_mixing_supply_3 = {}
        #        self.constraints.T_mixing_return_3 = {}
        #
        #        for t in time:
        #                for n in node:
        #                    for p in self.data.pipe_supply_end[n]:
        #                        self.constraints.T_mixing_supply_3[p,t] = m.addConstr(
        #                            self.variables.w[p,t],
        #                            gb.GRB.LESS_EQUAL,
        #                            self.data.mf_pipe_max[p]*(self.variables.T_supply[n,t]-self.variables.T_out[p,t])+self.variables.mf_pipe[p,t]*(self.data.T_supply_min[n,t]-self.data.T_out_max[p,t])-self.data.mf_pipe_max[p]*(self.data.T_supply_min[n,t]-self.data.T_out_max[p,t]),name='T mixing supply 3({0},{1})'.format(p,t))
        #
        #                    for p in self.data.pipe_return_end[n]:
        #                        self.constraints.T_mixing_return_3[p,t] = m.addConstr(
        #                            self.variables.w[p,t],
        #                            gb.GRB.LESS_EQUAL,
        #                            self.data.mf_pipe_max[p]*(self.variables.T_return[n,t]-self.variables.T_out[p,t])+self.variables.mf_pipe[p,t]*(self.data.T_return_min[n,t]-self.data.T_out_max[p,t])-self.data.mf_pipe_max[p]*(self.data.T_return_min[n,t]-self.data.T_out_max[p,t]),name='T mixing return 3({0},{1})'.format(p,t))
        #
        #        self.constraints.T_mixing_supply_4 = {}
        #        self.constraints.T_mixing_return_4 = {}
        #
        #        for t in time:
        #                for n in node:
        #                    for p in self.data.pipe_supply_end[n]:
        #                        self.constraints.T_mixing_supply_4[p,t] = m.addConstr(
        #                            self.variables.w[p,t],
        #                            gb.GRB.LESS_EQUAL,
        #                            self.data.mf_pipe_min[p]*(self.variables.T_supply[n,t]-self.variables.T_out[p,t])+self.variables.mf_pipe[p,t]*(self.data.T_supply_max[n,t]-self.data.T_out_min[p,t])-self.data.mf_pipe_min[p]*(self.data.T_supply_max[n,t]-self.data.T_out_min[p,t]),name='T mixing supply 4({0},{1})'.format(p,t))
        #
        #                    for p in self.data.pipe_return_end[n]:
        #                        self.constraints.T_mixing_return_4[p,t] = m.addConstr(
        #                            self.variables.w[p,t],
        #                            gb.GRB.LESS_EQUAL,
        #                            self.data.mf_pipe_min[p]*(self.variables.T_return[n,t]-self.variables.T_out[p,t])+self.variables.mf_pipe[p,t]*(self.data.T_return_max[n,t]-self.data.T_out_min[p,t])-self.data.mf_pipe_min[p]*(self.data.T_return_max[n,t]-self.data.T_out_min[p,t]),name='T mixing return 4({0},{1})'.format(p,t))

        self.constraints.T_in_supply = {}
        self.constraints.T_in_return = {}

        for t in time:
            for n in node:

                for p in self.data.pipe_supply_start[n]:
                    self.constraints.T_in_supply[p, t] = m.addConstr(
                        self.variables.T_supply[n, t],
                        gb.GRB.EQUAL,
                        self.variables.T_in[p, t], name='T in supply({0},{1})'.format(p, t))

                for p in self.data.pipe_return_start[n]:
                    self.constraints.T_in_return[p, t] = m.addConstr(
                        self.variables.T_return[n, t],
                        gb.GRB.EQUAL,
                        self.variables.T_in[p, t], name='T in return({0},{1})'.format(p, t))

                #        #temp dynamics in pipes

        self.constraints.T_in_init = {}

        for p in pipe_supply + pipe_return:
            for t in time_init[p]:
                self.constraints.T_in_init[p, t] = m.addConstr(
                    self.variables.T_in[p, t],
                    gb.GRB.EQUAL,
                    self.data.T_in_init[p, t], name='T in init({0},{1})'.format(p, t))

        self.constraints.mf_pipe_init = {}
        for p in pipe_supply + pipe_return:
            for t in time_init[p]:
                self.constraints.mf_pipe_init[p, t] = m.addConstr(
                    self.variables.mf_pipe[p, t],
                    gb.GRB.EQUAL,
                    self.data.mf_pipe_init[p, t], name='mf pipe init({0},{1})'.format(p, t))

        self.constraints.time_delay_discretization = {}

        for t in time:
            for p in pipe_supply + pipe_return:
                self.constraints.time_delay_discretization[p, t] = m.addConstr(
                    self.variables.tau[p, t],
                    gb.GRB.EQUAL,
                    gb.quicksum(n1 * (self.variables.u[p, n1, t] - self.variables.u[p, n2, t]) for (n1, n2) in
                                zip(self.data.time_delay_range[p][1:], self.data.time_delay_range[p][:-1])),
                    name='time delay descretization({0},{1})'.format(p, t))

        self.constraints.time_delay_u1 = {}
        self.constraints.time_delay_u2 = {}

        for p in pipe_supply + pipe_return:
            for t in self.data.time_list:
                for n in self.data.time_delay_range[p]:
                    self.constraints.time_delay_u1[p, n, time[t]] = m.addConstr(
                        -self.data.big_M * (1 - self.variables.u[p, n, time[t]]),
                        gb.GRB.LESS_EQUAL,
                        gb.quicksum(self.variables.mf_pipe[p, k] for k in time_extended[p][
                                                                          t + self.data.time_delay_max[p] - n:t +
                                                                                                              self.data.time_delay_max[
                                                                                                                  p] + 1]) * self.data.Dt / (
                                    np.pi * self.data.radius_pipe[p] * self.data.radius_pipe[
                                p] * self.data.water_density) - self.data.length_pipe[p],
                        name='time delay u1({0},{1})'.format(p, n, time[t]))

                    self.constraints.time_delay_u2[p, n, time[t]] = m.addConstr(
                        self.data.big_M * self.variables.u[p, n, time[t]],
                        gb.GRB.GREATER_EQUAL,
                        self.data.epsilon + gb.quicksum(self.variables.mf_pipe[p, k] for k in time_extended[p][t +
                                                                                                               self.data.time_delay_max[
                                                                                                                   p] - n:t +
                                                                                                                          self.data.time_delay_max[
                                                                                                                              p] + 1]) * self.data.Dt / (
                                    np.pi * self.data.radius_pipe[p] * self.data.radius_pipe[
                                p] * self.data.water_density) - self.data.length_pipe[p],
                        name='time delay u2({0},{1})'.format(p, n, time[t]))

        self.constraints.time_delay_v1 = {}
        self.constraints.time_delay_v2 = {}
        self.constraints.time_delay_v3 = {}
        self.constraints.time_delay_v4 = {}
        self.constraints.time_delay_v5 = {}
        self.constraints.time_delay_v6 = {}
        self.constraints.time_delay_v7 = {}

        for p in pipe_supply + pipe_return:
            for t in self.data.time_list:
                for n in self.data.time_delay_range[p]:
                    # NEW

                    self.constraints.time_delay_v1[p, n, time[t]] = m.addConstr(
                        -self.data.big_M * self.variables.v[p, n, time[t]],
                        gb.GRB.LESS_EQUAL,
                        self.variables.T_in_new[p, n, time[t]], name='time delay v1({0},{1})'.format(p, n, time[t]))

                    self.constraints.time_delay_v2[p, n, time[t]] = m.addConstr(
                        self.variables.T_in_new[p, n, time[t]],
                        gb.GRB.LESS_EQUAL,
                        self.data.big_M * self.variables.v[p, n, time[t]],
                        name='time delay v2({0},{1})'.format(p, n, time[t]))

                    # NEW

                    self.constraints.time_delay_v3[p, n, time[t]] = m.addConstr(
                        -self.data.big_M * (1 - self.variables.v[p, n, time[t]]),
                        gb.GRB.LESS_EQUAL,
                        self.variables.T_in_new[p, n, time[t]] - self.variables.T_in[
                            p, time_extended[p][t + self.data.time_delay_max[p] - n]] * (
                                    1 - 2 * n * self.data.water_thermal_loss_coeff / (
                                        self.data.water_heat_capacity * self.data.water_density * self.data.radius_pipe[
                                    p])), name='time delay v3({0},{1})'.format(p, n, time[t]))

                    self.constraints.time_delay_v4[p, n, time[t]] = m.addConstr(
                        self.variables.T_in_new[p, n, time[t]] - self.variables.T_in[
                            p, time_extended[p][t + self.data.time_delay_max[p] - n]] * (
                                    1 - 2 * n * self.data.water_thermal_loss_coeff / (
                                        self.data.water_heat_capacity * self.data.water_density * self.data.radius_pipe[
                                    p])),
                        gb.GRB.LESS_EQUAL,
                        self.data.big_M * (1 - self.variables.v[p, n, time[t]]),
                        name='time delay v4({0},{1})'.format(p, n, time[t]))

        #                        # WITHOUT LOSSES!!!!!!!!!!!!
        #                        self.constraints.time_delay_v3[p,n,time[t]] = m.addConstr(
        #                            -self.data.big_M*(1-self.variables.v[p,n,time[t]]),
        #                            gb.GRB.LESS_EQUAL,
        #                            self.variables.T_in_new[p,n,time[t]]-self.variables.T_in[p,time_extended[p][t+self.data.time_delay_max[p]-n]],name='time delay v3({0},{1})'.format(p,n,time[t]))
        #
        #                        self.constraints.time_delay_v4[p,n,time[t]] = m.addConstr(
        #                            self.variables.T_in_new[p,n,time[t]]-self.variables.T_in[p,time_extended[p][t+self.data.time_delay_max[p]-n]],
        #                            gb.GRB.LESS_EQUAL,
        #                            self.data.big_M*(1-self.variables.v[p,n,time[t]]),name='time delay v4({0},{1})'.format(p,n,time[t]))
        #
        #
        #                        self.constraints.time_delay_v5[p,n,time[t]] = m.addConstr(
        #                            -self.data.big_M*(1-self.variables.v[p,n,time[t]]),
        #                            gb.GRB.LESS_EQUAL,
        #                            n-self.variables.tau[p,time[t]],name='time delay v5({0},{1})'.format(p,n,time[t]))
        #
        #                        self.constraints.time_delay_v6[p,n,time[t]] = m.addConstr(
        #                            n-self.variables.tau[p,time[t]],
        #                            gb.GRB.LESS_EQUAL,
        #                            self.data.big_M*(1-self.variables.v[p,n,time[t]]),name='time delay v6({0},{1})'.format(p,n,time[t]))

        for p in pipe_supply + pipe_return:
            for t in self.data.time_list:
                self.constraints.time_delay_v7[p, time[t]] = m.addConstr(
                    gb.quicksum(self.variables.v[p, n, time[t]] for n in self.data.time_delay_range[p]),
                    gb.GRB.EQUAL,
                    1, name='time delay v7({0},{1})'.format(p, time[t]))

        #        self.constraints.temp_out = {}
        #
        #        for p in pipe_supply+pipe_return:
        #            for t in self.data.time_list:
        #
        #                    self.constraints.temp_out[p,time[t]] = m.addConstr(
        #                        gb.quicksum(self.variables.T_in[p,time_extended[p][t+self.data.time_delay_max[p]-n]]*self.variables.tau_new[p,n,time[t]] for n in self.data.time_delay_range[p]),
        #                        gb.GRB.EQUAL,
        #                        self.variables.T_out[p,time[t]])

        self.constraints.temp_out = {}

        for p in pipe_supply + pipe_return:
            for t in time:
                self.constraints.temp_out[p, t] = m.addConstr(
                    gb.quicksum(self.variables.T_in_new[p, k, t] for k in self.data.time_delay_range[p]),
                    gb.GRB.EQUAL,
                    self.variables.T_out[p, t], name='temp out 0({0},{1})'.format(p, t))

        # pressure loss in pipes LP!!!!!!!!! CONVEX

        #        self.constraints.pressure_loss_supply = {}
        #
        #        for t in time:
        #                for p in pipe_supply:
        #
        #                    self.constraints.pressure_loss_supply[p,t] = m.addConstr(
        #                        self.variables.pr_supply[self.data.pipe_supply_connexion[p][0],t]-self.variables.pr_supply[self.data.pipe_supply_connexion[p][1],t],
        #                        gb.GRB.GREATER_EQUAL,
        #                        self.data.pressure_loss_coeff[p]*self.variables.mf_pipe[p,t]*self.variables.mf_pipe[p,t],name='pressure loss supply({0},{1})'.format(p,t))
        #
        #        self.constraints.pressure_loss_return = {}
        #
        #        for t in time:
        #                for p in pipe_return:
        #
        #                    self.constraints.pressure_loss_return[p,t] = m.addConstr(
        #                        self.variables.pr_return[self.data.pipe_return_connexion[p][0],t]-self.variables.pr_return[self.data.pipe_return_connexion[p][1],t],
        #                        gb.GRB.GREATER_EQUAL,
        #                        self.data.pressure_loss_coeff[p]*self.variables.mf_pipe[p,t]*self.variables.mf_pipe[p,t],name='pressure loss return({0},{1})'.format(p,t))

        self.constraints.pressure_loss_supply = {}

        for t in time:
            for p in pipe_supply:
                for l in range(LL):
                    self.constraints.pressure_loss_supply[p, t, l] = m.addConstr(
                        2 * np.sqrt(self.data.pressure_loss_coeff[p]) * self.variables.mf_pipe[p, t],
                        gb.GRB.LESS_EQUAL,
                        np.sqrt(self.data.pr_supply_slice[l]) + (
                                    self.variables.pr_supply[self.data.pipe_supply_connexion[p][0], t] -
                                    self.variables.pr_supply[self.data.pipe_supply_connexion[p][1], t]) / np.sqrt(
                            self.data.pr_supply_slice[l]), name='pressure loss supply({0},{1},{2})'.format(p, t, l))

        self.constraints.pressure_loss_return = {}

        for t in time:
            for p in pipe_return:
                for l in range(LL):
                    self.constraints.pressure_loss_return[p, t, l] = m.addConstr(
                        2 * np.sqrt(self.data.pressure_loss_coeff[p]) * self.variables.mf_pipe[p, t],
                        gb.GRB.LESS_EQUAL,
                        np.sqrt(self.data.pr_return_slice[l]) + (
                                    self.variables.pr_return[self.data.pipe_return_connexion[p][0], t] -
                                    self.variables.pr_return[self.data.pipe_return_connexion[p][1], t]) / np.sqrt(
                            self.data.pr_return_slice[l]), name='pressure loss supply({0},{1},{2})'.format(p, t, l))

                    # wind realization

        self.constraints.wind_scenario = {}

        for t in time:
            for g in wind:
                self.constraints.wind_scenario[g, t] = m.addConstr(
                    self.variables.P[g, t],
                    gb.GRB.LESS_EQUAL,
                    self.data.wind_scenario[g, t] * self.data.elec_maxprod[g],
                    name='wind scenario({0},{1})'.format(n, t))

        # elec balance

        self.constraints.elec_balance = {}

        for t in time:
            for n in node:
                self.constraints.elec_balance[n, t] = m.addConstr(
                    gb.quicksum(self.variables.storage_plus[g, t] - self.variables.storage_moins[g, t] for g in
                                self.data.elec_storage_node[n]) + gb.quicksum(
                        self.variables.P[g, t] for g in self.data.elec_station_node[n]) + gb.quicksum(
                        self.variables.flow_line[l, t] for l in self.data.line_end[n]) - gb.quicksum(
                        self.variables.flow_line[l, t] for l in self.data.line_start[n]) - self.data.elec_load[n, t],
                    gb.GRB.EQUAL,
                    0, name='elec balance({0},{1})'.format(n, t))

                # elec transmission constraints

        self.constraints.angle_ref = {}

        for t in time:
            self.constraints.angle_ref[t] = m.addConstr(
                self.variables.node_angle[node[0], t],
                gb.GRB.EQUAL,
                0, name='angle ref({0})'.format(t))

        self.constraints.elec_flow = {}

        for t in time:
            for l in line:
                self.constraints.elec_flow[l, t] = m.addConstr(
                    self.variables.flow_line[l, t],
                    gb.GRB.EQUAL,
                    self.data.B[l] * (self.variables.node_angle[self.data.line_connexion[l][0], t] -
                                      self.variables.node_angle[self.data.line_connexion[l][1], t]),
                    name='elec flow({0},{1})'.format(l, t))

            #        self.constraints.time_delay_t00 = {}
        #
        #        for p in pipe_supply+pipe_return:
        #
        #                self.constraints.time_delay_t00[p] = m.addConstr(
        #                    self.variables.tau[p,'t00'],
        #                    gb.GRB.LESS_EQUAL,
        #                    0,name='time delay t00({0})'.format(p))

        self.constraints.balance_total = m.addConstr(
            gb.quicksum(self.variables.Q[h, t] for t in time for h in heat_station + heat_storage),
            gb.GRB.GREATER_EQUAL,
            gb.quicksum(self.data.heat_load[h, t] for t in time for h in node), name='balance overday')

    # %% SOLVE


dispatch_MILP = integrated_dispatch_MILP()
dispatch_MILP.model.params.OutputFlag = 1


#tstart = datetime.now()
#tstart = tt.time()
dispatch_MILP.optimize()
#tend = tt.time()
#tend = datetime.now()
#exectime=tend-tstart
syst_cost_MILP = dispatch_MILP.model.ObjVal

pipe = ['P{0}'.format(p + 1) for p in range(P)]

pipe_connexion = {'P1': ['N1', 'N2'], 'P2': ['N2', 'N3']}
pipe_start = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
pipe_end = {'N1': [], 'N2': [], 'N3': [], 'N4': [], 'N5': [], 'N6': []}
for n in node:
    for p in pipe:
        if pipe_connexion[p][0] == n:
            pipe_start[n].append(p)
        if pipe_connexion[p][1] == n:
            pipe_end[n].append(p)

# %% DH network parameters

pipe_maxflow = {p: 1000 for p in pipe}



class expando(object):
    '''


        A small class which can have attributes set
    '''
    pass


class integrated_dispatch_LP:
    def __init__(self):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self._load_data()
        self._build_model()

    def optimize(self):
        self.model.optimize()

    def computeIIS(self):
        self.model.computeIIS()

    def _load_data(self):

        # indexes
        self.data.time = time
        self.data.time_list = time_list
        self.data.node = node
        self.data.line = line
        self.data.pipe = pipe
        self.data.gen = gen
        self.data.heat_storage = heat_storage
        self.data.elec_storage = elec_storage
        self.data.heat_pump = heat_pump
        self.data.wind = wind
        self.data.heat_only = heat_only
        self.data.CHP_sorted = CHP_sorted
        self.data.CHP = CHP
        self.data.producers = producers
        self.data.heat_station = heat_station
        self.data.elec_station = elec_station
        self.data.heat_exchanger_station = heat_exchanger_station

        # producers sorted per node
        self.data.producers_node = producers_node
        self.data.heat_station_node = heat_station_node
        self.data.elec_station_node = elec_station_node
        self.data.CHP_node = CHP_node
        self.data.CHP_sorted_node = CHP_sorted_node
        self.data.heat_pump_node = heat_pump_node
        self.data.heat_only_node = heat_only_node
        self.data.heat_storage_node = heat_storage_node
        self.data.gen_node = gen_node
        self.data.wind_node = wind_node
        self.data.elec_storage_node = elec_storage_node
        self.data.heat_exchanger_station_node = heat_exchanger_station_node

        # connexions between nodes
        self.data.pipe_connexion = pipe_connexion
        self.data.pipe_start = pipe_start
        self.data.pipe_end = pipe_end
        self.data.line_connexion = line_connexion
        self.data.line_start = line_start
        self.data.line_end = line_end

        # LOADS
        self.data.heat_load = heat_load
        self.data.elec_load = elec_load

        # Heat station parameters
        self.data.CHP_maxprod = CHP_maxprod
        self.data.heat_maxprod = heat_maxprod
        self.data.rho_elec = rho_elec
        self.data.rho_heat = rho_heat
        self.data.r_min = r_min
        self.data.storage_rho_plus = storage_rho_plus
        self.data.storage_rho_moins = storage_rho_moins
        self.data.storage_maxcapacity = storage_maxcapacity
        self.data.storage_loss = storage_loss
        self.data.storage_init = storage_init
        self.data.COP = COP

        # DHN parameters
        self.data.pipe_maxflow = pipe_maxflow

        # Elec station parameters
        self.data.elec_maxprod = elec_maxprod
        self.data.wind_scenario = wind_scenario

        # Cost parameters
        self.data.alpha = alpha

        # elec transmission system
        self.data.B = B
        self.data.line_maxflow = line_maxflow

    def _build_model(self):

        self.model = gb.Model()
        self._build_variables()
        self._build_objective()
        self._build_constraints()

    def _build_variables(self):

        # indexes shortcuts
        time = self.data.time
        node = self.data.node
        line = self.data.line
        pipe = self.data.pipe
        gen = self.data.gen
        heat_storage = self.data.heat_storage
        elec_storage = self.data.elec_storage
        heat_pump = self.data.heat_pump
        wind = self.data.wind
        heat_only = self.data.heat_only
        CHP_sorted = self.data.CHP_sorted
        CHP = self.data.CHP
        producers = self.data.producers
        heat_station = self.data.heat_station
        elec_station = self.data.elec_station
        heat_exchanger_station=self.data.heat_exchanger_station
        m = self.model

        # heat market optimization variables

        self.variables.Q = {}  # heat production from CHPs and HO units (first satge)
        for t in time:
            for h in heat_station:
                self.variables.Q[h, t] = m.addVar(lb=0, ub=self.data.heat_maxprod[h], name='Q({0},{1})'.format(h, t))

        self.variables.flow_pipe = {}
        for t in time:
            for p in pipe:
                self.variables.flow_pipe[p, t] = m.addVar(lb=0, ub=self.data.pipe_maxflow[p],
                                                          name='flow pipe({0},{1})'.format(p, t))

        self.variables.storage_plus = {}  # heat storage: heat discharged (first stage)
        for t in time:
            for h in heat_storage + elec_storage:
                self.variables.storage_plus[h, t] = m.addVar(lb=0, ub=self.data.storage_maxprod[h],
                                                             name='storage plus({0},{1})'.format(h, t))

        self.variables.storage_moins = {}  # heat storage: heat charged (first stage)
        for t in time:
            for h in heat_storage + elec_storage:
                self.variables.storage_moins[h, t] = m.addVar(lb=0, ub=self.data.storage_maxprod[h],
                                                              name='storage moins({0},{1})'.format(h, t))

        self.variables.storage_energy = {}  # heat stored in heat storage h at end of time period t
        for t in time:
            for h in heat_storage + elec_storage:
                self.variables.storage_energy[h, t] = m.addVar(lb=0, ub=self.data.storage_maxcapacity[h],
                                                               name='storage energy({0},{1})'.format(h, t))

                # electricity market optimization variables : primal variables

        self.variables.P = {}  # electricity production from electricity generators, CHPs and wind producers
        for t in time:
            for g in CHP + gen + wind:
                self.variables.P[g, t] = m.addVar(lb=0, ub=self.data.elec_maxprod[g],
                                                  name='P({0},{1})'.format(g, t))  # dispatch of electricity generators

            for g in heat_pump:
                self.variables.P[g, t] = m.addVar(lb=-gb.GRB.INFINITY, ub=0,
                                                  name='P({0},{1})'.format(g, t))  # elec consumption of HP

        # electricity transmission system variables

        self.variables.node_angle = {}
        for t in time:
            for n in node:
                self.variables.node_angle[n, t] = m.addVar(lb=-gb.GRB.INFINITY, name='node angle({0},{1})'.format(n, t))

        self.variables.flow_line = {}
        for t in time:
            for l in line:
                self.variables.flow_line[l, t] = m.addVar(lb=-self.data.line_maxflow[l], ub=self.data.line_maxflow[l],
                                                          name='flow line({0},{1})'.format(l, t))

        m.update()

    def _build_objective(self):  # building the objective function for the heat maret clearing

        # indexes shortcuts
        time = self.data.time
        node = self.data.node
        line = self.data.line
        pipe = self.data.pipe
        gen = self.data.gen
        heat_storage = self.data.heat_storage
        elec_storage = self.data.elec_storage
        heat_pump = self.data.heat_pump
        wind = self.data.wind
        heat_only = self.data.heat_only
        CHP_sorted = self.data.CHP_sorted
        CHP = self.data.CHP
        producers = self.data.producers
        heat_station = self.data.heat_station
        elec_station = self.data.elec_station
        heat_exchanger_station=self.data.heat_exchanger_station
        m = self.model

        m.setObjective(
            gb.quicksum(self.data.alpha[g, t] * self.variables.P[g, t] for t in time for g in gen + wind) + gb.quicksum(
                self.data.alpha[g, t] * self.variables.Q[g, t] for t in time for g in heat_only) + gb.quicksum(
                self.data.alpha[g, t] * (
                            self.data.rho_elec[g] * self.variables.P[g, t] + self.data.rho_heat[g] * self.variables.Q[
                        g, t]) for t in time for g in CHP),
            gb.GRB.MINIMIZE)

    def _build_constraints(self):

        # indexes shortcuts
        time = self.data.time
        node = self.data.node
        line = self.data.line
        pipe = self.data.pipe
        gen = self.data.gen
        heat_storage = self.data.heat_storage
        elec_storage = self.data.elec_storage
        heat_pump = self.data.heat_pump
        wind = self.data.wind
        heat_only = self.data.heat_only
        CHP_sorted = self.data.CHP_sorted
        CHP = self.data.CHP
        producers = self.data.producers
        heat_station = self.data.heat_station
        elec_station = self.data.elec_station
        heat_exchanger_station=self.data.heat_exchanger_station
        m = self.model

        # heat balance!!

        self.constraints.heat_balance = {}

        for t in time:
            for n in node:
                self.constraints.heat_balance[n, t] = m.addConstr(
                    gb.quicksum(self.variables.storage_plus[i, t] - self.variables.storage_moins[i, t] for i in
                                self.data.heat_storage_node[n]) + gb.quicksum(
                        self.variables.Q[i, t] for i in self.data.heat_station_node[n]) + gb.quicksum(
                        self.variables.flow_pipe[p, t] for p in self.data.pipe_end[n]) - gb.quicksum(
                        self.variables.flow_pipe[p, t] for p in self.data.pipe_start[n]),
                    gb.GRB.EQUAL,
                    self.data.heat_load[n, t], name='heat balance({0},{1})'.format(n, t))

                ##1) CHPs

        self.constraints.CHP_maxprod = {}
        self.constraints.CHP_ratio = {}

        for t in time:
            for h in self.data.CHP_sorted['ex']:
                self.constraints.CHP_maxprod[h, t] = m.addConstr(
                    self.data.rho_heat[h] * self.variables.Q[h, t] + self.data.rho_elec[h] * self.variables.P[h, t],
                    gb.GRB.LESS_EQUAL,
                    self.data.CHP_maxprod[h], name='CHP maxprod({0},{1})'.format(h, t))

                self.constraints.CHP_ratio[h, t] = m.addConstr(
                    self.variables.P[h, t],
                    gb.GRB.GREATER_EQUAL,
                    self.data.r_min[h] * self.variables.Q[h, t], name='CHP ratio({0},{1})'.format(h, t))

            for h in self.data.CHP_sorted['bp']:
                self.constraints.CHP_ratio[h, t] = m.addConstr(
                    self.variables.P[h, t],
                    gb.GRB.EQUAL,
                    self.data.r_min[h] * self.variables.Q[h, t], name='CHP ratio({0},{1})'.format(h, t))

                ##2) storage

        self.constraints.storage_update = {}

        for (t1, t2) in zip(time[:-1], time[1:]):
            for h in heat_storage + elec_storage:
                self.constraints.storage_update[h, t2] = m.addConstr(
                    self.variables.storage_energy[h, t2],
                    gb.GRB.EQUAL,
                    self.variables.storage_energy[h, t1] - self.data.storage_rho_plus[h] * self.variables.storage_plus[
                        h, t2] + self.data.storage_rho_moins[h] * self.variables.storage_moins[h, t2] -
                    self.data.storage_loss[h], name='storage update({0},{1})'.format(h, t2))

        self.constraints.storage_init = {}

        for h in heat_storage + elec_storage:
            self.constraints.storage_init[h] = m.addConstr(
                self.variables.storage_energy[h, time[0]],
                gb.GRB.EQUAL,
                self.data.storage_init[h] - self.data.storage_rho_plus[h] * self.variables.storage_plus[h, time[0]] +
                self.data.storage_rho_moins[h] * self.variables.storage_moins[h, time[0]] - self.data.storage_loss[h],
                name='storage init({0})'.format(h))

        self.constraints.storage_final = {}

        for h in heat_storage + elec_storage:
            self.constraints.storage_final[h] = m.addConstr(
                self.variables.storage_energy[h, time[-1]],
                gb.GRB.GREATER_EQUAL,
                self.data.storage_init[h], name='storage final({0})'.format(h))

        # 3) heat pumps

        self.constraints.heat_pump_COP = {}

        for t in time:
            for h in heat_pump:
                self.constraints.heat_pump_COP[h, t] = m.addConstr(
                    -self.data.COP[h] * self.variables.P[h, t],
                    gb.GRB.EQUAL,
                    self.variables.Q[h, t], name='heat pump COP({0},{1})'.format(h, t))

                # wind realization

        self.constraints.wind_scenario = {}

        for t in time:
            for g in wind:
                self.constraints.wind_scenario[g, t] = m.addConstr(
                    self.variables.P[g, t],
                    gb.GRB.LESS_EQUAL,
                    self.data.wind_scenario[g, t] * self.data.elec_maxprod[g],
                    name='wind scenario({0},{1})'.format(n, t))

        # elec balance

        self.constraints.elec_balance = {}

        for t in time:
            for n in node:
                self.constraints.elec_balance[n, t] = m.addConstr(
                    gb.quicksum(self.variables.storage_plus[g, t] - self.variables.storage_moins[g, t] for g in
                                self.data.elec_storage_node[n]) + gb.quicksum(
                        self.variables.P[g, t] for g in self.data.elec_station_node[n]) + gb.quicksum(
                        self.variables.flow_line[l, t] for l in self.data.line_end[n]) - gb.quicksum(
                        self.variables.flow_line[l, t] for l in self.data.line_start[n]) - self.data.elec_load[n, t],
                    gb.GRB.EQUAL,
                    0, name='elec balance({0},{1})'.format(n, t))

                # elec transmission constraints

        self.constraints.angle_ref = {}

        for t in time:
            self.constraints.angle_ref[t] = m.addConstr(
                self.variables.node_angle[node[0], t],
                gb.GRB.EQUAL,
                0, name='angle ref({0})'.format(t))

        self.constraints.elec_flow = {}

        for t in time:
            for l in line:
                self.constraints.elec_flow[l, t] = m.addConstr(
                    self.variables.flow_line[l, t],
                    gb.GRB.EQUAL,
                    self.data.B[l] * (self.variables.node_angle[self.data.line_connexion[l][0], t] -
                                      self.variables.node_angle[self.data.line_connexion[l][1], t]),
                    name='elec flow({0},{1})'.format(l, t))

            # %% SOLVE


dispatch_LP = integrated_dispatch_LP()
dispatch_LP.model.params.OutputFlag = 1
# tstart = datetime.now()
# tstart = tt.time()
dispatch_LP.optimize()
# tend = tt.time()
# tend = datetime.now()
# exectime=tend-tstart
syst_cost_LP = dispatch_LP.model.objval

#-------------- РЕЗУЛЬТАТЫ ЧАСТИЧНО-ЦЕЛОЧИСЛЕННОГО ЛИНЕЙНОГОГО ПРОГРАММИРОВАНИЯ  ---------------------



P_LP = {}
for g in elec_station:
    for t in time:
        P_LP[g, t] = dispatch_LP.variables.P[g, t].x
for g in elec_storage:
    for t in time:
        P_LP[g, t] = dispatch_LP.variables.storage_plus[g, t].x - dispatch_LP.variables.storage_moins[g, t].x

Q_LP = {}
for g in heat_station:
    for t in time:
        Q_LP[g, t] = dispatch_LP.variables.Q[g, t].x
for g in heat_storage:
    for t in time:
        Q_LP[g, t] = dispatch_LP.variables.storage_plus[g, t].x - dispatch_LP.variables.storage_moins[g, t].x

    # %%


# %%

print('Затраты для полноценной интегрированной ЭЭС',syst_cost_MILP)
print('Затраты для ЭЭС с упрощённой тепловой частью', syst_cost_LP)
print('Затраты разница в % ',(syst_cost_LP-syst_cost_MILP)/syst_cost_LP*100)
# print('Cost % difference feasible',(syst_cost_LP-syst_cost_MILP_feasible)/syst_cost_LP*100)


P_LP = {}
for g in elec_station:
    for t in time:
        P_LP[g,t]=dispatch_LP.variables.P[g,t].x
for g in elec_storage:
    for t in time:
        P_LP[g,t]=dispatch_LP.variables.storage_plus[g,t].x - dispatch_LP.variables.storage_moins[g,t].x

Q_LP = {}
for g in heat_station:
    for t in time:
        Q_LP[g,t]=dispatch_LP.variables.Q[g,t].x
for g in heat_storage:
    for t in time:
        Q_LP[g,t]=dispatch_LP.variables.storage_plus[g,t].x - dispatch_LP.variables.storage_moins[g,t].x

wind_utilization_LP = sum(P_LP[w, t] for w in wind for t in time) / (
            sum(wind_scenario[w, t] * elec_maxprod[w] for w in wind for t in time))
print('wind_utilization_LP', wind_utilization_LP)

P_MILP = {}
for g in elec_station:
    for t in time:
        P_MILP[g,t]=dispatch_MILP.variables.P[g,t].x
for g in elec_storage:
    for t in time:
        P_MILP[g,t]=dispatch_MILP.variables.storage_plus[g,t].x - dispatch_MILP.variables.storage_moins[g,t].x

Q_MILP = {}
for g in heat_station:
    for t in time:
        Q_MILP[g,t]=dispatch_MILP.variables.Q[g,t].x
for g in heat_storage:
    for t in time:
        Q_MILP[g,t]=dispatch_MILP.variables.storage_plus[g,t].x - dispatch_MILP.variables.storage_moins[g,t].x

wind_utilization_MILP = sum(P_MILP[w, t] for w in wind for t in time) / (
            sum(wind_scenario[w, t] * elec_maxprod[w] for w in wind for t in time))
print('wind_utilization_MILP', wind_utilization_MILP)

import seaborn as sns
# Make a fake dataset:
height = [syst_cost_MILP*80/1000, syst_cost_LP*80/1000]
bars = ('IHPD model', 'CED model')
y_pos = [0,1]

#y_pos = np.arange(len(bars))

# Create bars
plt.bar(y_pos, height,  color=(0.1, 0.1, 0.1, 0.1),  edgecolor='blue')
plt.xticks(y_pos, bars)
plt.ylabel('Production costs, RUB * 10e3 ')
