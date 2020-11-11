# -*- coding: utf-8 -*-
"""
"""
import pandas as pd
import pypower as pp
import pyomo.core as pyomo
import numpy as np

#from pyomo.core.base.set import FiniteSimpleRangeSet
import numpy.matlib
from pypower.api import case30,makeYbus
from numpy.linalg import inv
from pyomo.environ import *
import pyomo.environ
from pyomo.bilevel import *
import pyomo.environ as pyo
from pyomo.opt import SolverFactory


#%% ПАРАМЕТРЫ ПЕРЕДАЮЩЕЙ СЕТИ

ppc = case30() # схема сети из pypower

nt = 24 # период моделирования
gennum = ppc["gencost"].shape[0] # считаем количество генераторов передающей сети (для case30 их 6)
busnum = ppc["bus"].shape[0] # считаем количество узлов передающей сети (для case30 их 30)
linenum = ppc["branch"].shape[0] # считаем количество узлов передающей сети (для case30 их 41)
 
# определим затраты
cg1 = ppc["gencost"][:, [5]].T # % линейный  коэффициент затрат на генерацию
cg2 = ppc["gencost"][:, [4]].T # квадратичный коэффициент затрат на генерацию
crgs = np.array([[6, 6.75, 7, 5.25, 5, 5]]) # коэффициент цены резерва генератора
conoff = 2*np.ones((1,gennum)) # коэффициент цены запуска агрегата

sampleL = np.array([[120.6, 115.8, 114.8, 112.6, 114.0, 113.4, 
                  117.1, 126.3, 130.7, 132.5, 135.6, 134.8, 
                  136.5, 137.7, 137.1, 138.0, 136.3, 133.3, 
                  131.7, 129.3, 128.2, 127.4, 125.6, 124.2]]) # образец профиля суммарной нагрузки системы передачи
sf = sampleL[0,:nt]/np.mean(sampleL[0,:nt])  # моделируем типа весовых коэффициентов для каждого узла схема
sf = np.matlib.repmat(sf,busnum,1)
loads = np.matlib.repmat(ppc["bus"][:, [2]],1,nt)
loads = loads*sf # получаем матрицу нагрузок на интервал моделирования nt для всех нагрузочных узлов схемы

# сделаем Y/B узловую матрицу и матрицу A в Pline=A*Pinj
Ybus, Yf, Yt = makeYbus(100,ppc["bus"],ppc["branch"])
B = (1j*Ybus).real.todense()
NB = -B
Bred = B[0:busnum-1,0:busnum-1] # редуцированная матрица В
# извлекаем топологию сети
frm = ppc["branch"][:,0]
to = ppc["branch"][:,1]
M = np.zeros((linenum,busnum)) # от/к матрице
lineC = np.zeros((linenum,1))

for i in range(linenum):
    M[i,int(frm[i]-1)] = 1
    M[i,int(to[i]-1)]= -1
    lineC[i]=1*ppc["branch"][i,6]
    
m = M[:,0:busnum-1] # матрица отражающая топологию сети, какие линии соединены с какими узлами схемы
x = ppc["branch"][:,[3]] # реактивное сопротивление линий
r = ppc["branch"][:,[2]] # активное соротивление линий
b = x/(x*x+r*r) # проводимость линии
D = np.diag(b[:,0]) # диагональная матрица проводимости


#%% Параметры оптимизации передающей сети и её ограничения

# Генераторы
one = np.ones(nt)
# Задаём максимальную и минимальную мощность генераторов
Gmax = np.array([80*one,80*one,40*one,50*one,30*one,55*one]) # задаём верхние органичения по вырабатываемой мощности для генераторов
Gmin = np.array(np.zeros((6,nt))) # задаём нижние органичения по вырабатываемой мощности для генераторов (фактически это нули)

# Задаём модель передающей сети с ограниченными диапазонами линий, генераторов и пр.
model = ConcreteModel()
model.bn = RangeSet(0,busnum-1)
model.bn1 = RangeSet(0,busnum-2)
model.nt = RangeSet(0,nt-1)
model.gn = RangeSet(0,gennum-1)
model.ln = RangeSet(0,linenum-1)

# Задаём переменные с помощью Pyomo Var function, которые мы хотим оптимизировать на уровне передающей сети 

# 1. добавляем в модель ограничения на максимальную и минимальную мощность генераторов, а также их статус работы
def b1(model,i,j):
    return (Gmin[i,j],Gmax[i,j])
model.pg = Var(model.gn, model.nt, bounds = b1) # задаём переменную генераторов с их ограничениями
model.onoff = Var(model.gn, model.nt, within = Binary) # задаём переменную статуса этих генераторов (включен/выключен)

# 2. добавляем в модель максимальный и минимальный рампинг генераторов (т.е. верхние и нижние предельно допустимые величины резерва мощности генератора)
Rgsmin = -0.5*Gmax
Rgsmax = 0.5*Gmax
def b2(model,i,j):
    return (Rgsmin[i,j],Rgsmax[i,j])
model.rgs = Var(model.gn, model.nt, bounds = b2)  # задаём переменную резерва мощности генераторов rgs с ограничением рампинга 
model.z = Var(model.gn, model.nt) # переменная мощности генераторов без ограничения на резерв

# 3. добавляем в модель органичения по пропускной способности ЛЭП
def b3(model,i,j):
    return (-lineC[i,0], lineC[i,0])
model.pflow =  Var(model.ln, model.nt, bounds=b3) # задаём переменную перетока мощности по линиям с ограничениями

model.Pinj = Var(model.bn, model.nt) #%  задаём переменную инъекции узловой полезной мощности для генерации и потребления (с учётом прогноза ветровой генерации)


model.theta =  Var(model.bn1, model.nt) # задаём переменную фазы напряжения в узлах схемы
Bredinv = inv(Bred)  #  берём обратную матрицу проводимостей B, т.е. получаем матрицу сопротивлений Х  

prod = np.dot(D,m) #  перемножаем диаг. матрицу проводимостей D и матрицу отражающая топологию сети m 
mbus = 7 # узел подключения микросети к передающей сети

#%% ПАРАМЕТРЫ РАСПРЕДЕЛИТЕЛЬНОЙ МИКРОСЕТИ 

#


lmp = 10.3*np.ones(nt)
cg = np.array([[3.3, 3.3, 3.3]])  # коэффициенты затрат на использование генераторов распредельной микросети, [евро/МВт]
cd = 1.01*np.array([[3.3, 3.3, 3.3, 3.3, 3.3, 3.3]]) # коэффициенты выгода микросети для потребления гибкой нагрузки, [евро/МВт] (немного больше чем cg, чтобы поощрить затраты) 
cb = 0.1 # коэффициент затрат на поддержание накопленной энергии, [евро/МВт] (должен быть небольшим, чтобы микросеть запасала энергию, а не экспортировала энергию)  
ci = lmp # цена импортируемой микросетью электроэнергии, [евро/МВт]
ce = 0.8*lmp #цена экспортируемой микросетью электроэнергии, [евро/МВт] (цена экспорта понятно меньше импорта)

ndl = 6 # количество гибких (активных) нагрузок
ng = 3 # количество генераторов в микросети
nb = 10 # количество накопителей электроэнергии

# Задаём модель распределительной микросети с ограниченными диапазонами линий, генераторов и пр.
model.sub = SubModel()
model.sub.nb = RangeSet(0,nb-1)
model.sub.ndl = RangeSet(0,ndl-1)
model.sub.ntt = RangeSet(0,nt-1)
model.sub.ng = RangeSet(0,ng-1)

# Задаём переменные с помощью Pyomo Var function, которые мы хотим оптимизировать на уровне распределитеной микросетии 
# 1. Активные нагрузки с их ограничениями 
one = np.ones(nt)
Pdmin = 0.1*np.array([0.5*one,4*one,2*one,5.5*one,1*one,7*one])
Pdmax = 0.03*np.array([10*one,16*one,15*one,20*one,27*one,32*one])

# добавляем в модель ограничения на регулирования активных нагрузок (ограничения управлением спросом)  
def b1(model,i,j):
    return (Pdmin[i,j],Pdmax[i,j])
model.sub.pd=Var(model.sub.ndl, model.sub.ntt, bounds = b1) # задаём переменную активных нагрузок с заданными ограничениями (ограничения управлением спросом)  

# 2. Пассивные нагрузки 
nload = 0.03 * np.array([[120.6, 115.8, 114.8, 112.6, 114.0, 113.4, 
                  117.1, 126.3, 130.7, 132.5, 135.6, 134.8, 
                  136.5, 137.7, 137.1, 138.0, 136.3, 133.3, 
                  131.7, 129.3, 128.2, 127.4, 125.6, 124.2]]) # профиль нагрузки для распределительной микросети 
nload = nload[:,0:nt]

# 3. Генераторы распредсети с ограничениями 
Gmaxs = 0.3 * np.array([5*one,4.5*one,7*one])
Gmins = np.array([1*one, 0.8*one, 1.5*one])
# добавляем в модель ограничения на генераторы распредсети   
def b2(model,i,j):
    return (Gmins[i,j],Gmaxs[i,j])
model.sub.pg = Var(model.sub.ng, model.sub.ntt, bounds=b2) # задаём переменную генераторов распредсети с ограничениями 
# добавляем в модель максимальный и минимальный рампинг генераторов 
Rgsmin = -0.3*Gmaxs
Rgsup = 0.3*Gmaxs
def b8(model,i,j):
    return (0,0) # однако пока ограничения нулевые 
model.sub.rgs = Var(model.sub.ng, model.sub.ntt, bounds=b8) # задаём переменную rgs генераторов распредсети с ограничением рампинга 
model.sub.z = Var(model.sub.ng, model.sub.ntt) # резервная переменная для rgs без ограничения рампинга 

# 4. Накопители
# задаём границы заряда/разряда накопителя 
Pbmin = -3*np.ones((nb,nt))
Pbmax = 3*np.ones((nb,nt))
# добавляем в модель ограничения на заряд/разряд накопителей (±3 МВт)
def b3(model,i,j):
    return (Pbmin[i,j],Pbmax[i,j])
model.sub.pb = Var(model.sub.nb, model.sub.ntt, bounds=b3) # задаём переменную для накопителя с ограничениями 

# задаём ограничения на состояния заряда накопителя
Bmin = np.zeros((nb,nt))
Bmax = 10*np.ones((nb,nt))
# добавляем в модель ограничения на состояния заряда накопителей (0 и 10 МВт)
def b4(model,i,j):
    return (Bmin[i,j],Bmax[i,j])
model.sub.b = Var(model.sub.nb, model.sub.ntt, bounds=b4)  # задаём переменную состояния заряда накопителя с ограничениями 

# 4. Импорт и экспорт мощности с ограничениям 
# добавляем в модель ограничения на импорт и экспорт мощности для распредсети 
def b5(model,i):
    return (0, 1000)
model.sub.ex = Var(model.sub.ntt, bounds=b5)   # задаём переменную импорта с ограничениями (0 и 1000 МВт) 

def b6(model,i):
    return (0,1000)
model.sub.im = Var(model.sub.ntt, bounds=b6)  # задаём переменную экспорта с ограничениями (0 и 1000 МВт)

Nmin = -20*one
Nmax = 20*one
def b7(model,i):
    return (Nmin[i],Nmax[i])
model.sub.net = Var(model.sub.ntt, bounds=b7)

# ОГРАНИЧЕНИЯ ОПТИМИЗАЦИОННОЙ МОДЕЛИ 

# Задаём список ограничений передающей сети 
model.cons = ConstraintList()

# Ограничения перетока мощности для генерирующих узлов Pinj = -L + Pg 
for i in range(busnum):
    for j in range(nt):
        if i==0:
            model.cons.add(model.Pinj[i,j]==-loads[i,j]+model.pg[0,j])
        elif i==1:
            model.cons.add(model.Pinj[i,j]==-loads[i,j]+model.pg[1,j])
        elif i==12:
            model.cons.add(model.Pinj[i,j]==-loads[i,j]+model.pg[2,j])
        elif i==21:
            model.cons.add(model.Pinj[i,j]==-loads[i,j]+model.pg[3,j])
        elif i==22:
            model.cons.add(model.Pinj[i,j]==-loads[i,j]+model.pg[4,j])
        elif i==26:
            model.cons.add(model.Pinj[i,j]==-loads[i,j]+model.pg[5,j])
        elif i==mbus:
            model.cons.add(model.Pinj[i,j]==-loads[i,j]+model.sub.net[j])
      #  elif i==mbus
       #   constraints=[constraints,Pinj(i,:)==-loads(i,:)+mg];
        else:
            model.cons.add(model.Pinj[i,j]==-loads[i,j]) # остальные чисто нагрузочные узлы с балансом Pinj = -L

 # задаём ограничения на расчёт потокраспределения передающей сети 
for r in range(busnum-1):
    for c in range(nt):
        model.cons.add(model.theta[r,c] == sum(Bredinv[r,i]*model.Pinj[i,c] for i in range(busnum-1)))  # задаём уравнение theta = Х*Pinj

for r in range(linenum):
    for c in range(nt):
        model.cons.add(model.pflow[r,c] == sum(prod[r,i]*model.theta[i,c] for i in range(busnum-1))) # задаём уравнение на переток мощности 
  
for h in range(nt):
    model.cons.add(sum(model.rgs[r,h] for r in range(gennum)) == sum(0.03*loads[r,h] for r in range(busnum))) 
    model.cons.add(sum(model.pg[r,h] for r in range(gennum)) == sum(loads[r,h] for r in range(busnum))) #задаём уравнение баланса мощности pg = l

# задаём ограничения на работу генераторов передающей сети с учётом их резерва
Rdn = 0.3*Gmax;
Rup = 0.3*Gmax;
for g in range(gennum):
    for h in range(nt-1):
        model.cons.add(-Rdn[g,h+1] <= model.pg[g,h+1]-model.pg[g,h] <= Rup[g,h+1]) # ограничения рампинга
        model.cons.add((model.pg[g,h+1]+model.rgs[g,h+1])-(model.pg[g,h]-model.rgs[g,h]) <= float(Rup[g,h+1])) # Ограничения по максимальному уровню резерва мощности 
        model.cons.add(float(-Rdn[g,h+1]) <= (model.pg[g,h+1]-model.rgs[g,h+1])-(model.pg[g,h]+model.rgs[g,h])) # Ограничения по минимальному уровню резерва мощности 

for g in range(gennum):
    for h in range(nt):
        model.cons.add(Gmin[g,h] * model.onoff[g,h] <= model.pg[g,h] + model.rgs[g,h])# ограничения для генератора с учётом резерва мощности (нижняя граница) 
        model.cons.add(model.pg[g,h] + model.rgs[g,h] <= Gmax[g,h] * model.onoff[g,h])# ограничения для генератора с учётом резерва мощности (верхняя граница)
        model.cons.add(model.rgs[g, h] <= model.z[g, h])# ограничения для текущего резерва мощности rgs<=z (т.е. резерв должен быть меньше или равен мощности генератора без ограничений)
        model.cons.add(-model.rgs[g,h] <= model.z[g,h])# ограничения для текущего резерва мощности -rgs<=z
        
#%% Целевая функция для передающей сети F
#def exps(model):
a = sum(crgs[0,i] * model.z[i,j] for i in range(gennum) for j in range(nt)) # первый член - это затраты на резерв генератора 
b = sum(model.onoff[i,j]*conoff[0,i] for i in range(gennum) for j in range(nt))# второй член - это затраты на запуск генератора 
c = sum(model.pg[i,j]*cg1[0,i] for i in range(gennum) for j in range(nt)) # третий член - это затраты на генерацию мощности (линейные) 
d = sum(model.pg[i,j]*model.pg[i,j]*cg2[0,i] for i in range(gennum) for j in range(nt)) # четвёртный член - это затраты на генерацию мощности (квадратичные) 
model.obj = Objective(expr = a+b+c+d, sense = minimize) # целевая функция F

# Задаём список ограничений микросети 
model.sub.cons = ConstraintList()

# задаём ограничения на работу генераторов микросети с учётом их резерва
Rdn = 0.3*Gmaxs
Rup = 0.3*Gmaxs
model.sub.cons.add(model.sub.b[0,0] <= 5) # задаём нижние ограничения состояния накопителя
model.sub.cons.add(model.sub.b[0,0] >= 5) # задаём верхние ограничения состояния накопителя

for g in range(ng):
    for h in range(nt-1):
        model.sub.cons.add(-Rdn[g,h+1] <= model.sub.pg[g,h+1]-model.sub.pg[g,h] <= Rup[g,h+1]) # ограничения рампинга
        model.sub.cons.add((model.sub.pg[g,h+1]+model.sub.rgs[g,h+1])-(model.sub.pg[g,h]-model.sub.rgs[g,h]) <= float(Rup[g,h+1])) # Ограничения по максимальному уровню резерва мощности 
        model.sub.cons.add(float(-Rdn[g,h+1]) <= (model.sub.pg[g,h+1]-model.sub.rgs[g,h+1])-(model.sub.pg[g,h]+model.sub.rgs[g,h])) # Ограничения по минимальному уровню резерва мощности 

for g in range(ng):
    for h in range(nt):
#        model.sub.cons.add(Gmin[g,h] <= model.sub.pg[g,h] + model.sub.rgs[g,h])# ограничения для генератора с учётом резерва мощности (нижняя граница)
        model.sub.cons.add(model.sub.pg[g,h] + model.sub.rgs[g,h] <= Gmaxs[g,h]) #ограничения для генератора с учётом резерва мощности (верхняя граница )
        model.sub.cons.add(model.sub.rgs[g,h] <= model.sub.z[g,h]) # ограничения для текущего резерва мощности rgs<=z (т.е. резерв должен быть меньше или равен мощности генератора без ограничений)
        model.sub.cons.add(-model.sub.rgs[g,h] <= model.sub.z[g,h]) # ограничения для текущего резерва мощности -rgs<=z

# ограничения на баланс мощности микросети (генерация - накопители - пассивная нагрузка - активная нагрузка = экспорт - импорт)
# также задаём ограничение что (экспорт - импорт = net)
for h in range(nt):
    model.sub.cons.add(sum(model.sub.pg[r,h] for r in range(ng)) >= sum(nload[r1,h] + model.sub.pb[r1,h] for r1 in range(1))+ model.sub.net[h] + sum(model.sub.pd[r2,h] for r2 in range(ndl)))
    model.sub.cons.add(sum(model.sub.pg[r,h] for r in range(ng)) <= sum(nload[r1,h] + model.sub.pb[r1,h] for r1 in range(1))+ model.sub.net[h] + sum(model.sub.pd[r2,h] for r2 in range(ndl)))
    model.sub.cons.add(model.sub.ex[h] - model.sub.im[h] <= model.sub.net[h])
    model.sub.cons.add(model.sub.ex[h] - model.sub.im[h] >= model.sub.net[h])

# ограничения на минимальный и максимальный уровень заряда накопителя (0 и 10 МВт)    
for i in range(nb):
    for j in range(nt-1):
        model.sub.cons.add(model.sub.b[i,j+1] >= model.sub.b[i,j]+model.sub.pb[i,j])
        model.sub.cons.add(model.sub.b[i,j+1] <= model.sub.b[i,j]+model.sub.pb[i,j])

# ограничения на текущую мощность заряда накопителя  
for i in range(nb):
    for j in range(nt):
        model.sub.cons.add(-model.sub.pb[i,j] <= model.sub.b[i,j])
        
#%% Целевая функция для микросети f
aa = sum(cb * model.sub.b[i,j] for i in range(nb) for j in range(nt))  # первый член - это затраты на поддержание накопленной энергии
bb = -sum(model.sub.pd[i,j]*cd[0,i] for i in range(ndl) for j in range(nt))  # второй член - это выгода от управления спросом в микросети (получается мы её максимизируем, там стоит минус) 
cc = sum(model.sub.pg[i,j]*cg[0,i] for i in range(ng) for j in range(nt)) # третий член - это затраты на генерацию микросети (от генераторов) 
dd = -sum(model.sub.ex[i]*ce[i] for i in range(nt)) # четвёртый член - это затраты на экспорт мощности в передающую сеть
ee = sum(model.sub.im[i]*ci[i] for i in range(nt)) # пятый член - это затраты на импорт мощности из передающей сети
model.sub.obj = Objective(expr = aa+bb+cc+dd+ee, sense = minimize) # Целевая функция f


#%% РЕШЕНИЕ НА БАЗЕ ДВУХУРОВНЕВОГО ПРОГРАММИРОВАНИЯ 
opt = SolverFactory("bilevel_blp_global")
opt.options["solver"] = 'gurobi'
results = opt.solve(model, tee=True)

#%% РЕЗУЛЬТАТЫ ОПТИМИЗАЦИИ 
print('--------------------------------------')
print('Оптимальная генерация передающей сети')
model.pg.display()
print('--------------------------------------')
print('Оптимальный  резерв мощности передающей сети')
model.rgs.display()

print('--------------------------------------')
print('Оптимальный  статус загрузки передающей сети')
model.onoff.display()

print('--------------------------------------')
print('Оптимальная  генерация микрогрида')
model.sub.pg.display()
print('--------------------------------------')
print('Оптимальный  резерв генерации микрогрида')
model.sub.rgs.display() #
print('--------------------------------------')
print('Оптимальная  выходная мощность накопителя в микрогриде')
model.sub.pb.display()
print('--------------------------------------')
print('Оптимальная состояния заряда/разряда накопителя в микрогриде')
model.sub.b.display()
print('--------------------------------------')
print('Оптимальная  мощность управления спросом')
model.sub.pd.display()
print('--------------------------------------')
print('Оптимальный экспорт мощности')
model.sub.ex.display()
print('--------------------------------------')
print('Оптимальный  импорт мощности')
model.sub.im.display() #
print('--------------------------------------')
print('Оптимальный net = экспорт - импорт')
model.sub.net.display()


print('--------------------------------------')
print('ОПТИМИМАЛЬНЫЕ ЗАТРАТЫ ПЕРЕДАЮЩЕЙ СЕТИ')
model.obj.display() #
print('--------------------------------------')
print('ОПТИМАЛЬНЫЕ ЗАТРАТЫ РАСПРЕДЕЛИТЕЛЬНОЙ МИКРОСЕТИ')

model.sub.obj.display()  #

#
#

#model.Pinj.display()
#model.pflow.display()


