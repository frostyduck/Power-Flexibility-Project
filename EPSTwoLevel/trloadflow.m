function [V,ang,IsConverged,bus_sol,line,BusStruct,Y] = trloadflow(bus,line,tol,iter_max,acc, display)
%В БУДУЩЕМ НЕОБХОДИМО УДАЛИТЬ СЛИШКОМ МАЛЕНЬКИЕ ВЕЛИЧИНЫ СОПРОТИВЛЕНИЙ, ЧТОБЫ
%ОБЕСПЕЧИТЬ СХОДИМОСТЬ!!! ПРИМЕР:
if strcmp(display,'on')
    disp('ПОМНИ ПРО РЕФЕРЕНСНЫЕ ЗНАЧЕНИЯ НАПРЯЖЕНИЯ В PU УЗЛАХ, ИХ МОЖНО И НУЖНО ВЫСТАВЛЯТЬ!!!');
    disp('ПОМНИ ПРО ВОЗМОЖНОСТЬ НАСТРОЙКИ МНОЖЕСТВЕННЫХ ПЕРЕКЛЮЧЕНИЙ PU->PQ->PU!!!');
    disp('ПОМНИ ПРО ВОЗМОЖНОСТЬ КОРРЕКТИРОВКИ КОЭФ.ТРАНСФ. ДЛЯ УЛУЧШЕНИЯ СХОДИМОСТИ!!!');
end
BusStruct=struct();%Структура, в которой описываются номера узлов, которые тносятся к PQ и PU, и тоько PQ узлов
BusStruct.CheckTapChangers=0;%ПЛОХО РАБОТАЕТ НЕ ИСПОЛЬЗОВАТЬ! Настраиваем коэффициенты трансформации для улучшения сходимости 1-да, 0-нет
BusStruct.MultiplePUPQ=1;%Множественные переключения PU->PQ->PU 1-да, 0-нет.
BusStruct.ConsiderUrefNewVersion=1;%1-используем в качестве референсных значений 16 столбец матрицы 0-используем в качестве референсных значений исходные знаения напряжений
if BusStruct.MultiplePUPQ%Множественные переключения PU->PQ->PU
    if BusStruct.ConsiderUrefNewVersion
        BusStruct.Uref0=bus(:,16);%используем в качестве референсных значений 16 столбец матрицы        
        BusStruct.NumOfPUPQSwitches=zeros(size(bus,1),1);%Число переключений PU->PQ, ЗАЩИТА ОТ ЗАЛИПАНИЯ ПЕРЕКЛЮЧЕНИЙ
    else        
        BusStruct.Uref0=bus(:,2);%Запоминаем исходное значение напряжений по всем узлам, для PU узлов это будут Uref
    end
end

[nbus] = size(bus,1);     % number of buses

CheckFlag=1;%Специальный флаг, который следит, чтобы использовались только решатели, которые проверяют проверку сходимости и с которых можно получить чувствительность решения.
IsConverged=0;%По умолчанию считаем, что УР разошелся!

ExtIter=0;%Outer iteration PU->PQ tripping
% Устаноква параметров расчета по умолчанию
if isempty(tol);tol = 1e-11;end
if isempty(iter_max);iter_max = 30;end

Y=y_sparseNew(bus,line);

bus_type = bus(:,10);
%Сначала задаем всем узлам напряжение номинала bus(:,2), затем, при
%необходимости (if BusStruct.ConsiderUrefNewVersion), выставляем для PU
%узлов требуемые значения
V = ones(nbus,1);%Сначала задаем всем узлам напряжение номинала bus(:,2);%НЕ удаляй здесь ничего!!!
%Функция коррекции КТ не работает, не использовать
if BusStruct.CheckTapChangers%bus_int нужна для корреткировки коэф. трансформации для улучшения сходимости
    busmax = max(bus(:,1));
    bus_int = zeros(busmax,1);
    ibus = (1:nbus)';
    bus_int(bus(:,1)) = ibus;    
    Ndx1=bus_int(line(:,1));
    Ndx2=bus_int(line(:,2));
end
%при необходимости (if BusStruct.ConsiderUrefNewVersion), выставляем для PU узлов требуемые значения
BusStruct.PU_no=find(bus_type==2);%Positions of PU nodes
if BusStruct.ConsiderUrefNewVersion%Новая версия с Uref для PU узлов в 16 столбце 
    V(BusStruct.PU_no) = bus(BusStruct.PU_no,16);%Выставляем в PU узлах Uref
end

ang = bus(:,3)*0;
%ang = bus(:,3);
Pg = bus(:,4);
Qg = bus(:,5);
Pl = bus(:,6);
Ql = bus(:,7);
qg_max = bus(:,11);
qg_min = bus(:,12);
sw_bno=ones(nbus,1);
g_bno=sw_bno;
% set up index for Jacobian calculation
%% form PQV_no and PQ_no
bus_zeros=zeros(nbus,1);
swing_index=find(bus_type==1);%swing bus index
sw_bno(swing_index)=bus_zeros(swing_index);%
BusStruct.PQV_no=find(bus_type >=2);%Positions of all nodes except swing bus %Узлы, которые соответствуют всем узлам кроме БУ
BusStruct.PQ_no=find(bus_type==3);%Positions of PQ nodes %Узлы, которые соответствуют PQ узлам
BusStruct.InitialPU=BusStruct.PU_no;%В исходном состоянии данные узлы были PU, только они могут переключаться PU->PQ
BusStruct.SW_no=swing_index;%Индекс БУ
BusStruct.SWPU_no=[swing_index; BusStruct.PU_no];%Индексы SW и PU узлов
g_bno(BusStruct.PU_no)=bus_zeros(BusStruct.PU_no);

st = clock;     % start the iteration time clock
if strcmp(display,'on')
    options = optimoptions('fsolve','Jacobian','on','TolFun',tol,'TolX',tol,'Display','iter-detailed','MaxIterations',iter_max);%Включаем вывод
else
    options = optimoptions('fsolve','Jacobian','on','TolFun',tol,'TolX',tol,'Display','off','MaxIterations',iter_max);%Выключаем вывод результатов
end
%% start iteration process for main Newton_Raphson solution
%Решение УР формируется в виде двойного цикла, на каждом цикле
%происходит решение уравнений без перефиксации PU->PQ узлов, после
%каждого решения проверяется необходимость перефиксации, в случае
%необходимости перефиксируются либо все разом, либо один узел с максимальным небалансом по Q (условие выбора внутри цикла)
while 1%Внешний цикл, проверка и перефиксакция PU->PQ
    if (ExtIter>=1)
                
        [lim_flag,Qg,bus_type,BusStruct] = chq_lim(qg_max,qg_min,bus_type,Qg,display,BusStruct,bus,V);
        
        if lim_flag %Производились какие-то модификации узлов в chq_lim
            BusStruct.PQV_no=find(bus_type >=2);
            BusStruct.PQ_no=find(bus_type==3);
            BusStruct.PU_no=find(bus_type==2);%Positions of PU nodes
            BusStruct.SWPU_no=[swing_index; BusStruct.PU_no];%Индексы SW и PU узлов
        else
            if strcmp(display,'on')
                disp(['Load flow converged after ',num2str(ExtIter),' external iterations']);
            end
            break%No PU->PQ switching, power flow convrged! Переключений узлов не было, режим сошелся            
        end
    end
    %В процессе расчета может произойти обратное переключение с PQ на PU,
    %поэтому необходимо вновь выставить заданное напряжение в PU узле
    if BusStruct.ConsiderUrefNewVersion%Новая версия с Uref для PU узлов в 16 столбце
        V(BusStruct.PU_no) = bus(BusStruct.PU_no,16);%Выставляем в PU узлах Uref
    end
    %Корректировка КТ не работает! Не используй!
    if BusStruct.CheckTapChangers%если нужна корреткировка коэф. трансформации Криво работает, не использую данную функцию, надо с ней разбираться
        V(Ndx1)=V(Ndx2).*line(:,6);
    end    
    V0=V(BusStruct.PQ_no);%Напряжения для PQ узлов берем из исходных, количенство PQ узлов может меняться вследствие перефикисрования PU-PQ-PU    
    ang0=ang(BusStruct.PQV_no);
    x0=[ang0;V0];
    
    f=@(x)myfunwithJac(x,V,ang,Y,Pg,Qg,Pl,Ql,BusStruct);
    [x,fval,~] = fsolve(f,x0,options);
    ang(BusStruct.PQV_no)=x(1:length(BusStruct.PQV_no));%Исправил!
    V(BusStruct.PQ_no)=x(length(BusStruct.PQV_no)+1:length(x));%Исправил!
    %Проверяем сходимость УР
    ndx=find(abs(fval)>options.TolFun,1);%If at least one element bigger than TolFun then loadflow diverged
    if ~isempty(ndx)%load flow diverged
        IsConverged=0;
        bus_sol=bus;
        line;
        return
    else
        IsConverged=1;
        PFTolFun=options.TolFun;%указываем чувстительность уравнений, необходимо для исключения ненужных величин по генерации в конце
    end    
    Q = calcQ(V,ang,Y);
    %Суммарная генерация в узле равна генерации во внешнюю сеть Q плюс
    %генерация, потребляемая нагрузкой Ql
    %Qg(BusStruct.PU_no) = Q(BusStruct.PU_no) + Ql(BusStruct.PU_no);
    Qg = Q + Ql;
    if strcmp(display,'on')
        disp(['Current external load flow iteration - ',num2str(ExtIter)]);
    end
    ExtIter=ExtIter+1;    
end
% voltage in rectangular coordinate
V_rect = V.*exp(1i*ang);
% bus current injection
cur_inj = Y*V_rect;
% power output based on voltages
S = V_rect.*conj(cur_inj);
P = real(S); Q = imag(S);
Pg=P+Pl;%Активная Генерация по всем узлам, включая балансирующий
notGenNdx=find(abs(Pg)<PFTolFun);%ищем все элементы, которые меньше чувствительности
Pg(notGenNdx)=0;
Qg=Q+Ql;%Реактивная Генерация по всем узлам, включая балансирующий
notGenNdx=find(abs(Qg)<PFTolFun);%ищем все элементы, которые меньше чувствительности
Qg(notGenNdx)=0;
bus_sol=bus;
bus_sol(:,2)=V;
bus_sol(:,3)=ang*(180/pi);
bus_sol(:,4)=Pg;%Тупо впечатываем генерацию по всем узлам
bus_sol(:,5)=Qg;%Тупо впечатываем генерацию по всем узлам
return
end
function F=myfun(x,V,ang,Y,Pg,Qg,Pl,Ql,BusStruct)%Finite Differences Jacobian Численный якобиан
ang(BusStruct.PQV_no)=x(1:length(BusStruct.PQV_no));%Исправил!
V(BusStruct.PQ_no)=x(length(BusStruct.PQV_no)+1:length(x));%Исправил!
V_rect = V.*exp(1i*ang);
cur_inj = Y*V_rect;
S = V_rect.*conj(cur_inj);
P = real(S); Q = imag(S);
delP = Pg - Pl - P;
delQ = Qg - Ql - Q;
F=[delP(BusStruct.PQV_no);delQ(BusStruct.PQ_no)];
end
function [F,Jac]=myfunwithJac(x,V,ang,Y,Pg,Qg,Pl,Ql,BusStruct)%Analytic Jacobian Аналитический Якобиан
ang(BusStruct.PQV_no)=x(1:length(BusStruct.PQV_no));%Исправил!
V(BusStruct.PQ_no)=x(length(BusStruct.PQV_no)+1:length(x));%Исправил!
Jac=-(form_jacMine(V,ang,Y,BusStruct));%Analytic Jacobian Аналитический Якобиан
V_rect = V.*exp(1i*ang);
cur_inj = Y*V_rect;
S = V_rect.*conj(cur_inj);
P = real(S); Q = imag(S);
delP = Pg - Pl - P;
delQ = Qg - Ql - Q;
F=[delP(BusStruct.PQV_no);delQ(BusStruct.PQ_no)];
end

