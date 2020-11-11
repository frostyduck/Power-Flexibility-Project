    clc;
    clear all;
    close all;
%% Параметры передающей сети
    nt=24; % временной горизонт моделирования 
    mpc = case30;  % загрузка файлы схемы передающей сети (ieee 30)
    [ng ~] = size(mpc.gen); % определяем количество генераторов
    [bus,~] = size(mpc.bus);% определяем количество узлов нагрузки
    cg1 = mpc.gencost(:,6); % линейный коэффциент цены (рубль/МВт)
    cg2 = mpc.gencost(:,5); % квадратичный коэффциент цены (рубль/МВт)
    crg=[6 6.75 7 5.25 5 5]'; % коэфициент цены резервов генератора(рубль/МВт)
    % Вектор нагрузки передающей сети для исследдуемого временного горизонта 
    sampleL=[120.6 115.8 114.8 112.6 114.0 113.4 117.1 126.3 130.7 132.5 135.6 134.8 136.5 137.7 137.1 138.0 136.3 133.3 131.7 129.3 128.2 127.4 125.6 124.2];
    sf=sampleL(:,1:nt)./mean(sampleL(:,1:nt));
    sf=repmat(sf,bus,1);
    loads=repmat(mpc.bus(:,3),1,nt); 
    loads=loads.*sf;
    loads(2,:) = 0;
    loads(8,:) = 0;
    wbus1 = 6; % определяем шину подключения ветряка 1
    wbus2 = 4; % определяем шину подключения ветряка 2
    dbus1 = 5; % определяем шину подключения распредсети 1 
    dbus2 = 15; % определяем шину подключения распредсети 2 
    dbus3 = 25; % определяем шину подключения распредсети 3 
    wc = 1000; % стоимость ветроэнергии
    lc = 1000; % стоимость отключения нагрузки 
    GSF = makePTDF(mpc.baseMVA, mpc.bus, mpc.branch, bus);% матрица коэффициентов сдвига генерации для системы передачи
    lineC = 2*repmat(mpc.branch(:,6),1,nt); % пропускная способность ЛЭП
    Conoff=2*ones(1,ng);    
 
%% Параметры распределительной микросети 1
   dbusA = 6;  % количество узлов распредсети
   totalA = [0, 4, 5, 3, 4, 2]; % нагрузка по узлам 
   totalA = repmat(totalA',1,nt);
   dssizeA = 1; % коэффицент масштабирования мощности распредсети (если 0 система отключена)
   cd2A = 0.02; % квадратичный коэффциент цены генерации микросети
   cd1A = 5;  % линейный коэффциент цены генерации микросети
   PA = sdpvar(dbusA,nt,'full'); % график диспетчеризации генерации микросети 
   pdA1 = sdpvar(dbusA,nt,'full');
   pdA1up = totalA*dssizeA;
   pdA1dn = 0.00001*totalA*dssizeA;
   pdA = 0*(totalA-pdA1up);
   CpdA1 = 3; % выгода микросети от испрользования активной нагрузки 
   PAup = repmat(0.8*[14.0000   8.2500    4.5000    6.5000    3.7500    0.0000]',1,24); %[14.0000   10.2500    5.5000    5.5000    1.7500    0.0000]
   PAdn = 0;
   dgA = sdpvar(1,nt,'full');
   dgAup = 0;
   dgAdn = 0;
   % Верхние и нижние границы управления спросом микросети 
   drupA = sdpvar(dbusA,nt,'full');
   drdnA = sdpvar(dbusA,nt,'full');
   drA1 = 1; % линейный коэффициент управления спросом
   drA2 = 1; % квадратичный коэффициент управления спросом
   drscale = 0.5; %коэффициент масштабирования активной нагрузки (сейчас это 50%)
 %% Параметры распределительной микросети 2
   dbusB = 6;
   totalB = [0, 4, 5, 3, 4, 2];
   totalB = repmat(totalB',1,nt);
   dssizeB = 1;
   cd2B = 0.02;
   cd1B = 5;  
   PB = sdpvar(dbusB,nt,'full');
   pdB1 = sdpvar(dbusB,nt,'full');
   pdB1up = totalB*dssizeB;
   pdB1dn = 0.00001*totalB*dssizeB;
   pdB = 0*(totalB-pdB1up);
   CpdB1 = 3;
   PBup = repmat(0.8*[14.0000   8.2500    4.5000    6.5000    3.7500    0.0000]',1,24); %[14.0000   10.2500    5.5000    5.5000    1.7500    0.0000]
   PBdn = 0;
   dgB = sdpvar(1,nt,'full');
   dgBup = 0;
   dgBdn = 0;
   drupB = sdpvar(dbusB,nt,'full');
   drdnB = sdpvar(dbusB,nt,'full');
   drB1 = 1;
   drB2 = 1;
   drscaleB = 0.5;
%% Параметры распределительной микросети 3
   dbusE = 6;
   totalE = [0, 4, 5, 3, 4, 2];
   totalE = repmat(totalE',1,nt);
   dssizeE = 1;
   cd2E = 0.02;
   cd1E = 5;  
   PE = sdpvar(dbusE,nt,'full');
   pdE1 = sdpvar(dbusE,nt,'full');
   pdE1up = totalE*dssizeE;
   pdE1dn = 0.00001*totalE*dssizeE;
   pdE = 0*(totalE-pdE1up);
   CpdE1 = 3;
   PEup = repmat(0.8*[14.0000   8.2500    4.5000    6.5000    3.7500    0.0000]',1,24); %[14.0000   10.2500    5.5000    5.5000    1.7500    0.0000]
   PEdn = 0;
   dgE = sdpvar(1,nt,'full');
   dgEup = 0;
   dgEdn = 0;
   drupE = sdpvar(dbusE,nt,'full');
   drdnE = sdpvar(dbusE,nt,'full');
   drE1 = 1;
   drE2 = 1;
   drscaleE = 0.5;
%% Прогноз ветромощности и моделирование сценариев
    %Прогноз ветромощности
    %Ветряк 1
    wl = 1; % коэффициент масштабирования проникновения ветровой энергии 
    load('ObsDays.mat'); %  загружаем исторический данные изменения ветромощности на основе NREL-Eastern Wind Integration Study
    k = 54; % набора из 54 траекторий, представляющих возможные реализации ветра
    wf1=wl*ObsDays(k,:,11)+2;
    corr1=corrcoef(ObsDays(:,:,11));
    std_dev = [0.12,0.15,0.18,0.5,0.6,0.67,0.72,0.76,0.79,0.82,0.83,0.8315,0.833,0.835,0.836,0.838,0.839,0.841,0.842,0.844,0.845,0.847,0.848,0.85];
    for j=1:nt
        for k=1:nt
        covarr1(j,k)=corr1(j,k)*std_dev(k)*std_dev(j);
        end
    end
    beta=0.01;
    eulernum=exp(1);
    epsilon=0.01;
    Ndelta=1*nt;%количество неопределённостей, одна ветростанция для каждого периода
    Nneed=ceil((1/epsilon)*(eulernum/(eulernum-1))*(log(1/beta)+4*(Ndelta+1)-1));% количество сценариев
    
    %создаём ошибки в наборе данных ветромощности 
    avg=zeros(nt,1)';
    for i=1:Nneed
        err1(i,:)=mvnrnd(avg,covarr1);
    end
    for j=1:Nneed
        wn1(j,:)=wf1+wf1.*err1(j,:);
        pos=wn1(j,:)>=0;
        wst(j,:)=wn1(j,:).*pos;
    end
    winderror=wst-repmat(wf1,Nneed,1); 
    
    %создаём ошибки в наборе данных ветромощности для ветряка 1
    pr1=100;%верхний процентиль (на всем протяжении)
    pr2=0;%нижний процентиль (на всем протяжении)
    windup1=prctile(winderror,pr1); % верхняя граница отклонений прогноза ветромощности 
    winddown1=prctile(winderror,pr2); % нижняя граница отклонений прогноза ветромощности 

    % ветряк 2
    wl=1; % коэффициент масштабирования проникновения ветровой энергии
    wf2=wl*ObsDays(54,:,11);
    for j=1:Nneed
        wn1(j,:)=wf2+wf2.*err1(j,:);
        pos=wn1(j,:)>=0;
        wst(j,:)=wn1(j,:).*pos;
    end
    winderror=wst-repmat(wf2,Nneed,1);   
    
    %создаём ошибки в наборе данных ветромощности для ветряка 2
    windup2=prctile(winderror,pr1);
    winddown2=prctile(winderror,pr2);
    
    % ветряк 1
    scale = 0.9;
    wf1 = scale*(wf1);
    windup1 = scale*(windup1);
    winddown1 = scale*(winddown1);
    
    windup1 = 1*[2.14850277468970,2.38176933526874,1.74677853896453,2.60511542849926,2.35355880939503,2.41625640922597,9.26333481098347,4.72603747192732,3.64443832323432,3.54370589728489,3.77247031861537,5.32507128646877,3.89024919339414,3.43366439452552,3.79975232769689,5.91542122062703,18.0071565379186,20.3629232486021,27.8562658759252,16.9300600310706,14.6166013530343,12.8505562939878,9.60886741793157,16.5816874966563];
    winddown1 = 1*[-2.29653290221721,-2.29980619610469,-1.86455639333191,-1.35000000000000,-1,-1,-3.40000000000000,-1.45000000000000,-1.05000000000000,-1,-1.10000000000000,-1.55000000000000,-1.15000000000000,-1,-1.10000000000000,-1.95000000000000,-5.70000000000000,-6.30000000000000,-9.10000000000000,-5.65000000000000,-4.55000000000000,-3.90000000000000,-3.15000000000000,-5.15000000000000];

    % ветряк 2 
    scale2 = 0.5;
    wf2 = scale2*(wf2);
    windup2 = scale2*(windup2);
    winddown2 = scale2*(winddown2);
    
      % Ветряк 1 + Ветряк 2
    wf = wf1+wf2;
    windup = windup1+windup2;
    winddown = winddown1+winddown2;
%% Переменные оптимизации передающей системы и их ограничения 
    % Генераторы 
    Gmax = [];
    for i = 1:ng
        Gmax = 2*[Gmax;mpc.gen(i,9)*ones(1,nt)];
    end
    Gmin=[zeros(ng,nt)];
    pg=sdpvar(ng,nt,'full');
    onoff=ones(ng,nt);
%     onoff=binvar(ng,nt,'full');


    Rgmax=0.1*Gmax;
    rgup=sdpvar(ng,nt,'full'); % верхние резервы мощности генерации 
    rgdn=sdpvar(ng,nt,'full'); % нижние резервы мощности генерации 
    Rdn=0.3*Gmax; % нижний рампинг
    Rup=0.3*Gmax; % верхний рампинг
    
    % Цены импорта и экспорта 
    Pim = sdpvar(1,nt,'full');
    Pimmax = 10*ones(1,nt);
    Pimmin = 0*ones(1,nt);
    
    % Цена управления спросом 
    Pdr=sdpvar(1,nt,'full');
    Pdrmax = 10;
    
    cimA = sdpvar(1,nt,'full'); % импортируемая мощность микросети 1
    cimB = sdpvar(1,nt,'full'); % импортируемая мощность микросети 2
    cimE = sdpvar(1,nt,'full'); % импортируемая мощность микросети 3
    drpA = sdpvar(1,nt,'full'); % цена управления спросом микросети 1
    drpB = sdpvar(1,nt,'full'); % цена управления спросом микросети 2
    drpE = sdpvar(1,nt,'full'); % цена управления спросом микросети 3
    % Нижние и верхние границы
    CO = [Gmin<=pg<=Gmax,0<=Pdr<=Pdrmax,0<=rgup<=Rgmax,0<=rgdn<=Rgmax,Pimmin<=Pim<=Pimmax,0<=cimA<=30,0<=drpA<=10,0<=cimB<=30,0<=drpB<=10,0<=cimE<=30,0<=drpE<=10];      
    %% Ограничения передающей системы 
    % Переток мощности с ошибками прогноза ветромощности 
    genbus = mpc.gen(:,1);
    genvec = zeros(1,bus);
    for i = 1:length(genbus)
        genvec(genbus(i))=1;
    end
    
     Pinj=sdpvar(bus,nt,'full'); % узловая матрица с ошибками прогноза ветромощности        
    count = 1;
    for i=1:bus
      if genvec(i) ~= 1
          if i == wbus1
              CO = [CO,Pinj(i,:)==-loads(i,:)+wf1];
          elseif i == wbus2
              CO = [CO,Pinj(i,:)==-loads(i,:)+wf2];
          elseif i == dbus1

          elseif i == dbus2

          elseif i == dbus3

          else
              CO=[CO,Pinj(i,:)==-loads(i,:)];
          end
      else
          CO=[CO,Pinj(i,:)==-loads(i,:)+pg(count,:)];
          count = count+1;
      end
    end

    CP = [-lineC<=GSF*Pinj<=lineC];
    
    Pinj1=sdpvar(bus,nt,'full'); % узловая матрица с ошибками прогноза ветромощности (ветряк 1)            
    count = 1;
    for i=1:bus
      if genvec(i) ~= 1
          if i == wbus1
              CO = [CO,Pinj1(i,:)==-loads(i,:)+wf1+winddown1];
          elseif i == wbus2
              CO = [CO,Pinj1(i,:)==-loads(i,:)+wf2+winddown2];
          elseif i == dbus1
              
          elseif i == dbus2

          elseif i == dbus3

          else
              CO=[CO,Pinj1(i,:)==-loads(i,:)];
          end
      else
          CO=[CO,Pinj1(i,:)==-loads(i,:)+pg(count,:)+rgup(count,:)];
          count = count+1;
      end
    end

    Pinj2=sdpvar(bus,nt,'full'); % узловая матрица с ошибками прогноза ветромощности (ветряк 2)        
    count = 1;
    for i=1:bus
      if genvec(i) ~= 1
          if i == wbus1
              CO = [CO,Pinj2(i,:)==-loads(i,:)+wf1+windup1];
          elseif i == wbus2
              CO = [CO,Pinj2(i,:)==-loads(i,:)+wf2+windup2];
          elseif i == dbus1
              
          elseif i == dbus2

          elseif i == dbus3

          else
              CO=[CO,Pinj2(i,:)==-loads(i,:)];
          end
      else
          CO=[CO,Pinj2(i,:)==-loads(i,:)+pg(count,:)-rgdn(count,:)];
          count = count+1;
      end
    end
    
    CO = [CO,-lineC<=GSF*Pinj2<=lineC];
    
   % Ограничения резервов генератора 
   CO = [CO,Gmin.*onoff<=pg<=Gmax.*onoff]; % Генератор граничит с коэффицентом рампинга Rgs
    CO = [CO,pg+rgup<=Gmax.*onoff,Gmin.*onoff<=pg-rgdn]; 
    CO = [CO,(pg(:,2:nt)+rgup(:,2:nt))-(pg(:,1:nt-1)-rgdn(:,1:nt-1))<=Rup(:,2:nt)]; % верхний предел коэффицента рампинга генератора 
    CO = [CO,-Rdn(:,2:nt)<=(pg(:,2:nt)-rgdn(:,2:nt))-(pg(:,1:nt-1)+rgup(:,1:nt-1))]; % нижний предел коэффицента рампинга генератора
    CO = [CO,windup-sum(rgdn)-sum(drupE)-sum(drupB)-sum(drupA)==0,-winddown-sum(rgup)-sum(drdnE)-sum(drdnB)-sum(drdnA)==0]; % note
%     CO = [CO,windup-sum(rgdn)-sum(drupA)-sum(drupB)-sum(drupE)==0,-winddown-sum(rgup)-sum(drdnA)-sum(drdnB)-sum(drdnE)==0]; % note
    % Ограничения генератора
    CO = [CO,-Rdn(:,2:nt)<=pg(:,2:nt)-pg(:,1:nt-1)<=Rup(:,2:nt)]; % рампинг добавляем в CO

%   EP=[sum(pg)-sum(loads)+wf-PA(1,:)-PB(1,:)-PE(1,:)>=0];   %note

    EP = [sum(pg)-sum(loads)+wf-PE(1,:)-PE(2,:)-PB(1,:)-PB(2,:)-PA(1,:)-PA(2,:)>=0];   
%% Ограничения микросети 1
   CDA = [PAdn<=PA<=PAup, pdA1dn<=pdA1<=pdA1up, 0<=drupA<=drscale*pdA1, 0<=drdnA<=drscale*pdA1];
   for i = 1:dbusA-1
       if i == 1
           CDA = [CDA, PA(4,:) == PA(1,:) - pdA(2,:) - pdA1(2,:)];
       elseif i == 2
           CDA = [CDA, PA(3,:) == PA(2,:) - pdA(5,:) - pdA1(5,:)];
       elseif i == 3
           CDA = [CDA, 0 == PA(3,:) - pdA(6,:) - pdA1(6,:)];
       elseif i == 4
           CDA = [CDA, PA(5,:) == PA(4,:) - pdA(3,:) - pdA1(3,:)];
       elseif i == 5
           CDA = [CDA, 0 == PA(5,:) - pdA(4,:) - pdA1(4,:)];
       end
   end 
   CDA = [CDA,PA(1,:)+PA(2,:) >= sum(pdA + pdA1)]; 
 %% Ограничения микросети 2
   CDB = [PBdn<=PB<=PBup, pdB1dn<=pdB1<=pdB1up, 0<=drupB<=drscale*pdB1, 0<=drdnB<=drscale*pdB1];
   for i = 1:dbusB-1
       if i == 1
           CDB = [CDB, PB(4,:) == PB(1,:) - pdB(2,:) - pdB1(2,:)];
       elseif i == 2
           CDB = [CDB, PB(3,:) == PB(2,:) - pdB(5,:) - pdB1(5,:)];
       elseif i == 3
           CDB = [CDB, 0 == PB(3,:) - pdB(6,:) - pdB1(6,:)];
       elseif i == 4
           CDB = [CDB, PB(5,:) == PB(4,:) - pdB(3,:) - pdB1(3,:)];
       elseif i == 5
           CDB = [CDB, 0 == PB(5,:) - pdB(4,:) - pdB1(4,:)];
       end
   end 
   CDB = [CDB,PB(1,:)+PB(2,:) >= sum(pdB + pdB1)];
 %% Ограничения микросети 3
   CDE = [PEdn<=PE<=PEup, pdE1dn<=pdE1<=pdE1up, 0<=drupE<=drscale*pdE1, 0<=drdnE<=drscale*pdE1];
   for i = 1:dbusE-1
       if i == 1
           CDE = [CDE, PE(4,:) == PE(1,:) - pdE(2,:) - pdE1(2,:)];
       elseif i == 2
           CDE = [CDE, PE(3,:) == PE(2,:) - pdE(5,:) - pdE1(5,:)];
       elseif i == 3
           CDE = [CDE, 0 == PE(3,:) - pdE(6,:) - pdE1(6,:)];
       elseif i == 4
           CDE = [CDE, PE(5,:) == PE(4,:) - pdE(3,:) - pdE1(3,:)];
       elseif i == 5
           CDE = [CDE, 0 == PE(5,:) - pdE(4,:) - pdE1(4,:)];
       end
   end 
   CDE = [CDE,PE(1,:)+PE(2,:) >= sum(pdE + pdE1)];
%% Общая целевая функция для двухуровневой ЭЭС  
%     OO = sum(onoff')*Conoff'+sum(pg')*cg1+sum(rgup' + rgdn')*crg+sum(windup-sum(rgdn)-sum(drupA)-sum(drupB)-sum(drupE))*wc+...
%         sum(-winddown-sum(rgup)-sum(drdnA)-sum(drdnB)-sum(drdnE))*lc;%note
    OO = sum(onoff')*Conoff'+sum(pg')*cg1+sum(rgup' + rgdn')*crg+sum(windup-sum(rgdn)-sum(drupE)-sum(drupB)-sum(drupA))*wc+...
    sum(-winddown-sum(rgup)-sum(drdnE)-sum(drdnB)-sum(drdnA))*lc;%note
    O1 = sum(pg'.*pg')*cg2;
    O2 = sum(sum((pdA1up-pdA1).*(pdA1up-pdA1)))*CpdA1 + sum(sum(drupA'+drdnA'))*drA1 + sum(sum(drupA'.*drupA'+drdnA'.*drdnA'))*drA2
    O3 = sum(sum((pdB1up-pdB1).*(pdB1up-pdB1)))*CpdB1 + sum(sum(drupB'+drdnB'))*drB1 + sum(sum(drupB'.*drupB'+drdnB'.*drdnB'))*drB2
    O4 = sum(sum((pdE1up-pdE1).*(pdE1up-pdE1)))*CpdE1 + sum(sum(drupE'+drdnE'))*drE1 + sum(sum(drupE'.*drupE'+drdnE'.*drdnE'))*drE2
%% Решение оптимизационной проблемы квадратичным решателем Yalmip 
%    optimize([CDA,CDB,CDE,CO,EP,CP],OO+O1+O2+03+04) % note -10000*scale +O3+O4
   optimize([CDE,CDB,CDA,CO,EP,CP],OO+O1+O4+O3+O2) % note -10000*scale +O3+O4
%    value(cimA)
%    value(drpA)
%    value(dgA)
%    value(drupA)
%    value(drdnA)
   EnergyPrice = dual(EP)';
   ConPrice = dual(CP);
%% Отдельной целевые функции микросетей с оптимизированными параметрами 
   OD1 = sum(sum(drupA'+drdnA'))*drA1 + sum(sum(drupA'.*drupA'+drdnA'.*drdnA'))*drA2 + sum(sum((pdA1up-pdA1).*(pdA1up-pdA1)))*CpdA1 - sum(drupA+drdnA)*EnergyPrice'+ (PA(1,:)+ PA(2,:))*EnergyPrice' 
   OD2 = sum(sum(drupB'+drdnB'))*drB1 + sum(sum(drupB'.*drupB'+drdnB'.*drdnB'))*drB2 + sum(sum((pdB1up-pdB1).*(pdB1up-pdB1)))*CpdB1 -...
         sum(drupB+drdnB)*EnergyPrice'+ (PB(1,:)+ PB(2,:))*EnergyPrice'
   OD3 = sum(sum(drupE'+drdnE'))*drE1 + sum(sum(drupE'.*drupE'+drdnE'.*drdnE'))*drE2 + sum(sum((pdE1up-pdE1).*(pdE1up-pdE1)))*CpdE1 -...
         sum(drupE+drdnE)*EnergyPrice'+ (PE(1,:)+ PE(2,:))*EnergyPrice' 

%% Полученные затраты
disp('Затраты передающей сети')
value(OO+O1+O4+O3+O2) 
disp('Затраты распределенных микросетей')
value(OD1+OD2+OD3) 
disp('Суммарные затраты по всей двухуровневой ЭЭС')
value(OO+O1+O4+O3+O2+OD1+OD2+OD3) 

Pene1 = max(value(wf))/149 % исходный уровень проникновения ВЭС
Pene2 = (max(value(windup))+max(value(wf)))/149 % повышенный уровень проникновения ВЭС после оптимизации 

   disp('Суммарно верхняя граница отклонений от прогноза ветромощности')
   sum(windup)
   disp('Суммарно нижняя граница отклонений от прогноза ветромощности')
   sum(winddown)
   disp('Цена импортируемой электроэнергии')
   EnergyPrice = dual(EP)'

   
