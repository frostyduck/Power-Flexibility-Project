    clc;
    clear all;
    close all;
    parpool; % включаем параллельное вычисления для ускорения расчёта
%% Параметры передающей сети
    nt=24; % time horizon
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
   drscaleA = 0.5; %коэффициент масштабирования активной нагрузки (сейчас это 50%)
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
    pr2=0;%lower percentile(throughout)
    windup1=prctile(winderror,pr1); % верхняя граница отклонений прогноза ветромощности 
    winddown1=prctile(winderror,pr2); % нижняя граница отклонений прогноза ветромощности 

    % ветряк 2
    wl=0; % коэффициент масштабирования проникновения ветровой энергии
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
    scale = 0.5;
    wf1 = scale*(wf1);
    windup1 = scale*(windup1);
    winddown1 = scale*(winddown1);
    
    windup1 = 1*[2.14850277468970,2.38176933526874,1.74677853896453,2.60511542849926,2.35355880939503,2.41625640922597,9.26333481098347,4.72603747192732,3.64443832323432,3.54370589728489,3.77247031861537,5.32507128646877,3.89024919339414,3.43366439452552,3.79975232769689,5.91542122062703,18.0071565379186,20.3629232486021,27.8562658759252,16.9300600310706,14.6166013530343,12.8505562939878,9.60886741793157,16.5816874966563];
    winddown1 = 1*[-2.29653290221721,-2.29980619610469,-1.86455639333191,-1.35000000000000,-1,-1,-3.40000000000000,-1.45000000000000,-1.05000000000000,-1,-1.10000000000000,-1.55000000000000,-1.15000000000000,-1,-1.10000000000000,-1.95000000000000,-5.70000000000000,-6.30000000000000,-9.10000000000000,-5.65000000000000,-4.55000000000000,-3.90000000000000,-3.15000000000000,-5.15000000000000];

    % ветряк 2 
    scale2 = 0;
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
    
     % Ограничения резервов генератора 
    CO = [CO,Gmin.*onoff<=pg<=Gmax.*onoff]; % Генератор граничит с коэффицентом рампинга Rgs
    CO = [CO,pg+rgup<=Gmax.*onoff,Gmin.*onoff<=pg-rgdn]; 
    CO = [CO,(pg(:,2:nt)+rgup(:,2:nt))-(pg(:,1:nt-1)-rgdn(:,1:nt-1))<=Rup(:,2:nt)]; % верхний предел коэффицента рампинга генератора 
    CO = [CO,-Rdn(:,2:nt)<=(pg(:,2:nt)-rgdn(:,2:nt))-(pg(:,1:nt-1)+rgup(:,1:nt-1))]; % нижний предел коэффицента рампинга генератора
    % Ограничения генератора
    CO=[CO,-Rdn(:,2:nt)<=pg(:,2:nt)-pg(:,1:nt-1)<=Rup(:,2:nt)]; % рампинг добавляем в CO
%% Процесс редукции двухуровневой задачи в одноуровневую  
% Ограничения для редукции  
   CD = [];
   bigM = [];
   ST = [];
   DF = [];
   O2 = 0;
% Ограничения микросети 1 
  CO = [CO,Pinj(dbus1,:)==-loads(dbus1,:)-PA(1,:)];%
  CO = [CO,Pinj1(dbus1,:)==-loads(dbus1,:)-PA(1,:)+sum(drdnA)];%
  CO = [CO,Pinj2(dbus1,:)==-loads(dbus1,:)-PA(1,:)-sum(drupA)];%
  CD = [CD, dgAdn<=dgA<=dgAup, PAdn<=PA<=PAup, pdA1dn<=pdA1<=pdA1up, 0<=drupA<=drscaleA*pdA1,...
      0<=drdnA<=drscaleA*pdA1];
   for i = 1:dbusA-1
       if i == 1
           CD = [CD, PA(4,:) == PA(1,:) - pdA(2,:) - pdA1(2,:)]; %mu1A(1,:)
       elseif i == 2
           CD = [CD, PA(3,:) == PA(2,:) - pdA(5,:) - pdA1(5,:)];
       elseif i == 3
           CD = [CD, 0 == PA(3,:) - pdA(6,:) - pdA1(6,:)];
       elseif i == 4
           CD = [CD, PA(5,:) == PA(4,:) - pdA(3,:) - pdA1(3,:)];
       elseif i == 5
           CD = [CD, 0 == PA(5,:) - pdA(4,:) - pdA1(4,:)]; %mu1A(5,:)
       end
   end  
   
   mu1A = sdpvar(5,nt,'full'); 
   l1A = sdpvar(1,nt,'full'); 
   b1A = binvar(1,nt,'full'); 
   m1A = 10000000;
%    DC = [DC,l1A.*(dgA-dgAdn) == 0];
   bigM = [bigM, l1A<=m1A*b1A, dgA-dgAdn<=m1A*(1-b1A)];
   
   l2A = sdpvar(1,nt,'full');  
   b2A = binvar(1,nt,'full');
%    DC = [DC,l2A.*(dgA-dgAup) == 0];
   bigM = [bigM, l2A<=m1A*b2A, -dgA+dgAup<=m1A*(1-b2A)];
   
   l3A = sdpvar(dbusA,nt,'full'); 
   b3A = binvar(dbusA,nt,'full'); 
%    DC = [DC,l3A.*(PA-PAdn) == 0];
   bigM = [bigM, l3A<=m1A*b3A, PA-PAdn<=m1A*(1-b3A)];
   
   l4A = sdpvar(dbusA,nt,'full');  
   b4A = binvar(dbusA,nt,'full');
%    DC = [DC,l4A.*(PA-PAup) == 0];
   bigM = [bigM, l4A<=m1A*b4A, -PA+PAup<=m1A*(1-b4A)];
   
   l5A = sdpvar(dbusA,nt,'full'); 
   b5A = binvar(dbusA,nt,'full'); 
%    DC = [DC,l5A.*(pdA1-pdA1dn) == 0];
   bigM = [bigM, l5A<=m1A*b5A, pdA1-pdA1dn<=m1A*(1-b5A)];
   
   l6A = sdpvar(dbusA,nt,'full');  
   b6A = binvar(dbusA,nt,'full');
%    DC = [DC,l6A.*(pdA1-pdA1up) == 0];
   bigM = [bigM, l6A<=m1A*b6A, -pdA1+pdA1up<=m1A*(1-b6A)];
   
   l7A = sdpvar(dbusA,nt,'full'); 
   b7A = binvar(dbusA,nt,'full'); 
%    DC = [DC,l7A.*(drupA) == 0];
   bigM = [bigM, l7A<=m1A*b7A, drupA<=m1A*(1-b7A)];
   
   l8A = sdpvar(dbusA,nt,'full');  
   b8A = binvar(dbusA,nt,'full');
%    DC = [DC,l8A.*(drupA-drscale*pdA1) == 0];
   bigM = [bigM, l8A<=m1A*b8A, -drupA+drscaleA*pdA1<=m1A*(1-b8A)];
   
   l9A = sdpvar(dbusA,nt,'full'); 
   b9A = binvar(dbusA,nt,'full'); 
%    DC = [DC,l9A.*(drdnA) == 0];
   bigM = [bigM, l9A<=m1A*b9A, drdnA<=m1A*(1-b9A)];
   
   l10A = sdpvar(dbusA,nt,'full');  
   b10A = binvar(dbusA,nt,'full');
%    DC = [DC,l10A.*(drdnA-drscale*pdA1) == 0];
   bigM = [bigM, l10A<=m1A*b10A, -drdnA+drscaleA*pdA1<=m1A*(1-b10A)];
   
   
  %% Задаём условия Каруша-Куна-Такера для микросетей  
   % Стационарность описывает набор условий оптимальности первого порядка, 
   % то есть первая производная функции Лагранжа по каждой переменной решения равна 0 (микросеть 1).
   
   DF = [DF,l1A>=0,l2A>=0,l3A>=0,l4A>=0,l5A>=0,l6A>=0,l7A>=0,l8A>=0,l9A>=0,l10A>=0];
   ST = [ST,2*cd2A*dgA+cd1A*ones(1,nt)+l2A-l1A-mu1A(5,:) == 0];%dgA
   ST = [ST,cimA-mu1A(1,:)+l4A(1,:)-l3A(1,:)==0, cimA-mu1A(2,:)+l4A(2,:)-l3A(2,:) == 0, mu1A(2,:)-mu1A(3,:)+l4A(3,:)-l3A(3,:) == 0, ...
       mu1A(1,:)-mu1A(4,:)+l4A(4,:)-l3A(4,:) == 0, mu1A(4,:)-mu1A(5,:)+l4A(5,:)-l3A(5,:) == 0, mu1A(5,:)+l4A(6,:)-l3A(6,:) == 0];%PA
   ST = [ST,2*CpdA1*pdA1(1,:)-2*CpdA1*pdA1up(1,:)+l6A(1,:)-l5A(1,:)-drscaleA*l8A(1,:)-drscaleA*l10A(1,:) == 0,...
       2*CpdA1*pdA1(2,:)-2*CpdA1*pdA1up(2,:)+mu1A(1,:)+l6A(2,:)-l5A(2,:)-drscaleA*l8A(2,:)-drscaleA*l10A(2,:) == 0,...
       2*CpdA1*pdA1(3,:)-2*CpdA1*pdA1up(3,:)+mu1A(4,:)+l6A(3,:)-l5A(3,:)-drscaleA*l8A(3,:)-drscaleA*l10A(3,:) == 0,...
       2*CpdA1*pdA1(4,:)-2*CpdA1*pdA1up(4,:)+mu1A(5,:)+l6A(4,:)-l5A(4,:)-drscaleA*l8A(4,:)-drscaleA*l10A(4,:) == 0,...
       2*CpdA1*pdA1(5,:)-2*CpdA1*pdA1up(5,:)+mu1A(2,:)+l6A(5,:)-l5A(5,:)-drscaleA*l8A(5,:)-drscaleA*l10A(5,:) == 0,...
       2*CpdA1*pdA1(6,:)-2*CpdA1*pdA1up(6,:)+mu1A(3,:)+l6A(6,:)-l5A(6,:)-drscaleA*l8A(6,:)-drscaleA*l10A(6,:) == 0];%pdA1      
   ST = [ST,2*drA2*drupA+(drA1)*ones(dbusA,nt)-repmat(drpA,dbusA,1)+l8A-l7A == 0,2*drA2*drdnA+(drA1)*ones(dbusA,nt)-repmat(drpA,dbusA,1)+l10A-l9A == 0];%drup,drdn
   ST = [ST,cimA == drpA];%dgA
   % Заменяем целевую функцию распределительной сети через двойственные оценки
   dualA = -sum(dgA.*dgA)*cd2A- sum(sum(drupA'.*drupA'+drdnA'.*drdnA'))*drA2-sum(sum(pdA1.*pdA1)*CpdA1)+...
       sum(sum(pdA(2:6,:).*mu1A)+sum(l3A*PAdn-l4A.*PAup)+sum(l5A.*pdA1dn-l6A.*pdA1up)-l2A*dgAup) + sum(sum(pdA1up.*pdA1up))*CpdA1
   OD1 = sum(sum(drupA'+drdnA'))*drA1 + sum(sum(drupA'.*drupA'+drdnA'.*drdnA'))*drA2 +sum(dgA')*cd1A +...
   sum(dgA'.*dgA')*cd2A  + sum(sum((pdA1up-pdA1).*(pdA1up-pdA1)))*CpdA1-sum(drpA.*sum(drupA+drdnA))+ sum((PA(1,:)+PA(2,:)).*cimA) 
   
   O2 = O2 + dualA-(sum(sum(drupA'+drdnA'))*drA1 + sum(sum(drupA'.*drupA'+drdnA'.*drdnA'))*drA2 +sum(dgA')*cd1A...
       + sum(dgA'.*dgA')*cd2A + sum(sum((pdA1up-pdA1).*(pdA1up-pdA1)))*CpdA1)

 % Ограничения микросети 2 
  CO = [CO,Pinj(dbus2,:)==-loads(dbus2,:)-PB(1,:)];%note
  CO = [CO,Pinj1(dbus2,:)==-loads(dbus2,:)-PB(1,:)+sum(drdnB)];%note
  CO = [CO,Pinj2(dbus2,:)==-loads(dbus2,:)-PB(1,:)-sum(drupB)];%note
  CD = [CD, dgBdn<=dgB<=dgBup, PBdn<=PB<=PBup, pdB1dn<=pdB1<=pdB1up, 0<=drupB<=drscaleB*pdB1,...
      0<=drdnB<=drscaleB*pdB1];
   for i = 1:dbusB-1
       if i == 1
           CD = [CD, PB(4,:) == PB(1,:) - pdB(2,:) - pdB1(2,:)]; %mu1B(1,:)
       elseif i == 2
           CD = [CD, PB(3,:) == PB(2,:) - pdB(5,:) - pdB1(5,:)];
       elseif i == 3
           CD = [CD, 0 == PB(3,:) - pdB(6,:) - pdB1(6,:)];
       elseif i == 4
           CD = [CD, PB(5,:) == PB(4,:) - pdB(3,:) - pdB1(3,:)];
       elseif i == 5
           CD = [CD, 0 == PB(5,:) - pdB(4,:) - pdB1(4,:)]; %mu1B(5,:)
       end
   end  
   
   mu1B = sdpvar(5,nt,'full'); 
   l1B = sdpvar(1,nt,'full'); 
   b1B = binvar(1,nt,'full'); 
   m1B = 10000000;
%    DC = [DC,l1B.*(dgB-dgBdn) == 0];
   bigM = [bigM, l1B<=m1B*b1B, dgB-dgBdn<=m1B*(1-b1B)];
   
   l2B = sdpvar(1,nt,'full');  
   b2B = binvar(1,nt,'full');
%    DC = [DC,l2B.*(dgB-dgBup) == 0];
   bigM = [bigM, l2B<=m1B*b2B, -dgB+dgBup<=m1B*(1-b2B)];
   
   l3B = sdpvar(dbusB,nt,'full'); 
   b3B = binvar(dbusB,nt,'full'); 
%    DC = [DC,l3B.*(PB-PBdn) == 0];
   bigM = [bigM, l3B<=m1B*b3B, PB-PBdn<=m1B*(1-b3B)];
   
   l4B = sdpvar(dbusB,nt,'full');  
   b4B = binvar(dbusB,nt,'full');
%    DC = [DC,l4B.*(PB-PBup) == 0];
   bigM = [bigM, l4B<=m1B*b4B, -PB+PBup<=m1B*(1-b4B)];
   
   l5B = sdpvar(dbusB,nt,'full'); 
   b5B = binvar(dbusB,nt,'full'); 
%    DC = [DC,l5B.*(pdB1-pdB1dn) == 0];
   bigM = [bigM, l5B<=m1B*b5B, pdB1-pdB1dn<=m1B*(1-b5B)];
   
   l6B = sdpvar(dbusB,nt,'full');  
   b6B = binvar(dbusB,nt,'full');
%    DC = [DC,l6B.*(pdB1-pdB1up) == 0];
   bigM = [bigM, l6B<=m1B*b6B, -pdB1+pdB1up<=m1B*(1-b6B)];
   
   l7B = sdpvar(dbusB,nt,'full'); 
   b7B = binvar(dbusB,nt,'full'); 
%    DC = [DC,l7B.*(drupB) == 0];
   bigM = [bigM, l7B<=m1B*b7B, drupB<=m1B*(1-b7B)];
   
   l8B = sdpvar(dbusB,nt,'full');  
   b8B = binvar(dbusB,nt,'full');
%    DC = [DC,l8B.*(drupB-drscale*pdB1) == 0];
   bigM = [bigM, l8B<=m1B*b8B, -drupB+drscaleB*pdB1<=m1B*(1-b8B)];
   
   l9B = sdpvar(dbusB,nt,'full'); 
   b9B = binvar(dbusB,nt,'full'); 
%    DC = [DC,l9B.*(drdnB) == 0];
   bigM = [bigM, l9B<=m1B*b9B, drdnB<=m1B*(1-b9B)];
   
   l10B = sdpvar(dbusB,nt,'full');  
   b10B = binvar(dbusB,nt,'full');
%    DC = [DC,l10B.*(drdnB-drscale*pdB1) == 0];
   bigM = [bigM, l10B<=m1B*b10B, -drdnB+drscaleB*pdB1<=m1B*(1-b10B)];

   % Стационарность описывает набор условий оптимальности первого порядка, 
   % то есть первая производная функции Лагранжа по каждой переменной решения равна 0 (микросеть 2).
   DF = [DF,l1B>=0,l2B>=0,l3B>=0,l4B>=0,l5B>=0,l6B>=0,l7B>=0,l8B>=0,l9B>=0,l10B>=0];
   ST = [ST,2*cd2B*dgB+cd1B*ones(1,nt)+l2B-l1B-mu1B(5,:) == 0];%dgB
   ST = [ST,cimB-mu1B(1,:)+l4B(1,:)-l3B(1,:)==0, cimB-mu1B(2,:)+l4B(2,:)-l3B(2,:) == 0, mu1B(2,:)-mu1B(3,:)+l4B(3,:)-l3B(3,:) == 0, ...
       mu1B(1,:)-mu1B(4,:)+l4B(4,:)-l3B(4,:) == 0, mu1B(4,:)-mu1B(5,:)+l4B(5,:)-l3B(5,:) == 0, mu1B(5,:)+l4B(6,:)-l3B(6,:) == 0];%PB
   ST = [ST,2*CpdB1*pdB1(1,:)-2*CpdB1*pdB1up(1,:)+l6B(1,:)-l5B(1,:)-drscaleB*l8B(1,:)-drscaleB*l10B(1,:) == 0,...
       2*CpdB1*pdB1(2,:)-2*CpdB1*pdB1up(2,:)+mu1B(1,:)+l6B(2,:)-l5B(2,:)-drscaleB*l8B(2,:)-drscaleB*l10B(2,:) == 0,...
       2*CpdB1*pdB1(3,:)-2*CpdB1*pdB1up(3,:)+mu1B(4,:)+l6B(3,:)-l5B(3,:)-drscaleB*l8B(3,:)-drscaleB*l10B(3,:) == 0,...
       2*CpdB1*pdB1(4,:)-2*CpdB1*pdB1up(4,:)+mu1B(5,:)+l6B(4,:)-l5B(4,:)-drscaleB*l8B(4,:)-drscaleB*l10B(4,:) == 0,...
       2*CpdB1*pdB1(5,:)-2*CpdB1*pdB1up(5,:)+mu1B(2,:)+l6B(5,:)-l5B(5,:)-drscaleB*l8B(5,:)-drscaleB*l10B(5,:) == 0,...
       2*CpdB1*pdB1(6,:)-2*CpdB1*pdB1up(6,:)+mu1B(3,:)+l6B(6,:)-l5B(6,:)-drscaleB*l8B(6,:)-drscaleB*l10B(6,:) == 0];%pdB1      
   ST = [ST,2*drB2*drupB+(drB1)*ones(dbusB,nt)-repmat(drpB,dbusB,1)+l8B-l7B == 0,2*drB2*drdnB+(drB1)*ones(dbusB,nt)-repmat(drpB,dbusB,1)+l10B-l9B == 0];%drup,drdn
   ST = [ST,cimB == drpB];%dgB
   % Заменяем целевую функцию распределительной сети через двойственные оценки
   dualB = -sum(dgB.*dgB)*cd2B- sum(sum(drupB'.*drupB'+drdnB'.*drdnB'))*drB2-sum(sum(pdB1.*pdB1)*CpdB1)+...
       sum(sum(pdB(2:6,:).*mu1B)+sum(l3B*PBdn-l4B.*PBup)+sum(l5B.*pdB1dn-l6B.*pdB1up)-l2B*dgBup) + sum(sum(pdB1up.*pdB1up))*CpdB1
   OD1B = sum(sum(drupB'+drdnB'))*drB1 + sum(sum(drupB'.*drupB'+drdnB'.*drdnB'))*drB2 +sum(dgB')*cd1B +...
   sum(dgB'.*dgB')*cd2B  + sum(sum((pdB1up-pdB1).*(pdB1up-pdB1)))*CpdB1-sum(drpB.*sum(drupB+drdnB))+ sum((PB(1,:)+PB(2,:)).*cimB) 
   
   O2 = O2 + dualB-(sum(sum(drupB'+drdnB'))*drB1 + sum(sum(drupB'.*drupB'+drdnB'.*drdnB'))*drB2 +sum(dgB')*cd1B...
       + sum(dgB'.*dgB')*cd2B + sum(sum((pdB1up-pdB1).*(pdB1up-pdB1)))*CpdB1)
  
   % Ограничения микросети 3
  CO = [CO,Pinj(dbus3,:)==-loads(dbus3,:)-PE(1,:)];%note
  CO = [CO,Pinj1(dbus3,:)==-loads(dbus3,:)-PE(1,:)+sum(drdnE)];%note
  CO = [CO,Pinj2(dbus3,:)==-loads(dbus3,:)-PE(1,:)-sum(drupE)];%note
  CD = [CD, dgEdn<=dgE<=dgEup, PEdn<=PE<=PEup, pdE1dn<=pdE1<=pdE1up, 0<=drupE<=drscaleE*pdE1,...
      0<=drdnE<=drscaleE*pdE1];
   for i = 1:dbusE-1
       if i == 1
           CD = [CD, PE(4,:) == PE(1,:) - pdE(2,:) - pdE1(2,:)]; %mu1E(1,:)
       elseif i == 2
           CD = [CD, PE(3,:) == PE(2,:) - pdE(5,:) - pdE1(5,:)];
       elseif i == 3
           CD = [CD, 0 == PE(3,:) - pdE(6,:) - pdE1(6,:)];
       elseif i == 4
           CD = [CD, PE(5,:) == PE(4,:) - pdE(3,:) - pdE1(3,:)];
       elseif i == 5
           CD = [CD, 0 == PE(5,:) - pdE(4,:) - pdE1(4,:)]; %mu1E(5,:)
       end
   end  
   
   mu1E = sdpvar(5,nt,'full'); 
   l1E = sdpvar(1,nt,'full'); 
   b1E = binvar(1,nt,'full'); 
   m1E = 10000000;
%    DC = [DC,l1E.*(dgE-dgEdn) == 0];
   bigM = [bigM, l1E<=m1E*b1E, dgE-dgEdn<=m1E*(1-b1E)];
   
   l2E = sdpvar(1,nt,'full');  
   b2E = binvar(1,nt,'full');
%    DC = [DC,l2E.*(dgE-dgEup) == 0];
   bigM = [bigM, l2E<=m1E*b2E, -dgE+dgEup<=m1E*(1-b2E)];
   
   l3E = sdpvar(dbusE,nt,'full'); 
   b3E = binvar(dbusE,nt,'full'); 
%    DC = [DC,l3E.*(PE-PEdn) == 0];
   bigM = [bigM, l3E<=m1E*b3E, PE-PEdn<=m1E*(1-b3E)];
   
   l4E = sdpvar(dbusE,nt,'full');  
   b4E = binvar(dbusE,nt,'full');
%    DC = [DC,l4E.*(PE-PEup) == 0];
   bigM = [bigM, l4E<=m1E*b4E, -PE+PEup<=m1E*(1-b4E)];
   
   l5E = sdpvar(dbusE,nt,'full'); 
   b5E = binvar(dbusE,nt,'full'); 
%    DC = [DC,l5E.*(pdE1-pdE1dn) == 0];
   bigM = [bigM, l5E<=m1E*b5E, pdE1-pdE1dn<=m1E*(1-b5E)];
   
   l6E = sdpvar(dbusE,nt,'full');  
   b6E = binvar(dbusE,nt,'full');
%    DC = [DC,l6E.*(pdE1-pdE1up) == 0];
   bigM = [bigM, l6E<=m1E*b6E, -pdE1+pdE1up<=m1E*(1-b6E)];
   
   l7E = sdpvar(dbusE,nt,'full'); 
   b7E = binvar(dbusE,nt,'full'); 
%    DC = [DC,l7E.*(drupE) == 0];
   bigM = [bigM, l7E<=m1E*b7E, drupE<=m1E*(1-b7E)];
   
   l8E = sdpvar(dbusE,nt,'full');  
   b8E = binvar(dbusE,nt,'full');
%    DC = [DC,l8E.*(drupE-drscale*pdE1) == 0];
   bigM = [bigM, l8E<=m1E*b8E, -drupE+drscaleE*pdE1<=m1E*(1-b8E)];
   
   l9E = sdpvar(dbusE,nt,'full'); 
   b9E = binvar(dbusE,nt,'full'); 
%    DC = [DC,l9E.*(drdnE) == 0];
   bigM = [bigM, l9E<=m1E*b9E, drdnE<=m1E*(1-b9E)];
   
   l10E = sdpvar(dbusE,nt,'full');  
   b10E = binvar(dbusE,nt,'full');
%    DC = [DC,l10E.*(drdnE-drscale*pdE1) == 0];
   bigM = [bigM, l10E<=m1E*b10E, -drdnE+drscaleE*pdE1<=m1E*(1-b10E)];
 % Стационарность описывает набор условий оптимальности первого порядка, 
   % то есть первая производная функции Лагранжа по каждой переменной решения равна 0 (микросеть 3).
   DF = [DF,l1E>=0,l2E>=0,l3E>=0,l4E>=0,l5E>=0,l6E>=0,l7E>=0,l8E>=0,l9E>=0,l10E>=0];
   ST = [ST,2*cd2E*dgE+cd1E*ones(1,nt)+l2E-l1E-mu1E(5,:) == 0];%dgE
   ST = [ST,cimE-mu1E(1,:)+l4E(1,:)-l3E(1,:)==0, cimE-mu1E(2,:)+l4E(2,:)-l3E(2,:) == 0, mu1E(2,:)-mu1E(3,:)+l4E(3,:)-l3E(3,:) == 0, ...
       mu1E(1,:)-mu1E(4,:)+l4E(4,:)-l3E(4,:) == 0, mu1E(4,:)-mu1E(5,:)+l4E(5,:)-l3E(5,:) == 0, mu1E(5,:)+l4E(6,:)-l3E(6,:) == 0];%PE
   ST = [ST,2*CpdE1*pdE1(1,:)-2*CpdE1*pdE1up(1,:)+l6E(1,:)-l5E(1,:)-drscaleE*l8E(1,:)-drscaleE*l10E(1,:) == 0,...
       2*CpdE1*pdE1(2,:)-2*CpdE1*pdE1up(2,:)+mu1E(1,:)+l6E(2,:)-l5E(2,:)-drscaleE*l8E(2,:)-drscaleE*l10E(2,:) == 0,...
       2*CpdE1*pdE1(3,:)-2*CpdE1*pdE1up(3,:)+mu1E(4,:)+l6E(3,:)-l5E(3,:)-drscaleE*l8E(3,:)-drscaleE*l10E(3,:) == 0,...
       2*CpdE1*pdE1(4,:)-2*CpdE1*pdE1up(4,:)+mu1E(5,:)+l6E(4,:)-l5E(4,:)-drscaleE*l8E(4,:)-drscaleE*l10E(4,:) == 0,...
       2*CpdE1*pdE1(5,:)-2*CpdE1*pdE1up(5,:)+mu1E(2,:)+l6E(5,:)-l5E(5,:)-drscaleE*l8E(5,:)-drscaleE*l10E(5,:) == 0,...
       2*CpdE1*pdE1(6,:)-2*CpdE1*pdE1up(6,:)+mu1E(3,:)+l6E(6,:)-l5E(6,:)-drscaleE*l8E(6,:)-drscaleE*l10E(6,:) == 0];%pdE1      
   ST = [ST,2*drE2*drupE+(drE1)*ones(dbusE,nt)-repmat(drpE,dbusE,1)+l8E-l7E == 0,2*drE2*drdnE+(drE1)*ones(dbusE,nt)-repmat(drpE,dbusE,1)+l10E-l9E == 0];%drup,drdn
   ST = [ST,cimE == drpE];%dgE
% Заменяем целевую функцию распределительной сети через двойственные оценки
   dualE = -sum(dgE.*dgE)*cd2E- sum(sum(drupE'.*drupE'+drdnE'.*drdnE'))*drE2-sum(sum(pdE1.*pdE1)*CpdE1)+...
       sum(sum(pdE(2:6,:).*mu1E)+sum(l3E*PEdn-l4E.*PEup)+sum(l5E.*pdE1dn-l6E.*pdE1up)-l2E*dgEup) + sum(sum(pdE1up.*pdE1up))*CpdE1
   OD1E = sum(sum(drupE'+drdnE'))*drE1 + sum(sum(drupE'.*drupE'+drdnE'.*drdnE'))*drE2 +sum(dgE')*cd1E +...
   sum(dgE'.*dgE')*cd2E  + sum(sum((pdE1up-pdE1).*(pdE1up-pdE1)))*CpdE1-sum(drpE.*sum(drupE+drdnE))+ sum((PE(1,:)+PE(2,:)).*cimE) 
   
   O2 = O2 + dualE-(sum(sum(drupE'+drdnE'))*drE1 + sum(sum(drupE'.*drupE'+drdnE'.*drdnE'))*drE2 +sum(dgE')*cd1E...
       + sum(dgE'.*dgE')*cd2E + sum(sum((pdE1up-pdE1).*(pdE1up-pdE1)))*CpdE1)
%%
    CO = [CO,-lineC<=GSF*Pinj2<=lineC];
    CO = [CO,-lineC<=GSF*Pinj1<=lineC];
    CO = [CO,-lineC<=GSF*Pinj<=lineC];
    CO = [CO,sum(pg)-sum(loads)+wf-PB(1,:)-PB(2,:)-PA(1,:)-PA(2,:)-PE(1,:)-PE(2,:)==0];   %note
    CO = [CO,windup-sum(rgdn)-sum(drupB)-sum(drupA)-sum(drupE)==0,-winddown-sum(rgup)-sum(drdnB)-sum(drdnA)-sum(drdnE)==0]; % note
%% Целевая функция системы передачи 
    OO = sum(onoff')*Conoff'+sum(pg')*cg1+sum(rgup' + rgdn')*crg+sum(windup-sum(rgdn)-sum(drupB)-sum(drupA)-sum(drupE))*wc...
        +sum(-winddown-sum(rgup)-sum(drdnB)-sum(drdnA)-sum(drdnE))*lc;%
    O1 = sum(pg'.*pg')*cg2;
%% Решение оптимизационной проблемы квадратичным решателем Yalmip 
   optimize([CD,ST,bigM,DF,CO],OO+O1-O2) % -10000*scale
%    value(cimA)
%    value(drpA)
%    value(dgA)
%    value(drupA)
%    value(drdnA)

%% Полученные затраты
disp('Затраты передающей сети')
value((OO+O1-O2)/1000000) 
disp('Затраты распределенных микросетей')
value(OD1+OD1B+OD1E) 
disp('Суммарные затраты по всей двухуровневой ЭЭС')
value((OO+O1-O2)/1000000 +OD1+OD1B+OD1E) 

Pene1 = max(value(wf))/149 % исходный уровень проникновения ВЭС
Pene2 = (max(value(windup))+max(value(wf)))/149 % повышенный уровень проникновения ВЭС после оптимизации 

   disp('Суммарно верхняя граница отклонений от прогноза ветромощности')
   sum(windup)
   disp('Суммарно нижняя граница отклонений от прогноза ветромощности')
   sum(winddown)
   disp('Цена импортируемой электроэнергии')
   EnergyPrice = value(cimA)
   

 


