clc;
clear all;
close all;
    
%% Параметры передающей сети 
    nt=24; % временной горизонт
    % определим затраты 
    cg1=[2 1.75 1 3.25 3 3]'; % линейный  коэффициент затрат на генерацию 
    cg2=[0.02 0.0175 0.0625 0.00834 0.025 0.025]'; % квадратичный коэффициент затрат на генерацию 
    crgs=[6 6.75 7 5.25 5 5]'; % коэффициент цены запуска резерва генератора
    mpc=case30; % схема передающей сети 
    sampleL=[120.6 115.8 114.8 112.6 114.0 113.4 117.1 126.3 130.7 132.5 135.6 134.8 136.5 137.7 137.1 138.0 136.3 133.3 131.7 129.3 128.2 127.4 125.6 124.2]; % вектор нагрузки
    sf=sampleL(:,1:nt)./mean(sampleL(:,1:nt));
    sf=repmat(sf,30,1);
    loads=repmat(mpc.bus(:,3),1,nt); 
    loads=loads.*sf;
    mbus=8; % узел для подключения микросети
    wbus=8; % узел для подключения ветростанции 
    
    % сделаем Y/B узловую матрицу и матрицу A в Pline=A*Pinj
    [Ybus, Yf, Yt] = makeYbus(mpc);
    B=full(real(i.*Ybus));
    NB=-B;
    Bred=B(1:29,1:29); % уменьшенная матрица B 

    % извлекаем топологию сети 
    from = mpc.branch(:,1);
    to = mpc.branch(:,2);
    M = zeros(41,30); % от/к матрице
    lineC = zeros(41,1);
    for i=1:41
        M(i,from(i))=1;
        M(i,to(i))=-1;
        lineC(i)=1*mpc.branch(i,7);
    end
    m = M(:,1:29);

    x=mpc.branch(:,4); % реактиное сопротивление
    r=mpc.branch(:,3); % активное сопротивление
    b=x./(x.^2+r.^2); % реактивная проводимость
    D=diag(b); % диагональная матрица проводимости   
%% Параметры микросети 
    lmp=5.3*ones(1,nt); %вектор нагрузки микросети
    cg=[3.3 3.3 3.3]'; 
    cd=1.01*[3.3 3.3 3.3 3.3 3.3 3.3]'; % немного больше чем cg, чтобы поощрить затраты
    cb=0.1; % должен быть небольшим, чтобы микросеть запасала энергию, а не экспортировала энергию
    ci=lmp;% цена импортируемой микросетью электроэнергии
    ce=0.8*lmp;% цена экспортируемой микросетью электроэнергии
    pdr = 1*ones(1,24);% штраф управление спросом 
%% Прогноз ветромощности и моделирование сценариев
    %прогноз ветра 
    wl=0.01; % фактор масштабирования проникновения ветрогенерации 
    load('ObsDays.mat'); % датасет с ретроспективными данными ветрогенерации за 3 года 
    k=54; % количество траекторий ветрогенерации 
    wf=wl*ObsDays(k,:,11);
    corr1=corrcoef(ObsDays(:,:,11));
    std_dev = [0.12,0.15,0.18,0.5,0.6,0.67,0.72,0.76,0.79,0.82,0.83,0.8315,0.833,0.835,0.836,0.838,0.839,0.841,0.842,0.844,0.845,0.847,0.848,0.85];
    for j=1:nt
        for k=1:nt
        covarr1(j,k)=corr1(j,k)*std_dev(k)*std_dev(j);
        end
    end

    beta=0.01;
    eulernum=exp(1);
    epsilon=0.05;
    Ndelta=1*nt;%количество неопределенностей, один ветер в каждом периоде
    Nneed=ceil((1/epsilon)*(eulernum/(eulernum-1))*(log(1/beta)+4*(Ndelta+1)-1));% количество сценариев 

    %создаём набор ошибок ветра 
    avg=zeros(nt,1)';
    for i=1:Nneed
        err1(i,:)=mvnrnd(avg,covarr1);
    end

    for j=1:Nneed
        wn1(j,:)=wf+wf.*err1(j,:);
        pos=wn1(j,:)>=0;
        wst(j,:)=wn1(j,:).*pos;
    end
    winderror=wst-repmat(wf,Nneed,1);
    %создаём набор ошибок ветра 
    pr1=90;%верхний процентиль (на всем протяжении)
    pr2=10;%нижний процентиль (на всем протяжении)
    % Задаём верхнее и нижнее отклонения мощности ветростанции 
    windup=prctile(winderror,pr1); 
    winddown=prctile(winderror,pr2);
%% Параметры оптимизации передающей сети и её ограничения 
    % Генераторы
    ng=6; % количество генераторов
    Gmax=[80*ones(1,nt);80*ones(1,nt);40*ones(1,nt);50*ones(1,nt);30*ones(1,nt);55*ones(1,nt)];
    Gmin=[zeros(6,nt)];
    pg=sdpvar(ng,nt,'full');

    Rgsmin=-0.5*Gmax;
    Rgsmax=0.5*Gmax;
    rgs=sdpvar(ng,nt,'full');
    % Цены импорта и экспорта 
    Pim=sdpvar(1,nt,'full');
    Pex=sdpvar(1,nt,'full');
    % Цена управления спросом 
    Pdr=sdpvar(1,nt,'full');
    % Верхние и нижние пределы задаём 
    CO = [Gmin<=pg<=Gmax,Rgsmin<=rgs<=Rgsmax,0<=Pdr<=20,0<=Pex<=10,9<=Pim<=10];      
%% Параметры оптимизации микросети и её ограничения 
    % 1. Активные нагрузки 
    nd=1; % количество активных нагрузок 
    Pdmin=0.1*[0.5*ones(1,nt)+4*ones(1,nt)+2*ones(1,nt)+5.5*ones(1,nt)+1*ones(1,nt)+7*ones(1,nt)];
    Pdmax=0.03*[10*ones(1,nt)+16*ones(1,nt)+15*ones(1,nt)+20*ones(1,nt)+27*ones(1,nt)+32*ones(1,nt)];
    Pdb = 0.5*(Pdmin+Pdmax);
    pdup=sdpvar(nd,nt,'full');
    pddn=sdpvar(nd,nt,'full');

    % 2. Пассивные нагрузки 
    NCL=0.03*[120.6 115.8 114.8 112.6 114.0 113.4 117.1 126.3...
             130.7 132.5 135.6 136.5 137.7 133.7 137.1 138.0...
             136.3 133.3 131.7 129.3 128.2 127.4 125.6 124.2];     
    NCL=NCL(:,1:nt);

    % 3. Генераторы в микросети 
    ng=1; % количество генераторов 
    Gmmax=0.3*[5*ones(1,nt)+4.5*ones(1,nt)+7*ones(1,nt)];
    Gmmin=1*[1*ones(1,nt)+0.8*ones(1,nt)+1.5*ones(1,nt)];
    mpg=sdpvar(ng,nt,'full'); % генерация микросети 
%     mrgs=sdpvar(ng,nt,'full');

    % 4. Накопители энергии 
    nb=1; % количество накопителей 
    % мощность в и из накопителя 
    Pbmin=-3*ones(nb,nt); % граница разрядки накопителя
    Pbmax=3*ones(nb,nt); % граница зарядки накопителя 
    pb=sdpvar(nb,nt,'full'); % уровень зарядки накопителя

    % состояние разряда 
    Bmin=zeros(nb,nt);
    Bmax=10*ones(nb,nt);
    mb=sdpvar(nb,nt,'full');
    mb(1)=5;% определяем состояние разряда

    % 4. Импорт и экспорт мощности 
    ex=sdpvar(1,nt,'full');
    im=sdpvar(1,nt,'full');
    net=sdpvar(1,nt,'full');
    Nmin=-20*ones(1,nt);
    Nmax=20*ones(1,nt);

    % Верхние и нижние границы задаём 
    CI = [Pdmin<=pdup<=Pdmax,Pdmin<=pddn<=Pdmax,Gmmin<=mpg<=Gmmax,Pbmin<=pb<=Pbmax,Bmin<=mb<=Bmax,0<=ex,0<=im,Nmin<=net<=Nmax];  
%% Ограничения сети передачи 
    % Потокраспределение с прогнозом ветра 
    Pinj=sdpvar(30,nt,'full'); % узловая матрциа с прогнозом ветра 
    for i=1:30
      if i==1
      CO=[CO,Pinj(i,:)==-loads(i,:)+pg(1,:)];
      elseif i==2
      CO=[CO,Pinj(i,:)==-loads(i,:)+pg(2,:)];
      elseif i==13
      CO=[CO,Pinj(i,:)==-loads(i,:)+pg(3,:)];
      elseif i==22
      CO=[CO,Pinj(i,:)==-loads(i,:)+pg(4,:)];
      elseif i==23
      CO=[CO,Pinj(i,:)==-loads(i,:)+pg(5,:)];
      elseif i==27
      CO=[CO,Pinj(i,:)==-loads(i,:)+pg(6,:)];
%       elseif i==wbus
%       CO=[CO,Pinj(i,:)==-loads(i,:)+wst(1,:)];
      elseif i==mbus
      CO=[CO,Pinj(i,:)==-loads(i,:)+net+wf];
      else
      CO=[CO,Pinj(i,:)==-loads(i,:)];
      end
    end

    theta=sdpvar(29,nt,'full'); % угол для узла 
    CO = [CO,(theta==inv(Bred)*Pinj(1:29,:)):'lmp'];
    pflow=sdpvar(41,nt,'full'); % поток по линии 
    CO = [CO,pflow==D*m*theta];
    CO = [CO,-repmat(lineC,1,nt)<=pflow<=repmat(lineC,1,nt)];%  линия в ограничения CO   
 
    % Ограничения по резерву 
    CO=[CO,pdup-Pdb>=windup,pddn-Pdb<=winddown];    

    %Ограничения генератора 
    Rdn=0.3*Gmax;
    Rup=0.3*Gmax;
    CO=[CO,-Rdn(:,2:nt)<=pg(:,2:nt)-pg(:,1:nt-1)<=Rup(:,2:nt)]; % рампинг добавляем в ограничения CO
    CO=[CO,Gmin<=pg+rgs<=Gmax]; % Границы генератора с рампингом Rgs
    CO=[CO,(pg(:,2:nt)+rgs(:,2:nt))-(pg(:,1:nt-1)-rgs(:,1:nt-1))<=Rup(:,2:nt)]; % Ограничения по рампинг вверх Rgs
    CO=[CO,-Rdn(:,2:nt)<=(pg(:,2:nt)-rgs(:,2:nt))-(pg(:,1:nt-1)+rgs(:,1:nt-1))]; % Ограничения по рампингу вниз Rgs

    %CO=[CO,sum(pg)-sum(loads)+mg==0];  % баланс мощности 
%     CO=[CO,sum(pg)-sum(loads)+sum(-im+ex)==0];   
    CO=[CO,sum(pg)-sum(loads)+net+wf==0];  % баланс мощности с ветром  
    
%% Ограничения микросети
    % Ограничения распредгенерации 
    Rdn=0.3*Gmmax;
    Rup=0.3*Gmmax;
    CI=[CI,-Rdn(:,2:nt)<=mpg(:,2:nt)-mpg(:,1:nt-1)<=Rup(:,2:nt)]; % рампинг генератора в ограничения CI
%     CI=[CI,Gmmin<=mpg+mrgs<=Gmmax]; % Границы генератора с рампингом Rgs
%     CI=[CI,(mpg(:,2:nt)+mrgs(:,2:nt))-(mpg(:,1:nt-1)-mrgs(:,1:nt-1))<=Rup(:,2:nt)]; % Ограничения по рампинг вверх Rgs
%     CI=[CI,-Rdn(:,2:nt)<=(mpg(:,2:nt)-mrgs(:,2:nt))-(mpg(:,1:nt-1)+mrgs(:,1:nt-1))]; % Ограничения по рампингу вниз Rgs

    % Ограничения для накопитеоей 
    CI=[CI,mb(:,2:nt)==mb(:,1:nt-1)+pb(:,1:nt-1)]; % переход состояния зарядки накопителя в ограничения CI
    CI=[CI,-pb<=mb]; % Батарея не может разряжаться больше чем она имеет

    % Ограничения импорта и экспорта
    %    CI=[CI,ex==max(net,0),im==max(-net,0)];
    CI=[CI,ex-im == net];

    % Ограничения баланса мощности
    CI=[CI,mpg-pb-NCL-Pdb==net];% прогноз 
    
%% Целевая функция передающей сети 
%     OO=sum(pg')*cg1+sum(abs(rgs'))*crgs+(Pdb-pddn)*Pdr'+(pdup-Pdb)*Pdr'+Pex*ex'-Pim*im';
    OO=sum(pg')*cg1+sum(abs(rgs'))*crgs+(Pdb-pddn)*Pdr'+(pdup-Pdb)*Pdr';
    for i=1:nt
        OO=OO+pg(:,i)'*diag(cg2)*pg(:,i);
    end   
    
%% Целевая функция распределительной микросети 
%     OI=sum(mpg')*cg-sum(pd')*cd+sum(mb')*cb-Pex*ex'+Pim*im'-(Pdb-pddn)*(Pdr-pdr)'-(pdup-Pdb)*(Pdr-pdr)';
  
     OI=sum(mpg')*cg(1)+sum(mb')*cb-(Pdb-pddn)*(Pdr-pdr)'-(pdup-Pdb)*(Pdr-pdr)'; 
     CCI = [Pdmin<=pdup<=Pdmax,Pdmin<=pddn<=Pdmax,Gmmin<=mpg<=Gmmax,Pbmin<=pb<=Pbmax,Bmin<=mb<=Bmax,0<=ex,0<=im,Nmin<=net<=Nmax];  ,'bilevel.outersolver','bmibnb'
     
%% Используем решатель solvebilevel - он для выпуклой, но чисто двухуровневой оптимизации
     [K,details] = kkt([CI,Pdmin<=pdup<=Pdmax,Pdmin<=pddn<=Pdmax,Pbmin<=pb<=Pbmax,Bmin<=mb<=Bmax,0<=im],OI,[pg; rgs; Pim; Pex; Pdr; Pinj; theta; pflow])
     ops = sdpsettings('bilevel.algorithm','external','kkt.dualpresolve.lplift',0,'bilevel.outersolver','bmibnb','verbose',1,'debug',1);
     sol = solvebilevel(CO,OO,CI,OI,[mpg; pdup;pddn; mb; ex; im; net;pb],ops);
     
%% Используем решатель bmibnb - он для невыпуклой оптимизации (здесь нужно,чтобы на компьтере была лицензия gurobi)
     %ops = sdpsettings('kkt.dualpresolve.lplift',0);
     %[K,details] = kkt(CI,OI,[pg; rgs; Pim; Pex; Pdr; Pinj; theta; pflow],ops)
     %ops = sdpsettings('solver','bmibnb','bmibnb.lowersolver', 'gurobi','bmibnb.uppersolver', 'scip','bmibnb.lpsolver', 'gurobi','verbose',2,'debug',1);
     
%% Используем решатель scip - он для невыпуклой оптимизации
     %ops = sdpsettings('solver','scip');
     %optimize([K,CO, details.dual<=5],OO,ops)
     
%% Полученные параметры оптимизации (то есть эксплуатационные затраты)
    OO=value(OO) % для передающей системы
    OI=value(OI) % для распределительной микросети


