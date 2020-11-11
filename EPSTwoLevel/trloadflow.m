function [V,ang,IsConverged,bus_sol,line,BusStruct,Y] = trloadflow(bus,line,tol,iter_max,acc, display)
%� ������� ���������� ������� ������� ��������� �������� �������������, �����
%���������� ����������!!! ������:
if strcmp(display,'on')
    disp('����� ��� ����������� �������� ���������� � PU �����, �� ����� � ����� ����������!!!');
    disp('����� ��� ����������� ��������� ������������� ������������ PU->PQ->PU!!!');
    disp('����� ��� ����������� ������������� ����.������. ��� ��������� ����������!!!');
end
BusStruct=struct();%���������, � ������� ����������� ������ �����, ������� �������� � PQ � PU, � ����� PQ �����
BusStruct.CheckTapChangers=0;%����� �������� �� ������������! ����������� ������������ ������������� ��� ��������� ���������� 1-��, 0-���
BusStruct.MultiplePUPQ=1;%������������� ������������ PU->PQ->PU 1-��, 0-���.
BusStruct.ConsiderUrefNewVersion=1;%1-���������� � �������� ����������� �������� 16 ������� ������� 0-���������� � �������� ����������� �������� �������� ������� ����������
if BusStruct.MultiplePUPQ%������������� ������������ PU->PQ->PU
    if BusStruct.ConsiderUrefNewVersion
        BusStruct.Uref0=bus(:,16);%���������� � �������� ����������� �������� 16 ������� �������        
        BusStruct.NumOfPUPQSwitches=zeros(size(bus,1),1);%����� ������������ PU->PQ, ������ �� ��������� ������������
    else        
        BusStruct.Uref0=bus(:,2);%���������� �������� �������� ���������� �� ���� �����, ��� PU ����� ��� ����� Uref
    end
end

[nbus] = size(bus,1);     % number of buses

CheckFlag=1;%����������� ����, ������� ������, ����� �������������� ������ ��������, ������� ��������� �������� ���������� � � ������� ����� �������� ���������������� �������.
IsConverged=0;%�� ��������� �������, ��� �� ���������!

ExtIter=0;%Outer iteration PU->PQ tripping
% ��������� ���������� ������� �� ���������
if isempty(tol);tol = 1e-11;end
if isempty(iter_max);iter_max = 30;end

Y=y_sparseNew(bus,line);

bus_type = bus(:,10);
%������� ������ ���� ����� ���������� �������� bus(:,2), �����, ���
%������������� (if BusStruct.ConsiderUrefNewVersion), ���������� ��� PU
%����� ��������� ��������
V = ones(nbus,1);%������� ������ ���� ����� ���������� �������� bus(:,2);%�� ������ ����� ������!!!
%������� ��������� �� �� ��������, �� ������������
if BusStruct.CheckTapChangers%bus_int ����� ��� ������������� ����. ������������� ��� ��������� ����������
    busmax = max(bus(:,1));
    bus_int = zeros(busmax,1);
    ibus = (1:nbus)';
    bus_int(bus(:,1)) = ibus;    
    Ndx1=bus_int(line(:,1));
    Ndx2=bus_int(line(:,2));
end
%��� ������������� (if BusStruct.ConsiderUrefNewVersion), ���������� ��� PU ����� ��������� ��������
BusStruct.PU_no=find(bus_type==2);%Positions of PU nodes
if BusStruct.ConsiderUrefNewVersion%����� ������ � Uref ��� PU ����� � 16 ������� 
    V(BusStruct.PU_no) = bus(BusStruct.PU_no,16);%���������� � PU ����� Uref
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
BusStruct.PQV_no=find(bus_type >=2);%Positions of all nodes except swing bus %����, ������� ������������� ���� ����� ����� ��
BusStruct.PQ_no=find(bus_type==3);%Positions of PQ nodes %����, ������� ������������� PQ �����
BusStruct.InitialPU=BusStruct.PU_no;%� �������� ��������� ������ ���� ���� PU, ������ ��� ����� ������������� PU->PQ
BusStruct.SW_no=swing_index;%������ ��
BusStruct.SWPU_no=[swing_index; BusStruct.PU_no];%������� SW � PU �����
g_bno(BusStruct.PU_no)=bus_zeros(BusStruct.PU_no);

st = clock;     % start the iteration time clock
if strcmp(display,'on')
    options = optimoptions('fsolve','Jacobian','on','TolFun',tol,'TolX',tol,'Display','iter-detailed','MaxIterations',iter_max);%�������� �����
else
    options = optimoptions('fsolve','Jacobian','on','TolFun',tol,'TolX',tol,'Display','off','MaxIterations',iter_max);%��������� ����� �����������
end
%% start iteration process for main Newton_Raphson solution
%������� �� ����������� � ���� �������� �����, �� ������ �����
%���������� ������� ��������� ��� ������������ PU->PQ �����, �����
%������� ������� ����������� ������������� ������������, � ������
%������������� ��������������� ���� ��� �����, ���� ���� ���� � ������������ ���������� �� Q (������� ������ ������ �����)
while 1%������� ����, �������� � ������������� PU->PQ
    if (ExtIter>=1)
                
        [lim_flag,Qg,bus_type,BusStruct] = chq_lim(qg_max,qg_min,bus_type,Qg,display,BusStruct,bus,V);
        
        if lim_flag %������������� �����-�� ����������� ����� � chq_lim
            BusStruct.PQV_no=find(bus_type >=2);
            BusStruct.PQ_no=find(bus_type==3);
            BusStruct.PU_no=find(bus_type==2);%Positions of PU nodes
            BusStruct.SWPU_no=[swing_index; BusStruct.PU_no];%������� SW � PU �����
        else
            if strcmp(display,'on')
                disp(['Load flow converged after ',num2str(ExtIter),' external iterations']);
            end
            break%No PU->PQ switching, power flow convrged! ������������ ����� �� ����, ����� �������            
        end
    end
    %� �������� ������� ����� ��������� �������� ������������ � PQ �� PU,
    %������� ���������� ����� ��������� �������� ���������� � PU ����
    if BusStruct.ConsiderUrefNewVersion%����� ������ � Uref ��� PU ����� � 16 �������
        V(BusStruct.PU_no) = bus(BusStruct.PU_no,16);%���������� � PU ����� Uref
    end
    %������������� �� �� ��������! �� ���������!
    if BusStruct.CheckTapChangers%���� ����� ������������� ����. ������������� ����� ��������, �� ��������� ������ �������, ���� � ��� �����������
        V(Ndx1)=V(Ndx2).*line(:,6);
    end    
    V0=V(BusStruct.PQ_no);%���������� ��� PQ ����� ����� �� ��������, ����������� PQ ����� ����� �������� ���������� ���������������� PU-PQ-PU    
    ang0=ang(BusStruct.PQV_no);
    x0=[ang0;V0];
    
    f=@(x)myfunwithJac(x,V,ang,Y,Pg,Qg,Pl,Ql,BusStruct);
    [x,fval,~] = fsolve(f,x0,options);
    ang(BusStruct.PQV_no)=x(1:length(BusStruct.PQV_no));%��������!
    V(BusStruct.PQ_no)=x(length(BusStruct.PQV_no)+1:length(x));%��������!
    %��������� ���������� ��
    ndx=find(abs(fval)>options.TolFun,1);%If at least one element bigger than TolFun then loadflow diverged
    if ~isempty(ndx)%load flow diverged
        IsConverged=0;
        bus_sol=bus;
        line;
        return
    else
        IsConverged=1;
        PFTolFun=options.TolFun;%��������� ��������������� ���������, ���������� ��� ���������� �������� ������� �� ��������� � �����
    end    
    Q = calcQ(V,ang,Y);
    %��������� ��������� � ���� ����� ��������� �� ������� ���� Q ����
    %���������, ������������ ��������� Ql
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
Pg=P+Pl;%�������� ��������� �� ���� �����, ������� �������������
notGenNdx=find(abs(Pg)<PFTolFun);%���� ��� ��������, ������� ������ ����������������
Pg(notGenNdx)=0;
Qg=Q+Ql;%���������� ��������� �� ���� �����, ������� �������������
notGenNdx=find(abs(Qg)<PFTolFun);%���� ��� ��������, ������� ������ ����������������
Qg(notGenNdx)=0;
bus_sol=bus;
bus_sol(:,2)=V;
bus_sol(:,3)=ang*(180/pi);
bus_sol(:,4)=Pg;%���� ����������� ��������� �� ���� �����
bus_sol(:,5)=Qg;%���� ����������� ��������� �� ���� �����
return
end
function F=myfun(x,V,ang,Y,Pg,Qg,Pl,Ql,BusStruct)%Finite Differences Jacobian ��������� �������
ang(BusStruct.PQV_no)=x(1:length(BusStruct.PQV_no));%��������!
V(BusStruct.PQ_no)=x(length(BusStruct.PQV_no)+1:length(x));%��������!
V_rect = V.*exp(1i*ang);
cur_inj = Y*V_rect;
S = V_rect.*conj(cur_inj);
P = real(S); Q = imag(S);
delP = Pg - Pl - P;
delQ = Qg - Ql - Q;
F=[delP(BusStruct.PQV_no);delQ(BusStruct.PQ_no)];
end
function [F,Jac]=myfunwithJac(x,V,ang,Y,Pg,Qg,Pl,Ql,BusStruct)%Analytic Jacobian ������������� �������
ang(BusStruct.PQV_no)=x(1:length(BusStruct.PQV_no));%��������!
V(BusStruct.PQ_no)=x(length(BusStruct.PQV_no)+1:length(x));%��������!
Jac=-(form_jacMine(V,ang,Y,BusStruct));%Analytic Jacobian ������������� �������
V_rect = V.*exp(1i*ang);
cur_inj = Y*V_rect;
S = V_rect.*conj(cur_inj);
P = real(S); Q = imag(S);
delP = Pg - Pl - P;
delQ = Qg - Ql - Q;
F=[delP(BusStruct.PQV_no);delQ(BusStruct.PQ_no)];
end

