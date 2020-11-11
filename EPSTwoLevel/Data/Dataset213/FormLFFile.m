function [] = FormLFFile()
%Form Load Flow Data File 
load('TestSystemData');%Load Bus and Lines data
%Make k times from peak load k<1
Sbase=100;
QtoP=0.3;%Qload/Pload
mnth=1;%Number of month for external load
k=0.7;%Load koefficient коэффициент загрузки от пиковой
kgen=1.05;%Generation+losses koefficient коэффициент загрузки по генерации+потери
A1SumLoad=Areas.Area1.PeakLocalLoad*k/Sbase;
A2SumLoad=Areas.Area2.PeakLocalLoad*k/Sbase;
A3SumLoad=Areas.Area3.PeakLocalLoad*k/Sbase;
A4SumLoad=Areas.Area4.PeakLocalLoad*k/Sbase;
A5SumLoad=Areas.Area5.PeakLocalLoad*k/Sbase;
BusNames=fieldnames(Bus);
BusNumName=struct();
NumOfBuses=length(BusNames);
LineNames=fieldnames(Lines);
NumOfLines=length(LineNames);
bus=zeros(NumOfBuses,16);%bus structure
acline=zeros(NumOfLines,12);%lines structure
PLA1=0;PLA2=0;PLA3=0;PLA4=0;PLA5=0;%sum load of areas
PGA1=0;PGA2=0;PGA3=0;PGA4=0;PGA5=0;%sum generation of areas
for count=1:NumOfBuses 
    CurrBus=Bus.(BusNames{count});%Current considered bus
    BusNumName.(['B',num2str(count)])=BusNames{count};
    Bus.(BusNames{count}).BusNum=count;%Set number for the current bus
    if strcmp(CurrBus.AreaName,'Area1')
        PLA1=PLA1+CurrBus.LocalLoadShare*A1SumLoad+...
            CurrBus.ExtLoad(mnth)/Sbase;%External load of the mnth-th month taken        
        PGA1=PGA1+CurrBus.GenCap/Sbase;
    elseif strcmp(CurrBus.AreaName,'Area2')
        PLA2=PLA2+CurrBus.LocalLoadShare*A2SumLoad+...
            CurrBus.ExtLoad(mnth)/Sbase;%External load of the mnth-th month taken
        PGA2=PGA2+CurrBus.GenCap/Sbase;
    elseif strcmp(CurrBus.AreaName,'Area3')
        PLA3=PLA3+CurrBus.LocalLoadShare*A3SumLoad+...
            CurrBus.ExtLoad(mnth)/Sbase;%External load of the mnth-th month taken
        PGA3=PGA3+CurrBus.GenCap/Sbase;        
    elseif strcmp(CurrBus.AreaName,'Area4')
        PLA4=PLA4+CurrBus.LocalLoadShare*A4SumLoad+...
            CurrBus.ExtLoad(mnth)/Sbase;%External load of the mnth-th month taken
        PGA4=PGA4+CurrBus.GenCap/Sbase;        
    elseif strcmp(CurrBus.AreaName,'Area5')
        PLA5=PLA5+CurrBus.LocalLoadShare*A5SumLoad+...
            CurrBus.ExtLoad(mnth)/Sbase;%External load of the mnth-th month taken
        PGA5=PGA5+CurrBus.GenCap/Sbase;        
    end        
end
for count=1:NumOfBuses% fill in bus matrix
    CurrBus=Bus.(BusNames{count});
    bus(count,1)=count;%Bus number
    bus(count,2)=1;%voltage magnitude pu 
    bus(count,13)=CurrBus.Unom;%nominal voltage     
    if strcmp(CurrBus.AreaName,'Area1')
        bus(count,4)=kgen*CurrBus.GenCap*PLA1/PGA1/Sbase;
        bus(count,6)=CurrBus.LocalLoadShare*A1SumLoad+...
            CurrBus.ExtLoad(mnth)/Sbase;%Pload
        bus(count,7)=QtoP*CurrBus.LocalLoadShare*PLA1;%Qload
    elseif strcmp(CurrBus.AreaName,'Area2')
        bus(count,4)=kgen*CurrBus.GenCap*PLA2/PGA2/Sbase;
        bus(count,6)=CurrBus.LocalLoadShare*A2SumLoad+...
            CurrBus.ExtLoad(mnth)/Sbase;%Pload
        bus(count,7)=QtoP*CurrBus.LocalLoadShare*PLA2;%Qload
    elseif strcmp(CurrBus.AreaName,'Area3')
        bus(count,4)=kgen*CurrBus.GenCap*PLA3/PGA3/Sbase;
        bus(count,6)=CurrBus.LocalLoadShare*A3SumLoad+...
            CurrBus.ExtLoad(mnth)/Sbase;%Pload
        bus(count,7)=QtoP*CurrBus.LocalLoadShare*PLA3;%Qload
    elseif strcmp(CurrBus.AreaName,'Area4')
        bus(count,4)=kgen*CurrBus.GenCap*PLA4/PGA4/Sbase;
        bus(count,6)=CurrBus.LocalLoadShare*A4SumLoad+...
            CurrBus.ExtLoad(mnth)/Sbase;%Pload
        bus(count,7)=QtoP*CurrBus.LocalLoadShare*PLA4;%Qload
    elseif strcmp(CurrBus.AreaName,'Area5')
        bus(count,4)=kgen*CurrBus.GenCap*PLA5/PGA5/Sbase;
        bus(count,6)=CurrBus.LocalLoadShare*A5SumLoad+...
            CurrBus.ExtLoad(mnth)/Sbase;%Pload
        bus(count,7)=QtoP*CurrBus.LocalLoadShare*PLA5;%Qload
    end
    if CurrBus.GenCap~=0%Add PV bus with reactive power
        bus(count,10)=2;%Add PV bus
        bus(count,11)=0.8*bus(count,4);%Add Qmax
        bus(count,12)=-0.8*bus(count,4);%Add Qmin 
        bus(count,16)=bus(count,2);%Reference voltage for PU bus         
    else%Add PQ bus
        bus(count,10)=3;%Add PQ bus
    end
end

bus(1,10)=1;%Swing bus


for count=1:NumOfLines
    CurLine=Lines.(LineNames{count});%Current considered bus
    FromBus=Bus.(CurLine.Fr);
    ToBus=Bus.(CurLine.To);
    acline(count,1)=FromBus.BusNum;
    acline(count,2)=ToBus.BusNum;
    acline(count,3)=CurLine.R;
    acline(count,4)=CurLine.X;
    acline(count,5)=CurLine.Bc;
    acline(count,6)=CurLine.Kt;
    acline(count,11)=1;%Num of parallel
    acline(count,12)=1;%State 1-on, 0-off
end

%Add distribution network to the model
MaxBusNum=bus(size(bus,1),1)+1;
%Add 33 node transmission system to C18 bus
acline=[acline; (Bus.C18.BusNum)  MaxBusNum  0.00082645/10    0.018512/10           0      1.0463           0           0           0           0           1           1];
bus=[bus; MaxBusNum 1 0 0 0 0 0 0 0 3 0 0 12 0 0 0];
MaxBusNum=MaxBusNum+1;
%Add 33 node transmission system to C23 bus
acline=[acline; (Bus.C23.BusNum)  MaxBusNum  0.00082645/10    0.018512/10           0      1.0463           0           0           0           0           1           1];
bus=[bus; MaxBusNum 1 0 0 0 0 0 0 0 3 0 0 12 0 0 0];
MaxBusNum=MaxBusNum+1;
%Add 33 node transmission system to C24 bus
acline=[acline; (Bus.C24.BusNum)  MaxBusNum  0.00082645/10    0.018512/10           0      1.0463           0           0           0           0           1           1];
bus=[bus; MaxBusNum 1 0 0 0 0 0 0 0 3 0 0 12 0 0 0];
MaxBusNum=MaxBusNum+1;
%Add 33 node transmission system to C27 bus
acline=[acline; (Bus.C27.BusNum)  MaxBusNum  0.00082645/10    0.018512/10           0      1.0463           0           0           0           0           1           1];
bus=[bus; MaxBusNum 1 0 0 0 0 0 0 0 3 0 0 12 0 0 0];
MaxBusNum=MaxBusNum+1;
%Add 33 node transmission system to C26 bus
acline=[acline; (Bus.C26.BusNum)  MaxBusNum  0.00082645/10    0.018512/10           0      1.0463           0           0           0           0           1           1];
bus=[bus; MaxBusNum 1 0 0 0 0 0 0 0 3 0 0 12 0 0 0];
MaxBusNum=MaxBusNum+1;
%Add 33 node transmission system to C31 bus
acline=[acline; (Bus.C31.BusNum)  MaxBusNum  0.00082645/10    0.018512/10           0      1.0463           0           0           0           0           1           1];
bus=[bus; MaxBusNum 1 0 0 0 0 0 0 0 3 0 0 12 0 0 0];
MaxBusNum=MaxBusNum+1;


disp('Write data to file');
% if ~isempty(folder)
%     data_file=[folder,'\',FileNameWithoutCDU,'.m'];
% else
%     data_file=[FileNameWithoutCDU,'.m'];
% end
[fid, message]=fopen('Test213RAW.m','wt');
if ~isempty(message)
    disp('ERROR:')
    disp(message);
    return
end
%Zapis informacii Bus.con
fprintf(fid,['function Test213RAW=Test213RAW()\n']);
name='bus';
if ~isempty(bus)
    fprintf(fid,'%s = [ ... \n',name);
    for counter=1:size(bus,1)
        fprintf(fid,[num2str(bus(counter,:)),';\n']);
    end
    fprintf(fid,' ];\n\n');
end
fprintf(fid,'Test213RAW.bus=bus;\n\n');
%Zapis informacii Line.con
name='acline';
if ~isempty(acline)
    fprintf(fid,'%s = [ ... \n',name);
    for counter=1:size(acline,1)
        fprintf(fid,[num2str(acline(counter,:)),';\n']);
    end
    fprintf(fid,' ];\n\n');
end
fprintf(fid,'Test213RAW.acline=acline;\n\n');
name=num2str(100);
fprintf(fid,'Test213RAW.basmva = %s ;\n\n',name);
fprintf(fid,'end');
fclose(fid);%Закрываем файл с информацией по УР
   
% tic
%[V,ang,IsConverged,bus_sol,line,BusStruct,Y] = loadflow(bus,acline,1e-6,100,1,'on');
% toc
% disp('123');
end

