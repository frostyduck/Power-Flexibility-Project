function [] = ReadAndSaveData()
%Read data from 'TestSystem.xls' and save data to 'TestSystemData.m' 
[~,~,ExtLoad] = xlsread('TestSystem','External load (MW)');%External load
[~,~,Buses] = xlsread('TestSystem','Bus');%Bus data
[~,~,Gen] = xlsread('TestSystem','Generation units');%Generation units data
Bus=struct();%fill Bus struct
BusNum=size(Buses,1);
for count=2:BusNum
    if ~isnan(Buses{count,1})
        Bus.(Buses{count,2}).Unom=Buses{count,3};%Voltage level
        Bus.(Buses{count,2}).LocalLoadShare=Buses{count,4};%Local Load Share
        Bus.(Buses{count,2}).ExtLoad=zeros(1,12);%External Load (MW)
        Bus.(Buses{count,2}).AreaName=Buses{count,5};%Area name 
        ndx=find(strcmp(Gen(:,1),Buses{count,2}));%Find generators in the Bus
        if ~isempty(ndx)
            GenSum=0;
            for c=1:length(ndx)            
                Bus.(Buses{count,2}).Gen.(['Gen',num2str(c)]).Cap=Gen{ndx(c),5};%GenCapacity MW
                GenSum=GenSum+Gen{ndx(c),5};
            end
            Bus.(Buses{count,2}).GenCap=GenSum;%Generation capacity MW  
        else
           Bus.(Buses{count,2}).Gen=[];
           Bus.(Buses{count,2}).GenCap=0;%Generation capacity MW    
        end
    end
end
for count=2:size(ExtLoad)%Add external load data
    if ~isnan(ExtLoad{count,1})
        Bus.(ExtLoad{count,2}).ExtLoad=cell2mat(ExtLoad(count,3:14));
    end
end
% for count=2:size(Gen)%Add generation units capacity data      
%     Bus.(Gen{count,1}).GenCap=Gen{count,5};%Generation capacity
% end
[~,~,Line] = xlsread('TestSystem','Lines');%Считываем данные по линиям
Lines=struct();%LinesData
LineNum=size(Line,1);
for count=2:LineNum
    if ~isnan(Line{count,1})        
        Lines.(Line{count,2}).Fr=Line{count,5};%From Bus
        Lines.(Line{count,2}).To=Line{count,6};%To Bus
        Lines.(Line{count,2}).X=Line{count,7};%Reactance pu
        Lines.(Line{count,2}).R=Line{count,8};%Resistance pu
        Lines.(Line{count,2}).Bc=2*Line{count,9};%Capacitance pu
        if ~strcmp(Line{count,12},'-')
            %Нашел косяки в коэффициентах трансформации, принял просто
            %отношения номинальных напряжений, можно корректировать при
            %желании
            %Lines.(Line{count,2}).Kt=(Line{count,12}/Bus.(Line{count,5}).Unom)/(Line{count,13}/Bus.(Line{count,6}).Unom);
            %Lines.(Line{count,2}).Kt=(Line{count,13}/Bus.(Line{count,6}).Unom)/(Line{count,12}/Bus.(Line{count,5}).Unom);
            %if Lines.(Line{count,2}).Kt>1.2 || Lines.(Line{count,2}).Kt<0.8
            %    disp('123');
            %end
            Lines.(Line{count,2}).Kt=1;%Transf koef
        else
            Lines.(Line{count,2}).Kt=1;%Transf koef
        end
    end
end
Areas=struct();%Areas data
Areas.Area1.PeakLocalLoad=33950;%MW
Areas.Area1.PeakExtLoad=11000;%MW

Areas.Area2.PeakLocalLoad=20960;%MW
Areas.Area2.PeakExtLoad=4000;%MW

Areas.Area3.PeakLocalLoad=10390;%MW
Areas.Area3.PeakExtLoad=430;%MW

Areas.Area4.PeakLocalLoad=15100;%MW
Areas.Area4.PeakExtLoad=16500;%MW

Areas.Area5.PeakLocalLoad=42070;%MW
Areas.Area5.PeakExtLoad=13500;%MW

save('TestSystemData','Bus','Lines','Areas');
end

