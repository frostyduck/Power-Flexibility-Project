function [Vdist_sol,convflag,Pfeeder,Qfeeder]=distloadflow(DistNet,Vsw,angsw,Settings)
%%
%Функция решения УР системы BF методом
%L, G, Zv - матрицы, участвующие в расчете BF-методом.
%bus - исходная информация по узлам схемы.
%IntToNode - матрица соответствия нумерации новый узел (последовательная нумерация) -> исходный узел (bus);
%ndxBusToInt - соответствие индексов узлов в матрице bus индексам узлов при последовательной нумерации
%eps - минимальный небаланс для проверки сходимости
%Vdist_sol - решение для распредсети
%convflag - 1 converged, 0 - diverged
%Pfeeder,Qfeeder - active and reactive power of the feeder (Transmission network connection)
%Vsw - модуль напряжения балансирующего узла, о.е.
%angsw - угол напряжения балансирующего узла, градусы IN RADIANS!!!
%%
Pfeeder=0;Qfeeder=0;%by default values
eps=Settings.epsDistLF;
convflag=0;%По умолчанию режим не сошелся
Vsw=Vsw*exp(1i*(angsw));%Напряжение балансирующего узла
%L=DistNet.L;
G=DistNet.G;
Zv=DistNet.Zv;
bus=DistNet.bus;
%IntToNode=DistNet.IntToNode;
ndxBusToInt=DistNet.ndxBusToInt;

%Vsw=bus(swndx,2)*exp(1i*deg2rad( bus(swndx,3) ));
V0=Vsw*ones(size(bus,1)-1,1);

V=bus(:,2).*exp(1i*deg2rad( bus(:,3) ));%Исходное приближение напряжений в узлах
Gt=G';
for N=1:Settings.IterMaxDistLF%Итерации BF-метода    
    I=conj(  (bus(:,6)+1i*bus(:,7))./V );%Инъекции тока от нагрузки в узлах
    I=I+conj(  (bus(:,4)+1i*bus(:,5))./V );%Инъекции тока от генерации в узлах
    I=I+conj(  (bus(:,8)+1i*bus(:,9)).*V );%Инъекции тока от шунтов в узлах
    Is=I(ndxBusToInt);
    Is(1)=[];%Удаляем инъекцию тока балансирующего узла
    Iv=Gt*Is;
    V1=[Vsw;V0-G*Zv*Iv];
    dV=V-V1;
    dVmax=max(abs(dV));
    if Settings.Display
        disp(['Номер итерации -',num2str(N),'; dVmax=',num2str(dVmax)])
    end
    if dVmax<eps
        V=V1;
        convflag=1;
        Sfeeder=V(1)*conj(Iv(1));
        Pfeeder=real(Sfeeder);
        Qfeeder=imag(Sfeeder);        
        break
    else
        %max(abs(dV))
        V=V1;
    end
end
%Vdist_sol=[IntToNode abs(V) rad2deg(angle(V))];
Vdist_sol=V;
end

