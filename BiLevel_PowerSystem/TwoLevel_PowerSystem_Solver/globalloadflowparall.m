function [NetData,IsConverged]=globalloadflowparall(NetData)
%%
%Функция совместного решения УР передающей и распределительной сетей
%Распределительные сети считаются параллельно
%%
DData=struct2cell(NetData.DData);
bus=NetData.TrData.bus;
TrSbase=NetData.TrData.basmva;%Base power for transmission side
Settings=NetData.Settings;
NumOfDNets=NetData.NumOfDNets;
VDsol=cell(NumOfDNets,1);%Distr network solutions
IsDdiverged=zeros(NumOfDNets,1);%Distr network solution divergence
TrBusNdx=zeros(NumOfDNets,1);%Indexes of common nodes for transmission and distr networks
PDfeeder=zeros(NumOfDNets,1);%Distr network feeder active power
QDfeeder=zeros(NumOfDNets,1);%Distr network feeder REactive power

DDataNames=fieldnames(NetData.DData);
ItNum=1;%Number of iterations
IsConverged=0;%Convergency flag, zero by default
while(ItNum<NetData.Settings.IterMaxLF)
    if ItNum>1
        bus_sol=NetData.TrData.bus_sol;
    else
        bus_sol=[];
    end
    parfor count=1:NumOfDNets%PF calculations for every distr net        
        CurDistNet=DData{count};%.(DDataNames{count});        
        if ItNum>1
            CurDistNet.Vsw=bus_sol(CurDistNet.TrBusNdx,2);
            CurDistNet.angsw=deg2rad(bus_sol(CurDistNet.TrBusNdx,3));
        else
            CurDistNet.Vsw=1;            
            CurDistNet.angsw=deg2rad(bus(CurDistNet.TrBusNdx,3));%Initial angle
        end
        [Vdist_sol,convflag,Pfeeder,Qfeeder]=distloadflow(CurDistNet,CurDistNet.Vsw,CurDistNet.angsw,Settings);%angsw in radians PF for current Dist net
        if ~convflag
            IsDconverged(count)=1;
            warning(['Solution divergency of the Distribution subsystem: ',DDataNames{count}]);
        else            
            VDsol{count}=Vdist_sol;
            PDfeeder(count)=Pfeeder*CurDistNet.baseMVA/TrSbase;%feeder active power in pu
            QDfeeder(count)=Qfeeder*CurDistNet.baseMVA/TrSbase;%feeder reactive power in pu
            TrBusNdx(count)=CurDistNet.TrBusNdx;
        end
%         %[DPowerFlows,P_loss,Q_loss]=DistLinePQ(Vdist_sol,CurDistNet);
%         if ~convflag
%             IsDistConverged=0;%Divergency of one Distr Subsystem
%             warning(['Solution divergency of the Distribution subsystem: ',DDataNames{count}]);
%             break%Divergency if the divergency of any Distr net solution
%         else
%             NetData.TrData.bus(CurDistNet.TrBusNdx,4)=Pfeeder;%Load is a negative generation
%             NetData.TrData.bus(CurDistNet.TrBusNdx,5)=Qfeeder;%Load is a negative generation                        
%         end        
    end 
    if isempty(find(IsDdiverged,1))
        NetData.TrData.bus(TrBusNdx,4)=PDfeeder;%Load is a negative generation
        NetData.TrData.bus(TrBusNdx,5)=QDfeeder;%Load is a negative generation
        for count=1:NumOfDNets
            NetData.DData.(DDataNames{count}).Vdist_sol=VDsol{count};
        end        
        [~,~,convflag,bus_sol,~,~,~] = trloadflow(NetData.TrData.bus,NetData.TrData.acline, NetData.Settings.epsTrLF, NetData.Settings.IterMaxTrLF,[],NetData.Settings.Display);
        if convflag%Check global convergence
            ItNum=ItNum+1;
            NetData.TrData.bus_sol=bus_sol;
            NetData.Un1=bus_sol(NetData.AllTrBusNdx,2);
            dV=NetData.Un1-NetData.Un0;    
            dVmax=max(abs(dV));
            if dVmax<NetData.Settings.epsLF
                IsConverged=1;
                break
            else                
                NetData.Un0=NetData.Un1;
            end
        else
            warning(['Solution divergency of the Transmission System']);
            break%Divergency if the divergency of any Distr net solution                        
        end                 
    else
        break%%Divergency of Distr Subsystem
    end
    
    
    
%     disp('123');
%     if ~IsDistConverged%%Divergency of Distr Subsystem
%         break
%     else
%         %ПЕРЕПИСАТЬ ЭТОТ РЕШАТЕЛЬ, СДЕЛАТЬ ЧТО-НИБУДЬ С DISPLAY!!!!!
%         [~,~,convflag,bus_sol,~,~,~] = trloadflow(NetData.TrData.bus,NetData.TrData.acline, NetData.Settings.epsTrLF, NetData.Settings.IterMaxTrLF,[],NetData.Settings.Display);
%         if convflag%Check global convergence
%             ItNum=ItNum+1;
%             NetData.TrData.bus_sol=bus_sol;
%             NetData.Un1=bus_sol(NetData.AllTrBusNdx,2);
%             dV=NetData.Un1-NetData.Un0;    
%             dVmax=max(abs(dV));
%             if dVmax<NetData.Settings.epsLF
%                 IsConverged=1;
%                 break
%             else                
%                 NetData.Un0=NetData.Un1;
%             end
%         else
%             warning(['Solution divergency of the Transmission System']);
%             break%Divergency if the divergency of any Distr net solution                        
%         end        
%     end
end
end

