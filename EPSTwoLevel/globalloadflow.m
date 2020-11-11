function [NetData,IsConverged]=globalloadflow(NetData)
%%
%Функция совместного решения УР передающей и распределительной сетей
%%
DDataNames=fieldnames(NetData.DData);
ItNum=1;%Number of iterations
IsConverged=0;%Convergency flag, zero by default
while(ItNum<NetData.Settings.IterMaxLF)
    IsDistConverged=1;%Does all Distr subsystems converged?
    for count=1:NetData.NumOfDNets%PF calculations for every distr net        
        CurDistNet=NetData.DData.(DDataNames{count});        
        if ItNum>1
            CurDistNet.Vsw=NetData.TrData.bus_sol(CurDistNet.TrBusNdx,2);
            CurDistNet.angsw=deg2rad(NetData.TrData.bus_sol(CurDistNet.TrBusNdx,3));
        else
            CurDistNet.Vsw=1;            
            CurDistNet.angsw=deg2rad(NetData.TrData.bus(CurDistNet.TrBusNdx,3));%Initial angle
        end                                     
        [Vdist_sol,convflag,Pfeeder,Qfeeder]=distloadflow(CurDistNet,CurDistNet.Vsw,CurDistNet.angsw,NetData.Settings);%angsw in radians PF for current Dist net        
        NetData.DData.(DDataNames{count}).Vdist_sol=Vdist_sol;
        %[DPowerFlows,P_loss,Q_loss]=DistLinePQ(Vdist_sol,CurDistNet);
        if ~convflag
            IsDistConverged=0;%Divergency of one Distr Subsystem
            warning(['Solution divergency of the Distribution subsystem: ',DDataNames{count}]);
            break%Divergency if the divergency of any Distr net solution
        else
            NetData.TrData.bus(CurDistNet.TrBusNdx,4)=Pfeeder;%Load is a negative generation
            NetData.TrData.bus(CurDistNet.TrBusNdx,5)=Qfeeder;%Load is a negative generation                        
        end        
    end 
    if ~IsDistConverged%%Divergency of Distr Subsystem
        break
    else
        %ПЕРЕПИСАТЬ ЭТОТ РЕШАТЕЛЬ, СДЕЛАТЬ ЧТО-НИБУДЬ С DISPLAY!!!!!
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
    end
end
end

