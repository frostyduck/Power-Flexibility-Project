function NetData = InitData(NetData)
%INITDATA Summary of this function goes here
%   Detailed explanation goes here
NetData.TrData=feval(NetData.TDataFileName);%Transmission Network Data
NumOfDNets=size(NetData.DnameTnode,1);%Number of dist nets
NetData.NumOfDNets=NumOfDNets;
NetData.AllTrBusNdx=zeros(NumOfDNets,1);%Number of all ndxes in Transm network
NetData.Un0=10*ones(NumOfDNets,1);
NetData.Un1=100*ones(NumOfDNets,1);
for count=1:NumOfDNets
    CurDistNet=feval(NetData.DnameTnode{count,1});%Current distribution system  
    CurDistNet.TrBusNdx=find(NetData.TrData.bus(:,1)==NetData.DnameTnode{count,2},1);%индекс узла в передающей сети 
    NetData.TrData.bus(CurDistNet.TrBusNdx,4:7)=0;%Delete any load or generation from distribution system node!!!!!
    NetData.AllTrBusNdx(count)=CurDistNet.TrBusNdx;
    %CurDistNet.Vsw=1;%Initial voltage
    %CurDistNet.angsw=deg2rad(NetData.TrData.bus(CurDistNet.TrBusNdx,3));%Initial angle    
    swndx=find(CurDistNet.bus(:,10)==1,1);%swing bus index Индекс балансирующего узла
    CurDistNet.swnum=CurDistNet.bus(swndx,1);%swing bus number Номер балансирующего узла        
    CurDistNet=reordernodes(CurDistNet,CurDistNet.swnum);%Процедура обработки графа и вычисления матриц    
    NetData.DData.(NetData.DnameTnode{count,1})=CurDistNet;        
end
end

