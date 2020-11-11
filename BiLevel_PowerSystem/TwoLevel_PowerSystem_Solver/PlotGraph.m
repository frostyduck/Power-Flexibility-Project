function PlotGraph()
%Визуализация графа сети
[baseMVA,acline,bus]=case33;%Подгружаемая схема
nbus=size(bus,1);%num of buses
busmax = max(bus(:,1));
bus_int = zeros(busmax,1);
ibus = (1:nbus)';
bus_int(bus(:,1)) = ibus;

nbranch=size(acline,1);%num of branches
br=zeros(nbranch,2);
br(:,1)=bus_int(acline(:,1));
br(:,2)=bus_int(acline(:,2));
weight=abs(acline(:,3)+1i*acline(:,4));
EdgeTable = table(br,weight,'VariableNames',{'EndNodes','Weight'})
G = graph(EdgeTable)
BusNames=cell(nbus,1);
for count=1:nbus
    BusNames{count,1}=['N',num2str(bus(count,1))];
    %NodeProps = table({['Bus',num2str(bus(count,1))]}', [1]', ...
    %    'VariableNames', {'NodeName' 'Voltage'});
    %G=addnode(G,NodeProps);
end
G.Nodes.Names=BusNames;
plot(G,'NodeLabel',G.Nodes.Names);

end