%������������ ����� ����
[baseMVA,acline,bus]=case33();%������������ �����
swndx=find(bus(:,10)==1,1);%������ �������������� ����
swnum=bus(swndx,1);%swing bus number ����� �������������� ����
s=acline(:,1);
t=acline(:,2);
G=graph(s,t);
subplot(2,1,1)
plot(G)
hold on
[nodetoint,L,G,Zv]=reordernodes(acline,swnum,bus);
ss=nodetoint(s);
tt=nodetoint(t);
G1=graph(ss,tt);
subplot(2,1,2)
plot(G1)