function [DistPowerFlows,P_loss,Q_loss]= DistLinePQ(V,DistNet)

%VV = V.*exp(1i*ang);  %Напряжение в точке решения
%tps = acline(:,6).*exp(1i*acline(:,7)*pi/180);
%tpj=conj(tps);
acline=DistNet.acline;
r = acline(:,3); 
x = acline(:,4);
%chrg = acline(:,12).*acline(:,5)/2;

z = r + 1i*x;
y = acline(:,12)./z;
Vf = V(DistNet.NodeToInt(acline(:,1)));
Vt = V(DistNet.NodeToInt(acline(:,2)));

%Is=Vf.*(y + i*chrg)./(tps.*tpj) - Vt.*y./tpj;
%Ir=Vt.*(y + i*chrg) - Vf.*y./tps;
dV=Vf - Vt;
Is=dV.*y;
Ir=-dV.*y;

Ssend  = Vf.*conj(Is);
Srec  = Vt.*conj(Ir);

P_loss = sum(real(Ssend)) + sum(real(Srec));
Q_loss = sum(imag(Ssend)) + sum(imag(Srec));

DistPowerFlows=zeros(size(acline,1),6);
DistPowerFlows(:,1)=acline(:,1);%Node1
DistPowerFlows(:,2)=acline(:,2);%Node2
%DistPowerFlows(:,3)=acline(:,11);%Parall
DistPowerFlows(:,3)=real(Ssend);%Pfrom pu
DistPowerFlows(:,4)=imag(Ssend);%Qfrom pu
DistPowerFlows(:,5)=real(Srec);%Pto pu
DistPowerFlows(:,6)=imag(Srec);%Qto pu
return
end