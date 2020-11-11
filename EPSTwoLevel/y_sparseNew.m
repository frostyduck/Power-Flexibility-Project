function     [Y] = y_sparseNew(bus,acline)

nbus = length(bus(:,1));     % number of buses
ibus = (1:nbus)';

busmax = max(bus(:,1));
bus_int = zeros(busmax,1);
bus_int(bus(:,1)) = ibus;

ndxfrom=bus_int(acline(:,1));    
ndxto=bus_int(acline(:,2)); 
chrg = acline(:,12).*acline(:,5)/2;
r = acline(:,3);
x = acline(:,4);
z = r + i*x;
y = acline(:,12)./z;
ts = acline(:,6).*exp(1i*acline(:,7)*pi/180);
ts2= ts.*conj(ts);
Y = sparse(ndxfrom,ndxto,-y./conj(ts),nbus,nbus) + ...
      sparse(ndxto,ndxfrom,-y./ts,nbus,nbus) + ...
      sparse(ndxfrom,ndxfrom,(y+j*chrg)./ts2,nbus,nbus)+ ...
      sparse(ndxto,ndxto,y+j*chrg,nbus,nbus); 

b = find(diag(Y) == 0);
if ~isempty(b)
  Y = Y - sparse(b,b,j*1e-6,nbus,nbus);
end
Gb = bus(:,8);     % bus conductance
Bb = bus(:,9);     % bus susceptance  
Y = Y + sparse(ibus,ibus,Gb+1i*Bb,nbus,nbus);  

return

