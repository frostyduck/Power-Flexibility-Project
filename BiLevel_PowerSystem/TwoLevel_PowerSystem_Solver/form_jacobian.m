function [Jac,J11,J12,J21,J22]=form_jacobian(V,ang,Y,BusStruct)
exp_ang = exp(1i*ang);
V_rect = V.*exp_ang;
CV_rect=conj(V_rect);
Y_con = conj(Y);
i_c=Y_con*CV_rect;
S=V_rect.*i_c;%S=U*conj(I)
S=sparse(diag(S));
Vdia=sparse(diag(V_rect));
CVdia=conj(Vdia);
Vmag=sparse(diag(abs(V)));
S1=Vdia*Y_con*CVdia;%S=U*conj(I)=U*conj(Y)*conj(U) written as diag
%Мой вариант записи исключения--------------------------------
tt2=S-S1;
tt1=(S+S1)/Vmag;
%dP/ddelta
J11=-imag(tt2);
J11(BusStruct.SW_no,:)=[];
J11(:,BusStruct.SW_no)=[];
%dP/dU
J12=real(tt1);
J12(:,BusStruct.SWPU_no)=[];
J12(BusStruct.SW_no,:)=[];
%dQ/ddelta
J21=real(tt2);
J21(BusStruct.SWPU_no,:)=[];
J21(:,BusStruct.SW_no)=[];
%dQ/dU
J22=imag(tt1);
J22(BusStruct.SWPU_no,:)=[];
J22(:,BusStruct.SWPU_no)=[];

Jac = [J11 J12;
    J21 J22];
end

