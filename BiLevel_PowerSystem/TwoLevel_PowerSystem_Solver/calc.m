function [delP,delQ,P,Q,conv_flag] = ...
                 calc(V,ang,Y,Pg,Qg,Pl,Ql,sw_bno,g_bno,tol)
% Syntax:  [delP,delQ,P,Q,conv_flag] = 
%                calc(V,ang,Y,Pg,Qg,Pl,Ql,sw_bno,g_bno,tol)
%
% Purpose: calculates power mismatch and checks convergence
%          also determines the values of P and Q based on the 
%          supplied values of voltage magnitude and angle
% ************************************************************
% swing_bus = 1;
% gen_bus = 2;
% load_bus = 3;
% voltage in rectangular coordinate
V_rect = V.*exp(1i*ang);  
% bus current injection
cur_inj = Y*V_rect;
% power output based on voltages 
S = V_rect.*conj(cur_inj);
P = real(S); Q = imag(S);
delP = Pg - Pl - P;
delQ = Qg - Ql - Q;
% zero out mismatches on swing bus and generation bus
delP=delP.*sw_bno;%Небалансы по активнйо мощности учитываются за всех узлов, кроме балансирующего
delQ=delQ.*sw_bno;
delQ=delQ.*g_bno;%Небалансы по реактивной мощности учитываются только для PQ узлов
%  total mismatch
[pmis]=max(abs(delP));
[qmis]=max(abs(delQ));
mism = pmis+qmis;
if mism > tol,
    conv_flag = 1;
  else
    conv_flag = 0;
end
return

