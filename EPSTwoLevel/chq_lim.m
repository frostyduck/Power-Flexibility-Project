function [lim_flag,Qg,bus_type,BusStruct] = chq_lim(qg_max,qg_min,bus_type,Qg,display,BusStruct,bus,V)
lim_flag = 0;% indicates whether limit has been reached ����, ������� ��������� �� ��������� �����������
%gen_idx = find(bus_type ==2);
%qg_max_idx = find(Qg(gen_idx)>qg_max(gen_idx));%Qg>Qgmax

%Qg>Qgmax and node is PU node
%qg_max_idx = find( (Qg>qg_max) & bus_type ==2);
qg_max_idx = BusStruct.PU_no( find( Qg(BusStruct.PU_no)>qg_max(BusStruct.PU_no) ) );
%Qg<Qgmin and node is PU node
qg_min_idx = BusStruct.PU_no( find( Qg(BusStruct.PU_no)<qg_min(BusStruct.PU_no) ) );
% %Choose the node with the biggest mismatch �������� ���� � ������������ ����������
% [dQ,idx]=max(Qg(qg_max_idx)-qg_max(qg_max_idx));
% qg_max_idx = qg_max_idx (idx);
% qg_min_idx = find(Qg(gen_idx)<qg_min(gen_idx));%Qg<Qgmin
if BusStruct.MultiplePUPQ%������������� ������������ PU->PQ->PU 1-��, 0-���.
    PQtoPUndx=BusStruct.InitialPU( find( bus_type(BusStruct.InitialPU)-bus(BusStruct.InitialPU,10)==1) );%������� �����, ������� ���������� ��������� PU->PQ
    if ~isempty(PQtoPUndx)
        ndxQmax=PQtoPUndx( find( abs(Qg(PQtoPUndx)-bus(PQtoPUndx,11))<1e-5 ));%������� �����, ������������� PU->PQmax
        ndxQmin=setdiff(PQtoPUndx,ndxQmax);%������� �����, ������������� PU->PQmin
        if ~isempty(ndxQmax)
            ndxQmaxSwitch=ndxQmax(find(V(ndxQmax)>BusStruct.Uref0(ndxQmax)));
            if ~isempty(ndxQmaxSwitch)
                ndxToSwitch=ndxQmaxSwitch( find(BusStruct.NumOfPUPQSwitches(ndxQmaxSwitch)<=1) );%����������� ������ �� ����, ������� ������������� �� ����� ������ 2 ���
                if ~isempty(ndxToSwitch)                    
                    if strcmp(display,'on')
                        disp(['Generator at Bus ',mat2str(bus(ndxQmaxSwitch)),' PQ->PU switching']);
                    end
                    bus_type(ndxQmaxSwitch) = 2*ones(length(ndxQmaxSwitch),1);
                    BusStruct.NumOfPUPQSwitches(ndxToSwitch)=BusStruct.NumOfPUPQSwitches(ndxToSwitch)+ones(length(ndxToSwitch),1);%����������� ������� ������������ �� �������
                end
            end
        end
        if ~isempty(ndxQmin)
            ndxQminSwitch=ndxQmin(find(V(ndxQmin)<BusStruct.Uref0(ndxQmin)));
            if ~isempty(ndxQminSwitch)
                ndxToSwitch=ndxQminSwitch( find(BusStruct.NumOfPUPQSwitches(ndxQminSwitch)<=1) );%����������� ������ �� ����, ������� ������������� �� ����� ������ 2 ���
                if ~isempty(ndxToSwitch)
                    if strcmp(display,'on')
                        disp(['Generator at Bus ',mat2str(bus(ndxQminSwitch)),' PQ->PU switching']);
                    end
                    bus_type(ndxQminSwitch) = 2*ones(length(ndxQminSwitch),1);
                    BusStruct.NumOfPUPQSwitches(ndxToSwitch)=BusStruct.NumOfPUPQSwitches(ndxToSwitch)+ones(length(ndxToSwitch),1);%����������� ������� ������������ �� �������
                end
            end
        end
    end
end



if ~isempty(qg_max_idx)%Qg>Qgmax    
    if strcmp(display,'on')            
        disp(['Generator at Bus ',mat2str(bus(qg_max_idx)),' reached Qmax, PU->PQ switching']);        
    end
    bus_type(qg_max_idx) = 3*ones(length(qg_max_idx),1);        
    %Qg(qg_max_idx) = zeros(length(qg_max_idx),1);    
    Qg(qg_max_idx) = qg_max(qg_max_idx);%��������� �������� ��������� PU ���� �� ���������, ���� ���������� � ��������� PQ                  
    lim_flag = 1;
end
if ~isempty(qg_min_idx)       
    if strcmp(display,'on')            
        disp(['Generator at Bus ',mat2str(bus(qg_min_idx)),' reached Qmin, PU->PQ switching']);       
    end
    bus_type(qg_min_idx) = 3*ones(length(qg_min_idx),1);
    %set Qg to zero �������� ��������� ���������� �������� � ����
    Qg(qg_min_idx) = qg_min(qg_min_idx);             
    lim_flag = 1;
end
return
