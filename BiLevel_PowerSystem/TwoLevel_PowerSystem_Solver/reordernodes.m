function CurDistNet=reordernodes(CurDistNet,mainnodenum)
branch=CurDistNet.acline;
bus=CurDistNet.bus;
%%
%������� ������������� ����� ����������� �������������� ���� (mainnodenum)
%branch - �������� ���������� �� ������ �����; mainnodenum - �����
%�������������� ���� (����� � ����� ��); bus - �������� ���������� �� �����
%�����.
%NodeToInt - ������� ������������ ��������� �������� ���� (bus) -> ����� ���� (���������������� ���������);
%IntToNode - ������� ������������ ��������� ����� ���� (���������������� ���������) -> �������� ���� (bus); 
%ndxBusToInt - ������������ �������� ����� � ������� bus �������� ����� ��� ���������������� ���������
%ndxBranchToInt - ������������ �������� ������ � ������� branch ������� ������ ��� ���������������� ��������� 
%L, G, Zv - �������, ����������� � ������� BF-�������.
%%
NumOfBr=size(branch,1);
L=sparse(1:1:NumOfBr,1:1:NumOfBr,-1*ones(NumOfBr,1),NumOfBr,NumOfBr);%
Zv=L;%������� Zv
NodeToInt=zeros(max(bus(:,1)),1);%������� ������������ ��������� �������� ���� -> ����� ����
NodeToInt(mainnodenum)=1;
IntToNode=zeros(size(bus,1),1);%������� ������������ ��������� ����� ���� -> �������� ����
IntToNode(1)=mainnodenum;%������� ������������ ��������� ����� ���� -> �������� ����
ndxBranchToInt=zeros(size(branch,1),1);%������� ������������ �������� ������ � ������� branch ������� ������ ��� ���������������� ���������
%�������� ������� ����� � ������
nodes=mainnodenum;
counter=2;
ndxsToSearch=1:size(branch,1);%������� �� ������� ��������� ����� 
while ~isempty(ndxsToSearch)
    nodesnew=[];
    for count=1:length(nodes)
        ndx11=find(branch(ndxsToSearch,1)==nodes(count));
        ndx1=ndxsToSearch(ndx11);%������� ��� �����, � ������� ���� ���� node �� ������ �����
        if ~isempty(ndx1)
            ndx1sze=length(ndx1);
            NodeToInt(branch(ndx1,2))=counter:1:counter+ndx1sze-1;
            IntToNode(counter:1:counter+ndx1sze-1)=branch(ndx1,2);
            Zv(counter-1:1:counter+ndx1sze-2,counter-1:1:counter+ndx1sze-2)=diag(branch(ndx1,3)+1i*branch(ndx1,4));
            nodesnew=[nodesnew;branch(ndx1,2)];          
            if nodes(count)~=mainnodenum%������ ���� ��������� �� ������� L
                L(counter-1:1:counter+ndx1sze-2,NodeToInt(nodes(count))-1)=ones(ndx1sze,1);%��������� ������� L
            end
            counter=counter+ndx1sze;
            ndxsToSearch(ndx11)=[];%������� ������� ��� ��������� ���������, ����� �������� �����
        end
        ndx22=find(branch(ndxsToSearch,2)==nodes(count));
        ndx2=ndxsToSearch(ndx22);%������� ��� �����, � ������� ���� ���� node �� ������ �����
        if ~isempty(ndx2)
            ndx2sze=length(ndx2);
            NodeToInt(branch(ndx2,1))=counter:1:counter+ndx2sze-1;
            IntToNode(counter:1:counter+ndx2sze-1)=branch(ndx2,1);
            Zv(counter-1:1:counter+ndx2sze-2,counter-1:1:counter+ndx2sze-2)=diag(branch(ndx2,3)+1i*branch(ndx2,4));
            nodesnew=[nodesnew;branch(ndx2,1)];
            if nodes(count)~=mainnodenum%������ ���� ��������� �� ������� L
                L(counter-1:1:counter+ndx2sze-2,NodeToInt(nodes(count))-1)=ones(ndx2sze,1);%��������� ������� L
            end
            counter=counter+ndx2sze;            
            ndxsToSearch(ndx22)=[];%������� ������� ��� ��������� ���������, ����� �������� �����
        end
    end
    nodes=nodesnew;%������ ����� ��� ������ �� ��������� ������
end
G=inv(L);
ndx2=NodeToInt(bus(:,1));
[~, ndxBusToInt]=sort(ndx2);

%ndx3=NodeToInt(bus(:,1));


CurDistNet.NodeToInt=NodeToInt;
CurDistNet.L=L;
CurDistNet.G=G;
CurDistNet.Zv=Zv;
CurDistNet.IntToNode=IntToNode;
CurDistNet.ndxBusToInt=ndxBusToInt;
end