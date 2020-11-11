function CurDistNet=reordernodes(CurDistNet,mainnodenum)
branch=CurDistNet.acline;
bus=CurDistNet.bus;
%%
%Функция перенумерации графа тносительно балансирующего узла (mainnodenum)
%branch - исходная информация по ветвям схемы; mainnodenum - номер
%балансирующего узла (связь с сетью ВН); bus - исходная информация по узлам
%схемы.
%NodeToInt - матрица соответствия нумерации исходный узел (bus) -> новый узел (последовательная нумерация);
%IntToNode - матрица соответствия нумерации новый узел (последовательная нумерация) -> исходный узел (bus); 
%ndxBusToInt - соответствие индексов узлов в матрице bus индексам узлов при последовательной нумерации
%ndxBranchToInt - соответствие индексов ветвей в матрице branch инексам ветвей при воследовательной нумерации 
%L, G, Zv - матрицы, участвующие в расчете BF-методом.
%%
NumOfBr=size(branch,1);
L=sparse(1:1:NumOfBr,1:1:NumOfBr,-1*ones(NumOfBr,1),NumOfBr,NumOfBr);%
Zv=L;%Матрица Zv
NodeToInt=zeros(max(bus(:,1)),1);%матрица соответствия нумерации исходный узел -> новый узел
NodeToInt(mainnodenum)=1;
IntToNode=zeros(size(bus,1),1);%матрица соответствия нумерации новый узел -> исходный узел
IntToNode(1)=mainnodenum;%матрица соответствия нумерации новый узел -> исходный узел
ndxBranchToInt=zeros(size(branch,1),1);%матрица соответствия индексов ветвей в матрице branch инексам ветвей при последовательной нумерации
%Изменяем порядок узлов в дереве
nodes=mainnodenum;
counter=2;
ndxsToSearch=1:size(branch,1);%Индексы на которых выполняем поиск 
while ~isempty(ndxsToSearch)
    nodesnew=[];
    for count=1:length(nodes)
        ndx11=find(branch(ndxsToSearch,1)==nodes(count));
        ndx1=ndxsToSearch(ndx11);%Находим все связи, в которых есть узел node на первом месте
        if ~isempty(ndx1)
            ndx1sze=length(ndx1);
            NodeToInt(branch(ndx1,2))=counter:1:counter+ndx1sze-1;
            IntToNode(counter:1:counter+ndx1sze-1)=branch(ndx1,2);
            Zv(counter-1:1:counter+ndx1sze-2,counter-1:1:counter+ndx1sze-2)=diag(branch(ndx1,3)+1i*branch(ndx1,4));
            nodesnew=[nodesnew;branch(ndx1,2)];          
            if nodes(count)~=mainnodenum%Первый узел исключаем из матрицы L
                L(counter-1:1:counter+ndx1sze-2,NodeToInt(nodes(count))-1)=ones(ndx1sze,1);%Заполняем матрицу L
            end
            counter=counter+ndx1sze;
            ndxsToSearch(ndx11)=[];%Удаляем индексы уже найденных элементов, чтобы ускорить поиск
        end
        ndx22=find(branch(ndxsToSearch,2)==nodes(count));
        ndx2=ndxsToSearch(ndx22);%Находим все связи, в которых есть узел node на втором месте
        if ~isempty(ndx2)
            ndx2sze=length(ndx2);
            NodeToInt(branch(ndx2,1))=counter:1:counter+ndx2sze-1;
            IntToNode(counter:1:counter+ndx2sze-1)=branch(ndx2,1);
            Zv(counter-1:1:counter+ndx2sze-2,counter-1:1:counter+ndx2sze-2)=diag(branch(ndx2,3)+1i*branch(ndx2,4));
            nodesnew=[nodesnew;branch(ndx2,1)];
            if nodes(count)~=mainnodenum%Первый узел исключаем из матрицы L
                L(counter-1:1:counter+ndx2sze-2,NodeToInt(nodes(count))-1)=ones(ndx2sze,1);%Заполняем матрицу L
            end
            counter=counter+ndx2sze;            
            ndxsToSearch(ndx22)=[];%Удаляем индексы уже найденных элементов, чтобы ускорить поиск
        end
    end
    nodes=nodesnew;%Номера узлов для поиска на следующем уровне
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