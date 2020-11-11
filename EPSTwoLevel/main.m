%Example
DataPath=[pwd,'\Data'];%Add a folder with data
addpath(DataPath);
NetData=feval('NetData');%Network Data

%NetData initialization
NetData=InitData(NetData);

% PF solution
% t=zeros(20,1);
% for count=1:length(t)
% tic
[NetData,IsConverged]=globalloadflowparall(NetData);%with parallel computing of distr subsystems
%[NetData,IsConverged]=globalloadflowparall(NetData);%without parallel computing of distr subsystems
% t(count)=toc;
% end
% disp(['tmean=',num2str(mean(t(2:length(t)))),'']);

% Net1=NetData;%Copy NetData structure
% ndx=find(Net1.TrData.acline(:,1)==3 & Net1.TrData.acline(:,2)==4);
% Net1.TrData.acline(ndx,12)=0;%Trip line node5-node4
% [Net1,IsConverged]=globalloadflow(Net1);

[DistPowerFlows,P_loss,Q_loss]= DistLinePQ(NetData.DData.case33_2.Vdist_sol,NetData.DData.case33_2);

