function NetData=NetData()
%Settings Œ¡ﬂ«¿“≈À‹ÕŒ «¿ƒ¿“‹ ¬—≈ «Õ¿◊≈Õ»ﬂ
NetData.Settings.Display=0;%1/0 - display / not display iterations
NetData.Settings.epsDistLF=1e-4;%distribution system PF epsilon
NetData.Settings.IterMaxDistLF=30;%Maximum number of iterations for distribution system PF
NetData.Settings.epsTrLF=1e-4;%transmission system PF epsilon
NetData.Settings.IterMaxTrLF=30;%Maximum number of iterations for transmission system PF
NetData.Settings.epsLF=1e-3;%Global PF solution epsilon
NetData.Settings.IterMaxLF=30;%Maximum number of iterations for global PF solution

NetData.TDataFileName='Test213RAW';%Transmission data file name (Newton algorithm)
NetData.DnameTnode={'case33_1',214;...%every row contains : 1.Distribution network file name 2.Node number in transm. network
    'case33_2',215;
    'case33_3',216;
    'case33_4',217;
    'case33_5',218;
    'case33_6',219;};
