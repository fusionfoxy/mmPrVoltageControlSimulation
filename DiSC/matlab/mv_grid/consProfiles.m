% Consumption and production profiles for MV grid

% P and Q values for 24 hours sampled every hour (Summer period)
% Industry
MV_Pindu = [41 42 45.72 76.6 58.2 97.96 107.44 117.84 111.64 115 111.84 115.76 117.68 121.2 112.2 105.04 91.84 ...
                60.2 58 39.88 38.24 38.12 37.88 41.08]'.*1000*2;
% Interpolate to 15 min
t = 0:length(MV_Pindu)-1;
ti = 0:(1/4):length(MV_Pindu)-1;
MV_Pindu = interp1(t,MV_Pindu(:,:),ti);
MV_Pindu(94:96) = MV_Pindu(end); 
MV_Pindu = MV_Pindu';
MV_Qindu = MV_Pindu*tan(acos(pfIndu));
%%
% Agriculture
MV_Pagri = [4.18 3.55 3.59 3.66 4.02 4.06 4.37 4.09 5.22 6.05 5.67 5.03 3.47 5.47 3.41 3.6 5.18 7.1 6.7 5.61 ...
                5.28 5.41 6.09 4.05]'.*1000*3;
t = 0:length(MV_Pagri)-1;
ti = 0:(1/4):length(MV_Pagri)-1;
MV_Pagri = interp1(t,MV_Pagri(:,:),ti);
MV_Pagri(94:96) = MV_Pagri(end);
MV_Pagri = MV_Pagri';
MV_Qagri = MV_Pagri*tan(acos(pfAgri));
% Residential
% PresiSummer = [61.29 57.14 55.07 54.85 63.21 80.43 101.23 104.58 104.49 86.58 85.61 77.35 81.45 83.26 85.52 98.8 ...
%                 142.74 187.92 168.89 136.94 117.25 114.58 100.53 86.58]*1000;
% QresiSummer = PresiSummer*tan(acos(pfResi)); 

% Commercial
MV_Pcomm = [13.62 12.88 12.73 12.95 12.59 12.88 13.77 14.74 21.92 20.51 20.13 18.18 24.17 24.34 25.76 26.27 31.17 ...
                22.31 16.22 17.12 15.52 15.88 15.16 15.44]'.*1000*4;
t = 0:length(MV_Pcomm)-1;
ti = 0:(1/4):length(MV_Pcomm)-1;
MV_Pcomm = interp1(t,MV_Pcomm(:,:),ti);
MV_Pcomm(94:96) = MV_Pcomm(end);
MV_Pcomm = MV_Pcomm';
MV_Qcomm = MV_Pcomm*tan(acos(pfComm));

% Wind power plant
MV_PWP = (8e3+3e3.*randn(1,length(MV_Pindu)))'*0;
MV_QWP = MV_PWP*tan(acos(0.9));

% PV power plant
MV_PPV = (8e3+4e3.*randn(1,length(MV_Pindu)))'*0;
MV_QPV = MV_PPV*tan(acos(0.9));

