% Script for running simulation examples of the test scenarios implemented
% in the AAU testbed
% clc; clear;


addpath(genpath('DiSC'));
addpath('data');

%%% REMEMBER TO COMMENT THESE OUT !!!! %%%%
% reactiveComm = false;
% DelayMultiplier = 0.1;
% numDays = 1;
% seed = 2;



%% General setup
Ts = 5;                        % Sampling time system [s].
N = floor((1439*60/Ts)*numDays); % Number of iterations [-].

Ts_ctrl = 120;                 % Sampling time controller [s].
Ts_AssetData = Ts_ctrl;             % How often does the Asset transmit data [s]
startDay = 160;                 % Julian start day of the year (1-365) [-].
onPU = false;                   % To indicate if the system should be simulated on a per unit basis [-].
rng(seed);                         % Seed random generator [-].

% Choose scenario (1-4)
% - Scenario 1: No control
% - Scenario 2: LVGC active
scenario = 2;
withRefrigerators = true;
withHeatPumps = true;
penetration = 4; % 2 = 50%, 3=33%, 4 = 25%, ...

%% COMMUNICATION SETTINGS
withCom = true;
comLossPr = 0;

%% Setup electrical grid
% Base voltages
MV_vBase = 20e3;
LV_vBase = 0.4e3;

% Grid admittance matrix
Z = LVlineImp(true);

%% Setup power flow solver
% Setup the busses
numBus = max(max(Z(:,2)));
param.type = [0 ones(1,numBus-1)]; % 0 = slack, 1 = PQ
param.vBase = [MV_vBase LV_vBase*ones(1,numBus-1)];
pFlow = powerFlow(param);
% Construct grid admittance matrix
Y = pFlow.lDataToY(Z);
Yorg = Y;   % Make copy for trafo and switch implementation

%% Setup on-load tap-changing transformers
% Low voltage grid
param.vBase = LV_vBase;
param.nomTapRatio = 50;
param.Z = 0.004+1i*0.04;
param.Ts = Ts;
param.onPU = onPU;
% Create LV transformer
LV_tc1 = tcAsset(param);
% Placement of transformer
LV_tc1.busFrom = 1;
LV_tc1.busTo = 2;

% Input transformer objects into admittance matrix
% LV_tc1
[y1,y2,y3] = LV_tc1.sample(LV_vBase,LV_vBase,0,1,1);
Y(LV_tc1.busFrom,LV_tc1.busFrom) = Yorg(LV_tc1.busFrom,LV_tc1.busFrom) + y1+y2;
Y(LV_tc1.busTo,LV_tc1.busTo) = Yorg(LV_tc1.busTo,LV_tc1.busTo) +y1+y3;
Y(LV_tc1.busFrom,LV_tc1.busTo) = -y1;
Y(LV_tc1.busTo,LV_tc1.busFrom) = -y1;

%% Load consumption data
% Power factors
pfResi = 0.97;

% Low voltage grid (15 min. sampling)
intVal = 1:ceil((N*Ts)/(60*4))+1;
HouseData = load('data/house1to200days14startSample17280.mat');  % Mat file containing consumption data
HouseData.Data.HouseP = repmat(HouseData.Data.HouseP,numDays+3,1);
t = 0:size(HouseData.Data.HouseP(intVal,:),1)-1;
ti = 0:(size(HouseData.Data.HouseP(intVal,:),1)/(60*60*(1/Ts)*size(HouseData.Data.HouseP(intVal,:),1)/4)):size(HouseData.Data.HouseP(intVal,:),1)-1;
HouseData.pTs = interp1(t,HouseData.Data.HouseP(intVal,:),ti);
HouseData.pTs = HouseData.pTs.*1000;

% Refrigerators (5 min. sampling)
refrigData = load('data/refrigeratorData.mat');
t = 0:size(refrigData.refrigerator.realPower,1)-1;
ti = 0:(size(refrigData.refrigerator.realPower,1)/(60*60*(1/Ts)*size(refrigData.refrigerator.realPower,1)/12)):size(refrigData.refrigerator.realPower,1)-1;
refrigData.pTs = interp1(t,refrigData.refrigerator.realPower(:,:),ti);
refrigData.pTs = repmat(refrigData.pTs,numDays+3,1); %ceil(numDays/6)

% Heat pumps (5 min. sampling)
heatPumpData = load('data/heatPumpRealPower.mat');
t = 0:size(heatPumpData.hpPower,1)-1;
ti = 0:(size(heatPumpData.hpPower,1)/(60*60*(1/Ts)*size(heatPumpData.hpPower,1)/12)):size(heatPumpData.hpPower,1)-1;
heatPumpData.pTs = interp1(t,heatPumpData.hpPower(:,:),ti);
heatPumpData.pTs = repmat(heatPumpData.pTs,numDays+3,1);%ceil(numDays/6)

%% Set active and reactive power on each bus (Pin,Qin)
Pin = zeros(N,numBus);
Qin = zeros(N,numBus);
% Low voltage grid
if withRefrigerators == true
    temp = 1;
else
    temp = 0;
end

if withHeatPumps == true
    temp2 = 1;
else
    temp2 = 0;
end
% Active Power
Pin(:,3) = -sum(HouseData.pTs(1:N,1:8),2)-sum(refrigData.pTs(1:N,1:8),2)*temp;
Pin(:,4) = -sum(HouseData.pTs(1:N,9:12),2)-sum(refrigData.pTs(1:N,9:12),2)*temp;
Pin(:,5) = -sum(HouseData.pTs(1:N,13:15),2)-sum(refrigData.pTs(1:N,13:15),2)*temp;
Pin(:,6) = -sum(HouseData.pTs(1:N,16:19),2)-sum(refrigData.pTs(1:N,16:19),2)*temp;
Pin(:,7) = -sum(HouseData.pTs(1:N,20:22),2)-sum(refrigData.pTs(1:N,20:22),2)*temp;
Pin(:,8) = -sum(HouseData.pTs(1:N,23:27),2)-sum(refrigData.pTs(1:N,23:27),2)*temp;
Pin(:,10) = -sum(HouseData.pTs(1:N,28:31),2)-sum(refrigData.pTs(1:N,28:31),2)*temp;
Pin(:,11) = -sum(HouseData.pTs(1:N,32:35),2)-sum(refrigData.pTs(1:N,32:35),2)*temp;
Pin(:,12) = -sum(HouseData.pTs(1:N,36:37),2)-sum(refrigData.pTs(1:N,36:37),2)*temp;
Pin(:,13) = -sum(HouseData.pTs(1:N,38:38),2)-sum(refrigData.pTs(1:N,38:38),2)*temp;
Pin(:,14) = -sum(HouseData.pTs(1:N,39:40),2)-sum(refrigData.pTs(1:N,39:40),2)*temp;
Pin(:,15) = -sum(HouseData.pTs(1:N,41:45),2)-sum(refrigData.pTs(1:N,41:45),2)*temp;
Pin(:,16) = -sum(HouseData.pTs(1:N,46:49),2)-sum(refrigData.pTs(1:N,46:49),2)*temp;
Pin(:,17) = -sum(HouseData.pTs(1:N,50:52),2)-sum(refrigData.pTs(1:N,50:52),2)*temp;
Pin(:,18) = -sum(HouseData.pTs(1:N,53:54),2)-sum(refrigData.pTs(1:N,53:54),2)*temp;
Pin(:,19) = -sum(HouseData.pTs(1:N,55:58),2)-sum(refrigData.pTs(1:N,55:58),2)*temp;
Pin(:,20) = -sum(HouseData.pTs(1:N,59:61),2)-sum(refrigData.pTs(1:N,59:61),2)*temp;
Pin(:,21) = -sum(HouseData.pTs(1:N,62:63),2)-sum(refrigData.pTs(1:N,62:63),2)*temp;
Pin(:,22) = -sum(HouseData.pTs(1:N,64:67),2)-sum(refrigData.pTs(1:N,64:67),2)*temp;
Pin(:,23) = -sum(HouseData.pTs(1:N,68:70),2)-sum(refrigData.pTs(1:N,68:70),2)*temp;
Pin(:,24) = -sum(HouseData.pTs(1:N,71:72),2)-sum(refrigData.pTs(1:N,71:72),2)*temp;
Pin(:,25) = -sum(HouseData.pTs(1:N,73:75),2)-sum(refrigData.pTs(1:N,73:75),2)*temp;
Pin(:,26) = -sum(HouseData.pTs(1:N,76:77),2)-sum(refrigData.pTs(1:N,76:77),2)*temp;
Pin(:,27) = -sum(HouseData.pTs(1:N,78:80),2)-sum(refrigData.pTs(1:N,78:80),2)*temp;
Pin(:,28) = -sum(HouseData.pTs(1:N,81:81),2)-sum(refrigData.pTs(1:N,81:81),2)*temp;
Pin(:,29) = -sum(HouseData.pTs(1:N,82:83),2)-sum(refrigData.pTs(1:N,82:83),2)*temp;
Pin(:,30) = -sum(HouseData.pTs(1:N,84:85),2)-sum(refrigData.pTs(1:N,84:85),2)*temp;
Pin(:,31) = -sum(HouseData.pTs(1:N,86:87),2)-sum(refrigData.pTs(1:N,86:87),2)*temp;
Pin(:,32) = -sum(HouseData.pTs(1:N,88:90),2)-sum(refrigData.pTs(1:N,88:90),2)*temp;
Pin(:,34) = -sum(HouseData.pTs(1:N,91:94),2)-sum(refrigData.pTs(1:N,91:94),2)*temp;
Pin(:,35) = -sum(HouseData.pTs(1:N,95:97),2)-sum(refrigData.pTs(1:N,95:97),2)*temp;
Pin(:,36) = -sum(HouseData.pTs(1:N,98:104),2)-sum(refrigData.pTs(1:N,98:104),2)*temp;
Pin(:,37) = -sum(HouseData.pTs(1:N,105:106),2)-sum(refrigData.pTs(1:N,105:106),2)*temp;
Pin(:,38) = -sum(HouseData.pTs(1:N,107:108),2)-sum(refrigData.pTs(1:N,107:108),2)*temp;
Pin(:,39) = -sum(HouseData.pTs(1:N,109:109),2)-sum(refrigData.pTs(1:N,109:109),2)*temp;
Pin(:,40) = -sum(HouseData.pTs(1:N,110:111),2)-sum(refrigData.pTs(1:N,110:111),2)*temp;
Pin(:,41) = -sum(HouseData.pTs(1:N,112:113),2)-sum(refrigData.pTs(1:N,112:113),2)*temp;
Pin(:,42) = -sum(HouseData.pTs(1:N,114:116),2)-sum(refrigData.pTs(1:N,114:116),2)*temp;
Pin(:,1:2:end) = Pin(:,1:2:end) - heatPumpData.pTs(1:N,1:21);

PinCons = Pin*1.1;
% Reactive power
Qin(:,:) = Pin.*tan(acos(pfResi));
Qin_cons = Qin;

%% Setup assets
% Low voltage grid
% PV systems
% Solar PV systems
LV_PVlocation = [3 16 17 20 24 25 34 37];%13:52;%[17 27 28 31 35 36 48 49];        % Set to empty for no PV system
LV_numPV = length(LV_PVlocation);
param.onPU = onPU;
param.vBase = LV_vBase;
param.pRated = 6e3;
param.sMax = 6e3;
param.eta = 0.22;
param.A = 28;
% Create LV PV objects
for i=1:LV_numPV
    LV_pv(i) = pvAsset(param);
end
% Solar irradiance
param.Ts = Ts;
param.lat = 56.889;         
param.t = 0.75;           
param.p = 100;
LV_si = solarIrradiance(param);

% Energy storages
LV_ESlocation = [3 20 34];
LV_numES = length(LV_ESlocation);
param.sBase = 0;
param.vBase = LV_vBase;
param.sMax = 10e3;
param.pRatedMax = 10e3;
param.pRatedMin = -10e3;
param.pRate = 100*1000; 
param.eRated = 65e3*60*60;
param.Ts = Ts;
param.onPU = onPU;
% Create LV ES objects
for i=1:LV_numES
    LV_es(i) = esAsset(param);
    LV_es(i).eOld = param.eRated/2; % Initialize energy storages to half capacity
    LV_es(i).setQmode(1);
end

%% Setup Controllers
pRef = load('data/pRefDataLVonly.mat');
pRef = -smooth(pRef.pSlack,60);%+5e3;
pRef(50:250)= pRef(50:250)+30e3;%+20 - no constraints met.
pRef(300:800) = pRef(300:800)-18e3;
pRef = zeros(N,1)-20e3;
pRef(400:850) = 10e3;
pRef(851:1050) = 20e3;
pRef(1051:end) = -10e3;

% Low voltage controller
lvgc_assetBusses = [3;20;34];
lvgc_numAssets = length(lvgc_assetBusses);
xAsset = zeros(N,length(lvgc_assetBusses)+1);
lvgc_L = [0 0 0 0;          % Laplacian of communication graph
         -1 1 0 0;
         -1 0 1 0;
         -1 0 0 1];
 
lvgc_sMaxAssets = zeros(lvgc_numAssets,1);
lvgc_sMinAssets = zeros(lvgc_numAssets,1);
for i=1:lvgc_numAssets
    lvgc_sMaxAssets(i) = LV_es(i).pRatedMax;
    lvgc_sMinAssets(i) = LV_es(i).pRatedMin;
end

% Insert LVGC object
param.numAssets = length(lvgc_assetBusses);
param.L = lvgc_L;
param.Ts = Ts_ctrl;
param.pAssetMax = lvgc_sMaxAssets;
param.pAssetMin = lvgc_sMinAssets;
param.vDroopGains = ones(lvgc_numAssets,1)*200000;

lvgc = lvgcM4_tfk(param);
% Set P and I gain of internal PI controller
lvgc.setGains(0.7,0.6); % 0.3 , 0.2; 0.7, 0.6
% Set hysterisis limits for voltage control
vMinStart = 0.98*LV_vBase;
vMinStop = 0.985*LV_vBase;
vMaxStart = 1.02*LV_vBase;
vMaxStop = 1.01*LV_vBase;
lvgc.setHystParam(vMaxStart,vMaxStop,vMinStart,vMinStop);

% LVGC PV indexing
lvgc_PVidx = zeros(lvgc_numAssets,1);
for i=1:lvgc_numAssets
    lvgc_PVidx(i) = find(LV_PVlocation==lvgc_assetBusses(i));
end

% Low pass filters
% lvgc_lpf_input = firstOrderLowPassFilter(Ts,10*Ts,lvgc_numAssets);

%% Setup communication links
% Setup communication link

numLinksIn = 6;
numLinksOut = 5;

param.Ts = Ts;                          % Sampling time
param.maxDelay = DelayMultiplier*Ts;                  % Maximum delay
param.numLinksIn = numLinksIn;        % Number of communication links
param.numLinksOut = numLinksOut;
param.mu = DelayMultiplier*Ts_AssetData;
param.sigma = 1*Ts;
% Create communication link object
for i=1:lvgc_numAssets;
    lvgc_cl(i) = comLink(param);
    lvgc_cl(i).setInverseCDFdelay('exp',param);%dataTrace
    lvgc_cl(i).setPrLoss(comLossPr);
end
% Delay data
% Delay from assets to LVGC
delayDataAssetsToLvgc = ones(lvgc_numAssets,N)*0;

% Delay from LVGC to assets
delayDataLvgcToAssets = ones(lvgc_numAssets,N)*2*Ts;

evDel = NetworkDelay(param);

%% Input to assets
% Low voltage
LV_tc1_vRef = LV_vBase;         % MV trafo voltage rereference
LV_tc1_uRef = 0;                % MV trafo tap position reference
LV_tc1.setMode(0);              % MV trafo set mode (0 = follow tap reference, 1 = automatic voltage control)
% PV systems
LV_pv_vRef = LV_vBase*ones(LV_numPV,1);
LV_pv_dP = zeros(LV_numPV,1);
LV_pv_dPlim = zeros(LV_numPV,1);
LV_pv_qRef = zeros(LV_numPV,1);
% Energy storages
LV_es_vRef = LV_vBase*ones(LV_numES,1);
LV_es_pRef = zeros(LV_numES,1);
LV_es_qRef = zeros(LV_numES,1);

%% Run simulation
% Allocate memory
% Grid
vOut = ones(N+1,numBus);        % Voltages at each bus
pSlack = zeros(N,1);            % Active power at slack bus
qSlack = zeros(N,1);            % Reactive power at slack bus

% Low voltage
% PV systems
LV_PV_p = zeros(N,LV_numPV);
LV_PV_q = zeros(N,LV_numPV);
LV_PV_pAvb = zeros(N,LV_numPV);
% ES systems
LV_ES_p = zeros(N,LV_numES);
LV_ES_q = zeros(N,LV_numES);
LV_ES_e = zeros(N,LV_numES);

% Control
% Low voltage
lvgc_uRefP = zeros(N,lvgc_numAssets);
lvgc_uRefQ = zeros(N,lvgc_numAssets);
lvgc_pFlexUpAgg = zeros(N,1);
lvgc_pFlexDownAgg = zeros(N,1);
lvgc_qFlexUpAgg = zeros(N,1);
lvgc_qFlexDownAgg = zeros(N,1);
lvgc_pMax = zeros(lvgc_numAssets,1);
lvgc_pMin = zeros(lvgc_numAssets,1);
lvgc_qMin = zeros(lvgc_numAssets,1);
lvgc_qMax = zeros(lvgc_numAssets,1);
lvgc_pMeas = zeros(N,1);
lvgc_qMeas = zeros(N,1);
lvgc_pFlexMax = zeros(N,1);
lvgc_pFlexMin = zeros(N,1);
lvgc_pv_lpf = zeros(N,lvgc_numAssets);

RunCtrlFlag = false;
reaction = false;
voltageEventStart = 0;
reactionDelay = 0;

if reactiveComm == true
    AssetTXevent = inf;
else
    AssetTXevent = 1;
end

ApplyCtrlEvent = inf;
ReqAssetDataEvent = inf;
mm = -1; %ones(1,3)*-1;
% lvgc_ComOut2 = zeros(6,3,1);
mmTimes = ones(1,3)*-1; %0;
mmSampleCounter = 0;
tic
for i=1:N
    %% Itterate day
    if ~mod(i,24*60*60/Ts)
       if startDay == 365
           startDay = 1;
       else
           startDay = startDay +1;
       end
    end

    if ~mod(i,1000)
        disp(100*i/N)
    end
    
%     if sum(lvgc_qMax-10000)
%         disp('woot');
%         disp(lvgc_qMax);
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set invironmental data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wMean = (9+3*sin(i*0.0004));
    cc = 0.1 + 0.05*sin(i*0.08);

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Insert assets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i>1
        % Low voltage grid
        % Tap changing transformers
        % LV_tc1
        [y1,y2,y3] = LV_tc1.sample(abs(vOut(i-1,LV_tc1.busTo)),LV_tc1_vRef,LV_tc1_uRef,i,startDay);
        Y(LV_tc1.busFrom,LV_tc1.busFrom) = Yorg(LV_tc1.busFrom,LV_tc1.busFrom) + y1+y2;
        Y(LV_tc1.busTo,LV_tc1.busTo) = Yorg(LV_tc1.busTo,LV_tc1.busTo) +y1+y3;
        Y(LV_tc1.busFrom,LV_tc1.busTo) = -y1;
        Y(LV_tc1.busTo,LV_tc1.busFrom) = -y1;
        
        % PV systems
        for j=1:LV_numPV
            SI = LV_si.sample(i,startDay,cc);
            [LV_PV_p(i,j),LV_PV_q(i,j),LV_PV_pAvb(i,j)] = LV_pv(j).sample(SI,abs(vOut(i-1,LV_PVlocation(j))),LV_pv_dP(j),LV_pv_dPlim(j),LV_pv_qRef(j),LV_pv_vRef(j));
            Pin(i,LV_PVlocation(j)) = Pin(i,LV_PVlocation(j)) + LV_PV_p(i,j);
            Qin(i,LV_PVlocation(j)) = Qin(i,LV_PVlocation(j)) + LV_PV_q(i,j);
        end
        % Energi storages
        for j=1:LV_numES
            [LV_ES_p(i,j),LV_ES_q(i,j),LV_ES_e(i,j)] = LV_es(j).sample(abs(vOut(i-1,LV_ESlocation(j))),LV_es_pRef(j),LV_es_qRef(j),LV_es_vRef(j));
            Pin(i,LV_ESlocation(j)) = Pin(i,LV_ESlocation(j)) + LV_ES_p(i,j);
            Qin(i,LV_ESlocation(j)) = Qin(i,LV_ESlocation(j)) + LV_ES_q(i,j);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulate Electrical Grid  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [vOut(i,:),pSlack(i),qSlack(i),~]=pFlow.nrLoadFlow(Y,Pin(i,:)',Qin(i,:)');
    if pFlow.nIte == 100
        disp('ERROR nrIte = 100')
        break
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply  Control 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Low Volotage Grid Control 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i > 1
        if scenario > 1
            % Measurement at substation
            lvgc_pMeas = -pSlack(i);
            % State of assets
            lvgc_PV_pAvb = LV_PV_pAvb(i,lvgc_PVidx);
            lvgc_pAssetState = LV_ES_p(i,:)'+lvgc_PV_pAvb';
            
            % Flexibility of assets (pMin, pMax)
            % Filter power available from PV systems
%             lvgc_pv_lpf(i,:) = lvgc_lpf_input.sampleFilter(LV_PV_pAvb(i,lvgc_PVidx(:))');
            for j=1:lvgc_numAssets
                % Active power flexibility
                if LV_ES_e(i,j) <= 0
                    lvgc_pMax(j) = 0+LV_PV_pAvb(i,lvgc_PVidx(j));
                    lvgc_pMin(j) = LV_es(j).pRatedMin;
                elseif LV_ES_e(i,j) >= LV_es(j).eRated
                    lvgc_pMax(j) = LV_es(j).pRatedMax+LV_PV_pAvb(i,lvgc_PVidx(j));% - LV_ES_p(i,j);
                    lvgc_pMin(j) = 0;
                else 
                    lvgc_pMax(j) = LV_es(j).pRatedMax+LV_PV_pAvb(i,lvgc_PVidx(j));
                    lvgc_pMin(j) = LV_es(j).pRatedMin;
                end
                
                % Reactive power flexibility
                lvgc_qMin(j) = sqrt(LV_es(j).sMax^2-LV_ES_p(i,j)^2);
                lvgc_qMax(j) = sqrt(LV_es(j).sMax^2-LV_ES_p(i,j)^2);
            end
            % Voltage measurement
            lvgc_vMeas = vOut(i,lvgc_assetBusses);
            

            %%% CHECK FOR OUT OF BOUND VOLTAGE!
            if reactiveComm == true
%                 disp(min(abs(lvgc_vMeas)))
                if RunCtrlFlag == true;
%                     disp('RunCtrlFlag true');
%                     disp(min(abs(lvgc_vMeas)));
%                     disp(max(abs(lvgc_vMeas)));
                    if max(abs(lvgc_vMeas))<=vMaxStop && min(abs(lvgc_vMeas))>=vMinStop
                        RunCtrlFlag = false;
                        ReqAssetDataEvent = inf;
                        ApplyCtrlEvent = inf;
                        AssetTXevent = inf;
%                         disp('ApplyCtrlEvent is inf')
                    end
                elseif max(abs(lvgc_vMeas))>=vMaxStart
                    RunCtrlFlag = true;
                    ReqAssetDataEvent = i+evDel.draw();
                elseif min(abs(lvgc_vMeas))<=vMinStart
%                     disp('Voltage out of bound')
                    RunCtrlFlag = true;
                    ReqAssetDataEvent = i+evDel.draw();
%                     disp(['voltage low event, reqAssetDataEvent-i: ' num2str(ReqAssetDataEvent-i)]);
                    
                end
            else
                if RunCtrlFlag == true;
                    if max(lvgc_vMeas)<=vMaxStop && min(lvgc_vMeas)>=vMinStop
                        RunCtrlFlag = false;
                        ApplyCtrlEvent = inf;
                    end
                elseif max(lvgc_vMeas)>=vMaxStart
                    RunCtrlFlag = true;
                    ApplyCtrlEvent = i+evDel.draw();
                elseif min(lvgc_vMeas)<=vMinStart && reaction == false
                    RunCtrlFlag = true;
                    ApplyCtrlEvent = i+evDel.draw();
                end                
            end
            if reactiveComm == true && i >= ReqAssetDataEvent
                AssetTXevent = i + evDel.draw();
                ReqAssetDataEvent = inf;
%                 disp(['setting AssetTXevent-i:' num2str(AssetTXevent-i)]);
                    
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Sample Communication from asset to LVGC
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if withCom == true && i >= AssetTXevent
                if reactiveComm == false
                    AssetTXevent = i+exprnd(Ts_AssetData/Ts);
                end
                for j=lvgc_numAssets:-1:1
                    % Form data vector to be communicated
                    in = [lvgc_pAssetState(j);lvgc_pMin(j);lvgc_pMax(j);lvgc_qMin(j);lvgc_qMax(j);lvgc_vMeas(j)];
                    b(j) = lvgc_cl(j).TxIn(i,in,delayDataAssetsToLvgc(j,i));
                end
                if reactiveComm == true
                    ApplyCtrlEvent = max(b);
%                     disp(['setting ApplyCtrlEvent-i:' num2str(ApplyCtrlEvent-i)]);
                    AssetTXevent = inf;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Sample LVGC with LVGC sampling time
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             RunCtrlFlag = true;
            if i >= ApplyCtrlEvent% && RunCtrlFlag%RunCtrlFlag && ~mod(Ts*i-voltageEventStart,Ts_ctrl)
%                 disp('i and apply')
%                 disp(i)
%                 disp(ApplyCtrlEvent);
                if reactiveComm == false
                    ApplyCtrlEvent = i + Ts_ctrl/Ts;
                else
                    ReqAssetDataEvent = i + Ts_ctrl/Ts+evDel.draw();
                    ApplyCtrlEvent = inf;
%                     disp(['setting ReqAssetDataEvent-i:' num2str(ReqAssetDataEvent-i)]);
                end
                    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Sample Communication from asset to LVGC
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if withCom == true
                    for j=1:lvgc_numAssets
                        % Form data vector to be communicated
                        lvgc_ComOut1(:,j,i+1) = lvgc_cl(j).RxIn(i);
                        for n=1:numLinksIn
                            if isnan(lvgc_ComOut1(n,j,i+1))
                                lvgc_ComOut1(n,j,i+1) = lvgc_ComOut1(n,j,i);
                            end
                        end
                    end
                    % Collect data after communication
                    input_lvgc_pAssetState = lvgc_ComOut1(1,:,i+1)';
                    input_lvgc_pMin = lvgc_ComOut1(2,:,i+1)';
                    input_lvgc_pMax = lvgc_ComOut1(3,:,i+1)';
                    input_lvgc_qMin = lvgc_ComOut1(4,:,i+1)';
                    input_lvgc_qMax = lvgc_ComOut1(5,:,i+1)';
                    input_lvgc_vMeas = lvgc_ComOut1(6,:,i+1);
                else
                    input_lvgc_pAssetState = lvgc_pAssetState;
                    input_lvgc_pMin = lvgc_pMin;
                    input_lvgc_pMax = lvgc_pMax;
                    input_lvgc_qMin = lvgc_qMin;
                    input_lvgc_qMax = lvgc_qMax;
                    input_lvgc_vMeas = lvgc_vMeas;
                end
                mmSampleCounter = mmSampleCounter + 1;
%                 disp('');
%                 disp(['mmSamplesCounter:' num2str(mmSampleCounter)]);
                
%                 mm(end+1,:) = abs(input_lvgc_vMeas-lvgc_vMeas);
%                 disp(['input_lvgc_vMeas:' num2str(input_lvgc_vMeas)]);
%                 disp(['lvgc_vMeas:' num2str(lvgc_vMeas)]);
%                 disp(['setting mm(end,:):' num2str(mm(end,:))]);
%                 disp('');
%                 size(mm)
%                 disp(['Applying control at: ' num2str(i)]);
                [lvgc_uRefP(i,:),lvgc_uRefQ(i,:),lvgc_pFlexUpAgg(i),lvgc_pFlexDownAgg(i),lvgc_qFlexUpAgg(i),lvgc_qFlexDownAgg(i)] =...
                    lvgc.sample(pRef(i),lvgc_pMeas,input_lvgc_pAssetState,input_lvgc_pMin,input_lvgc_pMax,input_lvgc_qMax,input_lvgc_qMin,input_lvgc_vMeas);
                lvgc_uRefP(i,:) = 0;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Sample Communication from LVGC to assets
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if withCom == true
%                     disp('dispatch');
                    for j=1:lvgc_numAssets
                        % Form data vector to be communicated
                        in = [lvgc_uRefP(i,j);lvgc_uRefQ(i,j);mmSampleCounter;input_lvgc_vMeas(j);i];%;input_lvgc_qMax
                        B = lvgc_cl(j).TxOut(i,in,delayDataLvgcToAssets(j,i));
%                         disp(['B: ' num2str(B)])
                    end
                end
            else
                lvgc_uRefP(i,:) = lvgc_uRefP(i-1,:);
                lvgc_uRefQ(i,:) = lvgc_uRefQ(i-1,:);
                lvgc_pFlexUpAgg(i) = lvgc_pFlexUpAgg(i-1);
                lvgc_pFlexDownAgg(i) = lvgc_pFlexDownAgg(i-1);
                lvgc_qFlexUpAgg(i) = lvgc_qFlexUpAgg(i-1);
                lvgc_qFlexDownAgg(i) = lvgc_qFlexDownAgg(i-1);
                
                input_lvgc_uRefP = lvgc_uRefP(i,:);
                input_lvgc_uRefQ = lvgc_uRefQ(i,:);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Sample Communication from LVGC to assets
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if withCom == true
                for j=1:lvgc_numAssets
                    % Form data vector to be communicated
                    in = [lvgc_uRefP(i,j);lvgc_uRefQ(i,j)];
                    lvgc_ComOut2(:,j,i+1) = lvgc_cl(j).RxOut(i);
                    for n=1:numLinksOut
                        if isnan(lvgc_ComOut2(n,j,i+1))
                            lvgc_ComOut2(n,j,i+1) = lvgc_ComOut2(n,j,i);
                        end
                    end
                end
                % Collect data after communication
                input_lvgc_uRefP = lvgc_ComOut2(1,:,i+1)';
                input_lvgc_uRefQ = lvgc_ComOut2(2,:,i+1)';
                input_lvgc_mmSampleCounter = lvgc_ComOut2(3,:,i+1)';
                input_asset_vMeas = lvgc_ComOut2(4,:,i+1);
%                 input_asset_qMax =  lvgc_ComOut2(6:8,:,i+1);
                lvgc_uRefP(i,:) = input_lvgc_uRefP;
                lvgc_uRefQ(i,:) = input_lvgc_uRefQ;
%                 disp('input_lvgc_mmSampleCounter: ');
%                 disp(input_lvgc_mmSampleCounter');
%                 disp('lvgc_ComOut2(3,:,i): ');
%                 disp(lvgc_ComOut2(3,:,i));
%                 temp = length(mm);
                for n = 1:lvgc_numAssets
                    if input_lvgc_mmSampleCounter(n) > lvgc_ComOut2(3,n,i)
                        mm(end+1) = abs(input_asset_vMeas(n)-lvgc_vMeas(n));
%                         mm(end+1) = abs(input_asset_qMax(n)-lvgc_qMax(n));
                        mmTimes(end+1,:) = [i lvgc_ComOut2(5,n,i+1) n];
                    end
                end
            else
                input_lvgc_uRefP = lvgc_uRefP(i,:);
                input_lvgc_uRefQ = lvgc_uRefQ(i,:);
            end

            % Dispatch reference to energy storages and PV systems
            for j=1:lvgc_numAssets
                if input_lvgc_uRefP(j)>=lvgc_PV_pAvb(j) && LV_ES_e(i,j)>0
                    LV_es_pRef(j) = input_lvgc_uRefP(j)-lvgc_PV_pAvb(j);
                    LV_pv_dPlim(lvgc_PVidx(j)) = 0;
                elseif input_lvgc_uRefP(j)>=lvgc_PV_pAvb(j) && LV_ES_e(i,j)>=LV_es(j).eRated
                    LV_es_pRef(j) = 0;
                    LV_pv_dPlim(lvgc_PVidx(j)) = 0;
                elseif input_lvgc_uRefP(j)>=lvgc_PV_pAvb(j) && LV_ES_e(i,j)<=0
                    LV_es_pRef(j) = 0;
                    LV_pv_dPlim(lvgc_PVidx(j)) = 0;
                elseif input_lvgc_uRefP(j) <= lvgc_PV_pAvb(j) && LV_ES_e(i,j)<LV_es(j).eRated
                    LV_es_pRef(j) = input_lvgc_uRefP(j)-lvgc_PV_pAvb(j);
                    LV_pv_dPlim(lvgc_PVidx(j)) = 0;
                elseif input_lvgc_uRefP(j) <= lvgc_PV_pAvb(j) && LV_ES_e(i,j)>=LV_es(j).eRated
                    LV_pv_dP(lvgc_PVidx(j)) = -(lvgc_PV_pAvb(j)-input_lvgc_uRefP(j));
                    LV_es_pRef(j) = 0;
                else
                    disp('hmmm')
                end

            end

            LV_es_qRef = input_lvgc_uRefQ;
            
            % Aggregated power flexibility to MVGC
            lvgc_pFlexMax(i+1) = lvgc_pMeas + lvgc_pFlexUpAgg(i);
            lvgc_pFlexMin(i+1) = lvgc_pMeas + lvgc_pFlexDownAgg(i);
             
        end
    end
end
mm = mm(5:end);
mmTimes = mmTimes(5:end,:);
toc

%% Plotting
close all;
tvec = (0:N-1);%/60*60/Ts;

% Voltages
figure
plot(tvec,abs(vOut(1:N,lvgc_assetBusses))/LV_vBase)
ylabel('Voltage [p.u.]')
title('LV Grid')
xlabel('Time [hrs]')

% PQ slack
figure
subplot(2,1,1)
plot(tvec,pRef(1:N)/1e3,tvec,-pSlack(1:N)/1e3)
ylabel('Power [kW]')
legend('P Ref','P')
title('Active Power')
subplot(2,1,2)
plot(tvec,qSlack/1e3)
ylabel('Power [kVar]')
xlabel('Time [hrs]')
title('Reactive Power')

% Low Voltage Assets
% Energy storages
figure
subplot(3,1,1)
plot(tvec,LV_ES_p/1e3)
ylabel('Power [kW]')
title('Energy Storages')
subplot(3,1,2)
plot(tvec,LV_ES_q/1e3)
ylabel('Power [kVar]')
subplot(3,1,3)
plot(tvec,LV_ES_e/param.eRated)
ylabel('SoC [-]')
xlabel('Time [hrs]')

% PV systems
figure
subplot(3,1,1)
plot(tvec,LV_PV_pAvb(:,lvgc_PVidx(1)),tvec,LV_PV_p(:,lvgc_PVidx(1)),tvec,LV_PV_q(:,lvgc_PVidx(1)))
title('PV Asset 1')
legend('P Available','P Produced','Q Produced')
ylabel('Power [kW/kVar]')
subplot(3,1,2)
plot(tvec,LV_PV_pAvb(:,lvgc_PVidx(2)),tvec,LV_PV_p(:,lvgc_PVidx(2)),tvec,LV_PV_q(:,lvgc_PVidx(2)))
title('PV Asset 2')
ylabel('Power [kW/kVar]')
subplot(3,1,3)
plot(tvec,LV_PV_pAvb(:,lvgc_PVidx(3)),tvec,LV_PV_p(:,lvgc_PVidx(3)),tvec,LV_PV_q(:,lvgc_PVidx(3)))
title('PV Asset 3')
ylabel('Power [kW/kVar]')
xlabel('Time [hrs]')
