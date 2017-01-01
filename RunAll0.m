close all; clear all; clc;

numDays = 1;

for gg = 2:2
    for ggg = [1:50]%11:30
        for g = [0:0.05:0.2 0.3:0.1:1 1.3 1.6]%[0:0.05:0.2 0.3:0.1:1 1.5 2]
            clearvars -except g gg ggg numDays
            if gg == 1
                reactiveComm = true;
            elseif gg == 2
                reactiveComm = false;
            end
            seed = ggg;
            DelayMultiplier = g;
            filename = char(strcat('results/5sec_0985lim_delay',num2str(DelayMultiplier),'_reactive',num2str(reactiveComm),'_seed',num2str(seed),'.mat'));
            disp(['Running: ' filename]) 
	    if exist(filename,'file')
		disp('skipping file')
		continue
	    end
	    fid = fopen( filename, 'wt' );
	    fclose(fid);
	    LVgridWithComm;
	    save(filename,'vOut','N','lvgc_assetBusses','LV_vBase','mm','tvec','Ts_AssetData')
%             save(filename)
            disp(filename)
        end
    end
end
