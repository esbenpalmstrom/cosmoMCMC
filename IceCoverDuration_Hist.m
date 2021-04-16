%IceCoverDuration_Hist.m

%Todo: Normalize, so it shows number of models with ice compared to total
%number of models.

clear; close all;
load('data/d18Ocurves.mat');
%choose which models to load:
models = [4];
eem = 116;
eemt = zeros(1,eem);
LGM = 26;
LGMt = zeros(1,LGM);

d18O_eem = d18O_10ky(1:eem);


LGM_iceduration = [];
eem_iceduration = [];

for i = 1:length(models)
    path = ['models/' num2str(models(i)) '/gausta_data_2.mat'];
    load(path);
    
    for j = 1:length(model.walker)
        
        u = model.walker{1,j}.u;
        I = find(model.walker{1,j}.status==1);
        d18O_T = u(1,I);
        Tdg = uint8(u(2,I)'*1000);
        
        
        %1 means ice in this timeslice, 0 is no ice
        LGM_ice = d18O_T' < d18O_10ky(1:LGM);
        
        for ii = 1:length(LGM_ice)
            LGM_ice(ii,1:Tdg(ii)+1) = 0;
        end
        
        
        LGMt = LGMt + sum(LGM_ice);
        LGM_ID = sum(LGM_ice,2);
        LGM_iceduration = cat(1,LGM_iceduration,LGM_ID);
        
        eem_ice = d18O_T' < d18O_10ky(1:eem);
        
        for ii = 1:length(eem_ice)
            eem_ice(ii,1:Tdg(ii)+1) = 0;
        end
        
        eemt = eemt + sum(eem_ice);
        eem_ID = sum(eem_ice,2);
        eem_iceduration = cat(1,eem_iceduration,eem_ID);        
        
    end
    
    figure; hold on; grid on;
    bar(eemt)
    xlabel('time'), ylabel('no. of models with ice')
    title('Number of models with ice at given time')
    
    figure; hold on; grid on;
    bar(eemt)
    xlabel('time'), ylabel('no. of models with ice')
    xlim([0 50]);
    title('Number of models with ice at given time')
    
    figure; hold on; grid on;
    bar(LGMt)
    xlabel('time, kyr'), ylabel('no. of models with ice')
    title('Number of models with ice at given time')
    
    %make plot with total ice cover duration for given time intervals
    figure; hold on; grid on;
    histogram(LGM_iceduration);
    xlabel('kyrs with ice cover since LGM')
    title('total ice cover duration frequency')
    
    figure; hold on; grid on;
    histogram(eem_iceduration);
    xlabel('kyrs with ice cover since eem');
    title('total ice cover duration frequency')
    
end
distFig
