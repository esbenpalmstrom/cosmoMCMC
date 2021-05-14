%IceCoverDuration_Hist.m

%Todo: Normalize, so it shows number of models with ice compared to total
%number of models.

clear; close all;
load('data/d18Ocurves.mat');
%choose which models to load:
models = [17];
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
    
    figure; hold on;
    sbplt1 = subplot(2,2,1);
    box on;
    %set(gcf,'Position',[100 100 900 750]);
    eemt_norm = (eemt./length(eem_iceduration));
    bar(eemt_norm,'BarWidth',1.0)
    xlabel('time'), ylabel('Likelihood of models to have ice at given time')
    title('')
    grid on;
    
    sbplt2 = subplot(2,2,2);
    box on;
    hold on; grid on;
    %set(gcf,'Position',[100 100 900 750]);
    bar(eemt_norm,'BarWidth',1.0)
    xlabel('time'), ylabel('Likelihood of models to have ice at given time')
    xlim([0 50]);
    title('')
    
    sbplt3 = subplot(2,2,3);
    box on;
    hold on; grid on;
    %set(gcf,'Position',[100 100 900 750]);
    LGMt_norm = (LGMt./length(LGM_iceduration));
    bar(LGMt_norm,'BarWidth',1.0)
    xlabel('time, kyr'), ylabel('Likelihood of models to have ice at given time')
    title('')
    
    
    %Ice cover duration freq. for models since LGM. Maybe irrelevant?
%     figure; hold on; grid on;
%     set(gcf,'Position',[100 100 900 750]);
%     %[LGMCount,LGMValues] = hist(LGM_iceduration(:),26);
%     %LGMCount = LGMCount./sum(LGMCount);
%     %bar(LGMCount,LGMValues,'BarWidth',1.0)
%     histogram(LGM_iceduration,'Normalization','probability','FaceAlpha',1);
%     xlabel('kyrs with ice cover since LGM'), ylabel('Likelihood of total ice cover duration of all accepted models')
%     title('total ice cover duration frequency')
%     %set(gca, 'YScale', 'log')
%     set(gcf,'color','white')
    
    
    
    %ice cover duration freq. for all models
    sbplt4 = subplot(2,2,4);
    box on;
    hold on; grid on;
    set(gcf,'Position',[100 100 1200 750]);
    [vCount, vValues] = hist(eem_iceduration(:),62);
    vCount = vCount./sum(vCount);
    bar(vValues,vCount,'BarWidth',1.0)
    r = rectangle('Position',[1.2 0 sbplt4.XLim(2) 0.04],'LineWidth',1.5);
    %histogram(eem_iceduration,'Normalization','probability');
    xlabel('kyrs with ice cover since eem'), ylabel('Likelihood of total ice cover duration of all accepted models')
    title('total ice cover duration frequency')
    set(gcf,'color','white')
    
    %sumcurve = cumsum(vCount,'reverse');
    %plot(vValues,sumcurve,'-b')
    
    
    ax = axes('Position',[sbplt4.Position(1)+0.08 sbplt4.Position(2)+0.1 0.23 0.20]);
    vCount(1) = 0;
    bar(vValues,vCount,'BarWidth',1.0)
    grid on;
    set(ax,'Xlim',[ax.XLim(1)+1 ax.XLim(2)]);
    set(ax,'XTick',[1 ax.XTick]);
    
    figurename = 'IceCoverDurationHistograms';
    section = 'discussion';
    export_fig(['/Users/esben/OneDrive - Aarhus Universitet/Speciale/skriv/latex/figures/' section '/' figurename],'-jpg','-r300');
    
end
%distFig
