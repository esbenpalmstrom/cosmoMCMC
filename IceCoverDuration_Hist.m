%IceCoverDuration_Hist.m

%Todo: Normalize, so it shows number of models with ice compared to total
%number of models.

clear; close all;
load('data/d18Ocurves.mat');
%choose which models to load:
models = [4];

figurename = {'IceCoverDurHist004_a_fraction'
    'IceCoverDurHist004_b_total'};
titles = {''};

set(groot','defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');


% figurename = {'IceCoverDurHist002'
%     'IceCoverDurHist004'
%     'IceCoverDurHist005'
%     'IceCoverDurHist006'
%     'IceCoverDurHist016'
%     'IceCoverDurHist017'};



% 
% titles = {'OG, all data'
%     'Revised, all data'
%     'Revised, only surface'
%     'OG, only surface'
%     'Revised, Peltier'
%     'OG, Peltier'};

eem = 116;
eemt = zeros(1,eem);
LGM = 26;
LGMt = zeros(1,LGM);

d18O_eem = d18O_10ky(1:eem);


LGM_iceduration = [];
eem_iceduration = [];

for i = 1:length(models)
    close all;
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
    
    %calculate fraction of models having ice since eem:
    %sum(eem_iceduration == 0)/length(eem_iceduration);
    
    fig1 = figure(); 
    hold on;
    sbplt1 = subplot(1,1,1);
    box on;
    %set(gcf,'Position',[100 100 900 750]);
    eemt_norm = (eemt./length(eem_iceduration));
    bar(eemt_norm,'BarWidth',1.0)
    xlabel('Time [kyr]'), ylabel('Fraction of models having ice at given time')
    title('')
    grid on;
    sbplt1.YLim = [0 0.6];
    r = rectangle('Position',[10 0 8 0.07],'LineWidth',1.5);
    
    ax1 = axes('Position',[sbplt1.Position(1)+0.5 sbplt1.Position(2)+0.4 0.23 0.30]);
    ax1.Position = [sbplt1.Position(1)+0.5 sbplt1.Position(2)+0.35 0.23 0.40];
    box on;
    hold on; grid on;
    LGMt_norm = (LGMt./length(LGM_iceduration));
    bar(LGMt_norm,'BarWidth',1.0)
    xlabel('Time [kyr]'), ylabel('Fraction of models having ice at given time')
    title('')
    ax1.YLim = [0 0.07];
    ax1.XLim = [10 18];
    ax1.FontSize = 8;
    set(gcf,'Position',[100 100 1200 375]);
    set(gcf,'color','white')
    export_fig(['/Users/esben/OneDrive - Aarhus Universitet/Speciale/skriv/latex/figures/' figurename{1}],'-jpg','-r300');
    
    
%     sbplt2 = subplot(2,2,2);
%     box on;
%     hold on; grid on;
%     %set(gcf,'Position',[100 100 900 750]);
%     bar(eemt_norm,'BarWidth',1.0)
%     xlabel('Time [kyr]'), ylabel('Likelihood of models to have ice at given time')
%     xlim([0 50]);
%     title('')
%     sbplt2.YLim = [0 0.6];
%     
%     sbplt3 = subplot(2,2,3);
%     box on;
%     hold on; grid on;
%     %set(gcf,'Position',[100 100 900 750]);
%     LGMt_norm = (LGMt./length(LGM_iceduration));
%     bar(LGMt_norm,'BarWidth',1.0)
%     xlabel('Time [kyr]'), ylabel('Likelihood of models to have ice at given time')
%     title('')
%     sbplt3.YLim = [0 0.6];
    
    
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
    
    
    fig2 = figure();
    %ice cover duration freq. for all models
    sbplt4 = subplot(1,1,1);
    box on;
    hold on; grid on;
    set(gcf,'Position',[100 100 1200 375]);
    [vCount, vValues] = hist(eem_iceduration(:),62);
    vCount = vCount./sum(vCount);
    bar(vValues,vCount,'BarWidth',1.0)
    %r = rectangle('Position',[1.2 0 sbplt4.XLim(2) 0.05],'LineWidth',1.5);
    %histogram(eem_iceduration,'Normalization','probability');
    xlabel('Duration of ice cover since the Eemian [kyr]'), ylabel('Fraction of models')
    %title('total ice cover duration frequency')
    set(gcf,'color','white')
    sbplt4.YLim = [0 0.6];
    sbplt4.YScale = 'log';
    sbplt4.YLim = [10e-6 10e-1];
    sbplt4.YTickLabel = {'0.0001'
        '0.001'
        '0.01'
        '0.1'
        '1'};
    
    
    %sumcurve = cumsum(vCount,'reverse');
    %plot(vValues,sumcurve,'-b')
    
    
%     ax2 = axes('Position',[sbplt4.Position(1)+0.56 sbplt4.Position(2)+0.2 0.23 0.20]);
%     vCount(1) = 0;
%     bar(vValues,vCount,'BarWidth',1.0)
%     grid on;
%     set(ax2,'Xlim',[ax2.XLim(1)+1 ax2.XLim(2)]);
%     %set(ax,'XTick',[1 ax.XTick]);
%     ax2.XTick = [1 10 20 30 40 50 60 70 80];
%     ax2.YLim = [0 0.05];
%     ax2.YLabel.String = 'Fraction of models';
%     ax2.XLabel.String = 'Duration of ice cover since the Eemian [kyr]';
%     ax2.FontSize = 8;
    
    sgtitle(titles{i});
    
    sbplt1.FontSize = 14;
    sbplt4.FontSize = 14;
    ax1.FontSize = 9;
    
    
    
    %export_fig(['/Users/esben/OneDrive - Aarhus Universitet/Speciale/skriv/latex/figures/' figurename{2}],'-jpg','-r300');
    
end
%distFig
