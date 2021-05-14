%crossplot_Tdg_D18OT

%make sure to load a model by using the GUI in matlab or something

%load('/Users/esben/OneDrive - Aarhus Universitet/Speciale/Inversionsmodellering/varDP/models/TESTG2_gausta_19_Feb_2021_16_47_47/gausta_data_2.mat')

clear; close all;

load models/2/gausta_data_2.mat



%figure;
%hold on;
%grid on;
%xlabel('d18O threshold'); ylabel('Tdg');


dint = linspace(0,5.3,100);
Tint = linspace(0.0093,0.0248,100);
[dbin,Tbin] = meshgrid(dint,Tint);
histgrid = zeros(size(Tbin));


for nw = 1:length(model.walker)
    I = find(model.walker{nw}.status == 1);
    d18O_t = model.walker{nw}.u(1,I);
    Tdg = model.walker{nw}.u(2,I);
    
    %scatter(d18O_t,Tdg,0.25)
    
    [~,Id] = min(abs(d18O_t-dint'));
    d18O_t_new = dint(Id);
    ddisp = [d18O_t;d18O_t_new];
    
    [~,It] = min(abs(Tdg-Tint'));
    Tdg_new = Tint(It);
    tdisp = [Tdg;Tdg_new];
    
    
    for i = 1:length(I)
        histgrid(It(i),Id(i)) = histgrid(It(i),Id(i)) + 1;
    end
end


histnorm = (histgrid/max(max(histgrid)));
logbins = logspace(-2,0,100);


figure;
hold on;
xlabel('d18O threshold'); ylabel('Tdg');

contourf(dint,Tint,histnorm,logbins,'linestyle','none')
colorbar('Ticks',[0.001,0.01,0.1,1],'TickLabels',{'0.001','0.01','0.1','1'})
set(gca,'ColorScale','log')
set(gcf,'Color','White')


%map = colormap(parula(1000));
%map(1,:) = [1,1,1];
%colormap(map);
%colorbar

%xlim([3;5.3])
%set(gca,'ColorScale','log')