%CNMP_depth_profiles
%esben v1
%plot production rate of cosmogenic nuclides from fast muons and negative
%muon capture
clear; close all;

set(groot','defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

addpath Functions

%create depth profile
d = linspace(0,10,201); %[m]
%ds = linspace(0,10,201);

CNprop = getCNprop;

%calculate production rate at each depth point

%old produciton
load data/Gausta/gausta_data_2.mat
P10spal = model.data{1}.P10spal;
P10tot = model.data{1}.P10spal + model.data{1}.P10muon;
P10fm = CNprop.pr_fm_Be*P10tot;
P10nmc = CNprop.pr_nmc_Be*P10tot;

P26spal = model.data{1}.P26spal;
P26tot = model.data{1}.P26spal + model.data{1}.P26muon;
P26fm = CNprop.pr_fm_Al*P26tot;
P26nmc = CNprop.pr_nmc_Al*P26tot;

rho = CNprop.rho;
Lspal = CNprop.Lspal;
P10Lnmc = CNprop.Lnmc;
P26Lnmc = CNprop.Lnmc;
P10Lfm = CNprop.Lfm;
P26Lfm = CNprop.Lfm;

%nmc
P10nmcd = P10nmc*exp(-CNprop.rho*100*d/CNprop.Lnmc);
P26nmcd = P26nmc*exp(-CNprop.rho*100*d/CNprop.Lnmc);
%fm
P10fmd = P10fm*exp(-CNprop.rho*100*d/CNprop.Lfm);
P26fmd = P26fm*exp(-CNprop.rho*100*d/CNprop.Lfm);
%summed
P10m = P10nmcd + P10fmd;
P26m = P26nmcd + P26fmd;

%spallation
P10spald = P10spal*exp((-CNprop.rho*100*d)/CNprop.Lspal);
P26spald = P26spal*exp((-CNprop.rho*100*d)/CNprop.Lspal);


% ************* figure 1 ***************
figure; hold on; grid on; box on;

plot(P10m,d,'--b'); %old muon production rate
P10tot = P10m + P10spald;
plot(P10tot,d,'-b'); %old total production rate
%plot(P26m,d,'-.r');

%new production
load data/gausta_v2/gausta_data_2.mat

rho = model.data{1}.density;
Lspal = model.data{1}.production.Lspal;

P10spal = model.data{1}.production.P10spal;
P10fm = model.data{1}.production.P10_fm;
P10Lfm = model.data{1}.production.P10_Lfm;
P10nmc = model.data{1}.production.P10_nmc;
P10Lnmc = model.data{1}.production.P10_Lnmc;

P26spal = model.data{1}.production.P26spal;
P26fm = model.data{1}.production.P26_fm;
P26Lfm = model.data{1}.production.P26_Lfm;
P26nmc = model.data{1}.production.P26_nmc;
P26Lnmc = model.data{1}.production.P26_Lnmc;

% disp('new')
% Lspal
% P10Lfm
% P10Lnmc
% P26Lfm
% P26Lnmc

%nmc
P10nmcd = P10nmc*exp(-rho*100*d/P10Lnmc);
P26nmcd = P26nmc*exp(-rho*100*d/P26Lnmc);
%fm
P10fmd = P10fm*exp(-rho*100*d/P26Lfm);
P26fmd = P26fm*exp(-rho*100*d/P26Lfm);
%summed
P10m = P10nmcd + P10fmd;
P26m = P26nmcd + P26fmd;

%spal
P10spald = P10spal*exp(-rho*100*d/Lspal);
P26spald = P26spal*exp(-rho*100*d/Lspal);

P10tot = P10spald + P10m;
P26tot = P26spald + P26m;

plot(P10spald,d,'-r','LineWidth',1.2); % spallation CNP
%plot(P26spald,d,':r');

plot(P10m,d,'--k','LineWidth',1.2); % new muon CNP
%plot(P26m,d,'--r');

plot(P10tot,d,'-k','LineWidth',1.2) % new total CNP
%plot(P26tot,d,'-r')



xlim([10e-4 10e2])
%ylim([0 500]);
set(gca,'Ydir','reverse')
set(gca, 'XScale', 'log')
legend('high muogenic Beryllium-10 production',...
'high total Beryllium-10 production',...
'spallation Be10 production',...
'low muogenic Be10 production','low total Beryllium-10 production','Location','northwest')
%legend('old Be10 P','old Al26 P','new Be10 P','new Al26 P','Be10 spallation','Al26 spallation','total Be10 production','total Al26 production');
ylabel('depth [m]')
xlabel('production rate [atoms/g/yr]')
xlim([10e-4 10e1])
set(gca,'FontSize',16);
set(gcf,'color','white')
set(gcf,'Position',[100 100 1200 600]);
figurename = 'CNP_newvsold';
%export_fig(['/Users/esben/OneDrive - Aarhus Universitet/Speciale/skriv/latex/figures/' figurename],'-jpg','-r300');

% ************* figure 2 ***************
figure; hold on; grid on; box on;
plot(P10tot,d,'k','LineWidth',1.2)
plot(P26tot,d,'g','LineWidth',1.2)
xlim([10e-4 10e2])
%ylim([0 500]);
set(gca,'Ydir','reverse')
set(gca, 'XScale', 'log')
ylabel('depth below surface [m]')
xlabel('production rate [atoms/g/yr]')
legend('Total Beryllium-10 production','Total Aluminium-26 production','Location','northwest')

set(gca,'FontSize',16);
set(gcf,'color','white')
set(gcf,'Position',[100 100 800 700]);
figurename = 'CNP_tot';
%export_fig(['/Users/esben/OneDrive - Aarhus Universitet/Speciale/skriv/latex/figures/' figurename],'-jpg','-r300');


% ************* figure 3 ***************
figure; hold on; grid on; box on;
plot(P10tot,d,'-k','LineWidth',1.2)
plot(P10spald,d,':k','LineWidth',1.2)
plot(P10m,d,'--k','LineWidth',1.2)

set(gca,'Ydir','reverse')
set(gca, 'XScale', 'log')
ylabel('depth below surface [m]')
xlabel('production rate [atoms/g/yr]')
legend('Total Beryllium-10 production','Beryllium-10 spallation production','Beryllium-10 muon production','Location','northwest')
xlim([10e-4 10e1])
ylim([0 10])
set(gca,'FontSize',16);
set(gcf,'color','white')
set(gcf,'Position',[100 100 1200 600]);

%figurename = 'defense_CNP_spalvsmuon';
%export_fig(['/Users/esben/OneDrive - Aarhus Universitet/Speciale/skriv/latex/figures/' figurename],'-jpg','-r300');


% % ******** figure 4, for defense presentation ********
% figure; hold on;
% plot(P10spald,d,'b','LineWidth',1.2)
% plot(P10m,d,'r','LineWidth',1.2)
% set(gca,'Ydir','reverse')
% ylim([0 2])
% set(gcf,'color','white')
% set(gcf,'Position',[100 100 500 700]);
% set(gca,'XTick',[]);
% set(gca,'YTick',[]);
% legend('Spallation','Muogenic','Location','southeast')
% set(gca,'FontSize',16);
