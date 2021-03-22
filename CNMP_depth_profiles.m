%CNMP_depth_profiles
%esben v1
%plot production rate of cosmogenic nuclides from fast muons and negative
%muon capture
clear; close all;

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

figure; hold on; grid on;
plot(P10m,d,'-r');
plot(P26m,d,'-b');

% format bank
% format compact
% disp('old')
% Lspal
% P10Lfm
% P10Lnmc
% P26Lfm
% P26Lnmc


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
%P10spald = P10spal*exp(-rho*100*d/Lspal);
%P26spald = P26spal*exp(-rho*100*d/Lspal);

plot(P10m,d,'--r');
plot(P26m,d,'--b');

plot(P10spald,d,'-.r');
plot(P26spald,d,'-.b');

xlim([10e-8 10e3])
%ylim([0 500]);
set(gca,'Ydir','reverse')
set(gca, 'XScale', 'log')
legend('old Be10 P','old Al26 P','new Be10 P','new Al26 P','Be10 spallation','Al26 spallation');
ylabel('depth [m]')
xlabel('production rate [atom/g/yr]')