function compile_GaustaData_vE2()

%script to read and compile data specifically for gaustatoppen from a
%specific excel file.
%using new muon production method

%{
WIP:
Reduce total production by the difference in the old and new muon
production rates.
%}

addpath('MuonP_Jane/Functions','MuonP_Jane/Functions/cl36','MuonP_Jane/Functions/cl36/scaling','MuonP_Jane/Functions/CronusCalc')


close all;

%read excel file
ns = 1; %number of samples
noc = 2; %number of nuclides, if 1, it is assumed that the nuclide is Be

excelfile = 'data/Gausta/GaustaData_2.xlsx';


[num,text,~] = xlsread(excelfile);
ij=1; %row start
for i=1:ns
    
    
    sample{i}.Ndp = size(num,1); %Number of data points at depth for sample
    sample{i}.name = 'gausta';
    sample{i}.type = 'bedrock';
    sample{i}.batchid = text{ij+1,1};
    %Read depth specific data
    sample{i}.depths = (num(ij:ij+sample{i}.Ndp-1,7))/100; %sample depths
    sample{i}.N10 = num(ij:ij+sample{i}.Ndp-1,1); %N10
    sample{i}.dN10 = num(ij:ij+sample{i}.Ndp-1,3); %err10
    sample{i}.N26 = num(ij:ij+sample{i}.Ndp-1,4); %N26
    sample{i}.dN26 = num(ij:ij+sample{i}.Ndp-1,6); %err26
    sample{i}.elev = 1715; %elevation (m)
    
    sample{i}.lat = 59.8482;
    sample{i}.lon = 8.6600;
    sample{i}.thick = 1; %what are the thickness for the gausta samples?
    
    sample{i}.r2610 = sample{i}.N26./sample{i}.N10; %calculate N26/N10 ratio
    sample{i}.dr2610 = sample{i}.r2610.*sqrt((sample{i}.dN10./sample{i}.N10).^2+(sample{i}.dN26./sample{i}.N26).^2);
    sample{i}.density = 2.65;
    sample{i}.T10 = 74000; %apparent 10Be age (yrs), use cronus earth
    
    
    rho = sample{i}.density; %Density of rock sample (g/cm3)
    
    %Define depths below surface z/rho cm/(g/cm3) [g/cm^2] for fitting of production profiles
    
    
    D_m = 100; %Depth, changed below (rho)
    %D_m = 352.45;
    
    z_m = linspace(0,10,100);
    z_D = D_m*z_m.^3/10*rho; %denser depth-grid near surface
    %Convert elevation to atm pressure (hPa)
    p = ERA40atm(sample{i}.lat,sample{i}.lon,sample{i}.elev);
    sample{i}.pressure = p;
    %Spallation attenuation and thickness correction
    Lspal=attenuationlength(sample{i}.lat,sample{i}.lon,sample{i}.elev,p); %Calculated from CronusCalc functions based on site cutoff rigidity, not considering terrain shielding %[g/cm2]
    sample{i}.production.Lspal=Lspal; %[g/cm2]
    sf_spal = exp(-sample{i}.thick/2*rho/Lspal); %Factor to correct production rate for thickness of sample, sets surface production =production midway through sample. Make sure this is not already factored in to site-specific production rates
    
    maxZ = 1200; %maxdepth (g/cm2) for muon-production profile used for fitting
    % of exponentials below. If this depth is very large, exponential terms will
    % be dominated by fast muon production, which isn't ideal. 1200g/cm2=4.5m
    % with rho ~2.65-2.7, Test effect of this choice
    
    
    %Be10 production
    %sample{i}.P10total = 19.3360; %Surface Be10 production rate
    %sample{i}.production.P10total = sample{i}.P10total;
    %sample{i}.production.P10spal = 0.98*sample{i}.production.P10total; %later combined, then split with other factors
    sample{i}.production.P10spal = 18.88;
    %muons
    sample{i}.P10spal = sample{i}.production.P10spal;
    mc10.k_neg = 0.00191 .* 0.704 .* 0.1828; % From BCO fit, Model 1A, alpha=1; f_star*f_C*f_D
    mc10.sigma0 = 0.280e-30; % From BCO fit, model 1A, alpha=1;
    mc10.Natoms = 2.006e22; %Oxygen atoms pr gram Quartz
    p10_muons = fit_Pmu_with_exp_1A(p,min(z_D),maxZ,2,mc10,0);
    sample{i}.production.P10_Lnmc=p10_muons.L(2);
    sample{i}.production.P10_Lfm=p10_muons.L(1);
    % out = P_mu_total_alpha1(z_D,p,mc10,'yes'); Split into fast/negative muons
    % P10_top_nm_= out.P_neg(1); % P10_top_fm_= out.P_fast(1);
    shield_fac10_nm = exp(-sample{i}.thick/2*rho/p10_muons.L(2));
    shield_fac10_fm = exp(-sample{i}.thick/2*rho/p10_muons.L(1));
    sample{i}.production.P10_nmc = p10_muons.P(2)*shield_fac10_nm;
    sample{i}.production.P10_fm = p10_muons.P(1)*shield_fac10_fm;
    sample{i}.P10muon = sample{i}.production.P10_nmc + sample{i}.production.P10_fm;
    
    sample{i}.P10total = sample{i}.production.P10spal + sample{i}.P10muon; %Surface Be10 production rate
    sample{i}.production.P10total = sample{i}.P10total;
    
    %Al26 production
%     sample{i}.P26total = 132.2438; %Surface Al26 production rate
%     sample{i}.production.P26total = sample{i}.P26total;
%     sample{i}.production.P26spal = 0.98*sample{i}.production.P26total; %later combined, then split with other factors
    sample{i}.production.P26spal = 127.39;
    %muons
    sample{i}.P26spal = sample{i}.production.P26spal;
    mc26.k_neg = 0.0133 .* 0.296 .* 0.6559; % From BCO fit, Model 1A, alpha=1; f_star*f_C*f_D
    mc26.sigma0 = 3.89e-30; % From BCO fit, model 1A, alpha=1;
    mc26.Natoms = 1.003e22; %Si atoms pr gram Quartz
    p26_muons = fit_Pmu_with_exp_1A(p,min(z_D),maxZ,2,mc26,0);
    sample{i}.production.P26_Lnmc=p26_muons.L(2);
    sample{i}.production.P26_Lfm=p26_muons.L(1);
    shield_fac26_nm = exp(-sample{i}.thick/2*rho/p26_muons.L(2));
    shield_fac26_fm = exp(-sample{i}.thick/2*rho/p26_muons.L(1));
    sample{i}.production.P26_nmc = p26_muons.P(2)*shield_fac26_nm;
    sample{i}.production.P26_fm = p26_muons.P(1)*shield_fac26_fm;
    sample{i}.P26muon = sample{i}.production.P26_nmc + sample{i}.production.P26_fm;
    
    sample{i}.P26total = sample{i}.production.P26spal + sample{i}.P26muon; %Surface Al26 production rate
    sample{i}.production.P26total = sample{i}.P26total;
    
    %%% set model parameters
    model.Nsnr = 1; %number of samples
    model.Nfree = 3; %Number of free depth points, changed from 2 to 3
    model.Nsmp = 2*model.Nfree + 2; %number of sample specific parameters
    model.data{i}=sample{i}; %put sample into model
    model.Ndp = length(model.data{i}.depths); %Number of data points in depth profile
    model.Nnc = noc; %number of nuclides
    model.Nds = model.Nnc*model.Ndp; %number of data per sample (nuclides*depths)
    model.age = 3.0; %Max time (Myr)
    model.z0 = 20; %Max depth
    model.Temp = 1.0; %adjusted in run file.
    model.Mmp = 2; %number of generic model parameters
    
    
    
    
    ij=ij+sample{i}.Ndp; %Look for next sample in this row
end

save ./data/gausta_v2/gausta_data_2.mat sample model excelfile