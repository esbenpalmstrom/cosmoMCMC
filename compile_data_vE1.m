close all; clear;
%{
modified from compile_data_vJ2.m by JLA

Should be able to include any number of real depth points and nuclides and/or a number
of synthetic depth points and nuclides.

WIP notes:

%}

% *****settings******
nuclideSetting = 'BeAl'; %'BeAl' or 'Al'
depthSetting = [0;0.5;1]; %add depths, ex: [0;0.5;1], use ; for separation
savefolder = strcat('data/FS');



depthstr = regexprep(num2str(depthSetting'), ' +', '_');
addpath Functions
addpath data
addpath data/FS
load fsamples.mat %focus samples made in compile_samplesv2.m
ns = length(fieldnames(fsamples))-1; %number of samples, -1 since ID cell is not a sample.
CNprop = getCNprop;

% *****compile real sample data******
for j = 1:ns
    
    excelfile = strcat('data/FS/',fsamples.IDs{j},'.xlsx');
    [num,text,~] = xlsread(excelfile);
    nd = size(text,1)-1; %number of samples from same location
    ij = 1;
    %compile data from xlsx into sample cell
    for i = 1:nd
        sample{i}.name = text{ij+1,1};
        sample{i}.type = text{ij+1,2};
        sample{i}.site = text{ij+1,3};
        sample{i}.Ndp = num(ij,1);
        sample{i}.depths = num(ij:ij+sample{i}.Ndp-1,2);
        sample{i}.Nnc = num(ij,3);
        
        sample{i}.N10 = num(ij:ij+sample{i}.Ndp-1,4); %N10
        sample{i}.dN10 = num(ij:ij+sample{i}.Ndp-1,5); %err10
        sample{i}.N26 = num(ij:ij+sample{i}.Ndp-1,6); %N26
        sample{i}.dN26 = num(ij:ij+sample{i}.Ndp-1,7); %err26
        
        sample{i}.P10total = num(ij,8); %Surface N10 production rate
        sample{i}.P10spal = num(ij,12);
        sample{i}.P10muon = num(ij,13);
        sample{i}.P26total = num(ij,9); %Surface N26 production rate
        sample{i}.P26spal = num(ij,14);
        sample{i}.P26muon = num(ij,15);
        
        sample{i}.elevation = num(ij,10);
        sample{i}.T10 = num(ij,11);
        sample{i}.r2610 = sample{i}.N26./sample{i}.N10;
        sample{i}.dr2610 = sample{i}.r2610.*sqrt((sample{i}.dN10./sample{i}.N10).^2+(sample{i}.dN26./sample{i}.N26).^2);
        sample{i}.Lspal = CNprop.Lspal;
        sample{i}.rho = CNprop.rho;
        
        %model parameters
        model.Nsnr = 1; %number of samples
        model.Nfree = 2; %Number of free depth points
        model.Nsmp = 2*model.Nfree + 2; %number of sample specific parameters
        
        %model.data{i}.depths=sample{i}.depths; %depths of datapoints (m)
        model.data{i}.depths = depthSetting;
        model.Ndp = length(model.data{i}.depths); %Number of data points in depth profile
        model.Nnc = num(ij,3); %number of nuclides
        model.Nds = model.Nnc*model.Ndp; %number of data per sample (nuclides*depths)
        model.age = 3.0; %Max time (Myr)
        model.z0 = 20; %Max depth
        model.Temp = 1.0;
        model.Mmp = 2; %number of generic model parameters
        
        model.data{i}.P10total = sample{i}.P10total;
        model.data{i}.P10muon = sample{i}.P10muon;
        model.data{i}.P10spal = sample{i}.P10spal;
        
        model.data{i}.P26total = sample{i}.P26total;
        model.data{i}.P26muon = sample{i}.P26muon;
        model.data{i}.P26spal = sample{i}.P26spal;
        
        
        sname = fsamples.IDs{j};
        d18Op = fsamples.(sname).d18Ot;  % d18O threshold (3.5-5)
        Tdgla = fsamples.(sname).Tdgla; % Time of last deglaciation [Myr] (2-15e-3)
        Z1    = fsamples.(sname).Z1;    % Depth at time of deglaciation [m] (0-0.05)
        dT2   = fsamples.(sname).dT2;  % Time change from Tdgla [Myr] (0-2 Myr)
        dZ2   = fsamples.(sname).dZ2;   % Depth change during T2 [log10(m)] (-1 to 1)
        dT3   = fsamples.(sname).dT3;  % Time change from Tdgla+dT2 [Myr] (0-2 Myr)
        dZ3   = fsamples.(sname).dZ3;  % Depth change during T3 [log10(m)] (-1 to 1)
        E0    = fsamples.(sname).E0;    %Erosion rate exponent, %10^(-2) = 0.01 m/Myr [log10(m/Myr)] (-2 to 2)
        up = [d18Op,Tdgla,Z1,dT2,dZ2,dT3,dZ3,E0]; %Pack model parameters to vector structure
        zp=[0.0,Z1,Z1+10^(dZ2),Z1+10^(dZ2)+10^(dZ3),Z1+10^(dZ2)+10^(dZ3)+(3 - Tdgla+dT2+dT3)*10^(E0)];
        tp=[0.0,Tdgla*10^3,(Tdgla+dT2)*10^3,(Tdgla+dT2+dT3)*10^3,3000];
        model.synpmval=up;
        
        
        switch nuclideSetting
            case 'BeAl'
                if model.Nnc == 1
                    model.Nnc = 2;
                    model.Nds = model.Nnc*model.Ndp;
                    model.synthAl = 1;
                elseif model.Nnc == 2
                    model.synthAl = 0;
                end
                
                [gm,~,~] = forward_bedrockvJ2(up,model,CNprop);
                
                sample{i}.depths = depthSetting;
                
                %if length(depthSetting) >= 2
                N10_temp = sample{i}.N10;
                N10unc_temp = sample{i}.dN10;
                N26_temp = sample{i}.N26;
                N26unc_temp = sample{i}.dN26;
                
                sample{i}.N10(1:(length(gm)/2)) = gm(1:length(gm)/2);
                sample{i}.N26(1:(length(gm)/2)) = gm((length(gm)/2)+1:end);
                
                N10unc = getCosmoUnc(sample{i}.N10,10,'Be10',0); %calculates likely error in % on 10Be based on 'historical data' in Aarhus lab + AMS
                err10 = sample{i}.N10.*N10unc/100; %N10spal.*N10unc/100; Corresponding error in at/g on synthetic samples
                sample{i}.dN10 = err10;
                
                N26unc = getCosmoUnc(sample{i}.N26,10,'Al26',0); %calculates likely error in % on 10Be based on 'historical data' in Aarhus lab + AMS
                err26 = sample{i}.N26.*N26unc/100; %N10spal.*N10unc/100; Corresponding error in at/g on synthetic samples
                sample{i}.dN26 = err26;
                
                
                %sample{i}.N10(1) = N10_temp; %uncomment to replace surface conc. with real data.
                %sample{i}.dN10(1) = N10unc_temp;
                
                if model.synthAl ~= 1
                    %sample{i}.N26(1) = N26_temp; %uncomment to replace with real data
                    %sample{i}.d26(1) = N26unc_temp;
                    
                end
                
                %end
                
                sample{i}.r2610 = sample{i}.N26./sample{i}.N10;
                sample{i}.dr2610 = sample{i}.r2610.*sqrt((sample{i}.dN10./sample{i}.N10).^2+(sample{i}.dN26./sample{i}.N26).^2);
                
            case 'Be'
                model.Nnc = 1;
                model.Nds = model.Nnc*model.Ndp;
                [gm,~,~] = forward_bedrockvJ2(up,model,CNprop);
                
                sample{i}.depths = depthSetting;
                
                if length(depthSetting) >= 2
                    N10_temp = sample{i}.N10;
                    N10unc_temp = sample{i}.dN10;
                    
                    sample{i}.N10(1:(length(gm)/2)) = gm(1:length(gm)/2);
                    N10unc = getCosmoUnc(sample{i}.N10,10,'Be10',0); %calculates likely error in % on 10Be based on 'historical data' in Aarhus lab + AMS
                    err10 = sample{i}.N10.*N10unc/100; %N10spal.*N10unc/100; Corresponding error in at/g on synthetic samples
                    sample{i}.dN10 = err10;
                    
                    %sample{i}.N10(1) = N10_temp; %uncomment to replace
                    %with real data
                    %sample{i}.dN10(1) = N10unc_temp; % 
                    
                end
                sample{i}.N26(1:(length(gm)/2)) = NaN;
                sample{i}.dN26(1:(length(gm)/2)) = NaN;
                sample{i}.r2610(1:(length(gm)/2)) = NaN;
                sample{i}.dr2610(1:(length(gm)/2)) = NaN;
                
        end
        
        ij=ij+sample{i}.Ndp; %Look for next sample in this row
    end
    
    foldername = strcat(nuclideSetting,depthstr);
    mkdir(savefolder,foldername);
    savefile = strcat(savefolder,'/',foldername,'/',fsamples.IDs{j},'.mat');
    save(savefile,'sample','model');
    clear sample
    clear model
end