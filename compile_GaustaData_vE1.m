function compile_GaustaData_vE1()
%script to read and compile data specifically for gaustatoppen from a
%specific excel file.
%using old muon production method

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
     
    %Read sample specific data
    sample{i}.P10total = 19.3360; %Surface N10 production rate
    sample{i}.P10spal = 0.98*sample{i}.P10total;
    sample{i}.P10muon = 0.02*sample{i}.P10total;
    
    sample{i}.P26total = 132.2438; %Surface N26 production rate
    sample{i}.P26spal = 0.98*sample{i}.P26total;
    sample{i}.P26muon = 0.02*sample{i}.P26total;
    
    sample{i}.elevation = 1715; %elevation (m)
    sample{i}.T10 = 74000; %apparent 10Be age (yrs), use cronus earth
    
    sample{i}.r2610 = sample{i}.N26./sample{i}.N10; %calculate N26/N10 ratio
    sample{i}.dr2610 = sample{i}.r2610.*sqrt((sample{i}.dN10./sample{i}.N10).^2+(sample{i}.dN26./sample{i}.N26).^2);
    
    
    %%% set model parameters
    model.Nsnr = 1; %number of samples
    model.Nfree = 3; %Number of free depth points, changed from 2 to 3
    model.Nsmp = 2*model.Nfree + 2; %number of sample specific parameters
    model.data{i}.depths=sample{i}.depths; %depths of datapoints (m)
    model.Ndp = length(model.data{i}.depths); %Number of data points in depth profile
    model.Nnc = noc; %number of nuclides
    model.Nds = model.Nnc*model.Ndp; %number of data per sample (nuclides*depths)
    model.age = 3.0; %Max time (Myr)
    model.z0 = 20; %Max depth
    model.Temp = 1.0; %adjusted in run file.
    model.Mmp = 2; %number of generic model parameters
    P10total = sample{i}.P10total;
    model.data{i}.P10muon = 0.02*P10total;
    model.data{i}.P10spal = 0.98*P10total;
    
    
    P26total = sample{i}.P26total;
    model.data{i}.P26muon = 0.02*P26total;
    model.data{i}.P26spal = 0.98*P26total;
    
    ij=ij+sample{i}.Ndp; %Look for next sample in this row
end


save ./data/Gausta/gausta_data_2.mat sample model excelfile