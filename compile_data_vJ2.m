function compile_data_vJ2()

% vJ2 updated Feb. 2020 by JLA: includes the possibility of extra data points at depth

%{
esben changelog
sept. 27: compiled data from stroeven_test.xlsx
october: added lines to compile data from focus_samples.xlsx
%}

close all;

%Read data info from Excelfile
ns = 1; %number of samples
noc = 2; %number of nuclides, if 1, it is assumed that the nuclide is Be

%excelfile = 'data/dptest.xlsx';
%excelfile = 'data/stroeven_test.xlsx'; 
%excelfile = 'data/focus_samplesBe.xlsx';
%excelfile = 'data/focus_samplesBeAl.xlsx';
excelfile = 'data/FS_test/stroe2002a_1.xlsx';


[num,text,~] = xlsread(excelfile);
ij=1; %row start
for i=1:ns
    sample{i}.Ndp = num(ij,1); %Number of data points at depth for sample i
    
    sample{i}.name = text{ij+1,1};
    sample{i}.type = text{ij+1,2};
    sample{i}.batchid = text{ij+1,3};
    %Read depth specific data
    sample{i}.depths = num(ij:ij+sample{i}.Ndp-1,2); %sample depths
    sample{i}.N10 = num(ij:ij+sample{i}.Ndp-1,3); %N10
    sample{i}.dN10 = num(ij:ij+sample{i}.Ndp-1,4); %err10
    sample{i}.N26 = num(ij:ij+sample{i}.Ndp-1,5); %N26
    sample{i}.dN26 = num(ij:ij+sample{i}.Ndp-1,6); %err26
     
    %Read sample specific data
    sample{i}.P10total = num(ij,7); %Surface N10 production rate
    sample{i}.P10spal = 0.98*sample{i}.P10total;
    sample{i}.P10muon = 0.02*sample{i}.P10total;
    
    sample{i}.P26total = num(ij,8); %Surface N26 production rate
    sample{i}.P26spal = 0.98*sample{i}.P26total;
    sample{i}.P26muon = 0.02*sample{i}.P26total;
    
    sample{i}.elevation = num(ij,9); %elevation (m)
    sample{i}.T10 = num(ij,10); %apparent 10Be age (yrs)
    
    sample{i}.r2610 = sample{i}.N26./sample{i}.N10; %calculate N26/N10 ratio
    sample{i}.dr2610 = sample{i}.r2610.*sqrt((sample{i}.dN10./sample{i}.N10).^2+(sample{i}.dN26./sample{i}.N26).^2);
    
    
    %%% set model parameters
    model.Nsnr = 1; %number of samples
    model.Nfree = 2; %Number of free depth points %2?
    model.Nsmp = 2*model.Nfree + 2; %number of sample specific parameters
    model.data{i}.depths=sample{i}.depths; %depths of datapoints (m)
    model.Ndp = length(model.data{i}.depths); %Number of data points in depth profile
    model.Nnc = noc; %number of nuclides
    model.Nds = model.Nnc*model.Ndp; %number of data per sample (nuclides*depths)
    model.age = 3.0; %Max time (Myr)
    model.z0 = 20; %Max depth
    model.Temp = 1.0;
    model.Mmp = 2; %number of generic model parameters
    P10total = sample{i}.P10total;
    model.data{i}.P10muon = 0.02*P10total;
    model.data{i}.P10spal = 0.98*P10total;
    
    
    P26total = sample{i}.P26total;
    model.data{i}.P26muon = 0.02*P26total;
    model.data{i}.P26spal = 0.98*P26total;
    
    ij=ij+sample{i}.Ndp; %Look for next sample in this row
end


%save ./data/dptest.mat sample model excelfile %call & save CNprop here?
%save ./data/stroeven_test.mat sample model excelfile
%save ./data/focus_samplesBe.mat sample model excelfile
%save ./data/focus_samplesBeAl.mat sample model excelfile
save ./data/FS_test/stroe2002a_1_test.mat sample model excelfile