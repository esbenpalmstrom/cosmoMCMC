%run inversions of real sample data.
close all;
load fsamples.mat

%proposed standard setting: 4 walkers, 40k accepted models, 4k accepted models burnin.
%(maybe max 500k proposed models).
snr = 1; %number of samples from same excel file.
numdp = 4; %including glaciation time and depth, right now only made for 3 or 4
CNprod = 'new'; %'new' or 'old', remember to correct the data file used.

nwalkers = 16; %number of walkers
nmodmax = 1000e4; %max number of proposed models
nmodacc = 80e3; %max number of accepted model
lburnin = 10e3; %number of accepted models before burnin is finished
k_meanlength = 100; %default: 100



%for manual Tdg
Tdglac = 0.015; %Boundaries are quite wide for Gausta data, see the bedrock file.
%From 28. march, Tglac is phased out. Set Tdg boundaries manually in
%bedrock file.

acctarget = 0.30; %default: 0.3

model_d18O = 1; %Use d18O threshold from peltier models or not. 1 will use, 0 will not.

starttime = datetime('now');
starttimer = now;
starttimestr = datestr(now,0);
starttimestr = strrep(starttimestr,'-','_');
starttimestr = strrep(starttimestr,' ','_');
starttimestr = strrep(starttimestr,':','_');

%*****base inversion loop******
%correct length if working with all fsamples, or just gaustadata
% for i = 1:length(fsamples.IDs)


foldername = '1';
while exist(['models/' foldername], 'dir') == 7
    foldernumber = str2num(foldername);
    foldername = num2str(foldernumber+1);
end
    
for i = 1:1 %length(fsamples.IDs)
    
    sampleID = 'gausta_data_2';
    %samplepath = ['./data/Gausta/' sampleID '.mat'];
    samplepath = ['./data/gausta_v2/' sampleID '.mat']; %Gausta or gausta_v2
    
    %sampleID = fsamples.IDs{i};
    %samplepath = ['./data/FS/' sampleID '.mat'];
    
    savepath = ['models/' foldername '/' sampleID '.mat'];
    reportpath = ['models/' foldername '/reports/' sampleID '/'];
    mkdir (reportpath);
    bedrockMCvJ2_E2(snr,nwalkers,nmodmax,nmodacc,samplepath,lburnin,savepath,sampleID,Tdglac,k_meanlength,acctarget,numdp,CNprod,model_d18O) %use vE1 or E2
    makereportEv2(snr,savepath,reportpath,sampleID,numdp)
end


datetime('now');
finishtime = datestr(now);
endtimer = now;
timer = datestr(endtimer-starttimer,'HH:MM:SS.FFF');

C = {'nr. of walkers',nwalkers;
    'nr. of accepted models pr. walker',nmodacc;
    'length of burnin',lburnin;
    'start time',starttimestr;
    'finish time',finishtime;
    'time elapsed',timer;
    'Notes:',''
    };

writecell(C,['models/' foldername '/log.txt'],'Delimiter',';');