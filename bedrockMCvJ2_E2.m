function bedrockMCvJ2_E2(snr,nwalkers,nmodmax,nmodacc,samplepath,lburnin,savepath,sampleID,Tdglac,k_meanlength,acctarget,numdp,CNprod)
%MCMC model modified from bedrockMCv7 (DLE) by JLA, June 2019
%
% vJ2 updated Feb. 2020 by JLA includes the possibility of extra data points at depth

%{
esben changelog

sept. 27: ran data from stroeven_test.mat (from stroeven_test.xlsx)
Only changed data loading and saving.

early october: ran with focus_samples.

late october: changed 'obs' and 'sigd' calcuation to be dependent on number
of nuclides. 1 nuclide assumes Be.

added input parameters in order to make compatible with run_inversions.m


todo:
Make a switch for choosing between old and new muon production scheme. This
should happen in the forward model function, forward_bedrock_vE1.m

parallelize the workflow



%}

load fsamples.mat
addpath Functions data models export_fig;
addpath data/FS

%graphics switch
gg = 0;
close all;

if gg
    set(gcf,'units','normalized','position',[.1,.2,.4,.7]);
    set(gca,'position',[0,0,1,1],'visible','off');
    set(gca,'xlim',[0,1],'ylim',[0,1]);
    ax0 = gca;
    col1 = 0.9*[1,1,1];
    col2 = 0.3*[1,1,1];
    col3 = 0.6*[1,1,1];
end


%Cosmogenic properties
CNprop = getCNprop;


%**** load data *****
%load ./data/synthetic_data_fw_dp_Temp.mat 'sample' 'model';

load(samplepath,'sample','model')

%initialize
%models = struct();b

%below is adjusted in the data compilation. Consider changing?
% model.age = 3.0; %Max time (Myr)
% model.z0 = 20; %Max depth
temp = 1.0;
model.Temp = temp; %standard: 1.0

model.tc = 1.0; %time constraint for each data point except Tdg, Myr.

% model.Mmp = 2; %number of generic model parameters
model.mp{1}.name = ['d18O threshold'];
model.mp{1}.vmin = 3.5;
model.mp{1}.vmax = 5.3; %changed from 5.0 to 5.3


model.mp{2}.name = ['Time of deglaciation (Tdg)']; % use uncertainty of +/- 1ka
if isfield(fsamples,sampleID)
    
    model.mp{2}.vmin = fsamples.(sampleID).Tdgla - 1e-3; %Myr
    model.mp{2}.vmax = fsamples.(sampleID).Tdgla + 1e-3;
else
    
    %If not analyzing a sample in fsamples, use a manually inserted Tdglac
    
    %for example the gaustatoppen deglaciation.
    model.mp{2}.vmin = Tdglac - 1e-3;
    model.mp{2}.vmax = Tdglac + 14.5e-3;
    %model.mp{2}.vmax = Tdglac + 34.5e-3; %test
end

%loop number of samples and for model parameters
for i=1:model.Nsnr
    
    switch numdp
        case 3
            model.mp{model.Mmp+(i-1)*model.Nsmp+1}.name =  ['Z at Tdg, sample ',num2str(i)];
            model.mp{model.Mmp+(i-1)*model.Nsmp+1}.vmin = 0;
            model.mp{model.Mmp+(i-1)*model.Nsmp+1}.vmax = 0.05; %m
            model.mp{model.Mmp+(i-1)*model.Nsmp+2}.name =  ['dT2, sample ',num2str(i)];
            model.mp{model.Mmp+(i-1)*model.Nsmp+2}.vmin = 0;
            model.mp{model.Mmp+(i-1)*model.Nsmp+2}.vmax = model.tc; %Myr
            model.mp{model.Mmp+(i-1)*model.Nsmp+3}.name =  ['dz2, sample ',num2str(i)];
            model.mp{model.Mmp+(i-1)*model.Nsmp+3}.vmin = -1; %0.1 m
            model.mp{model.Mmp+(i-1)*model.Nsmp+3}.vmax = 1.0; %1; %log 10 m
            model.mp{model.Mmp+(i-1)*model.Nsmp+4}.name =  ['dT3, sample ',num2str(i)];
            model.mp{model.Mmp+(i-1)*model.Nsmp+4}.vmin = 0;
            model.mp{model.Mmp+(i-1)*model.Nsmp+4}.vmax = model.tc; %Myr
            model.mp{model.Mmp+(i-1)*model.Nsmp+5}.name =  ['dz3, sample ',num2str(i)];
            model.mp{model.Mmp+(i-1)*model.Nsmp+5}.vmin = -1;%0.1 m
            model.mp{model.Mmp+(i-1)*model.Nsmp+5}.vmax = 1.0;%1;%log 10 m
            model.mp{model.Mmp+(i-1)*model.Nsmp+6}.name =  ['E4, sample ',num2str(i)];
            model.mp{model.Mmp+(i-1)*model.Nsmp+6}.vmin = -2; %log10 -2 =0.01m/Myr
            model.mp{model.Mmp+(i-1)*model.Nsmp+6}.vmax = 2; %2=100 m/Myr
            
            %add sample info to models
            model.data{i} = sample{snr(i)};
        case 4
            model.mp{model.Mmp+(i-1)*model.Nsmp+1}.name =  ['Z at Tdg, sample ',num2str(i)];
            model.mp{model.Mmp+(i-1)*model.Nsmp+1}.vmin = 0;
            model.mp{model.Mmp+(i-1)*model.Nsmp+1}.vmax = 0.05; %m
            model.mp{model.Mmp+(i-1)*model.Nsmp+2}.name =  ['dT2, sample ',num2str(i)];
            model.mp{model.Mmp+(i-1)*model.Nsmp+2}.vmin = 0;
            model.mp{model.Mmp+(i-1)*model.Nsmp+2}.vmax = model.tc; %Myr
            model.mp{model.Mmp+(i-1)*model.Nsmp+3}.name =  ['dz2, sample ',num2str(i)];
            model.mp{model.Mmp+(i-1)*model.Nsmp+3}.vmin = -1; %0.1 m
            model.mp{model.Mmp+(i-1)*model.Nsmp+3}.vmax = 1.0; %1; %log 10 m
            model.mp{model.Mmp+(i-1)*model.Nsmp+4}.name =  ['dT3, sample ',num2str(i)];
            model.mp{model.Mmp+(i-1)*model.Nsmp+4}.vmin = 0;
            model.mp{model.Mmp+(i-1)*model.Nsmp+4}.vmax = model.tc; %Myr
            model.mp{model.Mmp+(i-1)*model.Nsmp+5}.name =  ['dz3, sample ',num2str(i)];
            model.mp{model.Mmp+(i-1)*model.Nsmp+5}.vmin = -1;%0.1 m
            model.mp{model.Mmp+(i-1)*model.Nsmp+5}.vmax = 1.0;%1;%log 10 m
            model.mp{model.Mmp+(i-1)*model.Nsmp+6}.name =  ['dT4, sample ',num2str(i)];
            model.mp{model.Mmp+(i-1)*model.Nsmp+6}.vmin = 0;
            model.mp{model.Mmp+(i-1)*model.Nsmp+6}.vmax = model.tc;
            model.mp{model.Mmp+(i-1)*model.Nsmp+7}.name =  ['dz4, sample ',num2str(i)];
            model.mp{model.Mmp+(i-1)*model.Nsmp+7}.vmin = -1; %log10 -2 =0.01m/Myr
            model.mp{model.Mmp+(i-1)*model.Nsmp+7}.vmax = 1.0; %2=100 m/Myr
            model.mp{model.Mmp+(i-1)*model.Nsmp+8}.name =  ['E5, sample ',num2str(i)];
            model.mp{model.Mmp+(i-1)*model.Nsmp+8}.vmin = -2; %log10 -2 =0.01m/Myr
            model.mp{model.Mmp+(i-1)*model.Nsmp+8}.vmax = 2; %2=100 m/Myr
            
            %add sample info to models
            model.data{i} = sample{snr(i)};
    end
end

%save some general MCMC parameters
model.Nwalk = nwalkers; % Number of walkers
model.burnin = lburnin; % Length of burnin
model.Nmod = nmodacc; % Number of models
model.Nmax = nmodmax; % Maximum number of models

%number of total model parameters
model.Nmp = model.Mmp + model.Nsnr*model.Nsmp;


for i=1:model.Nmp
    umin(i) = model.mp{i}.vmin;
    umax(i) = model.mp{i}.vmax;
    du0(i) = model.mp{i}.vmax - model.mp{i}.vmin;
end
umin = umin(:);
umax = umax(:);

%data and covariance

% dobs=zeros(1,model.Nds*model.Nsnr); sigd=zeros(1,model.Nds*model.Nsnr); %initiate data observation and error vectors
% for i=1:model.Nsnr %check
%     for j=1:model.Ndp
%         dobs((i-1)*model.Nds+j) = model.data{i}.N10(j);
%         dobs((i-1)*model.Nds+j+model.Ndp) = model.data{i}.N26(j);
%         sigd((i-1)*model.Nds+j) = model.data{i}.dN10(j);
%         sigd((i-1)*model.Nds+j+model.Ndp) = model.data{i}.dN26(j);
%     end
% end

%possibility of only one nuclide
dobs=zeros(1,model.Nds*model.Nsnr); sigd=zeros(1,model.Nds*model.Nsnr); %initiate data observation and error vectors
if model.Nnc == 2
    for i=1:model.Nsnr %check
        for j=1:model.Ndp
            dobs((i-1)*model.Nds+j) = model.data{i}.N10(j);
            dobs((i-1)*model.Nds+j+model.Ndp) = model.data{i}.N26(j);
            sigd((i-1)*model.Nds+j) = model.data{i}.dN10(j);
            sigd((i-1)*model.Nds+j+model.Ndp) = model.data{i}.dN26(j);
        end
    end
elseif model.Nnc == 1
    for i=1:model.Nsnr
        for j=1:model.Ndp
            dobs((i-1)*model.Nds+j) = model.data{i}.N10(j);
            sigd((i-1)*model.Nds+j) = model.data{i}.dN10(j);
        end
    end
end

Cobs = model.Temp*diag(sigd.^2);
Cobsinv = inv(Cobs);

%%%%%% Initiate random models generation %%%%%%%

%Initialize random sequence
rng('default');

%walker starting points
if model.Nwalk > 1
    for i=1:model.Nmp
        wini(i,:) = 0.8*(randperm(model.Nwalk)-1)/(model.Nwalk-1)+0.1;
    end
else
    wini = randi([0,4],model.Nmp,model.Nwalk)/4;
end

%wini(:,1) = 0.5;
%wini(1,:) = 1;
%wini(4,:) = 0;

%loop walkers
for nw = 1:model.Nwalk
    %walker starting point - for initial parameter vector
    for i = 1:model.Nmp
        u(i) = (1-wini(i,nw))*model.mp{i}.vmin + wini(i,nw)*model.mp{i}.vmax;
    end
    
    switch numdp
        case 3
            Tdg = u(2);
            dT2 = u(4);
            dT3 = u(6);
            T4 = Tdg+dT2+dT3;
            while T4 > model.age %Only use the starting model, if the total age doesnt exceed model.age
                if model.Nwalk > 1
                    for i=1:model.Nmp
                        wini(i,:) = 0.8*(randperm(model.Nwalk)-1)/(model.Nwalk-1)+0.1;
                    end
                else
                    wini = randi([0,4],model.Nmp,model.Nwalk)/4;
                end
                for i = 1:model.Nmp
                    u(i) = (1-wini(i,nw))*model.mp{i}.vmin + wini(i,nw)*model.mp{i}.vmax;
                end
                Tdg = u(2);
                dT2 = u(4);
                dT3 = u(6);
                T4 = Tdg + dT2 + dT3;
            end
        case 4
            Tdg = u(2);
            dT2 = u(4);
            dT3 = u(6);
            dT4 = u(8);
            T5 = Tdg + dT2 + dT3 + dT4;
            while T5 > model.age %Only use the starting model, if the total age doesnt exceed model.age
                if model.Nwalk > 1
                    for i=1:model.Nmp
                        wini(i,:) = 0.8*(randperm(model.Nwalk)-1)/(model.Nwalk-1)+0.1;
                    end
                else
                    wini = randi([0,4],model.Nmp,model.Nwalk)/4;
                end
                for i = 1:model.Nmp
                    u(i) = (1-wini(i,nw))*model.mp{i}.vmin + wini(i,nw)*model.mp{i}.vmax;
                end
                Tdg = u(2);
                dT2 = u(4);
                dT3 = u(6);
                dT4 = u(8);
                T5 = Tdg + dT2 + dT3 + dT4;
            end
    end
    
    %initialize
    minres = 1e20; %not used
    res_current = 1e20; %current residual - first model run
    restot = 0;
    acount = 0;
    bcount = 0;
    rcount = 0;
    accrat = 0;
    status = zeros(model.Nmax,1);
    erosion_rec = zeros(model.Nmax,1);
    up_rec = zeros(model.Nmp,model.Nmax);
    u_rec = zeros(model.Nmp,model.Nmax);
    N10_rec = zeros(model.Nmax,1);
    N26_rec = zeros(model.Nmax,1);
    restot_rec = zeros(model.Nmax,1);
    accrat_rec = zeros(model.Nmax,1);
    k_rec = zeros(model.Nmax,1);
    duR = zeros(model.Nmp,1);
    
    
    accfac = 1e-2; % how does changing this one change step length (k)?
    
    ktarget = 0.025; %not used
    
    %run models
    mi = 0; %model iteration
    k = 0.01; %initial step length
    
    while ((mi < model.Nmax)&&(acount < model.Nmod))
        
        mi = mi + 1;
        
        if rem(mi,100) == 0 %only display output in terminal for every 100th model
            disp(['nw = ',num2str(nw),'/',num2str(model.Nwalk),' mi = ',num2str(mi),'/',num2str(model.Nmax),' bcount = ',num2str(bcount),' acount = ',num2str(acount),' accrat = ',num2str(accrat),' k = ',num2str(k)]);
        end
        %***** step length ******
        
        
        % acceptance ratio
        if (mi > k_meanlength)
            accrat = (sum(abs(status((mi-k_meanlength):(mi-1))))+1)/k_meanlength;
        elseif (bcount < model.burnin)
            accrat = 0.1; %target acceptance ratio during burnin /ELP: not used? see below
        else %it never enters this part of the loop?
            accrat = 0.30; %target acceptance ratio for the inversion /ELP: not used? see below
        end
        
%         if (mi > 100)
%             accrat = (sum(abs(status((mi-100):(mi-1))))+1)/100;
%         elseif (bcount < model.burnin)
%             accrat = 0.1; %target acceptance ratio during burnin
%         else
%             accrat = 0.3; %target acceptance ratio for the inversion
%         end
        
        %step length
        if (bcount < 0.5*model.burnin)
            
            k = k*((1-accfac) + accfac*accrat/0.1); %is target acceptance ratio defined here?
            
            if ((mi > 100)&&(bcount < 2)) k = 1.0;
            elseif ((mi > 10)&&(bcount < 2)) k = 0.01*mi;
            end
            
        elseif (bcount < model.burnin)
            
            k = k*((1-accfac) + accfac*accrat/0.2);
            
            if ((mi > 100)&&(bcount < 2)) k = 1.0;
            elseif ((mi > 10)&&(bcount < 2)) k = 0.01*mi;
            end
            
            %k = 2*ktarget;
            
        elseif (acount < model.Nmod)
            
            %k = k*((1-accfac) + accfac*accrat/0.3);
            k = k*((1-accfac) + accfac*accrat/acctarget); %target acceptance ratio
            
            %k = k*(accrat + 0.1)/(0.3 + 0.1)
            
            %k = ktarget;
            
        end
        
        if (k > 0.5)
            k = 0.5;
        end
        
        if (bcount < model.burnin)
            model.Temp = 1.0 + 10.0*(model.burnin-bcount)/model.burnin;
        else
            model.Temp = 1.0;
        end
        
        %********* propose new parameters ************
        
        %random step
        du = 0.5*randn(model.Nmp,1).*du0(:);
        
        %proposed model
        up = u(:) + k*du(:);
        
        
        %retake if model age and depth limits are exceeded
        switch numdp
            case 3
                while (any(up(:) < umin(:))||any(up(:) > umax(:))||((up(2)+up(4)+up(6))>model.age) ) %added final clause to make sure time intervals dont exceed model.age
                    
                    %random step
                    du = 0.5*randn(model.Nmp,1).*du0(:);
                    
                    %proposed model
                    up = u(:) + k*du(:);
                    
                end
            case 4
                while (any(up(:) < umin(:))||any(up(:) > umax(:))||((up(2)+up(4)+up(6)+up(8))>model.age) ) %added final clause to make sure time intervals dont exceed model.age
                    
                    %random step
                    du = 0.5*randn(model.Nmp,1).*du0(:);
                    
                    %proposed model
                    up = u(:) + k*du(:);
                    
                end
        end
        
        %********** Forward model *****************
        [gmp,time,burial] = forward_bedrock_vE1(up,model,CNprop,numdp,CNprod); %gmp: predicted data vector.
        
        gmp = gmp(1:model.Nds);
        
        %Acceptance criteria
        restot = (dobs(:)-gmp(:))'*Cobsinv*(dobs(:)-gmp(:));
        rfrac = exp(-0.5*restot)/exp(-0.5*res_current);
        alpha = rand(1);
        
        %if model is accepted
        if ((alpha < rfrac)||(mi == 1))
            
            u = up;
            gm = gmp;
            res_current = restot;
            
            %accepted after burnin
            if (bcount > model.burnin)
                
                status(mi) = 1;
                acount = acount + 1;
                
                %if accepted during burnin
            else
                
                status(mi) = -1;
                bcount = bcount + 1;
                
            end
            
            %rejected
        else
            
            status(mi) = 0;
            
            %rejected
            rcount = rcount + 1;
            
        end
        
        %save things
        up_rec(:,mi) = up(:);
        u_rec(:,mi) = u(:);
        gm_rec(:,mi) = gm(:);
        restot_rec(mi) = res_current;
        accrat_rec(mi) = accrat;
        k_rec(mi) = k;
        
    end
    
    %change status flag for models that are rejected during burnin
    Imin=find(status == 1,1);  % find index of first accepted model after burnin
    I = find(status(1:Imin) == 0); %find indices of rejected models in burnin
    status(I) = -2; %set status of these models to -2
    
    model.walker{nw}.status = status(1:mi);
    model.walker{nw}.up = up_rec(:,1:mi);
    model.walker{nw}.u = u_rec(:,1:mi);
    model.walker{nw}.gm = gm_rec(:,1:mi);
    model.walker{nw}.restot = restot_rec(1:mi);
    model.walker{nw}.acount = acount;
    model.walker{nw}.bcount = bcount;
    model.walker{nw}.rcount = rcount;
    model.walker{nw}.accrate = accrat_rec(1:mi);
    model.walker{nw}.kstep = k_rec(1:mi);
    
end

%save output

str = num2str(snr(1));
if length(snr) > 1
    for i=2:length(snr)
        str = [str,['-',num2str(snr(i))]];
    end
end

%savefile = ['models/Synthetic/Syn_sample_',str,'_BeAl_dp_Temp.mat'];
save(savepath,'model','-v7.3');