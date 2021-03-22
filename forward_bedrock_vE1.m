function [gm,time,burial] = forward_bedrock_vE1(up,model,CNprop,numdp,CNprod)
% function [gm] = forward_bedrock(up,models)
%
% Forward model for transient CN integration
% of bedrock sample,
% Modified from forward_bedrockv7 (DLE) by JLA in June 2019
%
% Updated July 2019, DLE: fixed 1e6 error on erate time, tested full?
% integrated solution in-stead of using midpoint - no noticeable difference
%
% vJ2 updated Feb. 2020 by JLA includes the possibility of extra data points at depth
% up = model parameter vector
% gm = predicted data vector
% CNprop = Cosmogenic Nuclide properties
% numdp = number of depth points in the geological model, currently 3 or 4
% are possible
% CNprod = method of CN production, 'old' or 'new'

%model parameters

%{
Esben notes, sep. 2020

INPUT:
up: model parameter vector containing timings and dephts.
model: parameter vector containing various other info. Also contains "up"?
CNprop: Cosmogenic Nuclide properties. See getCNprop.m

OUTPUT:
gm: Vector of data predicted from the given proposed model.
time:
burial:
%}

%generic glaciation parameters
d18Op = up(1); %d18O threshold
T1 = up(2); %time of deglaciation

%load and prepare d18O data
load('d18Ocurves.mat');
I = find(Age < model.age*1e6);
dO = d18O_10ky(I);
dOt = Age(I);

%time steps
model.dt = 1000; %yr


%loop samples
for i=1:model.Nsnr
    
    n0 = (i-1)*model.Nsmp+model.Mmp; %parameter number start
    
    switch numdp
        case 3
            z1 = up(n0+1);
            dT2 = up(n0+2); %Myr
            dz2 = 10^up(n0+3); %m
            dT3 = up(n0+4);
            dz3 = 10^up(n0+5); %m
            E4 = 10^up(n0+6); %m/Myr
            
            T2 = T1 + dT2;
            z2 = z1 + dz2;
            
            T3 = T2 + dT3;
            z3 = z2 + dz3;
            
            
            T4 = model.age; %this requires age > Tdg+dT2+dT3
            z4 = z3 + (model.age - T3)*E4;
            
            mT = [0,T1,T2,T3,T4]*1e6; %Myr to yr
            mz = [0,z1,z2,z3,z4];
            
            % Check starting condition of model
            if (z4 > model.z0) % Depth at start of model greater than max depth
                maxtime = interp1(mz,mT,model.z0);
                ssbc = 0;
            else
                maxtime = model.age*1e6;
                ssbc = 1;
            end
        case 4
            z1 = up(n0+1);
            dT2 = up(n0+2); %Myr
            dz2 = 10^up(n0+3); %m
            dT3 = up(n0+4);
            dz3 = 10^up(n0+5); %m
            dT4 = up(n0+6); %Myr
            dz4 = 10^up(n0+7); %m
            E5 = 10^up(n0+8); %m/Myr
            
            T2 = T1 + dT2;
            z2 = z1 + dz2;
            
            T3 = T2 + dT3;
            z3 = z2 + dz3;
            
            T4 = T3 + dT4;
            z4 = z3 + dz4;
            
            T5 = model.age; %this requires age > Tdg+dT2+dT3+dT4
            z5 = z4 + (model.age - T4)*E5;
            
            mT = [0,T1,T2,T3,T4,T5]*1e6; %Myr to yr
            mz = [0,z1,z2,z3,z4,z5];
            
            % Check starting condition of model
            if (z5 > model.z0) % Depth at start of model greater than max depth
                maxtime = interp1(mz,mT,model.z0);
                ssbc = 0;
            else
                maxtime = model.age*1e6;
                ssbc = 1;
            end
    end
    
    nt = ceil(maxtime/model.dt);
    time = maxtime*linspace(0,1,nt);
    
    
    %**************************************
    % Exhumation history
    %**************************************
    %     burial = interp1(mT,mz,time); %Depths at times in 'time' %old
    
    Ndp=model.Ndp; %Number of data points in depth profile
    depths=model.data{i}.depths; %depths of datapoints
    
    burial = zeros(Ndp,nt); %initiate 'burial' matrix with depths at times in 'time', for all data sample points in depth profile
    burial(1,:) = interp1(mT,mz,time); %Surface depths at times in 'time'
    if Ndp > 1 % if there are samples at depth
        %bug - if the top point, is also a depth point, this loop needs to
        %take account of that as well.
        burial(2:end,:) = burial(1,:)+depths(2:end).*ones(Ndp-1,nt); %add their depths below surface to 'burial'
        
        
    end
    %below loop added by ELP, 29.01.202
    %if top depth is not in the surface, it should be able to take
    %account of that.
    if model.data{1}.depths(1) ~= 0
        burial(1,:) = burial(1,:) + depths(1).*ones(1,nt);
    end
    dOn = interp1(dOt,dO,time,'linear',dO(end)); %d18O values at 'time'
    pfac = ones(nt,1); %controls exposure, modified below
%     I = find(dOn > d18Op); pfac(I) = 0; %No exposure when d18O
    I = dOn > d18Op; pfac(I) = 0; %No exposure when d18O
    %values above threshold (d18Op)
%     I = find(time < 25e3); pfac(I) = 0; %correct exposure around
    I = time < 25e3; pfac(I) = 0; %correct exposure around
    %deglaciation, set no
    %exposure last 25 kyr
%     I = find(time < T1*1e6); pfac(I) = 1; %then add exposure from
    I = time < T1*1e6; pfac(I) = 1; %then add exposure from
    %time of deglaciation
    
    
    %*****CN production*****
    
    switch CNprod
        case 'old'
            P10spal = model.data{i}.P10spal;
            P10tot = model.data{i}.P10spal + model.data{i}.P10muon;
            P10fm = CNprop.pr_fm_Be*P10tot;
            P10nmc = CNprop.pr_nmc_Be*P10tot;
            
            P26spal = model.data{i}.P26spal;
            P26tot = model.data{i}.P26spal + model.data{i}.P26muon;
            P26fm = CNprop.pr_fm_Al*P26tot;
            P26nmc = CNprop.pr_nmc_Al*P26tot;
            
            rho = CNprop.rho;
            Lspal = CNprop.Lspal;
            P10Lnmc = CNprop.Lnmc;
            P26Lnmc = CNprop.Lnmc;
            P10Lfm = CNprop.Lfm;
            P26Lfm = CNprop.Lfm;
            
        case 'new'
            rho = model.data{i}.density;
            Lspal = model.data{i}.production.Lspal;
            
            P10spal = model.data{i}.production.P10spal;
            P10fm = model.data{i}.production.P10_fm;
            P10Lfm = model.data{i}.production.P10_Lfm;
            P10nmc = model.data{i}.production.P10_nmc;
            P10Lnmc = model.data{i}.production.P10_Lnmc;
            
            
            P26spal = model.data{i}.production.P26spal;
            P26fm = model.data{i}.production.P26_fm;
            P26Lfm = model.data{i}.production.P26_Lfm;
            P26nmc = model.data{i}.production.P26_nmc;
            P26Lnmc = model.data{i}.production.P26_Lnmc;
    end
    
    N10 = zeros(Ndp,nt); % changed '1' to 'Ndp'
    N26 = zeros(Ndp,nt); % --
    
    
    if (ssbc == 0) %depth of sample at model start greater than maxdepth
        
        N10(:,nt) = 0.0; % changed '(nt)' to '(:,nt)'
        N26(:,nt) = 0.0; % --
        
    else
        
        %assume steady state concentration at starting point
        switch numdp
            case 3
                erate = E4*1e-6;
            case 4
                erate = E5*1e-6;
        end
        
        switch CNprod
            case 'old'
                %erate = E5*1e-6;
                %spallation
                fBe = CNprop.lambda_Be + CNprop.rho*erate*100/CNprop.Lspal;
                fAl = CNprop.lambda_Al + CNprop.rho*erate*100/CNprop.Lspal;
                N10(:,nt) = P10spal*exp(-CNprop.rho*100*burial(:,nt)/CNprop.Lspal)/fBe; % changed '(nt)' to '(:,nt)'
                N26(:,nt) = P26spal*exp(-CNprop.rho*100*burial(:,nt)/CNprop.Lspal)/fAl; % changed '(nt)' to '(:,nt)'
                %negative muon capture
                fBe = CNprop.lambda_Be + CNprop.rho*erate*100/CNprop.Lnmc;
                fAl = CNprop.lambda_Al + CNprop.rho*erate*100/CNprop.Lnmc;
                N10(:,nt) = N10(:,nt) + P10nmc*exp(-CNprop.rho*100*burial(:,nt)/CNprop.Lnmc)/fBe; % changed '(nt)' to '(:,nt)'
                N26(:,nt) = N26(:,nt) + P26nmc*exp(-CNprop.rho*100*burial(:,nt)/CNprop.Lnmc)/fAl; % changed '(nt)' to '(:,nt)'
                %fast muon capture
                fBe = CNprop.lambda_Be + CNprop.rho*erate*100/CNprop.Lfm;
                fAl = CNprop.lambda_Al + CNprop.rho*erate*100/CNprop.Lfm;
                N10(:,nt) = N10(:,nt) + P10fm*exp(-CNprop.rho*100*burial(:,nt)/CNprop.Lfm)/fBe; % changed '(nt)' to '(nt,:)'
                N26(:,nt) = N26(:,nt) + P26fm*exp(-CNprop.rho*100*burial(:,nt)/CNprop.Lfm)/fAl; % changed '(nt)' to '(nt,:)'
                keyboard
                
            case 'new'
                %erate = E5*1e-6;
                %spallation
                fBe = CNprop.lambda_Be + rho*erate*100/Lspal;
                fAl = CNprop.lambda_Al + rho*erate*100/Lspal;
                N10(:,nt) = P10spal*exp(-rho*100*burial(:,nt)/Lspal)/fBe; % changed '(nt)' to '(:,nt)'
                N26(:,nt) = P26spal*exp(-rho*100*burial(:,nt)/Lspal)/fAl; % changed '(nt)' to '(:,nt)'
                %negative muon capture
                fBe = CNprop.lambda_Be + rho*erate*100/P10Lnmc;
                fAl = CNprop.lambda_Al + rho*erate*100/P26Lnmc;
                N10(:,nt) = N10(:,nt) + P10nmc*exp(-rho*100*burial(:,nt)/P10Lnmc)/fBe; % changed '(nt)' to '(:,nt)'
                N26(:,nt) = N26(:,nt) + P26nmc*exp(-rho*100*burial(:,nt)/P26Lnmc)/fAl; % changed '(nt)' to '(:,nt)'
                %fast muon capture
                fBe = CNprop.lambda_Be + rho*erate*100/P10Lfm;
                fAl = CNprop.lambda_Al + CNprop.rho*erate*100/P26Lfm;
                N10(:,nt) = N10(:,nt) + P10fm*exp(-rho*100*burial(:,nt)/P10Lfm)/fBe; % changed '(nt)' to '(nt,:)'
                N26(:,nt) = N26(:,nt) + P26fm*exp(-rho*100*burial(:,nt)/P26Lfm)/fAl; % changed '(nt)' to '(nt,:)'
        end
    end
    
    %integrate time
    for kk=(nt-1):-1:1
        %dt = (time(kk+1)-time(kk));
        %pf = .5*(pfac(kk+1)+pfac(kk));
        %bz = .5*(burial(kk+1)+burial(kk));
        %P10z = pf*(P10spal*exp(-CNprop.rho*100*bz/CNprop.Lspal)+P10nmc*exp(-CNprop.rho*100*bz/CNprop.Lnmc)+P10fm*exp(-CNprop.rho*100*bz/CNprop.Lfm));
        %P26z = pf*(P26spal*exp(-CNprop.rho*100*bz/CNprop.Lspal)+P26nmc*exp(-CNprop.rho*100*bz/CNprop.Lnmc)+P26fm*exp(-CNprop.rho*100*bz/CNprop.Lfm));
        %N10(kk) = N10(kk+1)*exp(-dt*CNprop.lambda_Be) + P10z*(1-exp(-dt*CNprop.lambda_Be))/CNprop.lambda_Be;
        %N26(kk) = N26(kk+1)*exp(-dt*CNprop.lambda_Al) + P26z*(1-exp(-dt*CNprop.lambda_Al))/CNprop.lambda_Al;
        
        
        dt = (time(kk+1)-time(kk));
        pf = .5*(pfac(kk+1)+pfac(kk));
        bz = burial(:,kk); % changed '(kk)' to '(:,kk)', Depths at time step kk
        erate = 100*(burial(1,kk+1)-burial(1,kk))/dt; % changed '(kk)' to '(1,kk)'
        
        P10z_spal = pf*P10spal*exp(-rho*100*bz/Lspal);
        P10z_nmc = pf*P10nmc*exp(-rho*100*bz/P10Lnmc);
        P10z_fm = pf*P10fm*exp(-rho*100*bz/P10Lfm);
        
        P26z_spal = pf*P26spal*exp(-rho*100*bz/Lspal);
        P26z_nmc = pf*P26nmc*exp(-rho*100*bz/P26Lnmc);
        P26z_fm = pf*P26fm*exp(-rho*100*bz/P26Lfm);
        
        
        N10(:,kk) = N10(:,kk+1)*exp(-dt*CNprop.lambda_Be); % changed '(kk)' to '(:,kk+1)'
        
        ff = CNprop.lambda_Be+rho*erate/Lspal;
        N10(:,kk) = N10(:,kk) + P10z_spal*(1.0-exp(-ff*dt))/ff; % changed '(kk)' to '(:,kk)'
        
        ff = CNprop.lambda_Be+rho*erate/P10Lnmc;
        N10(:,kk) = N10(:,kk) + P10z_nmc*(1.0-exp(-ff*dt))/ff; % changed '(kk)' to '(:,kk)'
        
        ff = CNprop.lambda_Be+rho*erate/P10Lfm;
        N10(:,kk) = N10(:,kk) + P10z_fm*(1.0-exp(-ff*dt))/ff; % changed '(kk)' to '(:,kk)'
        
        N26(:,kk) = N26(:,kk+1)*exp(-dt*CNprop.lambda_Al); % changed '(kk)' to '(:,kk+1)'
        
        ff = CNprop.lambda_Al+rho*erate/Lspal;
        N26(:,kk) = N26(:,kk) + P26z_spal*(1.0-exp(-ff*dt))/ff; % changed '(kk)' to '(:,kk)'
        
        ff = CNprop.lambda_Al+rho*erate/P26Lnmc;
        N26(:,kk) = N26(:,kk) + P26z_nmc*(1.0-exp(-ff*dt))/ff; % changed '(kk)' to '(:,kk)'
        
        ff = CNprop.lambda_Al+rho*erate/P26Lfm;
        N26(:,kk) = N26(:,kk) + P26z_fm*(1.0-exp(-ff*dt))/ff; % changed '(kk)' to '(:,kk)'
        
    end
    
    %modeled data
    %     gm((i-1)*model.Nds+1) = N10(1); %old
    %     gm((i-1)*model.Nds+2) = N26(1); %old
    
    gm=zeros(1,model.Nds*model.Nsnr);
    for j=1:Ndp
        gm((i-1)*model.Nds+j) = N10(j,1);
        gm((i-1)*model.Nds+Ndp+j) = N26(j,1);
    end
end