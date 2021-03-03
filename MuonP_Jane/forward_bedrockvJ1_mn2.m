function [gm,time,burial] = forward_bedrockvJ1_mn2(up,model,CNprop)

% function [gm] = forward_bedrock(up,models)
%
% Forward model for transient CN integration
% of bedrock sample, 

% Modified from forward_bedrockv7 (DLE) by JLA in June 2019
% Updated July 2019, DLE: fixed 1e6 error on erate time, tested full?
% integrated solution in-stead of using midpoint - no noticeable difference
% Multi-nuclide handling March 2020, JLA
% Be-Al production parameters updated, April 2020, JLA

% up = model parameter vector
% gm = predicted data vector

%model parameters

%generic glaciation parameters
d18Op = up(1);
T1 = up(2);

%load and prepare d18O data
load('zachos_.mat');
I = find(t < model.age*1e6);
dO = dO_4ky(I);
dOt = t(I);

% load('d18Ocurves.mat');
% I = find(Age < model.age*1e6);
% dO = d18O_4ky(I);
% dOt = Age(I);

%time steps 
model.dt = 1000; %yr


%loop samples
for i=1:model.Nsnr
    n0 = (i-1)*model.Nsmp+model.Mmp; %parameter number start
    
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
    if (z5 > model.z0); % Depth at start of model greater than max depth
        maxtime = interp1(mz,mT,model.z0);
        ssbc = 0;
    else
        maxtime = model.age*1e6;
        ssbc = 1;
    end;
    
    nt = ceil(maxtime/model.dt); 
    time = maxtime*linspace(0,1,nt);
    
    
    %**************************************
    % Exhumation history
    %**************************************
    burial = interp1(mT,mz,time); %Depths at times in 'time'
    dOn = interp1(dOt,dO,time,'linear',dO(end)); %d18O values at 'time'
    pfac = ones(nt,1); %controls exposure, modified below
    I = find(dOn > d18Op); pfac(I) = 0; %No exposure when d18O
                                        %values above threshold (d18Op)
    I = find(time < 25e3); pfac(I) = 0; %correct exposure around
                                        %deglaciation, set no
                                        %exposure last 25 kyr
    I = find(time < T1*1e6); pfac(I) = 1; %then add exposure from
                                          %time of deglaciation
                                      
   
    %CN production
    rho = model.data{i}.density;
    Lspal = model.data{i}.production.Lspal;
    
    if ismember(1,model.data{i}.nuclides) %10Be
        P10spal = model.data{i}.production.P10spal;
        P10fm = model.data{i}.production.P10_fm;
        P10Lfm = model.data{i}.production.P10_Lfm;
        P10nmc = model.data{i}.production.P10_nmc;
        P10Lnmc = model.data{i}.production.P10_Lnmc;
    end
    
    
    if ismember(2,model.data{i}.nuclides) %26Al
        P26spal = model.data{i}.production.P26spal;
        P26fm = model.data{i}.production.P26_fm;
        P26Lfm = model.data{i}.production.P26_Lfm;
        P26nmc = model.data{i}.production.P26_nmc;
        P26Lnmc = model.data{i}.production.P26_Lnmc;
    end
    
    
    if ismember(3,model.data{i}.nuclides) %36Cl
        P36spal = model.data{i}.production.P36spal;
        P36nmc = model.data{i}.production.P36_nmc;
        P36Lnmc = model.data{i}.production.P36_Lnmc;
        P36fm = model.data{i}.production.P36_fm;
        P36Lfm = model.data{i}.production.P36_Lfm;
        P36th1 = model.data{i}.production.P36_th1;
        P36Lth1 = model.data{i}.production.P36_Lth1;
        P36th2 = model.data{i}.production.P36_th2;
        P36Lth2 = model.data{i}.production.P36_Lth2;
        P36th3 = model.data{i}.production.P36_th3;
        P36Lth3 = model.data{i}.production.P36_Lth3;
        P36eth1 = model.data{i}.production.P36_eth1;
        P36Leth1 = model.data{i}.production.P36_Leth1;
        P36eth2 = model.data{i}.production.P36_eth2;
        P36Leth2 = model.data{i}.production.P36_Leth2;
        P36theth1 = model.data{i}.production.P36_theth1;
        P36Ltheth1 = model.data{i}.production.P36_Ltheth1;
        P36theth2 = model.data{i}.production.P36_theth2;
        P36Ltheth2 = model.data{i}.production.P36_Ltheth2;
    end
    
    
    if ismember(4,model.data{i}.nuclides) %21Ne
        P21spal = model.data{i}.production.P21spal;
        P21fm = model.data{i}.production.P21_fm;
        P21Lfm = model.data{i}.production.P21_Lfm;
        P21nmc = model.data{i}.production.P21_nmc;
        P21Lnmc = model.data{i}.production.P21_Lnmc;
    end
    
    
    N10 = zeros(nt,1);
    N26 = zeros(nt,1);
    N36 = zeros(nt,1);
    N21 = zeros(nt,1);
    
    if (ssbc == 0) %depth of sample at model start greater than maxdepth
        
        N10(nt) = 0.0;
        N26(nt) = 0.0;
        N36(nt) = 0.0;
        N21(nt) = 0.0;

    else
        %assume steady state concentration at starting point
        
        erate = E5*1e-6;
        
        if ismember(1,model.data{i}.nuclides) %10Be
            fBe = CNprop.lambda_Be + rho*erate*100/Lspal; %spallation
            N10(nt) = P10spal*exp(-rho*100*burial(nt)/Lspal)/fBe;
            fBe = CNprop.lambda_Be + rho*erate*100/P10Lnmc; %negative muon capture
            N10(nt) = N10(nt) + P10nmc*exp(-rho*100*burial(nt)/P10Lnmc)/fBe;
            fBe = CNprop.lambda_Be + rho*erate*100/P10Lfm; %fast muons
            N10(nt) = N10(nt) + P10fm*exp(-rho*100*burial(nt)/P10Lfm)/fBe;
        end
        
        if ismember(2,model.data{i}.nuclides) %26Al
            fAl = CNprop.lambda_Al + rho*erate*100/Lspal; %spallation
            N26(nt) = P26spal*exp(-rho*100*burial(nt)/Lspal)/fAl;
            fAl = CNprop.lambda_Al + rho*erate*100/P26Lnmc; %negative muon capture
            N26(nt) = N26(nt) + P26nmc*exp(-rho*100*burial(nt)/P26Lnmc)/fAl;
            fAl = CNprop.lambda_Al + rho*erate*100/P26Lfm; %fast muons
            N26(nt) = N26(nt) + P26fm*exp(-rho*100*burial(nt)/P26Lfm)/fAl;
        end
        
        if ismember(3,model.data{i}.nuclides) %36Cl
            fCl = CNprop.lambda_Cl + rho*erate*100/Lspal; %spallation
            N36(nt) = P36spal*exp(-rho*100*burial(nt)/Lspal)/fCl;
            fCl = CNprop.lambda_Cl + rho*erate*100/P36Lnmc; %negative muon capture
            N36(nt) = N36(nt) + P36nmc*exp(-rho*100*burial(nt)/P36Lnmc)/fCl;
            fCl = CNprop.lambda_Cl + rho*erate*100/P36Lfm; %fast muons
            N36(nt) = N36(nt) + P36fm*exp(-rho*100*burial(nt)/P36Lfm)/fCl;
            %thermal and epithermal neutrons (36Cl)
            fCl = CNprop.lambda_Cl + rho*erate*100/P36Lth1;
            N36(nt) = N36(nt) + P36th1*exp(-rho*100*burial(nt)/P36Lth1)/fCl;
            fCl = CNprop.lambda_Cl + rho*erate*100/P36Lth2;
            N36(nt) = N36(nt) + P36th2*exp(-rho*100*burial(nt)/P36Lth2)/fCl;
            fCl = CNprop.lambda_Cl + rho*erate*100/P36Lth3;
            N36(nt) = N36(nt) + P36th3*exp(-rho*100*burial(nt)/P36Lth3)/fCl;
            fCl = CNprop.lambda_Cl + rho*erate*100/P36Leth1;
            N36(nt) = N36(nt) + P36eth1*exp(-rho*100*burial(nt)/P36Leth1)/fCl;
            fCl = CNprop.lambda_Cl + rho*erate*100/P36Leth2;
            N36(nt) = N36(nt) + P36eth2*exp(-rho*100*burial(nt)/P36Leth2)/fCl;
            fCl = CNprop.lambda_Cl + rho*erate*100/P36Ltheth1;
            N36(nt) = N36(nt) + P36theth1*exp(-rho*100*burial(nt)/P36Ltheth1)/fCl;
            fCl = CNprop.lambda_Cl + rho*erate*100/P36Ltheth2;
            N36(nt) = N36(nt) + P36theth2*exp(-rho*100*burial(nt)/P36Ltheth2)/fCl;
        end
    
        if ismember(4,model.data{i}.nuclides) %21Ne
            fNe = rho*erate*100/Lspal; %spallation
            N21(nt) = P21spal*exp(-rho*100*burial(nt)/Lspal)/fNe;
            fNe = rho*erate*100/P21Lnmc; %negative muon capture
            N21(nt) = N21(nt) + P21nmc*exp(-rho*100*burial(nt)/P21Lnmc)/fNe;
            fNe = rho*erate*100/P21Lfm; %fast muons
            N21(nt) = N21(nt) + P21fm*exp(-rho*100*burial(nt)/P21Lfm)/fNe;
        end
    end
    
        
    %integrate time
    for kk=(nt-1):-1:1
    
        %dt = (time(kk+1)-time(kk));
        %pf = .5*(pfac(kk+1)+pfac(kk));
        %bz = .5*(burial(kk+1)+burial(kk));
        %P10z = pf*(P10spal*exp(-CNprop.rho*100*bz/Lspal)+P10nmc*exp(-CNprop.rho*100*bz/CNprop.Lnmc)+P10fm*exp(-CNprop.rho*100*bz/CNprop.Lfm));
        %P26z = pf*(P26spal*exp(-CNprop.rho*100*bz/Lspal)+P26nmc*exp(-CNprop.rho*100*bz/CNprop.Lnmc)+P26fm*exp(-CNprop.rho*100*bz/CNprop.Lfm));
        %N10(kk) = N10(kk+1)*exp(-dt*CNprop.lambda_Be) + P10z*(1-exp(-dt*CNprop.lambda_Be))/CNprop.lambda_Be;
        %N26(kk) = N26(kk+1)*exp(-dt*CNprop.lambda_Al) + P26z*(1-exp(-dt*CNprop.lambda_Al))/CNprop.lambda_Al;
        
        dt = (time(kk+1)-time(kk)); %length of timestep
        pf = .5*(pfac(kk+1)+pfac(kk)); %exposure: 1=full exposure, 0=no exposure (glacial cover)
        bz = burial(kk); %depth at time kk
        erate = 100*(burial(kk+1)-burial(kk))/dt; %erosion rate within timestep
        
        % Be-10
        if ismember(1,model.data{i}.nuclides) %10Be
            %Production depth profiles, spallation, fast and negative muon pathways
            P10z_spal = pf*P10spal*exp(-rho*100*bz/Lspal);
            P10z_nmc = pf*P10nmc*exp(-rho*100*bz/P10Lnmc);
            P10z_fm = pf*P10fm*exp(-rho*100*bz/P10Lfm);

            N10(kk) = N10(kk+1)*exp(-dt*CNprop.lambda_Be); %decay from previous step

            ff = CNprop.lambda_Be+rho*erate/Lspal; 
            N10(kk) = N10(kk) + P10z_spal*(1.0-exp(-ff*dt))/ff; %add spallation production (Dunai 2010 eq. 4.10)

            ff = CNprop.lambda_Be+rho*erate/P10Lnmc;
            N10(kk) = N10(kk) + P10z_nmc*(1.0-exp(-ff*dt))/ff; %add negative muon capture production

            ff = CNprop.lambda_Be+rho*erate/P10Lfm;
            N10(kk) = N10(kk) + P10z_fm*(1.0-exp(-ff*dt))/ff; %add fast muon production
        end
        
        % Al-26
        if ismember(2,model.data{i}.nuclides) %26Al
            %Production depth profiles, spallation, fast and negative muon pathways
            P26z_spal = pf*P26spal*exp(-rho*100*bz/Lspal);
            P26z_nmc = pf*P26nmc*exp(-rho*100*bz/P26Lnmc);
            P26z_fm = pf*P26fm*exp(-rho*100*bz/P26Lfm);

            N26(kk) = N26(kk+1)*exp(-dt*CNprop.lambda_Al); %decay from previous step

            ff = CNprop.lambda_Al+rho*erate/Lspal;
            N26(kk) = N26(kk) + P26z_spal*(1.0-exp(-ff*dt))/ff; %add spallation production

            ff = CNprop.lambda_Al+rho*erate/P26Lnmc;
            N26(kk) = N26(kk) + P26z_nmc*(1.0-exp(-ff*dt))/ff; %add negative muon capture production

            ff = CNprop.lambda_Al+rho*erate/P26Lfm;
            N26(kk) = N26(kk) + P26z_fm*(1.0-exp(-ff*dt))/ff; %add fast muon production
        end
        
        % Cl-36   
        if ismember(3,model.data{i}.nuclides) %36Cl
            %Production depth profiles, spallation, fast and negative muon, thermal and epithermal pathways(36Cl)
            P36z_spal = pf*P36spal*exp(-rho*100*bz/Lspal);
            P36z_nmc = pf*P36nmc*exp(-rho*100*bz/P36Lnmc);
            P36z_fm = pf*P36fm*exp(-rho*100*bz/P36Lfm);
            P36z_th1 = pf*P36th1*exp(-rho*100*bz/P36Lth1);
            P36z_th2 = pf*P36th2*exp(-rho*100*bz/P36Lth2);
            P36z_th3 = pf*P36th3*exp(-rho*100*bz/P36Lth3);
            P36z_eth1 = pf*P36eth1*exp(-rho*100*bz/P36Leth1);
            P36z_eth2 = pf*P36eth2*exp(-rho*100*bz/P36Leth2);
            P36z_theth1 = pf*P36theth1*exp(-rho*100*bz/P36Ltheth1);
            P36z_theth2 = pf*P36theth2*exp(-rho*100*bz/P36Ltheth2);

            N36(kk) = N36(kk+1)*exp(-dt*CNprop.lambda_Cl); %decay from previous step

            ff = CNprop.lambda_Cl+rho*erate/Lspal;
            N36(kk) = N36(kk) + P36z_spal*(1.0-exp(-ff*dt))/ff; %add spallation production

            ff = CNprop.lambda_Cl+rho*erate/P36Lnmc;
            N36(kk) = N36(kk) + P36z_nmc*(1.0-exp(-ff*dt))/ff; %add negative muon capture production

            ff = CNprop.lambda_Cl+rho*erate/P36Lfm;
            N36(kk) = N36(kk) + P36z_fm*(1.0-exp(-ff*dt))/ff; %add fast muon production

            ff = CNprop.lambda_Cl+rho*erate/P36Lth1;
            N36(kk) = N36(kk) + P36z_th1*(1.0-exp(-ff*dt))/ff; %add thermal 1 production

            ff = CNprop.lambda_Cl+rho*erate/P36Lth2;
            N36(kk) = N36(kk) + P36z_th2*(1.0-exp(-ff*dt))/ff; %add thermal 2 production

            ff = CNprop.lambda_Cl+rho*erate/P36Lth3;
            N36(kk) = N36(kk) + P36z_th3*(1.0-exp(-ff*dt))/ff; %add thermal 3 production

            ff = CNprop.lambda_Cl+rho*erate/P36Leth1;
            N36(kk) = N36(kk) + P36z_eth1*(1.0-exp(-ff*dt))/ff; %add epithermal 1 production

            ff = CNprop.lambda_Cl+rho*erate/P36Leth2;
            N36(kk) = N36(kk) + P36z_eth2*(1.0-exp(-ff*dt))/ff; %add epithermal 2 production

            ff = CNprop.lambda_Cl+rho*erate/P36Ltheth1;
            N36(kk) = N36(kk) + P36z_theth1*(1.0-exp(-ff*dt))/ff; %add thermal-epithermal 1 production

            ff = CNprop.lambda_Cl+rho*erate/P36Ltheth2;
            N36(kk) = N36(kk) + P36z_theth2*(1.0-exp(-ff*dt))/ff; %add thermal-epithermal 2 production
        end
        
        % Ne-21  
        if ismember(4,model.data{i}.nuclides) %21Ne
            %Production depth profiles, spallation, fast and negative muon pathways
            P21z_spal = pf*P21spal*exp(-rho*100*bz/Lspal);
            P21z_nmc = pf*P21nmc*exp(-rho*100*bz/P21Lnmc);
            P21z_fm = pf*P21fm*exp(-rho*100*bz/P21Lfm);

            N21(kk) = N21(kk+1); %No decay - 21Ne is stable

            ff = rho*erate/Lspal;
            N21(kk) = N21(kk) + P21z_spal*(1.0-exp(-ff*dt))/ff; %add spallation production

            ff = rho*erate/P21Lnmc;
            N21(kk) = N21(kk) + P21z_nmc*(1.0-exp(-ff*dt))/ff; %add negative muon capture production

            ff = rho*erate/P21Lfm;
            N21(kk) = N21(kk) + P21z_fm*(1.0-exp(-ff*dt))/ff; %add fast muon production
        end
    end

    %modeled data    
    for j=1:model.data{i}.Nnuc
        nuclide=model.data{i}.nuclides(j); %get nuclide identification
            if nuclide == 1 %10Be
                gm((i-1)*model.Nds+j) = N10(1);
            elseif nuclide == 2 %26Al
                gm((i-1)*model.Nds+j) = N26(1);
            elseif nuclide == 3 %36Cl
                gm((i-1)*model.Nds+j) = N36(1);
            elseif nuclide == 4 %21Ne
                gm((i-1)*model.Nds+j) = N21(1);
            end
    end

end

