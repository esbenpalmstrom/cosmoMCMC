function compile_data_mn()
 %Reads data info from excelfile and calculate production parameters
close all;

addpath('Functions','Functions/cl36','Functions/cl36/scaling','Functions/CronusCalc') 

%Read data info from Excelfile
ns = 1; %number of samples in data file

excelfile = ['data/InputTest.xlsx'];
sheet = 'MDML2';
[num,text,~] = xlsread(excelfile,sheet);

ij=1; %row start
for i=1:ns
    sample{i}.Nnuc = num(ij,9); %Number of nuclides for sample i
    sample{i}.nuclides = num(ij:ij+sample{i}.Nnuc-1,10); %nuclide identifications
    
    %Sample specific data
    sample{i}.name = text(ij+1,1);
    sample{i}.type = text(ij+1,2);
    sample{i}.site = text(ij+1,3);
    sample{i}.lat=num(ij,1);
    sample{i}.lon=num(ij,2);
    sample{i}.elev=num(ij,3);
    sample{i}.shield=num(ij,4);
    sample{i}.thick=num(ij,5);
    sample{i}.density = num(ij,6);
    sample{i}.depth=num(ij,7);
    sample{i}.sampleyr = num(ij,8);
    
    %Site-specific production parameters
    rho = sample{i}.density; %Density of rock sample (g/cm3)
    %Define depths below surface z/rho cm/(g/cm3) [g/cm^2] for fitting of production profiles
    D_m = 100; %Depth, changed below (rho)
    z_m = linspace(0,10,100);
    z_D = D_m*z_m.^3/10*rho; %denser depth-grid near surface
    keyboard
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


%     nuclide specific data
        for j=1:sample{i}.Nnuc
            nuclide=num(ij+j-1,10); %get nuclide identification
            switch nuclide
                case 1 %10Be
                sample{i}.N10 = num(ij+j-1,11);
                sample{i}.dN10 = num(ij+j-1,12);
                sample{i}.production.P10total = num(ij+j-1,13);
                sample{i}.production.P10spal = 0.98*sample{i}.production.P10total; %later combined, then split with other factors
                %sample{i}.production.P10muon = 0.02*sample{i}.production.P10total; %later combined, then split with other factors
                sample{i}.T10 = num(ij+j-1,14); 
                
                %Muons
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
                
                case 2 %26Al
                sample{i}.N26 = num(ij+j-1,11);
                sample{i}.dN26 = num(ij+j-1,12);
                sample{i}.production.P26total = num(ij+j-1,13);
                sample{i}.production.P26spal = 0.98*sample{i}.production.P26total; %later combined, then split with other factors
                %sample{i}.production.P26muon = 0.02*sample{i}.production.P26total; %later combined, then split with other factors   

                %Muons
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
                
                case 3 %36Cl
                sample{i}.N36=num(ij+j-1,11);
                sample{i}.dN36=num(ij+j-1,12);
                sample{i}.chem=num(ij+j-1,15:end);
                
                %Assemble Cl36 input data and calculate production parameters
                Cl36chem = [sample{i}.N36 zeros(1,2) sample{i}.chem(1) sample{i}.density sample{i}.thick sample{i}.lat,sample{i}.lon,sample{i}.elev p sample{i}.shield Lspal sample{i}.chem(2:end) sample{i}.depth sample{i}.sampleyr];
                Cl36prod=get36clProd(Cl36chem,z_D(z_D<1000)); %calc production rate from spallation, muons, thermal and epithermal neutrons

                %Spallation 36Cl (at/g/yr)
                sample{i}.production.P36spal=Cl36prod.P36s(1)*sf_spal;

                %Muons 36Cl (at/kg/yr)
                p36_muons = fit_P_exp(z_D(z_D<1000),Cl36prod.P36nm,Cl36prod.P36fm,2,0); %Fit muon production profile with two exponentials
                sample{i}.production.P36_Lnmc=p36_muons.L(2);
                sample{i}.production.P36_Lfm=p36_muons.L(1);
                shield_fac36_nm = exp(-sample{i}.thick/2*rho/p36_muons.L(2));
                shield_fac36_fm = exp(-sample{i}.thick/2*rho/p36_muons.L(1));
                sample{i}.production.P36_nmc = p36_muons.P(2)*shield_fac36_nm;
                sample{i}.production.P36_fm = p36_muons.P(1)*shield_fac36_fm;
            %     plot(Cl36prod.P36nm,z_D(z_D<1000)/rho,'.-',Cl36prod.P36fm,z_D(z_D<1000)/rho,'.-'), set(gca,'ydir','reverse'),legend('Negative muons','Fast muons')

                % Thermal neutron production 36Cl (at/kg/yr) - 4 terms 
                % Term 1-3
                sample{i}.production.P36_Lth1=Cl36prod.Lth1;
                sample{i}.production.P36_Lth2=Cl36prod.Lth2; 
                sample{i}.production.P36_Lth3=Cl36prod.Lth3; 
                shield_fac36_th1 = exp(-sample{i}.thick/2*rho/Cl36prod.Lth1);
                shield_fac36_th2 = exp(-sample{i}.thick/2*rho/Cl36prod.Lth2);
                shield_fac36_th3 = exp(-sample{i}.thick/2*rho/Cl36prod.Lth3);
                sample{i}.production.P36_th1 = Cl36prod.P36th1(1)*shield_fac36_th1;
                sample{i}.production.P36_th2 = Cl36prod.P36th2(1)*shield_fac36_th2;
                sample{i}.production.P36_th3 = Cl36prod.P36th3(1)*shield_fac36_th3;
            %     plot(Cl36prod.P36th3,z_D(z_D<1000)/rho,'.-'), set(gca,'ydir','reverse','xscale','log'),
            %     hold on, plot(Cl36prod.P36th3(1).*exp(-z_D(z_D<1000)./Cl36prod.Lth3),z_D(z_D<1000)/rho,'.-')

                %Epithermal 36Cl (at/kg/yr) - 3 terms% Term 1-2
                sample{i}.production.P36_Leth1=Cl36prod.Lth1; %same lengthscale as thermal1
                sample{i}.production.P36_Leth2=Cl36prod.Lth2;%same lengthscale as thermal2
                sample{i}.production.P36_eth1 = Cl36prod.P36eth1(1)*shield_fac36_th1;
                sample{i}.production.P36_eth2 = Cl36prod.P36eth2(1)*shield_fac36_th2;
            %     plot(Cl36prod.P36eth1,z_D(z_D<1000)/rho,'.-',Cl36prod.P36eth2,z_D(z_D<1000)/rho,'.-',Cl36prod.P36eth3,z_D(z_D<1000)/rho,'.-',Cl36prod.P36eth,z_D(z_D<1000)/rho,'k-'), set(gca,'ydir','reverse'),legend('Epithermal1','Epithermal2','Epithermal3','Epithermal'), grid on
            %     hold on, plot(Cl36prod.P36eth2(1).*exp(-z_D(z_D<1000)./Cl36prod.Lth2),z_D(z_D<1000)/rho,'.-')        

                % Non-exponential thermal+epithermal production terms, fit with two exponentials
                p36_theth = fit_P_exp(z_D(z_D<1000),Cl36prod.P36eth3+Cl36prod.P36th4,0,2,0);
                sample{i}.production.P36_Ltheth1=p36_theth.L(1);
                sample{i}.production.P36_Ltheth2=p36_theth.L(2);
                shield_fac36_theth1 = exp(-sample{i}.thick/2*rho/p36_theth.L(1));
                shield_fac36_theth2 = exp(-sample{i}.thick/2*rho/p36_theth.L(2));
                sample{i}.production.P36_theth1 = p36_theth.P(1)*shield_fac36_theth1;
                sample{i}.production.P36_theth2 = p36_theth.P(2)*shield_fac36_theth2;
            %     plot(Cl36prod.P36th1,z_D(z_D<1000)/rho,'.-',Cl36prod.P36th2,z_D(z_D<1000)/rho,'.-',Cl36prod.P36th3,z_D(z_D<1000)/rho,'.-',Cl36prod.P36th4,z_D(z_D<1000)/rho,'.-',Cl36prod.P36th,z_D(z_D<1000)/rho,'k-'), set(gca,'ydir','reverse'),legend('Thermal1','Thermal2','Thermal3','Thermal4','Thermal')
            %     plot(Cl36prod.P36th3,z_D(z_D<1000)/rho,'.-'), set(gca,'ydir','reverse','xscale','log'),%legend('Thermal1','Thermal2','Thermal3','Thermal4','Thermal')
 
                case 4 %21Ne
                sample{i}.N21=num(ij+j-1,11);
                sample{i}.dN21=num(ij+j-1,12);
                sample{i}.production.P21total=num(ij+j-1,13);
                sample{i}.production.P21spal = 0.964*sample{i}.production.P21total;
                %sample{i}.production.P21muon = 0.036*sample{i}.production.P21total;
                
                % Muons 21Ne. These are not calibrated, taken from Fernandez-Mosquera 2010
                mc21.Natoms = 1.0003e22; %Si atoms pr gram Quartz
                mc21.k_neg = 0.296.*0.6559.*0.0029; 
                mc21.sigma190 = 0.79e-27; 
                p21_muons = fit_Pmu_with_exp_1A(p,min(z_D),maxZ,2,mc21,0);
                sample{i}.production.P21_Lfm=p21_muons.L(1);
                sample{i}.production.P21_Lnmc=p21_muons.L(2);
                shield_fac21_nm = exp(-sample{i}.thick/2*rho/p21_muons.L(2));
                shield_fac21_fm = exp(-sample{i}.thick/2*rho/p21_muons.L(1));
                sample{i}.production.P21_nmc = p21_muons.P(2)*shield_fac21_nm;
                sample{i}.production.P21_fm = p21_muons.P(1)*shield_fac21_fm;

            end
        end
        
        %Calculate nuclide ratios
        if isfield(sample{i},'N26')
            sample{i}.r2610 = sample{i}.N26./sample{i}.N10;
            sample{i}.dr2610 = sample{i}.r2610.*sqrt((sample{i}.dN10./sample{i}.N10).^2+(sample{i}.dN26./sample{i}.N26).^2);
        end       
        if isfield(sample{i},'N21')
            sample{i}.r2110 = sample{i}.N21./sample{i}.N10;
            sample{i}.dr2110 = sample{i}.r2110.*sqrt((sample{i}.dN10./sample{i}.N10).^2+(sample{i}.dN21./sample{i}.N21).^2);
        end
        if isfield(sample{i},'N36')
            sample{i}.r3610 = sample{i}.N36./sample{i}.N10;
            sample{i}.dr3610 = sample{i}.r3610.*sqrt((sample{i}.dN10./sample{i}.N10).^2+(sample{i}.dN36./sample{i}.N36).^2);
        end
        
    ij=ij+sample{i}.Nnuc; %Look for next sample in this row
end

save ./data/InputTest.mat sample excelfile