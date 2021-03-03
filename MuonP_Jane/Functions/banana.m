
function banana(CNprop)

%
% banana() produces a cosmo banana plot
% 
% DLE 31/1 2019
%

%close all;
%set(gcf,'units','centimeters','paperunits','centimeters');
%set(gcf,'position',[20,20,30,20]);



%************* Cosmo parameters **************


Ptot_Be = CNprop.PBe0;
Ptot_Al = Ptot_Be*CNprop.pratio;
Pfm_Be = CNprop.pr_fm_Be*Ptot_Be;
Pnmc_Be = CNprop.pr_nmc_Be*Ptot_Be;
Pspal_Be = Ptot_Be - Pfm_Be - Pnmc_Be;
Pfm_Al = CNprop.pr_fm_Al*Ptot_Al;
Pnmc_Al = CNprop.pr_nmc_Al*Ptot_Al;
Pspal_Al = Ptot_Al - Pfm_Al - Pnmc_Al;
Lspal = CNprop.Lspal; %g/cm^2
Lnmc = CNprop.Lnmc;
Lfm = CNprop.Lfm;
lambda_Be = CNprop.lambda_Be;
lambda_Al = CNprop.lambda_Al;
rho = CNprop.rho; %density g/cm3

%********************************************

time = linspace(0,1e7,2000)';
time(2)=1;

%erates = [0,1e-7,1e-6,1e-5,1e-4,1e-3];
erates = [0,1e-7,5e-6,1e-6,1e-5,2e-5,1e-3,1];

NBe = zeros(length(time),length(erates));
NAl = zeros(length(time),length(erates));

for i=1:length(erates),

    %Full exposure with erosion
    fBe = lambda_Be + rho*erates(i)*100/Lspal;
    fAl = lambda_Al + rho*erates(i)*100/Lspal;
    NBe(:,i) = Pspal_Be/fBe*(1-exp(-fBe*time));
    NAl(:,i) = Pspal_Al/fAl*(1-exp(-fAl*time));

    fBe = lambda_Be + rho*erates(i)*100/Lnmc;
    fAl = lambda_Al + rho*erates(i)*100/Lnmc;
    NBe(:,i) = NBe(:,i) + Pnmc_Be/fBe*(1-exp(-fBe*time));
    NAl(:,i) = NAl(:,i) + Pnmc_Al/fAl*(1-exp(-fAl*time));

    fBe = lambda_Be + rho*erates(i)*100/Lfm;
    fAl = lambda_Al + rho*erates(i)*100/Lfm;
    NBe(:,i) = NBe(:,i) + Pfm_Be/fBe*(1-exp(-fBe*time));
    NAl(:,i) = NAl(:,i) + Pfm_Al/fAl*(1-exp(-fAl*time));
    
    
end;
    
%No erosion    
line(NBe(:,1),NAl(:,1)./NBe(:,1),'color','k');


%with erosion
for i=2:length(erates),
    line(NBe(:,i),NAl(:,i)./NBe(:,i),'color','k','linestyle','--');
end;


%infinite time with erosion
fero = logspace(-8,2,400);
NBe_f = Pspal_Be./(lambda_Be + rho*fero/Lspal)+Pnmc_Be./(lambda_Be + rho*fero/Lnmc)+Pfm_Be./(lambda_Be + rho*fero/Lfm);
NAl_f = Pspal_Al./(lambda_Al + rho*fero/Lspal)+Pnmc_Al./(lambda_Al + rho*fero/Lnmc)+Pfm_Al./(lambda_Al + rho*fero/Lfm);;
line(NBe_f,NAl_f./NBe_f,'color','k');


%burial isolines
tburial = [0.5e6,1e6,1.5e6,2e6];
NBe_b = zeros(length(time),length(tburial));
NAl_b = zeros(length(time),length(tburial));
for i=1:length(tburial),
   NBe_b(:,i) = NBe(:,1).*exp(-lambda_Be*tburial(i)); 
   NAl_b(:,i) = NAl(:,1).*exp(-lambda_Al*tburial(i)); 
   line(NBe_b(:,i),NAl_b(:,i)./NBe_b(:,i),'color','k','linestyle',':');
end;

%burial paths
tb = linspace(0,3e6,400);
etime = [1e3,1e4,1e5,1e6,1e7]; %preburial exposure times
NBe_p = zeros(length(tb),length(etime));
NAl_p = zeros(length(tb),length(etime));
for i=1:length(etime),
    NBe_p(:,i) = Ptot_Be/lambda_Be*(1-exp(-lambda_Be*etime(i)))*exp(-lambda_Be*tb);
    NAl_p(:,i) = Ptot_Al/lambda_Al*(1-exp(-lambda_Al*etime(i)))*exp(-lambda_Al*tb);

    line(NBe_p(:,i),NAl_p(:,i)./NBe_p(:,i),'color','k','linestyle',':');    
end;



