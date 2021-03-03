function [CNprop] = getCNprop();

%Cosmogenic props
TBe = 1.387e6;
TAl = 0.705e6;
CNprop.lambda_Be = log(2)/TBe;
CNprop.lambda_Al = log(2)/TAl;
CNprop.rho = 2.65; %density
CNprop.PBe0 = 4.0; %normalization production
CNprop.pratio = 6.75; %26Al/10Be surface production ratio
CNprop.pr_fm_Be = 0.005; %fast muons (of total p)
CNprop.pr_fm_Al = 0.006;
CNprop.pr_nmc_Be = 0.015;
CNprop.pr_nmc_Al = 0.018;
CNprop.Lspal = 150; %g/cm^2
CNprop.Lnmc = 1500;
CNprop.Lfm = 4320;
