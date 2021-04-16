function showhist(snr)

close all;
set(groot','defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
addpath data

%Cosmogenic props
% TBe = 1.387e6;
% TAl = 0.705e6;
% CNprop.lambda_Be = log(2)/TBe;
% CNprop.lambda_Al = log(2)/TAl;
% CNprop.rho = 2.65; %density
% CNprop.PBe0 = 4.0; %normalization production
% CNprop.pratio = 6.75; %26Al/10Be surface production ratio
% CNprop.pr_fm_Be = 0.005; %fast muons (of total p)
% CNprop.pr_fm_Al = 0.006;
% CNprop.pr_nmc_Be = 0.015;
% CNprop.pr_nmc_Al = 0.018;
% CNprop.Lspal = 150; %g/cm^2
% CNprop.Lnmc = 1500;
% CNprop.Lfm = 4320;

%load MC results
str = num2str(snr(1));
if length(snr) > 1
    for i=2:length(snr)
        str = [str,['-',num2str(snr(i))]];
    end
end

%mname = ['models/focus_samples_SyntheticDepths/focus_sample_',str,'_BeAl_surface.mat'];
%mname = ['models/FSbase_gausta27_Jan_2021_20_45_31_noThree/gausta_data_noThree.mat'];
modelfolder = ['4'];
mname = ['models/' modelfolder '/gausta_data_2.mat'];

load(mname);

col1 = 0.9*[1,1,1];
col2 = 0.6*[1,1,1];
col3 = 0.3*[1,1,1];

%colors
cg = [.6,.6,.9];
cig = [.6,.9,.6];

map = colormap;
nc = length(map);
xc = linspace(0,1,nc);
rc = map(:,1);
gc = map(:,2);
bc = map(:,3);
for i=1:model.Nwalk
    ri = interp1(xc,rc,(i-1)/model.Nwalk);
    gi = interp1(xc,gc,(i-1)/model.Nwalk);
    bi = interp1(xc,bc,(i-1)/model.Nwalk);
    wcol(i,:) = [ri,gi,bi];
end



%******* figures showing walkers *************

set(gcf,'units','pixels','position',[1000,1000,800,800]);
set(gcf,'color',[1,1,1]);
set(gca,'position',[0,0,1,1],'visible','off');
ax0 = gca;
text(0.5,0.035,'Time before present (kyr)','horizontalalignment','center','fontsize',20);
text(0.05,0.825,'$\delta_{18}$O ($10^{-3}$)','rotation',90,'horizontalalignment','center','fontsize',20);
text(0.05,0.35,'Depth below surface (m)','rotation',90,'horizontalalignment','center','verticalalignment','middle','fontsize',20);
%text(0.925,0.1,'20','horizontalalignment','center','verticalalignment','middle','fontsize',24);
%text(0.915,0.6,'0','horizontalalignment','center','verticalalignment','middle','fontsize',24);
text(0.92,0.65,'P($\delta_{18}$O)','horizontalalignment','center','fontsize',20);
text(0.935,0.25,'Normalized P','rotation',90,'horizontalalignment','center','fontsize',2);


ax1 = axes('position',[0.1,0.675,0.75,0.28]);
hold on; box on; %grid on;
set(gca,'ydir','reverse','layer','top');
set(gca,'xtick',[0:500:5000],'xticklabel',[]);
set(gca,'ytick',[3:5]);
set(gca,'fontsize',20);
axis([0,model.age*1e3,3,5.3]);

ax2 = axes('position',[0.1,0.1,0.75,0.5]);
hold on; box on; %grid on;
set(gca,'ydir','reverse','layer','top');
set(gca,'xtick',[0:500:5000]);
set(gca,'ytick',[0:5:20]); %changed to 20 from 10
set(gca,'fontsize',20);
axis([0,model.age*1e3,0,20]);%changed to 20 from 10
title(model.data{1}.name);

ax3 = axes('position',[0.86,0.675,0.12,0.28]);
hold on; box on; %grid on;
axis([0,1,3,5.3]);
set(gca,'ydir','reverse','layer','top');
set(gca,'ytick',[3:5]);
set(gca,'xticklabel',[],'yticklabel',[]);


% ********* d18O history ********
axes(ax1);
load d18Ocurves.mat
tta = Age*1e-3;
dd = d18O_10ky;

%extract dO18 distribution

%loop walkers
uval = [];
for nw = 1:model.Nwalk
    I = find(model.walker{nw}.status == 1);
    if (length(I) > 0)
        uval = [uval(:);model.walker{nw}.up(1,I)'];
    end
end
if (length(uval) > 0)
    [f,xi] = ksdensity(uval);
end
fint = cumsum(f);
fint = fint/max(fint);

ddg = dd;
ddig = dd;
I = find(dd < xi(1)); ddig(I) = xi(1);
patch([tta(:);flipud(tta(:))],[dd(:);flipud(ddig(:))],cig,'linestyle','none');

Nxi = length(xi);
on = ones(size(dd));
for i=1:Nxi-1
    top = xi(i)*on;
    base = xi(i+1)*on;
    I = find((dd > xi(i))&(dd < xi(i+1))); top(I) = dd(I);
    I = find(dd > xi(i+1)); top(I) = xi(i+1);
    cf = fint(i);
    col = cf*[1,1,1]+(1-cf)*cig;
    patch([tta(:);flipud(tta(:))],[top(:);flipud(base(:))],col,'linestyle','none');
    
    
    top = xi(i)*on;
    base = xi(i+1)*on;
    I = find((dd < xi(i+1))&(dd > xi(i))); base(I) = dd(I);
    I = find(dd < xi(i)); base(I) = xi(i);
    cf = fint(i);
    col = (1-cf)*[1,1,1]+cf*cg;
    patch([tta(:);flipud(tta(:))],[top(:);flipud(base(:))],col,'linestyle','none');
end;

ddg = dd;
ddig = dd;
I = find(dd > xi(end)); ddig(I) = xi(end);
patch([tta(:);flipud(tta(:))],[dd(:);flipud(ddig(:))],cg,'linestyle','none');


%draw d18O curve
l2 = line(tta,dd,'color','k');

%plot dO18 distribution
axes(ax3);
fval = 0.9*f/max(f);
patch([fval(:);zeros(size(fval(:)))],[xi(:);flipud(xi(:))],cg,'linestyle','none');
line(fval,xi,'color','k');


axes(ax2);

Maxz = model.z0;
Maxt = model.age*1e3;
z0 = model.z0;
Nz = 100;
Nt = 200;
zint = linspace(0,Maxz,Nz);
tint = linspace(0,Maxt,Nt);
[tbin,zbin]=meshgrid(tint,zint);

%[np,nn] = numSubplots(model.Nsnr);

for ns = 1:model.Nsnr
    
    %subplot(np(1),np(2),ns);
    %hold on; box on;
    
    %set(gca,'ydir','reverse');
    %xlabel('Time (Ma)');
    %ylabel('Burial depth');
    %title(['Sample ',num2str(ns)]);
    
    histgrid = zeros(size(tbin));
    dzi = zint(2)-zint(1);
    
    n0 = model.Mmp + (ns-1)*model.Nsmp;
    
    for nw=1:model.Nwalk
        I = find(model.walker{nw}.status == 1);
        for i=1:length(I)
            
            switch model.Nfree
                case 2
                    
                    T1 = model.walker{nw}.u(2,I(i));
                    z1 = model.walker{nw}.u(n0+1,I(i));
                    dT2 = model.walker{nw}.u(n0+2,I(i));
                    dz2 = 10^model.walker{nw}.u(n0+3,I(i));
                    dT3 = model.walker{nw}.u(n0+4,I(i));
                    dz3 = 10^model.walker{nw}.u(n0+5,I(i));
                    E4 = 10^model.walker{nw}.u(n0+6,I(i));
                    
                    T2 = T1 + dT2;
                    z2 = z1 + dz2;
                    T3 = T2 + dT3;
                    z3 = z2 + dz3;
                    T4 = model.age; %this requires age > Tdg+dT2+dT3
                    z4 = z3 + (model.age - T3)*E4;
                    
                    
                    Tm = [0,T1,T2,T3,T4]*1e3;
                    zm = [0,z1,z2,z3,z4];
                    
                    zinterp = interp1(Tm,zm,tint,'linear',2*z0);
                    izs = 1+floor(zinterp/dzi); %bin index at each tsfine
                    for itsfine=1:Nt
                        if izs(itsfine)<=Nz
                            histgrid(izs(itsfine),itsfine)=histgrid(izs(itsfine),itsfine)+ 1;
                        end
                    end
                    
                case 3
                    %keyboard
                    T1 = model.walker{nw}.u(2,I(i));
                    z1 = model.walker{nw}.u(n0+1,I(i));
                    dT2 = model.walker{nw}.u(n0+2,I(i));
                    dz2 = 10^model.walker{nw}.u(n0+3,I(i));
                    dT3 = model.walker{nw}.u(n0+4,I(i));
                    dz3 = 10^model.walker{nw}.u(n0+5,I(i));
                    dT4 = model.walker{nw}.u(n0+6,I(i));
                    dz4 = 10^model.walker{nw}.u(n0+7,I(i));
                    E5 = 10^model.walker{nw}.u(n0+8,I(i));
                    
                    
                    T2 = T1 + dT2;
                    z2 = z1 + dz2;
                    T3 = T2 + dT3;
                    z3 = z2 + dz3;
                    T4 = T3 + dT4; %this requires age > Tdg+dT2+dT3
                    z4 = z3 + dz4; %hvis T3 bliver større end model.age kan der opstå problemer
                    z5 = z4 + (model.age - T4)*E5;
                    T5 = model.age;
                    
                    Tm = [0,T1,T2,T3,T4,T5]*1e3;
                    zm = [0,z1,z2,z3,z4,z5];
                    
                    zinterp = interp1(Tm,zm,tint,'linear',2*z0);
                    izs = 1+floor(zinterp/dzi); %bin index at each itsfine
                    for itsfine=1:Nt
                        if izs(itsfine)<=Nz %JLA added second clause 11.03.20(quickfix)
                            histgrid(izs(itsfine),itsfine)=histgrid(izs(itsfine),itsfine)+ 1;
                        end
                    end
                    
            end
        end
    end
    
    histgrid = histgrid/max(histgrid(:));
    
    cvals = linspace(-3,0,100);
    [C,h] = contourf(tbin,zbin,log10(histgrid+1e-6),cvals,'linestyle','none');
    
    map = colormap;
    map(1,:) = [1,1,1];
    colormap(map);
    caxis([-3,0]);
    set(gca,'layer','top');
    cb = colorbar;
    set(cb,'fontsize',16);
    set(cb,'position',[0.87,0.1,0.03,0.3],'ytick',[-2,-1,0],'xticklabel',[0.01,0.10,1.00]);
    
    
    
end;


% mname = ['models/SGreenland/histplots/hist_sample_',str,'.gif'];
% mname = ['models/MDML/histplots/hist_sample_',str,'BeAlCl.gif'];
%mname = ['models/Synthetic/hist_sample_',str,'BeAl_dp_Temp.gif'];
%mname = ['models/stroeven/hist_sample_',str,'BeAl_dp_Temp.gif'];
%mname = ['models/focus_samples_SyntheticDepths/hist_sample_',str,'BeAl.gif'];
%mname = ['models/FSbase_gausta27_Jan_2021_20_45_31_noThree/gausta.gif'];

mname = ['models/' modelfolder '/hist.gif'];

frame = getframe(gcf);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,mname,'gif');


set(groot','defaulttextinterpreter','default');
set(groot, 'defaultAxesTickLabelInterpreter','default');
set(groot, 'defaultLegendInterpreter','default');