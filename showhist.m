%function showhist(snr,mname,figurename,figtitle)
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

%modelfolder = ['16'];
%mname = ['models/' modelfolder '/gausta_data_2.mat'];

% %synthetic base models
% fignr = 1;
% 
% mnames = {'Darm2008_1.mat','Ling2006a_11.mat','Nesj2007_9.mat','Stroe2002a_1.mat','Sven2015_12.mat'};
% figtitles = {'Pyha-Nattanen - Darmody et al. 2008'
%     'Gaustatoppen - Linge et al. 2006'
%     'Andoya - Nesje et al. 2007'
%     'Naakakarhakka - Stroeven et al. 2002'
%     'Karmoy - Svendsen et al. 2015'};
% figurenames = {'SyntheticBaseModelsDarmody2008'
%     'SyntheticBaseModelsLinge2006'
%     'SyntheticBaseModelsNesje2007'
%     'SyntheticBaseModelsStroeven2002'
%     'SyntheticBaseModelsSvendsen2015'};
% 
% 
% mname = ['models/FSbase_12_Dec_2020_13_04_14/' mnames{fignr}];
% figtitle = figtitles{fignr};
% figurename = figurenames{fignr};


mname = 'models/4/gausta_data_2.mat'; %path to .mat file
figurename = 'GaustaExhum004_extracted_models'; %name of file saved to figures folder
figtitle = 'Low production rate scheme'; %title of figure in plot


isSynth = 0; % 1 for on
load parametermeans.mat
load(mname);
dispmodels = 0;

showsmallexhum = 0; %1 for on with synth model, 2 for on without synth model. 0 for no small insert in the exhumation plot

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
text(0.5,0.035,'Time before present [kyr]','horizontalalignment','center','fontsize',20);
text(0.05,0.825,'$\delta_{18}$O [$10^{-3}$]','rotation',90,'horizontalalignment','center','fontsize',20);
text(0.05,0.35,'Depth below surface [m]','rotation',90,'horizontalalignment','center','verticalalignment','middle','fontsize',20);
%text(0.925,0.1,'20','horizontalalignment','center','verticalalignment','middle','fontsize',24);
%text(0.915,0.6,'0','horizontalalignment','center','verticalalignment','middle','fontsize',24);
text(0.92,0.65,'P($\delta_{18}$O)','horizontalalignment','center','fontsize',20);
text(0.96,0.25,'Normalized P','rotation',90,'horizontalalignment','center','fontsize',20);


ax1 = axes('position',[0.1,0.675,0.75,0.28]);
hold on; box on; %grid on;
set(gca,'ydir','reverse','layer','top');
set(gca,'xtick',[0:500:5000],'xticklabel',[]);
set(gca,'ytick',[3:5]);
set(gca,'fontsize',20);
axis([0,model.age*1e3,3,5.3]);
title(figtitle);

ax2 = axes('position',[0.1,0.1,0.75,0.5]); %exhumation plot axes
hold on; box on; %grid on;
set(gca,'ydir','reverse','layer','top');
set(gca,'xtick',[0:500:5000]);
set(gca,'ytick',[0:5:20]); %changed to 20 from 10
set(gca,'fontsize',20);
axis([0,model.age*1e3,0,20]);%changed to 20 from 10
%title(model.data{1}.name);

ax3 = axes('position',[0.86,0.675,0.12,0.28]);
hold on; box on; %grid on;
axis([0,1,3,5.3]);
set(gca,'ydir','reverse','layer','top');
set(gca,'ytick',[3:5]);
set(gca,'xticklabel',[],'yticklabel',[]);

if showsmallexhum == 1
    ax4 = axes('position',[0.15,0.15,0.25,0.2]); %small exhumation plot inside exhumation plot
    hold on; box on;
    set(gca,'ydir','reverse');
    set(gca,'xtick',[0:10:50]);
    set(gca,'ytick',[0:0.2:0.5]);
    set(gca,'fontsize',20);
    axis([0,25,0,0.5]);
end

if showsmallexhum == 2
    ax5 = axes('position',[0.45,0.15,0.25,0.2]);
    %ax5.Position = [0.55,0.35,0.25,0.2]; % for all data, high production
    ax5.Position = [0.61,0.42,0.20,0.16]; % for all data, low production
    hold on; box on;
    set(gca,'ydir','reverse');
    ax5.XLim = [0 150];
    ax5.YLim = [0 2];
    %set(gca,'xtick',[0:10:50]);
    %set(gca,'ytick',[0:0.2:0.5]);
    %set(gca,'fontsize',20);
    %axis([0,150,0,2]);

end
    


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
end

ddg = dd;
ddig = dd;
I = find(dd > xi(end)); ddig(I) = xi(end);
patch([tta(:);flipud(tta(:))],[dd(:);flipud(ddig(:))],cg,'linestyle','none');


%draw d18O curve
l2 = line(tta,dd,'color','k');
%plot synthetic data if relevant
if isSynth ==1
    synpmval = extractfield(parameterMeans,model.data{1}.name);
    plot([0 3000],[synpmval(1) synpmval(1)],'r--','LineWidth',2)
end

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
                    %line(tp(1:3),zp(1:3),'color',cg,'linewidth',4);
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
    
    if isSynth ==1
        syn_exhu(synpmval,[1,0,0],[.6,.9,.6])
    end
    
    if showsmallexhum == 1
        axes(ax4); hold on;
        contourf(tbin,zbin,log10(histgrid+1e-6),cvals,'linestyle','none');
        zp = [0.0,synpmval(3),synpmval(3)+10^(synpmval(5))];
        tp = [0.0,synpmval(2),synpmval(2)+synpmval(4)].*10e2;
        line(tp(1:3),zp(1:3),'color','r','linewidth',4,'LineStyle','--');
        mz = 15;
        plot(tp(1:2),zp(1:2),'o','markeredgecolor','k','markersize',mz,'markerfacecolor',cig);
        pos = get(gca,'position');
        xlim = get(gca,'xlim');
        ylim = get(gca,'ylim');
        
    end
    
    if showsmallexhum == 2
        axes(ax5); hold on;
        contourf(tbin,zbin,log10(histgrid+1e-6),cvals,'linestyle','none');
        ax5.FontSize = 12;
        
        ax5.XLabel.String = 'Time before present [kyr]';
        ax5.YLabel.String = 'Depth below surface [m]';
        %ax5.XLabel.FontSize = 14;
        %ax5.YLabel.FontSize = 14;
        
        axes(ax2)
        r = rectangle('Position',[0 0 150 2]);
        r.LineWidth = 2;
        axes(ax5)
    end
    
    if dispmodels == 1
        wlkn = 1;
        F = find(model.walker{wlkn}.u(1,:)<4.4);
        RF = F(randperm(length(F),1)); % get random number in f
        d18ORF = model.walker{wlkn}.u(1,RF);
        mdlRF = model.walker{wlkn}.u(:,RF);
        
        FF = find(model.walker{wlkn}.u(1,:)>5.25);
        RFF = FF(randperm(length(FF),1)); % get random number in ff
        d18ORFF = model.walker{wlkn}.u(1,RFF);
        mdlRFF = model.walker{wlkn}.u(:,RFF);
        
        
        axes(ax1);
        plot([0 3000],[d18ORF d18ORF],'r--','LineWidth',2)
        plot([0 3000],[d18ORFF d18ORFF],'m--','LineWidth',2)
        axes(ax2);
        syn_exhu(mdlRF,[1,0,0],[.6,.9,.6])
        syn_exhu(mdlRFF,[1,0,1],[.6,.9,.6])
        
    end
    
end



% mname = ['models/SGreenland/histplots/hist_sample_',str,'.gif'];
% mname = ['models/MDML/histplots/hist_sample_',str,'BeAlCl.gif'];
%mname = ['models/Synthetic/hist_sample_',str,'BeAl_dp_Temp.gif'];
%mname = ['models/stroeven/hist_sample_',str,'BeAl_dp_Temp.gif'];
%mname = ['models/focus_samples_SyntheticDepths/hist_sample_',str,'BeAl.gif'];
%mname = ['models/FSbase_gausta27_Jan_2021_20_45_31_noThree/gausta.gif'];

%mname = ['models/' modelfolder '/hist.gif'];

%frame = getframe(gcf);
%im = frame2im(frame);
%[imind,cm] = rgb2ind(im,256);
%imwrite(imind,cm,mname,'gif');
keyboard
export_fig(['/Users/esben/OneDrive - Aarhus Universitet/Speciale/skriv/latex/figures/' figurename],'-jpg','-r300');

set(groot','defaulttextinterpreter','default');
set(groot, 'defaultAxesTickLabelInterpreter','default');
set(groot, 'defaultLegendInterpreter','default');
