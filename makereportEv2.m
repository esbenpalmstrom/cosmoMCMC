function makereportEv2(snr,savepath,reportpath,sampleID,numdp)
%{
%based on makereportJv2.m
%input: snr,savepath,reportpath,sampleID
%example savepath: ['models/FS_' starttimestr '/' nuclidesetting depthsetting '/' sampleID '.mat']
%example reportpath: reportpath = ['models/FS_' starttimestr '/' nuclidesetting depthsetting '/reports/' sampleID '/']

WIP notes:
Make the ylim for exhumation figure dependent on total erosion, so that low
erosion scenarios are better displayed.
%}
close all

addpath Functions
addpath Functions/export_fig

CNprop = getCNprop;


str = num2str(snr(1));
if length(snr) > 1
    for i=2:length(snr)
        str = [str,['-',num2str(snr(i))]];
    end
end


mname = savepath;

sname = reportpath;

load(mname);

col1 = 0.9*[1,1,1];
col2 = 0.6*[1,1,1];
col3 = 0.3*[1,1,1];

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

%plot margins
lpm = 0.1;
ddm = 0.02;

%******* figures showing walkers *************

set(gcf,'units','normalized','position',[.1,.3,.4,.6]);
set(gca,'position',[0,0,1,1],'visible','off');
set(gca,'xlim',[0,1],'ylim',[0,1]);
set(gcf,'Name','Walker information');
ax0 = gca;

[np,nn] = numSubplots(model.Nmp);

%plot models params
for i = 1:model.Nmp
    subplot(np(1),np(2),i); hold on; box on; grid on;
    xlabel('model nr.');
    ylabel(model.mp{i}.name);
    %set(gca,'xlim',[0,model.Nmod],'ylim',[model.mp{i}.vmin,model.mp{i}.vmax]);
    set(gca,'xlim',[0,length(model.walker{1,1}.status)],'ylim',[model.mp{i}.vmin,model.mp{i}.vmax]);
    for nw = 1:model.Nwalk
        %         yyaxis left
        I = find(model.walker{nw}.status == 0);
        plot(I,model.walker{nw}.up(i,I),'.','color',col1);
        I = find(model.walker{nw}.status == -1);
        plot(I,model.walker{nw}.up(i,I),'.','color',col2);
        I = find(model.walker{nw}.status == 1);
        plot(I,model.walker{nw}.up(i,I),'.','color',wcol(nw,:));
        
        %plot development of residual
        %         yyaxis right
        %         resmean = movmean(model.walker{nw}.restot,499);
        %         plot(resmean,'-r');
        %         ylim([0 10])
        
    end
    if isfield(model,'synpmval') == 1
        plot([0,length(model.walker{1}.status)],[model.synpmval(i),model.synpmval(i)],'r-','LineWidth',2); %plot synthetic pm values
    end
end


print('temp3.pdf','-dpdf','-fillpage');

%set(gcf, 'Position', get(0, 'Screensize'));
filename = strcat(sname,'walkers.png');
saveas(gcf,filename)

%********************** Parameter distributions ************************
figure()
set(gcf,'units','normalized','position',[.3,.2,.4,.6]);
set(gcf,'Name','Parameter distributions');

[np,nn] = numSubplots(model.Nmp);

%loop parameters
for i=1:model.Nmp
    
    subplot(np(1),np(2),i);
    hold on; box on; grid on;
    ylabel('Frequency');
    xlabel(model.mp{i}.name);
    set(gca,'xlim',[model.mp{i}.vmin, ...
        model.mp{i}.vmax]);
    uval = [];
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status == 1);
        if (~isempty(I))
            [f,xi] = ksdensity(model.walker{nw}.up(i,I));
            line(xi,f,'color',wcol(nw,:));
            uval = [uval(:);model.walker{nw}.up(i,I)'];
        end
    end
    if (~isempty(uval))
        [f,xi] = ksdensity(uval);
    end
    line(xi,f,'color','k','linewidth',2);
    ylims = get(gca,'ylim');
    
    if isfield(model,'synpmval') == 1
        plot([model.synpmval(i),model.synpmval(i)],[ylims(1),ylims(2)],'r-','LineWidth',2); %plot synthetic pm values
    end
    
    %add median, 25% and 75% quartile lines
    med = median(uval);
    quart = quantile(uval,[0.25 0.75]);
    xline(med,'--b')
    xline(quart(1),'--r')
    xline(quart(2),'--r')
    
end



print('temp2.pdf','-dpdf','-fillpage');

%set(gcf, 'Position', get(0, 'Screensize'));
filename = strcat(sname,'parameters.jpg');
saveas(gcf,filename)

%******* report figure *************
figure;
set(gcf,'papertype','a4');
set(gcf,'units','centimeters','position',[5,5,21,29.7]);
set(gcf,'Name','Report');


ax0 = gca;
set(ax0,'position',[0,0,1,1]);
set(gca,'visible','off');
text(0.05,0.97,['File: ',mname],'HorizontalAlignment','left','fontsize',14);

dpx = (1-2*lpm-(model.Nsnr-1)*ddm)/model.Nsnr;

%Al/Be vs Be plot. Or alternatively if no Al, just a Be histogram
if model.Nnc == 2
    for ns = 1:model.Nsnr
        
        axes('position',[lpm+(ns-1)*(ddm+dpx),0.675,dpx,0.25]);
        hold on; box on; grid on;
        
        title(model.data{ns}.name);
        
        n0 = (ns-1)*model.Nds;
        Ndp = model.Ndp;
        
        Be = [];
        Al = [];
        
        for nw = 1:model.Nwalk
            I = find(model.walker{nw}.status > -1);
            Be = [Be(:,:);model.walker{nw}.gm(n0+1:Ndp+n0,I)'];
            if model.Nnc == 2
                Al = [Al(:,:);model.walker{nw}.gm(n0+Ndp+1:n0+2*Ndp,I)'];
            elseif model.Nnc == 1
                Al = NaN(size(Be));
            end
        end
        AlBe = Al./Be;
        
        Beint = linspace(0.5*min(min(Be)),max(max(Be)),40)';
        Alint = linspace(min(min(Al)),max(max(Al)),40)';
        AlBeint = linspace(0.5*min(min(AlBe)),1.5*max(max(AlBe)),40)';
        
        %banana(CNprop);
        xlabel(['NBe (atoms/g), sample ',num2str(ns)]);
        if (ns == 1), ylabel(['NAl/NBe']);
        else, set(gca,'yticklabel',[]);
        end
        
        if model.Nnc == 2
            for j=1:Ndp
                N = hist3([AlBe(:,j),Be(:,j)],{AlBeint Beint});
                N = N/sum(N(:));
                
                nfac = CNprop.PBe0/(model.data{ns}.P10spal + model.data{ns}.P10muon);
                [X,Y] = meshgrid(Beint*nfac,AlBeint);
                contour(X,Y,N,40);
            end
            
            errorbar(model.data{ns}.N10*nfac,model.data{ns}.r2610,model.data{ns}.dN10*nfac,'horizontal','.k');
            errorbar(model.data{ns}.N10*nfac,model.data{ns}.r2610,model.data{ns}.dr2610,'vertical','.k');
        end
        
    end
end

if model.Nnc == 1
    for ns = 1:model.Nsnr
        axes('position',[lpm+(ns-1)*(ddm+dpx),0.675,dpx,0.25]);
        hold on; box on; grid on;
        title(model.data{ns}.name);
        n0 = (ns-1)*model.Nds;
        Ndp = model.Ndp;
        Be = [];
        
        for nw = 1:model.Nwalk
            I = find(model.walker{nw}.status > -1);
            Be = [Be(:,:);model.walker{nw}.gm(n0+1:Ndp+n0,I)'];
        end
        Beint = linspace(0.5*min(min(Be)),max(max(Be)),40)';
        xlabel(['NBe (atoms/g), sample ',num2str(ns)]);
        for j =1:Ndp
            histogram(Be(:,j));
            errorbar(model.data{ns}.N10(j),4000,model.data{ns}.dN10(j),'-or','horizontal','LineWidth',2);
        end
    end
end


%****** glaciation history *********

dpx = (1-2*lpm-ddm)/2;

%loop parameters
for i=1:2
    
    axes('position',[lpm+(i-1)*(ddm+dpx),0.37,dpx,0.25]);
    hold on; box on; grid on;
    
    if (i == 1), ylabel('Frequency'); end
    xlabel(model.mp{i}.name);
    set(gca,'xlim',[model.mp{i}.vmin, ...
        model.mp{i}.vmax]);
    set(gca,'yticklabel',[]);
    uval = [];
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status == 1);
        if (~isempty(I))
            [f,xi] = ksdensity(model.walker{nw}.up(i,I));
            line(xi,f,'color',wcol(nw,:));
            uval = [uval(:);model.walker{nw}.up(i,I)'];
        end
    end
    if (~isempty(uval))
        [f,xi] = ksdensity(uval);
    end
    line(xi,f,'color','k','linewidth',2);
    
    ylims = get(gca,'ylim');
    if isfield(model,'synpmval') == 1
        plot([model.synpmval(i),model.synpmval(i)],[ylims(1),ylims(2)],'r-','LineWidth',2); %plot synthetic pm values
    end
    
    med = median(uval);
    quart = quantile(uval,[0.25 0.75]);
    xline(med,'--b')
    xline(quart(1),'--r')
    xline(quart(2),'--r')
    
end



%********exhumation plots*********
Maxz = model.z0;
Maxt = model.age;
z0 = model.z0;
Nz = 100;
Nt = 200;
zint = linspace(0,Maxz,Nz);
tint = linspace(0,Maxt,Nt);
[tbin,zbin]=meshgrid(tint,zint);

dpx = (1-2*lpm-(model.Nsnr-1)*ddm)/model.Nsnr;

map = colormap;
map(1,:) = [1,1,1];
colormap(map);

for ns = 1:model.Nsnr
    
    %subplot(np(1),np(2),ns);
    axes('position',[lpm+(ns-1)*(ddm+dpx),0.05,dpx,0.25]);
    hold on; box on; set(gca,'layer','top');
    
    set(gca,'ydir','reverse');
    xlabel('Time (Ma)');
    if (ns == 1)
        ylabel('Burial depth');
    else
        set(gca,'yticklabel',[]);
    end
    
    histgrid = zeros(size(tbin));
    dzi = zint(2)-zint(1);
    title(model.data{ns}.name);
    
    n0 = model.Mmp + (ns-1)*model.Nsmp;
    
    for nw=1:model.Nwalk
        I = find(model.walker{nw}.status == 1);
        for i=1:length(I)
            
            switch numdp
                case 3
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
                    z4 = z3 + (model.age - T3)*E4; %hvis T3 bliver større end model.age kan der opstå problemer
                    if z4 < 0 %test if problems occur.
                        error('error: some models are not being properly exhumed')
                    end
                    
                    Tm = [0,T1,T2,T3,T4];
                    zm = [0,z1,z2,z3,z4];
                    
                    zinterp = interp1(Tm,zm,tint,'linear',2*z0);
                    izs = 1+floor(zinterp/dzi); %bin index at each itsfine
                    for itsfine=1:Nt
                        if (izs(itsfine)<=Nz && izs(itsfine)>0) %JLA added second clause 11.03.20(quickfix)
                            histgrid(izs(itsfine),itsfine)=histgrid(izs(itsfine),itsfine)+ 1;
                        end
                    end
                    
                case 4
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
                    if z5 < 0 %test if problems occur.
                        error('error: some models are not being properly exhumed')
                    end
                    
                    Tm = [0,T1,T2,T3,T4,T5];
                    zm = [0,z1,z2,z3,z4,z5];
                    
                    zinterp = interp1(Tm,zm,tint,'linear',2*z0);
                    izs = 1+floor(zinterp/dzi); %bin index at each itsfine
                    for itsfine=1:Nt
                        if (izs(itsfine)<=Nz && izs(itsfine)>0) %JLA added second clause 11.03.20(quickfix)
                            histgrid(izs(itsfine),itsfine)=histgrid(izs(itsfine),itsfine)+ 1;
                        end
                    end
            end
        end
        
    end
    
    [C,h] = contourf(tbin,zbin,(histgrid+0.5).^0.25,50, ...
        'linestyle','none');
    
    
end
if isfield(model,'synpmval') == 1
    syn_exhu(model.synpmval)
end
ylim([0 20]) %make variable depending on total erosion?

mname = strcat(reportpath,sampleID,'.pdf');

print(mname,'-dpdf','-fillpage');
append_pdfs(mname, 'temp2.pdf', 'temp3.pdf')
delete('temp2.pdf','temp3.pdf')

set(groot','defaulttextinterpreter','default');
set(groot, 'defaultAxesTickLabelInterpreter','default');
set(groot, 'defaultLegendInterpreter','default');

%set(gcf, 'Position', get(0, 'Screensize'));
filename = strcat(sname,'report.jpg');
saveas(gcf,filename)

end