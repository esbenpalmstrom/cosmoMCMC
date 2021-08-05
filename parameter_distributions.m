%file: parameter_distributions.m

%{
Minimum residual parameters for each walker in each base inversion
output parametermeans.mat is a mean of minimum residual model parameters
from all walkers. 100 per walker.

remember to uncomment/comment the save command at bottom of script.

WIP notes:
Add plots of exhumation based on real data inversion where the different
proposed synthetic models are displayed. Mark out the chosen model.
Consider displaying only for one walker, or for all walkers.
%}
clear; close all;
%set to 'on' if you want
exhumationplots = 'on';
frequencyplots = 'on';
set(0, 'DefaultTextInterpreter', 'none')



load fsamples.mat
addpath Functions
cd models/FSbase_12_Dec_2020_13_04_14/ %pick models!!!
files = dir('*.mat');
nsamples = 6; %number of samples in files

nsamp = 100; %number of points to use for mean

walksamp = [
    4
    1
    3
    4
    1
    2
    ]; %define which walker to use for sample parameters

exhumfig = figure();
numsubplot = 1;

sitenames = {'Pyha-Nättanen','lol','Gaustatoppen','Andøya','Naakakarhakka','Karmøy'};


for ff = 1:nsamples
    load (files(ff).name);
    
    if ff == 2
        continue
    end
    
    figure;
    
    for j = 1:model.Nwalk
        [sortedValues,sortedIndex] = sort(model.walker{j}.restot);
        min_res_index{j} = sortedIndex(1:nsamp);
        min_res{j} = model.walker{j}.restot(min_res_index{j});
        min_u{j} = model.walker{j}.up(:,min_res_index{j});
        wlk_mean(:,j) = mean(min_u{j},2);
    end
    
    walker_totmean = mean(wlk_mean,2);
    
    
    
    
    %plot of parameter distributions vs frequency with sampled parameters
    %included
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
    
    
    set(gcf,'units','normalized','position',[.3,.2,.4,.6]);
    set(gcf,'Name','Parameter distributions with sampled parameters');
    [np,nn] = numSubplots(model.Nmp);
    
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
                wlkline{nw} = line(xi,f,'color',wcol(nw,:));
                uval = [uval(:);model.walker{nw}.up(i,I)'];
                
                fi = interp1(xi,f,min_u{nw}(i,:));
                plot(min_u{nw}(i,:),fi,'o','color',wcol(nw,:));
                
                [ffi,xii] = ksdensity(min_u{nw}(i,:));
                yyaxis right
                line(xii,ffi,'color',wcol(nw,:),'LineStyle','--');
                yyaxis left
                
                if nw == walksamp(ff) && i == 1 %i==1 if model picked based only on d18Ot
                    [maxF,maxI] = max(fi);
                    parameterMeans.(model.data{1}.name) = min_u{walksamp(ff)}(:,maxI);
                    %plot(min_u{walksamp(ff)}(i,maxI),fi(maxI),'ro')
                end
                if i == 1 %save best models from all walkers in separate variable
                    [maxsF,maxsI] = max(fi);
                    walkerModels(:,nw) = min_u{nw}(:,maxsI);
                end
                
                if nw == walksamp(ff)
                    plot(min_u{walksamp(ff)}(i,maxI),fi(maxI),'ro','MarkerSize',10,'LineWidth',3)
                end
                
            end
            
        end
        
        lgd = legend([wlkline{1} wlkline{2} wlkline{3} wlkline{4}],{'wlk 1','wlk 2','wlk 3','wlk 4'});
        
        if (~isempty(uval))
            [f,xi] = ksdensity(uval);
            
        end
        line(xi,f,'color','k','linewidth',2);
        fg = interp1(xi,f,walker_totmean(i));
        plot(walker_totmean(i),fg,'ko','MarkerSize',15,'LineWidth',5)
        ylims = get(gca,'ylim');
        
        keyboard
        if isfield(model,'synpmval') == 1
            plot([model.synpmval(i),model.synpmval(i)],[ylims(1),ylims(2)],'r-','LineWidth',2); %plot synthetic pm values
        end
        
    end
    sgtitle(files(ff).name);
    parameterMeans.IDs{ff} = model.data{1}.name;
    for mpn = 1:length(model.mp)
        parameterMeans.mpname{mpn} = model.mp{1,mpn}.name;
    end
    
    
    
    %make exhumation plot with proposed synthetic models displayed.
    %figure()
    set(0,'CurrentFigure',exhumfig);
    
    %plot margins
    Maxz = model.z0;
    Maxt = model.age;
    z0 = model.z0;
    Nz = 100;
    Nt = 200;
    zint = linspace(0,Maxz,Nz);
    tint = linspace(0,Maxt,Nt);
    [tbin,zbin]=meshgrid(tint,zint);
    map = colormap;
    map(1,:) = [1,1,1];
    colormap(map);
    for ns = 1:model.Nsnr
        sbplt = subplot(3,2,numsubplot);
        
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
        %title(model.data{ns}.site);
        %title(sitenames{ff});
        
        n0 = model.Mmp + (ns-1)*model.Nsmp;
        
        for nw=1:model.Nwalk
            I = find(model.walker{nw}.status == 1);
            for i=1:length(I)
                
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
                izs = 1+floor(zinterp/dzi); %bin index at each tsfine
                for itsfine=1:Nt
                    if (izs(itsfine)<=Nz && izs(itsfine)>0) %JLA added second clause 11.03.20(quickfix)
                        histgrid(izs(itsfine),itsfine)=histgrid(izs(itsfine),itsfine)+ 1;
                    end
                end
            end
        end
        
        [C,h] = contourf(tbin,zbin,(histgrid+0.5).^0.25,50, ...
            'linestyle','none');
        
        %plot best models for each walker.
        for nw = 1:model.Nwalk
            
            %             for k = 1:nsamp %plot a cloud of "good" models
            %                 up = min_u{nw}(:,k);
            %                 syn_exhu(up);
            %             end
            
            %             up = walkerModels(:,nw);
            %             syn_exhu(up); %plot the best model from each walker.
            
            
        end
        
        up = parameterMeans.(parameterMeans.IDs{ff});
        syn_exhu(up,[1,0,0],[.6,.9,.6]); %plot the "chosen one", make it red
        
        ylim([0 20])
    end
    keyboard
    
    
    
    
    
    %draw box for d18O curve
    sbplt.Position(4) = sbplt.Position(4)-0.05;
    
    ax = axes;
    hold on; box on;
    ax.Position = sbplt.Position;
    ax.Position(2) = sbplt.Position(2)+sbplt.Position(4)+0.01;
    ax.Position(4) = 0.5*sbplt.Position(4);
    
    ax.XTick = sbplt.XTick;
    ax.XTickLabel = '';
    
    %*******load and draw the d18O curve*******
    cg = [.6,.6,.9];
    cig = [.6,.9,.6];

    load d18Ocurves.mat
    tta = Age*1e-3;
    dd = d18O_10ky;
    uval = [];
    for nw = 1:model.Nwalk
        I = find(model.walker{nw}.status == 1);
        if (~isempty(I))
            uval = [uval(:);model.walker{nw}.up(1,I)'];
        end
    end
    
    if (~isempty(uval))
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
    
    l2 = line(tta,dd,'color','k');
    keyboard
    fval = 0.9*f/max(f);
    patch([fval(:);zeros(size(fval(:)))],[xi(:);flipud(xi(:))],cg,'linestyle','none');
    line(fval,xi,'color','k');

    
    
    
    title(sitenames{ff});
    
    
    
    
    
    
    numsubplot = numsubplot +1;
    
    
end
set(gcf,'color','white')
set(gcf,'Position',[100 100 500 750]);

figurename = 'SyntheticModels';
export_fig(['/Users/esben/OneDrive - Aarhus Universitet/Speciale/skriv/latex/figures/' figurename],'-jpg','-r300',exhumfig);


cd ../..
%save('parametermeans.mat','parameterMeans')