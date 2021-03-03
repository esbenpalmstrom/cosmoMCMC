function model_fig(doc,zp,tp)

set(groot','defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');

% ***** figure props ********
% close all;
set(gcf,'units','pixels','position',[1000,1000,800,800]);
set(gcf,'color',[1,1,1]);
set(gca,'position',[0,0,1,1],'visible','off');
ax0 = gca;
set(ax0,'xlim',[0,1],'ylim',[0,1]);
text(0.5,0.035,'Time before present (kyr)','horizontalalignment','center','fontsize',24);
text(0.05,0.825,'$\delta_{18}$O','rotation',90,'horizontalalignment','center','fontsize',24);
text(0.04,0.35,'Depth below surface (m)','rotation',90,'horizontalalignment','center','verticalalignment','middle','fontsize',24);
%text(0.075,0.1,'20','horizontalalignment','right','verticalalignment','middle','fontsize',24);
%text(0.075,0.6,'0','horizontalalignment','right','verticalalignment','middle','fontsize',24);

ax1 = axes('position',[0.1,0.675,0.75,0.28]);
hold on; box on; %grid on;
set(gca,'ydir','reverse','layer','top');
set(gca,'xtick',[0:1000:3000],'xticklabel',[]);
set(gca,'ytick',[3:5]);
set(gca,'fontsize',24);
axis([0,3000,3,5.3]);

ax2 = axes('position',[0.1,0.1,0.75,0.5]);
hold on; box on; %grid on;
set(gca,'ydir','reverse');
set(gca,'xtick',[0:1000:3000]);
set(gca,'ytick',[0:5:20]);
set(gca,'fontsize',24);
axis([0,3000,0,20]);

ax3 = axes('position',[0.20,0.20,0.25,0.2]);
hold on; box on; %grid on;
set(gca,'ydir','reverse');
set(gca,'xtick',[0:10:50]);
set(gca,'ytick',[0:0.2:0.5]);
set(gca,'fontsize',20);
axis([0,25,0,0.5]);

%colors
cg = [.6,.6,.9];
cig = [.6,.9,.6];

axes(ax1);
load d18Ocurves.mat
tta = Age*1e-3;
dd = d18O_6ky;

%critical d18O
% doc = 3.9;

%patch glacial/interglacial domains
ddg = dd;
ddig = dd;
I = find(dd > doc); ddg(I) = doc;
I = find(dd < doc); ddig(I) = doc;
p1 = patch([tta(:);flipud(tta(:))],[dd(:);flipud(ddg(:))],cg,'linestyle','none');
p2 = patch([tta(:);flipud(tta(:))],[dd(:);flipud(ddig(:))],cig,'linestyle','none');

xlim = get(gca,'xlim');
l1 = line(xlim,doc*[1,1],'color','k','linestyle','--');

%draw d18O curve
l2 = line(tta,dd,'color','k');

%time of deglaciation
Imin = find(dd > doc);
% Tdg = tta(Imin(1));

pos = get(gca,'position');
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');
xOt = pos(1)+pos(3)+0.015;
yOt = pos(2)+pos(4)*(ylim(2)-doc)/(ylim(2)-ylim(1));

axes(ax2); cla;
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');

nnx = xlim(2)-xlim(1);
nny = ylim(2)-ylim(1);

%depth points
% zp = [0.0,0.2,3.5,6.0,25];
% tp = [0.0,Tdg,150,1500,2500];
zend = zp(end) + (xlim(2)-tp(end))*(zp(end)-zp(end-1))/(tp(end)-tp(end-1));

x0 = 0; dx = 100;
y0 = 0; dy = 1;
line([x0,x0+dx,x0+dx,x0,x0],[y0,y0,y0+dy,y0+dy,y0],'color','k','linewidth',2);


%draw line segment
line(tp,zp,'color',cg,'linewidth',4);
        
%draw dashed lines    
for j=3:length(tp)-1,
    lh(j) = line([tp(j),tp(j)+(xlim(2)-tp(j))],[zp(j),zp(j)],'linestyle','--','color','k');
    lv(j) = line([tp(j),tp(j)],zp(j)*[1,0],'linestyle','--','color','k');
end;

%draw points
mz = 15;
plot(tp(2:(end-1)),zp(2:(end-1)),'o','markeredgecolor','k','markersize',mz,'markerfacecolor',cig);

pos = get(gca,'position');
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');


for j=3:length(tp)-1,
    xtp(j-1) = pos(1) + pos(3)*(tp(j)-xlim(1))/(xlim(2)-xlim(1));
    ytp(j-1) = pos(2) + pos(4) + 0.02;
    tlabel{j-1} = ['$T_',num2str(j-1),'$'];
    xzp(j-1) = pos(1) + pos(3) + 0.02;
    yzp(j-1) = pos(2) + pos(4)*(ylim(2)-zp(j))/(ylim(2)-ylim(1));
    zlabel{j-1} = ['$Z_',num2str(j-1),'$'];
end;

axes(ax3); hold on;


%draw dashed lines    
for j=2,
    lh(j) = line([tp(j),tp(j)+(xlim(2)-tp(j))],[zp(j),zp(j)],'linestyle','--','color','k');
    lv(j) = line([tp(j),tp(j)],zp(j)*[1,0],'linestyle','--','color','k');
end;

%draw line segment
line(tp(1:3),zp(1:3),'color',cg,'linewidth',4);

%draw points
mz = 15;
plot(tp(1:2),zp(1:2),'o','markeredgecolor','k','markersize',mz,'markerfacecolor',cig);

pos = get(gca,'position');
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');

for j=2,
    xtp(j-1) = pos(1) + pos(3)*(tp(j)-xlim(1))/(xlim(2)-xlim(1));
    ytp(j-1) = pos(2) + pos(4) + 0.03;
    tlabel{j-1} = ['$T_{\textrm{dg}}$'];
    xzp(j-1) = pos(1) + pos(3) + 0.02;
    yzp(j-1) = pos(2) + pos(4)*(ylim(2)-zp(j))/(ylim(2)-ylim(1));
    zlabel{j-1} = ['$Z_{\textrm{dg}}$'];
end;

axes(ax0); hold on;
for j=2:length(tp)-1,
    ttp(j-1)=text(xtp(j-1),ytp(j-1),tlabel{j-1},'fontsize',24,'horizontalalignment','center');
    ttz(j-1)=text(xzp(j-1),yzp(j-1),zlabel{j-1},'fontsize',24,'verticalalignment','middle');
end;
ttO=text(xOt,yOt,'$\delta_{18}$O$_c$','fontsize',24,'verticalalignment','middle');

text(0.325,0.15,'Time before present (kyr)','horizontalalignment','center','fontsize',20);
text(0.14,0.3,'Depth below surface (m)','rotation',90,'horizontalalignment','center','verticalalignment','middle','fontsize',20);
text(0.6,0.25,'$E_{0}$','fontsize',24,'verticalalignment','middle');


print -depsc model_setup.eps;


set(groot','defaulttextinterpreter','default');
set(groot, 'defaultAxesTickLabelInterpreter','default'); 
set(groot, 'defaultLegendInterpreter','default');