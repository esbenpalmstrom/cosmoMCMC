function analyse_walker(snr,nw)

close all;

str = num2str(snr(1));
if length(snr) > 1,
    for i=2:length(snr),
        str = [str,['-',num2str(snr(i))]];
    end;
end;
%mname = ['models/NEGreenland/Koldewey/MC_sample_',str,'.mat'];
% mname = ['models/NEGreenland/Hellefjorden/MC_sample_',str,'.mat'];
mname = ['models/MDML/MC_sample_',str,'.mat'];
load(mname);

figure()
set(gcf,'units','normalized','position',[.1,.3,.4,.6]);
set(gca,'position',[0,0,1,1],'visible','off');
set(gca,'xlim',[0,1],'ylim',[0,1]);
set(gcf,'Name','Walker information');
ax0 = gca;


ax1 = axes('position',[.1,.55,.8,.4]);
hold on; box on;
set(gca,'yscale','log');
I = find(model.walker{nw}.status > -2);
plot(model.walker{nw}.restot(I),'.k');
%plot(model.walker{nw}.kstep(I),model.walker{nw}.restot(I),'.k');
xlim = get(gca,'xlim');
line(xlim,[1,1],'color','r');
ylabel('Residual');

ax2 = axes('position',[.1,.1,.8,.4]);
hold on; box on;
plot(model.walker{nw}.accrate,'.r');
plot(model.walker{nw}.kstep,'.k');
ylabel('Acceptance ratio and kstep');

