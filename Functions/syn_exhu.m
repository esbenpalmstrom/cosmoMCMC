function syn_exhu(up,cg,cig)
%{
modified by ELP, dec. 2020
Added option to define colors in function input. Default parameters will be
used if only 'up' parameter is passed to the function
%}

%colors
if nargin == 1
    cg = [.6,.6,.9];
    cig = [.6,.9,.6];
end
xlim = get(gca,'xlim');

%depth points
zp=[0.0,up(3),up(3)+10^(up(5)),up(3)+10^(up(5))+10^(up(7)),up(3)+10^(up(5))+10^(up(7))+(3 - up(2)+up(4)+up(6))*10^(up(8))];
tp=[0.0,up(2),(up(2)+up(4)),(up(2)+up(4)+up(6)),3];
% zend = zp(end) + (xlim(2)-tp(end))*(zp(end)-zp(end-1))/(tp(end)-tp(end-1));

%draw line segment
line(tp,zp,'color',cg,'linewidth',4);
        
%draw points
mz = 15;
plot(tp(2:(end-1)),zp(2:(end-1)),'o','markeredgecolor','k','markersize',mz,'markerfacecolor',cig);

pos = get(gca,'position');
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');


for j=3:length(tp)-1
    xtp(j-1) = pos(1) + pos(3)*(tp(j)-xlim(1))/(xlim(2)-xlim(1));
    ytp(j-1) = pos(2) + pos(4) + 0.02;
    tlabel{j-1} = ['$T_',num2str(j-1),'$'];
    xzp(j-1) = pos(1) + pos(3) + 0.02;
    yzp(j-1) = pos(2) + pos(4)*(ylim(2)-zp(j))/(ylim(2)-ylim(1));
    zlabel{j-1} = ['$Z_',num2str(j-1),'$'];
end