%file: makeStatistics_vE1.m
%{
Purpose of script is to make figures diplaying fx. erosion, glaciation and d18O
in various plots for all inversions that have been run.

WIP:
plot a frequency band as in Skov et al (2020)
Boxplot
Violinplot?

make sure the script is expandable - so that any additional models can be
added easily

youll have to add a special case for the walker variations, since all the
.mat files will be in the same folder in that case.

Add subplots where you can display d18O threshold and maybe other erosion
timings fx. erosion at 200k years.

Consider calculating cosmogenic nuclide memory, and then using that memory
limit as timing to calculate total erosion. Wait until having a meeting
with MFK regarding this?
%}
close all;

load fsamples.mat

%path variables for the different inversions. Remember to check if youre
%using the latest model!
path1 = 'models/FSbase_12_Dec_2020_13_04_14/';
path2 = 'models/FS_Be0/13_Dec_2020_17_38_44/';
path3 = 'models/FS_BeAl0/15_Dec_2020_14_54_48/';
path4 = 'models/FS_BeAl0_0.5_1/12_Dec_2020_23_53_04/';
path5 = 'models/FS_BeAl0_1_2/15_Dec_2020_22_25_00/';
%path5 = 'models/FS_BeAl0_1_2/15_Dec_2020_22_25_00/walkervariations/';
paths = {path1, path2, path3, path4, path5};

ypos = linspace(0.05,0.95,length(paths));

for i = 1:length(fsamples.IDs)
    
    figure()
    set(gcf,'Name','model frequencies');
    
    for j = 1:length(paths)
        
        subplot(1,2,1)
        hold on; box on; grid on;
        set(gca,'YTickLabel',[]);
        xlim([0 20])
        ylim([0 1])
        
        path = [paths{j} fsamples.IDs{i} '.mat'];
        load(path);
        
        
        uval = [];
        for nw = 1:model.Nwalk
            I = find(model.walker{nw}.status == 1); %use only accepted models.
            if (~isempty(I))
                uval = [uval(:,:) model.walker{nw}.up(:,I)];
            end
        end
        
        d18OT = uval(1,:);
        T1 = uval(2,:);
        z1 = uval(3,:);
        dT2 = uval(4,:);
        dz2 = 10.^uval(5,:);
        dT3 = uval(6,:);
        dz3 = 10.^uval(7,:);
        E4 = 10.^uval(8,:);
        
        T2 = T1 + dT2;
        z2 = z1 + dz2;
        T3 = T2 + dT3;
        z3 = z2 + dz3;
        T4 = ones(1,length(T3)).*model.age;
        z4 = z3 + (model.age - T3).*E4;
        
        
        %erosion at 1 ma
%         eroeins = zeros(1,length(T4));
%         for nnmod = 1:length(uval)
%             %interpolate erosion at 1 ma. for all models.
%             eroeins(nnmod) = interp1([0,T1(nnmod),T2(nnmod),T3(nnmod),T4(nnmod)],[0,z1(nnmod),z2(nnmod),z3(nnmod),z4(nnmod)],1);
%             
% %             if eroeins(nnmod) > 20 %anything more than 20 m erosion is gonna be outside the xlim anyways
% %                 eroeins(nnmod) = NaN;
% %             end
%         end
        
        
        load([paths{j} 'param_stats/' fsamples.IDs{i} 'Syn_stats.mat']);
        
        
        xpos = linspace(0,20,100); %frequency boxes and position of patches
        f = histc(eroeins,xpos); %count erosion into boxes
        normf = f./sum(f); %normalize to get frequency
        
        ywidth = 0.03; %width of frequency bars
        
        %make boxes with frequency color
        for nbox = 1:length(f)
            if normf(nbox) ~= 0
                x = [xpos(nbox) xpos(nbox+1) xpos(nbox+1) xpos(nbox)];
                y = [ypos(j)-ywidth ypos(j)-ywidth ypos(j)+ywidth ypos(j)+ywidth];
                patch(x,y,normf(nbox),'EdgeColor','none')
            end
        end
        
        colorbar
        
        textoffset = 3;

        if model.Nnc == 2 && ~isfield(model,'synpmval')
            text((xpos(1)-textoffset),ypos(j),['Real BeAl' num2str(model.data{1}.depths','% 5.1f')])
        elseif model.Nnc == 1 && ~isfield(model,'synpmval')
            text((xpos(1)-textoffset),ypos(j),['Real Be' num2str(model.data{1}.depths','% 5.1f')])
        end
        
        if  isfield(model,'synpmval') %add indicators at the synthetic model used.
            
            synT1 = model.synpmval(2);
            synz1 = model.synpmval(3);
            synT2 = model.synpmval(2) + model.synpmval(4);
            synz2 = model.synpmval(3) + 10.^model.synpmval(5);
            synT3 = synT2 + model.synpmval(6);
            synz3 = synz2 + 10.^model.synpmval(7);
            synT4 = ones(1,length(synT3)).*model.age;
            synz4 = synz3 + (model.age - synT3).*10.^model.synpmval(8);
            
            synero = interp1([0,synT1,synT2,synT3,synT4],[0,synz1,synz2,synz3,synz4],1);
            line([synero synero],[(ypos(j)-ywidth-0.02) (ypos(j)+ywidth+0.02)],'Color','red','LineWidth',1.5)
            
            if model.Nnc == 2
                text((xpos(1)-textoffset),ypos(j),['BeAl' num2str(model.data{1}.depths','% 5.1f')])
            elseif model.Nnc == 1
                text((xpos(1)-textoffset),ypos(j),['Be' num2str(model.data{1}.depths','% 5.1f')])
            end
            
        end
        
        
        
        %make plot of d18Ot distribution
        subplot(1,2,2)
        set(gca,'YTickLabel',[]);
        
        xpos = linspace(3.5,5,100); %frequency boxes and position of patches
        f = histc(d18OT,xpos); %count erosion into boxes
        normf = f./sum(f); %normalize to get frequency
        
        for nbox = 1:length(f)
            if normf(nbox) ~= 0
                x = [xpos(nbox) xpos(nbox+1) xpos(nbox+1) xpos(nbox)];
                y = [ypos(j)-ywidth ypos(j)-ywidth ypos(j)+ywidth ypos(j)+ywidth];
                patch(x,y,normf(nbox),'EdgeColor','none')
            end
        end
        
        
        if isfield(model,'synpmval')
            synd18 = model.synpmval(1);
            line([synd18 synd18],[(ypos(j)-ywidth-0.02) (ypos(j)+ywidth+0.02)],'Color','red','LineWidth',1.5)
        end
        
        
        
    end
    sgtitle((model.data{1}.name))

    
end
distFig