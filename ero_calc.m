%{
erosion calculator.
short script to calculate total erosion frequency at different points in
time from the models from MCMC inversion.
Save erosion to mat file
%}

load fsamples.mat

path1 = 'models/FSbase_12_Dec_2020_13_04_14/';
path2 = 'models/FS_Be0/13_Dec_2020_17_38_44/';
path3 = 'models/FS_BeAl0/15_Dec_2020_14_54_48/';
path4 = 'models/FS_BeAl0_0.5_1/12_Dec_2020_23_53_04/';
path5 = 'models/FS_BeAl0_1_2/15_Dec_2020_22_25_00/';


paths = {path1, path2, path3, path4, path5};

ypos = linspace(0.05,0.95,length(paths));

for i = 1:length(fsamples.IDs)
    
    for j = 1:length(paths)
        
        
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
        eroeins = zeros(1,length(T4));
        for nnmod = 1:length(uval)
            %interpolate erosion at 1 ma. for all models.
            eroeins(nnmod) = interp1([0,T1(nnmod),T2(nnmod),T3(nnmod),T4(nnmod)],[0,z1(nnmod),z2(nnmod),z3(nnmod),z4(nnmod)],1);
            
%             if eroeins(nnmod) > 20 %anything more than 20 m erosion is gonna be outside the xlim anyways
%                 eroeins(nnmod) = NaN;
%             end
        end
        
        mkdir(paths{j},'param_stats')
        
        savepath = [paths{j} 'param_stats/' fsamples.IDs{i} 'Syn_stats.mat'];
        save(savepath,'eroeins')
    end
end