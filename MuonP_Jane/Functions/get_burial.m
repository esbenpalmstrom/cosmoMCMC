function [dO_i,bt_i] = get_burial()

%load and prepare d18O data
load('d18Ocurves.mat');
I = find(Age < 1e6);
dO = d18O_4ky(I);
dOt = Age(I);
dt = diff(dOt);
dO = dO(1:(end-1));

dO_i = linspace(3.5,5.0,200);
bt_i = zeros(size(dO_i));

for i=1:length(dO_i),
    
    I = find(dO > dO_i(i));
    
    bt_i(i) = sum(dt(I));

end;



