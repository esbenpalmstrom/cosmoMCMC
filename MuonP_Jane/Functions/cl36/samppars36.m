%
% sp=samppars36(sampledata)
%
% Extracts sample parameters from a sample vector and puts them
% into sp.
%
function sp=samppars36(sampledata)
%
% Make sure that sampledata is a vector rather than an array!
%
if (min([size(sampledata,1); size(sampledata,2)]) > 1)
  error('sampledata must be a vector!');
end
%
% Make sampledata a column vectors if it isn't already.
%
if (size(sampledata,1)==1)
  sampledata=sampledata';
end
%
% First, check that the input data is reasonable.
%
if (length(sampledata) ~= 38)
  error('sampledata has wrong size!');
end
%
% Setup the values of sample parameters.
%
sp.concentration36=sampledata(1);
sp.inheritance36=sampledata(2);
sp.epsilon=sampledata(3);
sp.qavg=sampledata(4);
sp.rb=sampledata(5);

%calculate the normalization factor for the oxides after converting the
%fractional water content to weight % water
sp.oxidenormfactor=1-(sp.qavg/sp.rb);

sp.ls=sampledata(6);
sp.latitude=sampledata(7);
sp.longitude=sampledata(8);
sp.ST=sampledata(11);
sp.Lambdafe=sampledata(12);
sp.originaloxideinput=sampledata(13:23);

%renormalize the original oxide info to account for water content of rock
newoxidepercent=sp.originaloxideinput*sp.oxidenormfactor;
sp.ci=[NaN; NaN; newoxidepercent;  sampledata(24:31)];

%no renormalization for the target element composition
sp.tci=sampledata(32:36);

%
%Note: sp.ci and sp.tci are for chlorine only and do not have an indicator
%after them for this reason. 
%
%
% Extract the pressure and elevation from the sampledata vector.
%
sp.elevation=sampledata(9);
sp.P=sampledata(10);
%
% Extra the depth to top of sample.
%
sp.depthtotop=sampledata(37);
%
% The year the sample was collected determines the final point in
% time for integration of nuclide production.  This also depends on
% the time axis of the geomagnetic history being used.  Currently,
% "t=0" means "2010."  However, this will eventually change when a
% new geomagnetic history is provided.  Here, we computed tfinal
% based on the current geomag history, which ends in 2010.
%
sp.tfinal=(sampledata(38)-2010)/1000;
