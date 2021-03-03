% Calculate and plot current 36Cl production depth profiles for a sample pf
% composition specified in input sampledata vector. This function and
% subfunctions are based on CRONUSCalc m-files downloaded from 
% https://bitbucket.org/cronusearth/cronus-calc in March 2019 (last update
% of bitbucket code was August 2018 at this point). JLA

%  The sampledata vector contains the following information.
%
%1.     Sample 36-Cl concentration (atoms of 36-Cl/g of target)
%2.     Inheritance (atoms 36-Cl/g of target)  
%3.     erosion-rate epsilon (g/(cm^2*kyr))
%4.     fractional volumetric water-content (unitless) 
%5.     bulk density (g/cm^3)
%6.     sample thickness (cm)
%7.     Latitude (decimal degrees, -90(S) to +90(N))
%8.     Longitude (decimal degrees, 0-360 degrees east)
%9.     Elevation (meters)
%10.    Pressure (hPa)                Both 9 and 10 must be present!
%11.    Shielding factor for terrain, snow, etc. (unitless)
%12.    Effective attenuation length -Lambdafe (g/cm^2)
%13.    % CO2                        Rock
%14.    % Na2O                       Rock
%15.    % MgO                        Rock
%16.    % Al2O3                      Rock
%17.    % SiO2                       Rock
%18.    % P2O5                       Rock
%19.    % K2O                        Rock
%20.    % CaO                        Rock
%21.    % TiO2                       Rock
%22.    % MnO                        Rock
%23.    % Fe2O3                      Rock
%24.    Cl (ppm)                     Rock
%25.    B (ppm)                      Rock
%26.    Sm (ppm)                     Rock
%27.    Gd (ppm)                     Rock
%28.    U (ppm)                      Rock
%29.    Th (ppm)                     Rock
%30.    Cr (ppm)                     Rock
%31.    Li (ppm)                     Rock
%32.	Target element %K2O          Target
%33.    Target element %CaO          Target
%34.    Target element %TiO2         Target
%35.    Target element %Fe2O3        Target
%36.    Target element Cl (ppm)      Target
%37.    Depth to top of sample (g/cm^2)
%38.    Year sampled (e.g. 2010)
%
% A second input vector, sampleuncertainties, contains 1-sigma
% uncertainties for all 36 inputs.  In general, we assume that
% these 36 inputs are uncorrelated.  The one exception is that 
% we allow for correlation between input parameter 1 (36-Cl
% concentration in the target rock in atoms/g) and parameter 36
% (Concentration of Cl in the target rock in ppm.)  The covariance
% between these paramters is given as the third input parameter, covar.
%
% scaling_model is one of 'DE','DU','LI','LM','SA','SF','ST' and 
% informs which scaling model is being used


function cl36_JLA2(bulkcomp,K,Ca,Cl,scaling_model)

  close all
  addpath scaling
  
  load Basen.mat
  sampledata=data36(1,:);
  sampledata(1,36)=Cl;sampledata(1,20)=Cl; %Update target and bulk rock Cl composition
  sampledata(1,32)=K;sampledata(1,33)=Ca; %Update target Ca and K composition - note this approximation does not change bulk cross section values!
  switch bulkcomp 
      case 'basalt'
      case 'granite'
  end
  % Make sampledata and uncertainties column vectors if they aren't already.
  
  if (size(sampledata,1)==1)
    sampledata=sampledata';
  end
%   if (size(sampleuncertainties,1)==1)
%     sampleuncertainties=sampleuncertainties';
%   end
  
  % First, check that the input data is reasonable.
  if (length(sampledata) ~= 38)
    error('sampledata has wrong size!');
  end
%   if (length(sampleuncertainties) ~= 38)
%     error('sampleuncertainties has wrong size!');
%   end
 
%   % Give a warning about the composition not adding up to 100%.
%   if (abs(sum(sampledata(13:23))-100) > 2.0)
%     warning('major element composition does not add up to one hundred percent!');
%   end

  % Setup the physical parameters.
  pp=physpars(scaling_model);

  % Extract the sample parameters from the sampledatavector.
  sp=samppars36(sampledata);
  
  % Get the scale factors.
  sf=scalefacs36(sp,scaling_model);
  
  % We need an absolute maximum age for several purposes, including
  % detecting saturated samples and setting the maximum depth for comppars.
  maxage=2000;               % 2Ma > 6 half lives              
  
  % Figure out the maximum possible depth at which we'll ever need a
  % production rate.  This is depthtotop + maxage * erosion +
  % thickness * density + a safety factor.
  
  maxdepth=sp.depthtotop+maxage*sp.epsilon+sp.ls*sp.rb+1000;
  %check for maximum depth: muon formulation is only good down to 2e5 g/cm2.
  % This will likely never happen, but if it does, Matlab will error out saying we have not
  % supplied times and plotprod and we will receive that error and be able to help the user.
  if maxdepth > 2e5
    fprintf(1,'Maximum sample depth (%f) exceeds muon formulation of 2e5 g/cm2. \n Options: Lower the erosion rate, lower the maxage in file cl36age, or change muon formulation in muonfluxsato.m',[maxdepth])
    warning('This sample exceeds muon maximum depth. Try lowering the erosion rate of the sample or decreasing sample depth.');
% %     output=NaN*ones(expectedOutputs,1);
%     [~, plotprodca, plotprodk, plotprodcl ] = getPlotData(maxage, sf, 0, ...
%       0, 0, 0, 0, 0, 0, 'cl', scaling_model);
    return;
  end
  
  % Computed parameters.
  cp=comppars36(pp,sp,sf,maxdepth);
  
  % Get contemporary depth production rates in atoms/g 
  sf.currentsf=getcurrentsf(sf,0,scaling_model,'cl');
    D_m =3.33; %Max depth of profile (cm)
    z_m = linspace(0,10,100);
    z_D = D_m*z_m.^3/10*sp.rb; %denser depth-grid near surface
%     z=linspace(0,333,100); z=z*sp.rb; %Linear depth spacing
  [Prodtotal,Prods,ProdsCa,ProdsK,ProdsTi,ProdsFe,Prodth,Prodeth,Prodmu,...
      ~,~,~,~,Kpercent,Capercent,Clpercent]...
      =prodz36(z_D,pp,sf,cp);

%   Plot production depth profiles
  figure(),
  plot(Prodtotal,z_D,'k','Linewidth',2), hold on
  plot(Prods,z_D,'.-'),
%   plot(ProdsCa,z_D,ProdsK,z,ProdsTi,z,ProdsFe,z),
  plot(Prodth,z_D,Prodeth,z_D,'.-'),
  plot(Prodmu,z_D,'--'),
  legend('Total production','Spallation','Thermal','Epithermal','Muon','Location','southeast')
%   legend('Total production','Spallation','Ca spal','K spal','Ti spal',...
%       'Fe spal','Thermal','Epithermal','Muon','Location','southeast')
  xlabel('Production [at/g/yr]'),ylabel('Depth [g/cm^2]')
  set(gca,'ydir','reverse','xscale','log')
  set(gcf,'DefaultTextInterpreter','Latex');
  text(1.5e-2,50,['Ca: ' num2str(sampledata(33)) ' $\%$, K: ' ...
      num2str(sampledata(32)) ' $\%$, Cl: ' num2str(sampledata(36)) ' ppm'])
  
  figure(),
  plot(Capercent,z_D,Kpercent,z_D,Clpercent,z_D,'Linewidth',2)
  legend('Ca','K','Cl')
  xlabel('Production [%]'),ylabel('Depth [g/cm^2]')
  set(gca,'ydir','reverse')
  set(gcf,'DefaultTextInterpreter','Latex');
  text(30,50,['Ca: ' num2str(sampledata(33)) ' $\%$, K: ' ...
      num2str(sampledata(32)) ' $\%$, Cl: ' num2str(sampledata(36)) ' ppm'])