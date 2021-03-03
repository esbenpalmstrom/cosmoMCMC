clear;

load parametermeans.mat

%names
fsamples.IDs{1} = 'Sven2015_12';
fsamples.IDs{2} = 'Ling2006a_11';
fsamples.IDs{3} = 'Nesj2007_9';
fsamples.IDs{4} = 'Hall2019_12';
fsamples.IDs{5} = 'Stroe2002a_1';
fsamples.IDs{6} = 'Darm2008_1';

%Sven2015_12
fsamples.Sven2015_12.lat = 59.1843;
fsamples.Sven2015_12.lon = 5.195;
fsamples.Sven2015_12.N10 = 442731;
fsamples.Sven2015_12.N10unc = 8227;
fsamples.Sven2015_12.nnc = 1; %number of nuclides
fsamples.Sven2015_12.P10total = 4.588;
fsamples.Sven2015_12.P10spals = 4.51;
fsamples.Sven2015_12.P10muon = 0.078;
fsamples.Sven2015_12.P26spals = fsamples.Sven2015_12.P10spals*6.75;
fsamples.Sven2015_12.P26total = fsamples.Sven2015_12.P26spals*1.0204;
fsamples.Sven2015_12.P26muon = fsamples.Sven2015_12.P26total*0.02;


%Ling2006a_11
fsamples.Ling2006a_11.lat = 59.848157999999998;
fsamples.Ling2006a_11.lon = 8.660021000000000;
fsamples.Ling2006a_11.N10 = 1745206;
fsamples.Ling2006a_11.N10unc = 102780;
fsamples.Ling2006a_11.nnc = 1;
fsamples.Ling2006a_11.P10total = 19.336;
fsamples.Ling2006a_11.P10spals = 19.20;
fsamples.Ling2006a_11.P10muon = 0.136;
fsamples.Ling2006a_11.P26spals = fsamples.Ling2006a_11.P10spals*6.75;
fsamples.Ling2006a_11.P26total = fsamples.Ling2006a_11.P26spals*1.0204;
fsamples.Ling2006a_11.P26muon = fsamples.Ling2006a_11.P26total*0.02;


%Nesj2007_9
fsamples.Nesj2007_9.lat = 69.218000000000004;
fsamples.Nesj2007_9.lon = 15.909000000000001;
fsamples.Nesj2007_9.N10 = 210000;
fsamples.Nesj2007_9.N10unc = 14000;
fsamples.Nesj2007_9.nnc = 1;
fsamples.Nesj2007_9.P10total = 4.528;
fsamples.Nesj2007_9.P10spals = 4.45;
fsamples.Nesj2007_9.P10muon = 0.078;
fsamples.Nesj2007_9.P26spals = fsamples.Nesj2007_9.P10spals*6.75;
fsamples.Nesj2007_9.P26total = fsamples.Nesj2007_9.P26spals*1.0204;
fsamples.Nesj2007_9.P26muon = fsamples.Nesj2007_9.P26total*0.02;

%Hall2019_12
fsamples.Hall2019_12.lat = 60.379240000000003;
fsamples.Hall2019_12.lon = 18.235260000000000;
fsamples.Hall2019_12.N10 = 315194;
fsamples.Hall2019_12.N26 = 2086003;
fsamples.Hall2019_12.N10unc = 8751;
fsamples.Hall2019_12.N26unc = 78972;
fsamples.Hall2019_12.nnc = 2;
fsamples.Hall2019_12.P10total = 4.186;
fsamples.Hall2019_12.P10spals = 4.11;
fsamples.Hall2019_12.P10muon = 0.076;
fsamples.Hall2019_12.P26tot = 28.447;
fsamples.Hall2019_12.P26spals = 27.75;
fsamples.Hall2019_12.P26muon = 0.697;


%Stroe2002a_1
fsamples.Stroe2002a_1.lat = 67.677199999999999;
fsamples.Stroe2002a_1.lon = 23.072500000000002;
fsamples.Stroe2002a_1.N10 = 470742;
fsamples.Stroe2002a_1.N26 = 2965550;
fsamples.Stroe2002a_1.N10unc = 24710;
fsamples.Stroe2002a_1.N26unc = 175728;
fsamples.Stroe2002a_1.nnc = 2;
fsamples.Stroe2002a_1.P10total = 6.067;
fsamples.Stroe2002a_1.P10spals = 5.98;
fsamples.Stroe2002a_1.P10muon = 0.087;
fsamples.Stroe2002a_1.P26tot = 41.180;
fsamples.Stroe2002a_1.P26spals = 40.38;
fsamples.Stroe2002a_1.P26muon = 0.8;

%Darm2008_1
fsamples.Darm2008_1.lat = 68.122249999999994;
fsamples.Darm2008_1.lon = 27.370120000000000;
fsamples.Darm2008_1.N10 = 389474;
fsamples.Darm2008_1.N26 = 2094000;
fsamples.Darm2008_1.N10unc = 17544;
fsamples.Darm2008_1.N26unc = 134000;
fsamples.Darm2008_1.nnc = 2;
fsamples.Darm2008_1.P10total = 7.001;
fsamples.Darm2008_1.P10spals = 6.91;
fsamples.Darm2008_1.P10muon = 0.091;
fsamples.Darm2008_1.P26tot = 47.486;
fsamples.Darm2008_1.P26spals = 46.64;
fsamples.Darm2008_1.P26muon = 0.846;


%constrain Tdgla from Stroeven and Hughes

%from stroeven (2016) isochrones
fsamples.Sven2015_12.Tdgla = 20.0e-3;
fsamples.Ling2006a_11.Tdgla = 10.25e-3;
fsamples.Nesj2007_9.Tdgla = 18.0e-3;
fsamples.Hall2019_12.Tdgla = 10.75e-3;
fsamples.Stroe2002a_1.Tdgla = 10.1e-3;
fsamples.Darm2008_1.Tdgla = 10.5e-3;



for i = 1:length(fsamples.IDs)
    fsamples.(fsamples.IDs{i}).d18Ot = parameterMeans.(fsamples.IDs{i})(1);
    %Tdgla are set manually above.
    fsamples.(fsamples.IDs{i}).Z1 = parameterMeans.(fsamples.IDs{i})(3);
    fsamples.(fsamples.IDs{i}).dT2 = parameterMeans.(fsamples.IDs{i})(4);
    fsamples.(fsamples.IDs{i}).dZ2 = parameterMeans.(fsamples.IDs{i})(5);
    fsamples.(fsamples.IDs{i}).dT3 = parameterMeans.(fsamples.IDs{i})(6);
    fsamples.(fsamples.IDs{i}).dZ3 = parameterMeans.(fsamples.IDs{i})(7);
    fsamples.(fsamples.IDs{i}).E0 = parameterMeans.(fsamples.IDs{i})(8);
end


save fsamples.mat fsamples