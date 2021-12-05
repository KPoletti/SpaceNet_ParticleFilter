clc;clear;close all;
%% Parse data
c = 299792458e-03; %(km/s)
sig_t = 100e-09;
addpath('./MatFiles');
addpath('./Particle Filter Testing Files')
sat = 'AIM';
% create  True data and range data filenames then load in
Mat1 = strcat(sat,'_Range.mat');
Txt1 = strcat(sat,'_True_Position_1day.txt');
load(Mat1, 'finalCell')
Truth = readmatrix(Txt1);


% sets up times in Year,Month,day,hour,minute,sec
AIM.Times = cell2mat(finalCell(:,1));
% sets up ranges from Sat to Boulder, Wyoming border, Kremmling, Pueblo
AIM.Ranges = cell2mat(finalCell(:,2));

% find truth data
AIMTrue.R = Truth(:,5:7);
AIMTrue.V = Truth(:,9:11);

% Loads ECI frame sensor locations
Rs = [lla2ecef([40.00888, -105.24774,  1612]); % Boulder
      lla2ecef([40.935583, -105.380917,1868]); % Wymoing  Border
	  lla2ecef([39.95121, -106.34978,  2331]); % Kremmmling
      lla2ecef([38.24047, -104.57511,  2000])];% Pueblo

% save ECEF positions
Sensors(1).ECEF = Rs(1,:)'*10^(-3);
Sensors(2).ECEF = Rs(2,:)'*10^(-3);
Sensors(3).ECEF = Rs(3,:)'*10^(-3);
Sensors(4).ECEF = Rs(4,:)'*10^(-3);
longs = [-105.24774, -105.380917, -106.34978,-104.57511];
lats = [40.00888,40.935583,39.95121,38.24047];

% Find azimuth elvation and range
aer1 = eci2aer(AIMTrue.R*10^3,AIM.Times,[40.00888, -105.24774,  1612].*ones(length(AIMTrue.R),1),'IAU-2000/2006');
aer2 = eci2aer(AIMTrue.R*10^3,AIM.Times,[40.935583, -105.380917,1868].*ones(length(AIMTrue.R),1),'IAU-2000/2006');
aer3 = eci2aer(AIMTrue.R*10^3,AIM.Times,[39.95121, -106.34978,  2331].*ones(length(AIMTrue.R),1),'IAU-2000/2006');
aer4 = eci2aer(AIMTrue.R*10^3,AIM.Times,[38.24047, -104.57511,  2000].*ones(length(AIMTrue.R),1),'IAU-2000/2006');
%%
los = 15;
% Find the indices were all sats are flying by limiting range to 2000 km
Range_Mask1 = aer1(:,2) >= los;
Range_Mask2 = aer2(:,2) >= los;
Range_Mask3 = aer3(:,2) >= los;
Range_Mask4 = aer4(:,2) >= los;

Flyby_Mask = Range_Mask1 & Range_Mask2 & Range_Mask3 & Range_Mask4;
ind_Mask = find(Flyby_Mask == 1);

% add the fly by indices to the mask.
AIM.ind = ind_Mask;
L = length(ind_Mask);
% Find ECI Locatoins of the sensors in ECI 
for j =1:length(Rs)
    for i = 1:length(ind_Mask)
        % Creates a Direction Cosine Matrix of eci2ecef
        dcm = dcmeci2ecef('IAU-2000/2006',AIM.Times(ind_Mask(i),:));
        AIMTrue.RECI(i,:) = dcm*AIMTrue.R(ind_Mask(i),:)';
        DCM(i).ECI2ECEF = dcm;
        R_ECI(i,:) = dcm\Rs(j,:)';
    end
    % Saves to struct and converts into km
    Sensors(j).ECI = R_ECI'*10^(-3);
end
%%

% Create TDoA Data
t1 = AIM.Ranges(ind_Mask,1)/c + sig_t*randn(L,1);
t2 = AIM.Ranges(ind_Mask,2)/c + sig_t*randn(L,1);
t3 = AIM.Ranges(ind_Mask,3)/c + sig_t*randn(L,1);
t4 = AIM.Ranges(ind_Mask,4)/c + sig_t*randn(L,1);

TDoA.SN1 = [t1 - t1, t2 - t1, t3 - t1, t4 - t1];
TDoA.SN2 = [t1 - t2, t2 - t2, t3 - t2, t4 - t2];
TDoA.SN3 = [t1 - t3, t2 - t3, t3 - t3, t4 - t3];
TDoA.SN4 = [t1 - t4, t2 - t4, t3 - t4, t4 - t4];



% later I need to add in the ECEF coordinates of the A

%% Save to .mat files
filename1 = strcat(sat,'Range.mat');
filename2 = strcat(sat,'_Sensors.mat');
filename3 = strcat(sat,'TDoA.mat');
filename4 = strcat(sat,'True.mat');
save(['./MatFiles\' filename1],'AIM')
save(['./MatFiles\' filename2], 'Sensors')
save(['./MatFiles\' filename3],'TDoA')
save(['./MatFiles\' filename4],'AIMTrue')
save(['./MatFiles\' 'ECI2ECEF.mat'],'DCM')

function [rhoNew] = topHorz(time, long,lat, rho)

% converts refernce frames from local to Topological Horizon.
year = time(1);
m = time(2);
d = time(3);
h = time(4);
min = time(5);
s = time(6);


J0 = 367*year-fix(7*(year+fix((m+9)/12))/4)+fix(275*m/9)+d+1721013.5;
T0=(J0-2451545)/36525;
Th0=100.4606184+36000.77004*T0+.000387933*T0^2-2.583*10^(-8)*T0^3;
% OMdot=0.985647*pi/180/(24*3600);
% period =1.62 *3600;
% a=(sqrt(mu)/(2*pi)*period)^(2/3);
UT = h + min/60 + s/3600;

Th0 = Th0-fix(Th0/360)*360;

thetaG = Th0 + UT/24 *360.98564724;

theta = thetaG + long;
theta = theta - fix(theta/360)*360;
st = sind(theta);
ct = cosd(theta);
sp = sind(lat);
cp = sind(lat);
dcm = [ -st,     ct,     0;
       -sp*ct, -sp*st, cp;
       cp*ct,  cp*st, sp ];
rhoNew = dcm*rho;
       

end