clc;clear;close all;
%% Parse data
c = 299792458e-03; %(km/s)
sig_t = 100e-09;
addpath('./MatFiles');
addpath('./Particle Filter Testing Files')
load AIM_Range.mat
Truth = readmatrix('AIM_True_Position_1day.txt');
% sets up times in Year,Month,day,hour,minute,sec
AIM.Times = cell2mat(finalCell(:,1));
% sets up ranges from Sat to Boulder, Wyoming border, Kremmling, Pueblo
AIM.Ranges = cell2mat(finalCell(:,2));

los = 1500;
% Find the indices were all sats are flying by limiting range to 2000 km
Range_Mask1 = AIM.Ranges(:,1) <= los;
Range_Mask2 = AIM.Ranges(:,2) <= los;
Range_Mask3 = AIM.Ranges(:,3) <= los;
Range_Mask4 = AIM.Ranges(:,4) <= los;

Flyby_Mask = Range_Mask1 & Range_Mask2 & Range_Mask3 & Range_Mask4;
ind_Mask = find(Flyby_Mask == 1);

% add the fly by indices to the mask.
AIM.ind = ind_Mask;
% Loads ECI frame sensor locations
Rs = [lla2ecef([40.00888, -105.24774,  1612]); % Boulder
      lla2ecef([40.935583, -105.380917,1868]); % Wymoing  Border
	  lla2ecef([39.95121, -106.34978,  2331]); % Kremmmling
      lla2ecef([38.24047, -104.57511,  2000])];% Pueblo
% Rs = [lla2ecef([40.00888, -105.24774,  0]); % Boulder
%       lla2ecef([40.935583, -105.380917,0]); % Wymoing  Border
% 	  lla2ecef([39.95121, -106.34978,  2331]); % Kremmmling
%       lla2ecef([38.24047, -104.57511,  0])];% Pueblo 
% Find ECI Locatoins of the sensors in ECI 
for j =1:length(Rs)
    for i = 1:length(ind_Mask)
        % Creates a Direction Cosine Matrix of eci2ecef
        dcm = dcmeci2ecef('IAU-2000/2006',AIM.Times(ind_Mask(i),:));
        R_ECI(i,:) = dcm\Rs(j,:)';
    end
    Sensors(j).ECEF = Rs(j,:)'*10^(-3);
    % Saves to struct and converts into km
    Sensors(j).ECI = R_ECI'*10^(-3);
end
% Add in the ecef to Sensor

L = length(ind_Mask);
% Create TDoA Data
t1 = AIM.Ranges(ind_Mask,1)/c + sig_t*randn(L,1);
t2 = AIM.Ranges(ind_Mask,2)/c + sig_t*randn(L,1);
t3 = AIM.Ranges(ind_Mask,3)/c + sig_t*randn(L,1);
t4 = AIM.Ranges(ind_Mask,4)/c + sig_t*randn(L,1);

TDoA.SN1 = [t1 - t1, t2 - t1, t3 - t1, t4 - t1];
TDoA.SN2 = [t1 - t2, t2 - t2, t3 - t2, t4 - t2];
TDoA.SN3 = [t1 - t3, t2 - t3, t3 - t3, t4 - t3];
TDoA.SN4 = [t1 - t4, t2 - t4, t3 - t4, t4 - t4];


AIMTrue.R = Truth(:,5:7);
AIMTrue.V = Truth(:,9:11);

% later I need to add in the ECEF coordinates of the A

%% Save to .mat files
save('AIMRange.mat','AIM')
save('AIM_Sensors.mat', 'Sensors')
save('AIMTDoA.mat','TDoA')
save('AIMTrue.mat','AIMTrue')