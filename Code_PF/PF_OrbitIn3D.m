%%
clc;clear;close all;
%% 3D Orbit Particle filter.
%{
Will input Raw Time delay of arrival data, ie the time between the signal
arriving. With the data, it will calculate an initial position with the
TDoA algorithm and an initial velocity with Lambert's method.

Then the PF will be initialized and plot everything
%}
% 1/19/2021
%By Keith Poletti
%% Plot options
DoPlot = 1;
plotlims = 1;
% record mp4 of the particles
DoFilm = 0;
%% Pick a flyby and the satellite name
ind1 = 3;
sat = 'AIM';
%% Load in data
addpath('./MatFiles');


% AIM Data
filename1 = strcat(sat,'Range.mat');
filename2 = strcat(sat,'_Sensors.mat');
filename3 = strcat(sat,'TDoA.mat');
filename4 = strcat(sat,'True.mat');
load(filename1)
load(filename2)
load(filename3)
load(filename4)
load('ECI2ECEF.mat')
%% Constants
R_E = 6378; %[km]
R_sat = R_E + 575; %[km]
mu = 398600; %[km^3/s^2]
c = 299792458 * 10^(-3); %(km/s)
sig_t = 100 *10^(-9); % s
sig_r = 1e-03; %km
%Line of sight max distance for a unit to see the satellite
% Indices of flybys
Ind_TDoA = AIM.ind;
% find the indice where each fly by starts
Starts = [1;find(diff(Ind_TDoA)>1)+1];
%% Initial Conditions

% P must be in ECEF for good convergrance
P = [Sensors(1).ECEF, Sensors(2).ECEF,...
     Sensors(3).ECEF, Sensors(4).ECEF];

% pick an arbitrary flyby to analyze
ind1 = Starts(ind1);
ind2 = ind1+7;

% calculate ECEF position vectors for two points in flyby
XYZ1 = TDOA_calc(P,c,sig_r,TDoA.SN1(ind1,:));
XYZ2 = TDOA_calc(P,c,sig_r,TDoA.SN1(ind2,:));
% load TimesOfFlyby.mat

Times = AIM.Times;
dt = abs(Times(2,end)-Times(1,end));
% found locations convert from ECEF to ECI 
dcm = dcmeci2ecef('IAU-2000/2006',Times(Ind_TDoA(ind1),:));

test = dcm*AIMTrue.R(Ind_TDoA(ind1),:)';
tdoa = vecnorm((P-test)/c,2,1);
tdoa = tdoa - tdoa(1);

XYZ1 = dcm\XYZ1;
az = atan2(XYZ1(2),XYZ1(1))*180/pi + 90;

dcm = dcmeci2ecef('IAU-2000/2006',Times(Ind_TDoA(ind2),:));
XYZ2 = dcm\XYZ2;

% Lamberts does a good job it will probably need to have the data spread
% over a minute for best accuracy with 0.5km/s
[V1,V2] = lambert(XYZ1,XYZ2,(ind2-ind1)*dt,'retro');
V1 = real(V1);
V2 = real(V2);


%% PF Time
% Create initial PF state first position from TDoA and Velocity from
% Lambert's
Initial=[XYZ1; V1];

% start up PF
pf = particleFilter(@Orbit3DStateFcn,@Orbit3DPFMeasurementLikelihoodFcn);

% initialize PF with 10,000 particles, 10 km cov and 1.5km/s covariance
% change options on the particle filter.
initialize(pf, 10000,Initial,diag([10^2*ones(3,1); 1.5^2*ones(3,1)]));
pf.ResamplingMethod = 'systematic';
pf.ResamplingPolicy.MinEffectiveParticleRatio = 0.050;
% pf.ResamplingMethod = 'residual';
% pf.ResamplingMethod = 'stratified';

%%
% define the true values 
xTrueMeas = [AIMTrue.R(Ind_TDoA,:),AIMTrue.V(Ind_TDoA,:)];
%Allocate

xCorrectedPF = zeros(size(xTrueMeas));
COV=zeros(size(xTrueMeas));

% Sensor unit 1 or LASP TDoA measurements
yMeas = TDoA.SN1';

%Make a sphere for Earth
[N,E,W]=sphere;
if ind1 == Starts(end)
    endOfFB = length(Ind_TDoA);
else 
    endOfFB = Starts(Starts>ind1)-1; %end of flyby
    endOfFB = endOfFB(1);
end
set(0,'defaultfigurecolor',[1 1 1])
if DoPlot ==1
    if DoFilm ==1
        vidfile = VideoWriter('Orbital3DparticlesGoZoom.mp4','MPEG-4');
        vidfile.FrameRate=5;
        open(vidfile)
        vidPart = VideoWriter('ParticlesOrigin.mp4','MPEG-4');
        vidPart.FrameRate=5;
        open(vidPart)
    end
    % Plot a lot more stuff 
    neon_bl = [31, 81, 255]/255;
    % subplot(2,1,1)
    figure(1)
    k = ind1;
    p = DCM(k).ECI2ECEF*pf.Particles(1:3,:);
    state = getStateEstimate(pf);
    state_ECEF = DCM(k).ECI2ECEF*state(1:3);
    plot3(state_ECEF(1),state_ECEF(2),state_ECEF(3),'ob','linewidth',2)
    hold on;
    plot3(AIMTrue.RECI(k,1),AIMTrue.RECI(k,2), AIMTrue.RECI(k,3),'*m','linewidth',2)
    %plot orbit
    plot3(AIMTrue.RECI(ind1:endOfFB,1),AIMTrue.RECI(ind1:endOfFB,2),AIMTrue.RECI(ind1:endOfFB,3),'m','linewidth',2)
    plot3(p(1),p(2),p(3),'.c','MarkerIndices',floor(linspace(1,length(pf.Particles),10)),'linewidth',0.5)

    %plot Measurement location
    plot3(Sensors(1).ECEF(1),Sensors(1).ECEF(2),Sensors(1).ECEF(3),'*g')
    plot3(Sensors(2).ECEF(1),Sensors(2).ECEF(2),Sensors(2).ECEF(3),'*g')
    plot3(Sensors(3).ECEF(1),Sensors(3).ECEF(2),Sensors(3).ECEF(3),'*g')
    plot3(Sensors(4).ECEF(1),Sensors(4).ECEF(2),Sensors(4).ECEF(3),'*g')
    earth_sphere(50,'km')
    xlabel('Position [x]')
    ylabel('Position [y]')
    zlabel('Position [z]')
    legend('Particles','PF Estimate','STK Satellite Position','STK Orbit','Sensors','location','northeast')
    titlestr=strcat('Iteration=',num2str(0));
    title('Orbital View of Particle Filter')
    view(az,15)
%             xlim([0 5100])
%         ylim([-6050 -1700])
%         zlim([2000 5500])
%     xlim([0 4000])
%     ylim([-6000 -2500])
%     zlim([2000 5500])
    xlim([-3250 1100])
    ylim([-6000 -2500])
    zlim([2000 5300])
    drawnow
    hold off;
    if DoFilm ==1
        im=getframe(gcf);
        writeVideo(vidfile, im)
    end
    % subplot(2,1,2)
    figure(2);
    plot3(pf.Particles(1,:)-xTrueMeas(ind1,1),pf.Particles(2,:)-xTrueMeas(ind1,2),pf.Particles(3,:)-xTrueMeas(ind1,3),'.c','linewidth',0.5)
    hold on;
    plot3(pf.State(1)-xTrueMeas(ind1,1),pf.State(2)-xTrueMeas(ind1,2),pf.State(3)-xTrueMeas(ind1,3),'*b','linewidth',3)

    plot3(0,0,0,'*m','linewidth',3);
    grid on
    xlabel('X Position [km]')
    ylabel('Y Position [km]')
    zlabel('Z Position [km]')
    xlim([-50 50])
    ylim([-50 50])
    zlim([-50 50])
    title('Zoomed in with STK Satellite Position at the Origin')
    legend('Particles','PF Estimate','STK Satellite Position','location','northeast')
    % xlim([-2 2])
    % ylim([-2 2])

    drawnow
    hold off;
    if DoFilm ==1
        im1=getframe(gcf);
        writeVideo(vidPart, im1)
    end
end
% Allocate some more
xCorrectedPF(ind1,:) = getStateEstimate(pf);
XYZ = zeros(length(ind1:endOfFB),3);
TDoAres = zeros(size(XYZ));
timeVector = ind1:endOfFB;
TDoAPF =zeros(size(yMeas));
TDoATr =zeros(size(yMeas));
for k=ind1:endOfFB
    % Use measurement y[k] to correct the particles for time k

    xCorrectedPF(k,:) = correct(pf,[yMeas(:,k);k;xTrueMeas(k,1:3)'],filename2); % Filter updates and stores Particles[k|k], Weights[k|k]

    Neff = 1/sum(pf.Weights.^2);
    Nratio(k-ind1+1,1) = Neff / pf.NumParticles;
    
    % The result is x[k|k]: Estimate of states at time k, utilizing
    % measurements up to time k. This estimate is the mean of all particles
    % because StateEstimationMethod was 'mean'.
    %
    % Now, predict particles at next time step. These are utilized in the
    % next correct command
    %
    % Find the position if we only did TDoA this will  be in ECEF.
    XYZ(k -ind1+1,:) =  TDOA_calc(P,c,sig_r,yMeas(:,k)');
    TDoAres(k -ind1+1,:) = AIMTrue.RECI(k,1:3) -  XYZ(k -ind1+1,:);
    %Find  Covariance
    [state,tmp] = getStateEstimate(pf);
    COV(k+1,:)=diag(tmp);
    tmp = diag(tmp);
    if DoPlot ==1
    figure(1);
%     subplot(1,2,1)
    state_ECEF = DCM(k).ECI2ECEF*state(1:3);

    plot3(state_ECEF(1),state_ECEF(2),state_ECEF(3),'ob','linewidth',2)
    hold on;
    %plot Measurement
    p = DCM(k).ECI2ECEF*pf.Particles(1:3,:);
    px = p(1) + 3*sqrt(tmp(1))*randn(size(pf.Particles(1,:)));
    py = p(2) + 3*sqrt(tmp(2))*randn(size(pf.Particles(1,:)));
    pz = p(3) + 3*sqrt(tmp(3))*randn(size(pf.Particles(1,:)));
    %plot k true state
    

    plot3(AIMTrue.RECI(k,1),AIMTrue.RECI(k,2), AIMTrue.RECI(k,3),'*m','linewidth',2)
    %plot orbit
    plot3(AIMTrue.RECI(ind1:endOfFB,1),AIMTrue.RECI(ind1:endOfFB,2),AIMTrue.RECI(ind1:endOfFB,3),'m','linewidth',2)
    plot3(px,py,pz,'.c','linewidth',2)

    plot3(Sensors(1).ECEF(1),Sensors(1).ECEF(2),Sensors(1).ECEF(3),'*g')
    plot3(Sensors(2).ECEF(1),Sensors(2).ECEF(2),Sensors(2).ECEF(3),'*g')
    plot3(Sensors(3).ECEF(1),Sensors(3).ECEF(2),Sensors(3).ECEF(3),'*g')
    plot3(Sensors(4).ECEF(1),Sensors(4).ECEF(2),Sensors(4).ECEF(3),'*g')
    earth_sphere(50,'km')
    legend('PF Estimate','STK Satellite Position','STK Orbit','Particles','Sensors','location','northeast')

    view(az,15)
    xlabel('Position X [km]')
    ylabel('Position Y [km]')
    zlabel('Position Z [km]')
    titlestr=strcat('Iteration=',num2str(k));
    title('Orbital View of Particle Filter')
    if plotlims==1
%         xlim([0 5100])
%         ylim([-6050 -1700])
%         zlim([2000 5500])
%     xlim([0 4000])
%     ylim([-6000 -2500])
%     zlim([2000 5500])
    xlim([-3250 1100])
    ylim([-6000 -2500])
    zlim([2000 5300])
    end
    if DoFilm == 1
        im=getframe(gcf);
        writeVideo(vidfile, im);
    end
    hold off
    figure(2)
%     subplot(1,2,2)
    plot3(pf.Particles(1,:)-xTrueMeas(k,1),pf.Particles(2,:)-xTrueMeas(k,2),pf.Particles(3,:)-xTrueMeas(k,3),'.c','Markerindices',floor(linspace(1,length(pf.Particles),500)),'linewidth',0.5)
    hold on;
    grid on
    plot3(pf.State(1)-xTrueMeas(k,1),pf.State(2)-xTrueMeas(k,2),pf.State(3)-xTrueMeas(k,3),'*b','linewidth',3)
    plot3(0,0,0,'*m','linewidth',3);
    if plotlims==1
        xlim([-50 50])
        ylim([-50 50])
        zlim([-50 50])
    end
    legend('Particles','PF Estimate','STK Satellite Position','location','northeast')
    xlabel('X Position Error [km]')
    ylabel('Y Position Error [km]')
    zlabel('Z Position Error [km]')
    title('Zoomed in View with STK Satellite Position at the Origin')
    % xlim([-2 2])
    % ylim([-2 2])

    drawnow
    hold off;
    if DoFilm == 1
        im1=getframe(gcf);
        writeVideo(vidPart, im1);
    end
    end
%   Update and predict
    predict(pf, Nratio(k-ind1+1,1)); % Filter updates and stores Particles[k+1|k]
end
if DoFilm == 1
close(vidfile)
close(vidPart)
end
%% Plot the data
res = Plot_PF(timeVector,xTrueMeas,xCorrectedPF,k,COV,DoFilm,dt,TDoAres);
%%
Rres = res(:,1:3);
Rres_mag = vecnorm(Rres,2,2);
[Rres_min, Rres_ind] = min(Rres_mag);
Vres = res(:,4:6);
Vres_mag = vecnorm(Vres,2,2);
[Vres_min, Vres_ind] = min(Vres_mag);
[res_min, res_ind] = min(vecnorm([Rres_mag,Vres_mag],2,2));
fprintf('Residual in model: \n')
fprintf('Minimum Position residual: %0.3f [km] at %d \n',Rres_min, Rres_ind) 
fprintf('Minimum Velocity residual: %0.3f [km/s] at %d \n',Vres_min, Vres_ind) 
fprintf('Minimum Overall residual: R = %0.3f km V = %0.3f km/s  at %d \n \n',Rres_mag(res_ind),Vres_mag(res_ind), res_ind) 

fprintf('Residual in model: \n')
fprintf('Mean Position residual: %0.3f [km] \n',sqrt(mean((vecnorm(Rres,1,2)).^2))) 
fprintf('Mean Velocity residual: %0.3f [km/s] \n',sqrt(mean((vecnorm(Vres,1,2)).^2,'all'))) 
fprintf('Minimum Overall residual: R = %0.3f km V = %0.3f km/s  at %d \n \n',Rres_mag(res_ind),Vres_mag(res_ind), res_ind) 
indCOV = ind1+1;
Rres = sqrt(COV(indCOV:endOfFB,1:3));
Rres_mag = vecnorm(Rres,2,2);
[Rres_min, Rres_ind] = min(Rres_mag);
Vres = sqrt(COV(indCOV:endOfFB,4:6));
Vres_mag = vecnorm(Vres,2,2);
[Vres_min, Vres_ind] = min(Vres_mag);
[res_min, res_ind] = min(vecnorm([Rres_mag,Vres_mag],2,2));
fprintf('Error in model: \n')
fprintf('Minimum Position Uncertainty: %0.3f [km] at %d \n',Rres_min, Rres_ind+1) 
fprintf('Minimum Velocity Uncertainty: %0.3f [km/s] at %d \n',Vres_min, Vres_ind+1) 
fprintf('Minimum Overall Uncertainty: R = %0.3f km V = %0.3f km/s  at %d \n',Rres_mag(res_ind),Vres_mag(res_ind), res_ind+1) 
function [TDoAPF] = tdoa_times(predictedMeasurement, Sensors, ind)
    c = 299792458 * 10^(-3); %(km/s)  
    Dist2S1 = vecnorm( predictedMeasurement - Sensors(1).ECI(:,ind),2,1);
    Dist2S2 = vecnorm( predictedMeasurement - Sensors(2).ECI(:,ind),2,1);
    Dist2S3 = vecnorm( predictedMeasurement - Sensors(3).ECI(:,ind),2,1);
    Dist2S4 = vecnorm( predictedMeasurement - Sensors(4).ECI(:,ind),2,1);

    % calculate the time of arrivial
    t1 = round(Dist2S1 / c,7);
    t2 = round(Dist2S2 / c,7);
    t3 = round(Dist2S3 / c,7);
    t4 = round(Dist2S4 / c,7);
    % Round to 100ns resolution to be more representative
    t = [t1; t2; t3; t4]; 

    % calculate the T ime delay of arrival for each particle
    TDoAPF =  t - t(1,:);
end
