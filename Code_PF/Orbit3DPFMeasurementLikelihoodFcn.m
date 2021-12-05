function likelihood = Orbit3DPFMeasurementLikelihoodFcn(particles,measurement,filename)

addpath('./MatFiles');
load(filename,'Sensors')

c = 299792458 * 10^(-3); %(km/s)
R_E = 6378; %[km]
sig_t = 100 *10^(-9) ; % s

% pass in a indice for the R_eci sensor locations
ind = measurement(5);

xTrue = measurement(6:end);
% Round to 100ns resolution to be more representative
measurement = round(measurement(2:4),7);

% Validate the sensor measurement
numberOfMeasurements = 3; % Expected number of measurements
validateattributes(measurement, {'double'}, {'vector', 'numel', numberOfMeasurements}, ...
    'vdpExamplePFMeasurementLikelihoodFcn', 'measurement');
% predict the measurement
predictedMeasurement = particles(1:3,:);

% calculate the distance to each particle
Dist2S1 = vecnorm( predictedMeasurement - Sensors(1).ECI(:,ind),2,1);
Dist2S2 = vecnorm( predictedMeasurement - Sensors(2).ECI(:,ind),2,1);
Dist2S3 = vecnorm( predictedMeasurement - Sensors(3).ECI(:,ind),2,1);
Dist2S4 = vecnorm( predictedMeasurement - Sensors(4).ECI(:,ind),2,1);

% calculate the time of arrivial
t1 = round(Dist2S1 / c,7);
t2 = round(Dist2S2 / c,7);
t3 = round(Dist2S3 / c,7);
t4 = round(Dist2S4 / c,7);

t1 = Dist2S1 / c;
t2 = Dist2S2 / c;
t3 = Dist2S3 / c;
t4 = Dist2S4 / c;
% Round to 100ns resolution to be more representative
t = [t1; t2; t3; t4]; 

% calculate the T ime delay of arrival for each particle
tdoa1 =  t - t(1,:);
tdoa2 = [ t1 - t2; t2 - t2; t3 - t2; t4 - t2];
tdoa3 = [ t1 - t3; t2 - t3; t3 - t3; t4 - t3];
tdoa4 = [ t1 - t4; t2 - t4; t3 - t4; t4 - t4];

predictedMeasurement = tdoa1(2:end,:);
sigma = sig_t^2 * eye(numberOfMeasurements); % variance
numParticles = size(particles,2);
likelihood = zeros(numParticles,1);

% Coefficient multiplied by the exponenet
C = det(2*pi*sigma) ^ (-0.5);
for kk=1:numParticles
    if norm(particles(1:3,kk)) < R_E 
        likelihood(kk)=0;
    else
        v = -(predictedMeasurement(:,kk)-measurement);
        tmp = sigma\v;
        InTheExponent(kk) = -0.5 * ( v' * tmp); 
        likelihood(kk) =  C *exp(InTheExponent(kk));
    end
end
%weights = likelihood/(sum(likelihood));
% Neff = 1/sum(weights.^2);
% Nratio = Neff / numParticles;
%Plot_PF_weight(particles(1:3,:),xTrue',weights)
end