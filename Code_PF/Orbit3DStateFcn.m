function particles = Orbit3DStateFcn(particles,Nratio) 

addpath('./MatFiles');
[numberOfStates, numberOfParticles] = size(particles);
c = 299792458 * 10^(-3); %(km/s)
R_E = 6378; %[km]


dt = 5; % [s] Sample time
for kk=1:numberOfParticles
    particles(:,kk) = RK4(0,particles(:,kk),dt,@OrbitalODE3D);
    % calculate the distance to each particle
end
% Add Gaussian noise with variance 0.1 on each state variable
%was R= 0.75 V=0.15
processNoiseR = 0.75 * eye(3);
processNoiseV = 0.15 * eye(3);
if Nratio == 1
%     rng(42)
    particles(1:3,:) = particles(1:3,:) + processNoiseR * randn( 3, numberOfParticles);
    particles(4:6,:) = particles(4:6,:) + processNoiseV * randn( 3, numberOfParticles);
else
%     rng(42)
    processNoiseR = Nratio * processNoiseR;
    processNoiseV = Nratio * processNoiseV;
    particles(1:3,:) = particles(1:3,:) + processNoiseR * randn( 3, numberOfParticles);
    particles(4:6,:) = particles(4:6,:) + processNoiseV * randn( 3, numberOfParticles);
end
end
function sols = OrbitalODE3D(t,X)
    mu = 398600;
    R=[X(1) ; X(2); X(3)];
    V=[X(4) ; X(5); X(6)];
    % Accelertion based on 
    A=(-mu/(norm(R))^3) * R;
    
    sols=[V ; A ];

end

function[out] = RK4(t,X,h,ODE)
%Simple Runge-Kutta Method integrator
    k1 = ODE(t,X);
    k2 = ODE(t+h/2,X+h*k1/2);
    k3 = ODE(t+h/2,X+h*k2/2);
    k4 = ODE(t+h,X+h*k3);
    out = X+h/6*(k1+2*k2+2*k3+k4);

end