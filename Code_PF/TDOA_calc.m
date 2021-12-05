function[XYZ]=TDOA_calc(P,c,err_r,tdoa)
% This is based on TDoA file from MATLAB file exchange

method=1;


trials=100;
M=length(P);
in_est_error=0;
% Induce Noise into Sensor Position SEED For repeatablility
% rng(40)


og = find(tdoa==0);
% makes an initial guess for this case 575 km above LASP 
guess = P(1:3,1)+ P(1:3,1)/norm(P(1:3,1))*575;
for k=1:trials

% %finding TOAs 

% Add in a Perturbation of initial guess
Perr = P+err_r;
% tdoa = tdoa +err_t;



    
%%% Taylor Series Expansion Solution
p_T_0 = guess + in_est_error*randn(size(guess));    %initial estimate with some error (penalty term)
d = c*tdoa';
f = zeros(M,1);
del_f = zeros(M,2);

if method==1
    


for ii=1:M
   % 
   f(ii)       = norm(p_T_0-Perr(:,ii)) - norm(p_T_0-Perr(:,og)); 
   del_f(ii,1) = (p_T_0(1)-Perr(1,ii))*norm(p_T_0-Perr(:,ii))^-1 - (p_T_0(1)-Perr(1,og))*norm(p_T_0-Perr(:,og))^-1;
   del_f(ii,2) = (p_T_0(2)-Perr(2,ii))*norm(p_T_0-Perr(:,ii))^-1 - (p_T_0(2)-Perr(2,og))*norm(p_T_0-Perr(:,og))^-1;
   del_f(ii,3) = (p_T_0(3)-Perr(3,ii))*norm(p_T_0-Perr(:,ii))^-1 - (p_T_0(3)-Perr(3,og))*norm(p_T_0-Perr(:,og))^-1;    

end
% Solve in the direction of the gradient
   x_nonlin = pinv(del_f)*(d-f)+p_T_0;
   % set the next guess equal to the direction.
   guess    = x_nonlin;
   X        = x_nonlin;
else
    % Construct an A and b to solve Ax = b
    A= -2 * [P(:,2)'- P(:,1)', d(2);
             P(:,3)'- P(:,1)', d(3); 
             P(:,4)'- P(:,1)', d(4)];
    b=[ d(2)^2-norm(P(:,2))^2+norm(P(:,1))^2;
        d(3)^2-norm(P(:,3))^2+norm(P(:,1))^2;
        d(4)^2-norm(P(:,4))^2+norm(P(:,1))^2];
% Calculate position based on pseudoinverse
    X=pinv(A)*b;
    
end
end



XYZ=X;
end
