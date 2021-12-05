function [res ] = Plot_PF(timeVector,xTrueMeas,xCorrectedPF,k,COV,DoFilm,dt,TDoAres)
figure(2);
sgtitle('Position Measurement and predicitions')
subplot(3,1,1);
ind = timeVector;
timeVector = (timeVector - timeVector(1))* dt;
plot(timeVector,xTrueMeas(ind,1),'r',timeVector,xCorrectedPF(ind,1),'b');
hold on
grid on
plot(timeVector,xCorrectedPF(ind,1)-3*sqrt(COV(ind,1)),'-m',timeVector,xCorrectedPF(ind,1)+3*sqrt(COV(ind,1)),'-m')
legend('True','Particlte filter estimate','Measured','3 \sigma')
% ylim([-1 1]);
xlim([min(timeVector) max(timeVector)])
ylabel('X [km]');


subplot(3,1,2);
plot(timeVector,xTrueMeas(ind,2),'r',timeVector,xCorrectedPF(ind,2),'b');
hold on
grid on
plot(timeVector,xCorrectedPF(ind,2)-3*sqrt(COV(ind,2)),'-m',timeVector,xCorrectedPF(ind,2)+3*sqrt(COV(ind,2)),'-m')
% ylim([-1.5 1.5]);
xlim([min(timeVector) max(timeVector)])
xlabel('Time [s]');
ylabel('Y [km]');


subplot(3,1,3);
plot(timeVector,xTrueMeas(ind,3),'r',timeVector,xCorrectedPF(ind,3),'b');
hold on
grid on
plot(timeVector,xCorrectedPF(ind,3)-3*sqrt(COV(ind,3)),'-m',timeVector,xCorrectedPF(ind,3)+3*sqrt(COV(ind,3)),'-m')
% ylim([-1.5 1.5]);
xlim([min(timeVector) max(timeVector)])
xlabel('Time [s]');
ylabel('Z[km]');

figure(3);
sgtitle('Velocity Measurement and predicitions')
subplot(3,1,1);
plot(timeVector,xTrueMeas(ind,4),'r',timeVector,xCorrectedPF(ind,4),'b');
hold on
grid on
plot(timeVector,xCorrectedPF(ind,4)-3*sqrt(COV(ind,4)),'-m',timeVector,xCorrectedPF(ind,4)+3*sqrt(COV(ind,4)),'-m')
legend('True','Particlte filter estimate','3 \sigma')
% ylim([-1 1]);
xlim([min(timeVector) max(timeVector)])
ylabel('$\dot{X}[km/s]$', 'Interpreter','latex');

subplot(3,1,2);
plot(timeVector,xTrueMeas(ind,5),'r',timeVector,xCorrectedPF(ind,5),'b');
hold on
grid on
plot(timeVector,xCorrectedPF(ind,5)-3*sqrt(COV(ind,5)),'-m',timeVector,xCorrectedPF(ind,5)+3*sqrt(COV(ind,5)),'-m')
% ylim([-1.5 1.5]);
xlim([min(timeVector) max(timeVector)])
xlabel('Time [s]');
ylabel('$\dot{Y}[km/s]$', 'Interpreter','latex');


subplot(3,1,3);
plot(timeVector,xTrueMeas(ind,6),'r',timeVector,xCorrectedPF(ind,6),'b');
hold on
grid on
plot(timeVector,xCorrectedPF(ind,6)-3*sqrt(COV(ind,6)),'-m',timeVector,xCorrectedPF(ind,6)+3*sqrt(COV(ind,6)),'-m')
% ylim([-1.5 1.5]);
xlim([min(timeVector) max(timeVector)])
xlabel('Time [s]');
ylabel('$\dot{Z}[km/s]$', 'Interpreter','latex');

%%
res = xTrueMeas(ind,:)-xCorrectedPF(ind,:);
figure(4);

plot(timeVector,xTrueMeas(ind,1)-xCorrectedPF(ind,1),'r','linewidth',2);
grid on
hold on




plot(timeVector,xTrueMeas(ind,2)-xCorrectedPF(ind,2),'g','linewidth',2);
grid on



plot(timeVector,xTrueMeas(ind,3)-xCorrectedPF(ind,3),'b','linewidth',2);
grid on
xlabel('Time [s]')
ylabel('Residual Error [km]')

% yline(100,'--r','linewidth',2);
% yline(-100,'--r','linewidth',2);
% ylim([-150 150])
ylim([-50 50])
xlim([timeVector(1) timeVector(end)])
title('Position Residuals vs time')
legend('R_x','R_y','R_z')%,'Allowed Error')



figure(5)

plot(timeVector,xTrueMeas(ind,4)-xCorrectedPF(ind,4),'r','linewidth',2);
grid on
hold on



plot(timeVector,xTrueMeas(ind,5)-xCorrectedPF(ind,5),'b','linewidth',2);
grid on


plot(timeVector,xTrueMeas(ind,6)-xCorrectedPF(ind,6),'g','linewidth',2);
grid on
title('Velocity Residuals vs time')
xlabel('Time [s]')
ylabel('Velocity Residual [km/s]', 'Interpreter','latex')
% yline(3,'--r','linewidth',2);
% yline(-3,'--r','linewidth',2);
V_leg = legend('$\dot{R_x}$', '$\dot{R_y}$', '$\dot{R_z}$');%,'Allowed Error');
set(V_leg,'Interpreter','latex');
ylim([-5 5])
xlim([timeVector(1) timeVector(end)])
if DoFilm == 1

vidRes = VideoWriter('ResidualsVTime.mp4','MPEG-4');
vidRes.FrameRate=5;
open(vidRes)
vidRes3D = VideoWriter('Residuals3D.mp4','MPEG-4');
vidRes3D.FrameRate=5;
open(vidRes3D)

for i = 1:length(timeVector)
    
    figure(8);

    plot(timeVector(1:i),res(1:i,1),'r','linewidth',2);
    grid on
    hold on

    plot(timeVector(1:i),res(1:i,2),'g','linewidth',2);
    grid on

    plot(timeVector(1:i),res(1:i,3),'b','linewidth',2);
    grid on
    xlabel('Time [s]')
    ylabel('Residual Error [km]')
   xline(timeVector(i),'linewidth',2);
    yline(100,'--r','linewidth',2);
    yline(-100,'--r','linewidth',2);
    ylim([-150 150])
    xlim([timeVector(1) timeVector(end)])
    title('Position Residuals vs time')
    legend('R_x','R_y','R_z','Time','Allowed Error')
    im2 = getframe(gcf);
    writeVideo(vidRes, im2)
    drawnow
    hold off
    
    figure(9);
    plot3(res(i,1),res(i,2),res(i,3),'*c','linewidth',1)
    hold on;
    plot3(0,0,0,'*m','linewidth',3);
    quiver3(0,0,0,res(i,1),0,0,'r','linewidth',1)
    quiver3(res(i,1),0,0,0,res(i,2),0,'g','linewidth',1)
    quiver3(res(i,1),res(i,2),0,0,0,res(i,3),'b','linewidth',1)
    quiver3(0,0,0,res(i,1),res(i,2),res(i,3),'k','linewidth',1)
    % round to closest five
    big = ceil(max(abs(res(i,1:3)))/5)*5;
    xlabel('X Position Error [km]')
    ylabel('Y Position Error [km]')
    zlabel('Z Position Error [km]')
    legend('PF Estimate','STK Satellite Location','R_x','R_y','R_z','R_t_o_t')
    title('Residual Vectors')
    xlim([-big big])
    ylim([-big big])
    zlim([-big big])
    grid on
    hold off;
    im3 = getframe(gcf);
    writeVideo(vidRes3D, im3)
    drawnow
end
close(vidRes3D)
close(vidRes)
end
res_mag = vecnorm(res(:,1:3),2,2);
tdoa_mag = vecnorm(TDoAres,2,2);
figure;
plot(timeVector,res_mag,'b','linewidth',2)
hold on 
grid on
plot(timeVector,tdoa_mag,'r','linewidth',2)
xlabel('Time [s]')
ylabel('Magnitude of Residual Error [km]')
ylim([-1 50])
xlim([timeVector(1) timeVector(end)])
legend('PF Residual', 'TDoA Residual')
end