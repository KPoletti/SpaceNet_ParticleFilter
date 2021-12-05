function[] = Plot_PF_weight(pf,xTrueMeas,W)
    ind = find(W>1e-4);
    ind = sort([ind ; randi([1 10^4],2*length(ind),1)]);
    pf = pf(:,ind);
    W = W(ind);
    figure(2)
%     subplot(1,2,2)
    scatter3(pf(1,:)-xTrueMeas(:,1),pf(2,:)-xTrueMeas(:,2),pf(3,:)-xTrueMeas(:,3),[],W)
    grid on
    hold on
    plot3(0,0,0,'*m','linewidth',3);
    xlim([-50 50])
    ylim([-50 50])
    zlim([-50 50])
    legend('Particles','STK Satellite Position','location','northeast')
    xlabel('X Position Error [km]')
    ylabel('Y Position Error [km]')
    zlabel('Z Position Error [km]')
    title('Zoomed in View with STK Satellite Position at the Origin')
    hold off
    
end