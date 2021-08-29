function DSPMvisualizer(fname)
% DSPM VISUALISER %
% Reads data from DSPM output .mat file and plots
% (1) GEOMETRY
% (2) X-Z PRESSURE FIELD
% (3) PRESSURE-FIELD + GEOMETRY
% (4) RADIATION PRESSURE

load([fname, '.mat']); 
%%% (1) PLOT GEOMETRY %%%
% Plot sources and transducer surface points (target points)
figure
plot3(posS(:,1)*1000,posS(:,2)*1000,posS(:,3)*1000,'rO','LineWidth',1.5); hold on;
hold on;
plot3(posT(:,1)*1000,posT(:,2)*1000,posT(:,3)*1000,'bx','LineWidth',1.5); hold on;
legend('Point Source','Target Point');
xlabel('x-position, mm');
ylabel('y-position, mm');
zlabel('z-position, mm');
formatPlots();
axis equal;
%%% (2) PLOT X-Z PRESSURE FIELD %%%
% Plot sources and transducer surface points (target points)
figure
sz = 200;
scatter(fldPts(:,1)*10^3,fldPts(:,3)*10^3,...
    sz,abs(press_scaled),'s','filled');
xlabel('$x, mm$','Interpreter','latex');
ylabel('$z, mm$','Interpreter','latex');
colormap('hot'); cb = colorbar;
ylabel(cb,'{\it Pressure, Pa}','Fontname','Times New Roman','Fontsize',18);
formatPlots();
%ylim(1000*[0,L0]);xlim(1000*[-L0/2,L0/2]);
axis square;

%%% (3) ACOUSTIC PRESSURE FIELD + GEOMETRY %%%%%
% Plot sources and transducer surface points (target points)
figure
plot3(posS(:,1)*1000,posS(:,2)*1000,posS(:,3)*1000,'rO','LineWidth',1.5); hold on;
hold on;
plot3(posT(:,1)*1000,posT(:,2)*1000,posT(:,3)*1000,'bx','LineWidth',1.5); hold on;
scatter3(fldPts(:,1)*10^3,fldPts(:,2)*10^3,fldPts(:,3)*10^3,...
    sz,(abs(press_scaled)),'s','filled');
xlabel('$x, mm$','Interpreter','latex');
ylabel('$y, mm$','Interpreter','latex');
zlabel('$z, mm$','Interpreter','latex');
colormap('hot'); cb = colorbar;
ylabel(cb,'{\it Pressure, Pa}','Fontname','Times New Roman','Fontsize',18);
legend('Point Source','Target Point', 'Pressure Field','Location','northwest');
formatPlots();
axis equal;

if dropOn
%%% (4) ACOUSTIC RADIATION PRESSURE %%%%%
figure
scatter3(posD(:,1)*10^3,posD(:,2)*10^3,posD(:,3)*10^3,...
    2500,abs(p_rad),'filled');
xlabel('x-position, mm','Interpreter','latex');
ylabel('y-position, mm','Interpreter','latex');
zlabel('z-position, mm','Interpreter','latex'); hold on;
colormap('jet'); cb = colorbar;
ylabel(cb,'{\it Radiation Pressure, Pa}','Fontname','Times New Roman','Fontsize',18);
view(90,0);
axis('equal');
formatPlots();
end
end
