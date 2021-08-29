function  scaleFactorStudy(fname)
load(fname);
ff = [1, 0.9, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05]; 
vel_avg = []; N = [];
for ii = 1:length(ff)
[posTT,nrmTT] = descritize(R0,ff(ii)*pS,L0,dir0,des_typ); 
N = [N; size(posTT,1)];
[~,vel] = getPressureVelocity(posTT,nrmTT,posS,A,k,omega,rho);
vel_avg = [vel_avg; mean(abs(vel))]; % m/s, average normal surface velocity
end
[N,idx] = sort(N);
vel_avg = vel_avg(idx); 
save('scaleFactorStudy.mat'); 
figure
plot(log(N), vel_avg,'k^','LineWidth',1.5);
xlabel('log(Number of points)');
ylabel('Average velocity, m/s'); 
formatPlots();
end