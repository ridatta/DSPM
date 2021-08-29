clc; close all; clear;

f = 1;
wl = 1;
A = 1;
k = 2 * pi / wl;
omega = 2 * pi * f;
rho = 1;

r = linspace(0,k*wl,100);

t = linspace(0,1/f,10);

figure
for ii = 1:length(t)
    p = A .* exp(1i .* (k .* r - omega .* t(ii))) ./ r; 
    v = A ./ (rho .* 1i .* omega) .* (1i .* k ./ r - 1./ r.^2) .* ...
    (exp(1i .* (k .* r - omega .* t(ii))));
    plot(r,real(p),'r'); hold on;
    plot(r,real(v),'b');
    hold off;
    pause(1);
end


