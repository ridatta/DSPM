function [P,V] = getPressureVelocity(posT,nrmT,posS,A,k,omega,rho)
% RETURNS THE UNSCALED DSPM ACOUSTIC PRESSURE and NORMAL VELOCITY AT
% TARGET POINTS posT
% posT = [Nx3] matrix conatining target points
% nrmT = [Nx3] vector conatining normal vectors
% posS = [Mx3] matrix containing source point positions
% A = [Mx1] Acoustic strength matrix
% k = wave number
% omega = angular freqeuncy
% rho = medium density
N = size(posT,1); % Number of target points
M = size(posS,1); % number of source pts
P = zeros(N,1); % pressure matrix
V = zeros(N,1); % normal velocity matrix
%%% CALCULATE PRESSURE %%%
parfor n = 1:N
    for m = 1:M
        r_mn = posT(n,:) - posS(m,:); % Position vector from source m to target pt n
        nrm_n = nrmT(n,:); % normal tp the surface at pt. n
        n_mn = dot(r_mn'/norm(r_mn), nrm_n'); % Projection of r_mn on z-axis
        %%% NORMAL VELOCITY %%%
        V(n) = V(n) + A(m) * 1 / (1i * omega * rho) * n_mn *...
            exp(-1i * k * norm(r_mn)) ...
            * (1i * k / norm(r_mn) + 1 / norm(r_mn)^2); 
        %%% ACOUSTIC PRESSURE %%%
        P(n) = P(n) + A(m) * exp(-1i * k * norm(r_mn)) / norm(r_mn); 
    end
end
end