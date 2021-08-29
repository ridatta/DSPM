function [P,posT] = getPressureField(xrng,zrng,posS,A,k)
% RETURNS THE UNSCALED DSPM ACOUSTIC PRESSURE FIELD FOR THE XZ PLANE
% xrng, zrng = x and z field ranges
% posS = [Mx3] matrix containing source point positions
% A = [Mx1] Acoustic strength matrix
% k = wave number
%%% DESCRITIZE FIELD %%%
res = 300;
X = linspace(xrng(1),xrng(end),res);
Z = linspace(zrng(1),zrng(end),res);
[XX,ZZ] = meshgrid(X,Z);
posT = [XX(:), 0 .* XX(:), ZZ(:)];
posT = unique(posT,'rows'); % Target points
%%% CALCULATE PRESSURE %%%
N = size(posT,1); % Number of target points
M = size(posS,1); % number of source pts
P = zeros(N,1); % pressure matrix
parfor n = 1:N
    for m = 1:M
        r_mn = posT(n,:) - posS(m,:); % Position vector from source m to target pt n
        P(n) = P(n) + A(m) * exp(-1i * k * norm(r_mn)) / norm(r_mn); 
    end
end
end