clc; close all; clear;

% DSPM Implementation based on Cheng & Qin (2011)
% Geometry = Circular Transducer D = 12.7 mm

% Parameters
%%% ACOUSTIC PROPERTIES %%%
f = 27.44*10^3; % hz, frequency
rho = 1.3; % kg/m^3, medium density
c = 340; % m/s, speed of sound
v0 = 0.5; % m/s, amplitude of surface velocity
omega = 2 * pi * f; % angular frequency, rad/s
k = 2 * pi * f / c; % wave number, m^-1
wl = 2 * pi / k; % wavelength
%%% DROPLET PROPERTIES %%%
z0 = 2.2 * 10^-3; % drop position, [m]
rd = 0.79*10^-3; % droplet radius, m
rho_f = 998; % density of drop, kg/m^3
c_f = 1540; % drop sound velocity, m/s 
Zf = rho_f * c_f; % drop impedance
isRigid = false; % is the drop a rigid body
%%% DSPM/GEOMETRY PARAMETERS %%%
rs = wl / 8; % m, source radius
pS = 2 * rs; % m, source pitch distance
R0 = 0.5 * 30 * 10^-3; % m, Transducer radius
L0 = 12.7*10^-3; % Transducer location, m 
dir0 = [0;0;-1]; % transducer direction vector
R_ref = 30 * 10^-3; % Reflector radius, [m]
L_ref = 0; % reflector location, m 
dir_ref = [0;0;1]; % reflector direction vector
des_typ = 'hex'; % descretization scheme

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CIRCULAR TRANSDUCER GEOMETRY %

[posAT,nrmAT] = descritize(R0,pS,L0,dir0,des_typ); 
posAS = posAT; posAS(:,3) = L0 -1 * dir0(end) * rs; % transducer source points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% REFLECTOR GEOMETRY %%%
refOn = true;
if refOn 
    [posR,nrmR] = descritize(R_ref,pS,L_ref,dir_ref,des_typ);
    posRS = posR; posRS(:,3) = L_ref + -1 * dir_ref(end) * rs; % reflector source points
else
    posR = []; posRS = []; nrmR = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DROP GEOMETRY %%%
dropOn = false;
if dropOn
    res = 9;
    [xx,yy,zz] = sphere(res);
    xx = xx * rd; yy = yy * rd; zz = zz * rd + z0;
    posD = [xx(:), yy(:), zz(:)];
    posD = unique(posD,'rows');
    rds = rd - 0.2*rd; % sphere radius, m
    [xx,yy,zz] = sphere(res);
    xx = xx * rds; yy = yy * rds; zz = zz * rds + z0;
    posDS = [xx(:), yy(:), zz(:)];
    posDS = unique(posDS,'rows'); % source position
    % Normal vectors
    C = [0, 0, z0]; % position vector of droplet center, [m]
    nrmD = posD - C; % vector from droplet center to surface, [m]
    nrmD_len = sqrt(nrmD(:,1) .^ 2 + nrmD(:,2) .^ 2 + nrmD(:,3) .^ 2); % norm
    nrmD = nrmD ./ nrmD_len; % normalize
else
    posD = []; posDS = []; nrmD = [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatenate all source and target points
posS = [posAS; posRS; posDS];
posT = [posAT; posR; posD]; 
nrmT = [nrmAT; nrmR; nrmD]; 
numTR = size(posAT,1); % Number of transducer points
numR = size(posR,1); % Number of reflector points
numD = size(posD,1); % Number of drop points
M = size(posS,1); % Number of source points
N = size(posT,1); % Number of target points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
% Create elements of MM matrix
MM = zeros(N,M); % M matrix
parfor m = 1:M
    for n = 1:N
        nrm_n = nrmT(n,:); % normal tp the surface at pt. n
        r_mn = posT(n,:) - posS(m,:); % Position vector from source m to n
        n_mn = dot(r_mn'/norm(r_mn), nrm_n'); % Projection of r_mn on normal to surface
        if n <= (numTR + numR) % if target point does not lie on surface of droplet
        MM(n,m) = 1 / (1i * omega * rho) * n_mn * exp(-1i * k * norm(r_mn)) * ...
            (1i * k / norm(r_mn) + 1 / norm(r_mn)^2); % element in row n and column m
        else
        if ~(isRigid)
        MM(n,m) = 1 / (1i * omega * rho) * n_mn * exp(-1i * k * norm(r_mn)) ...
            * (1i * k / norm(r_mn)+ 1 / norm(r_mn)^2) - (1 / norm(r_mn)) * exp(-1i * k * norm(r_mn)) / Zf; % element in row n and column m
        else
            MM(n,m) = 1 / (1i * omega * rho) * n_mn * exp(-1i * k * norm(r_mn)) * ...
            (1i * k / norm(r_mn) + 1 / norm(r_mn)^2);  % element in row n and column m
        end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%% BOUNDARY CONDITIONS %%%%
V = zeros(N,1); % Velocity, m/s
V(1:numTR) = v0; % Rigid transducer surface; Velocity at transducer surface 
V(numTR+1:numTR+numR) = 0; % Rigid reflector
V(numTR+numR+1:end) = 0; % Droplet
%%% SOLVE FOR ACOUSTIC STRENGTHS %%%%
A = pinv(MM) * V; % Acoustic strength
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCALE-FACTOR %%%
% Find normal surface velocity and scaling factor
% Get coordinates of transducer surface points with 10x finer pitch
[posTT,nrmTT] = descritize(R0,0.1*pS,L0,dir0,des_typ); 
[~,vel] = getPressureVelocity(posTT,nrmTT,posS(1:numTR,:),A(1:numTR),k,omega,rho);
vel_avg = mean(abs(vel)); % m/s, average normal surface velocity
scale_factor = v0 / vel_avg; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ACOUSTIC PRESSURE FIELD %%%
% Define field points and calculate pressure
[press, fldPts] = getPressureField([-L0,L0],sort([L0,L_ref]),posS,A,k); 
press_scaled = scale_factor * press;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ANTINODAL PRESSURE %%%
AN = [0,0,L0+dir0(end)*0.5*wl];
[antiNodalPressure,~] = getPressureVelocity(AN,[0,1,0],posS,A,k,omega,rho);
antiNodalPressure = scale_factor * antiNodalPressure; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ACOUSTIC RADIATION PRESSURES %%%
% Get first-order pressures and velocities
[p1,u1] = getPressureVelocity(posD,nrmD,posS,A,k,omega,rho);
p1 = scale_factor * p1; % Scale pressure field
u1 = scale_factor * u1; % Scale velocity field
if isRigid
    err = mean(abs(u1));
else
err = mean(abs(u1) - abs(p1) / Zf);
end
% Get RADIATION PRESSURE
p_rad = (p1).^2 / (2 * rho * c^2) - rho * abs(u1) .^ 2 / 2; 
% Projection in z-direction
PZ = zeros(1,length(p_rad)); PX = zeros(1,length(p_rad));
PY = zeros(1,length(p_rad));
for ii = 1:length(p_rad)
    ez = [0;0;1];ex = [1;0;0];ey = [0;1;0];
    n_vect = nrmD(ii,:); 
    PZ(ii) = p_rad(ii) * dot(n_vect',ez);
    PX(ii) = p_rad(ii) * dot(n_vect',ex);
    PY(ii) = p_rad(ii) * dot(n_vect',ey);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUT LOG %%%     
fname = ['Wada-2011-DSPM-Output-' datestr(clock,29)];
save([fname, '.mat']); 

DSPMlogger(fname);
DSPMvisualizer(fname);



                    