function [pos,nrm] = descritize(R0,pS,L,n,des_typ)
% Descritizes a circular surface
% R0 = surface radius, [m]
% pS = pitch distance, [m]
% L = z-position of surface
% n = normal vector to surface [3x1]

if strcmpi(des_typ,'rect')
    x = 0:pS:R0;
    x = [-x, x];
    [XX,YY] = meshgrid(x,x);
    pos = [XX(:), YY(:)];
    pos = unique(pos,'rows');
    pos = [pos, L*ones(size(pos,1),1)]; % point position, m
    % REMOVE EXCESS POINTS
    d = sqrt((pos(:,1)).^2 + (pos(:,2)).^2); % distance from the origin (x0,y0), m
    % remove points where d > transducer radius
    idx = find(d >  R0);
    pos(idx,:) = []; % Posiution matrix of points [N x 3]
elseif strcmpi(des_typ,'hex')
    % Get co-ordinates of point sources
    R = 0; % Radial dist. of points in row, [m]
    count = 0; % Counts current row
    X = []; Y = [];
    while R <= R0
        if count == 0
            num = 1;
        else
            num = 6 * count; % Number of points in row
        end
        theta = 0:2*pi/num:2*pi*(1-1/num);% Angular placement of each transducer
        theta = theta';
        X = [X; R * cos(theta)];
        Y = [Y; R * sin(theta)]; 
        R = R + pS; count = count + 1; 
    end
    pos = [X, Y, L*ones(size(X))]; % Coordinates of point sources
    pos = unique(pos,'rows'); 
end
 nrm = repmat(n',[size(pos,1),1]); % Matrix of normal vectors [N x 3]
end