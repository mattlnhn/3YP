function [rho_p] = densdist(pos, n, L, np, f, sigma_x, sigma_y)
% DENSDIST proton density at each node
%   pos         position of beam centre relative to midpoint of left edge of
%               detector, 1x2 vector
%   n           no. of nodes
%   L           detector side length [m]
%   np          no. protons per bunch
%   f           bunches per second [s-1]
%   sigma_xy    beam core width (1 s.d.) [m]

    dl = L/n;
    xrange = (0:dl:L-dl) + dl/2;
    yrange = (0:dl:L-dl) + dl/2 - L/2;
    xbeam = pos(1);
    ybeam = pos(2);
    [xmesh, ymesh] = meshgrid(xrange, yrange);
    % node distance from beam centre squared
    %d2 = (xpos - xmesh).^2 + (ypos - ymesh).^2;
    xrel2 = (xbeam - xmesh).^2;
    yrel2 = (ybeam - ymesh).^2;

    % gaussian
    p = np*f; % avg protons s^-1
    rho_p = (p/(2*pi*sigma_x*sigma_y))*exp(-.5*((xrel2/sigma_x^2)+(yrel2/sigma_y^2)));

end

