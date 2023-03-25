function [rho_p] = densdist(pos, n, L, xfrac, np, nb, f, sigma_x, sigma_y)
% DENSDIST proton density at each node
%   pos         x position of beam centre relative to midpoint of left edge of
%               detector
%   n           no. of nodes
%   L           detector side length [m]
%   xfrac       fraction of detector in x-direction
%   np          no. protons per bunch
%   f           bunches per second [s-1]
%   sigma_xy    beam core width (1 s.d.) [m]

    dl = L/n;
    xrange = (0:dl:L-dl) + dl/2;
    xrange = xrange(:, 1:xfrac*n);
    yrange = ((L/2)-dl:-dl:0) + dl/2;
    xbeam = pos;
    ybeam = 0;
    [xmesh, ymesh] = meshgrid(xrange, yrange);
    % node distance from beam centre squared
    %d2 = (xpos - xmesh).^2 + (ypos - ymesh).^2;
    xrel2 = (xbeam - xmesh).^2;
    yrel2 = (ybeam - ymesh).^2;

    % gaussian
    p = np*f*nb; % avg protons s^-1
    rho_p = (p/(2*pi*sigma_x*sigma_y))*exp(-.5*((xrel2/sigma_x^2)+(yrel2/sigma_y^2)));

end

