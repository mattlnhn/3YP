function [rho_p] = densdist(beamxpos, dl, dx, dy, secxpos, secypos, phi, ...
                   np, nb, f, sigma_x, sigma_y)
%   DENSDIST proton density at each node
%   N.B. system origin at midpoint of leading (left) edge
%   dx, dy, secxpos, secypos MUST be multiples of dl
%   INPUTS ----------------------------------------------------------------
%   beamxpos    beam centre x pos wrt system origin [m]
%   dl          distance between nodes [m]
%   dx          size in x dimension [m]
%   dy          size in y dimension [m]
%   secxpos     x pos of bottom left corner of section wrt system origin
%               [m]
%   secypos     y pos of bottom left corner of section wrt system origin
%               [m]
%   np          protons per bunch
%   nb          bunches in machine
%   f           bunch revolution frequency [Hz]
%   sigma_x     beam core x width [m]
%   sigma_y     beam core y width [m]
%   OUTPUTS ---------------------------------------------------------------
%   rho_p
    xrange = (0:dl:dx-dl) + secxpos + dl/2;
    yrange = (dy-dl:-dl:0) + secypos + dl/2;
    xrange = xrange.*cosd(phi);

    xbeam = beamxpos;
    ybeam = 0;
    [xmesh, ymesh] = meshgrid(xrange, yrange);
    xrel2 = (xbeam - xmesh).^2;
    yrel2 = (ybeam - ymesh).^2;

    p = np*nb*f; % avg protons s^-1
    rho_p = (p/(2*pi*sigma_x*sigma_y))*...
            exp(-.5*((xrel2/sigma_x^2)+(yrel2/sigma_y^2)));
end