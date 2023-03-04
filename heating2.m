function [T] = heating2(n, L, theta, dt, nt, mat, beam)
% HEATING
%   INPUTS ----------------------------------------------------------------
%   n       no. of nodes
%   L       length of side of detector [m]
%   theta   detector thickness [m]
%   dt      time step [s]
%   nt      no. of time steps
%   mat     material properties structure
%           mat.A           atomic mass
%           mat.Z           atomic no.
%           mat.I           mean ionisation energy
%           mat.w           mass fraction
%           mat.rho         density
%           mat.epsilon     emissivity
%           mat.k           thermal conductivity
%           mat.c_p         specific heat capacity
%   beam    beam properties structure
%           beam.M          particle mass [MeV c-2]
%           beam.gamma      Lorentz factor
%           beam.f          bunches per second
%           beam.np         protons per bunch
%           beam.sigma      beam width
%           beam.pos        beam centre position w.r.t. midpoint of 
%                           detector left edge
%   OUTPUT ----------------------------------------------------------------
%   T       n*n*nt matrix of temp at each node at each time step
    
    %% CONSTANTS
    const.sigma = ; % stefan-boltzmann
    
    %% MESH
    dl = L/n; % volume cell side length
    T = zeros(n, n, nt+1); % n*n mesh at each t inc. t=0
    T_0 = 1.9; % K
    T(:, :, 1) = ones(n, n)*T_0; % T(0)
    E.L = [ones(n, 1) zeros(n, n-1)]; % left edge
    E.T = [ones(1, n); zeros(n-1, n)]; % top
    E.R = [zeros(n, n-1) ones(n, 1)]; % right
    E.B = [zeros(n-1, n); ones(1, n)]; % bottom
    

    %% ENERGY DEPOSITION
    beta = (1-beam.gamma^-2)^.5;
    % mass stopping power = <1/rho dE/dx>
    msp = bethe(beam.M, mat.A, mat.Z, mat.I, mat.w, mat.rho, beta, 1);
    rho_p = densdist(beam.pos, n, L, beam.np, beam.f, beam.sigma);
    
    %% TIME ITERATION
    i = 1;

    % coefficients
    C1 = (mat.k*dt)/(mat.rho*mat.c_p*dl*dl); % inner node
    C2 = (dt)/(mat.rho*mat.c_p*dl*dl*theta); % boundary node
    C3 = (dt*msp)/(mat.rho*mat.c_p*mat.c_p*dl*dl*theta); % protons
    C4 = (2*mat.epsilon*const.sigma*dt)/(mat.rho*mat.c_p*theta); % rad

    while i <= nt
        % CURRENT T is T(:, :, (i+1))
        % NEXT T is T(:, :, (i+1)+1) = T(:, :, i+2)
        % offset by one due to MATLAB indexing and t=0 being first frame
        
        % moving beam centre relative to edge
        % rho_p = ;

        % boundary conditions
        BQ.L = 0; % boundary Qdot, left
        BQ.T = 0; % top
        BQ.R = 0; % right
        BQ.B = 0; % bottom
        
        % temp matrices
        Ti.N = T(:, :, i+1); % node
        Ti.L = [zeros(n, 1) Ti.N(:, 1:n-1)]; % left
        Ti.T = [zeros(1, n); Ti.N(1:n-1, :)]; % top
        Ti.R = [Ti.N(:, 2:n) zeros(n, 1)]; % right
        Ti.B = [Ti.N(2:n, :); zeros(1, n)]; % bottom

        % next T
        T(:, :, i+2) = Ti.N + ...
            C1*(Ti.L-Ti.N).*not(E.L) + C2*BQ.L.*E.L + ...
            C1*(Ti.T-Ti.N).*not(E.T) + C2*BQ.T.*E.T + ...
            C1*(Ti.R-Ti.N).*not(E.R) + C2*BQ.R.*E.R + ...
            C1*(Ti.B-Ti.N).*not(E.B) + C2*BQ.B.*E.B + ...
            C3*rho_p + ...
            -C4*(Ti.N.^4-T_0^4);

        i = i + 1;

    end

end

%% FUNCTIONS

function [rho_p] = densdist(pos, n, L, np, f, sigma_xy)
% DENSDIST proton density at each node
%   pos         position of beam centre relative to midpoint of left edge of
%               detector, 1x2 vector
%   n           no. of nodes
%   L           detector side length [m]
%   np          no. protons per bunch
%   f           bunches per second [s-1]
%   sigma_xy    beam core width (1 s.d.) [m]

    dl = L/n;
    x = (0:dl:L-dl) + dl/2;
    y = (0:dl:L-dl) + dl/2 - L/2;
    xpos = pos(1);
    ypos = pos(2);
    [xmesh, ymesh] = meshgrid(x, y);
    % node distance from beam centre squared
    d2 = (xpos - xmesh).^2 + (ypos - ymesh).^2;
    
    % gaussian
    p = np*f; % avg protons s^-1
    rho_p = p*(2*pi*sigma_xy^2)^-1*exp(-.5*d2/sigma_xy^2);

end
