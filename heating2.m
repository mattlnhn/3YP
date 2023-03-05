function [dT] = heating2(n, L, theta, dt, nt, mat, beam)
% HEATING
%   INPUTS ----------------------------------------------------------------
%   n       no. of nodes
%   L       length of side of detector [m]
%   theta   detector thickness [m]
%   dt      time step [s]
%   nt      no. of time steps
%   mat     material properties structure
%           mat.A           atomic mass, g mol-1
%           mat.Z           atomic no.
%           mat.I           mean ionisation energy, MeV
%           mat.w           mass fraction
%           mat.rho         density, g cm-3
%           mat.epsilon     emissivity
%           mat.k           thermal conductivity, W m-1 K-1
%           mat.c_p         specific heat capacity, J kg-1 K-1
%   beam    beam properties structure
%           beam.M          particle mass, MeV c-2
%           beam.gamma      Lorentz factor
%           beam.f          bunches per second
%           beam.np         protons per bunch
%           beam.sigma      beam width, m
%           beam.pos        beam centre position w.r.t. midpoint of 
%                           detector left edge, 1x2, m
%   OUTPUT ----------------------------------------------------------------
%   T       n*n*nt matrix of temp at each node at each time step
%   NOTES -----------------------------------------------------------------
%   bethe function uses combination of CGS and eV for units; heat transfer
%   working in SI units hence unit conversion AFTER bethe function has been
%   called
    
    %% CONSTANTS
    const.sigma = 5.670374419e-8; % stefan-boltzmann, W m-2 K-4
    
    %% MESH
    dl = L/n; % volume cell side length
    dT = zeros(n, n, nt+1); % n*n mesh at each t inc. t=0
    T_0 = 1.9*ones(n, n); % K
    %T(:, :, 1) = ones(n, n)*T_0; % T(0)
    E.L = [ones(n, 1) zeros(n, n-1)]; % left edge
    E.T = [ones(1, n); zeros(n-1, n)]; % top
    E.R = [zeros(n, n-1) ones(n, 1)]; % right
    E.B = [zeros(n-1, n); ones(1, n)]; % bottom

    %% ENERGY DEPOSITION
    beta = (1-beam.gamma^-2)^.5;
    % mass stopping power = <1/rho dE/dx>
    msp = bethe(beam.M, mat.A, mat.Z, mat.I, mat.w, mat.rho, beta, 1);
    rho_p = densdist(beam.pos, n, L, beam.np, beam.f, beam.sigma, beam.sigma);

    %% CONVERT UNITS
    % mass stopping power MeV g-1 cm2 --> J g-1 m2
    % * joules per mev * (m per cm)^2
    msp = msp * 1.6021773e-13 * (0.1)^2;
    % density g cm-3 --> kg m-3
    mat.rho = mat.rho * 1000;

    
    %% TIME ITERATION
    i = 1;

    % coefficients
    C1 = (mat.k*dt)/(mat.rho*mat.c_p*dl*dl); % inner node
    C2 = (dt)/(mat.rho*mat.c_p*dl*dl*theta); % boundary node
    %C3 = (dt*msp)/(mat.rho*mat.c_p*mat.c_p*dl*dl*theta) % protons
    C3 = (dt*msp)/mat.c_p;
    C4 = (2*mat.epsilon*const.sigma*dt)/(mat.rho*mat.c_p*theta); % rad

    while i <= nt
        % PREVIOUS T is T(:, :, i)
        % NEW T is T(:, :, i+1)
        
        % moving beam centre relative to edge
        % rho_p = ;

        % boundary conditions
        BQ.L = 0; % boundary Qdot, left
        BQ.T = 0; % top
        BQ.R = 0; % right
        BQ.B = 0; % bottom
        
        % temp matrices
        Ti.N = T_0 + dT(:, :, i); % node
        Ti.L = [zeros(n, 1) Ti.N(:, 1:n-1)]; % left
        Ti.T = [zeros(1, n); Ti.N(1:n-1, :)]; % top
        Ti.R = [Ti.N(:, 2:n) zeros(n, 1)]; % right
        Ti.B = [Ti.N(2:n, :); zeros(1, n)]; % bottom

        % next T
        T = Ti.N + ...
            C1*(Ti.L-Ti.N).*not(E.L) + C2*BQ.L.*E.L + ...
            C1*(Ti.T-Ti.N).*not(E.T) + C2*BQ.T.*E.T + ...
            C1*(Ti.R-Ti.N).*not(E.R) + C2*BQ.R.*E.R + ...
            C1*(Ti.B-Ti.N).*not(E.B) + C2*BQ.B.*E.B + ...
            C3*rho_p + ...
            -C4*(Ti.N.^4-T_0^4);

        dT(:, :, i+1) = T - T_0;

        i = i + 1;

    end

end

