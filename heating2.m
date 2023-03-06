function [dTframes] = heating2(n, L, theta, dt, nt, mat, beam)
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

    %% PHYSICAL CONSTANTS
    const.sigma = 5.670374419e-8; % stefan-boltzmann, W m-2 K-4
    
    %% MESH
    dl = L/n; % volume cell side length
    dT = zeros(n, n); % preallocate
    T_0 = 1.9*ones(n, n); % K
    T_04 = T_0.^4; % used in loop

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

    %% ANIMATION OUTPUT PARAMS
    fps = 100; % frames per simulated second
    nf = round((dt*fps)^-1); % no. of time steps btwn frames
    frames = round(nt/nf); % no. of frames total
    dTframes = zeros(n, n, frames); % preallocate

    %% TIME ITERATION
    i = 1;

    % coefficients
    C1 = (mat.k*dt)/(mat.rho*mat.c_p*dl*dl); % inner node
    C2 = (dt)/(mat.rho*mat.c_p*dl*dl*theta); % boundary node
    % if rho_p constant w.r.t. t
    C3 = (dt*msp*rho_p)/mat.c_p; % protons
    % if not, include rho_p term in calculation later
    %C3 = (dt*msp)/mat.c_p; % protons
    C4 = (2*mat.epsilon*const.sigma*dt)/(mat.rho*mat.c_p*theta); % rad

    % if boundary conditions stay constant, calculate here
    BQ.L = 0; % boundary Qdot, left
    BQ.T = 0; % top
    BQ.R = 0; % right
    BQ.B = 0; % bottom

    % combine before loop to save computation time
    C1L = C1*not(E.L);
    C1T = C1*not(E.T);
    C1R = C1*not(E.R);
    C1B = C1*not(E.B);
    C2L = C2*E.L;
    C2T = C2*E.T;
    C2R = C2*E.R;
    C2B = C2*E.B;
    
    while i <= nt

        %fprintf('Current time is t+%d\n%.4f%% complete\n\n', i*dt, 100*i/nt)
        
        % PREVIOUS T is T(:, :, i)
        % NEW T is T(:, :, i+1)
        
        % moving beam centre relative to edge
        %rho_p = ;

        % if boundary conditions change at each step, calculate here
        %BQ.L = 0; % boundary Qdot, left
        %BQ.T = 0; % top
        %BQ.R = 0; % right
        %BQ.B = 0; % bottom
        
        % temp matrices
        Ti.N = T_0 + dT; % node

        Ti.L = [zeros(n, 1) Ti.N(:, 1:n-1)]; % left
        Ti.T = [zeros(1, n); Ti.N(1:n-1, :)]; % top
        Ti.R = [Ti.N(:, 2:n) zeros(n, 1)]; % right
        Ti.B = [Ti.N(2:n, :); zeros(1, n)]; % bottom

        % next T
        T = Ti.N + ...
            C1L.*(Ti.L-Ti.N) + C2L.*BQ.L + ...
            C1T.*(Ti.T-Ti.N) + C2T.*BQ.T + ...
            C1R.*(Ti.R-Ti.N) + C2R.*BQ.R + ...
            C1B.*(Ti.B-Ti.N) + C2B.*BQ.B + ...
            C3 + ...
            -C4*(Ti.N.^4-T_04);

        dT = T - T_0;

        if rem(i, nf) == 0
            fprintf('%.4f%% complete. Animation frame captured.\n', 100*i/nt)
            dTframes(:, :, i/nf) = dT;
        end

        i = i + 1;

    end

end

