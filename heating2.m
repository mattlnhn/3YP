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

    % adaptive grid layout
    %{

    x f n n c c c c     ]           x = xfine
    x f n n c c c c     |           f = fine
    x f n n c c c c     |           n = normal
    x f n n c c c c      >  L       c = coarse
    x f n n c c c c     |
    x f n n c c c c     |
    x f n n c c c c     |
    x f n n c c c c     ]

    %}

    nc = n
    nn = n*2
    nf = n*2*2
    nx = n*2*2*2

    wc = nc/2; % fractions must add to 1
    wn = nn/4;
    wf = nf/8;
    wx = nx/8;

    % xfine
    dl.x = L/nx;
    dT.x = zeros(nx, wx);
    T0.x = 1.9*ones(nx, wx);
    T04.x = 1.9^4*ones(nx, wx);
    E.x.l = dT.x;
    E.x.l(:, 1) = 1;
    E.x.t = dT.x;
    E.x.t(1, :) = 1;
    E.x.r = dT.x;
    E.x.r(:, end) = 1;
    E.x.b = dT.x;
    E.x.b(end, :) = 1;

    % fine
    dl.f = L/nf;
    dT.f = zeros(nf, wf);
    T0.f = 1.9*ones(nf, wf);
    T04.f = 1.9^4*ones(nf, wf);
    E.f.l = dT.f;
    E.f.l(:, 1) = 1;
    E.f.t = dT.f;
    E.f.t(1, :) = 1;
    E.f.r = dT.f;
    E.f.r(:, end) = 1;
    E.f.b = dT.f;
    E.f.b(end, :) = 1;

    % normal
    dl.n = L/nn;
    dT.n = zeros(nn, wn);
    T0.n = 1.9*ones(nn, wn);
    T04.n = 1.9^4*ones(nn, wn);
    E.n.l = dT.n;
    E.n.l(:, 1) = 1;
    E.n.t = dT.n;
    E.n.t(1, :) = 1;
    E.n.r = dT.n;
    E.n.r(:, end) = 1;
    E.n.b = dT.n;
    E.n.b(end, :) = 1;

    % coarse
    dl.c = L/nc;
    dT.c = zeros(nc, wc);
    T0.c = 1.9*ones(nc, wc);
    T04.c = 1.9^4*ones(nc, wc);
    E.c.l = dT.c;
    E.c.l(:, 1) = 1;
    E.c.t = dT.c;
    E.c.t(1, :) = 1;
    E.c.r = dT.c;
    E.c.r(:, end) = 1;
    E.c.b = dT.c;
    E.c.b(end, :) = 1;
    

    %% ENERGY DEPOSITION
    beta = (1-beam.gamma^-2)^.5;
    % mass stopping power = <1/rho dE/dx>
    msp = bethe(beam.M, mat.A, mat.Z, mat.I, mat.w, mat.rho, beta, 1);

    rho_p.x = densdist(beam.pos, nx, L, 1/8, beam.np, beam.f, beam.sigma, beam.sigma);
    rho_p.f = densdist(beam.pos-L/8, nf, L, 1/8, beam.np, beam.f, beam.sigma, beam.sigma);
    rho_p.n = densdist(beam.pos-L/4, nn, L, 1/4, beam.np, beam.f, beam.sigma, beam.sigma);
    rho_p.c = densdist(beam.pos-L/2, nc, L, 1/2, beam.np, beam.f, beam.sigma, beam.sigma);


    %% CONVERT UNITS
    % mass stopping power MeV g-1 cm2 --> J g-1 m2
    % * joules per mev * (m per cm)^2
    msp = msp * 1.6021773e-13 * (0.1)^2;
    % density g cm-3 --> kg m-3
    mat.rho = mat.rho * 1000;


    %% ANIMATION OUTPUT PARAMS
    fps = 10; % frames per simulated second
    nframes = round((dt*fps)^-1); % no. of time steps btwn frames
    frames = round(nt/nframes); % no. of frames total
    dTframes = zeros(nx, nx, frames); % preallocate


    %% TIME ITERATION
    i = 1;

    % coefficients
    c1.x.coef = (mat.k*dt)/(mat.rho*mat.c_p*dl.x*dl.x); % inner node
    c1.f.coef = (mat.k*dt)/(mat.rho*mat.c_p*dl.f*dl.f);
    c1.n.coef = (mat.k*dt)/(mat.rho*mat.c_p*dl.n*dl.n);
    c1.c.coef = (mat.k*dt)/(mat.rho*mat.c_p*dl.c*dl.c);
    c2.x.coef = (dt)/(mat.rho*mat.c_p*dl.x*dl.x*theta); % boundary node
    c2.f.coef = (dt)/(mat.rho*mat.c_p*dl.f*dl.f*theta);
    c2.n.coef = (dt)/(mat.rho*mat.c_p*dl.n*dl.n*theta);
    c2.c.coef = (dt)/(mat.rho*mat.c_p*dl.c*dl.c*theta);
    % if rho_p constant w.r.t. t
    c3.x = (dt*msp*rho_p.x)/mat.c_p; % protons
    c3.f = (dt*msp*rho_p.f)/mat.c_p;
    c3.n = (dt*msp*rho_p.n)/mat.c_p;
    c3.c = (dt*msp*rho_p.c)/mat.c_p;
    % if not, include rho_p term in calculation later
    %c3 = (dt*msp)/mat.c_p; % protons
    c4 = (2*mat.epsilon*const.sigma*dt)/(mat.rho*mat.c_p*theta); % rad


    c1.x.l = c1.x.coef*not(E.x.l);
    c1.x.t = c1.x.coef*not(E.x.t);
    c1.x.r = c1.x.coef*not(E.x.r);
    c1.x.b = c1.x.coef*not(E.x.b);
    c1.f.l = c1.f.coef*not(E.f.l);
    c1.f.t = c1.f.coef*not(E.f.t);
    c1.f.r = c1.f.coef*not(E.f.r);
    c1.f.b = c1.f.coef*not(E.f.b);
    c1.n.l = c1.n.coef*not(E.n.l);
    c1.n.t = c1.n.coef*not(E.n.t);
    c1.n.r = c1.n.coef*not(E.n.r);
    c1.n.b = c1.n.coef*not(E.n.b);
    c1.c.l = c1.c.coef*not(E.c.l);
    c1.c.t = c1.c.coef*not(E.c.t);
    c1.c.r = c1.c.coef*not(E.c.r);
    c1.c.b = c1.c.coef*not(E.c.b);

    Ti.x.l = dT.x;
    Ti.x.t = dT.x;
    Ti.x.r = dT.x;
    Ti.x.b = dT.x;
    Ti.f.l = dT.f;
    Ti.f.t = dT.f;
    Ti.f.r = dT.f;
    Ti.f.b = dT.f;
    Ti.n.l = dT.n;
    Ti.n.t = dT.n;
    Ti.n.r = dT.n;
    Ti.n.b = dT.n;
    Ti.c.l = dT.c;
    Ti.c.t = dT.c;
    Ti.c.r = dT.c;
    Ti.c.b = dT.c;

    Bi.x.l = dT.x;
    Bi.x.t = dT.x;
    Bi.x.r = dT.x;
    Bi.x.b = dT.x;
    Bi.f.l = dT.f;
    Bi.f.t = dT.f;
    Bi.f.r = dT.f;
    Bi.f.b = dT.f;
    Bi.n.l = dT.n;
    Bi.n.t = dT.n;
    Bi.n.r = dT.n;
    Bi.n.b = dT.n;
    Bi.c.l = dT.c;
    Bi.c.t = dT.c;
    Bi.c.r = dT.c;
    Bi.c.b = dT.c;

    % inner boundaries
    Tx2f = zeros(nx, 1); % xfine to fine
    Tf2x = zeros(nf, 1); % fine to xfine
    Tf2n = zeros(nf, 1); % fine to normal
    Tn2f = zeros(nn, 1); % etc.
    Tn2c = zeros(nn, 1);
    Tc2n = zeros(nc, 1);

    while i <= nt
        tic;

        % moving beam centre relative to edge
        %rho_p = ;

        % node temps
        Ti.x.n = dT.x + T0.x;
        Ti.f.n = dT.f + T0.f;
        Ti.n.n = dT.n + T0.n;
        Ti.c.n = dT.c + T0.c;
        Ti.x.n4 = Ti.x.n.*Ti.x.n.*Ti.x.n.*Ti.x.n;
        Ti.f.n4 = Ti.f.n.*Ti.f.n.*Ti.f.n.*Ti.f.n;
        Ti.n.n4 = Ti.n.n.*Ti.n.n.*Ti.n.n.*Ti.n.n;
        Ti.c.n4 = Ti.c.n.*Ti.c.n.*Ti.c.n.*Ti.c.n;

        % XFINE
        % left, top, right, bottom temps
        Ti.x.l(:, 2:end) = Ti.x.n(:, 1:end-1);
        Ti.x.t(2:end, :) = Ti.x.n(1:end-1, :);
        Ti.x.r(:, 1:end-1) = Ti.x.n(:, 2:end);
        Ti.x.b(1:end-1, :) = Ti.x.n(2:end, :);
        % boundary conditions
        %Bi.x.l = ;
        %Bi.x.t = ;
        Tx2f(1:2:end-1) = Ti.f.n(:, 1);
        Tx2f(2:2:end) = Ti.f.n(:, 1);
        %Tx2f = kron(Ti.f.n(:, 1), [1; 1]);
        Bi.x.r(:, end) = c1.x.coef*(Tx2f-Ti.x.n(:, end));
        %Bi.x.b = ;
        % update temp
        T.x = Ti.x.n + ... 
              c1.x.l.*(Ti.x.l-Ti.x.n) + Bi.x.l + ...
              c1.x.t.*(Ti.x.t-Ti.x.n) + Bi.x.t + ...
              c1.x.r.*(Ti.x.r-Ti.x.n) + Bi.x.r + ...
              c1.x.b.*(Ti.x.b-Ti.x.n) + Bi.x.b + ...
              c3.x + ...
              -c4*(Ti.x.n4-T04.x);
        dT.x = T.x - T0.x;

        % FINE
        % left, top, right, bottom temps
        Ti.f.l(:, 2:end) = Ti.f.n(:, 1:end-1);
        Ti.f.t(2:end, :) = Ti.f.n(1:end-1, :);
        Ti.f.r(:, 1:end-1) = Ti.f.n(:, 2:end);
        Ti.f.b(1:end-1, :) = Ti.f.n(2:end, :);
        % boundary conditions
        Tf2x = (Ti.x.n(1:2:end-1, end) + Ti.x.n(2:2:end, end))/2;
        Bi.f.l(:, 1) = c1.f.coef*(Tf2x-Ti.f.n(:, 1));
        %Bi.f.t = ;
        Tf2n(1:2:end-1) = Ti.n.n(:, 1);
        Tf2n(2:2:end) = Ti.n.n(:, 1);
        %Tf2n = kron(Ti.n.n(:, 1), [1; 1]);
        Bi.f.r(:, end) = c1.f.coef*(Tf2n-Ti.f.n(:, end));
        %Bi.f.b = ;
        % update temp
        T.f = Ti.f.n + ... 
              c1.f.l.*(Ti.f.l-Ti.f.n) + Bi.f.l + ...
              c1.f.t.*(Ti.f.t-Ti.f.n) + Bi.f.t + ...
              c1.f.r.*(Ti.f.r-Ti.f.n) + Bi.f.r + ...
              c1.f.b.*(Ti.f.b-Ti.f.n) + Bi.f.b + ...
              c3.f + ...
              -c4*(Ti.f.n4-T04.f);
        dT.f = T.f - T0.f;


        % NORMAL
        Ti.n.l(:, 2:end) = Ti.n.n(:, 1:end-1);
        Ti.n.t(2:end, :) = Ti.n.n(1:end-1, :);
        Ti.n.r(:, 1:end-1) = Ti.n.n(:, 2:end);
        Ti.n.b(1:end-1, :) = Ti.n.n(2:end, :);
        % boundary conditions
        Tn2f = (Ti.f.n(1:2:end-1, end) + Ti.f.n(2:2:end, end))/2;
        Bi.n.l(:, 1) = c1.n.coef*(Tn2f-Ti.n.n(:, 1));
        %Bi.n.t = ;
        Tn2c(1:2:end-1) = Ti.c.n(:, 1);
        Tn2c(2:2:end) = Ti.c.n(:, 1);
        %Tn2c = kron(Ti.c.n(:, 1), [1; 1]);
        Bi.n.r(:, end) = c1.n.coef*(Tn2c-Ti.n.n(:, end));
        %Bi.n.b = ;
        % update temp
        T.n = Ti.n.n + ... 
              c1.n.l.*(Ti.n.l-Ti.n.n) + Bi.n.l + ...
              c1.n.t.*(Ti.n.t-Ti.n.n) + Bi.n.t + ...
              c1.n.r.*(Ti.n.r-Ti.n.n) + Bi.n.r + ...
              c1.n.b.*(Ti.n.b-Ti.n.n) + Bi.n.b + ...
              c3.n + ...
              -c4*(Ti.n.n4-T04.n);
        dT.n = T.n - T0.n;

        % COARSE
        Ti.c.l(:, 2:end) = Ti.c.n(:, 1:end-1);
        Ti.c.t(2:end, :) = Ti.c.n(1:end-1, :);
        Ti.c.r(:, 1:end-1) = Ti.c.n(:, 2:end);
        Ti.c.b(1:end-1, :) = Ti.c.n(2:end, :);
        % boundary conditions
        Tc2n = (Ti.n.n(1:2:end-1, end) + Ti.n.n(2:2:end, end))/2;
        Bi.c.l(:, 1) = c1.c.coef*(Tc2n-Ti.c.n(:, 1));
        %Bi.c.t = ;
        %Bi.c.r = ;
        %Bi.c.b = ;
        % update temp
        T.c = Ti.c.n + ... 
              c1.c.l.*(Ti.c.l-Ti.c.n) + Bi.c.l + ...
              c1.c.t.*(Ti.c.t-Ti.c.n) + Bi.c.t + ...
              c1.c.r.*(Ti.c.r-Ti.c.n) + Bi.c.r + ...
              c1.c.b.*(Ti.c.b-Ti.c.n) + Bi.c.b + ...
              c3.c + ...
              -c4*(Ti.c.n4-T04.c);
        dT.c = T.c - T0.c;

        if rem(i, nframes) == 0
            ind = i/nframes;
            dTframes(:, 1:wx, ind) = dT.x;
            dTframes(:, wx+1:wx*2, ind) = kron(dT.f, ones(2));
            dTframes(:, (wx*2)+1:wx*4, ind) = kron(dT.n, ones(4));
            dTframes(:, (wx*4)+1:end, ind) = kron(dT.c, ones(8));
            fprintf('%.4f%% complete. Animation frame captured.\n', 100*i/nt)
        end

        i = i + 1;

        toc;
    end

end

