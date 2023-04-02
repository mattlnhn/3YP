function [dT1, dT2, dT3, dT4, dT5, dT6, dT7] = heating2(dlf, Lf, dl, L, dlc, Lc, theta, dt, nt, fps, mat, beam)
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
%   beampos immediately assigned to beam.pos, workaround for parallel
%           processing
%   OUTPUT ----------------------------------------------------------------
%   dTframes                temp matrix at 1000fps
%   NOTES -----------------------------------------------------------------
%   bethe function uses combination of CGS and eV for units; heat transfer
%   working in SI units hence unit conversion AFTER bethe function has been
%   called

    %% PHYSICAL CONSTANTS
    const_sigma = 5.670374419e-8; % stefan-boltzmann, W m-2 K-4
    T0 = 1.9;
    T04 = 1.9^4;

    %% MESH
    % N.B. system origin at midpoint of leading (left) edge
    %{
    5 5 5 5 5 7 7 7 7 7 7 7 7 7 7
    5 5 5 5 5 7 7 7 7 7 7 7 7 7 7
    5 5 5 5 5 7 7 7 7 7 7 7 7 7 7
    5 5 5 5 5 7 7 7 7 7 7 7 7 7 7
    5 5 5 5 5 7 7 7 7 7 7 7 7 7 7
    2 4 4 4 4 7 7 7 7 7 7 7 7 7 7
    2 4 4 4 4 7 7 7 7 7 7 7 7 7 7
    1 4 4 4 4 7 7 7 7 7 7 7 7 7 7
    3 4 4 4 4 7 7 7 7 7 7 7 7 7 7
    3 4 4 4 4 7 7 7 7 7 7 7 7 7 7
    6 6 6 6 6 7 7 7 7 7 7 7 7 7 7
    6 6 6 6 6 7 7 7 7 7 7 7 7 7 7
    6 6 6 6 6 7 7 7 7 7 7 7 7 7 7
    6 6 6 6 6 7 7 7 7 7 7 7 7 7 7
    6 6 6 6 6 7 7 7 7 7 7 7 7 7 7
    %}

    mf = round(dl/dlf)
    mc = round(dlc/dl)

    nf = round(Lf/dlf); % fine nodes
    n = round(L/dl); % normal nodes
    nc = round(Lc/dlc); % coarse nodes

    % section sizes in normal nodes
    n_h2 = round(n*(L-Lf)/(2*L));       % h of (2 || 3)
    n_h12 = round(n*(L+Lf)/(2*L));      % h of 1 && (2 || 3)
    n_1 = round(n*Lf/L);                % h & w of 1 = w of (2 || 3)
    n_w2 = n_1;                         % alias for above
    n_w4 = round(n*(L-Lf)/L);           % w of 4

    % section sizes in coarse nodes
    c_h5 = round(nc*(Lc-L)/(2*Lc));     % h of (5 || 6)
    c_h45 = round(nc*(Lc+L)/(2*Lc));    % h of 4 && (5 || 6)
    c_1 = round(nc*Lf/Lc);              % h & w of 1 = w of (2 || 3)
    c_w2 = c_1;                         % alias for above
    c_h4 = round(nc*L/Lc);              % h of 4 = w of (2 || 3) && 4
    c_w24 = c_h4;                       % alias for above
    c_w4 = round(nc*(L-Lf)/Lc);         % w of 4
    c_w7 = round(nc*(Lc-L)/Lc);         % w of 7

    % section 1--fine
    dT1 = zeros(nf, nf);
    E1L = dT1; E1L(:, 1) = 1;
    E1T = dT1; E1T(1, :) = 1;
    E1R = dT1; E1R(:, end) = 1;
    E1B = dT1; E1B(end, :) = 1;

    % section 2--normal
    dT2 = zeros(n_h2, n_1);
    E2L = dT2; E2L(:, 1) = 1;
    E2T = dT2; E2T(1, :) = 1;
    E2R = dT2; E2R(:, end) = 1;
    E2B = dT2; E2B(end, :) = 1;

    % section 3--normal
    dT3 = zeros(n_h2, n_1);
    E3L = dT3; E3L(:, 1) = 1;
    E3T = dT3; E3T(1, :) = 1;
    E3R = dT3; E3R(:, end) = 1;
    E3B = dT3; E3B(end, :) = 1;

    % section 4--normal
    dT4 = zeros(n, n_w4);
    E4L = dT4; E4L(:, 1) = 1;
    E4T = dT4; E4T(1, :) = 1;
    E4R = dT4; E4R(:, end) = 1;
    E4B = dT4; E4B(end, :) = 1;

    % section 5--coarse
    dT5 = zeros(c_h5, c_h4);
    E5L = dT5; E5L(:, 1) = 1;
    E5T = dT5; E5T(1, :) = 1;
    E5R = dT5; E5R(:, end) = 1;
    E5B = dT5; E5B(end, :) = 1;

    % section 6--coarse
    dT6 = zeros(c_h5, c_h4);
    E6L = dT6; E6L(:, 1) = 1;
    E6T = dT6; E6T(1, :) = 1;
    E6R = dT6; E6R(:, end) = 1;
    E6B = dT6; E6B(end, :) = 1;

    % section 7--coarse
    dT7 = zeros(nc, c_w7);
    E7L = dT7; E7L(:, 1) = 1;
    E7T = dT7; E7T(1, :) = 1;
    E7R = dT7; E7R(:, end) = 1;
    E7B = dT7; E7B(end, :) = 1;

    %% ENERGY DEPOSITION
    beta = (1-beam.gamma^-2)^.5;
    % mass stopping power = <1/rho dE/dx>
    msp = bethe(beam.M, mat.A, mat.Z, mat.I, mat.w, mat.rho, beta, 1);
    rho_p1 = densdist(beam.pos, dlf, Lf, Lf, 0, -Lf/2, beam.np, ...
             beam.nb, beam.f, beam.sigma, beam.sigma);
    rho_p2 = densdist(beam.pos, dl, Lf, (L-Lf)/2, 0, +Lf/2, beam.np, ...
             beam.nb, beam.f, beam.sigma, beam.sigma);
    rho_p3 = densdist(beam.pos, dl, Lf, (L-Lf)/2, 0, -L/2, beam.np, ...
             beam.nb, beam.f, beam.sigma, beam.sigma);
    rho_p4 = densdist(beam.pos, dl, L-Lf, L, +Lf, -L/2, beam.np, ...
             beam.nb, beam.f, beam.sigma, beam.sigma);
    % rho_p zero for sections 5-7

    %% CONVERT UNITS
    % mass stopping power MeV g-1 cm2 --> J g-1 m2
    % * joules per mev * (m per cm)^2
    msp = msp * 1.6021773e-13 * (0.1)^2;
    % density g cm-3 --> kg m-3
    mat.rho = mat.rho * 1000;

    %% ANIMATION OUTPUT PARAMS
    %{
    ntframes = round((dt*fps)^-1); % no. of time steps btwn frames
    frames = round(nt/ntframes); % no. of frames total
    dTframes = zeros(nn/2, nn, frames); % preallocate
    %}

    %% COEFFICIENTS
    % inner node conduction coefficients
    condf = (mat.k*dt)/(mat.rho*mat.c_p*dlf*dlf); % conduction, fine
    cond = (mat.k*dt)/(mat.rho*mat.c_p*dl*dl); % conduction, normal
    condc = (mat.k*dt)/(mat.rho*mat.c_p*dlc*dlc); % conduction, coarse

    % proton energy deposition coefficients
    pro1 = (dt*msp*rho_p1)/mat.c_p; % section 1
    pro2 = (dt*msp*rho_p2)/mat.c_p; % section 2
    pro3 = (dt*msp*rho_p3)/mat.c_p; % section 3
    pro4 = (dt*msp*rho_p4)/mat.c_p; % section 4
    % rho_p zero for sections 5-7

    % radiation
    % all cells front & back
    rad = (2*mat.epsilon*const_sigma*dt)/(mat.rho*mat.c_p*theta);
    % fine edge
    radEf = (mat.epsilon*const_sigma*dt)/(mat.rho*mat.c_p*dlf);
    % normal edge
    radE = (mat.epsilon*const_sigma*dt)/(mat.rho*mat.c_p*dl);
    % coarse edge
    radEc = (mat.epsilon*const_sigma*dt)/(mat.rho*mat.c_p*dlc);

    % exclude edges of each section from inner node conduction
    cond1L = condf*not(E1L); cond1T = condf*not(E1T);
    cond1R = condf*not(E1R); cond1B = condf*not(E1B);
    cond2L = cond*not(E2L); cond2T = cond*not(E2T);
    cond2R = cond*not(E2R); cond2B = cond*not(E2B);
    cond3L = cond*not(E3L); cond3T = cond*not(E3T);
    cond3R = cond*not(E3R); cond3B = cond*not(E3B);
    cond4L = cond*not(E4L); cond4T = cond*not(E4T);
    cond4R = cond*not(E4R); cond4B = cond*not(E4B);
    cond5L = condc*not(E5L); cond5T = condc*not(E5T);
    cond5R = condc*not(E5R); cond5B = condc*not(E5B);
    cond6L = condc*not(E6L); cond6T = condc*not(E6T);
    cond6R = condc*not(E6R); cond6B = condc*not(E6B);
    cond7L = condc*not(E7L); cond7T = condc*not(E7T);
    cond7R = condc*not(E7R); cond7B = condc*not(E7B);

    % preallocate matrices for boundary conditions
    B1L = dT1; B1T = dT1; B1R = dT1; B1B = dT1;
    B2L = dT2; B2T = dT2; B2R = dT2; B2B = dT2;
    B3L = dT3; B3T = dT3; B3R = dT3; B3B = dT3;
    B4L = dT4; B4T = dT4; B4R = dT4; B4B = dT4;
    B5L = dT5; B5T = dT5; B5R = dT5; B5B = dT5;
    B6L = dT6; B6T = dT6; B6R = dT6; B6B = dT6;
    B7L = dT7; B7T = dT7; B7R = dT7; B7B = dT7;

    % preallocate inner boundaries
    % conduction into sec 1
    T2to1 = zeros(1, mf*n_1); % n->f
    T3to1 = zeros(1, mf*n_1); % n->f
    T4to1 = zeros(mf*n_1, 1); % n->f
    % conduction into sec 2
    T1to2 = zeros(1, n_1); % f->n
    T4to2 = zeros(n_h2, 1); % n
    T5to2 = zeros(1, n_w2); % c->n
    % conduction into sec 3
    T1to3 = zeros(1, n_1); % f->n
    T4to3 = zeros(n_h2, 1); % n
    T6to3 = zeros(1, n_w2); % c->n
    % conduction into sec 4
    T1to4 = zeros(n_1, 1); % f->n
    T2to4 = zeros(n_h2, 1); % n
    T3to4 = zeros(n_h2, 1); % n
    T5to4 = zeros(1, n_w4); % c->n
    T6to4 = zeros(1, n_w4); % c->n
    T7to4 = zeros(n, 1); % c->n
    % conduction into sec 5
    T2to5 = zeros(1, c_w2); % n->c
    T4to5 = zeros(1, c_w4); % n->c
    T7to5 = zeros(c_h5, 1); % c
    % conduction into sec 6
    T3to6 = zeros(1, c_w2); % n->c
    T4to6 = zeros(1, c_w4); % n->c
    T7to6 = zeros(c_h5, 1); % c
    % conduction into sec 7
    T4to7 = zeros(c_h4, 1); % n->c
    T5to7 = zeros(c_h5, 1); % c
    T6to7 = zeros(c_h5, 1); % c

    %% TIME ITERATION
    iter = 1;
    %starttime = tic;
    
    while iter <= nt

        % node temps
        Ti1 = dT1 + T0;
        Ti2 = dT2 + T0;
        Ti3 = dT3 + T0;
        Ti4 = dT4 + T0;
        Ti5 = dT5 + T0;
        Ti6 = dT6 + T0;
        Ti7 = dT7 + T0;
        Ti14 = Ti1.*Ti1.*Ti1.*Ti1;
        Ti24 = Ti2.*Ti2.*Ti2.*Ti2;
        Ti34 = Ti3.*Ti3.*Ti3.*Ti3;
        Ti44 = Ti4.*Ti4.*Ti4.*Ti4;
        Ti54 = Ti5.*Ti5.*Ti5.*Ti5;
        Ti64 = Ti6.*Ti6.*Ti6.*Ti6;
        Ti74 = Ti7.*Ti7.*Ti7.*Ti7;

        % SECTION 1 #######################################################
        % boundary conditions
        % left edge radiation
        B1L(:, 1) = -radEf*(Ti14(:, 1)-T04);
        % top edge conduction 2->1
        for i = 1:n_1 % larger cells on 2->1 boundary
            T2to1(1, mf*(i-1)+1:mf*i) = Ti2(end, i);
        end
        B1T(1, :) = condf*(T2to1-Ti1(1, :));
        % right edge conduction 4->1
        for i = 1:n_1 % larger cells on 4->1 boundary
            T4to1(mf*(i-1)+1:mf*i, 1) = Ti4(n_h2+i, 1);
        end
        B1R(:, end) = condf*(T4to1-Ti1(:, end));
        % bottom edge conduction 3->1
        for i = 1:n_1 % larger cells on 3->1 boundary
            T3to1(1, mf*(i-1)+1:mf*i) = Ti3(1, i);
        end
        B1B(end, :) = condf*(T3to1-Ti1(end, :));

        % update temp
        Tu1 = Ti1 + ...
              cond1L.*(circshift(Ti1, [0 1])-Ti1) + B1L + ...
              cond1T.*(circshift(Ti1, [1 0])-Ti1) + B1T + ...
              cond1R.*(circshift(Ti1, [0 -1])-Ti1) + B1R + ...
              cond1B.*(circshift(Ti1, [-1 0])-Ti1) + B1B + ...
              pro1 + ...
              -rad*(Ti14-T04);
        dT1 = Tu1 - T0;

        % SECTION 2 #######################################################
        % boundary conditions
        % left edge radiation
        B2L(:, 1) = -radE*(Ti24(:, 1)-T04);
        % top edge conduction 5->2
        for i = 1:c_w2 % larger cells on 5->2 boundary
            T5to2(1, mc*(i-1)+1:mc*i) = Ti5(end, i);
        end
        B2T(1, :) = cond*(T5to2-Ti2(1, :));
        % right edge conduction sec 4->2
        B2R(:, end) = cond*(Ti4(1:n_h2, 1)-Ti2(:, end));
        % bottom edge conduction sec 1->2
        for i = 1:n_1 % no. of large cells on 1->2 boundary
            rc = mf*(i-1)+1:mf*i;
            T1to2(1, i) = sum(Ti1(1:mf, rc), 'all')/mf^2;
        end
        B2B(end, :) = cond*(T1to2-Ti2(end, :));

        % update temp
        Tu2 = Ti2 + ...
              cond2L.*(circshift(Ti2, [0 1])-Ti2) + B2L + ...
              cond2T.*(circshift(Ti2, [1 0])-Ti2) + B2T + ...
              cond2R.*(circshift(Ti2, [0 -1])-Ti2) + B2R + ...
              cond2B.*(circshift(Ti2, [-1 0])-Ti2) + B2B + ...
              pro2 + ...
              -rad*(Ti24-T04);
        dT2 = Tu2 - T0;

        % SECTION 3 #######################################################
        % boundary conditions
        % left edge radiation
        B3L(:, 1) = -radE*(Ti34(:, 1)-T04);
        % top edge conduction sec 1->3
        for i = 1:n_1 % larger cells on 1->3 boundary
            rc = mf*(i-1)+1:mf*i;
            T1to3(1, i) = sum(Ti1((end-mf)+1:end, rc), 'all')/mf^2;
        end
        B3T(1, :) = cond*(T1to3-Ti3(1, :));
        % right edge conduction sec 4->3
        B3R(:, end) = cond*(Ti4(n_h12+1:end, 1)-Ti3(:, end));
        % bottom edge conduction sec 6->3
        for i = 1:c_w2 % larger cells on 6->3 boundary
            T6to3(1, mc*(i-1)+1:mc*i) = Ti6(1, i);
        end
        B3B(end, :) = cond*(T6to3-Ti3(end, :));
        
        % update temp
        Tu3 = Ti3 + ...
              cond3L.*(circshift(Ti3, [0 1])-Ti3) + B3L + ...
              cond3T.*(circshift(Ti3, [1 0])-Ti3) + B3T + ...
              cond3R.*(circshift(Ti3, [0 -1])-Ti3) + B3R + ...
              cond3B.*(circshift(Ti3, [-1 0])-Ti3) + B3B + ...
              pro3 + ...
              -rad*(Ti34-T04);
        dT3 = Tu3 - T0;

        % SECTION 4 #######################################################
        % boundary conditions
        % left edge conduction 2->4
        B4L(1:n_h2, 1) = cond*(Ti2(:, end)-Ti4(1:n_h2, 1));
        % left edge conduction 1->4
        for i = 1:n_1 % no. of large cells on 1->4 boundary
            rr = mf*(i-1)+1:mf*i;
            T1to4(i, 1) = sum(Ti1(rr, (end-mf)+1:end), 'all')/mf^2;
        end
        B4L(n_h2+1:n_h12, 1) = cond*(T1to4-Ti4(n_h2+1:n_h12, 1));
        % left edge conduction 3->4
        B4L(n_h12+1:end, 1) = cond*(Ti3(:, end)-Ti4(n_h12+1:end, 1));
        % top edge conduction 5->4
        for i = 1:c_w4 % larger nodes on 5->4 boundary
            T5to4(1, mc*(i-1)+1:mc*i) = Ti5(end, c_w2+i);
        end
        B4T(1, :) = cond*(T5to4-Ti4(1, :));
        % right edge conduction 7->4
        for i = 1:c_h4 % larger nodes on 7->4 boundary
            T7to4(mc*(i-1)+1:mc*i, 1) = Ti7(c_h5+1:c_h45, 1);
        end
        B4R(:, end) = cond*(T7to4-Ti4(:, end));
        % bottom edge conduction 6->4
        for i = 1:c_w4 % larger nodes on 6->4 boundary
            T6to4(1, mc*(i-1)+1:mc*i) = Ti6(1, c_w2+i);
        end
        B4B(end, :) = cond*(T6to4-Ti4(end, :));

        % update temp
        Tu4 = Ti4 + ...
              cond4L.*(circshift(Ti4, [0 1])-Ti4) + B4L + ...
              cond4T.*(circshift(Ti4, [1 0])-Ti4) + B4T + ...
              cond4R.*(circshift(Ti4, [0 -1])-Ti4) + B4R + ...
              cond4B.*(circshift(Ti4, [-1 0])-Ti4) + B4B + ...
              pro4 + ...
              -rad*(Ti44-T04);
        dT4 = Tu4 - T0;

        % SECTION 5 #######################################################
        % boundary conditions
        % left edge radiation
        B5L(:, 1) = -radEc*(Ti54(:, 1)-T04);
        % top edge radiation
        B5T(1, :) = -radEc*(Ti54(1, :)-T04);
        % right edge conduction 7->5
        T7to5 = Ti7(1:c_h5, 1);
        B5R(:, end) = condc*(T7to5-Ti5(:, end));
        % bottom edge conduction 2->5
        for i = 1:c_w2 % larger nodes on 2->5 boundary
            rc = mc*(i-1)+1:mc*i;
            T2to5(1, i) = sum(Ti2(1:mc, rc), 'all')/mc^2;
        end
        B5B(end, 1:c_w2) = condc*(T2to5-Ti5(end, 1:c_w2));
        % bottom edge conduction 4->5
        for i = 1:c_w4 % larger nodes on 4->5 boundary
            rc = mc*(i-1)+1:mc*i;
            T4to5(1, i) = sum(Ti4(1:mc, rc), 'all')/mc^2;
        end
        B5B(end, c_w2+1:end) = condc*(T4to5-Ti5(end, c_w2+1:end));

        % update temp
        Tu5 = Ti5 + ...
              cond5L.*(circshift(Ti5, [0 1])-Ti5) + B5L + ...
              cond5T.*(circshift(Ti5, [1 0])-Ti5) + B5T + ...
              cond5R.*(circshift(Ti5, [0 -1])-Ti5) + B5R + ...
              cond5B.*(circshift(Ti5, [-1 0])-Ti5) + B5B + ...
              -rad*(Ti54-T04);
        dT5 = Tu5 - T0;

        % SECTION 6 #######################################################
        % boundary conditions
        % left edge radiation
        B6L(:, 1) = -radEc*(Ti64(:, 1)-T04);
        % top edge conduction 3->6
        for i = 1:c_w2 % larger nodes on 3->6 boundary
            rc = mc*(i-1)+1:mc*i;
            T3to6(1, i) = sum(Ti3((end-mc)+1:end, rc), 'all')/mc^2;
        end
        B6T(1, 1:c_w2) = condc*(T3to6-Ti6(1, 1:c_w2));
        % top edge conduction 4->6
        for i = 1:c_w4
            rc = mc*(i-1)+1:mc*i;
            T4to6(1, i) = sum(Ti4((end-mc)+1:end, rc), 'all')/mc^2;
        end
        B6T(1, c_w2+1:end) = condc*(T4to6-Ti6(1, c_w2+1:end));
        % right edge conduction 7->6
        T7to6 = Ti7(c_h45+1:end, 1);
        B6R(:, end) = condc*(T7to6-Ti6(:, end));
        % bottom edge radiation
        B6B(end, :) = -radEc*(Ti64(end, :)-T04);

        % update temp
        Tu6 = Ti6 + ...
              cond6L.*(circshift(Ti6, [0 1])-Ti6) + B6L + ...
              cond6T.*(circshift(Ti6, [1 0])-Ti6) + B6T + ...
              cond6R.*(circshift(Ti6, [0 -1])-Ti6) + B6R + ...
              cond6B.*(circshift(Ti6, [-1 0])-Ti6) + B6B + ...
              -rad*(Ti64-T04);
        dT6 = Tu6 - T0;

        % SECTION 7 #######################################################
        % boundary conditions
        % left edge conduction 5->7
        T5to7 = Ti5(:, end);
        B7L(1:c_h5, 1) = condc*(T5to7-Ti7(1:c_h5, 1));
        % left edge conduction 4->7
        for i = 1:c_h4 % larger nodes on 4->7 boundary
            rr = mc*(i-1)+1:mc*i;
            T4to7(i, 1) = sum(Ti4(rr, (end-mc)+1:end), 'all')/mc^2;
        end
        B7L(c_h5+1:c_h45, 1) = condc*(T4to7-Ti7(c_h5+1:c_h45, 1));
        % left edge conduction 6->7
        T6to7 = Ti6(:, end);
        B7L(c_h45+1:end, 1) = condc*(T6to7-Ti7(c_h45+1:end, 1));
        % top edge radiation
        B7T(1, :) = -radEc*(Ti74(1, :)-T04);
        % right edge radiation
        B7R(:, end) = -radEc*(Ti74(:, end)-T04);
        % bottom edge radiation
        B7B(end, :) = -radEc*(Ti74(end, :)-T04);

        % update temp
        Tu7 = Ti7 + ...
              cond7L.*(circshift(Ti7, [0 1])-Ti7) + B7L + ...
              cond7T.*(circshift(Ti7, [1 0])-Ti7) + B7T + ...
              cond7R.*(circshift(Ti7, [0 -1])-Ti7) + B7R + ...
              cond7B.*(circshift(Ti7, [-1 0])-Ti7) + B7B + ...
              -rad*(Ti74-T04);
        dT7 = Tu7 - T0;




        %{
        % SAVE ############################################################
        if rem(i, nframes) == 0
            ind = i/nframes;
            dTframes(:, 1:nn/100, ind) = dTx(1:4:end, 1:4:end);
            dTframes(:, (nn/100)+1:2*nn/100, ind) = dTf(1:2:end, 1:2:end);
            dTframes(:, (2*nn/100)+1:3*nn/100, ind) = dTn;
            dTframes(:, (3*nn/100)+1:end, ind) = kron(dTc, ones(2));
            endtime = toc(starttime);
            fprintf('%.4f%% complete in %.4f s. Frame captured.\n', 100*i/nt, endtime)
            starttime = tic;
        end
        %}
        fprintf('%.4f%% complete\n', 100*iter/nt)
        iter = iter + 1;
        

    end

end

