function [dT1out, dT2out, dT4out, dT5out, dT7out] = heating2A(dlf, Lf, dl, L, dlc, Lc, theta, dt, nt, frames, mat, beam, phi)
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
    5 5 5 5 5 5 7 7 7 7 7 7 7 7 7 7
    5 5 5 5 5 5 7 7 7 7 7 7 7 7 7 7
    5 5 5 5 5 5 7 7 7 7 7 7 7 7 7 7
    5 5 5 5 5 5 7 7 7 7 7 7 7 7 7 7
    5 5 5 5 5 5 7 7 7 7 7 7 7 7 7 7
    2 2 4 4 4 4 7 7 7 7 7 7 7 7 7 7
    2 2 4 4 4 4 7 7 7 7 7 7 7 7 7 7
    1 1 4 4 4 4 7 7 7 7 7 7 7 7 7 7
------------------------------------- Q = 0 (horizontal symmetry)
    1 1 4 4 4 4 7 7 7 7 7 7 7 7 7 7
    3 3 4 4 4 4 7 7 7 7 7 7 7 7 7 7
    3 3 4 4 4 4 7 7 7 7 7 7 7 7 7 7
    6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7
    6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7
    6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7
    6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7
    6 6 6 6 6 6 7 7 7 7 7 7 7 7 7 7
    %}

    mf = round(dl/dlf);
    mc = round(dlc/dl);

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
    dT1 = zeros(nf/2, nf);
    dT1out = zeros(nf/2, nf, frames);
    E1L = dT1; E1L(:, 1) = 1;
    E1T = dT1; E1T(1, :) = 1;
    E1R = dT1; E1R(:, end) = 1;
    E1B = dT1; E1B(end, :) = 1;

    % section 2--normal
    dT2 = zeros(n_h2, n_1);
    dT2out = zeros(n_h2, n_1, frames);
    E2L = dT2; E2L(:, 1) = 1;
    E2T = dT2; E2T(1, :) = 1;
    E2R = dT2; E2R(:, end) = 1;
    E2B = dT2; E2B(end, :) = 1;

    % section 4--normal
    dT4 = zeros(n/2, n_w4);
    dT4out = zeros(n/2, n_w4, frames);
    E4L = dT4; E4L(:, 1) = 1;
    E4T = dT4; E4T(1, :) = 1;
    E4R = dT4; E4R(:, end) = 1;
    E4B = dT4; E4B(end, :) = 1;

    % section 5--coarse
    dT5 = zeros(c_h5, c_h4);
    dT5out = zeros(c_h5, c_h4, frames);
    E5L = dT5; E5L(:, 1) = 1;
    E5T = dT5; E5T(1, :) = 1;
    E5R = dT5; E5R(:, end) = 1;
    E5B = dT5; E5B(end, :) = 1;

    % section 7--coarse
    dT7 = zeros(nc/2, c_w7);
    dT7out = zeros(nc/2, c_w7, frames);
    E7L = dT7; E7L(:, 1) = 1;
    E7T = dT7; E7T(1, :) = 1;
    E7R = dT7; E7R(:, end) = 1;
    E7B = dT7; E7B(end, :) = 1;

    %% ENERGY DEPOSITION
    beta = (1-beam.gamma^-2)^.5;
    % mass stopping power = <1/rho dE/dx>
    msp = bethe(beam.M, mat.A, mat.Z, mat.I, mat.w, mat.rho, beta, 1);
    rho_p1 = densdist(beam.pos, dlf, Lf, Lf/2, 0, 0, phi, beam.np, ...
             beam.nb, beam.f, beam.sigma, beam.sigma);
    rho_p2 = densdist(beam.pos, dl, Lf, (L-Lf)/2, 0, +Lf/2, phi, beam.np, ...
             beam.nb, beam.f, beam.sigma, beam.sigma);
    rho_p4 = densdist(beam.pos, dl, L-Lf, L/2, +Lf, 0, beam.np, phi, ...
             beam.nb, beam.f, beam.sigma, beam.sigma);
    % rho_p zero for sections 5-7

    %% CONVERT UNITS
    % mass stopping power MeV g-1 cm2 --> J g-1 m2
    % * joules per mev * (m per cm)^2
    msp = msp * 1.6021773e-13 * (0.1)^2;
    % density g cm-3 --> kg m-3
    mat.rho = mat.rho * 1000;

    %% ANIMATION OUTPUT PARAMS
    f2f = floor(nt/frames); % steps from frame to frame

    %% COEFFICIENTS
    % inner node conduction coefficients
    condfV = (mat.k*dt)/(mat.rho*mat.c_p*dlf*dlf); % conduction, fine
    condfH = (mat.k*dt)/(mat.rho*mat.c_p*dlf*cosd(phi)*dlf*cosd(phi)); % conduction, fine
    condV = (mat.k*dt)/(mat.rho*mat.c_p*dl*dl); % conduction, normal
    condH = (mat.k*dt)/(mat.rho*mat.c_p*dl*cosd(phi)*dl*cosd(phi)); % conduction, normal
    condcV = (mat.k*dt)/(mat.rho*mat.c_p*dlc*dlc); % conduction, coarse
    condcH = (mat.k*dt)/(mat.rho*mat.c_p*dlc*cosd(phi)*dlc*cosd(phi)); % conduction, coarse

    % proton energy deposition coefficients
    pro1 = (dt*msp*rho_p1)/(mat.c_p); % section 1
    pro2 = (dt*msp*rho_p2)/(mat.c_p); % section 2
    pro4 = (dt*msp*rho_p4)/(mat.c_p); % section 4
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
    cond1L = condfH*not(E1L); cond1T = condfV*not(E1T);
    cond1R = condfH*not(E1R); cond1B = condfV*not(E1B);
    cond2L = condH*not(E2L); cond2T = condV*not(E2T);
    cond2R = condH*not(E2R); cond2B = condV*not(E2B);
    cond4L = condH*not(E4L); cond4T = condV*not(E4T);
    cond4R = condH*not(E4R); cond4B = condV*not(E4B);
    cond5L = condcH*not(E5L); cond5T = condcV*not(E5T);
    cond5R = condcH*not(E5R); cond5B = condcV*not(E5B);
    cond7L = condcH*not(E7L); cond7T = condcV*not(E7T);
    cond7R = condcH*not(E7R); cond7B = condcV*not(E7B);

    % preallocate matrices for boundary conditions
    B1L = dT1; B1T = dT1; B1R = dT1; B1B = dT1;
    B2L = dT2; B2T = dT2; B2R = dT2; B2B = dT2;
    B4L = dT4; B4T = dT4; B4R = dT4; B4B = dT4;
    B5L = dT5; B5T = dT5; B5R = dT5; B5B = dT5;
    B7L = dT7; B7T = dT7; B7R = dT7; B7B = dT7;

    % preallocate inner boundaries
    % conduction into sec 1
    T2to1 = zeros(1, mf*n_1); % n->f
    T4to1 = zeros(mf*n_1/2, 1); % n->f
    % conduction into sec 2
    T1to2 = zeros(1, n_1); % f->n
    T4to2 = zeros(n_h2, 1); % n
    T5to2 = zeros(1, n_w2); % c->n
    % conduction into sec 4
    T1to4 = zeros(n_1/2, 1); % f->n
    T2to4 = zeros(n_h2, 1); % n
    T5to4 = zeros(1, n_w4); % c->n
    T7to4 = zeros(n/2, 1); % c->n
    % conduction into sec 5
    T2to5 = zeros(1, c_w2); % n->c
    T4to5 = zeros(1, c_w4); % n->c
    T7to5 = zeros(c_h5, 1); % c
    % conduction into sec 7
    T4to7 = zeros(c_h4/2, 1); % n->c
    T5to7 = zeros(c_h5, 1); % c


    %% TIME ITERATION
    iter = 1;
    %starttime = tic;
    
    while iter <= nt

        % node temps
        Ti1 = dT1 + T0;
        Ti2 = dT2 + T0;
        Ti4 = dT4 + T0;
        Ti5 = dT5 + T0;
        Ti7 = dT7 + T0;
        Ti14 = Ti1.*Ti1.*Ti1.*Ti1;
        Ti24 = Ti2.*Ti2.*Ti2.*Ti2;
        Ti44 = Ti4.*Ti4.*Ti4.*Ti4;
        Ti54 = Ti5.*Ti5.*Ti5.*Ti5;
        Ti74 = Ti7.*Ti7.*Ti7.*Ti7;

        % SECTION 1 #######################################################
        % boundary conditions
        % left edge radiation
        B1L(:, 1) = -radEf*(Ti14(:, 1)-T04);
        % top edge conduction 2->1
        for i = 1:n_1 % larger cells on 2->1 boundary
            T2to1(1, mf*(i-1)+1:mf*i) = Ti2(end, i);
        end
        B1T(1, :) = condfV*(T2to1-Ti1(1, :));
        % right edge conduction 4->1
        for i = 1:n_1/2 % larger cells on 4->1 boundary
            T4to1(mf*(i-1)+1:mf*i, 1) = Ti4(n_h2+i, 1);
        end
        B1R(:, end) = condfH*(T4to1-Ti1(:, end));
        % bottom edge adiabatic
        B1B(end, :) = 0;

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
        B2T(1, :) = condV*(T5to2-Ti2(1, :));
        % right edge conduction sec 4->2
        B2R(:, end) = condH*(Ti4(1:n_h2, 1)-Ti2(:, end));
        % bottom edge conduction sec 1->2
        for i = 1:n_1 % no. of large cells on 1->2 boundary
            rc = mf*(i-1)+1:mf*i;
            T1to2(1, i) = sum(Ti1(1:mf, rc), 'all')/mf^2;
        end
        B2B(end, :) = condV*(T1to2-Ti2(end, :));

        % update temp
        Tu2 = Ti2 + ...
              cond2L.*(circshift(Ti2, [0 1])-Ti2) + B2L + ...
              cond2T.*(circshift(Ti2, [1 0])-Ti2) + B2T + ...
              cond2R.*(circshift(Ti2, [0 -1])-Ti2) + B2R + ...
              cond2B.*(circshift(Ti2, [-1 0])-Ti2) + B2B + ...
              pro2 + ...
              -rad*(Ti24-T04);
        dT2 = Tu2 - T0;

        % SECTION 4 #######################################################
        % boundary conditions
        % left edge conduction 2->4
        B4L(1:n_h2, 1) = condH*(Ti2(:, end)-Ti4(1:n_h2, 1));
        % left edge conduction 1->4
        for i = 1:n_1/2 % no. of large cells on 1->4 boundary
            rr = mf*(i-1)+1:mf*i;
            T1to4(i, 1) = sum(Ti1(rr, (end-mf)+1:end), 'all')/mf^2;
        end
        B4L(n_h2+1:end, 1) = condH*(T1to4-Ti4(n_h2+1:end, 1));
        % top edge conduction 5->4
        for i = 1:c_w4 % larger nodes on 5->4 boundary
            T5to4(1, mc*(i-1)+1:mc*i) = Ti5(end, c_w2+i);
        end
        B4T(1, :) = condV*(T5to4-Ti4(1, :));
        % right edge conduction 7->4
        for i = 1:c_h4/2 % larger nodes on 7->4 boundary
            T7to4(mc*(i-1)+1:mc*i, 1) = Ti7(c_h5+i, 1);
        end
        B4R(:, end) = condH*(T7to4-Ti4(:, end));
        % bottom edge adiabatic
        B4B(end, :) = 0;

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
        B5R(:, end) = condcH*(T7to5-Ti5(:, end));
        % bottom edge conduction 2->5
        for i = 1:c_w2 % larger nodes on 2->5 boundary
            rc = mc*(i-1)+1:mc*i;
            T2to5(1, i) = sum(Ti2(1:mc, rc), 'all')/mc^2;
        end
        B5B(end, 1:c_w2) = condcV*(T2to5-Ti5(end, 1:c_w2));
        % bottom edge conduction 4->5
        for i = 1:c_w4 % larger nodes on 4->5 boundary
            rc = mc*(i-1)+1:mc*i;
            T4to5(1, i) = sum(Ti4(1:mc, rc), 'all')/mc^2;
        end
        B5B(end, c_w2+1:end) = condcV*(T4to5-Ti5(end, c_w2+1:end));

        % update temp
        Tu5 = Ti5 + ...
              cond5L.*(circshift(Ti5, [0 1])-Ti5) + B5L + ...
              cond5T.*(circshift(Ti5, [1 0])-Ti5) + B5T + ...
              cond5R.*(circshift(Ti5, [0 -1])-Ti5) + B5R + ...
              cond5B.*(circshift(Ti5, [-1 0])-Ti5) + B5B + ...
              -rad*(Ti54-T04);
        dT5 = Tu5 - T0;

        % SECTION 7 #######################################################
        % boundary conditions
        % left edge conduction 5->7
        T5to7 = Ti5(:, end);
        B7L(1:c_h5, 1) = condcH*(T5to7-Ti7(1:c_h5, 1));
        % left edge conduction 4->7
        for i = 1:c_h4/2 % larger nodes on 4->7 boundary
            rr = mc*(i-1)+1:mc*i;
            T4to7(i, 1) = sum(Ti4(rr, (end-mc)+1:end), 'all')/mc^2;
        end
        B7L(c_h5+1:end, 1) = condcH*(T4to7-Ti7(c_h5+1:end, 1));
        % top edge radiation
        B7T(1, :) = -radEc*(Ti74(1, :)-T04);
        % right edge radiation
        B7R(:, end) = -radEc*(Ti74(:, end)-T04);
        % bottom edge adiabatic
        B7B(end, :) = 0;

        % update temp
        Tu7 = Ti7 + ...
              cond7L.*(circshift(Ti7, [0 1])-Ti7) + B7L + ...
              cond7T.*(circshift(Ti7, [1 0])-Ti7) + B7T + ...
              cond7R.*(circshift(Ti7, [0 -1])-Ti7) + B7R + ...
              cond7B.*(circshift(Ti7, [-1 0])-Ti7) + B7B + ...
              -rad*(Ti74-T04);
        dT7 = Tu7 - T0;

        % SAVE ############################################################
        if rem(iter, f2f) == 0
            fn = iter/f2f; % frame number
            dT1out(:, :, fn) = dT1;
            dT2out(:, :, fn) = dT2;
            dT4out(:, :, fn) = dT4;
            dT5out(:, :, fn) = dT5;
            dT7out(:, :, fn) = dT7;
            fprintf('Frame %d captured (%.2f%%)\n', fn, 100*iter/nt)
        end

        iter = iter + 1;

    end

end

