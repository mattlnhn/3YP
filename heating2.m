function [dT1, dT2, dT3, dT4] = heating2(dlf, Lf, dl, L, theta, dt, nt, fps, mat, beam)
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

    m = round(dl/dlf);
    nf = round(Lf/dlf); % fine nodes
    n = round(L/dl); % normal nodes

    n_h2 = round(n*(L-Lf)/(2*L)); % height of sections 2 or 3 in normal nodes
    n_h12 = round(n*(L+Lf)/(2*L)); % height of section 1 + 2 or 3 in normal nodes
    n_1 = round(n*Lf/L); % height & width of section 1 in normal nodes
    n_w4 = round(n*(L-Lf)/L); % width of section 4 in normal nodes

    % section 1--fine
    dT1 = zeros(nf, nf);
    E1L = dT1; E1L(:, 1) = 1;
    E1T = dT1; E1T(1, :) = 1;
    E1R = dT1; E1R(:, end) = 1;
    E1B = dT1; E1B(end, :) = 1;

    % section 2--normal, above fine section

    dT2 = zeros(n_h2, n_1);
    E2L = dT2; E2L(:, 1) = 1;
    E2T = dT2; E2T(1, :) = 1;
    E2R = dT2; E2R(:, end) = 1;
    E2B = dT2; E2B(end, :) = 1;

    % section 3--normal, below fine section
    dT3 = zeros(n_h2, n_1);
    E3L = dT3; E3L(:, 1) = 1;
    E3T = dT3; E3T(1, :) = 1;
    E3R = dT3; E3R(:, end) = 1;
    E3B = dT3; E3B(end, :) = 1;

    % section 4--normal, right of fine section
    dT4 = zeros(n, n_w4);
    E4L = dT4; E4L(:, 1) = 1;
    E4T = dT4; E4T(1, :) = 1;
    E4R = dT4; E4R(:, end) = 1;
    E4B = dT4; E4B(end, :) = 1;

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

    % proton energy deposition coefficients
    pro1 = (dt*msp*rho_p1)/mat.c_p; % section 1
    pro2 = (dt*msp*rho_p2)/mat.c_p; % section 2
    pro3 = (dt*msp*rho_p3)/mat.c_p; % section 3
    pro4 = (dt*msp*rho_p4)/mat.c_p; % section 4

    % radiation
    % all cells front & back
    rad = (2*mat.epsilon*const_sigma*dt)/(mat.rho*mat.c_p*theta);
    % fine edge
    radEf = (mat.epsilon*const_sigma*dt)/(mat.rho*mat.c_p*dlf);
    % normal edge
    radE = (mat.epsilon*const_sigma*dt)/(mat.rho*mat.c_p*dl);

    % exclude edges of each section from inner node conduction
    cond1L = condf*not(E1L); cond1T = condf*not(E1T);
    cond1R = condf*not(E1R); cond1B = condf*not(E1B);
    cond2L = cond*not(E2L); cond2T = cond*not(E2T);
    cond2R = cond*not(E2R); cond2B = cond*not(E2B);
    cond3L = cond*not(E3L); cond3T = cond*not(E3T);
    cond3R = cond*not(E3R); cond3B = cond*not(E3B);
    cond4L = cond*not(E4L); cond4T = cond*not(E4T);
    cond4R = cond*not(E4R); cond4B = cond*not(E4B);

    % preallocate matrices for boundary conditions
    B1L = dT1; B1T = dT1; B1R = dT1; B1B = dT1;
    B2L = dT2; B2T = dT2; B2R = dT2; B2B = dT2;
    B3L = dT3; B3T = dT3; B3R = dT3; B3B = dT3;
    B4L = dT4; B4T = dT4; B4R = dT4; B4B = dT4;

    % preallocate inner boundaries
    % conduction into sec 1
    T2to1 = zeros(1, m*n_1); % lo->hi
    T3to1 = zeros(1, m*n_1); % lo->hi
    T4to1 = zeros(m*n_1, 1); % lo->hi
    % conduction into sec 2
    T1to2 = zeros(1, n_1); % hi->lo
    T4to2 = zeros(n_h2, 1); % lo
    % conduction into sec 3
    T1to3 = zeros(1, n_1); % hi->lo
    T4to3 = zeros(n_h2, 1); % lo
    % conduction into sec 4
    T1to4 = zeros(n_1, 1); % hi->lo
    T2to4 = zeros(n_h2, 1); % lo
    T3to4 = zeros(n_h2, 1); % lo

    %% TIME ITERATION
    iter = 1;
    %starttime = tic;
    
    while iter <= nt

        % node temps
        Ti1 = dT1 + T0;
        Ti2 = dT2 + T0;
        Ti3 = dT3 + T0;
        Ti4 = dT4 + T0;
        Ti14 = Ti1.*Ti1.*Ti1.*Ti1;
        Ti24 = Ti2.*Ti2.*Ti2.*Ti2;
        Ti34 = Ti3.*Ti3.*Ti3.*Ti3;
        Ti44 = Ti4.*Ti4.*Ti4.*Ti4;

        % SECTION 1 #######################################################
        % boundary conditions
        % left edge radiation
        B1L(:, 1) = -radEf*(Ti14(:, 1)-T04);
        % top edge conduction sec 2->1
        for i = 1:n_1 % no. of large cells on 2->1 boundary
            T2to1(1, m*(i-1)+1:m*i) = Ti2(end, i);
        end
        B1T(1, :) = condf*(T2to1-Ti1(1, :));
        % right edge conduction sec 4->1
        for i = 1:n_1
            T4to1(m*(i-1)+1:m*i, 1) = Ti4(n_h2+i, 1);
        end
        B1R(:, end) = condf*(T4to1-Ti1(:, end));
        % bottom edge conduction sec 3->1
        for i = 1:n_1
            T3to1(1, m*(i-1)+1:m*i) = Ti3(1, i);
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
        % top edge radiation
        B2T(1, :) = -radE*(Ti24(1, :)-T04);
        % right edge conduction sec 4->2
        B2R(:, end) = cond*(Ti4(1:n_h2, 1)-Ti2(:, end));
        % bottom edge conduction sec 1->2
        for i = 1:n_1 % no. of large cells on 1->2 boundary
            rc = m*(i-1)+1:m*i;
            T1to2(1, i) = sum(Ti1(1:m, rc), 'all')/m^2;
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
        for i = 1:n_1 % no. of large cells on 1->3 boundary
            rc = m*(i-1)+1:m*i;
            T1to3(1, i) = sum(Ti1((end-m)+1:end, rc), 'all')/m^2;
        end
        B3T(1, :) = cond*(T1to3-Ti3(1, :));
        % right edge conduction sec 4->3
        B3R(:, end) = cond*(Ti4(n_h12+1:end, 1)-Ti3(:, end));
        % bottom edge radiation
        B3B(end, :) = -radE*(Ti34(end, :)-T04);
        
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
            rr = m*(i-1)+1:m*i;
            T1to4(i, 1) = sum(Ti1(rr, (end-m)+1:end), 'all')/m^2;
        end
        B4L(n_h2+1:n_h12, 1) = cond*(T1to4-Ti4(n_h2+1:n_h12, 1));
        % left edge conduction 3->4
        B4L(n_h12+1:end, 1) = cond*(Ti3(:, end)-Ti4(n_h12+1:end, 1));
        % top edge radiation
        B4T(1, :) = -radE*(Ti44(1, :)-T04);
        % right edge radiation
        B4R(:, end) = -radE*(Ti44(:, end)-T04);
        % bottom edge radiation
        B4B(end, :) = -radE*(Ti44(end, :)-T04);

        % update temp
        Tu4 = Ti4 + ...
              cond4L.*(circshift(Ti4, [0 1])-Ti4) + B4L + ...
              cond4T.*(circshift(Ti4, [1 0])-Ti4) + B4T + ...
              cond4R.*(circshift(Ti4, [0 -1])-Ti4) + B4R + ...
              cond4B.*(circshift(Ti4, [-1 0])-Ti4) + B4B + ...
              pro4 + ...
              -rad*(Ti44-T04);
        dT4 = Tu4 - T0;

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
        iter = iter + 1

    end

end

