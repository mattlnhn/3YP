function [dTframes] = heating2(n, L, theta, dt, nt, mat, beam, beampos)
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
    
    %% FIX FOR PARFOR
    beam.pos = beampos;

    %% MESH

    nc = n;
    nn = n*2;
    nf = n*2*2;
    nx = n*2*2*2;
    
    % fractions must add to 1
    % n must be multiple of denominator
    wc = nc*97/100;
    wn = nn/100;
    wf = nf/100;
    wx = nx/100;

    % xfine
    dlx = L/nx;
    dTx = zeros(nx/2, wx);

    Exl = dTx;
    Exl(:, 1) = 1;
    Ext = dTx;
    Ext(1, :) = 1;
    Exr = dTx;
    Exr(:, end) = 1;
    Exb = dTx;
    Exb(end, :) = 1;

    % fine
    dlf = L/nf;
    dTf = zeros(nf/2, wf);

    Efl = dTf;
    Efl(:, 1) = 1;
    Eft = dTf;
    Eft(1, :) = 1;
    Efr = dTf;
    Efr(:, end) = 1;
    Efb = dTf;
    Efb(end, :) = 1;

    % normal
    dln = L/nn;
    dTn = zeros(nn/2, wn);

    Enl = dTn;
    Enl(:, 1) = 1;
    Ent = dTn;
    Ent(1, :) = 1;
    Enr = dTn;
    Enr(:, end) = 1;
    Enb = dTn;
    Enb(end, :) = 1;

    % coarse
    dlc = L/nc;
    dTc = zeros(nc/2, wc);

    Ecl = dTc;
    Ecl(:, 1) = 1;
    Ect = dTc;
    Ect(1, :) = 1;
    Ecr = dTc;
    Ecr(:, end) = 1;
    Ecb = dTc;
    Ecb(end, :) = 1;
    

    %% ENERGY DEPOSITION
    beta = (1-beam.gamma^-2)^.5;
    % mass stopping power = <1/rho dE/dx>
    msp = bethe(beam.M, mat.A, mat.Z, mat.I, mat.w, mat.rho, beta, 1);

    rho_px = densdist(beam.pos, nx, L, 1/100, beam.np, beam.nb, ...
        beam.f, beam.sigma, beam.sigma);
    rho_pf = densdist(beam.pos-L/100, nf, L, 1/100, beam.np, beam.nb, ...
        beam.f, beam.sigma, beam.sigma);
    rho_pn = densdist(beam.pos-2*L/100, nn, L, 1/100, beam.np, ...
        beam.nb, beam.f, beam.sigma, beam.sigma);
    rho_pc = densdist(beam.pos-3*L/100, nc, L, 97/100, beam.np, ...
        beam.nb, beam.f, beam.sigma, beam.sigma);


    %% CONVERT UNITS
    % mass stopping power MeV g-1 cm2 --> J g-1 m2
    % * joules per mev * (m per cm)^2
    msp = msp * 1.6021773e-13 * (0.1)^2;
    % density g cm-3 --> kg m-3
    mat.rho = mat.rho * 1000;


    %% ANIMATION OUTPUT PARAMS
    fps = 1000; % frames per simulated second
    nframes = round((dt*fps)^-1); % no. of time steps btwn frames
    frames = round(nt/nframes); % no. of frames total
    dTframes = zeros(nn/2, nn, frames); % preallocate


    %% TIME ITERATION
    i = 1;

    % inner node conduction coefficients
    c1x = (mat.k*dt)/(mat.rho*mat.c_p*dlx*dlx); % xfine
    c1f = (mat.k*dt)/(mat.rho*mat.c_p*dlf*dlf); % fine
    c1n = (mat.k*dt)/(mat.rho*mat.c_p*dln*dln); % normal
    c1c = (mat.k*dt)/(mat.rho*mat.c_p*dlc*dlc); % coarse

    % proton energy deposition coefficients
    c3x = (dt*msp*rho_px)/mat.c_p; % xfine
    c3f = (dt*msp*rho_pf)/mat.c_p; % fine
    c3n = (dt*msp*rho_pn)/mat.c_p; % normal
    c3c = (dt*msp*rho_pc)/mat.c_p; % coarse

    % radiation
    % all cells
    c4 = (2*mat.epsilon*const_sigma*dt)/(mat.rho*mat.c_p*theta);
    % xfine boundary
    c4Ex = (2*mat.epsilon*const_sigma*dt)/(mat.rho*mat.c_p*dlx);
    % fine boundary
    c4Ef = (2*mat.epsilon*const_sigma*dt)/(mat.rho*mat.c_p*dlf);
    % normal boundary
    c4En = (2*mat.epsilon*const_sigma*dt)/(mat.rho*mat.c_p*dln);
    % coarse boundary
    c4Ec = (2*mat.epsilon*const_sigma*dt)/(mat.rho*mat.c_p*dlc);

    % pre-calculate conduction coefficient matrices
    c1xl = c1x*not(Exl); % xfine excluding left edge
    c1xt = c1x*not(Ext); % xfine excluding top edge
    c1xr = c1x*not(Exr); % xfine excluding right edge
    c1xb = c1x*not(Exb); % xfine excluding bottom edge
    c1fl = c1f*not(Efl); % fine etc.
    c1ft = c1f*not(Eft);
    c1fr = c1f*not(Efr);
    c1fb = c1f*not(Efb);
    c1nl = c1n*not(Enl); % normal etc.
    c1nt = c1n*not(Ent);
    c1nr = c1n*not(Enr);
    c1nb = c1n*not(Enb);
    c1cl = c1c*not(Ecl); % coarse etc.
    c1ct = c1c*not(Ect);
    c1cr = c1c*not(Ecr);
    c1cb = c1c*not(Ecb);

    % preallocate matrices for boundary conditions
    Bixl = dTx; % xfine left
    Bixt = dTx; % xfine top
    Bixr = dTx; % xfine right
    Bixb = dTx; % xfine bottom
    Bifl = dTf; % fine etc.
    Bift = dTf;
    Bifr = dTf;
    Bifb = dTf;
    Binl = dTn; % normal etc.
    Bint = dTn;
    Binr = dTn;
    Binb = dTn;
    Bicl = dTc; % coarse etc.
    Bict = dTc;
    Bicr = dTc;
    Bicb = dTc;

    % inner boundaries
    Tx2f = zeros(nx/2, 1); % xfine to fine
    Tf2x = zeros(nf/2, 1); % fine to xfine
    Tf2n = zeros(nf/2, 1); % fine to normal
    Tn2f = zeros(nn/2, 1); % normal to fine
    Tn2c = zeros(nn/2, 1); % normal to coarse
    Tc2n = zeros(nc/2, 1); % coarse to normal

    starttime = tic;

    while i <= nt

        % node temps
        Tixn = dTx + T0; % T^i_node xfine
        Tifn = dTf + T0; % T^i_node fine
        Tinn = dTn + T0; % T^i_node normal
        Ticn = dTc + T0; % T^i_node coarse
        Tixn4 = Tixn.*Tixn.*Tixn.*Tixn; % (T^i_node)^4 xfine
        Tifn4 = Tifn.*Tifn.*Tifn.*Tifn; % (T^i_node)^4 fine
        Tinn4 = Tinn.*Tinn.*Tinn.*Tinn; % (T^i_node)^4 normal
        Ticn4 = Ticn.*Ticn.*Ticn.*Ticn; % (T^i_node)^4 coarse

        % XFINE
        % boundary conditions
        Bixl(:, 1) = -c4Ex*(Tixn4(:, 1)-T04); % rad
        Bixt(1, :) = -c4Ex*(Tixn4(1, :)-T04); % rad
        Tx2f(1:2:end-1) = Tifn(:, 1);
        Tx2f(2:2:end) = Tifn(:, 1);
        Bixr(:, end) = c1x*(Tx2f-Tixn(:, end)); % cond
        Bixb(end, :) = -c4Ex*(Tixn4(end, :)-T04); % rad
        % update temp
        Tx = Tixn + ...
              c1xl.*(circshift(Tixn, [0 1])-Tixn) + Bixl + ...
              c1xt.*(circshift(Tixn, [1 0])-Tixn) + Bixt + ...
              c1xr.*(circshift(Tixn, [0 -1])-Tixn) + Bixr + ...
              c1xb.*(circshift(Tixn, [-1 0])-Tixn) + Bixb + ...
              c3x + ...
              -c4*(Tixn4-T04);
        dTx = Tx - T0;


        % FINE
        % boundary conditions
        Tf2x = (Tixn(1:2:end-1, end) + Tixn(2:2:end, end) + ...
            Tixn(1:2:end-1, end-1) + Tixn(2:2:end, end-1))/4;
        Bifl(:, 1) = c1f*(Tf2x-Tifn(:, 1)); % cond
        Bift(1, :) = -c4Ef*(Tifn4(1, :)-T04); % rad
        Tf2n(1:2:end-1) = Tinn(:, 1);
        Tf2n(2:2:end) = Tinn(:, 1);
        Bifr(:, end) = c1f*(Tf2n-Tifn(:, end)); % cond
        Bifb(end, :) = -c4Ef*(Tifn4(end, :)-T04); % rad

        % update temp
        Tf = Tifn + ... 
              c1fl.*(circshift(Tifn, [0 1])-Tifn) + Bifl + ...
              c1ft.*(circshift(Tifn, [1 0])-Tifn) + Bift + ...
              c1fr.*(circshift(Tifn, [0 -1])-Tifn) + Bifr + ...
              c1fb.*(circshift(Tifn, [-1 0])-Tifn) + Bifb + ...
              c3f + ...
              -c4*(Tifn4-T04);
        dTf = Tf - T0;


        % NORMAL
        % boundary conditions
        Tn2f = (Tifn(1:2:end-1, end) + Tifn(2:2:end, end) + ...
            Tifn(1:2:end-1, end-1) + Tifn(2:2:end, end-1))/4;
        Binl(:, 1) = c1n*(Tn2f-Tinn(:, 1)); % cond
        Bint(1, :) = -c4En*(Tinn4(1, :)-T04); % rad
        Tn2c(1:2:end-1) = Ticn(:, 1);
        Tn2c(2:2:end) = Ticn(:, 1);
        Binr(:, end) = c1n*(Tn2c-Tinn(:, end)); % cond
        Binb(end, :) = -c4En*(Tinn4(end, :)-T04); % rad
        % update temp
        Tn = Tinn + ... 
              c1nl.*(circshift(Tinn, [0 1])-Tinn) + Binl + ...
              c1nt.*(circshift(Tinn, [1 0])-Tinn) + Bint + ...
              c1nr.*(circshift(Tinn, [0 -1])-Tinn) + Binr + ...
              c1nb.*(circshift(Tinn, [-1 0])-Tinn) + Binb + ...
              c3n + ...
              -c4*(Tinn4-T04);
        dTn = Tn - T0;

        % COARSE
        % boundary conditions
        Tc2n = (Tinn(1:2:end-1, end) + Tinn(2:2:end, end) + ...
            Tinn(1:2:end, end-1) + Tinn(2:2:end, end-1))/4;
        Bicl(:, 1) = c1c*(Tc2n-Ticn(:, 1)); % cond
        Bict(1, :) = -c4Ec*(Ticn4(1, :)-T04); % rad
        Bicr(:, end) = -c4Ec*(Ticn4(:, end)-T04); % rad
        Bicb(end, :) = -c4Ec*(Ticn4(end, :)-T04); % rad

        % update temp
        Tc = Ticn + ... 
              c1cl.*(circshift(Ticn, [0 1])-Ticn) + Bicl + ...
              c1ct.*(circshift(Ticn, [1 0])-Ticn) + Bict + ...
              c1cr.*(circshift(Ticn, [0 -1])-Ticn) + Bicr + ...
              c1cb.*(circshift(Ticn, [-1 0])-Ticn) + Bicb + ...
              c3c + ...
              -c4*(Ticn4-T04);
        dTc = Tc - T0;

        if rem(i, nframes) == 0
            ind = i/nframes;
            dTframes(:, 1:nn/100, ind) = dTx(1:4:end, 1:4:end);
            dTframes(:, (nn/100)+1:2*nn/100, ind) = dTf(1:2:end, 1:2:end);
            dTframes(:, (2*nn/100)+1:3*nn/100, ind) = dTn;
            dTframes(:, (3*nn/100)+1:end, ind) = kron(dTc, ones(2));
            endtime = toc(starttime);
            fprintf('%.4f%% complete in %.4f s. Animation frame captured.\n', 100*i/nt, endtime)
            starttime = tic;
        end

        i = i + 1;

    end

end

