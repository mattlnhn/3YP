clc; clear;

Lc = 1e-2;
dlc = 1e-4;
L = 5e-3;
dl = 5e-5;
Lf = 1e-3;
dlf = 1e-5;
mf = round(dl/dlf);
mc = round(dlc/dl);

theta = 100e-6; % thickness, m
phi = 45; % degrees

mat.A = [26.982, 51.996, 15.999]; % atomic mass, g mol-1
mat.Z = [13, 24, 8]; % atomic no.
mat.I = [166e-6, 315.8e-6, 95e-6]; % mean ionisation energy, MeV
mat.w = [0.52533, 0.00509, 0.46959]; % mass fraction
mat.rho = 3.85; % density, g cm-3
mat.epsilon = .8; % emissivity coefficient
mat.k = 30; % thermal conductivity, W m-1 K-1
mat.c_p = 900; % specific heat capacity J kg-1 K-1

beam.M = 938.27208816; % particle mass, MeV c-2
beam.gamma = 7460.69; % Lorentz factor
beam.f = 11.2e3;  % bunches per second, Hz
beam.np = 2.2e11; % protons per bunch (nominal 25ns spacing)
beam.nb = 2760; % no. of bunches
%beam.beta_star = 0.15;
%beam.epsilon_n = 2.50e-6;
%beam.sigma = sqrt(beam.beta_star*beam.epsilon_n/beam.gamma);
beam.sigma = 220e-6;
%beam.pos = -4*beam.sigma; % beam centre pos w.r.t. midpoint of left edge, m

tau = .25;
dt = tau*(mat.rho*1000*mat.c_p*dlf*dlf*cosd(phi)*cosd(phi))/mat.k; % time increment, s
T = 1e0;
nt = round(T/dt);
T = dt*nt;

bprange = 3.5:.1:5.5;

for b = 1:length(bprange)

    beam.pos = -bprange(b)*beam.sigma;

    fprintf('beam pos = %d sigma, dt = %d s with %d time steps for total simulation time %d s\n', bprange(b), dt, nt, T)

    total_time = tic;
    [dT1, dT2, dT4, dT5, dT7] = heating2A(dlf, Lf, dl, L, dlc, Lc, theta, dt, nt, 1000, mat, beam, phi);
    toc(total_time)

    filename = sprintf('Chromox100um45degIP4%.1fsigma1e0s.mat', bprange(b));
    save(filename, 'dT1', 'dT2', 'dT4', 'dT5', 'dT7')

end

T = 1e-3;
nt = round(T/dt);
T = dt*nt;

for b = 1:length(bprange)

    beam.pos = -bprange(b)*beam.sigma;

    fprintf('beam pos = %d sigma, dt = %d s with %d time steps for total simulation time %d s\n', bprange(b), dt, nt, T)

    total_time = tic;
    [dT1, dT2, dT4, dT5, dT7] = heating2A(dlf, Lf, dl, L, dlc, Lc, theta, dt, nt, 1000, mat, beam, phi);
    toc(total_time)

    filename = sprintf('Chromox100um45degIP4%.1fsigma1e-3s.mat', bprange(b));
    save(filename, 'dT1', 'dT2', 'dT4', 'dT5', 'dT7')

end

%{
dT2up = kron(dT2(:, :, end-1), ones(mf));
dT4up = kron(dT4(:, :, end-1), ones(mf));

dT5up = kron(dT5(:, :, end-1), ones(mf*mc));
dT7up = kron(dT7(:, :, end-1), ones(mf*mc));

dTs1 = [dT2up; dT1(:, :, end-1)];
dTs2 = [dTs1 dT4up];
dTs3 = [dT5up; dTs2];
dTs4 = [dTs3 dT7up];
dTfinal = [dTs4; flipud(dTs4)];

figure()
imagesc(dTfinal)
colorbar
axis square
%}