%% INCIDENT PARTICLE

M = 938.27208816; % incident proton mass MeV/c^2
gamma = 7460.69;
beta = sqrt(1-gamma^-2);

%% MATERIAL PROPERTIES

% Aluminium
Al.A = 26.9815384; % g mol-1
Al.Z = 13;
Al.I = 166e-6; % MeV
Al.w = 1;
Al.rho = 2.6989; %g cm-3
Al.dEdx = bethe(M, Al.A, Al.Z, Al.I, Al.w, Al.rho, beta, 1)

% Carbon
C.A = 12.011;
C.Z = 6;
C.I = 78e-6;
C.w = 1;
C.rho = 1.7;
C.dEdx = bethe(M, C.A, C.Z, C.I, C.w, C.rho, beta, 1)

% NaI(Tl)
NaITl.A = [22.98976928, 126.90447, 204.3833];
NaITl.Z = [11, 53, 81];
NaITl.I = [149e-6, 515.2e-6, 798.9e-6];
NaITl.w = [.064892, .358206, .576902];
NaITl.rho = 3.67;
NaITl.dEdx = bethe(M, NaITl.A, NaITl.Z, NaITl.I, NaITl.w, NaITl.rho,...
    beta, 1)

% Chromox (99.5 % mol Al2O3, .5% mol Cr2O3)
% [Al, Cr, O]
Chromox.A = [26.982, 51.996, 15.999];
Chromox.Z = [13, 24, 8];
Chromox.I = [166e-6, 315.8e-6, 95e-6];
Chromox.w = [0.52533, 0.00509, 0.46959];
Chromox.rho = 3.85;
Chromox.dEdx = bethe(M, Chromox.A, Chromox.Z, Chromox.I, Chromox.w,...
    Chromox.rho, beta, 1)


% YAG(Ce) (95 % mol Y3Al5O12, 5 % mol Ce)
% [Y, Al, O, Ce]
YAGCe.A = [88.906, 26.982, 15.999, 140.12];
YAGCe.Z = [39, 13, 8, 58];
YAGCe.I = [374.0e-6, 166e-6, 95e-6, 565.8e-6];
YAGCe.w = [.44380, .22448, .31945, .01227];
YAGCe.rho = 4.55;
YAGCe.dEdx = bethe(M, YAGCe.A, YAGCe.Z, YAGCe.I, YAGCe.w, YAGCe.rho,...
    beta, 1)

