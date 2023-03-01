% protons on Cu
% MATERIAL CONSTANTS
A = 63.546; % g mol-1
Z = 29;
I = 322.0e-6; % MeV
w = 1;
rho = 8.96; %g cm-3

% INCIDENT PARTICLE
M_p = 938.27208816; % incident proton mass MeV/c^2
M_m = 105.6583755; % incident muon mass MeV/c^2
beta = [0.2:0.01:0.8 0.8:0.001:0.999 0.999:1e-5:0.99999 0.99999:1e-8:0.99999999];
gamma = (1-beta).^-0.5;
betagamma = beta.*gamma;

l = length(beta);
p_dEdx = zeros(1, l);
p_dEdx_d = zeros(1, l);
m_dEdx = zeros(1, l);
m_dEdx_d = zeros(1, l);

for i = 1:l
    p_dEdx(i) = bethe(M_p, A, Z, I, w, rho, beta(i), 0);
    p_dEdx_d(i) = bethe(M_p, A, Z, I, w, rho, beta(i), 1);
    m_dEdx(i) = bethe(M_m, A, Z, I, w, rho, beta(i), 0);
    m_dEdx_d(i) = bethe(M_m, A, Z, I, w, rho, beta(i), 1);
end

figure()
axes('XScale', 'log', 'YScale', 'log')
hold on

%{
loglog(log10(betagamma), log10(p_dEdx), 'g')
loglog(log10(betagamma), log10(p_dEdx_d), 'r')
loglog(log10(betagamma), log10(m_dEdx), 'g:')
loglog(log10(betagamma), log10(m_dEdx_d), 'r:')
%}


loglog(betagamma, p_dEdx, 'g')
loglog(betagamma, p_dEdx_d, 'r')
loglog(betagamma, m_dEdx, 'g:')
loglog(betagamma, m_dEdx_d, 'r:')


legend('Proton, without delta', 'Proton, with delta', 'Muon, without delta', 'Muon, with delta')
xlabel('\beta\gamma')
ylabel('-\langledE/\rhodx\rangle [MeV g^{-1} cm^2]')
ylim([1 10])
grid on
grid minor