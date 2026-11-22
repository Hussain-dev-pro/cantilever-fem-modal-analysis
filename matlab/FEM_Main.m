%% ========================================================================
%  FEM_Main.m  -  Cantilever bracket: FEM modal analysis, verification,
%                 forced response, fatigue check and robust redesign
%  ------------------------------------------------------------------------
%  Case study: a cantilevered steel bracket (0.6 m, SHS 25x2.5) carries a
%  1.2 kg junction box on the 955 rpm pump skid. Its weld cracked.
%  Phases
%   1  Problem definition and baseline modal analysis (own FE code)
%   2  Verification: exact theory, static check, mesh convergence
%   3  Mode shapes and effective modal mass
%   4  Forced response to skid vibration, root stress, fatigue
%   5  Robust redesign: catalogue search with uncertainty
%   6  Final design check and summary
%  Needs only base MATLAB (eig, fzero). Uses beamFE.m and exactCantilever.m.
%  ========================================================================
clear; clc; close all;
figdir = 'figures'; if ~isfolder(figdir), mkdir(figdir); end

%% 1) PROBLEM DEFINITION ---------------------------------------------------
P.E = 210e9;  P.rho = 7850;  P.L = 0.6;  P.mt = 1.2;  P.zeta = 0.02;  P.ne = 40;
n_op = 955;  f1x = n_op/60;                       % running speed (Hz)
EXC  = [f1x 0.5; 2*f1x 1.5; 6*f1x 2.0];           % [Hz, base accel amplitude m/s^2]: 1x, 2x, vane pass (6 vanes)
excName = {'1x','2x','Vane pass'};
P.band = 0.10;  SEP = 0.20;  RSEP = 0.10;         % speed band, separation margins (nominal, robust)
FAT = 71;  Slim = FAT*(2e6/1e7)^(1/3)/2;          % IIW FAT71 weld: amplitude fatigue limit at 1e7 cycles
fprintf('==== PHASE 1: PROBLEM ====\n');
fprintf('Excitations: 1x = %.2f Hz, 2x = %.2f Hz, vane pass = %.1f Hz (each +/-%d%% speed band)\n', ...
        EXC(1,1), EXC(2,1), EXC(3,1), P.band*100);
fprintf('Weld fatigue limit (IIW FAT %d, 1e7 cycles): %.1f MPa stress amplitude\n', FAT, Slim);
[A0,I0,c0] = shs(0.025, 0.0025);                  % original bracket SHS 25x2.5
[f0,Phi0,K0,M0] = feModes(P,A0,I0,P.mt,P.ne);
fprintf('Original SHS 25x2.5 (%.2f kg): f = %s Hz\n', P.rho*A0*P.L, mat2str(round(f0(1:4)',1)));
fprintf('Separation margin to nearest excitation band: %.1f %%  (required >= %d %%)\n\n', ...
        100*sepMargin(f0, EXC, P.band), SEP*100);

% interference (Campbell-type) diagram
fig = newFig(); hold on;
rpm = linspace(0, 1300, 200);
for j = 1:3, plot(rpm, rpm/60*EXC(j,1)/f1x, 'LineWidth', 1.4); end
fill([0.9 1.1 1.1 0.9]*n_op, [0 0 120 120], [1 .88 .65], 'EdgeColor','none', 'FaceAlpha', .5);
for j = 1:2, yline(f0(j), 'k-', sprintf('f_%d = %.1f Hz', j, f0(j)), 'LineWidth', 1.6); end
xlabel('Pump speed (rpm)'); ylabel('Frequency (Hz)'); ylim([0 120]); grid on;
legend({'1x', '2x', 'Vane pass (6x)', 'Operating band'}, 'Location', 'northwest');
title('Interference diagram, original bracket: f_1 sits on the 2x line');
saveas(fig, fullfile(figdir,'fig1_interference.png'));

%% 2) VERIFICATION ---------------------------------------------------------
fprintf('==== PHASE 2: VERIFICATION ====\n');
for tip = [0 P.mt]
    fe = feModes(P,A0,I0,tip,P.ne);  fx = exactCantilever(P.L,P.E,P.rho,A0,I0,tip,4);
    fprintf('Tip mass %.1f kg  FE: %s  Exact: %s  max error %.1e %%\n', tip, ...
        mat2str(round(fe(1:4)',2)), mat2str(round(fx,2)), max(abs(fe(1:4)'-fx)./fx)*100);
end
[K,~] = beamFE(P.L,P.E,P.rho,A0,I0,P.ne,P.mt);  F = zeros(size(K,1),1);  F(end-1) = 1;
d = K\F;  fprintf('Static tip deflection per N: FE %.6e m, exact L^3/3EI %.6e m\n', d(end-1), P.L^3/(3*P.E*I0));
NE = [1 2 4 8 16 32 64];  fx = exactCantilever(P.L,P.E,P.rho,A0,I0,P.mt,4);
errC = nan(numel(NE),4);  errL = errC;
for i = 1:numel(NE)
    fc = feModes(P,A0,I0,P.mt,NE(i));               n = min(4,numel(fc));  errC(i,1:n) = (fc(1:n)'-fx(1:n))./fx(1:n)*100;
    fl = feModes(P,A0,I0,P.mt,NE(i),'lumped');      n = min(4,numel(fl));  errL(i,1:n) = (fl(1:n)'-fx(1:n))./fx(1:n)*100;
end
fprintf('Mesh convergence, error in f1..f4 (%%):\n  ne   consistent                         lumped\n');
for i = 1:numel(NE)
    fprintf('%4d  %s   %s\n', NE(i), sprintf('%9.2e ',errC(i,:)), sprintf('%9.2e ',errL(i,:)));
end
pC = polyfit(log(NE(3:6)), log(abs(errC(3:6,2)))', 1);  pL = polyfit(log(NE(3:6)), log(abs(errL(3:6,2)))', 1);
fprintf('Observed convergence order (mode 2): consistent %.1f, lumped %.1f\n\n', -pC(1), -pL(1));
fig = newFig();
loglog(NE, abs(errC), '-o', 'LineWidth', 1.4); hold on; set(gca,'ColorOrderIndex',1);
loglog(NE, abs(errL), '--s', 'LineWidth', 1.2); grid on;
xlabel('Number of elements'); ylabel('|Frequency error| (%)');
legend('Mode 1 consistent','Mode 2 consistent','Mode 3 consistent','Mode 4 consistent', ...
       'Mode 1 lumped','Mode 2 lumped','Mode 3 lumped','Mode 4 lumped','Location','southwest','NumColumns',2);
title(sprintf('Mesh convergence: consistent mass ~h^{%.1f}, lumped ~h^{%.1f}', -pC(1), -pL(1)));
saveas(fig, fullfile(figdir,'fig2_convergence.png'));

%% 3) MODE SHAPES AND EFFECTIVE MASS ---------------------------------------
fprintf('==== PHASE 3: MODES ====\n');
r = zeros(size(M0,1),1);  r(1:2:end) = 1;          % rigid-body vertical base motion
Lm = Phi0'*M0*r;  meff = Lm.^2/(r'*M0*r)*100;
fprintf('Effective modal mass (%%): %s, first 4 modes = %.1f %%\n\n', mat2str(round(meff(1:4)',1)), sum(meff(1:4)));
x = linspace(0, P.L, P.ne+1);
fig = newFig();
for j = 1:4
    subplot(2,2,j); w = [0; Phi0(1:2:end,j)];  w = w/max(abs(w));
    plot(x*1e3, w, 'LineWidth', 2); hold on; yline(0,'k:'); grid on; ylim([-1.1 1.1]);
    title(sprintf('Mode %d: %.1f Hz  (eff. mass %.1f%%)', j, f0(j), meff(j)));
    xlabel('Position along bracket (mm)'); ylabel('Normalised deflection');
end
saveas(fig, fullfile(figdir,'fig3_mode_shapes.png'));

%% 4) FORCED RESPONSE AND FATIGUE ------------------------------------------
fprintf('==== PHASE 4: FORCED RESPONSE (original) ====\n');
[S0,T0,parts0] = response(P,A0,I0,c0,EXC,1.0);
for j = 1:3
    fprintf('%-10s base %.1f m/s^2 -> tip %.1f m/s^2, root stress %.2f MPa\n', excName{j}, EXC(j,2), parts0(j,2), parts0(j,1));
end
[Sw0,Tw0,sw0] = worstBand(P,A0,I0,c0,EXC);
fprintf('Sum at 955 rpm: %.1f MPa; worst over speed band: %.1f MPa at %.0f rpm  (limit %.1f MPa) -> %s\n\n', ...
        S0, Sw0, sw0*n_op, Slim, passfail(Sw0 <= Slim));
fH = linspace(2, 150, 1500);  [Tr0,Sr0] = frf(P,A0,I0,c0,fH);

%% 5) ROBUST REDESIGN ------------------------------------------------------
fprintf('==== PHASE 5: REDESIGN (catalogue search) ====\n');
Bs = [20 25 30 35 40 50 60];  ts = [2 2.5 3 4 5];  C = [];
for B = Bs
    for t = ts
        if t > B/5, continue; end
        [A,I,c] = shs(B/1000, t/1000);
        f  = feModes(P,A,I,P.mt,P.ne);
        mg = sepMargin(f, EXC, P.band);   rm = robustMargin(P,A,I,EXC);
        Sw = worstBand(P,A,I,c,EXC);
        dfl = (P.mt*9.81*P.L^3/(3*P.E*I) + P.rho*A*9.81*P.L^4/(8*P.E*I))*1e3;
        ok = mg >= SEP && rm >= RSEP && Sw <= 0.5*Slim && t >= 2.5 && dfl <= P.L/300*1e3;
        C = [C; B t P.rho*A*P.L f(1) f(2) 100*mg 100*rm Sw dfl ok]; %#ok<AGROW>
    end
end
fprintf('  Section     mass kg  f1 Hz   f2 Hz  margin%%  robust%%  stress MPa  defl mm  OK\n');
for i = 1:size(C,1)
    fprintf('  SHS %2gx%-4g %6.2f %7.1f %7.0f %7.1f %8.1f %9.2f %8.2f   %s\n', C(i,1:9), passfail(C(i,10)));
end
okIdx = find(C(:,10));  [~,k] = min(C(okIdx,3) - 1e-6*C(okIdx,7));  best = C(okIdx(k),:);
fprintf('Lightest section passing ALL checks: SHS %gx%g (%.2f kg)\n', best(1), best(2), best(3));
nomOnly = C(C(:,6) >= SEP*100 & C(:,8) <= 0.5*Slim & C(:,2) >= 2.5, :);  [~,k2] = min(nomOnly(:,3));
fprintf('Note: SHS %gx%g is lighter and passes the nominal margin (%.0f%%) but FAILS the robust margin (%.0f%%)\n\n', ...
        nomOnly(k2,1), nomOnly(k2,2), nomOnly(k2,6), nomOnly(k2,7));
fig = newFig(); hold on;
scatter(C(~C(:,10),3), C(~C(:,10),7), 60, [0.75 0.2 0.2], 'x', 'LineWidth', 1.5);
scatter(C(C(:,10)==1,3), C(C(:,10)==1,7), 60, [0.1 0.55 0.2], 'filled');
plot(best(3), best(7), 'ko', 'MarkerSize', 14, 'LineWidth', 2);
yline(RSEP*100,'k--','Robust margin 10%'); grid on;
xlabel('Bracket mass (kg)'); ylabel('Robust separation margin (%)');
legend({'Fails a check','Passes all checks','Selected design'},'Location','southeast');
title(sprintf('Design space: selected SHS %gx%g', best(1), best(2)));
saveas(fig, fullfile(figdir,'fig4_design_space.png'));

%% 6) FINAL DESIGN CHECK ---------------------------------------------------
fprintf('==== PHASE 6: FINAL DESIGN ====\n');
[A1,I1,c1] = shs(best(1)/1000, best(2)/1000);
f1 = feModes(P,A1,I1,P.mt,P.ne);
[S1,T1,parts1] = response(P,A1,I1,c1,EXC,1.0);  [Sw1,Tw1] = worstBand(P,A1,I1,c1,EXC);
[Tr1,Sr1] = frf(P,A1,I1,c1,fH);
fprintf('                      Original SHS 25x2.5   Final SHS %gx%g\n', best(1), best(2));
fprintf('Mass (kg)             %8.2f              %8.2f\n', P.rho*A0*P.L, P.rho*A1*P.L);
fprintf('f1 / f2 (Hz)          %5.1f / %-5.0f         %5.1f / %-5.0f\n', f0(1), f0(2), f1(1), f1(2));
fprintf('Separation margin     %8.1f %%            %8.1f %%\n', 100*sepMargin(f0,EXC,P.band), 100*sepMargin(f1,EXC,P.band));
fprintf('Robust margin         %8.1f %%            %8.1f %%\n', 100*robustMargin(P,A0,I0,EXC), 100*robustMargin(P,A1,I1,EXC));
fprintf('Root stress, worst    %8.1f MPa          %8.2f MPa   (limit %.1f)\n', Sw0, Sw1, Slim);
fprintf('Tip accel, worst      %8.1f m/s^2        %8.1f m/s^2\n', Tw0, Tw1);
fprintf('Stress reduction: %.0f %%,  tip vibration reduction: %.0f %%\n', (1-Sw1/Sw0)*100, (1-Tw1/Tw0)*100);
fig = newFig();
subplot(2,1,1); semilogy(fH, Tr0, 'r', fH, Tr1, 'b', 'LineWidth', 1.5); hold on;
for j = 1:3, xline(EXC(j,1), 'k:'); end
ylabel('Tip / base acceleration'); grid on; legend('Original SHS 25x2.5', sprintf('Final SHS %gx%g', best(1), best(2)));
title('Transmissibility: dotted lines = 1x, 2x, vane pass');
subplot(2,1,2); semilogy(fH, Sr0, 'r', fH, Sr1, 'b', 'LineWidth', 1.5); grid on;
xlabel('Excitation frequency (Hz)'); ylabel('Root stress per 1 m/s^2 (MPa)');
saveas(fig, fullfile(figdir,'fig5_frf_compare.png'));
save('FEM_results.mat','P','EXC','A0','I0','c0','A1','I1','c1','best','f0','f1','C','Slim');
fprintf('\nSaved FEM_results.mat and figures/*.png\n');

%% ======================= LOCAL FUNCTIONS ================================
function [A,I,c] = shs(B,t)               % square hollow section (sharp corners)
    A = B^2 - (B-2*t)^2;  I = (B^4 - (B-2*t)^4)/12;  c = B/2;
end

function [f,Phi,K,M] = feModes(P,A,I,mtip,ne,massType)
    if nargin < 6, massType = 'consistent'; end
    [K,M] = beamFE(P.L,P.E,P.rho,A,I,ne,mtip,massType);
    if strcmpi(massType,'lumped')        % condense massless rotations (Guyan)
        w = 1:2:size(K,1);  t = 2:2:size(K,1);
        Kc = K(w,w) - K(w,t)*(K(t,t)\K(t,w));  Kc = (Kc + Kc')/2;
        lam = sort(real(eig(Kc, M(w,w))));  f = sqrt(lam)/(2*pi);  Phi = [];  return
    end
    [Phi,Lam] = eig(K, M);  [lam,ix] = sort(diag(Lam));  Phi = Phi(:,ix);
    Phi = Phi ./ sqrt(diag(Phi'*M*Phi))';  f = sqrt(lam)/(2*pi);
end

function m = sepMargin(f, EXC, band)      % distance of each mode (<300 Hz) from each band
    m = inf;
    for j = 1:size(EXC,1)
        lo = EXC(j,1)*(1-band);  hi = EXC(j,1)*(1+band);
        for fi = f(f < 300)'
            if fi >= lo && fi <= hi, d = 0; elseif fi < lo, d = (lo-fi)/lo; else, d = (fi-hi)/hi; end
            m = min(m, d);
        end
    end
end

function m = robustMargin(P,A,I,EXC)      % worst margin: stiffness +/-15 %, tip mass +/-20 %
    m = inf;  E0 = P.E;
    for es = [0.85 1 1.15]
        for ms = [0.8 1 1.2]
            P.E = E0*es;  f = feModes(P,A,I,P.mt*ms,P.ne);  m = min(m, sepMargin(f,EXC,P.band));
        end
    end
end

function [S,T,parts] = response(P,A,I,c,EXC,scale)
% steady harmonic base excitation, modal damping zeta; root bending stress
    [K,M,ke] = beamFE(P.L,P.E,P.rho,A,I,P.ne,P.mt);
    [Phi,Lam] = eig(K,M);  w = sqrt(diag(Lam));  Phi = Phi ./ sqrt(diag(Phi'*M*Phi))';
    Cd = M*Phi*diag(2*P.zeta*w)*Phi'*M;  r = zeros(size(K,1),1);  r(1:2:end) = 1;
    parts = zeros(size(EXC,1),2);
    for j = 1:size(EXC,1)
        W = 2*pi*EXC(j,1)*scale;  a = EXC(j,2);
        u = (K - W^2*M + 1i*W*Cd) \ (-M*r*a);           % displacement relative to base
        Fe = ke*[0;0;u(1:2)];                            % first element end forces
        parts(j,1) = abs(Fe(2))*c/I/1e6;                 % root moment -> stress (MPa)
        parts(j,2) = abs(a - W^2*u(end-1));              % absolute tip acceleration
    end
    S = sum(parts(:,1));  T = sum(parts(:,2));
end

function [Sw,Tw,sw] = worstBand(P,A,I,c,EXC)
    sc = linspace(1-P.band, 1+P.band, 41);  S = zeros(size(sc));  T = S;
    for i = 1:numel(sc), [S(i),T(i)] = response(P,A,I,c,EXC,sc(i)); end
    [Sw,i] = max(S);  Tw = max(T);  sw = sc(i);
end

function [Tr,Sr] = frf(P,A,I,c,fH)        % unit base acceleration sweep
    [K,M,ke] = beamFE(P.L,P.E,P.rho,A,I,P.ne,P.mt);
    [Phi,Lam] = eig(K,M);  w = sqrt(diag(Lam));  Phi = Phi ./ sqrt(diag(Phi'*M*Phi))';
    Cd = M*Phi*diag(2*P.zeta*w)*Phi'*M;  r = zeros(size(K,1),1);  r(1:2:end) = 1;
    Tr = zeros(size(fH));  Sr = Tr;
    for j = 1:numel(fH)
        W = 2*pi*fH(j);  u = (K - W^2*M + 1i*W*Cd) \ (-M*r);
        Tr(j) = abs(1 - W^2*u(end-1));  Fe = ke*[0;0;u(1:2)];  Sr(j) = abs(Fe(2))*c/I/1e6;
    end
end

function s = passfail(ok), if ok, s = 'PASS'; else, s = 'FAIL'; end, end

function fig = newFig()
    fig = figure('Color','w','Position',[80 80 900 560]);
    try, theme(fig,'light'); end %#ok<TRYNC>
end
