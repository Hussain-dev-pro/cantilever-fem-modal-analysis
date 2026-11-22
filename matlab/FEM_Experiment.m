%% FEM_Experiment.m  -  validate the FE code with a steel ruler and a phone
% Rig: 300 mm steel ruler (25 x 1 mm) clamped to a table edge with 250 mm
% free length; phone taped at the tip (tip mass AND sensor).
% Record with phyphox "Acceleration (without g)", ONE tap per file, export CSV.
% The script predicts f1..f3 with the same FE code (beamFE.m), checks it
% against exact theory, then finds the measured peaks and compares.
clear; clc; close all;
E = 200e9; rho = 7850; b = 0.025; h = 0.001; L = 0.25;
m_phone = 0.19;                 % kg  -> WEIGH YOUR PHONE
A = b*h;  I = b*h^3/12;
[K,M] = beamFE(L,E,rho,A,I,40,m_phone);
fFE = sort(sqrt(eig(K,M)))/(2*pi);  fFE = fFE(1:3)';
fEx = exactCantilever(L,E,rho,A,I,m_phone,3);
fprintf('Ruler + phone: FE f1..f3 = %s Hz, exact = %s Hz\n', mat2str(round(fFE,2)), mat2str(round(fEx,2)));
fprintf('Ruler alone (for reference): %s Hz\n\n', mat2str(round(exactCantilever(L,E,rho,A,I,0,3),1)));

file = 'ruler_tap.csv';  axisCol = 4;   % phyphox: 1 time, 2 x, 3 y, 4 z (phone flat -> z)
if exist(file,'file')
    D = readmatrix(file);  t = D(:,1);  a = D(:,axisCol);  src = 'MEASURED';
else
    fs = 400;  t = (0:1/fs:10)';  a = 0.02*randn(size(t));  on = t > 0.5;  tt = t(on) - 0.5;
    fd = fFE.*[0.97 0.96 0.95];   % demo: what a slightly flexible clamp typically does
    for k = 1:3, a(on) = a(on) + (1/k)*exp(-2*pi*fd(k)*0.01*tt).*sin(2*pi*fd(k)*tt); end
    src = 'DEMO (no csv found)';
end
ok = ~isnan(t) & ~isnan(a);  t = t(ok);  a = a(ok);
fs = 1/median(diff(t));  tu = (t(1):1/fs:t(end))';  au = interp1(t,a,tu);  au = au - mean(au);
N = numel(au);  win = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));
Af = abs(fft(au.*win))/N;  Af = Af(1:floor(N/2));  f = (0:numel(Af)-1)'*fs/N;
pk = zeros(1,3);
for k = 1:3                                % search +/-25 % around each prediction
    s = f > 0.75*fFE(k) & f < 1.25*fFE(k) & f < fs/2*0.9;
    if any(s), [~,i] = max(Af.*s); pk(k) = f(i); else, pk(k) = NaN; end
end
fprintf('%s data, fs = %.0f Hz\n  mode   FE (Hz)   measured (Hz)   error\n', src, fs);
for k = 1:3, fprintf('  %d     %7.2f   %9.2f       %+5.1f %%\n', k, fFE(k), pk(k), (pk(k)-fFE(k))/fFE(k)*100); end
fig = figure('Color','w');  try, theme(fig,'light'); end %#ok<TRYNC>
semilogy(f, Af, 'k'); hold on;  for k = 1:3, xline(fFE(k),'r--',sprintf('FE f_%d',k)); end
xlim([0 min(fs/2, 1.3*fFE(3))]); grid on; xlabel('Frequency (Hz)'); ylabel('|A|');
title(['Ruler tap test: ' src]);
if ~isfolder('figures'), mkdir('figures'); end
saveas(fig, fullfile('figures','fig6_experiment.png'));
