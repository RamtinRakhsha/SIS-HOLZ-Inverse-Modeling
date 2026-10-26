%% generalized_reanalysis_original_data.m
% Generalized reanalysis of the ORIGINAL SIS-HOLZ experimental phase data.
%
% This file deliberately keeps two calculations separate:
%
%   (1) HISTORICAL REFERENCE
%       Reproduces the previously used analytical bell-profile calculation
%       from NoBiprism_DBI_bell_split_USED4paper.m.
%       That old plotting script contains t = 200 nm.
%
%   (2) NEW GENERALIZED INVERSION
%       Fits the ORIGINAL ExPhase.txt data using the generalized numerical
%       forward model and the final published experimental thickness
%       t = 166 nm.
%
% The historical curve is used only for comparison. It is NOT used as
% input to the generalized inversion.
%
% Keep ExPhase.txt in the same folder as this M-file.

clear; close all; clc;

%% ============================================================
% 1. ORIGINAL EXPERIMENTAL DATA AND CALIBRATION
% =============================================================

ExPhase = load('ExPhase.txt');

pixel = ExPhase(:,1);
phaseExp = ExPhase(:,2);

a = 0.5472;                 % nm, value used in the original script
g = 8/a;                    % |g_008|, nm^-1

% Exact coordinate transformation used in the original journal M-file
pixelCentered = pixel - 108/2;
sExp = 0.8*pixelCentered/5400;
thetaExp = sExp/g;

fprintf('\nOriginal experimental dataset\n');
fprintf('-----------------------------\n');
fprintf('Number of points   : %d\n',numel(phaseExp));
fprintf('|g_008|             : %.9f nm^-1\n',g);
fprintf('Delta-theta range  : [%+.8e, %+.8e] rad\n', ...
    min(thetaExp),max(thetaExp));

%% ============================================================
% 2. HISTORICAL BELL-PROFILE REFERENCE
% =============================================================
%
% This reproduces the old journal plotting script as closely as possible.
% IMPORTANT: the old M-file contains t = 200 nm, whereas the final paper
% reports t = 166 nm for the experimental ROI. Therefore this curve is kept
% strictly as the historical plotted reference.

tHistorical = 200;          % nm, exactly as in the old M-file
alphaHistorical = 0.0325;   % nm
zeta = 250;
gamma = 50;
gamma1 = 0;
r = 1;
u = 1;
sBound = 201;

sigma = 2*pi*g*alphaHistorical;

phi_1 = zeros(1,sBound);
phi_2 = zeros(1,sBound);
phiOld = zeros(1,sBound);
phaseOldRaw = zeros(1,sBound);

for i = 1:sBound

    s_i = (i-1)/10000 - (sBound-1)/20000;
    OMEG = 2*r*u*gamma;

    phi_1(i) = besselj(0,sigma) * ...
        ( exp( 1i*pi*s_i*(tHistorical/2-(OMEG+gamma1))) + ...
          exp(-1i*pi*s_i*(tHistorical/2-(OMEG-gamma1))) ) * ...
        tHistorical*sinc(pi*s_i*tHistorical/2)/2;

    for n = 1:40
        phi_2(i) = phi_2(i) + ...
            (1i)^n*besselj(n,sigma) * ...
            ( cos(pi*(s_i-n/tHistorical) * ...
                  (tHistorical/2-(OMEG+gamma1))) * ...
              tHistorical*sinc(pi*(tHistorical/2)*(s_i-n/tHistorical)) ...
            + cos(pi*(s_i+n/tHistorical) * ...
                  (tHistorical/2-(OMEG-gamma1))) * ...
              tHistorical*sinc(pi*(tHistorical/2)*(s_i+n/tHistorical)) );
    end

    phiOld(i) = (1i*pi/zeta)*(phi_1(i)+phi_2(i));

    % Preserve the exact historical phase definition.
    phaseOldRaw(i) = atan(imag(phiOld(i))/real(phiOld(i)));
end

sHistorical = -(sBound-1)/20000 : 0.0001 : (sBound-1)/20000;

phaseOldA = phaseOldRaw;
phaseOldA(8:18)    = phaseOldA(8:18)    + pi;
phaseOldA(85:117)  = phaseOldA(85:117)  + pi;
phaseOldA(184:194) = phaseOldA(184:194) + pi;

phaseHistorical = pi + 0.1*unwrap(phaseOldA) + pi/2.6;
thetaHistorical = sHistorical/g;

%% ============================================================
% 3. NEW GENERALIZED INVERSION
% =============================================================
%
% Final published experimental condition:
%   t = 166 nm, 008 HOLZ line, R_old = 0.0325 nm.
%
% Generalized displacement:
%
%   u(z) = sum a_n sin[n*pi*(z/t + 1/2)]
%
% This basis vanishes at both surfaces and does not prescribe the overall
% profile shape.

t = 166;                    % nm, published experimental ROI thickness
N = 5;                      % low-order basis for experimental fit
Nz = 1601;
z = linspace(-t/2,t/2,Nz);

% Precompute propagation matrix at the experimental s values
P = exp(-2*pi*1i*(sExp(:)*z(:).'));

idxMinus = z <= 0;
idxPlus  = z >= 0;

% Fit parameters:
% p = [a1 ... aN, phi0, c1]
%
% a_n : pm
% phi0: arbitrary constant phase reference
% c1  : residual linear instrumental phase ramp versus theta/(1e-3 rad)
%
% No quadratic phase term is included in this first fit.

p0 = [32.5 zeros(1,N-1) 4.5 0];

lb = [-80*ones(1,N) -20 -20];
ub = [ 80*ones(1,N)  20  20];

% Deterministic multi-start initial guesses
startCoeffs = [
     32.5   0     0     0     0
     25     5     0     0     0
     35    -5     5     0     0
     20    10    -5     2     0
     40   -10     5    -2     0
];

bestCost = inf;
bestP = p0;

for kStart = 1:size(startCoeffs,1)

    pStart = [startCoeffs(kStart,1:N) 4.5 0];

    if exist('lsqnonlin','file') == 2

        opts = optimoptions('lsqnonlin', ...
            'Display','off', ...
            'MaxIterations',1000, ...
            'MaxFunctionEvaluations',10000, ...
            'FunctionTolerance',1e-12, ...
            'StepTolerance',1e-12);

        pTry = lsqnonlin( ...
            @(p) phaseResidualGeneralized( ...
                p,phaseExp,thetaExp,P,z,t,g,N,idxMinus,idxPlus), ...
            pStart,lb,ub,opts);

    else

        obj = @(p) sum(phaseResidualGeneralized( ...
            p,phaseExp,thetaExp,P,z,t,g,N,idxMinus,idxPlus).^2) ...
            + boundPenalty(p,lb,ub);

        opts = optimset('Display','off', ...
            'MaxIter',5000, ...
            'MaxFunEvals',20000, ...
            'TolX',1e-10, ...
            'TolFun',1e-12);

        pTry = fminsearch(obj,pStart,opts);
    end

    rTry = phaseResidualGeneralized( ...
        pTry,phaseExp,thetaExp,P,z,t,g,N,idxMinus,idxPlus);

    costTry = sum(rTry.^2);

    if costTry < bestCost
        bestCost = costTry;
        bestP = pTry;
    end
end

pFit = bestP;

[phaseFit,uFit_nm,PsiFit,Aminus,Aplus] = generalizedPhaseModel( ...
    pFit,thetaExp,P,z,t,g,N,idxMinus,idxPlus);

residual = phaseExp - phaseFit;
phaseRMSE = sqrt(mean(residual.^2));

aFit_pm = pFit(1:N);
uFit_pm = 1000*uFit_nm;
uPeak_pm = max(abs(uFit_pm));

%% ============================================================
% 4. PREVIOUSLY REPORTED FINAL BELL PROFILE
% =============================================================
%
% Final paper:
% u_old(z) = R [cos(2*pi*z/t)+1]/2
% R = 32.5 pm, t = 166 nm

Rold_pm = 32.5;
uOldFinal_pm = (Rold_pm/2) .* (cos(2*pi*z/t)+1);

%% ============================================================
% 5. REPORT
% =============================================================

fprintf('\nGeneralized experimental inversion\n');
fprintf('----------------------------------\n');
fprintf('Thickness used              : %.2f nm\n',t);
fprintf('Number of basis functions   : %d\n',N);
fprintf('Phase RMSE                  : %.6f rad\n',phaseRMSE);
fprintf('Recovered peak |u|          : %.4f pm\n',uPeak_pm);
fprintf('Published bell peak         : %.4f pm\n',Rold_pm);
fprintf('Phase offset                : %.6f rad\n',pFit(N+1));
fprintf('Linear phase-ramp parameter : %.6f\n',pFit(N+2));

fprintf('\nRecovered coefficients:\n');
for n = 1:N
    fprintf('a_%d = %+12.6f pm\n',n,aFit_pm(n));
end

%% ============================================================
% 6. PUBLICATION FIGURES: EXPORT PANELS SEPARATELY
% =============================================================
%
% The panel letters and short panel descriptions are added in LaTeX
% using the subfigure environment, matching the format used for Figure 4.
% Therefore, no MATLAB titles are added here.

%% -------------------------
% Panel (a): Experimental phase fit
% -------------------------
figA = figure('Color','w', ...
    'Units','inches', ...
    'Position',[0.5 0.5 5.0 5.0]);

plot(thetaExp,phaseExp,'or', ...
    'LineWidth',1.4, ...
    'MarkerSize',6.5, ...
    'MarkerFaceColor','none', ...
    'DisplayName','Original experiment');
hold on

plot(thetaExp,phaseFit,'-k', ...
    'LineWidth',2.0, ...
    'DisplayName','Generalized fit');

xlabel('$\Delta\theta$ (rad)', 'Interpreter','latex')
ylabel('Phase (rad)')

xlim([-7e-4 7e-4])
ylim([4.2 4.8])

xticks((-6:2:6)*1e-4)
yticks(4.2:0.1:4.8)

lgdA = legend('Location','none');
grid on
box on
pbaspect([1 1 1])

% Center the legend horizontally and place it slightly below the top border.
drawnow;
axA = gca;
lgdA.Units = 'normalized';
axPosA = axA.Position;
lgdPosA = lgdA.Position;
lgdPosA(1) = axPosA(1) + (axPosA(3)-lgdPosA(3))/2;
lgdPosA(2) = axPosA(2) + axPosA(4) - lgdPosA(4) - 0.025;
lgdA.Position = lgdPosA;

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',14, ...
    'LineWidth',1.2, ...
    'TickDir','in');

drawnow;
exportgraphics(figA,'ExpFit_a.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');


%% -------------------------
% Panel (b): Recovered displacement
% -------------------------
figB = figure('Color','w', ...
    'Units','inches', ...
    'Position',[0.5 0.5 5.0 5.0]);

plot(z/t,uFit_pm,'-k', ...
    'LineWidth',2.2, ...
    'DisplayName','Generalized recovery');
hold on

plot(z/t,uOldFinal_pm,'--', ...
    'LineWidth',2.0, ...
    'DisplayName','Published 32.5 pm bell profile');

xlabel('Normalized depth, $z/t$', 'Interpreter','latex')
ylabel('Projected displacement, $u$ (pm)', 'Interpreter','latex')

xlim([-0.5 0.5])
ylim([0 70])

lgdB = legend('Location','none');
grid on
box on
pbaspect([1 1 1])

% Center the legend horizontally and place it slightly below the top border.
drawnow;
axB = gca;
lgdB.Units = 'normalized';
axPosB = axB.Position;
lgdPosB = lgdB.Position;
lgdPosB(1) = axPosB(1) + (axPosB(3)-lgdPosB(3))/2;
lgdPosB(2) = axPosB(2) + axPosB(4) - lgdPosB(4) - 0.025;
lgdB.Position = lgdPosB;

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',14, ...
    'LineWidth',1.2, ...
    'TickDir','in');

drawnow;
exportgraphics(figB,'ExpFit_b.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');


%% -------------------------
% Panel (c): Phase residual
% -------------------------
figC = figure('Color','w', ...
    'Units','inches', ...
    'Position',[0.5 0.5 5.0 5.0]);

plot(thetaExp,residual,'-o', ...
    'LineWidth',1.4, ...
    'MarkerSize',5.5, ...
    'MarkerFaceColor','none');

hold on
yline(0,'--');

xlabel('$\Delta\theta$ (rad)', 'Interpreter','latex')
ylabel('Phase residual (rad)')

xlim([-7e-4 7e-4])
ylim([-pi/8-0.1 pi/8+0.1])

xticks((-6:2:6)*1e-4)

ax3 = gca;
ax3.YTick = [-pi/8 -pi/16 0 pi/16 pi/8];
ax3.YTickLabel = {'$-\pi/8$','$-\pi/16$','$0$','$\pi/16$','$\pi/8$'};
ax3.TickLabelInterpreter = 'latex';

grid on
box on
pbaspect([1 1 1])

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',14, ...
    'LineWidth',1.2, ...
    'TickDir','in');

drawnow;
exportgraphics(figC,'ExpFit_c.pdf', ...
    'ContentType','vector', ...
    'BackgroundColor','white');

% Optional PNG exports:
% exportgraphics(figA,'ExpFit_a.png','Resolution',600,'BackgroundColor','white');
% exportgraphics(figB,'ExpFit_b.png','Resolution',600,'BackgroundColor','white');
% exportgraphics(figC,'ExpFit_c.png','Resolution',600,'BackgroundColor','white');

%% Save result
result.thetaExp_rad = thetaExp;
result.phaseExp_rad = phaseExp;
result.phaseFit_rad = phaseFit;
result.residual_rad = residual;
result.z_nm = z;
result.uFit_pm = uFit_pm;
result.uOldBell_pm = uOldFinal_pm;
result.coefficients_pm = aFit_pm;
result.phaseRMSE_rad = phaseRMSE;
result.peakDisplacement_pm = uPeak_pm;
result.thickness_nm = t;
result.g008_nm_inv = g;
result.historicalTheta_rad = thetaHistorical;
result.historicalPhase_rad = phaseHistorical;

save('generalized_reanalysis_original_data_result.mat','result');

%% ============================================================
% LOCAL FUNCTIONS
% =============================================================

function r = phaseResidualGeneralized( ...
    p,phaseExp,thetaExp,P,z,t,g,N,idxMinus,idxPlus)

    phaseModel = generalizedPhaseModel( ...
        p,thetaExp,P,z,t,g,N,idxMinus,idxPlus);

    r = phaseModel - phaseExp;
end

function [phaseModel,u_nm,Psi,Aminus,Aplus] = generalizedPhaseModel( ...
    p,thetaExp,P,z,t,g,N,idxMinus,idxPlus)

    coeff_pm = p(1:N);
    phi0 = p(N+1);
    c1 = p(N+2);

    % Generalized through-thickness displacement profile
    u_nm = zeros(size(z));

    for n = 1:N
        basis = sin(n*pi*(z/t + 1/2));
        u_nm = u_nm + coeff_pm(n)*1e-3*basis;
    end

    % Displacement-dependent scattering phase from the new manuscript
    depthFactor = exp(-2*pi*1i*g*u_nm);

    integrand = P .* depthFactor;

    % Contributions from the two specimen halves
    Aminus = trapz(z(idxMinus),integrand(:,idxMinus),2);
    Aplus  = trapz(z(idxPlus), integrand(:,idxPlus), 2);

    % Ideal-overlap resultant wave for the first experimental reanalysis.
    % Any remaining constant phase origin and linear instrumental ramp are
    % fitted as nuisance parameters below.
    Psi = Aminus + Aplus;

    phiRaw = unwrap(angle(Psi));

    thetaScaled = thetaExp(:)/1e-3;

    phaseModel = phiRaw + phi0 + c1*thetaScaled;
end

function pen = boundPenalty(p,lb,ub)

    below = max(lb-p,0);
    above = max(p-ub,0);

    pen = 1e6*sum(below.^2 + above.^2);
end

