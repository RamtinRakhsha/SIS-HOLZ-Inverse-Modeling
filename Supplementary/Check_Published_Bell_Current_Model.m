%% Check_Published_Bell_Current_Model.m
% =========================================================================
% PRIORITY CHECK BEFORE SUBMISSION
%
% Purpose:
%   Compare the previously published 32.5 pm surface-relaxed bell profile
%   with the experimental phase data using EXACTLY the same numerical
%   forward-model convention currently used for the N=5 reanalysis.
%
% Fair comparison:
%   - The 32.5 pm bell displacement profile is held fixed.
%   - Only the nuisance phase offset (phi0) and linear phase ramp (c1)
%     are fitted to the experimental phase.
%
% Diagnostic:
%   The script also finds the best bell amplitude R under the same current
%   forward model. This is diagnostic only and does not replace the
%   published 32.5 pm comparison.
%
% The current N=5 coefficients reported in the manuscript are included
% only as a reference curve to verify that the current forward-model
% implementation reproduces the reported ~0.0145 rad RMSE.
%
% Outputs:
%   Published_Bell_Current_Model_Check.pdf
%   Published_Bell_Current_Model_Check.csv
% =========================================================================

clear;
close all;
clc;

%% ============================================================
% 1. LOAD EXPERIMENTAL DATA
% =============================================================

scriptPath = mfilename('fullpath');

if isempty(scriptPath)
    scriptFolder = pwd;
else
    scriptFolder = fileparts(scriptPath);
end

candidateFiles = {
    fullfile(scriptFolder,'ExPhase.txt')
    fullfile(scriptFolder,'..','Data','ExPhase.txt')
    fullfile(scriptFolder,'Data','ExPhase.txt')
    };

dataFile = '';

for k = 1:numel(candidateFiles)
    if isfile(candidateFiles{k})
        dataFile = candidateFiles{k};
        break
    end
end

if isempty(dataFile)
    error(['ExPhase.txt was not found. Put it beside this script or in ' ...
           '../Data/ExPhase.txt.']);
end

ExPhase = load(dataFile);

pixel = ExPhase(:,1);
phaseExp = ExPhase(:,2);

%% ============================================================
% 2. SAME EXPERIMENTAL CALIBRATION USED IN CURRENT REANALYSIS
% =============================================================

a = 0.5472;                 % nm
g = 8/a;                    % |g_008|, nm^-1

pixelCentered = pixel - 108/2;
sExp = 0.8*pixelCentered/5400;
thetaExp = sExp/g;

t = 166;                    % nm, published ROI thickness

Nz = 1601;
z = linspace(-t/2,t/2,Nz);

P = exp(-2*pi*1i*(sExp(:)*z(:).'));

idxMinus = z <= 0;
idxPlus  = z >= 0;

thetaScaled = thetaExp(:)/1e-3;

% Nuisance design matrix:
% phase = raw_phase + phi0 + c1*(theta/1e-3)
Xnuis = [ones(numel(thetaExp),1), thetaScaled];

%% ============================================================
% 3. PUBLISHED 32.5 pm SURFACE-RELAXED BELL
% =============================================================
%
% Published profile:
%
% u_old(z) = R/2 * [cos(2*pi*z/t) + 1]
%
% with R = 32.5 pm.

Rpublished_pm = 32.5;

uPublished_pm = ...
    (Rpublished_pm/2) .* ...
    (cos(2*pi*z/t)+1);

uPublished_nm = 1e-3*uPublished_pm;

phiRawPublished = ...
    currentForwardPhase( ...
        uPublished_nm,P,z,g,idxMinus,idxPlus);

% Fit ONLY the nuisance offset and linear ramp.
betaPublished = ...
    Xnuis \ (phaseExp(:)-phiRawPublished(:));

phiPublished = ...
    phiRawPublished(:) + ...
    Xnuis*betaPublished;

resPublished = ...
    phaseExp(:)-phiPublished;

rmsePublished = ...
    sqrt(mean(resPublished.^2));

%% ============================================================
% 4. CURRENT N=5 REFERENCE CURVE
% =============================================================
%
% Coefficients reported in the current manuscript.
% This is included to verify that this standalone script reproduces the
% current Figure 5 fit convention.

aN5_pm = ...
    [47.7217 -0.0061 14.5754 0.7070 -10.3383];

uN5_nm = zeros(size(z));

for n = 1:5
    uN5_nm = ...
        uN5_nm + ...
        aN5_pm(n)*1e-3 .* ...
        sin(n*pi*(z/t + 1/2));
end

phiRawN5 = ...
    currentForwardPhase( ...
        uN5_nm,P,z,g,idxMinus,idxPlus);

betaN5 = ...
    Xnuis \ (phaseExp(:)-phiRawN5(:));

phiN5 = ...
    phiRawN5(:) + ...
    Xnuis*betaN5;

resN5 = ...
    phaseExp(:)-phiN5;

rmseN5 = ...
    sqrt(mean(resN5.^2));

%% ============================================================
% 5. DIAGNOSTIC: BEST BELL AMPLITUDE UNDER CURRENT MODEL
% =============================================================
%
% This does NOT alter the published comparison.
% It answers:
%   "If the profile is forced to retain the published bell shape, what
%    amplitude gives the lowest full-profile phase RMSE under the current
%    forward model?"

objBell = @(Rpm) bellObjective( ...
    Rpm,phaseExp,Xnuis,P,z,t,g,idxMinus,idxPlus);

opts = optimset( ...
    'Display','off', ...
    'TolX',1e-8);

[Rbest_pm,~] = ...
    fminbnd(objBell,0,100,opts);

uBestBell_pm = ...
    (Rbest_pm/2) .* ...
    (cos(2*pi*z/t)+1);

phiRawBestBell = ...
    currentForwardPhase( ...
        1e-3*uBestBell_pm,P,z,g,idxMinus,idxPlus);

betaBestBell = ...
    Xnuis \ (phaseExp(:)-phiRawBestBell(:));

phiBestBell = ...
    phiRawBestBell(:) + ...
    Xnuis*betaBestBell;

resBestBell = ...
    phaseExp(:)-phiBestBell;

rmseBestBell = ...
    sqrt(mean(resBestBell.^2));

%% ============================================================
% 6. REPORT
% =============================================================

fprintf('\n');
fprintf('PUBLISHED BELL / CURRENT FORWARD-MODEL CHECK\n');
fprintf('============================================\n');
fprintf('Experimental points                  : %d\n',numel(phaseExp));
fprintf('Thickness                            : %.1f nm\n',t);
fprintf('|g_008|                              : %.8f nm^-1\n',g);

fprintf('\nPublished 32.5 pm bell:\n');
fprintf('  R                                   = %.4f pm\n',Rpublished_pm);
fprintf('  fitted phi0                         = %.8f rad\n',betaPublished(1));
fprintf('  fitted c1                           = %.8f\n',betaPublished(2));
fprintf('  full-profile phase RMSE             = %.8f rad\n',rmsePublished);

fprintf('\nCurrent N=5 reference:\n');
fprintf('  phase RMSE                          = %.8f rad\n',rmseN5);

fprintf('\nBest forced-bell diagnostic:\n');
fprintf('  best bell amplitude R               = %.4f pm\n',Rbest_pm);
fprintf('  full-profile phase RMSE             = %.8f rad\n',rmseBestBell);

fprintf('\nRMSE ratios:\n');
fprintf('  published bell / N=5                = %.2f\n', ...
    rmsePublished/rmseN5);
fprintf('  best forced bell / N=5              = %.2f\n', ...
    rmseBestBell/rmseN5);

%% ============================================================
% 7. OVERLAY FIGURE
% =============================================================

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[0.5 0.5 6.2 5.2]);

plot( ...
    thetaExp,phaseExp,'or', ...
    'LineWidth',1.3, ...
    'MarkerSize',6.0, ...
    'MarkerFaceColor','none', ...
    'DisplayName','Experiment');

hold on;

plot( ...
    thetaExp,phiN5,'-k', ...
    'LineWidth',2.0, ...
    'DisplayName',sprintf('$N=5$ fit, RMSE = %.4f rad',rmseN5));

plot( ...
    thetaExp,phiPublished,'--', ...
    'LineWidth',2.0, ...
    'DisplayName',sprintf('Published 32.5 pm bell, RMSE = %.4f rad', ...
    rmsePublished));

xlabel( ...
    '$\Delta\theta$ (rad)', ...
    'Interpreter','latex');

ylabel('Phase (rad)');

xlim([-7e-4 7e-4]);
ylim([4.2 4.8]);

xticks((-6:2:6)*1e-4);
yticks(4.2:0.1:4.8);

grid on;
box on;
pbaspect([1 1 1]);

legend( ...
    'Location','best', ...
    'Interpreter','latex');

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',13, ...
    'LineWidth',1.2, ...
    'TickDir','in');

exportgraphics( ...
    fig, ...
    fullfile(scriptFolder, ...
    'Published_Bell_Current_Model_Check.pdf'), ...
    'ContentType','vector', ...
    'BackgroundColor','white');

%% ============================================================
% 8. SAVE CURVES FOR RECORD
% =============================================================

T = table( ...
    thetaExp(:), ...
    phaseExp(:), ...
    phiPublished(:), ...
    phiBestBell(:), ...
    phiN5(:), ...
    resPublished(:), ...
    resBestBell(:), ...
    resN5(:), ...
    'VariableNames', { ...
    'DeltaTheta_rad', ...
    'ExperimentalPhase_rad', ...
    'Published32p5BellPhase_rad', ...
    'BestForcedBellPhase_rad', ...
    'N5Phase_rad', ...
    'Published32p5BellResidual_rad', ...
    'BestForcedBellResidual_rad', ...
    'N5Residual_rad'});

writetable( ...
    T, ...
    fullfile(scriptFolder, ...
    'Published_Bell_Current_Model_Check.csv'));

fprintf('\nSaved:\n');
fprintf('  Published_Bell_Current_Model_Check.pdf\n');
fprintf('  Published_Bell_Current_Model_Check.csv\n');

%% ============================================================
% LOCAL FUNCTIONS
% =============================================================

function phiRaw = currentForwardPhase( ...
    u_nm,P,z,g,idxMinus,idxPlus)

    depthFactor = ...
        exp(-2*pi*1i*g*u_nm);

    integrand = ...
        P .* depthFactor;

    Aminus = ...
        trapz( ...
            z(idxMinus), ...
            integrand(:,idxMinus), ...
            2);

    Aplus = ...
        trapz( ...
            z(idxPlus), ...
            integrand(:,idxPlus), ...
            2);

    Psi = Aminus + Aplus;

    phiRaw = unwrap(angle(Psi));
end


function rmse = bellObjective( ...
    Rpm,phaseExp,Xnuis,P,z,t,g,idxMinus,idxPlus)

    u_pm = ...
        (Rpm/2) .* ...
        (cos(2*pi*z/t)+1);

    phiRaw = ...
        currentForwardPhase( ...
            1e-3*u_pm,P,z,g,idxMinus,idxPlus);

    beta = ...
        Xnuis \ (phaseExp(:)-phiRaw(:));

    phiModel = ...
        phiRaw(:) + ...
        Xnuis*beta;

    rmse = ...
        sqrt(mean((phaseExp(:)-phiModel).^2));
end
