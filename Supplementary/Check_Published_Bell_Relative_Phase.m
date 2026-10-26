%% Check_Published_Bell_Relative_Phase.m
% =========================================================================
% Purpose
% -------
% Re-evaluate the published 32.5 pm surface-relaxed bell using the phase
% observable appropriate to a two-wave SIS-HOLZ hologram:
%
%       delta_phi(s) = arg[A_-(s)] - arg[A_+(s)]
%                    = arg[A_-(s) * conj(A_+(s))]
%
% rather than:
%
%       arg[A_-(s) + A_+(s)]
%
% Both signs of delta_phi are tested because the experimental phase sign
% depends on which split-line side is taken as S1 versus S2.
%
% Only the nuisance constant phase offset and linear phase ramp are fitted
% for the fixed published 32.5 pm bell.
%
% A current N=5 reference curve is evaluated using the SAME relative-phase
% observable for comparison.
%
% Outputs
% -------
% Published_Bell_Relative_Phase_Check.pdf
% Published_Bell_Relative_Phase_Check.csv
% =========================================================================

clear;
close all;
clc;

%% 1. Locate experimental data

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
    error(['ExPhase.txt was not found. Put it beside this script or ' ...
           'in ../Data/ExPhase.txt.']);
end

ExPhase = load(dataFile);

pixel = ExPhase(:,1);
phaseExp = ExPhase(:,2);

%% 2. Same current experimental calibration

a = 0.5472;                  % nm
g = 8/a;                     % |g_008|, nm^-1

pixelCentered = pixel - 108/2;
sExp = 0.8*pixelCentered/5400;
thetaExp = sExp/g;

t = 166;                     % nm
Nz = 1601;
z = linspace(-t/2,t/2,Nz);

P = exp(-2*pi*1i*(sExp(:)*z(:).'));

idxMinus = z <= 0;
idxPlus  = z >= 0;

thetaScaled = thetaExp(:)/1e-3;
Xnuis = [ones(numel(thetaExp),1), thetaScaled];

%% 3. Published 32.5 pm bell

Rpublished_pm = 32.5;

uPublished_pm = ...
    (Rpublished_pm/2) .* ...
    (cos(2*pi*z/t)+1);

[AminusPub,AplusPub] = ...
    splitWaves(1e-3*uPublished_pm,P,z,g,idxMinus,idxPlus);

phiRelPub = ...
    unwrap(angle(AminusPub .* conj(AplusPub)));

% Test both possible S1-S2 sign conventions
[phiPubPlus,betaPubPlus,rmsePubPlus] = ...
    fitNuisance(+phiRelPub,phaseExp,Xnuis);

[phiPubMinus,betaPubMinus,rmsePubMinus] = ...
    fitNuisance(-phiRelPub,phaseExp,Xnuis);

if rmsePubPlus <= rmsePubMinus
    phiPublished = phiPubPlus;
    betaPublished = betaPubPlus;
    rmsePublished = rmsePubPlus;
    publishedSign = +1;
else
    phiPublished = phiPubMinus;
    betaPublished = betaPubMinus;
    rmsePublished = rmsePubMinus;
    publishedSign = -1;
end

%% 4. Current N=5 displacement, evaluated with relative phase

aN5_pm = ...
    [47.7217 -0.0061 14.5754 0.7070 -10.3383];

uN5_nm = zeros(size(z));

for n = 1:5
    uN5_nm = ...
        uN5_nm + ...
        aN5_pm(n)*1e-3 .* ...
        sin(n*pi*(z/t + 1/2));
end

[AminusN5,AplusN5] = ...
    splitWaves(uN5_nm,P,z,g,idxMinus,idxPlus);

phiRelN5 = ...
    unwrap(angle(AminusN5 .* conj(AplusN5)));

[phiN5Plus,betaN5Plus,rmseN5Plus] = ...
    fitNuisance(+phiRelN5,phaseExp,Xnuis);

[phiN5Minus,betaN5Minus,rmseN5Minus] = ...
    fitNuisance(-phiRelN5,phaseExp,Xnuis);

if rmseN5Plus <= rmseN5Minus
    phiN5 = phiN5Plus;
    betaN5 = betaN5Plus;
    rmseN5 = rmseN5Plus;
    n5Sign = +1;
else
    phiN5 = phiN5Minus;
    betaN5 = betaN5Minus;
    rmseN5 = rmseN5Minus;
    n5Sign = -1;
end

%% 5. Diagnostic: optimize bell amplitude with relative-phase observable

objBell = @(Rpm) bellRelativeObjective( ...
    Rpm,phaseExp,Xnuis,P,z,t,g,idxMinus,idxPlus);

opts = optimset('Display','off','TolX',1e-8);

[Rbest_pm,rmseBestBell] = ...
    fminbnd(objBell,0,100,opts);

uBestBell_pm = ...
    (Rbest_pm/2) .* ...
    (cos(2*pi*z/t)+1);

[AminusBest,AplusBest] = ...
    splitWaves(1e-3*uBestBell_pm,P,z,g,idxMinus,idxPlus);

phiRelBest = ...
    unwrap(angle(AminusBest .* conj(AplusBest)));

[phiBestPlus,~,rmseBestPlus] = ...
    fitNuisance(+phiRelBest,phaseExp,Xnuis);

[phiBestMinus,~,rmseBestMinus] = ...
    fitNuisance(-phiRelBest,phaseExp,Xnuis);

if rmseBestPlus <= rmseBestMinus
    phiBestBell = phiBestPlus;
    bestSign = +1;
else
    phiBestBell = phiBestMinus;
    bestSign = -1;
end

%% 6. Report

fprintf('\n');
fprintf('PUBLISHED BELL - RELATIVE SIS-HOLZ PHASE CHECK\n');
fprintf('===============================================\n');
fprintf('Experimental points                 : %d\n',numel(phaseExp));
fprintf('Thickness                           : %.1f nm\n',t);
fprintf('|g_008|                             : %.8f nm^-1\n',g);

fprintf('\nPublished 32.5 pm bell:\n');
fprintf('  RMSE for +(arg A- - arg A+)       = %.8f rad\n',rmsePubPlus);
fprintf('  RMSE for -(arg A- - arg A+)       = %.8f rad\n',rmsePubMinus);
fprintf('  selected sign                     = %+d\n',publishedSign);
fprintf('  selected RMSE                     = %.8f rad\n',rmsePublished);
fprintf('  fitted phi0                       = %.8f rad\n',betaPublished(1));
fprintf('  fitted c1                         = %.8f\n',betaPublished(2));

fprintf('\nCurrent N=5 coefficients, relative-phase observable:\n');
fprintf('  selected sign                     = %+d\n',n5Sign);
fprintf('  selected RMSE                     = %.8f rad\n',rmseN5);

fprintf('\nBest forced-bell amplitude, relative-phase observable:\n');
fprintf('  R                                 = %.4f pm\n',Rbest_pm);
fprintf('  selected sign                     = %+d\n',bestSign);
fprintf('  RMSE                              = %.8f rad\n',rmseBestBell);

%% 7. Plot

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[0.5 0.5 6.3 5.2]);

plot( ...
    thetaExp,phaseExp,'or', ...
    'LineWidth',1.3, ...
    'MarkerSize',6, ...
    'MarkerFaceColor','none', ...
    'DisplayName','Experiment');

hold on;

plot( ...
    thetaExp,phiPublished,'--', ...
    'LineWidth',2.1, ...
    'DisplayName',sprintf( ...
    'Published 32.5 pm bell, RMSE = %.4f rad',rmsePublished));

plot( ...
    thetaExp,phiN5,'-k', ...
    'LineWidth',2.0, ...
    'DisplayName',sprintf( ...
    'Current $N=5$ profile, RMSE = %.4f rad',rmseN5));

xlabel('$\Delta\theta$ (rad)','Interpreter','latex');
ylabel('Phase (rad)');

xlim([-7e-4 7e-4]);
ylim([4.2 4.8]);

xticks((-6:2:6)*1e-4);
yticks(4.2:0.1:4.8);

grid on;
box on;
pbaspect([1 1 1]);

legend('Location','best','Interpreter','latex');

set(gca, ...
    'FontName','Times New Roman', ...
    'FontSize',13, ...
    'LineWidth',1.2, ...
    'TickDir','in');

exportgraphics( ...
    fig, ...
    fullfile(scriptFolder, ...
    'Published_Bell_Relative_Phase_Check.pdf'), ...
    'ContentType','vector', ...
    'BackgroundColor','white');

%% 8. Save curves

T = table( ...
    thetaExp(:), ...
    phaseExp(:), ...
    phiPublished(:), ...
    phiBestBell(:), ...
    phiN5(:), ...
    'VariableNames', { ...
    'DeltaTheta_rad', ...
    'ExperimentalPhase_rad', ...
    'Published32p5Bell_RelativePhase_rad', ...
    'BestForcedBell_RelativePhase_rad', ...
    'N5_RelativePhase_rad'});

writetable( ...
    T, ...
    fullfile(scriptFolder, ...
    'Published_Bell_Relative_Phase_Check.csv'));

fprintf('\nSaved:\n');
fprintf('  Published_Bell_Relative_Phase_Check.pdf\n');
fprintf('  Published_Bell_Relative_Phase_Check.csv\n');

%% ============================================================
% LOCAL FUNCTIONS
% =============================================================

function [Aminus,Aplus] = splitWaves( ...
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
end


function [phiFit,beta,rmse] = fitNuisance( ...
    phiRaw,phaseExp,Xnuis)

    beta = ...
        Xnuis \ (phaseExp(:)-phiRaw(:));

    phiFit = ...
        phiRaw(:) + Xnuis*beta;

    rmse = ...
        sqrt(mean((phaseExp(:)-phiFit).^2));
end


function rmse = bellRelativeObjective( ...
    Rpm,phaseExp,Xnuis,P,z,t,g,idxMinus,idxPlus)

    u_pm = ...
        (Rpm/2) .* ...
        (cos(2*pi*z/t)+1);

    [Aminus,Aplus] = ...
        splitWaves( ...
            1e-3*u_pm,P,z,g,idxMinus,idxPlus);

    phiRel = ...
        unwrap(angle(Aminus .* conj(Aplus)));

    [~,~,rmsePlus] = ...
        fitNuisance(+phiRel,phaseExp,Xnuis);

    [~,~,rmseMinus] = ...
        fitNuisance(-phiRel,phaseExp,Xnuis);

    rmse = min(rmsePlus,rmseMinus);
end
