%% Fast N=5 false-double-hump test
% =========================================================================
% PURPOSE
% -------
% Test one question only:
%
%   If the TRUE displacement profile is the published 32.5 pm
%   single-maximum bell, how often does the CURRENT N=5 experimental
%   inversion reconstruct a false double-hump when realistic phase noise
%   is added?
%
% SETTINGS
% --------
%   Truth                  : published 32.5 pm surface-relaxed bell
%   Experimental grid      : original ExPhase data grid
%   Thickness              : 166 nm
%   Reflection             : 008
%   Inverse basis          : N = 5 sine basis functions
%   Free nuisance params   : constant phase offset + linear phase ramp
%   Phase noise            : sigma_phi = 0.015 rad
%   Realizations           : 50
%   Multi-starts           : 3 per realization
%
% OUTPUTS
% -------
%   Fast_N5_false_double_hump_summary.csv
%   Fast_N5_false_double_hump_profiles.pdf
%   Fast_N5_false_double_hump_results.mat
%
% Requires Optimization Toolbox for lsqnonlin.
% No Parallel Computing Toolbox is required.
% =========================================================================

clear;
close all;
clc;

rng(12345);

%% ============================================================
% 1. USER CONTROLS
% =============================================================

N = 5;

sigmaPhi = 0.015;     % rad
nRealizations = 50;
nStarts = 3;

minProminence_pm = 1.0;

%% ============================================================
% 2. LOAD ORIGINAL EXPERIMENTAL GRID
% =============================================================

scriptPath = mfilename('fullpath');

if isempty(scriptPath)
    dataFolder = pwd;
else
    dataFolder = fileparts(scriptPath);
end

candidate1 = fullfile(dataFolder,'ExPhase.txt');
candidate2 = fullfile(dataFolder,'ExPhase(1).txt');

if isfile(candidate1)
    dataFile = candidate1;
elseif isfile(candidate2)
    dataFile = candidate2;
else
    error(['Put ExPhase.txt or ExPhase(1).txt in the same folder ' ...
           'as this M-file.']);
end

ExPhase = load(dataFile);

pixel = ExPhase(:,1);

%% ============================================================
% 3. EXPERIMENTAL CALIBRATION AND MODEL SETTINGS
% =============================================================

a = 0.5472;             % nm
g = 8/a;                % |g_008|, nm^-1

pixelCentered = pixel - 108/2;
sExp = 0.8*pixelCentered/5400;
thetaExp = sExp/g;

t = 166;                % nm
Nz = 1601;
z = linspace(-t/2,t/2,Nz);

P = exp(-2*pi*1i*(sExp(:)*z(:).'));

idxMinus = z <= 0;
idxPlus  = z >= 0;

thetaScaled = thetaExp(:)/1e-3;

%% ============================================================
% 4. N=5 SINE BASIS
% =============================================================

basis = zeros(N,Nz);

for n = 1:N
    basis(n,:) = sin(n*pi*(z/t + 1/2));
end

%% ============================================================
% 5. TRUE SINGLE-MAXIMUM BELL PROFILE
% =============================================================

Rbell_pm = 32.5;

uTrue_pm = ...
    (Rbell_pm/2) .* (cos(2*pi*z/t)+1);

uTrue_nm = uTrue_pm*1e-3;

trueNumMaxima = ...
    countMeaningfulMaxima( ...
        uTrue_pm,minProminence_pm);

fprintf('\nFAST N=5 FALSE-DOUBLE-HUMP TEST\n');
fprintf('==============================\n');
fprintf('Truth peak              : %.2f pm\n',Rbell_pm);
fprintf('True meaningful maxima  : %d\n',trueNumMaxima);
fprintf('Noise sigma             : %.3f rad\n',sigmaPhi);
fprintf('Realizations            : %d\n',nRealizations);
fprintf('Starts per realization  : %d\n\n',nStarts);

%% ============================================================
% 6. CREATE CLEAN SYNTHETIC PHASE
% =============================================================

phaseRawTrue = ...
    phaseFromProfile( ...
        uTrue_nm,P,z,g,idxMinus,idxPlus);

phi0Truth = 1.7;
c1Truth = -0.09;

phaseTrue = ...
    phaseRawTrue + ...
    phi0Truth + ...
    c1Truth*thetaScaled;

%% ============================================================
% 7. FIRST FIT THE NOISE-FREE BELL WITH N=5
% =============================================================

lb = [-80*ones(1,N) -20 -20];
ub = [ 80*ones(1,N)  20  20];

options = optimoptions( ...
    'lsqnonlin', ...
    'Display','off', ...
    'MaxIterations',1000, ...
    'MaxFunctionEvaluations',15000, ...
    'FunctionTolerance',1e-10, ...
    'StepTolerance',1e-10, ...
    'OptimalityTolerance',1e-10);

aBellProj_pm = ...
    (basis.' \ uTrue_pm.').';

noiseFreeStarts = [
    aBellProj_pm
    32.5 0 0 0 0
    zeros(1,N)
];

bestCost = inf;
pNoiseFree = [];

for q = 1:size(noiseFreeStarts,1)

    p0 = [noiseFreeStarts(q,:) phi0Truth c1Truth];

    [pTry,resnormTry] = ...
        lsqnonlin( ...
            @(p) phaseResidual( ...
                p,phaseTrue,thetaExp,P,z,t,g,N, ...
                idxMinus,idxPlus), ...
            p0,lb,ub,options);

    if resnormTry < bestCost
        bestCost = resnormTry;
        pNoiseFree = pTry;
    end
end

[phaseNoiseFreeFit,~] = ...
    phaseModel( ...
        pNoiseFree,thetaExp,P,z,t,g,N, ...
        idxMinus,idxPlus);

noiseFreeRMSE = ...
    sqrt(mean((phaseNoiseFreeFit-phaseTrue).^2));

fprintf('Noise-free N=5 reference fit\n');
fprintf('-----------------------------\n');
fprintf('Phase RMSE = %.6e rad\n',noiseFreeRMSE);

for n = 1:N
    fprintf('a_%d = %+10.5f pm\n',n,pNoiseFree(n));
end

fprintf('\n');

%% ============================================================
% 8. MONTE-CARLO LOOP
% =============================================================

peakRecovered_pm = nan(nRealizations,1);
nMaxRecovered = nan(nRealizations,1);
phaseRMSE = nan(nRealizations,1);

uRecovered_pm = nan(nRealizations,Nz);

aExperimental_pm = ...
    [47.7217 -0.0061 14.5754 0.7070 -10.3383];

for iRep = 1:nRealizations

    phaseObs = ...
        phaseTrue + ...
        sigmaPhi*randn(size(phaseTrue));

    perturbation = ...
        [2 2 2 2 2].*randn(1,N);

    starts = [
        pNoiseFree(1:N)
        aExperimental_pm
        pNoiseFree(1:N) + perturbation
    ];

    bestCost = inf;
    pBest = [];

    for q = 1:nStarts

        p0 = [
            starts(q,:) ...
            pNoiseFree(N+1) ...
            pNoiseFree(N+2)
        ];

        [pTry,resnormTry] = ...
            lsqnonlin( ...
                @(p) phaseResidual( ...
                    p,phaseObs,thetaExp,P,z,t,g,N, ...
                    idxMinus,idxPlus), ...
                p0,lb,ub,options);

        if resnormTry < bestCost
            bestCost = resnormTry;
            pBest = pTry;
        end
    end

    [phaseBest,uBest_nm] = ...
        phaseModel( ...
            pBest,thetaExp,P,z,t,g,N, ...
            idxMinus,idxPlus);

    uBest_pm = 1000*uBest_nm;

    uRecovered_pm(iRep,:) = uBest_pm;

    peakRecovered_pm(iRep) = ...
        max(abs(uBest_pm));

    nMaxRecovered(iRep) = ...
        countMeaningfulMaxima( ...
            uBest_pm,minProminence_pm);

    phaseRMSE(iRep) = ...
        sqrt(mean((phaseObs-phaseBest).^2));

    fprintf('Realization %3d/%3d: maxima = %d, peak = %.2f pm, RMSE = %.4f rad\n', ...
        iRep,nRealizations, ...
        nMaxRecovered(iRep), ...
        peakRecovered_pm(iRep), ...
        phaseRMSE(iRep));
end

%% ============================================================
% 9. PRIMARY RESULT
% =============================================================

falseDoubleRate = ...
    mean(nMaxRecovered >= 2);

fprintf('\nPRIMARY RESULT\n');
fprintf('==============\n');
fprintf('False double-hump rate = %.1f %%\n', ...
    100*falseDoubleRate);

fprintf('Recovered peak median  = %.2f pm\n', ...
    median(peakRecovered_pm));

fprintf('Recovered peak 16-84%%  = [%.2f, %.2f] pm\n', ...
    prctile(peakRecovered_pm,16), ...
    prctile(peakRecovered_pm,84));

fprintf('Median phase RMSE       = %.4f rad\n', ...
    median(phaseRMSE));

%% ============================================================
% 10. SUMMARY TABLE
% =============================================================

T = table( ...
    (1:nRealizations).', ...
    peakRecovered_pm, ...
    nMaxRecovered, ...
    phaseRMSE, ...
    'VariableNames', ...
    {'Realization','Peak_pm','NumMaxima','PhaseRMSE_rad'});

writetable( ...
    T, ...
    fullfile(dataFolder, ...
    'Fast_N5_false_double_hump_summary.csv'));

%% ============================================================
% 11. SIMPLE PROFILE FIGURE
% =============================================================

u16 = prctile(uRecovered_pm,16,1);
u50 = median(uRecovered_pm,1);
u84 = prctile(uRecovered_pm,84,1);

fig = figure( ...
    'Color','w', ...
    'Units','inches', ...
    'Position',[0.5 0.5 6 5]);

ax = axes(fig);

c = ax.ColorOrder(1,:);

fill( ...
    ax, ...
    [z/t fliplr(z/t)], ...
    [u16 fliplr(u84)], ...
    c, ...
    'FaceAlpha',0.18, ...
    'EdgeColor','none', ...
    'DisplayName','16--84\% interval');

hold(ax,'on');

plot( ...
    ax,z/t,uTrue_pm, ...
    '-k', ...
    'LineWidth',2, ...
    'DisplayName','True 32.5 pm bell');

plot( ...
    ax,z/t,u50, ...
    '-', ...
    'LineWidth',2, ...
    'DisplayName','Median N=5 recovery');

xlabel( ...
    ax, ...
    'Normalized depth, $z/t$', ...
    'Interpreter','latex');

ylabel( ...
    ax, ...
    'Projected displacement, $u$ (pm)', ...
    'Interpreter','latex');

title( ...
    ax, ...
    sprintf('False double-hump rate = %.1f\\%%', ...
    100*falseDoubleRate), ...
    'Interpreter','latex', ...
    'FontWeight','normal');

xlim(ax,[-0.5 0.5]);

grid(ax,'on');
box(ax,'on');

legend( ...
    ax, ...
    'Location','best', ...
    'Interpreter','latex');

set(ax, ...
    'FontName','Times New Roman', ...
    'FontSize',12, ...
    'LineWidth',1.0, ...
    'TickDir','in');

exportgraphics( ...
    fig, ...
    fullfile(dataFolder, ...
    'Fast_N5_false_double_hump_profiles.pdf'), ...
    'ContentType','vector', ...
    'BackgroundColor','white');

%% ============================================================
% 12. SAVE RESULTS
% =============================================================

result.N = N;
result.sigmaPhi_rad = sigmaPhi;
result.nRealizations = nRealizations;
result.nStarts = nStarts;

result.t_nm = t;
result.g_nm_inv = g;

result.z_nm = z;
result.normalizedDepth = z/t;

result.uTrue_pm = uTrue_pm;
result.trueNumMaxima = trueNumMaxima;

result.noiseFreeFitParameters = pNoiseFree;
result.noiseFreeRMSE_rad = noiseFreeRMSE;

result.uRecovered_pm = uRecovered_pm;
result.peakRecovered_pm = peakRecovered_pm;
result.nMaxRecovered = nMaxRecovered;
result.phaseRMSE_rad = phaseRMSE;

result.falseDoubleHumpRate = falseDoubleRate;

save( ...
    fullfile(dataFolder, ...
    'Fast_N5_false_double_hump_results.mat'), ...
    'result');

fprintf('\nSaved:\n');
fprintf('  Fast_N5_false_double_hump_summary.csv\n');
fprintf('  Fast_N5_false_double_hump_profiles.pdf\n');
fprintf('  Fast_N5_false_double_hump_results.mat\n');

%% ============================================================
% LOCAL FUNCTIONS
% =============================================================

function r = phaseResidual( ...
    p,phaseObs,thetaExp,P,z,t,g,N, ...
    idxMinus,idxPlus)

    phaseModelValue = ...
        phaseModel( ...
            p,thetaExp,P,z,t,g,N, ...
            idxMinus,idxPlus);

    r = ...
        phaseModelValue-phaseObs;
end


function [phaseModelValue,u_nm] = ...
    phaseModel( ...
        p,thetaExp,P,z,t,g,N, ...
        idxMinus,idxPlus)

    coeff_pm = p(1:N);
    phi0 = p(N+1);
    c1 = p(N+2);

    u_nm = zeros(size(z));

    for n = 1:N

        u_nm = ...
            u_nm + ...
            coeff_pm(n)*1e-3 .* ...
            sin(n*pi*(z/t + 1/2));
    end

    phaseRaw = ...
        phaseFromProfile( ...
            u_nm,P,z,g, ...
            idxMinus,idxPlus);

    thetaScaled = ...
        thetaExp(:)/1e-3;

    phaseModelValue = ...
        phaseRaw + ...
        phi0 + ...
        c1*thetaScaled;
end


function phaseRaw = ...
    phaseFromProfile( ...
        u_nm,P,z,g, ...
        idxMinus,idxPlus)

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

    Psi = ...
        Aminus + Aplus;

    phaseRaw = ...
        unwrap(angle(Psi));
end


function nMax = ...
    countMeaningfulMaxima( ...
        u_pm,minProminence_pm)

    nMax = 0;

    for i = 2:numel(u_pm)-1

        isLocalMax = ...
            u_pm(i) > u_pm(i-1) && ...
            u_pm(i) >= u_pm(i+1);

        if ~isLocalMax
            continue
        end

        leftMin = ...
            min(u_pm(1:i));

        rightMin = ...
            min(u_pm(i:end));

        approximateProminence = ...
            u_pm(i) - ...
            max(leftMin,rightMin);

        if approximateProminence >= minProminence_pm
            nMax = nMax + 1;
        end
    end
end
