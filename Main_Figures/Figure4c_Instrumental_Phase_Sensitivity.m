% FIGURE 3c
% Sensitivity of SIS-HOLZ inverse recovery to instrumental phase curvature
%
% Synthetic data are generated using the true phase curvature a_xx,0.
% During inversion, the assumed phase curvature is varied while the
% three displacement coefficients are refitted.
%
% Outputs:
%   fig3c_instrumental_phase_sensitivity.png
%   fig3c_instrumental_phase_sensitivity_results.csv
%
% Requires Optimization Toolbox for lsqnonlin.

clear;
close all;
clc;

rng(8);

%% ===============================================================
% USER CONTROLS
% ===============================================================

fastMode = true;

% Test assumed curvature relative to the true value
curvatureRatios = [0.50 0.75 1.00 1.25 1.50];

% Known synthetic coefficients in pm
aTruePm = [32; -4; 2];

% Bounds for fitted displacement coefficients in pm
coefficientBoundPm = 60;

% Keep zero to isolate instrumental phase sensitivity
noiseFraction = 0.0;

showProgress = true;

if fastMode

    z = linspace(-0.5, 0.5, 201);
    s = linspace(-8, 8, 121);
    x = linspace(-1, 1, 81);

else

    z = linspace(-0.5, 0.5, 401);
    s = linspace(-8, 8, 241);
    x = linspace(-1, 1, 241);

end

%% ===============================================================
% GLOBAL LATEX FORMATTING
% ===============================================================

set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");

set(groot, "defaultAxesFontSize", 10);
set(groot, "defaultTextFontSize", 10);
set(groot, "defaultLegendFontSize", 9);
set(groot, "defaultAxesLineWidth", 0.8);

%% ===============================================================
% TRUE FORWARD MODEL PARAMETERS
% ===============================================================

paramsTrue.g = 13;
paramsTrue.carrier = 8;
paramsTrue.curvature = 0.35;
paramsTrue.cross = 0.12;

trueCurvature = paramsTrue.curvature;

%% ===============================================================
% TRUE DISPLACEMENT PROFILE
% ===============================================================

Nmode = 3;

basis = zeros(Nmode, numel(z));

for n = 1:Nmode

    basis(n,:) = sin(n*pi*(z + 0.5));

end

uTrue = ...
    (reshape(aTruePm,1,[]) * 1e-3) * basis;

truePeakPm = max(abs(uTrue)) * 1e3;

%% ===============================================================
% GENERATE SYNTHETIC DATA USING TRUE INSTRUMENTAL PHASE
% ===============================================================

preTrue = make_precomputed_model(z, s, x, paramsTrue);

psiClean = forward_precomputed(uTrue, preTrue);

psiScale = sqrt(mean(abs(psiClean(:)).^2));

if noiseFraction > 0

    noiseSigma = noiseFraction * psiScale;

    complexNoise = ...
        noiseSigma/sqrt(2) .* ...
        (randn(size(psiClean)) + 1i*randn(size(psiClean)));

    psiObs = psiClean + complexNoise;

else

    psiObs = psiClean;

end

%% ===============================================================
% STORAGE
% ===============================================================

nCurvature = numel(curvatureRatios);

aFitAllPm = nan(Nmode, nCurvature);

peakRecoveredPm = nan(1, nCurvature);
peakBiasPercent = nan(1, nCurvature);
profileRMSEPm = nan(1, nCurvature);
complexRelativeError = nan(1, nCurvature);

%% ===============================================================
% OPTIMIZATION SETTINGS
% ===============================================================

lowerPm = -coefficientBoundPm * ones(Nmode,1);
upperPm =  coefficientBoundPm * ones(Nmode,1);

a0Pm = zeros(Nmode,1);

options = optimoptions( ...
    "lsqnonlin", ...
    "Display", "off", ...
    "FunctionTolerance", 1e-11, ...
    "StepTolerance", 1e-11, ...
    "OptimalityTolerance", 1e-11, ...
    "MaxIterations", 400, ...
    "MaxFunctionEvaluations", 4000);

%% ===============================================================
% INSTRUMENTAL PHASE SENSITIVITY STUDY
% ===============================================================

for iC = 1:nCurvature

    ratio = curvatureRatios(iC);

    paramsAssumed = paramsTrue;

    paramsAssumed.curvature = ...
        ratio * trueCurvature;

    preAssumed = make_precomputed_model( ...
        z, ...
        s, ...
        x, ...
        paramsAssumed);

    residualFun = @(aPm) residual_complex_pm( ...
        aPm, ...
        psiObs, ...
        preAssumed, ...
        basis);

    aFitPm = lsqnonlin( ...
        residualFun, ...
        a0Pm, ...
        lowerPm, ...
        upperPm, ...
        options);

    aFitAllPm(:,iC) = aFitPm;

    uFit = ...
        (reshape(aFitPm,1,[]) * 1e-3) * basis;

    peakRecoveredPm(iC) = ...
        max(abs(uFit)) * 1e3;

    peakBiasPercent(iC) = ...
        100 * ...
        (peakRecoveredPm(iC) - truePeakPm) / ...
        truePeakPm;

    profileErrorPm = ...
        (uFit - uTrue) * 1e3;

    profileRMSEPm(iC) = ...
        sqrt(mean(profileErrorPm.^2));

    psiFit = ...
        forward_precomputed(uFit, preAssumed);

    complexRelativeError(iC) = ...
        norm(psiFit(:) - psiObs(:)) / ...
        max(norm(psiObs(:)), eps);

    if showProgress

        fprintf( ...
            "a_xx/a_xx0 = %.3f, peak = %.5f pm, bias = %.5f percent, profile RMSE = %.5f pm, relative complex error = %.6e\n", ...
            ratio, ...
            peakRecoveredPm(iC), ...
            peakBiasPercent(iC), ...
            profileRMSEPm(iC), ...
            complexRelativeError(iC));

    end

end

%% ===============================================================
% PRINT COEFFICIENT SUMMARY
% ===============================================================

fprintf("\n");
fprintf("TRUE COEFFICIENTS, pm\n");
fprintf("---------------------\n");
disp(aTruePm);

fprintf("RECOVERED COEFFICIENTS, pm\n");
fprintf("--------------------------\n");

for iC = 1:nCurvature

    fprintf("a_xx/a_xx0 = %.3f : ",curvatureRatios(iC));
    fprintf("%10.5f ",aFitAllPm(:,iC));
    fprintf("\n");

end

%% ===============================================================
% FIGURE 3c
% SINGLE PANEL FOR COMBINED FIGURE 3
% ===============================================================

fig = figure( ...
    "Color","w", ...
    "Units","centimeters", ...
    "Position",[2 2 8.5 7.0]);

ax = axes(fig);

plot( ...
    ax, ...
    curvatureRatios, ...
    peakRecoveredPm, ...
    "-o", ...
    "LineWidth",1.8, ...
    "MarkerSize",6, ...
    "MarkerFaceColor","w");

hold(ax,"on");

yline( ...
    ax, ...
    truePeakPm, ...
    "--", ...
    "True peak", ...
    "Interpreter","latex", ...
    "LabelHorizontalAlignment","left", ...
    "HandleVisibility","off");

xline( ...
    ax, ...
    1.0, ...
    ":", ...
    "$a_{xx}/a_{xx,0}=1$", ...
    "Interpreter","latex", ...
    "LabelVerticalAlignment","bottom", ...
    "HandleVisibility","off");

grid(ax,"on");
box(ax,"on");

ax.GridAlpha = 0.15;
ax.TickDir = "out";
ax.TickLength = [0.015 0.015];

xlabel( ...
    ax, ...
    "Assumed curvature ratio, $a_{xx}/a_{xx,0}$");

ylabel( ...
    ax, ...
    "Recovered peak displacement (pm)");

%% ===============================================================
% EXPORT
% ===============================================================

if fastMode

    exportName = ...
        "fig3c_instrumental_phase_sensitivity_FAST.png";

else

    exportName = ...
        "fig3c_instrumental_phase_sensitivity.png";

end

exportgraphics( ...
    fig, ...
    exportName, ...
    "Resolution",600);

%% ===============================================================
% SAVE NUMERICAL RESULTS
% ===============================================================

CurvatureRatio = reshape(curvatureRatios,[],1);
AssumedCurvature = CurvatureRatio * trueCurvature;

RecoveredPeakPm = reshape(peakRecoveredPm,[],1);
PeakBiasPercent = reshape(peakBiasPercent,[],1);
ProfileRMSEPm = reshape(profileRMSEPm,[],1);
ComplexRelativeError = reshape(complexRelativeError,[],1);

A1RecoveredPm = reshape(aFitAllPm(1,:),[],1);
A2RecoveredPm = reshape(aFitAllPm(2,:),[],1);
A3RecoveredPm = reshape(aFitAllPm(3,:),[],1);

resultsTable = table( ...
    CurvatureRatio, ...
    AssumedCurvature, ...
    RecoveredPeakPm, ...
    PeakBiasPercent, ...
    ProfileRMSEPm, ...
    ComplexRelativeError, ...
    A1RecoveredPm, ...
    A2RecoveredPm, ...
    A3RecoveredPm);

writetable( ...
    resultsTable, ...
    "fig3c_instrumental_phase_sensitivity_results.csv");

%% ===============================================================
% LOCAL FUNCTIONS
% ===============================================================

function pre = make_precomputed_model(z,s,x,params)

    z = reshape(z,1,[]);
    s = reshape(s,[],1);
    x = reshape(x,1,[]);

    pre.z = z;

    pre.lowerMask = z <= 0;
    pre.upperMask = z >= 0;

    pre.kernel = exp(-2*pi*1i*(s*z));

    [X,S] = meshgrid(x,s);

    chi = ...
        2*pi*params.carrier*X ...
        + params.curvature*X.^2 ...
        + params.cross*X.*S;

    pre.phaseMinus = exp(1i*chi/2);
    pre.phasePlus = exp(-1i*chi/2);

    pre.g = params.g;

end

function psi = forward_precomputed(u,pre)

    phaseDepth = exp(-2*pi*1i*pre.g*u);

    lowerIntegrand = ...
        pre.kernel(:,pre.lowerMask) .* ...
        phaseDepth(pre.lowerMask);

    upperIntegrand = ...
        pre.kernel(:,pre.upperMask) .* ...
        phaseDepth(pre.upperMask);

    Aminus = trapz( ...
        pre.z(pre.lowerMask), ...
        lowerIntegrand, ...
        2);

    Aplus = trapz( ...
        pre.z(pre.upperMask), ...
        upperIntegrand, ...
        2);

    psi = ...
        Aminus .* pre.phaseMinus ...
        + Aplus .* pre.phasePlus;

end

function r = residual_complex_pm(aPm,psiObs,preAssumed,basis)

    uModel = ...
        (reshape(aPm,1,[]) * 1e-3) * basis;

    psiModel = ...
        forward_precomputed(uModel,preAssumed);

    difference = psiModel - psiObs;

    r = [
        real(difference(:))
        imag(difference(:))
        ];

end
