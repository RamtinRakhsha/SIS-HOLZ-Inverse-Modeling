% FIGURE 3b
% Sensitivity of SIS-HOLZ inverse recovery to specimen thickness
%
% Synthetic data are generated at the true reference thickness t0.
% During inversion, the assumed thickness is varied while the three
% displacement coefficients are refitted.
%
% Outputs:
%   fig3b_thickness_sensitivity.png
%   fig3b_thickness_sensitivity_results.csv
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

thicknessRatios = [0.95 0.975 1.00 1.025 1.05];

% Reference thickness is normalized to 1.
% Only the ratio t/t0 matters in this synthetic study.
t0 = 1.0;

% Known synthetic coefficients in pm
aTruePm = [32; -4; 2];

% Bounds for the fitted coefficients in pm
coefficientBoundPm = 60;

% Optional synthetic complex wave noise
% Keep at zero to isolate thickness sensitivity.
noiseFraction = 0.0;

showProgress = true;

if fastMode

    eta = linspace(-0.5, 0.5, 201);
    s = linspace(-8, 8, 121);
    x = linspace(-1, 1, 81);

else

    eta = linspace(-0.5, 0.5, 401);
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
% FIXED FORWARD MODEL PARAMETERS
% ===============================================================

params.g = 13;
params.carrier = 8;
params.curvature = 0.35;
params.cross = 0.12;

%% ===============================================================
% TRUE DISPLACEMENT PROFILE
%
% The basis is defined on normalized depth eta = z/t.
% This keeps the profile shape tied to relative depth while physical
% thickness enters the diffraction integral through z = t*eta.
% ===============================================================

Nmode = 3;

basis = zeros(Nmode, numel(eta));

for n = 1:Nmode
    basis(n,:) = sin(n*pi*(eta + 0.5));
end

uTrue = ...
    (reshape(aTruePm,1,[]) * 1e-3) * basis;

truePeakPm = max(abs(uTrue)) * 1e3;

%% ===============================================================
% GENERATE SYNTHETIC DATA AT TRUE THICKNESS
% ===============================================================

zTrue = t0 * eta;

preTrue = make_precomputed_model(zTrue, s, x, params);

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

nThickness = numel(thicknessRatios);

aFitAllPm = nan(Nmode, nThickness);

peakRecoveredPm = nan(1, nThickness);
profileRMSEPm = nan(1, nThickness);
complexRelativeError = nan(1, nThickness);

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
% THICKNESS SENSITIVITY STUDY
% ===============================================================

for iT = 1:nThickness

    ratio = thicknessRatios(iT);

    tAssumed = ratio * t0;

    zAssumed = tAssumed * eta;

    preAssumed = make_precomputed_model( ...
        zAssumed, ...
        s, ...
        x, ...
        params);

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

    aFitAllPm(:,iT) = aFitPm;

    uFit = ...
        (reshape(aFitPm,1,[]) * 1e-3) * basis;

    peakRecoveredPm(iT) = ...
        max(abs(uFit)) * 1e3;

    profileErrorPm = ...
        (uFit - uTrue) * 1e3;

    profileRMSEPm(iT) = ...
        sqrt(mean(profileErrorPm.^2));

    psiFit = ...
        forward_precomputed(uFit, preAssumed);

    complexRelativeError(iT) = ...
        norm(psiFit(:) - psiObs(:)) / ...
        max(norm(psiObs(:)), eps);

    if showProgress

        fprintf( ...
            "t/t0 = %.3f, peak = %.5f pm, profile RMSE = %.5f pm, relative complex error = %.6e\n", ...
            ratio, ...
            peakRecoveredPm(iT), ...
            profileRMSEPm(iT), ...
            complexRelativeError(iT));

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

for iT = 1:nThickness

    fprintf("t/t0 = %.3f : ",thicknessRatios(iT));
    fprintf("%10.5f ",aFitAllPm(:,iT));
    fprintf("\n");

end

%% ===============================================================
% FIGURE 3b
% PUBLICATION STYLE
% ===============================================================

fig = figure( ...
    "Color","w", ...
    "Units","centimeters", ...
    "Position",[2 2 9.0 11.0]);

tl = tiledlayout( ...
    fig, ...
    2,1, ...
    "TileSpacing","compact", ...
    "Padding","compact");

%% ---------------------------------------------------------------
% Panel a, recovered peak displacement
% ---------------------------------------------------------------

ax1 = nexttile(tl,1);

plot( ...
    ax1, ...
    thicknessRatios, ...
    peakRecoveredPm, ...
    "-o", ...
    "LineWidth",1.8, ...
    "MarkerSize",6, ...
    "MarkerFaceColor","w");

hold(ax1,"on");

yline( ...
    ax1, ...
    truePeakPm, ...
    "--", ...
    "True peak", ...
    "Interpreter","latex", ...
    "LabelHorizontalAlignment","left", ...
    "HandleVisibility","off");

xline( ...
    ax1, ...
    1.0, ...
    ":", ...
    "$t/t_0=1$", ...
    "Interpreter","latex", ...
    "LabelVerticalAlignment","bottom", ...
    "HandleVisibility","off");

grid(ax1,"on");
box(ax1,"on");

ax1.GridAlpha = 0.15;
ax1.TickDir = "out";

xlabel(ax1,"Assumed thickness ratio, $t/t_0$");
ylabel(ax1,"Recovered peak displacement (pm)");

title(ax1,"(a) Recovered displacement amplitude");

%% ---------------------------------------------------------------
% Panel b, displacement profile RMSE
% ---------------------------------------------------------------

ax2 = nexttile(tl,2);

plot( ...
    ax2, ...
    thicknessRatios, ...
    profileRMSEPm, ...
    "-o", ...
    "LineWidth",1.8, ...
    "MarkerSize",6, ...
    "MarkerFaceColor","w");

hold(ax2,"on");

xline( ...
    ax2, ...
    1.0, ...
    ":", ...
    "$t/t_0=1$", ...
    "Interpreter","latex", ...
    "LabelVerticalAlignment","bottom", ...
    "HandleVisibility","off");

grid(ax2,"on");
box(ax2,"on");

ax2.GridAlpha = 0.15;
ax2.TickDir = "out";

xlabel(ax2,"Assumed thickness ratio, $t/t_0$");
ylabel(ax2,"Profile RMSE (pm)");

title(ax2,"(b) Displacement profile error");

%% ===============================================================
% EXPORT
% ===============================================================

if fastMode

    exportName = ...
        "fig3b_thickness_sensitivity_FAST.png";

else

    exportName = ...
        "fig3b_thickness_sensitivity.png";

end

exportgraphics( ...
    fig, ...
    exportName, ...
    "Resolution",600);

%% ===============================================================
% SAVE NUMERICAL RESULTS
% ===============================================================

ThicknessRatio = reshape(thicknessRatios,[],1);
RecoveredPeakPm = reshape(peakRecoveredPm,[],1);
ProfileRMSEPm = reshape(profileRMSEPm,[],1);
ComplexRelativeError = reshape(complexRelativeError,[],1);

A1RecoveredPm = reshape(aFitAllPm(1,:),[],1);
A2RecoveredPm = reshape(aFitAllPm(2,:),[],1);
A3RecoveredPm = reshape(aFitAllPm(3,:),[],1);

resultsTable = table( ...
    ThicknessRatio, ...
    RecoveredPeakPm, ...
    ProfileRMSEPm, ...
    ComplexRelativeError, ...
    A1RecoveredPm, ...
    A2RecoveredPm, ...
    A3RecoveredPm);

writetable( ...
    resultsTable, ...
    "fig3b_thickness_sensitivity_results.csv");

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

    difference = ...
        psiModel - psiObs;

    r = [
        real(difference(:))
        imag(difference(:))
        ];

end
