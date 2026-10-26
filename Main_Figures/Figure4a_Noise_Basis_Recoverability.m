% FIGURE 3a
% Noise sensitivity and individual depth mode recoverability
%
% Corrected publication study:
% Each synthetic test contains one known depth mode with the same
% amplitude. The inverse fit recovers the coefficient of that mode.
% This directly tests whether finer depth modes become less recoverable.
%
% Fast mode is intended for development.
% Set fastMode = false once for the final publication calculation.

clear;
close all;
clc;

rng(8);

%% ===============================================================
% USER CONTROLS
% ===============================================================

fastMode = true;

noiseLevels = [0 0.005 0.01 0.02 0.05];

modeNumbers = 1:20;

trueAmplitudePm = 10;

coefficientLowerPm = 0;
coefficientUpperPm = 20;

useLogY = true;

if fastMode

    % Rapid development calculation
    nRealizations = 5;

    z = linspace(-0.5, 0.5, 201);
    s = linspace(-8, 8, 121);
    x = linspace(-1, 1, 81);

else

    % Final publication calculation
    nRealizations = 20;

    z = linspace(-0.5, 0.5, 401);
    s = linspace(-8, 8, 241);
    x = linspace(-1, 1, 241);

end

showProgress = true;

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
% PRECOMPUTE TERMS THAT DO NOT CHANGE
% ===============================================================

pre = make_precomputed_model(z, s, x, params);

%% ===============================================================
% STORAGE
% ===============================================================

nNoise = numel(noiseLevels);
nMode = numel(modeNumbers);

coefficientErrorPercent = nan(nNoise, nMode, nRealizations);

meanErrorPercent = nan(nNoise, nMode);
stdErrorPercent = nan(nNoise, nMode);

%% ===============================================================
% SCALAR OPTIMIZATION SETTINGS
% ===============================================================

optimOptions = optimset( ...
    "Display", "off", ...
    "TolX", 1e-7, ...
    "MaxIter", 100);

%% ===============================================================
% MAIN RECOVERABILITY STUDY
% ===============================================================

for iMode = 1:nMode

    modeNumber = modeNumbers(iMode);

    % Surface relaxed sine mode
    phi = sin(modeNumber*pi*(z + 0.5));

    % Known synthetic truth
    uTrue = trueAmplitudePm*1e-3*phi;

    psiClean = forward_precomputed(uTrue, pre);

    psiScale = sqrt(mean(abs(psiClean(:)).^2));

    for iNoise = 1:nNoise

        noiseFraction = noiseLevels(iNoise);

        if noiseFraction == 0

            realizationsToRun = 1;

        else

            realizationsToRun = nRealizations;

        end

        temporaryErrors = nan(1, realizationsToRun);

        for iReal = 1:realizationsToRun

            if noiseFraction == 0

                psiObs = psiClean;

            else

                noiseSigma = noiseFraction*psiScale;

                complexNoise = ...
                    noiseSigma/sqrt(2) .* ...
                    (randn(size(psiClean)) + 1i*randn(size(psiClean)));

                psiObs = psiClean + complexNoise;

            end

            costFun = @(aPm) scalar_complex_cost( ...
                aPm, ...
                phi, ...
                psiObs, ...
                pre);

            aFitPm = fminbnd( ...
                costFun, ...
                coefficientLowerPm, ...
                coefficientUpperPm, ...
                optimOptions);

            temporaryErrors(iReal) = ...
                100*abs(aFitPm - trueAmplitudePm)/trueAmplitudePm;

            if showProgress

                fprintf( ...
                    "Mode %02d, noise %4.1f percent, run %02d, fitted %.5f pm, error %.5f percent\n", ...
                    modeNumber, ...
                    100*noiseFraction, ...
                    iReal, ...
                    aFitPm, ...
                    temporaryErrors(iReal));

            end

        end

        if noiseFraction == 0

            coefficientErrorPercent(iNoise,iMode,:) = temporaryErrors(1);

        else

            coefficientErrorPercent(iNoise,iMode,1:realizationsToRun) = ...
                temporaryErrors;

        end

        meanErrorPercent(iNoise,iMode) = mean(temporaryErrors);

        if realizationsToRun > 1

            stdErrorPercent(iNoise,iMode) = std(temporaryErrors);

        else

            stdErrorPercent(iNoise,iMode) = 0;

        end

    end

end

%% ===============================================================
% PRINT SUMMARY
% ===============================================================

fprintf("\n");
fprintf("MEAN RELATIVE COEFFICIENT ERROR, percent\n");
fprintf("----------------------------------------\n");

for iNoise = 1:nNoise

    fprintf("Noise %4.1f percent : ",100*noiseLevels(iNoise));
    fprintf("%8.4f ",meanErrorPercent(iNoise,:));
    fprintf("\n");

end

fprintf("\n");
fprintf("STANDARD DEVIATION, percent\n");
fprintf("---------------------------\n");

for iNoise = 1:nNoise

    fprintf("Noise %4.1f percent : ",100*noiseLevels(iNoise));
    fprintf("%8.4f ",stdErrorPercent(iNoise,:));
    fprintf("\n");

end

%% ===============================================================
% FIGURE 3a
% PUBLICATION STYLE WITH ERROR BARS
% ===============================================================

fig = figure( ...
    "Color","w", ...
    "Units","centimeters", ...
    "Position",[2 2 11.5 8.0]);

ax = axes(fig);

hold(ax,"on");
box(ax,"on");

colorOrder = get(groot,"defaultAxesColorOrder");
nColors = size(colorOrder,1);

lineHandles = gobjects(nNoise,1);

%% Plotting floor for logarithmic axis

allPositive = meanErrorPercent(meanErrorPercent > 0);

if isempty(allPositive)

    plotFloor = 1e-6;

else

    plotFloor = max(min(allPositive)/10,1e-6);

end

%% Plot mean plus or minus one standard deviation

markerList = ["o","s","^","d","v"];

for iNoise = 1:nNoise

    thisColor = colorOrder(mod(iNoise-1,nColors)+1,:);

    yMean = meanErrorPercent(iNoise,:);
    yStd = stdErrorPercent(iNoise,:);

    if useLogY

        yPlot = max(yMean,plotFloor);

        lowerErr = min(yStd,0.95*yPlot);

        upperErr = yStd;

    else

        yPlot = yMean;

        lowerErr = min(yStd,yPlot);

        upperErr = yStd;

    end

    lineHandles(iNoise) = errorbar( ...
        ax, ...
        modeNumbers, ...
        yPlot, ...
        lowerErr, ...
        upperErr, ...
        "LineWidth",1.6, ...
        "Marker",markerList(iNoise), ...
        "MarkerSize",5.5, ...
        "MarkerFaceColor","w", ...
        "Color",thisColor, ...
        "CapSize",4);

end

%% ===============================================================
% AXIS FORMATTING
% ===============================================================

xlim(ax,[min(modeNumbers)-0.4 max(modeNumbers)+0.4]);

xticks(ax,modeNumbers);

xlabel(ax,"Depth mode number, $n$");

ylabel(ax,"Mean coefficient recovery error (\%)");

if useLogY

    set(ax,"YScale","log");

end

grid(ax,"on");

ax.GridAlpha = 0.16;
ax.MinorGridAlpha = 0.08;

if useLogY

    ax.YMinorGrid = "on";

end

ax.TickDir = "out";
ax.TickLength = [0.015 0.015];

%% ===============================================================
% LEGEND
% ===============================================================

legendText = strings(1,nNoise);

for iNoise = 1:nNoise

    legendText(iNoise) = ...
        sprintf("%.1f\\%%",100*noiseLevels(iNoise));

end

lgd = legend( ...
    ax, ...
    lineHandles, ...
    legendText, ...
    "Location","northwest", ...
    "Interpreter","latex", ...
    "Box","off");

title(lgd,"Complex-wave noise");

%% ===============================================================
% SMALL FIGURE NOTE
% ===============================================================

text( ...
    ax, ...
    0.98, ...
    0.04, ...
    "Error bars: $\pm 1$ SD", ...
    "Units","normalized", ...
    "HorizontalAlignment","right", ...
    "VerticalAlignment","bottom", ...
    "Interpreter","latex", ...
    "FontSize",9);

%% ===============================================================
% EXPORT
% ===============================================================

if fastMode

    exportName = ...
        "fig3a_noise_mode_recoverability_FAST.png";

else

    exportName = ...
        "fig3a_noise_mode_recoverability.png";

end

exportgraphics( ...
    fig, ...
    exportName, ...
    "Resolution",600);

%% ===============================================================
% SAVE RESULTS WITHOUT TABLE WARNINGS
% ===============================================================

nRows = nNoise*nMode;

NoisePercent = zeros(nRows,1);
ModeNumber = zeros(nRows,1);
MeanErrorPercent = zeros(nRows,1);
StdErrorPercent = zeros(nRows,1);

row = 0;

for iNoise = 1:nNoise

    for iMode = 1:nMode

        row = row + 1;

        NoisePercent(row) = 100*noiseLevels(iNoise);
        ModeNumber(row) = modeNumbers(iMode);
        MeanErrorPercent(row) = meanErrorPercent(iNoise,iMode);
        StdErrorPercent(row) = stdErrorPercent(iNoise,iMode);

    end

end

resultsTable = table( ...
    NoisePercent, ...
    ModeNumber, ...
    MeanErrorPercent, ...
    StdErrorPercent);

writetable( ...
    resultsTable, ...
    "fig3a_noise_mode_recoverability_results.csv");

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

function cost = scalar_complex_cost(aPm,phi,psiObs,pre)

    uModel = aPm*1e-3*phi;

    psiModel = forward_precomputed(uModel,pre);

    difference = psiModel - psiObs;

    cost = mean(abs(difference(:)).^2);

end
