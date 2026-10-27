% FIGURE 3
% Publication synthetic validation for generalized SIS-HOLZ inversion
%
% Standalone script. Local functions at the end implement the forward
% SIS-HOLZ model and phase extraction used in this validation.
%
% The test uses:
%   - a smooth asymmetric surface-relaxed displacement profile,
%   - 3 sine depth modes,
%   - 1% complex-wave RMS noise,
%   - fixed optical parameters,
%   - lsqnonlin fitting of the displacement coefficients only.
%
% Requires Optimization Toolbox for lsqnonlin.

clear;
close all;
clc;

rng(8);

%% USER CONTROLS

noiseFraction = 0.01;   % 0.01 = 1 percent complex-wave RMS noise
Nmode = 3;
nStarts = 8;

colorbarGap = 0.004;
colorbarWidth = 0.010;

%% Global LaTeX formatting

set(groot, "defaultTextInterpreter", "latex");
set(groot, "defaultAxesTickLabelInterpreter", "latex");
set(groot, "defaultLegendInterpreter", "latex");

set(groot, "defaultAxesFontSize", 10);
set(groot, "defaultTextFontSize", 10);
set(groot, "defaultLegendFontSize", 9);
set(groot, "defaultAxesLineWidth", 0.8);

%% Coordinates

z = linspace(-0.5, 0.5, 401);
s = linspace(-8, 8, 241);
x = linspace(-1, 1, 241);

%% Fixed forward-model parameters

params.g = 13;
params.carrier = 8;
params.curvature = 0.35;
params.cross = 0.12;

%% Surface-relaxed basis
% Every sine mode is zero at z/t = -0.5 and +0.5.

basis = zeros(Nmode, numel(z));

for n = 1:Nmode
    basis(n,:) = sin(n*pi*(z + 0.5));
end

%% Known synthetic truth, pm

aTruePm = [32; -4; 2];

uTrue = (reshape(aTruePm,1,[]) * 1e-3) * basis;

%% Generate clean synthetic data

[Iclean, psiClean] = sis_holz_forward(z, s, x, uTrue, params);

%% Add known complex Gaussian noise

psiScale = sqrt(mean(abs(psiClean(:)).^2));
noiseSigma = noiseFraction * psiScale;

complexNoise = ...
    noiseSigma/sqrt(2) .* ...
    (randn(size(psiClean)) + 1i*randn(size(psiClean)));

psiObs = psiClean + complexNoise;
Iobs = abs(psiObs).^2;

%% Check exact forward-model consistency before fitting

rConsistency = residual_complex_pm( ...
    aTruePm, z, s, x, psiClean, params, basis);

consistencyNorm = norm(rConsistency);

fprintf("\nMODEL CONSISTENCY TEST\n");
fprintf("----------------------\n");
fprintf("Residual norm at true coefficients = %.6e\n\n", consistencyNorm);

if consistencyNorm > 1e-8
    error("Forward generation and inverse fitting models are inconsistent.");
end

%% Fit coefficients only

lowerPm = -60 * ones(Nmode,1);
upperPm =  60 * ones(Nmode,1);

residualFun = @(aPm) residual_complex_pm( ...
    aPm, z, s, x, psiObs, params, basis);

options = optimoptions( ...
    "lsqnonlin", ...
    "Display", "off", ...
    "FunctionTolerance", 1e-12, ...
    "StepTolerance", 1e-12, ...
    "OptimalityTolerance", 1e-12, ...
    "MaxIterations", 1000, ...
    "MaxFunctionEvaluations", 10000);

% Deterministic multistart
startMatrix = zeros(Nmode,nStarts);
startMatrix(:,1) = zeros(Nmode,1);

if nStarts >= 2
    startMatrix(:,2) = [20;0;0];
end

if nStarts >= 3
    startMatrix(:,3) = [35;0;0];
end

for q = 4:nStarts
    startMatrix(:,q) = -20 + 40*rand(Nmode,1);
end

bestResnorm = inf;
aFitPm = zeros(Nmode,1);
bestExitflag = NaN;

for q = 1:nStarts

    [aCandidate,resnormCandidate,~,exitflagCandidate] = lsqnonlin( ...
        residualFun, ...
        startMatrix(:,q), ...
        lowerPm, ...
        upperPm, ...
        options);

    if resnormCandidate < bestResnorm
        bestResnorm = resnormCandidate;
        aFitPm = aCandidate;
        bestExitflag = exitflagCandidate;
    end

end

%% Reconstruct best fit

uFit = (reshape(aFitPm,1,[]) * 1e-3) * basis;

[Ifit, psiFit] = sis_holz_forward(z, s, x, uFit, params);

%% Metrics

profileRMSEPm = sqrt(mean(((uFit-uTrue)*1e3).^2));

coefficientErrorPm = aFitPm-aTruePm;
maxCoefficientErrorPm = max(abs(coefficientErrorPm));

phaseObs = sis_holz_phase_trace_local(psiObs, x);
phaseFit = sis_holz_phase_trace_local(psiFit, x);

phaseOffset = mean(phaseFit-phaseObs);
phaseFitAligned = phaseFit-phaseOffset;

phaseRMSE = sqrt(mean((phaseFitAligned-phaseObs).^2));

commonImax = max([Iobs(:); Ifit(:)]);

IobsNorm = Iobs/commonImax;
IfitNorm = Ifit/commonImax;

Iresidual = IobsNorm-IfitNorm;
intensityRMSE = sqrt(mean(Iresidual(:).^2));

fprintf("PUBLICATION SYNTHETIC VALIDATION\n");
fprintf("--------------------------------\n");
fprintf("Complex-wave RMS noise = %.2f percent\n\n",100*noiseFraction);

fprintf("True coefficients, pm\n");
disp(aTruePm);

fprintf("Recovered coefficients, pm\n");
disp(aFitPm);

fprintf("Coefficient errors, pm\n");
disp(coefficientErrorPm);

fprintf("Maximum coefficient error = %.6f pm\n",maxCoefficientErrorPm);
fprintf("Profile RMSE = %.6f pm\n",profileRMSEPm);
fprintf("Phase RMSE = %.6e rad\n",phaseRMSE);
fprintf("Normalized intensity RMSE = %.6e\n",intensityRMSE);
fprintf("Best residual norm = %.6e\n",bestResnorm);
fprintf("Best exit flag = %d\n\n",bestExitflag);

%% Publication figure

fig = figure( ...
    "Color","w", ...
    "Units","centimeters", ...
    "Position",[2 2 20.5 12.5]);

tl = tiledlayout( ...
    fig,2,3, ...
    "TileSpacing","compact", ...
    "Padding","compact");

% Reserve white space on the right for the panel-(c) color bar label.
% Without this, the tiled layout expands to the figure edge and MATLAB
% clips the vertical label even when the color bar itself remains visible.
tl.OuterPosition = [0.03 0.04 0.88 0.92];

ax = gobjects(1,6);

%% (a) Displacement recovery

ax(1) = nexttile(tl,1);

plot(ax(1),z,uTrue*1e3,"LineWidth",2);
hold(ax(1),"on");
plot(ax(1),z,uFit*1e3,"--","LineWidth",2);

grid(ax(1),"on");
box(ax(1),"on");

xlabel(ax(1),"Normalized depth, $z/t$");
ylabel(ax(1),"Projected displacement, $u$ (pm)");
title(ax(1),"(a) Displacement recovery", ...
    "Interpreter","latex", ...
    "FontWeight","normal");

xlim(ax(1),[-0.6 0.6]);
xticks(ax(1),[-0.5 0 0.5]);
ylim(ax(1),[0 35]);

lgd1 = legend(ax(1), ...
    {"True profile","Recovered profile"}, ...
    "Location","south");

%% (b) Synthetic data

ax(2) = nexttile(tl,2);

imagesc(ax(2),x,s,IobsNorm);
axis(ax(2),"xy");
clim(ax(2),[0 1]);

xlabel(ax(2),"Image coordinate, $x$");
ylabel(ax(2),"Diffraction coordinate, $s$");
title(ax(2),"(b) Synthetic data, 1\\% noise", ...
    "Interpreter","latex", ...
    "FontWeight","normal");

colormap(ax(2),hot(256));

%% (c) Best-fit model

ax(3) = nexttile(tl,3);

imagesc(ax(3),x,s,IfitNorm);
axis(ax(3),"xy");
clim(ax(3),[0 1]);

xlabel(ax(3),"Image coordinate, $x$");
ylabel(ax(3),"Diffraction coordinate, $s$");
title(ax(3),"(c) Best-fit model", ...
    "Interpreter","latex", ...
    "FontWeight","normal");

colormap(ax(3),hot(256));

%% (d) Phase agreement

ax(4) = nexttile(tl,4);

plot(ax(4),s,phaseObs,"LineWidth",1.7);
hold(ax(4),"on");
plot(ax(4),s,phaseFitAligned,"--","LineWidth",2);

grid(ax(4),"on");
box(ax(4),"on");

xlabel(ax(4),"Diffraction coordinate, $s$");
ylabel(ax(4),"Relative phase, $\phi$ (rad)");
title(ax(4),"(d) Phase agreement", ...
    "Interpreter","latex", ...
    "FontWeight","normal");

ylim(ax(4),[0 35]);

lgd4 = legend(ax(4), ...
    {"Noisy synthetic phase","Recovered model"}, ...
    "Location","north");

%% (e) Recovered coefficients

ax(5) = nexttile(tl,5);

bar(ax(5),1:Nmode,[aTruePm,aFitPm]);

grid(ax(5),"on");
box(ax(5),"on");

xlabel(ax(5),"Depth mode, $n$");
ylabel(ax(5),"Coefficient, $a_n$ (pm)");
title(ax(5),"(e) Recovered coefficients", ...
    "Interpreter","latex", ...
    "FontWeight","normal");

legend(ax(5),{"True","Recovered"},"Location","northeast");
xticks(ax(5),1:Nmode);

%% (f) Residual

ax(6) = nexttile(tl,6);

residualScale = 3*std(Iresidual(:));

if residualScale <= 0
    residualScale = 1e-12;
end

imagesc(ax(6),x,s,Iresidual);
axis(ax(6),"xy");
clim(ax(6),[-residualScale residualScale]);

xlabel(ax(6),"Image coordinate, $x$");
ylabel(ax(6),"Diffraction coordinate, $s$");
title(ax(6),"(f) Residual, data minus fit", ...
    "Interpreter","latex", ...
    "FontWeight","normal");

text(ax(6),0.04,0.07, ...
    sprintf("$\\mathrm{RMSE}_I=%.3g$",intensityRMSE), ...
    "Units","normalized", ...
    "Interpreter","latex", ...
    "Color","w");

%% Shared intensity color bar, aligned with panel c only

drawnow;

for j = 1:6
    ax(j).Units = "normalized";
end

savedPos = cell(1,6);

for j = 1:6
    savedPos{j} = ax(j).Position;
end

% Reserve a small amount of space inside panel (c) for the color bar and
% its vertical label. This prevents the label from being clipped at the
% right edge of the exported figure.
p = savedPos{3};
rightReserve = 0.055;

ax(3).Position = [ ...
    p(1), ...
    p(2), ...
    p(3)-rightReserve, ...
    p(4)];

p3 = ax(3).Position;

cb = colorbar(ax(3));
cb.Units = "normalized";

cb.Position = [ ...
    p3(1)+p3(3)+colorbarGap+0.003, ...
    p3(2), ...
    colorbarWidth, ...
    p3(4)];

cb.TickLabelInterpreter = "latex";
cb.FontSize = 9;
cb.Label.String = "Normalized intensity, $I/I_{\max}$";
cb.Label.Interpreter = "latex";
cb.Label.FontSize = 10;
cb.Label.Rotation = 270;
cb.Label.VerticalAlignment = "bottom";

% Restore all other panel positions exactly.
for j = [1 2 4 5 6]
    ax(j).Position = savedPos{j};
end

%% Export

exportgraphics( ...
    fig, ...
    "fig3_publication_synthetic_validation.png", ...
    "Resolution",600);

%% ===============================================================
% Local functions
% ===============================================================

function [I,psi] = sis_holz_forward(z,s,x,u,params)

    z = z(:).';
    s = s(:);
    x = x(:).';

    phaseDepth = exp(-2*pi*1i*params.g*u);

    kernel = exp(-2*pi*1i*(s*z));

    lowerMask = z <= 0;
    upperMask = z >= 0;

    lowerIntegrand = ...
        kernel(:,lowerMask) .* phaseDepth(lowerMask);

    upperIntegrand = ...
        kernel(:,upperMask) .* phaseDepth(upperMask);

    Aminus = trapz( ...
        z(lowerMask), ...
        lowerIntegrand, ...
        2);

    Aplus = trapz( ...
        z(upperMask), ...
        upperIntegrand, ...
        2);

    [X,S] = meshgrid(x,s);

    chi = ...
        2*pi*params.carrier*X ...
        + params.curvature*X.^2 ...
        + params.cross*X.*S;

    psi = ...
        Aminus .* exp(1i*chi/2) ...
        + Aplus .* exp(-1i*chi/2);

    I = abs(psi).^2;

end

function phaseTrace = sis_holz_phase_trace_local(psi,x)

    [~,ix0] = min(abs(x));

    phaseTrace = unwrap(angle(psi(:,ix0)));

end

function r = residual_complex_pm(aPm,z,s,x,psiObs,params,basis)

    uModel = ...
        (reshape(aPm,1,[]) * 1e-3) * basis;

    [~,psiModel] = ...
        sis_holz_forward(z,s,x,uModel,params);

    difference = psiModel-psiObs;

    r = [
        real(difference(:))
        imag(difference(:))
        ];

end
