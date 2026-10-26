%% Figure 1 - Generalized SIS-HOLZ forward simulations
% Standalone MATLAB script for the manuscript:
% "Generalized Inverse Modeling of SIS-HOLZ Strain Profiles"
%
% Top row    : prescribed projected displacement profiles
% Middle row : normalized SIS-HOLZ interferograms
% Bottom row : relative phase traces
%
% Each interferogram is normalized by its own maximum intensity.

clear;
close all;
clc;

%% Plot defaults
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultTextFontName','Times New Roman');
set(groot,'defaultLegendFontName','Times New Roman');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Coordinates
z = linspace(-0.5,0.5,401);
s = linspace(-8,8,241);
x = linspace(-1,1,241);

%% Forward-model parameters
params.g = 13;
params.carrier = 8;
params.curvature = 0.35;
params.cross = 0.12;

%% Prescribed displacement profiles (nm)
profiles = {
    0.035*cos(pi*z)
    0.035*(1-2*abs(z))
    0.035*cos(pi*(z-0.10)).*(1+0.35*z)
    0.040*(1-(2*z).^2).*(1+0.28*z)
    };

names = {
    'Symmetric bell'
    'V-shaped'
    'Asymmetric bell'
    'Surface-relaxed asymmetric'
    };

nProfile = numel(profiles);

Iall = cell(nProfile,1);
phaseAll = cell(nProfile,1);

for k = 1:nProfile
    [I,psi] = sis_holz_forward(z,s,x,profiles{k},params);
    Iall{k} = I./max(I(:));

    [~,ix0] = min(abs(x));
    ph = unwrap(angle(psi(:,ix0)));
    ph = ph - ph(round((numel(ph)+1)/2));
    phaseAll{k} = ph;
end

%% Figure
fig = figure('Color','w','Units','centimeters','Position',[1 1 23 16]);
tl = tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');

% Row 1: displacement
for k = 1:nProfile
    ax = nexttile(tl,k);
    plot(ax,z,1000*profiles{k},'LineWidth',1.6);
    grid(ax,'on');
    box(ax,'on');
    xlim(ax,[-0.5 0.5]);
    xticks(ax,[-0.5 0 0.5]);
    title(ax,names{k},'FontWeight','normal');
    if k == 1
        ylabel(ax,'Projected displacement, $u$ (pm)');
    end
    xlabel(ax,'Normalized depth, $z/t$');
end

% Row 2: interferograms
axMiddle = gobjects(1,nProfile);
for k = 1:nProfile
    ax = nexttile(tl,4+k);
    axMiddle(k) = ax;
    imagesc(ax,x,s,Iall{k});
    axis(ax,'xy');
    colormap(ax,hot);
    caxis(ax,[0 1]);
    if k == 1
        ylabel(ax,'Diffraction coordinate, $s$');
    end
    xlabel(ax,'Image coordinate, $x$');
end

% Row 3: relative phase
for k = 1:nProfile
    ax = nexttile(tl,8+k);
    plot(ax,s,phaseAll{k},'LineWidth',1.6);
    grid(ax,'on');
    box(ax,'on');
    xlim(ax,[min(s) max(s)]);
    if k == 1
        ylabel(ax,'Relative phase (rad)');
    end
    xlabel(ax,'Diffraction coordinate, $s$');
end

% Colorbar attached to the rightmost middle-row axes.
% This syntax is compatible with older MATLAB releases that do not support
% colorbar(tiledlayout,...).
cb = colorbar(axMiddle(end));
cb.Location = 'eastoutside';
cb.Label.String = 'Normalized intensity, $I/I_{\max}$';
cb.Label.Interpreter = 'latex';

exportgraphics(fig,'fig1_forward_profiles.png','Resolution',600);

%% Local function
function [I,psi] = sis_holz_forward(z,s,x,u,params)

    z = z(:).';
    s = s(:);
    x = x(:).';

    phaseDepth = exp(-2*pi*1i*params.g*u);
    kernel = exp(-2*pi*1i*(s*z));

    lowerMask = z <= 0;
    upperMask = z >= 0;

    lowerIntegrand = kernel(:,lowerMask).*phaseDepth(lowerMask);
    upperIntegrand = kernel(:,upperMask).*phaseDepth(upperMask);

    Aminus = trapz(z(lowerMask),lowerIntegrand,2);
    Aplus  = trapz(z(upperMask),upperIntegrand,2);

    [X,S] = meshgrid(x,s);

    chi = 2*pi*params.carrier*X ...
        + params.curvature*X.^2 ...
        + params.cross*X.*S;

    psi = Aminus.*exp(1i*chi/2) ...
        + Aplus.*exp(-1i*chi/2);

    I = abs(psi).^2;
end
