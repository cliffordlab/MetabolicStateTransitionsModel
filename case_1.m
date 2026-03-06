%% case_1.m
% Lambda–omega system with animated trajectory and time-varying stability.
% 2x3 tiled layout with spans:
%   - ONE big LEFT panel spanning tiles (1,2,4,5)   => nexttile(1,[2 2])
%   - ONE big RIGHT panel spanning tiles (3,6)      => nexttile(3,[2 1])
%
% LEFT  : Phase plane (x–y) with nullclines, vector field, trajectory, peaks, initial point
% RIGHT : Time-series of peak values of x(t) (gray, left y-axis)
%         and cumulative sum of absolute differences between consecutive peaks (black, right y-axis)

clear; close all; clc;

%% ==============================
%  DISPLAY & STYLE PARAMETERS
%  ==============================
fontSizeAxes   = 38;
fontSizeLabel  = 38;
fontSizeTitle  = 38;

% Line / marker styles
lineWidthMain       = 5;    % for nullclines and main curves
lineWidthTrajectory = 3.5;  % trajectory line width
lineWidthPeaks      = 4;    % connection lines between peaks
markerSizePeaks     = 10;

% Colors
stable_color   = [0.65 0.82 1.00];   % stronger blue
unstable_color = [1.00 0.70 0.65];   % stronger red

xnull_color    = [0.80 0.00 0.00];   % dark red   (x-nullcline)
ynull_color    = [0.05 0.55 0.05];   % dark green (y-nullcline)
traj_color     = [0.05 0.10 0.55];   % dark blue  (trajectory)

% Phase-plane domain
xy_lim      = 2.5;        % x,y in [-xy_lim, xy_lim]
nGrid_bg    = 200;      % resolution of background shading
nGrid_field = 20;       % resolution of vector field

% Animation controls (visual speed, not dynamics)
frame_skip  = 40;       % draw every 'frame_skip'-th time sample
pause_time  = 0.01;     % pause between frames (s)

%% ==============================
%  MOVIE SETTINGS
%  ==============================
make_movie      = true;
movie_filename  = 'Movie/case_1_movie.mp4';
movie_fps       = 10;

%% ==============================
%  MODEL PARAMETERS
%  ==============================
lambda_scale = 1.0;
omega_scale  = 1.5;

lambda_base  = @(r) lambda_scale * (1 - r.^2);
omega        = @(r) omega_scale  * 2*pi;

tspan        = [0 30];
x0           = [0.999; 0.0];

% --- ONE-TIME COSINE FLIP WITH DURATION CONTROL ---
flip_start    = 0.0;
flip_duration = 10.0;
flip_end      = flip_start + flip_duration;

lambda_sign = @(t) ...
    (t < flip_start) .* 1 + ...
    (t >= flip_start & t <= flip_end) .* cos(pi * (t - flip_start) / flip_duration) + ...
    (t > flip_end) .* (-1);

%% ==============================
%  NULLCLINE SNAPSHOTS EVERY 10 CYCLES
%  ==============================
osc_period         = 1.0;    % omega ~ 2*pi => period ~1
cycles_per_mark    = 10;
nullcline_interval = osc_period * cycles_per_mark;

%% ==============================
%  SOLVE ODE
%  ==============================
odefun = @(t, x) lambda_omega_rhs_timevarying(t, x, lambda_base, omega, lambda_sign);
opts   = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',0.01);

[t, X]  = ode45(odefun, tspan, x0, opts);
x_traj  = X(:,1);
y_traj  = X(:,2);
N       = numel(t);

%% ==============================
%  PEAK DETECTION (x(t))
%  ==============================
[peakVals, peakTimes] = findpeaks(x_traj, t);
t_peaks = peakTimes(:).';
x_peaks = peakVals(:).';

% Map peak times to indices in the full trajectory
peak_idx_in_traj = zeros(1, numel(t_peaks));
for i = 1:numel(t_peaks)
    [~, peak_idx_in_traj(i)] = min(abs(t - t_peaks(i)));
end
x_peaks_phase = x_traj(peak_idx_in_traj).';
y_peaks_phase = y_traj(peak_idx_in_traj).';

% Precompute max cumulative |diff| value for axis scaling (right axis)
if numel(x_peaks) >= 2
    diffsAll   = abs(diff(x_peaks));
    maxCumDiff = sum(diffsAll);
else
    maxCumDiff = 1;
end

%% ==============================
%  ADAPTATION WINDOW (BASED ON |Δpeak| > THRESHOLD)
%  ==============================
diffThreshold = 0.1;
leverage      = 4;

if numel(x_peaks) >= 2
    diffs = abs(diff(x_peaks));
    idx_above = find(diffs > diffThreshold);

    if ~isempty(idx_above)
        firstPairIdx = idx_above(1);
        lastPairIdx  = idx_above(end);

        idx_start_peak = max(1, firstPairIdx - leverage);
        idx_end_peak   = min(numel(t_peaks), lastPairIdx + 1 + leverage);

        t_adapt_start = t_peaks(idx_start_peak);
        t_adapt_end   = t_peaks(idx_end_peak);
    else
        t_adapt_start = NaN;
        t_adapt_end   = NaN;
    end
else
    t_adapt_start = NaN;
    t_adapt_end   = NaN;
end

%% ==============================
%  FIGURE (2x3 WITH TILE SPANS) + MARGINS (OPTIONAL)
%  ==============================
fig = figure('Color','w','Position',[100 80 1900 980]);

tl  = tiledlayout(fig, 2, 3, 'Padding','compact', 'TileSpacing','compact');

% OPTIONAL: explicit outer margins (uncomment to control)
mLeft=0.01; mRight=0.1; mBottom=0.2; mTop=0.2;
tl.OuterPosition = [mLeft, mBottom, 1-mLeft-mRight, 1-mTop-mBottom];
tl.InnerPosition = tl.OuterPosition;

% Big LEFT panel spanning tiles (1,2,4,5)
axLeft  = nexttile(tl, 1, [2 2]);

% Big RIGHT panel spanning tiles (3,6)
axRight = nexttile(tl, 3, [2 1]);

%% ==============================
%  INIT LEFT PANEL (PHASE PLANE)
%  ==============================
L = init_left_panel(axLeft, xy_lim, nGrid_bg, nGrid_field, ...
    fontSizeAxes, fontSizeLabel, ...
    lineWidthMain, lineWidthTrajectory, lineWidthPeaks, markerSizePeaks, ...
    stable_color, unstable_color, xnull_color, ynull_color, traj_color, ...
    lambda_base, omega, lambda_sign, t(1), x0);

%% ==============================
%  INIT RIGHT PANEL (PEAKS + CUMULATIVE |Δpeak|)
%  ==============================
R = init_right_panel(axRight, ...
    fontSizeAxes, fontSizeLabel, lineWidthPeaks, ...
    t, x_traj, maxCumDiff);

%% ==============================
%  LEGENDS — CENTERED ABOVE EACH BIG PANEL
%  ==============================
setup_legends(axLeft, axRight, L, R, ...
    xnull_color, ynull_color, lineWidthMain);

%% ==============================
%  VIDEO WRITER (OPTIONAL)
%  ==============================
if make_movie
    v = VideoWriter(movie_filename, 'MPEG-4');
    v.FrameRate = movie_fps;
    open(v);
end

%% ==============================
%  ANIMATION LOOP
%  ==============================
next_mark_time = t(1) + nullcline_interval;

for k = 2:frame_skip:N

    tk = t(k);
    s  = lambda_sign(tk);

    % ---- LEFT PANEL UPDATE ----
    set(L.bg, 'CData', s * lambda_base(L.R_bg));

    Lambda_grid = s * lambda_base(L.R0);
    Omega_grid  = omega(L.R0);

    Ux = Lambda_grid .* L.Xg2 - Omega_grid .* L.Yg2;
    Vy = Omega_grid .* L.Xg2 + Lambda_grid .* L.Yg2;

    Len = hypot(Ux, Vy);
    set(L.hQuiv, 'UData', Ux./(Len+eps), 'VData', Vy./(Len+eps));

    % Leave frozen nullclines every 10 cycles
    if tk >= next_mark_time
        contour(axLeft, L.Xg2, L.Yg2, Ux, [0 0], ...
            'Color', xnull_color, 'LineWidth', lineWidthMain, ...
            'HandleVisibility','off');
        contour(axLeft, L.Xg2, L.Yg2, Vy, [0 0], ...
            'Color', ynull_color, 'LineWidth', lineWidthMain, ...
            'HandleVisibility','off');

        next_mark_time = next_mark_time + nullcline_interval;
    end

    % Dynamic nullclines (latest only)
    delete(L.h1_nc); delete(L.h2_nc);
    [~, L.h1_nc] = contour(axLeft, L.Xg2, L.Yg2, Ux, [0 0], ...
        'Color', xnull_color, 'LineWidth', lineWidthMain, ...
        'HandleVisibility','off');
    [~, L.h2_nc] = contour(axLeft, L.Xg2, L.Yg2, Vy, [0 0], ...
        'Color', ynull_color, 'LineWidth', lineWidthMain, ...
        'HandleVisibility','off');

    % Trajectory + moving point
    set(L.hTrajLeft,  'XData', x_traj(1:k), 'YData', y_traj(1:k));
    set(L.hPointLeft, 'XData', x_traj(k),   'YData', y_traj(k));

    % Peaks + connections
    idx_last = find(t_peaks <= tk, 1, 'last');

    if ~isempty(idx_last)
        % Phase-plane peaks
        x_seg_phase = [x_traj(1) x_peaks_phase(1:idx_last)];
        y_seg_phase = [y_traj(1) y_peaks_phase(1:idx_last)];
        set(L.hPeakLinePhase, 'XData', x_seg_phase, 'YData', y_seg_phase);
        set(L.hPeakPtsPhase,  'XData', x_peaks_phase(1:idx_last), ...
                              'YData', y_peaks_phase(1:idx_last));

        % Right panel: peaks time-series
        t_seg = [t(1) t_peaks(1:idx_last)];
        x_seg = [x_traj(1) x_peaks(1:idx_last)];
        yyaxis(axRight,'left');
        set(R.hPeakLine, 'XData', t_seg, 'YData', x_seg);

        % Right panel: cumulative |diff|
        if idx_last >= 2
            diffs    = abs(diff(x_peaks(1:idx_last)));
            cumDiffs = cumsum(diffs);
            t_diff   = t_peaks(2:idx_last);
            yyaxis(axRight,'right');
            set(R.hCumLine, 'XData', t_diff, 'YData', cumDiffs);
        else
            yyaxis(axRight,'right');
            set(R.hCumLine, 'XData', NaN, 'YData', NaN);
        end
    else
        set(L.hPeakLinePhase,'XData',NaN,'YData',NaN);
        set(L.hPeakPtsPhase, 'XData',NaN,'YData',NaN);

        yyaxis(axRight,'left');
        set(R.hPeakLine, 'XData',NaN,'YData',NaN);
        yyaxis(axRight,'right');
        set(R.hCumLine,  'XData',NaN,'YData',NaN);
    end

    drawnow;

    if make_movie
        frame = getframe(fig);
        writeVideo(v, frame);
    end

    pause(pause_time);
end

%% ==============================
%  ADD SHADED ADAPTATION REGION AFTER ANIMATION
%  ==============================
if ~isnan(t_adapt_start)
    yyaxis(axRight,'right');
    yL = ylim(axRight);

    set(R.hAdapt, ...
        'XData', [t_adapt_start t_adapt_end t_adapt_end t_adapt_start], ...
        'YData', [yL(1)         yL(1)        yL(2)        yL(2)], ...
        'Visible','on');

    if make_movie
        drawnow;
        frame = getframe(fig);
        writeVideo(v, frame);
    end
end

%% ==============================
%  CLOSE VIDEO WRITER
%  ==============================
if make_movie
    close(v);
end

%% ==============================
%  LOCAL FUNCTIONS
%  ==============================
function L = init_left_panel(ax, xy_lim, nGrid_bg, nGrid_field, ...
    fontSizeAxes, fontSizeLabel, ...
    lineWidthMain, lineWidthTrajectory, lineWidthPeaks, markerSizePeaks, ...
    stable_color, unstable_color, xnull_color, ynull_color, traj_color, ...
    lambda_base, omega, lambda_sign, t0, x0)

    axes(ax); %#ok<LAXES>
    hold(ax,'on');
    set(ax,'FontSize',fontSizeAxes,'TickLabelInterpreter','latex');

    % Background shading grid
    x_range = linspace(-xy_lim, xy_lim, nGrid_bg);
    y_range = linspace(-xy_lim, xy_lim, nGrid_bg);
    [Xg, Yg] = meshgrid(x_range, y_range);
    R_bg     = hypot(Xg, Yg);

    s0         = lambda_sign(t0);
    Lambda_bg0 = s0 * lambda_base(R_bg);

    bg = imagesc(ax, x_range, y_range, Lambda_bg0);
    set(ax,'YDir','normal');
    colormap(ax, [stable_color; unstable_color]);
    caxis(ax, [-1 1]);
    set(bg,'AlphaData',0.5);

    % Vector field grid
    x_field = linspace(-xy_lim, xy_lim, nGrid_field);
    y_field = linspace(-xy_lim, xy_lim, nGrid_field);
    [Xg2, Yg2] = meshgrid(x_field, y_field);

    R0           = hypot(Xg2, Yg2);
    Lambda_grid0 = s0 * lambda_base(R0);
    Omega_grid0  = omega(R0);

    Ux0 = Lambda_grid0 .* Xg2 - Omega_grid0 .* Yg2;
    Vy0 = Omega_grid0 .* Xg2 + Lambda_grid0 .* Yg2;

    L0 = hypot(Ux0, Vy0);
    quivScale = 0.55;

    hQuiv = quiver(ax, Xg2, Yg2, ...
        Ux0./(L0+eps), Vy0./(L0+eps), quivScale, ...
        'Color',[0 0 0 0.6], ...
        'LineWidth',1.6, ...
        'MaxHeadSize',2.2, ...
        'HandleVisibility','off');

    % Dynamic nullclines (will be updated)
    [~, h1_nc] = contour(ax, Xg2, Yg2, Ux0, [0 0], ...
        'Color', xnull_color, 'LineWidth', lineWidthMain, 'HandleVisibility','off');
    [~, h2_nc] = contour(ax, Xg2, Yg2, Vy0, [0 0], ...
        'Color', ynull_color, 'LineWidth', lineWidthMain, 'HandleVisibility','off');

    % Frozen first nullclines
    contour(ax, Xg2, Yg2, Ux0, [0 0], ...
        'Color', xnull_color, 'LineWidth', lineWidthMain, 'HandleVisibility','off');
    contour(ax, Xg2, Yg2, Vy0, [0 0], ...
        'Color', ynull_color, 'LineWidth', lineWidthMain, 'HandleVisibility','off');

    % Trajectory + moving dot
    hTrajLeft = plot(ax, x0(1), x0(2), ...
        'Color', traj_color, 'LineWidth', lineWidthTrajectory, ...
        'DisplayName','Trajectory');

    movingDotSize = 10;
    hPointLeft = plot(ax, x0(1), x0(2), 'o', ...
        'MarkerFaceColor', traj_color, ...
        'Color', traj_color, ...
        'MarkerSize', movingDotSize, ...
        'HandleVisibility','off');

    % Peaks in phase plane
    hPeakLinePhase = plot(ax, NaN, NaN, 'k--', ...
        'LineWidth', lineWidthPeaks, 'HandleVisibility','off');
    hPeakPtsPhase  = plot(ax, NaN, NaN, 'ks', ...
        'MarkerFaceColor','k', 'MarkerSize', markerSizePeaks, ...
        'HandleVisibility','off');

    % Initial point marker
    plot(ax, x0(1), x0(2), 'ks', ...
        'MarkerFaceColor','k', 'MarkerSize', markerSizePeaks, ...
        'HandleVisibility','off');

    axis(ax,'equal');
    xlim(ax,[-xy_lim xy_lim]);
    ylim(ax,[-xy_lim xy_lim]);
    xlabel(ax,'$x$','FontSize',fontSizeLabel,'Interpreter','latex');
    ylabel(ax,'$y$','FontSize',fontSizeLabel,'Interpreter','latex');
    grid(ax,'on');
    box(ax,'on');

    % Return struct (IMPORTANT: bg is stored here)
    L.ax             = ax;
    L.bg             = bg;
    L.R_bg           = R_bg;
    L.Xg2            = Xg2;
    L.Yg2            = Yg2;
    L.R0             = R0;
    L.hQuiv          = hQuiv;
    L.h1_nc          = h1_nc;
    L.h2_nc          = h2_nc;
    L.hTrajLeft      = hTrajLeft;
    L.hPointLeft     = hPointLeft;
    L.hPeakLinePhase = hPeakLinePhase;
    L.hPeakPtsPhase  = hPeakPtsPhase;
end

function R = init_right_panel(ax, ...
    fontSizeAxes, fontSizeLabel, lineWidthPeaks, ...
    t, x_traj, maxCumDiff)

    axes(ax); %#ok<LAXES>
    hold(ax,'on');
    set(ax,'FontSize',fontSizeAxes,'TickLabelInterpreter','latex');

    % LEFT axis limits for peaks
    yyaxis(ax,'left');
    yPeaksMin = max(0, min(x_traj));
    yPeaksMax = 1.1*max(x_traj);
    ylim(ax, [yPeaksMin yPeaksMax]);
    ylabel(ax,'$x(t)$ at peaks','FontSize',fontSizeLabel,'Interpreter','latex');

    % RIGHT axis limits for cumulative |diff|
    yyaxis(ax,'right');
    ylim(ax, [0 1.1*maxCumDiff]);
    ylabel(ax,'Cumulative $|\Delta$peak$|$','FontSize',fontSizeLabel,'Interpreter','latex');

    % Shaded adaptation patch FIRST (hidden initially)
    yyaxis(ax,'right');
    yL = ylim(ax);
    adapt_color = [0.90 0.90 0.90];
    hAdapt = patch(ax, ...
        [t(1) t(1) t(1) t(1)], ...
        [yL(1) yL(1) yL(2) yL(2)], ...
        adapt_color, ...
        'FaceAlpha',0.3, ...
        'EdgeColor','none', ...
        'Visible','off', ...
        'HandleVisibility','off');

    % Curves
    yyaxis(ax,'left');
    hPeakLine = plot(ax, NaN, NaN, '-o', ...
        'Color',[0.5 0.5 0.5], ...
        'LineWidth', lineWidthPeaks, ...
        'DisplayName','Peaks');

    yyaxis(ax,'right');
    hCumLine = plot(ax, NaN, NaN, '-o', ...
        'Color','k', ...
        'LineWidth', lineWidthPeaks, ...
        'DisplayName','Cumulative $|\Delta$peak$|$');

    xlabel(ax,'$t$','FontSize',fontSizeLabel,'Interpreter','latex');
    grid(ax,'on');
    box(ax,'on');
    xlim(ax, [t(1) t(end)]);

    R.ax       = ax;
    R.hAdapt   = hAdapt;
    R.hPeakLine= hPeakLine;
    R.hCumLine = hCumLine;
end

function setup_legends(axLeft, axRight, L, R, xnull_color, ynull_color, lineWidthMain)
    % Make only these appear in legends
    hX_leg = plot(axLeft, NaN, NaN, ...
        'Color', xnull_color, 'LineWidth', lineWidthMain, ...
        'DisplayName','x-nullcline');
    hY_leg = plot(axLeft, NaN, NaN, ...
        'Color', ynull_color, 'LineWidth', lineWidthMain, ...
        'DisplayName','y-nullcline');

    leg1 = legend(axLeft, [hX_leg hY_leg L.hTrajLeft], ...
        'Interpreter','latex', ...
        'Orientation','horizontal', ...
        'Location','northoutside', ...
        'Box','on');

    % Center above left panel
    leg1.Position(1) = axLeft.Position(1) + (axLeft.Position(3) - leg1.Position(3)) / 2;
    leg1.Position(2) = leg1.Position(2) + 0.02;

    leg2 = legend(axRight, [R.hPeakLine R.hCumLine], ...
        'Interpreter','latex', ...
        'Orientation','horizontal', ...
        'Location','northoutside', ...
        'Box','on');

    % Center above right panel
    leg2.Position(1) = axRight.Position(1) + (axRight.Position(3) - leg2.Position(3)) / 2;
    leg2.Position(2) = leg2.Position(2) + 0.02;
end

function dxdt = lambda_omega_rhs_timevarying(t, x, lambda_base, omega, lambda_sign_handle)
    s   = lambda_sign_handle(t);
    r   = hypot(x(1), x(2));
    lam = s * lambda_base(r);
    om  = omega(r);

    dxdt = [lam*x(1) - om*x(2); ...
            om*x(1) + lam*x(2)];
end