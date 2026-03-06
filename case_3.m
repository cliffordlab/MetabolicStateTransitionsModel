%% case_3.m
% Order-2 lambda–omega system with animated trajectory and time-varying stability.
%
% LAYOUT (2x3 pattern with spans):
%   LEFT  panel spans tiles (2,3,1) (2,3,2) (2,3,4) (2,3,5)  => nexttile(1,[2 2])
%   RIGHT panel spans tiles (2,3,3) (2,3,6)                  => nexttile(3,[2 1])
%
% LEFT  : Phase plane (x–y) with nullclines, vector field, trajectory, peaks, initial point
% RIGHT : Time-series of peak values of x(t) (gray, left y-axis)
%         and cumulative sum of absolute differences between consecutive peaks (black, right y-axis)
%
% Model:
%   Λ(r,t) = λ(t) - b r^2
%   Ω(r,t) = ω0 + a r^2
% with λ(t) chosen so that the limit-cycle radius shrinks smoothly over time.

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
xy_lim      = 1.25;      % both x and y in [-xy_lim, xy_lim]
nGrid_bg    = 200;       % resolution of background shading
nGrid_field = 20;        % resolution of vector field

% Animation controls
frame_skip  = 40;        % draw every 'frame_skip'-th time sample
pause_time  = 0.01;      % pause between frames (s)

%% ==============================
%  MOVIE SETTINGS
%  ==============================
make_movie      = true;                               % set false for no video
movie_filename  = 'Movie/case_3_movie.mp4';
movie_fps       = 10;                                 % frames per second

%% ==============================
%  ORDER-2 LAMBDA–OMEGA PARAMETERS (SHRINKING LIMIT CYCLE)
%  ==============================
omega0 = 1.0;   % base frequency term
a      = 1.0;   % quadratic freq term
b      = 1.0;   % quadratic amp term

% Time-varying limit-cycle radius r̄(t) that SHRINKS smoothly
R_max      = 1.0;    % initial large radius
R_min      = 0.6;    % final small radius
tau_shrink = 10.0;   % time constant for exponential decay

% Desired LC radius (shrinking from R_max to R_min)
R_target = @(t) R_min + (R_max - R_min) * exp(-t / tau_shrink);

% Time-varying lambda(t) chosen so that LC radius ~ R_target(t)
lambda_t = @(t) b * (R_target(t)).^2;

% Define Λ(r,t) and Ω(r,t)
Lambda_fn = @(r,t) lambda_t(t) - b .* r.^2;
Omega_fn  = @(r,t) omega0 + a .* r.^2;

%% ==============================
%  TIME SPAN & INITIAL CONDITION
%  ==============================
tspan = [0 60];
x0    = [R_max; 0.0];   % start on the large (initial) limit cycle

%% ==============================
%  NULLCLINE SNAPSHOTS EVERY 10 CYCLES
%  ==============================
Omega_init        = Omega_fn(R_max, 0);
osc_period        = 2*pi / Omega_init;
cycles_per_mark   = 10;
nullcline_interval= osc_period * cycles_per_mark;

%% ==============================
%  SOLVE ODE
%  ==============================
odefun = @(t, x) lambda_omega_rhs_order2_shrinkingLC(t, x, Lambda_fn, Omega_fn);
opts   = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'MaxStep',0.01);

[t, X]  = ode45(odefun, tspan, x0, opts);
x_traj  = X(:,1);
y_traj  = X(:,2);
N       = numel(t);

%% ==============================
%  PEAK DETECTION
%  ==============================
[peakVals, peakTimes] = findpeaks(x_traj, t);
t_peaks = peakTimes(:).';
x_peaks = peakVals(:).';

peak_idx_in_traj = zeros(1,numel(t_peaks));
for i=1:numel(t_peaks)
    [~, peak_idx_in_traj(i)] = min(abs(t - t_peaks(i)));
end
x_peaks_phase = x_traj(peak_idx_in_traj).';
y_peaks_phase = y_traj(peak_idx_in_traj).';

if numel(x_peaks) >= 2
    diffsAll   = abs(diff(x_peaks));
    maxCumDiff = sum(diffsAll);
else
    maxCumDiff = 1;
end

%% ==============================
%  ADAPTATION WINDOW BASED ON >0.1 diff (SAFE BOUNDS)
%  ==============================
diffThreshold = 0.1;
leverage      = 4;

t_adapt_start = NaN;
t_adapt_end   = NaN;

if numel(x_peaks) >= 2
    diffs     = abs(diff(x_peaks));
    idx_above = find(diffs > diffThreshold);

    if ~isempty(idx_above)
        firstPairIdx = idx_above(1);
        lastPairIdx  = idx_above(end);

        idx_start_peak = max(1, firstPairIdx - leverage);
        idx_end_peak   = min(numel(t_peaks), lastPairIdx + 1 + leverage);

        t_adapt_start = t_peaks(idx_start_peak);
        t_adapt_end   = t_peaks(idx_end_peak);
    end
end

%% ==============================
%  FIGURE / AXES SETUP  (2x3 PATTERN WITH SPANS)
%  ==============================
fig = figure('Color','w','Position',[100 100 1900 980]);

tl = tiledlayout(fig, 2, 3);
tl.Padding     = 'compact';
tl.TileSpacing = 'compact';


% OPTIONAL: explicit outer margins (uncomment to control)
mLeft=0.01; mRight=0.1; mBottom=0.2; mTop=0.2;
tl.OuterPosition = [mLeft, mBottom, 1-mLeft-mRight, 1-mTop-mBottom];
tl.InnerPosition = tl.OuterPosition;


% LEFT panel spans tiles (1,2,4,5)
ax1 = nexttile(tl, 1, [2 2]);
hold(ax1,'on');
set(ax1,'FontSize',fontSizeAxes,'TickLabelInterpreter','latex');

% RIGHT panel spans tiles (3,6)
ax3 = nexttile(tl, 3, [2 1]);
hold(ax3,'on');
set(ax3,'FontSize',fontSizeAxes,'TickLabelInterpreter','latex');

%% ==============================
%  LEFT PANEL — PHASE PLANE
%  ==============================
x_range = linspace(-xy_lim, xy_lim, nGrid_bg);
y_range = linspace(-xy_lim, xy_lim, nGrid_bg);
[Xg, Yg] = meshgrid(x_range, y_range);
R_bg     = hypot(Xg,Yg);

Lambda_bg0 = Lambda_fn(R_bg, t(1));
bg = imagesc(ax1,x_range,y_range,Lambda_bg0);
set(ax1,'YDir','normal');
colormap(ax1,[stable_color; unstable_color]);
caxis(ax1,[-1 1]);
set(bg,'AlphaData',0.5);

% Vector field
x_field = linspace(-xy_lim, xy_lim, nGrid_field);
y_field = linspace(-xy_lim, xy_lim, nGrid_field);
[Xg2,Yg2] = meshgrid(x_field,y_field);

R0           = hypot(Xg2,Yg2);
Lambda_grid0 = Lambda_fn(R0,t(1));
Omega_grid0  = Omega_fn(R0,t(1));

Ux0 = Lambda_grid0 .* Xg2 - Omega_grid0 .* Yg2;
Vy0 = Omega_grid0 .* Xg2 + Lambda_grid0 .* Yg2;
L0  = hypot(Ux0,Vy0);

quivScale = 0.55;

hQuiv = quiver(ax1, Xg2, Yg2, ...
    Ux0./(L0+eps), Vy0./(L0+eps), quivScale, ...
    'Color',[0 0 0 0.6], ...
    'LineWidth', 1.6, ...
    'MaxHeadSize', 2.2, ...
    'HandleVisibility','off');

% Initial nullclines (dynamic)
[~,h1_nc] = contour(ax1,Xg2,Yg2,Ux0,[0 0], ...
    'Color',xnull_color,'LineWidth',lineWidthMain,'HandleVisibility','off');
[~,h2_nc] = contour(ax1,Xg2,Yg2,Vy0,[0 0], ...
    'Color',ynull_color,'LineWidth',lineWidthMain,'HandleVisibility','off');

% Frozen snapshot
contour(ax1,Xg2,Yg2,Ux0,[0 0],'Color',xnull_color,'LineWidth',lineWidthMain,'HandleVisibility','off');
contour(ax1,Xg2,Yg2,Vy0,[0 0],'Color',ynull_color,'LineWidth',lineWidthMain,'HandleVisibility','off');

% Trajectory + moving point
hTrajLeft = plot(ax1,x_traj(1),y_traj(1), ...
    'Color',traj_color,'LineWidth',lineWidthTrajectory, ...
    'DisplayName','Trajectory');

movingDotSize = 10;
hPointLeft = plot(ax1,x_traj(1),y_traj(1),'o', ...
    'MarkerFaceColor',traj_color,'Color',traj_color, ...
    'MarkerSize',movingDotSize, ...
    'HandleVisibility','off');

% Peaks on phase plane
hPeakLinePhase = plot(ax1,NaN,NaN,'k--', ...
    'LineWidth',lineWidthPeaks,'HandleVisibility','off');
hPeakPtsPhase  = plot(ax1,NaN,NaN,'ks', ...
    'MarkerFaceColor','k','MarkerSize',markerSizePeaks,'HandleVisibility','off');
hInitPtPhase   = plot(ax1,x_traj(1),y_traj(1),'ks', ...
    'MarkerFaceColor','k','MarkerSize',markerSizePeaks,'HandleVisibility','off');

axis(ax1,'equal');
xlim(ax1,[-xy_lim xy_lim]);
ylim(ax1,[-xy_lim xy_lim]);
xlabel(ax1,'$x$','FontSize',fontSizeLabel,'Interpreter','latex');
ylabel(ax1,'$y$','FontSize',fontSizeLabel,'Interpreter','latex');
grid(ax1,'on');
box(ax1,'on');

% Legend anchors
hX_leg = plot(ax1,NaN,NaN,'Color',xnull_color,'LineWidth',lineWidthMain,'DisplayName','x-nullcline');
hY_leg = plot(ax1,NaN,NaN,'Color',ynull_color,'LineWidth',lineWidthMain,'DisplayName','y-nullcline');

%% ==============================
%  RIGHT PANEL — PEAKS + CUMULATIVE DIFF
%  ==============================
yyaxis(ax3,'left');
yPeaksMin = max(0,1.5*min(x_traj));
yPeaksMax = 1.1*max(x_traj);
ylim(ax3,[yPeaksMin yPeaksMax]);
ylabel(ax3,'$x(t)$ at peaks','FontSize',fontSizeLabel,'Interpreter','latex');

yyaxis(ax3,'right');
yCumMin = 0.08;
yCumMax = 1.2*maxCumDiff;
ylim(ax3,[yCumMin yCumMax]);
ylabel(ax3,'Cumulative $|\Delta$peak$|$','FontSize',fontSizeLabel,'Interpreter','latex');

% Adaptation patch (create first, hidden)
yyaxis(ax3,'right');
yL = ylim(ax3);
adapt_color = [0.90 0.90 0.90];
hAdapt = patch(ax3, ...
    [t(1) t(1) t(1) t(1)], ...
    [yL(1) yL(1) yL(2) yL(2)], ...
    adapt_color, ...
    'FaceAlpha',0.3, ...
    'EdgeColor','none', ...
    'Visible','off', ...
    'HandleVisibility','off');

yyaxis(ax3,'left');
hPeakLine = plot(ax3,NaN,NaN,'-o', ...
    'Color',[0.5 0.5 0.5],'LineWidth',lineWidthPeaks, ...
    'DisplayName','Peaks');

yyaxis(ax3,'right');
hCumLine = plot(ax3,NaN,NaN,'-o', ...
    'Color','k','LineWidth',lineWidthPeaks, ...
    'DisplayName','Cumulative $|\Delta$peak$|$');

xlabel(ax3,'$t$','FontSize',fontSizeLabel,'Interpreter','latex');
grid(ax3,'on');
box(ax3,'on');
xlim(ax3,[t(1) t(end)]);

%% ==============================
%  LEGENDS — CENTERED ABOVE EACH PANEL
%  ==============================
set([h1_nc h2_nc hQuiv hPointLeft ...
     hPeakLinePhase hPeakPtsPhase hInitPtPhase], ...
    'HandleVisibility','off');

leg1 = legend(ax1, [hX_leg hY_leg hTrajLeft], ...
    'Interpreter','latex', ...
    'Orientation','horizontal', ...
    'Location','northoutside', ...
    'Box','on');
leg1.Position(1) = ax1.Position(1) + (ax1.Position(3) - leg1.Position(3))/2;
leg1.Position(2) = leg1.Position(2) + 0.05;

leg2 = legend(ax3, [hPeakLine hCumLine], ...
    'Interpreter','latex', ...
    'Orientation','horizontal', ...
    'Location','northoutside', ...
    'Box','on');
leg2.Position(1) = ax3.Position(1) + (ax3.Position(3) - leg2.Position(3))/2;
leg2.Position(2) = leg2.Position(2) + 0.05;

%% ==============================
%  SET UP VIDEO WRITER (OPTIONAL)
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

for k=2:frame_skip:N

    tk = t(k);

    % Background field update
    set(bg,'CData',Lambda_fn(R_bg,tk));

    Lambda_grid = Lambda_fn(R0,tk);
    Omega_grid  = Omega_fn(R0,tk);

    Ux = Lambda_grid .* Xg2 - Omega_grid .* Yg2;
    Vy = Omega_grid .* Xg2 + Lambda_grid .* Yg2;
    L  = hypot(Ux,Vy);

    set(hQuiv,'UData',Ux./(L+eps), 'VData',Vy./(L+eps));

    % Leave frozen nullclines every ~10 cycles
    if tk >= next_mark_time
        contour(ax1,Xg2,Yg2,Ux,[0 0], ...
            'Color',xnull_color,'LineWidth',lineWidthMain,'HandleVisibility','off');
        contour(ax1,Xg2,Yg2,Vy,[0 0], ...
            'Color',ynull_color,'LineWidth',lineWidthMain,'HandleVisibility','off');
        next_mark_time = next_mark_time + nullcline_interval;
    end

    delete(h1_nc); delete(h2_nc);
    [~,h1_nc] = contour(ax1,Xg2,Yg2,Ux,[0 0], ...
        'Color',xnull_color,'LineWidth',lineWidthMain,'HandleVisibility','off');
    [~,h2_nc] = contour(ax1,Xg2,Yg2,Vy,[0 0], ...
        'Color',ynull_color,'LineWidth',lineWidthMain,'HandleVisibility','off');

    % Trajectory updates
    set(hTrajLeft,'XData',x_traj(1:k),'YData',y_traj(1:k));
    set(hPointLeft,'XData',x_traj(k),'YData',y_traj(k));

    idx_last = find(t_peaks <= tk, 1,'last');

    if ~isempty(idx_last)
        % phase-plane peaks
        x_seg_phase = [x_traj(1) x_peaks_phase(1:idx_last)];
        y_seg_phase = [y_traj(1) y_peaks_phase(1:idx_last)];
        set(hPeakLinePhase,'XData',x_seg_phase,'YData',y_seg_phase);
        set(hPeakPtsPhase,'XData',x_peaks_phase(1:idx_last), ...
                          'YData',y_peaks_phase(1:idx_last));

        % time peaks
        t_seg = [t(1) t_peaks(1:idx_last)];
        x_seg = [x_traj(1) x_peaks(1:idx_last)];
        yyaxis(ax3,'left');
        set(hPeakLine,'XData',t_seg,'YData',x_seg);

        % cumulative |diff|
        if idx_last >= 2
            diffs    = abs(diff(x_peaks(1:idx_last)));
            cumDiffs = cumsum(diffs);
            t_diff   = t_peaks(2:idx_last);
            yyaxis(ax3,'right');
            set(hCumLine,'XData',t_diff,'YData',cumDiffs);
        else
            yyaxis(ax3,'right');
            set(hCumLine,'XData',NaN,'YData',NaN);
        end
    else
        set(hPeakLinePhase,'XData',NaN,'YData',NaN);
        set(hPeakPtsPhase,'XData',NaN,'YData',NaN);

        yyaxis(ax3,'left');  set(hPeakLine,'XData',NaN,'YData',NaN);
        yyaxis(ax3,'right'); set(hCumLine,'XData',NaN,'YData',NaN);
    end

    drawnow;

    if make_movie
        frame = getframe(fig);
        writeVideo(v, frame);
    end

    pause(pause_time);
end

%% ==============================
%  ADD SHADED ADAPTATION REGION AFTER ANIMATION (OPTIONAL)
%  ==============================
if ~isnan(t_adapt_start)
    yyaxis(ax3,'right');
    yL = ylim(ax3);
    set(hAdapt, ...
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
%  LOCAL FUNCTION (SHRINKING LC)
%  ==============================
function dxdt = lambda_omega_rhs_order2_shrinkingLC(t, x, Lambda_fn, Omega_fn)
    r   = hypot(x(1), x(2));
    lam = Lambda_fn(r, t);
    om  = Omega_fn(r, t);

    dxdt = [lam*x(1) - om*x(2);
            om*x(1) + lam*x(2)];
end