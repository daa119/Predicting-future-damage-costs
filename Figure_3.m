%% Logistic growth panels: 2x3 subplot layout (older MATLAB compatible) - styled to match other figure
clear; clc; close all;

% ---- style knobs (matching your other figure) ----
MS_PLUS = 9;      % marker size for r+
MS_STAR = 10;     % marker size for stars (both black & blue)
LW      = 2;      % line/marker linewidth
FS_X    = 12;     % x-axis tick label size
FS_Y    = 12;     % y-axis tick label size
FS_TIT  = 16;     % panel title size
FS_LEG  = 16;     % legend font size

% ---------- Data for the five species ----------
panels = struct( ...
  'letter', {'a','b','c','d','e'}, ...
  'species', {'Callosciurus erythraeus', 'Herpestes javanicus', 'Myocastor coypus', 'Paguma larvata', 'Procyon lotor'}, ...
  'gamma',   {1.9792, 1.1525, 2.0487, 1.9258, 2.9258}, ...
  'alpha',   {0.2172, 0.0848, 0.1264, 0.1513, 0.1400}, ...
  'g_ci',    {[1.80 2.16], [1.07 1.24], [2.01 2.09], [1.88 1.97], [2.80 3.05]}, ...
  'a_ci',    {[0.185 0.249], [-0.032 0.201], [0.118 0.134], [0.142 0.161], [0.130 0.150]} ...
);

figure('Color','w');

% legend handle for blue star (created later)
h_blue = [];

for k = 1:5
    subplot(2,3,k); hold on;

    % Panel params
    species   = panels(k).species;
    letter    = panels(k).letter;
    gamma     = panels(k).gamma;
    alpha     = panels(k).alpha;
    gamma_ci  = panels(k).g_ci;
    alpha_ci  = panels(k).a_ci;

    % Time grid and central curve
    t = -30:0.01:20;
    z = 1 ./ (1 + (gamma - 1) .* exp(-alpha .* t));

    % Confidence bounds (ensure non-negative lower alpha)
    a_low  = max(alpha_ci(1), 0);
    a_high = alpha_ci(2);
    g_low  = gamma_ci(1);
    g_high = gamma_ci(2);

    z_low  = 1 ./ (1 + (g_low  - 1) .* exp(-a_low  .* t));
    z_high = 1 ./ (1 + (g_high - 1) .* exp(-a_high .* t));

    % Density-time inflection (where z = 0.5)
    tinf = (1/alpha) * log(gamma - 1);

    % --- Confidence band ---
    c = [0.3 0.3 0.3];
    fill([t, fliplr(t)], [z_low, fliplr(z_high)], c, ...
        'EdgeColor','none', 'FaceAlpha', 0.10);

    % --- Main curve ---
    plot(t, z, 'k-', 'LineWidth', 1.8);

    % --- Key points (styled to match) ---
    plot(0,    1/gamma, 'r+', 'MarkerSize', MS_PLUS, 'LineWidth', LW);      % initial z(0)
    plot(tinf, 0.5,     'k*', 'MarkerSize', MS_STAR, 'LineWidth', LW);      % density-time inflection (z=0.5)

    % --- Cost–density inflection (blue star), only for Myocastor & Procyon ---
    if ismember(k, [3 5])            % panels (c) and (e)
        s1_cd = 0.75;
        t_cd  = -(1/alpha) * log((1/(gamma-1))*(1/s1_cd - 1));  % time where z = 0.75
        % blue star marker on the curve
        plot(t_cd, s1_cd, 'b*', 'MarkerSize', MS_STAR, 'LineWidth', LW);
        if isempty(h_blue)
            % dummy handle for legend with blue star:
            h_blue = plot(NaN, NaN, 'b*', 'MarkerSize', MS_STAR, 'LineWidth', LW);
        end
    end

    % --- Axes & ticks ---
    axis([t(1) t(end) 0 1]);

    if k == 3
        % Panel (c): show full x-axis labels from -30 to 20
        set(gca, 'XTick', -30:10:20, 'XTickLabel', -30:10:20);
    elseif ismember(k, [4 5])
        % Bottom row, other panels
        set(gca, 'XTick', -30:10:20);
    else
        set(gca, 'XTick', []);
    end

    if ismember(k, [1 4])   % first column
        set(gca, 'YTick', 0:0.2:1);
    else
        set(gca, 'YTick', []);
    end

    set(gca, 'TickLength', [0.008 0.008]);   % shorter tick marks

    % --- Axis font sizes (match other figure) ---
    axh = gca;
    try
        axh.XAxis.FontSize = FS_X;
        axh.YAxis.FontSize = FS_Y;
    catch
        set(axh, 'FontSize', max(FS_X,FS_Y)); % fallback for older MATLAB
    end

    % --- Title (panel letter NOT italic; species italic, Arial) ---
    % Abbreviate species to "C. erythraeus", etc.
    toks = strsplit(species, ' ');
    speciesShort = sprintf('%c. %s', toks{1}(1), toks{2});
    title(sprintf('(%s) \\it %s', letter, speciesShort), ...
          'FontSize', FS_TIT, 'FontWeight', 'normal', ...
          'FontAngle', 'normal', 'Interpreter', 'tex', 'FontName', 'Arial');

    axis square;

    % --- In-panel text (LaTeX) ---
    xbr = t(end) - 0.5;   % inset from right edge
    ybr = 0 + 0.045;      % inset from bottom
    txtLines = {sprintf('$z(0)=%.3f$', 1/gamma); sprintf('$z(%.2f)=0.5$', tinf)};
    if k == 3
        % (c) Myocastor coypus: show z(9.07)=0.75
        txtLines{end+1} = '$z(9.07)=0.75$';
    elseif k == 5
        % (e) Procyon lotor: show z(12.53)=0.75
        txtLines{end+1} = '$z(12.53)=0.75$';
    end
    text(xbr, ybr, txtLines, 'Interpreter','latex', ...
         'HorizontalAlignment','right', 'VerticalAlignment','bottom', ...
         'FontSize', 14);

    hold off;
end

% ---------- Blank 6th subplot: shared legend (Arial) ----------
subplot(2,3,6); cla; axis off; hold on;
h_plus = plot(NaN, NaN, 'r+', 'MarkerSize', MS_PLUS, 'LineWidth', LW);  % z(0)
h_star = plot(NaN, NaN, 'k*', 'MarkerSize', MS_STAR, 'LineWidth', LW);  % density-time inflection

% Always include the blue star entry (handle exists if plotted)
if isempty(h_blue)
    lg_handles = [h_plus, h_star];
else
    lg_handles = [h_plus, h_star, h_blue];
end

% TeX-safe labels so they render correctly with Arial
lg_labels  = { ...
    'Initial density', ...
    'Density-time inflection', ...
    'Cost-density inflection' ...
};

% Legend in Arial (Interpreter='tex' so FontName is honored)
lg6 = legend(lg_handles, lg_labels(1:numel(lg_handles)), ...
    'Interpreter','tex', 'Location','northwest', 'FontName','Arial');
set(lg6, 'Box','off', 'FontSize', FS_LEG, 'FontWeight','normal');
hold off;

%% ================== MAKE SUBPANELS CLOSER (layout only) ==================
% Adjust these for tighter/looser layout (normalized units)
left   = 0.08;   right  = 0.05;   % outer left/right margins
bottom = 0.08;   top    = 0.08;   % outer bottom/top margins
hgap   = -0.2;   vgap   = 0.1;    % gaps between columns/rows

% Compute cell sizes for a 2x3 grid
W = (1 - left - right - 2*hgap) / 3;   % width of each column
H = (1 - top  - bottom - vgap ) / 2;   % height of each row

% Target positions [x y w h]; row 1 is the top row (subplots 1..3)
pos_target = [ ...
  left                 , bottom+H+vgap, W, H;  % subplot 1 (row1,col1)
  left + (W+hgap)      , bottom+H+vgap, W, H;  % subplot 2 (row1,col2)
  left + 2*(W+hgap)    , bottom+H+vgap, W, H;  % subplot 3 (row1,col3)
  left                 , bottom       , W, H;  % subplot 4 (row2,col1)
  left + (W+hgap)      , bottom       , W, H;  % subplot 5 (row2,col2)
  left + 2*(W+hgap)    , bottom       , W, H   % subplot 6 (row2,col3)
];

% Get all axes in the figure and normalize units
ax = findobj(gcf,'Type','axes');
set(ax, 'Units','normalized');

% Collect current positions into an N-by-4 array
posCells = get(ax, 'Position');
if iscell(posCells)
    P = vertcat(posCells{:});  % N x 4
else
    P = posCells;              % 1 x 4 (rare)
end

% Sort axes top-to-bottom (y DESC) and left-to-right (x ASC)
[~, order] = sortrows([ -P(:,2), P(:,1) ]);
ax = ax(order);

% Apply the tight target positions
nApply = min(6, numel(ax));    % safeguard
for i = 1:nApply
    set(ax(i), 'Position', pos_target(i,:));
end
