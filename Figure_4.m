%% 2x3 figure: unified milestone formulas for threshold & high-density curves
clear; clc; close all;

figure('Color','w');

% Shared steepness in time transform
s2 = 0.1;  % appears in B/W normalization and t(z) mapping

% Panel specs (A-E)
P(1).label = '(a) \it{C. erythraeus}';
P(1).time  = [2001 2002 2003 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017]';
P(1).cost  = [0.0000227820; 0.0000336500; 0.0000451435; 0.0000679191; ...
              0.0002142905; 0.0002745281; 0.0003406172; 0.0004383024; ...
              0.0005151388; 0.0006380282; 0.0007316905; 0.0007748157; ...
              0.0008334545; 0.0009359811; 0.0009735314; 0.0009913621]*1e3;
P(1).ctype = 'threshold';   % (s1,aa) = (1, 2)

P(2).label = '(b) \it{H. javanicus}';
P(2).time  = [2000 2001 2002 2003 2008 2010 2011 2016 2017]';
P(2).cost  = [0.0001056708; 0.0001398438; 0.0001507118; 0.0001622053; ...
              0.0001732201; 0.0001988323; 0.0002261410; 0.0002355286; ...
              0.0002444439]*1e3;
P(2).ctype = 'threshold';

P(3).label = '(c) \it{M. coypus}';
P(3).time  = (2000:2017)';
P(3).cost  = [0.0011227525; 0.0019315128; 0.0027574837; 0.0037918926; ...
              0.0049553719; 0.0060599913; 0.0072205077; 0.0084654167; ...
              0.0098312580; 0.0112110614; 0.0124788627; 0.0140081528; ...
              0.0153328050; 0.0162492165; 0.0168551508; 0.0173421520; ...
              0.0179523441; 0.0184694342]*1e3;
P(3).ctype = 'highdensity'; % (s1,aa) = (0.75, 1+exp(-2.5))

P(4).label = '(d) \it{P. larvata}';
P(4).time  = (2000:2017)';
P(4).cost  = [0.0005943984; 0.0015740235; 0.0023347862; 0.0042082156; ...
              0.0063072763; 0.0083912490; 0.0107959225; 0.0133158593; ...
              0.0168295962; 0.0207247934; 0.0255526833; 0.0301678621; ...
              0.0359615432; 0.0406945389; 0.0451999535; 0.0489165414; ...
              0.0529344219; 0.0566521213]*1e3;
P(4).ctype = 'threshold';

P(5).label = '(e) \it{P. lotor}';
P(5).time  = (2000:2017)';
P(5).cost  = [0.0004755187; 0.0008855944; 0.0017333013; 0.0026412825; ...
              0.0041885900; 0.0059537035; 0.0076683403; 0.0097866935; ...
              0.0119456040; 0.0153523752; 0.0198601133; 0.0250897390; ...
              0.0295453875; 0.0332002521; 0.0364644787; 0.0394035735; ...
              0.0425577974; 0.0454731157]*1e3;
P(5).ctype = 'highdensity';

% Legend placeholders (grab from first panel encountered)
legendCaptured = false;
hThresh = []; hMid = []; hSat = [];
hLast = []; h2050 = []; hUpper = [];

% jbfill presence
has_jbfill = exist('jbfill','file')==2;

% Loop over panels
for k = 1:5
    subplot(2,3,k); hold on;

    % (s1, aa) per curve type
    if strcmp(P(k).ctype,'threshold')
        s1 = 1;              % high-threshold curve
        aa = 2;
    else
        s1 = 0.75;           % high-density curve
        aa = 1 + exp(-2.5);
    end

    % Time / data
    time_year = P(k).time(:);
    cost_data = P(k).cost(:);
    time = time_year - time_year(1);
    n    = 2050 - time_year(1);
    t    = 0:0.1:n;

    % Transform constants (B, W)
    d = (1 - s1)/s2;
    B = 1/(1 + exp(1/s2 - d));
    W = (1 + exp(-d)) / (1 - B*(1 + exp(-d)));

    % Model
    modelfun1 = @(b1,tv) (b1(1) .* W .* ...
        (1 ./ (1 + ((1 - B)/B) .* ...
        exp(-(1 ./ (1 + (b1(2) - 1).*exp(-b1(3).*tv))) / s2)) - B));

    % Fit
    beta0 = [1 1 1];
    mdl = fitnlm(time, cost_data, modelfun1, beta0);

    % Params & prediction
    params = mdl.Coefficients.Estimate;
    Cmax = params(1); gamma = params(2); alp = params(3);
    tcol = t(:);
    [~, yci] = predict(mdl, tcol);
    CIu = yci(:,2);
    cost_curve = modelfun1(params, t);

    % Unified milestone formulas (work for both curve types)
    zthresh = s1 - 0.1 * log(5*aa - 1 + sqrt(5*aa*(5*aa - 2)));
    zmid    = s1 - 0.1 * log(aa/0.5 - 1);
    zsat    = s1 - 0.1 * log(aa/0.9 - 1);

    % Times from z-values
    t_crit = -(1/alp) * log((1/(gamma-1))*(1/zthresh - 1));
    t_mid  = -(1/alp) * log((1/(gamma-1))*(1/zmid    - 1));
    t_sat  = -(1/alp) * log((1/(gamma-1))*(1/zsat    - 1));

    % Costs at those times
    cost_crit = modelfun1(params, t_crit);
    cost_mid  = modelfun1(params, t_mid);
    cost_sat  = modelfun1(params, t_sat);

    % Inflection (z = s1); C = Cmax * (aa/2)
    t_inflect = -(1/alp) * log((1/(gamma-1))*(1/s1 - 1));
    cost_infl = Cmax * (aa/2);

    % Points at last report and 2050
    est_cost_last = cost_curve(10*time(end)+1);
    cost_2050     = cost_curve(10*n+1);
    CIu_2050      = CIu(end);

    % Plot shaded region to upper CI (fallback if jbfill not found)
    if has_jbfill
        jbfill(t, cost_curve, CIu', [0.3 0.3 0.3], 'none', 0, 0.1);
    else
        Xp = [t, fliplr(t)];
        Yp = [cost_curve, fliplr(CIu')];
        patch('XData',Xp,'YData',Yp,'FaceColor',[0.8 0.8 0.8], ...
              'EdgeColor','none','FaceAlpha',0.4);
    end

    % Main curve + observed data
    plot(t, cost_curve, 'k-', 'LineWidth', 2);
    plot(time, cost_data, 'kx', 'LineWidth', 2);

    % --- Smaller markers (all reduced slightly) ---
    hL  = plot(time(end), est_cost_last, 'k^', 'MarkerFaceColor','k','MarkerSize',8,'LineWidth',2);
    hP  = plot(n,         cost_2050,     'kp', 'MarkerFaceColor','k','MarkerSize',8,'LineWidth',2);
    hU  = plot(n,         CIu_2050,      'kh', 'MarkerFaceColor','k','MarkerSize',8,'LineWidth',2);

    hT  = plot(t_crit, cost_crit, 'ko', 'MarkerFaceColor','k','MarkerSize',8,'LineWidth',2);
    hM  = plot(t_mid,  cost_mid,  'ks', 'MarkerFaceColor','k','MarkerSize',9,'LineWidth',2);
    hS  = plot(t_sat,  cost_sat,  'kd', 'MarkerFaceColor','k','MarkerSize',8,'LineWidth',2);

    % Inflection marker (red star) — reduced as well
    %plot(t_inflect, cost_infl, 'r*', 'MarkerSize',10,'LineWidth',2);

    % Axes & title formatting
    if k == 2
        axis([0 n+1 0 0.45]);      % Force panel (b) Y-limit
    else
        axis([0 n+1 0 1.1*Cmax]);
    end
    axis square;

    % Reduce tick label fonts (make X a bit smaller to avoid 2040/2050 overlap)
    axh = gca;
    try
        axh.XAxis.FontSize = 9;    % changed from 10 -> 9
        axh.YAxis.FontSize = 12;
    catch
        % older MATLAB fallback
        set(axh,'FontSize',13);
    end

    % Titles in Arial, italic only (no bold)
    title(P(k).label, ...
          'FontSize',16, 'FontWeight','normal', ...
          'Interpreter','tex', 'FontName','Arial');

    % Calendar year ticks
    tick_years = unique([1995:10:2050, 2050]);
    set(axh,'XTick', tick_years - time_year(1), 'XTickLabel', tick_years);

    % Capture legend handles onceupper
    if ~legendCaptured
        hThresh = hT; hMid = hM; hSat = hS;
        hLast = hL; h2050 = hP; hUpper = hU;
        legendCaptured = true;
    end
    hold off;
end

% ---------- Legend in the vacant subpanel (row 2, col 3 -> index 6) ----------
subplot(2,3,6); cla; axis off; hold on;
L = legend([hThresh, hMid, hSat, hLast, h2050, hUpper], ...
    {'Threshold point','Midpoint', 'Near saturation point', ...
     'Estimated cost at last cost report', ...
     'Projected damage cost in 2050', ...
     'Upper projected damage cost in 2050'}, ...
    'Location','northwest', 'Interpreter','tex', 'FontName','Arial');
set(L,'Box','off','FontSize',16);
hold off;

% -------- Layout & spacing (match your other figure) --------
left   = 0.08;   right  = 0.05;   % outer left/right margins
bottom = 0.08;   top    = 0.08;   % outer bottom/top margins
hgap   = -0.2;   vgap   = 0.1;    % gaps between columns/rows

% Compute sizes for a 2x3 grid
W = (1 - left - right - 2*hgap) / 3;
H = (1 - top  - bottom - vgap ) / 2;

% Target positions [x y w h] (row1 = top row)
pos_target = [ ...
  left                 , bottom+H+vgap, W, H;  % (1)
  left + (W+hgap)      , bottom+H+vgap, W, H;  % (2)
  left + 2*(W+hgap)    , bottom+H+vgap, W, H;  % (3)
  left                 , bottom       , W, H;  % (4)
  left + (W+hgap)      , bottom       , W, H;  % (5)
  left + 2*(W+hgap)    , bottom       , W, H   % (6) legend panel
];

% Normalize & reorder axes to assign positions
ax = findobj(gcf,'Type','axes'); set(ax,'Units','normalized');
posCells = get(ax,'Position');
if iscell(posCells)
    Ppos = vertcat(posCells{:});
else
    Ppos = posCells;
end
[~, order] = sortrows([ -Ppos(:,2), Ppos(:,1) ]);
ax = ax(order);
for i = 1:min(6, numel(ax))
    set(ax(i), 'Position', pos_target(i,:));
end
