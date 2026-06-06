%% Display radar plot per defect category

% Input:    values: Numeric values to be displayed in the radar plot
%                   - rows = pelvis defects
%                   - columns = segments / radar axes
%           plotTitle: Title of the radar plot
%           ids: Pelvis IDs corresponding to rows of values
%           labels: Cell array/string array with axis labels
%                   - number of labels must match number of value columns
%
% Output:   None
%                   - creates radar plot in current axes
%                   - plots one curve per pelvis defect and the category mean

% Developed by C.Micheler,
% Department of Orthopaedics and Sports Orthopaedics, TUM School of Medicine and Health, Klinikum rechts der Isar, Technical University of Munich
% Institute for Machine Tools and Industrial Management, TUM School of Engineering and Design, Technical University of Munich


function plotRadarCategory(values, plotTitle, ids, labels)
% values: rows = pelvis defects, cols = segments
% ids: pelvis IDs corresponding to rows
% labels: labels corresponding to columns of values

rMax = 100;
rTicks = [0 25 50 75 100];

% Check dimensions
nAxes = size(values, 2);

if numel(labels) ~= nAxes
    error('plotRadarCategory:SizeMismatch', ...
        'Number of labels (%d) must match number of value columns (%d).', ...
        numel(labels), nAxes);
end

if numel(ids) ~= size(values, 1)
    error('plotRadarCategory:IDMismatch', ...
        'Number of IDs (%d) must match number of value rows (%d).', ...
        numel(ids), size(values, 1));
end

labels = labels(:)';

% Dynamic order around border:
% first label at top, then clockwise
thetaDeg = 90 - (0:nAxes-1) * 360/nAxes;
theta = deg2rad(thetaDeg);
thetaClosed = [theta theta(1)];

hold on
axis equal
axis off

% Grid polygons
for r = rTicks(2:end)
    xg = r * cos(thetaClosed);
    yg = r * sin(thetaClosed);
    plot(xg, yg, '-', ...
        'Color', [0.75 0.75 0.75], ...
        'LineWidth', 1);
end

% Axes
for k = 1:nAxes
    plot([0 rMax*cos(theta(k))], ...
        [0 rMax*sin(theta(k))], '-', ...
        'Color', [0.6 0.6 0.6], ...
        'LineWidth', 1);
end

% Plot all pelvis defects of this category
nCurves = size(values, 1);
hLines = gobjects(nCurves, 1);

for iCurve = 1:nCurves
    valuesCur = values(iCurve, :);
    valuesClosed = [valuesCur valuesCur(1)];

    xPoly = valuesClosed .* cos(thetaClosed);
    yPoly = valuesClosed .* sin(thetaClosed);

    hLines(iCurve) = plot(xPoly, yPoly, '-', ...
        'LineWidth', 1.5);

    plot(xPoly, yPoly, '.', ...
        'Color', hLines(iCurve).Color, ...
        'MarkerSize', 10);
end

% Mean curve of category
meanValues = mean(values, 1, 'omitnan');
meanClosed = [meanValues meanValues(1)];

xMean = meanClosed .* cos(thetaClosed);
yMean = meanClosed .* sin(thetaClosed);

hMean = plot(xMean, yMean, 'k-', ...
    'LineWidth', 3);

plot(xMean, yMean, 'ko', ...
    'MarkerFaceColor', 'k', ...
    'MarkerSize', 5);

% Labels
labelRadius = rMax * 1.12;

for k = 1:nAxes
    xl = labelRadius * cos(theta(k));
    yl = labelRadius * sin(theta(k));

    text(xl, yl, labels{k}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', ...
        'FontSize', 12);
end

% Tick labels
for r = rTicks
    text(-8, r, num2str(r), ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FontWeight', 'bold', ...
        'FontSize', 10, ...
        'Color', [0.1 0.1 0.1]);
end

xlim([-1.25*rMax, 1.25*rMax]);
ylim([-1.15*rMax, 1.15*rMax]);

title(plotTitle, 'FontWeight', 'bold');

% Legend with pelvis IDs
legendStrings = cell(nCurves + 1, 1);

for iCurve = 1:nCurves
    legendStrings{iCurve} = ['Pelvis ', num2str(ids(iCurve))];
end

legendStrings{end} = 'Mean';

legend([hLines; hMean], legendStrings, ...
    'Location', 'bestoutside');

hold off
end