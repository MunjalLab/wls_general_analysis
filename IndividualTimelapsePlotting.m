%% Script to plot pillar angle and OV volume data in time from individually timelapsed Zebrafish OVs
% Kira L. Heikes 20260504
% Made with the help of Copilot GPT-5
% For use in manuscript: Hydrostatic pressure shapes and canalizes
% semicircular canal morphology to ensure vestibular function.
% doi.org/10.64898/2026.07.22.740136

clear; clc;


%% I. USER INPUT FIELDS
% Time conversion settings
T_ref   = 9;     % reference time index (e.g. Timelapse starts at timepoint 9)
hpf_ref = 54.5;   % hpf at reference index (e.g. timepoint 9 corresponds to 54.5 hours post fertilization (hpf)
dt_hpf  = 0.5;    % hpf per index step (stepsize of timepoints e.g. 0.5 hours)

% Set true to force single combined plot or false to keep plots separated by detected base names in csv entries (from 'Embryo' ID column)
forceSingleFigure = false;

% Choose plotting colors
angleColor  = [0, 0.620, 0.451]; % green
volumeColor = [0.8980, 0.6275, 0.1412] ; % orange

% User prompted selection of csv file
[fileName, filePath] = uigetfile({'*.csv','CSV files (*.csv)'; '*.*','All files'}, ...
                                 'Select timelapse CSV');

% check if file was selected
if isequal(fileName,0)
    error('No file selected.');
end % script ends

%% II. Select CSV file and read data
csvFile = fullfile(filePath, fileName);

% Read table (preserve headers exactly as written in CSV)
opts = detectImportOptions(csvFile, 'Delimiter', ','); % assume it is comma-delimited
opts.VariableNamingRule = 'preserve';   % preserve the table column headers as they are
T = readtable(csvFile, opts); % Reads CSVs where each row is a timepoint of data from an individual timelapse containing a name (Embryo), a pillar angle, and an OV volume

% Diagnostics to print (comment out if preferred)
disp("VariableNames MATLAB sees:");
disp(string(T.Properties.VariableNames));

% Identify which columns contain each data type expected
getCol = @(tbl, candidates) findColumn(tbl, candidates); % Required CSV columns (or close variants): "Embryo", "Angle", "Volume"

% Identify names of columns for each data type
colEmbryo = getCol(T, ["Embryo"]);
colAngle  = getCol(T, ["Angle", "Angle (deg)", "Angle(deg)"]);
colVolume = getCol(T, ["Volume", "Volume (pl)", "Volume(pl)", "Volume (pL)", "Volume(pL)"]);

% Pull raw column data into tables using the matched names found above
embryoRaw = T.(colEmbryo);
angleRaw  = T.(colAngle);
volumeRaw = T.(colVolume);

% convert 'Embryo' column entries to strings
embryoStr = string(embryoRaw);

% Parse time index from Embryo strings (e.g. "..._T29" --> timeIndex = 29)
arrayOfTimepoints = regexp(embryoStr, '\d{2}(?=\D*)$', 'match'); % Time index is parsed from the "Embryo" column strings which end with "..._T29".
timeIndex = nan(height(T),1); % initialize array of time indexes

% loop through each entry in the data to find its timepoint and store in timeIndex array
for i = 1:height(T)
    if ~isempty(arrayOfTimepoints{i})
        timeIndex(i) = str2double(arrayOfTimepoints{i}{1});
    end
end

% find the base name for each entry (ID) by removing the timepoint identifier at the end (assuming this is 4 characters long on each row)
baseName = regexprep(embryoStr, '.{4}$', '', 'lineanchors');

% if any entries do not have an identifiable time index, report this as the number of entries with this issue
badTime = isnan(timeIndex);
if any(badTime)
    warning('Could not parse time index from %d row(s). They will be ignored.', sum(badTime));
end


%% III. Build cleaned table and group data for plotting
Tc = table;
Tc.Base      = baseName(~badTime);
Tc.TimeIndex = timeIndex(~badTime);

% Converts time index to absolute time in hpf:
Tc.Time_hpf = hpf_ref + (Tc.TimeIndex - T_ref) * dt_hpf;

% Pull all angle and volume data from rows with identifiable timepoints
ang = angleRaw(~badTime);
vol = volumeRaw(~badTime);

% Convert angle and volume data to numeric data type if not
if ~isnumeric(ang), ang = str2double(string(ang)); end
if ~isnumeric(vol), vol = str2double(string(vol)); end

% if any entries do not have a numeric angle and/or volume value after the above, report this as the number of entries with this issue
badAV = isnan(ang) | isnan(vol);
if any(badAV)
    warning('Dropping %d row(s) with non-numeric Angle or Volume.', sum(badAV));
end

% pull all angle and volume data from rows with numeric values for both
Tc = Tc(~badAV,:);
ang = ang(~badAV);
vol = vol(~badAV);
Tc.Angle  = ang;
Tc.Volume = vol;

% Grouping by base name
bases = unique(Tc.Base);

% based on user input above, decide if splitting plots by unique entry base names or plotting all together in one plot
if forceSingleFigure
    basesToPlot = "All series";
else
    basesToPlot = bases;
end


%% IV. Plot
% Plots Angle (left y-axis) and Volume (right y-axis) vs time (hpf) with markers and lines connecting points in temporal order.
for b = 1:numel(basesToPlot)

    if forceSingleFigure
        idx = true(height(Tc),1);
        plotTitleBase = sprintf('All series (%s)', fileName);
    else
        thisBase = basesToPlot(b);
        idx = (Tc.Base == thisBase);
        plotTitleBase = char(thisBase);
    end
    
    % pull data for either all bases or the unique bases depending if split figures or not
    t   = Tc.Time_hpf(idx);
    ang = Tc.Angle(idx);
    vol = Tc.Volume(idx);

    % Sort data by timepoint
    [t, order] = sort(t);
    ang = ang(order);
    vol = vol(order);

    % Initialize figure
    figure('Color','w');
    hold on;

    % Define left axis: Angle
    yyaxis left
    h1 = plot(t, ang, '-o', ...
        'LineWidth', 5, ...
        'MarkerSize', 30, ...
        'Color', angleColor, ...
        'MarkerFaceColor', angleColor, ...
        'MarkerEdgeColor', angleColor);
    ylabel('Angle (degrees)');
    ylim([80,200])
    yticks(100 : 20 : 180);
    ax.LineWidth = 5;

    % Optional: make left y-axis tick labels match angle color defined above
    ax = gca;
    ax.YAxis(1).Color = angleColor;

    % Define right axis: Volume
    yyaxis right
    h2 = plot(t, vol, '-o', ...
        'LineWidth', 5, ...
        'MarkerSize', 30, ...
        'Color', volumeColor, ...
        'MarkerFaceColor', volumeColor, ...
        'MarkerEdgeColor', volumeColor);
    ylabel('Volume (pL)');
    ylim([400,1010])
    %yticks(100 : 20 : 200);
    ax.LineWidth = 5;

    % Optional: make right y-axis tick labels match volume color defined above
    ax = gca;
    ax.YAxis(2).Color = volumeColor;

    % X axis: absolute timepoint in hpf
    xlabel('Time (hpf)');
    xlim([min(t-0.5),max(t+0.5)])

    % define plot title
    title(sprintf('Angle and Volume vs Time (hpf): %s', plotTitleBase), 'Interpreter', 'none');
    
    % figure formatting settings
    grid off; box on;
    set(ax,'FontSize',90) %and other properties
    ax.LineWidth = 5;

    % Legend that works with yyaxis to define both colors
    legend([h1 h2], {'Angle','Volume'}, 'Location','best');

    % save figures as matlab figure files for the user to choose output settings in the same folder as the csv file
    figureName = sprintf('figure_%s.fig', plotTitleBase);
    figurePath = fullfile(filePath, figureName);
    saveas(gcf, figurePath);
    % can choose to have figures saved as other formats besides .fig in script, but I prefer manually saving figures using Matlab figure export settings

    hold off;
end

% print diagnostics once finished plotting
fprintf('\nTime conversion used:\n');
fprintf('  T_ref   = %g\n', T_ref);
fprintf('  hpf_ref = %g\n', hpf_ref);
fprintf('  dt_hpf  = %g hpf per index\n\n', dt_hpf);


%% V. Define functions for script
function colName = findColumn(T, candidates) % find column names flexible to capitalization
    vars = string(T.Properties.VariableNames);
    varsNorm = normalizeNames(vars);
    candNorm = normalizeNames(string(candidates));

    % 1) Exact match (case-insensitive)
    for i = 1:numel(candNorm)
        hit = find(strcmpi(varsNorm, candNorm(i)), 1);
        if ~isempty(hit)
            colName = vars(hit);
            return;
        end
    end

    % 2) Fallback: contains match
    for i = 1:numel(candNorm)
        hit = find(contains(lower(varsNorm), lower(candNorm(i))), 1);
        if ~isempty(hit)
            colName = vars(hit);
            return;
        end
    end

    error("Could not find a column matching any of: %s\nColumns MATLAB sees: %s", ...
          strjoin(candidates, ", "), strjoin(vars, ", "));
end

function out = normalizeNames(in) % find the name without any zero break spacers or spaces, tabs, or new lines or quotes
    out = strtrim(in);

    % Remove UTF-8 BOM if present (common on first column)
    out = replace(out, char(65279), "");

    % Collapse whitespace
    out = regexprep(out, '\s+', ' ');

    % Remove surrounding quotes if present
    out = strip(out, '"');
end