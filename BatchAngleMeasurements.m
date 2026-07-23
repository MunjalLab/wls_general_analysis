%% Batch Measuring Pillar Angles from 3D Coordinates Chosen in Mastodon ('spots')
% Kira L. Heikes 20251204
% Made with the help of Copilot GPT-5
% For use in manuscript: Hydrostatic pressure shapes and canalizes
% semicircular canal morphology to ensure vestibular function.
% doi.org/10.64898/2026.07.22.740136


%% I. USER INPUT FIELDS
fileDir = "path\to\your\csv\files"; % USER INPUT: please enter the path to 
% your files starting with drive letter - should contain your Mastodon csv 
% files with 3d coordinate data and no other csv files
saveSuffix = "NameForOutputWithoutCSVsuffix"; % USER INPUT: please enter 
% the name you'd like the computed angles saved in - this will be used to 
% create separate files for anterior, posterior, and ventral files, with a
% prefix indicating which of these three added to each file


%% II. Pull files and initialize tables
% Preparing directory from which to extract csv files
fileNames = fullfile(fileDir, '*.csv'); % builds structure for full paths
% to the csv files in the directory
csvFiles = dir(fileNames); % list of all csv files in the directory

% define headers and data types for  empty tables in which to store data
varNames = ["Timepoint","Angle","File"]; % headers for table
varTypes = {'double', 'double','string'}; % types of data to populate

% initiate tables to populate with anterior, posterior, and ventral angles
anteriorTable=table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);
posteriorTable=table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);
ventralTable=table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);

saveName='l'; % initiate saveName as string variable type


%% III. Compute angles with AngleMeasurement.m script
% loop through all files to compute angles from 3d coordinates in each row
for k = 1:length(csvFiles)
    baseFileName = csvFiles(k).name; % pull csv file name each iteration
    fullFileName = fullfile(fileDir, baseFileName); % define full directory
    [~, saveName, ~] = fileparts(baseFileName); % define row label

    % Run AngleMeasurement function script 
    [anteriorEntry,posteriorEntry,ventralEntry] = AngleMeasurement(fullFileName); % populate outputs into anterior, posterior, and ventral data tables respectively
    anteriorEntry.File = repmat({saveName},height(anteriorEntry),1); % supply saveName as entry(ies) to File column of anterior data table based on number of angles calculated from this csv file
    posteriorEntry.File = repmat({saveName},height(posteriorEntry),1); % supply saveName as entry(ies) to File column of posterior data table based on number of angles calculated from this csv file
    ventralEntry.File = repmat({saveName},height(ventralEntry),1); % supply saveName as entry(ies) to File column of ventral data table based on number of angles calculated from this csv file

    if ~isempty(anteriorEntry), anteriorTable(k,:) = anteriorEntry; end % check if anterior data is not empty, then append to anteriorTable of data from all csv files
    if ~isempty(posteriorEntry), posteriorTable(k,:) = posteriorEntry; end % check if posterior data is not empty, then append to posteriorTable of data from all csv files
    if ~isempty(ventralEntry), ventralTable(k,:) = ventralEntry; end % check if ventral data is not empty, then append to ventralTable of data from all csv files

end


%% IV. Save computed angles to csv files in file directory
% for each anterior, posterior, and ventral tables, save to csv output with
% name defined by user above
if ~isempty(anteriorTable) % check if anteriorTable is not empty
    antFile = ['anteriorAngles_', saveSuffix, '.csv']; % define anterior csv file save name
    fullAntPath = fullfile(fileDir, antFile); % define anterior save file with directory
    writetable(anteriorTable,fullAntPath) % save data to csv
end

if ~isempty(posteriorTable) % check if posteriorTable is not empty
    postFile = ['posteriorAngles_', saveSuffix, '.csv']; % define posterior csv file save name
    fullPostPath = fullfile(fileDir, postFile); % define posterior save file with directory
    writetable(posteriorTable,fullPostPath) % save data to csv
end

if ~isempty(ventralTable) % check if ventralTable is not empty
    ventFile = ['ventralAngles_', saveSuffix, '.csv']; % define anterior csv file save name
    fullVentPath = fullfile(fileDir, ventFile); % define ventral save file with directory
    writetable(ventralTable,fullVentPath) % save data to csv
end

