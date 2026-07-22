%% Function to MeasurePillar Angles from 3D Coordinates Chosen in Mastodon ('spots')
% Kira L. Heikes 20251204
% Made with the help of Copilot GPT-5
% For use in manuscript: Hydrostatic pressure shapes and canalizes
% semicircular canal morphology to ensure vestibular function.


function [anteriorData,posteriorData,ventralData] = AngleMeasurement(filename)

    %% I. Load 3D coordinate data for calculating angles
    % read data into table
    fullTable = readtable(filename);
    
    % pull subset of data - only certain columns based on mastodon csv output of spots from table view
    dataTable = fullTable(3:end,[1,6,7,8,9]);
    
    % define variable names in data table
    varNames = ["Label","Timepoint","Xpoint","Ypoint","Zpoint"];
    dataTable.Properties.VariableNames = varNames;

    % initialize empty tables with specified variable names and types
    varTypes = {'string', 'double', 'double', 'double', 'double'};
    AbTable = table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);
    AtTable = table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);
    ALbTable = table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);
    PbTable = table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);
    PtTable = table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);
    PLbTable = table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);
    VbTable = table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);
    VtTable = table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);
    VLbTable = table('Size', [0, numel(varNames)], 'VariableNames', varNames, 'VariableTypes', varTypes);

    dataRows = height(dataTable); % Get the number of rows in the table
    

    %% II. Loop through each row of data table to extract data by point Label
    for rowIndex = 1:dataRows
        currentRow = dataTable(rowIndex, :); % Access current row of data table
        currentLabel = currentRow.Label; % Access current label and assign data to new table by label
        % populate 3d coordinates for each type of coordinate
        % assumes labeling matches the below labels to check for (must be
        % set manually by user in mastodon before exporting csv
        if currentLabel == "Ab"
            AbTable = [AbTable;currentRow];
        elseif currentLabel == "At"
            AtTable = [AtTable;currentRow];
        elseif currentLabel == "ALb"
            ALbTable = [ALbTable;currentRow];
        elseif currentLabel == "Pb"
            PbTable = [PbTable;currentRow];
        elseif currentLabel == "Pt"
            PtTable = [PtTable;currentRow];
        elseif currentLabel == "PLb"
            PLbTable = [PLbTable;currentRow];
        elseif currentLabel == "Vb"
            VbTable = [VbTable;currentRow];
        elseif currentLabel == "Vt"
            VtTable = [VtTable;currentRow];
        elseif currentLabel == "VLb"
            VLbTable = [VLbTable;currentRow];
        else
            print('No label found in row '+rowIndex+' of file '+filename+'.')
        end
    end
    % sort data in ascending order by timepoint
    AbTable = sortrows(AbTable, "Timepoint");
    AtTable = sortrows(AtTable, "Timepoint");
    ALbTable = sortrows(ALbTable, "Timepoint");
    PbTable = sortrows(PbTable, "Timepoint");
    PtTable = sortrows(PtTable, "Timepoint");
    PLbTable = sortrows(PLbTable, "Timepoint");
    VbTable = sortrows(VbTable, "Timepoint");
    VtTable = sortrows(VtTable, "Timepoint");
    VLbTable = sortrows(VLbTable, "Timepoint");
    

    %% III. Loop  through timepoints to calculate pillar angles at each time
    %% a) Anterior "A" Pillars
    numATimes = height(AbTable); % get number A timepoints
    ATimes = AbTable(:,"Timepoint"); % store timepoints in array
    anteriorAngle = zeros(numATimes,1); % pre-allocate array to store angles
    for m = 1:numATimes
        currAb = AbTable(m,:); % bud base coordinates
        currAt = AtTable(m,:); % pillar/bud contact site coordinates
        currALb = ALbTable(m,:); % lateral bud base coordinates
        xarray = [currAb.Xpoint,currAt.Xpoint,currALb.Xpoint];
        yarray = [currAb.Ypoint,currAt.Ypoint,currALb.Ypoint];
        zarray = [currAb.Zpoint,currAt.Zpoint,currALb.Zpoint];
        anteriorAngle(m) = findAngle(xarray,yarray,zarray);
    end
    
    %% b) Posterior "P" Pillars
    numPTimes = height(PbTable); % get number P timepoints
    PTimes = PbTable(:,"Timepoint"); % store timepoints in array
    posteriorAngle = zeros(numPTimes,1); % pre-allocate array to store angles
    for n = 1:numPTimes
        currPb = PbTable(n,:); % bud base coordinates
        currPt = PtTable(n,:); % pillar/bud contact site coordinates
        currPLb = PLbTable(n,:); % lateral bud base coordinates
        xarray = [currPb.Xpoint,currPt.Xpoint,currPLb.Xpoint];
        yarray = [currPb.Ypoint,currPt.Ypoint,currPLb.Ypoint];
        zarray = [currPb.Zpoint,currPt.Zpoint,currPLb.Zpoint];
        posteriorAngle(n) = findAngle(xarray,yarray,zarray);
    end

    %% c) Ventral "V" Pillars
    numVTimes = height(VbTable); % get number P timepoints
    VTimes = VbTable(:,"Timepoint"); % store timepoints in array
    ventralAngle = zeros(numVTimes,1); % pre-allocate array to store angles
    for n = 1:numVTimes
        currVb = VbTable(n,:); % bud base coordinates
        currVt = VtTable(n,:); % pillar/bud contact site coordinates
        currVLb = VLbTable(n,:); % lateral bud base coordinates
        xarray = [currVb.Xpoint,currVt.Xpoint,currVLb.Xpoint];
        yarray = [currVb.Ypoint,currVt.Ypoint,currVLb.Ypoint];
        zarray = [currVb.Zpoint,currVt.Zpoint,currVLb.Zpoint];
        ventralAngle(n) = findAngle(xarray,yarray,zarray);
    end
    

    %% IV. Save data in outputs.
    anteriorData = addvars(ATimes,anteriorAngle,'NewVariableNames', "Angle"); % save output for AngleMeasurement function
    posteriorData = addvars(PTimes,posteriorAngle,'NewVariableNames', "Angle"); % save output for AngleMeasurement function
    ventralData = addvars(VTimes,ventralAngle,'NewVariableNames',"Angle"); % save output for AngleMeasurement function
    

    %% V. Define angle computation function.
    function angleoutput = findAngle(x,y,z)
        P1 = [x(1),y(1),z(1)]; % bud base coordinates
        P2 = [x(2),y(2),z(2)]; % pillar/bud contact site coordinates
        P3 = [x(3),y(3),z(3)]; % lateral bud base coordinates
        A = P1-P2; % define vector A from middle of pillar/bud contact site to bud base
        B = P3-P2; % define vector A from middle of pillar/bud contact site to lateral bud base
        dotprod = dot(A,B); % find vector dot product
        mags = norm(A)*norm(B); % find product of vector normals
        anglerad = acos(dotprod/mags); % compute angle between vectors
        angleoutput = rad2deg(anglerad); % convert to degrees
    end

end
