function [dataset, metadata] = imarisread(filename, channels)
% imarisread loads an Imaris5.5 file
%
% SYNOPSIS  [dataset, metadata] = imarisread(filename, channel)
%
% INPUT filename: (optional) String containing the full name with path of
%                 the Imaris5.5 stack to be imported.
%                 If no argument is passed, or if it is set to [], the 
%                 user is asked to pick a file from a dialog.
%       channels: (optional) Indices of the channels to load (0-based!);
%                  - omit or set it to [] to load all channels.
%                  - set it to -1 to return only the metadata.
%
% OUTPUT  dataset : cell array with size (nTimepoints, nChannels). Each
%                   position contains a 3D stack. 
%         metadata: associated metadata
%
% Aaron Ponti, 2008/08/27
%              2015/03/23: Added channels input parameter.

%--------------------------------------------------------------------------
%
%   The contents of this file are subject to the Mozilla Public License
%   Version 1.1 (the "License"); you may not use this file except in
%   compliance with the License. You may obtain a copy of the License at 
%   http://www.mozilla.org/MPL/
% 
%   Software distributed under the License is distributed on an "AS IS"
%   basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
%   License for the specific language governing rights and limitations 
%   under the License.
% 
%   The Original Code is "Qu for MATLAB".
% 
%   The Initial Developer of the Original Code is Aaron Ponti.
%   All Rights Reserved.
%
%--------------------------------------------------------------------------

if nargin == 0
    filename = [];
    channels = [];
elseif nargin == 1
    channels = [];
elseif nargin == 2
    % Nothing to do
else
    error('0, 1, or 2 input parameters expected.');
end

% Initialize outputs
dataset  = {};
metadata = [];

% If no file was passed, ask the user to pick one
if isempty(filename)
    
    [fName, dirName] = uigetfile('*.ims', 'Select an Imaris file');
    if isequal(fName, 0)
        return;
    end
    filename = fullfile(dirName, fName);
end

% Get the info about the file
try

    % Use hdf5info to check wheteher the file is a valid HDF5 file...
    fileStructure = hdf5info(filename);
    
    % ... and a valid Imaris file
    if strcmp(fileStructure.GroupHierarchy.Attributes(2).Name, ...
            '/ImarisVersion') == 0
        error('Not a valid Imaris file.');
    end
    
catch

    return;

end

% Simplify notation - alias to the /DataSet/ResolutionLevel 0 substructure
datasetResLevel0 = fileStructure.GroupHierarchy.Groups(1).Groups(1);

% Open the file
fileID = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');

% Get the number of channels
readMetadataOnly = 0;
nChannelsInFile = numel(datasetResLevel0.Groups(1).Groups);
if isempty(channels)
    nChannels = nChannelsInFile;
elseif isequal(channels, -1)
    readMetadataOnly = 1;
    nChannels = nChannelsInFile;
else
    nChannels = numel(channels);
end

% Get the number of time points
nTimepoints = numel(datasetResLevel0.Groups);

% Open the group Image
imageID = H5G.open(fileID, '/DataSetInfo/Image/');

% Get some attributes
attributeXID = H5A.open_name(imageID, 'X');
metadata.width = toNumeric(H5A.read(attributeXID, 'H5ML_DEFAULT'));
H5A.close(attributeXID);

attributeYID = H5A.open_name(imageID, 'Y');
metadata.height = toNumeric(H5A.read(attributeYID, 'H5ML_DEFAULT'));
H5A.close(attributeYID);

attributeZID = H5A.open_name(imageID, 'Z');
metadata.nPlanes = toNumeric(H5A.read(attributeZID, 'H5ML_DEFAULT'));
H5A.close(attributeZID);

% Voxel size X
attributeExtMin0ID = H5A.open_name(imageID, 'ExtMin0');
extMinX = toNumeric(H5A.read(attributeExtMin0ID, 'H5ML_DEFAULT'));
H5A.close(attributeExtMin0ID);
attributeExtMax0ID = H5A.open_name(imageID, 'ExtMax0');
extMaxX = toNumeric(H5A.read(attributeExtMax0ID, 'H5ML_DEFAULT'));
H5A.close(attributeExtMax0ID);
metadata.voxelX = (extMaxX - extMinX) / metadata.width;

% Voxel size Y
attributeExtMin1ID = H5A.open_name(imageID, 'ExtMin1');
extMinY = toNumeric(H5A.read(attributeExtMin1ID, 'H5ML_DEFAULT'));
H5A.close(attributeExtMin1ID);
attributeExtMax1ID = H5A.open_name(imageID, 'ExtMax1');
extMaxY = toNumeric(H5A.read(attributeExtMax1ID, 'H5ML_DEFAULT'));
H5A.close(attributeExtMax1ID);
metadata.voxelY = (extMaxY - extMinY) / metadata.height;

% Voxel size Z
attributeExtMin2ID = H5A.open_name(imageID, 'ExtMin2');
extMinZ = toNumeric(H5A.read(attributeExtMin2ID, 'H5ML_DEFAULT'));
H5A.close(attributeExtMin2ID);
attributeExtMax2ID = H5A.open_name(imageID, 'ExtMax2');
extMaxZ = toNumeric(H5A.read(attributeExtMax2ID, 'H5ML_DEFAULT'));
H5A.close(attributeExtMax2ID);
metadata.voxelZ = (extMaxZ - extMinZ) / metadata.nPlanes;

% Close the group Image
H5G.close(imageID);

% Give the dataset array the correct dimensions
dataset = cell(nTimepoints, nChannels);

% Build the channels array
if isempty(channels) || channels == -1
    channelsToRead = 0 : nChannels - 1;
else
    % 0-based
    channelsToRead = channels;
end

% Store 
% Fill the metadata information (partially)
metadata.nChannelsInFile = nChannelsInFile;
metadata.channelsRead = channelsToRead;
metadata.nTimepoints = nTimepoints;
metadata.stagePosition = []; % TODO Fix
metadata.timeInterval = []; % TODO Fix
metadata.color = zeros(nChannels, 3);

nStacks = 0;
currentIndx = 0;
for n = channelsToRead

    % Keep track of index
    currentIndx = currentIndx + 1;

    % Get the channel ID
    channelStringID  = ['/DataSetInfo/Channel ', num2str(n), '/'];
    channelID        = H5G.open(fileID, channelStringID);
    
    % Try to read the channel color
    try
        attributeColorID = H5A.open_name(channelID, 'Color');
        metadata.color(currentIndx, 1 : 3) = toNumericVector(...
            H5A.read(attributeColorID, 'H5ML_DEFAULT'));
        H5A.close(attributeColorID);
    catch
        if nChannels == 1
            metadata.color(currentIndx, 1 : 3) = [1 1 1];
        else
            switch (currentIndx)
                case 1,
                    metadata.color(currentIndx, 1 : 3) = [1 0 0];
                case 2,
                    metadata.color(currentIndx, 1 : 3) = [0 1 0];
                case 3,
                    metadata.color(currentIndx, 1 : 3) = [0 0 1];
                otherwise,
                    metadata.color(currentIndx, 1 : 3) = rand(1, 3);
            end
        end
    end
    
    % Try to read the channel name
    try
        attributeNameID = H5A.open_name(channelID, 'Name');
        metadata.channelNames{currentIndx} = ...
            H5A.read(attributeNameID, 'H5ML_DEFAULT');
        H5A.close(attributeNameID);
    catch
       metadata.channelNames{currentIndx} = ''; 
    end
    
    H5G.close(channelID);
    
    % Get the list of time point names - they will have to be sorted
    timepointStrings = cell(1, nTimepoints);
    numericalIndices = zeros(1, nTimepoints);
    for m = 1 : nTimepoints

        % Timepoint string
        timepointString = datasetResLevel0.Groups(m).Name;
        index = strfind(timepointString, '/');
        if isempty(index)
            return
        end
        timepointString = timepointString(index(end) + 1 : end);
        
        timepointStrings{m} = timepointString;
        [~, ~, num] = getFilenameBody(timepointString);
        numericalIndices(m) = str2double(num);
    end
    
    % Now we sort them
    [~, ind] = sort(numericalIndices);
    timepointStrings = timepointStrings(ind);
    
    % Read the data if requested
    if readMetadataOnly == 0
        for m = 1 : nTimepoints
            
            nStacks = nStacks + 1;
            
            % Dataset string id
            datasetStringID = [...
                '/DataSet/ResolutionLevel 0/', ...
                timepointStrings{m}, '/Channel ', ...
                num2str(n), '/Data'];
            
            % Open the dataset
            datasetID = H5D.open (fileID, datasetStringID);
            
            try
                
                % Get the whole dataset at once
                data = ...
                    H5D.read(datasetID, 'H5ML_DEFAULT', 'H5S_ALL', ...
                    'H5S_ALL','H5P_DEFAULT');
                
                % Permute the dimensions
                data = permute(data, [2 1 3]);
                
                % Make sure the dimensions are correct
                data = data(...
                    1 : metadata.height, ...
                    1 : metadata.width, ...
                    1 : metadata.nPlanes);
                
                dataset{m, currentIndx} = data;
                
            catch
                
                dataset = {};
                metadata = [];
                return
                
            end
        end
    end
end

% Close the dataset
if readMetadataOnly == 0
    H5D.close(datasetID);
end

% Close the file
H5F.close(fileID);

% =========================================================================
% =========================================================================

function n = toNumeric(attr)

n = str2double(attr);

% -------------------------------------------------------------------------

function channelColor = toNumericVector(attr)

% Make sure to order the char array as a row vector
if size(attr, 1) > size(attr, 2)
    attr = attr';
end

% Find the indices of blank spaces that divide elements from each other
indices   = strfind(attr, ' ');
nElements = 1 + numel(indices); 

% If only one element (indices is empty), convert to double and return
if nElements == 1
    channelColor = str2double(attr);
    return
end

% If spaces were found, we split and convert one value after the other 
starts = [1, (indices + 1)];
ends   = [indices, numel(attr)];

channelColor = zeros(1, nElements);
for i = 1 : nElements
    part = attr(starts(i) : ends(i));
    channelColor(i) = str2double(part);
end


