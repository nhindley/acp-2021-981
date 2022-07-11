

% nph_mpd2mat.m, 20171114 - NPH University of Bath.
%
% Convert .mpd files into a MATLAB structure.
%
% Method is similar to the method of Mossy/Robin/Dave's MPD_Read.m, but
% subtley improved in parts. The output is better formatted and much
% easier to understand, as each variable is clearly named in its own field.
%
% Further, I also export all the metadata contained in the header of the
% file, such as site, location, PRF, etc. and the number of usable meteors.
%
%
% filename = '/Users/neil/data/MeteorRadar/KingEdwardPoint/mpd/mp20160708.kep.mpd';
%
%

function MPD = nph_mpd2mat(filename)

%==========================================================================
%% INITIALISE

MPD = struct;

%% READ FILE ==============================================================
% Open and get each line of data
% doing it this way as textscan can't cope with the metadata and data
% together. This is just easier. Might be faster to split, but it's hard to
% tell when the metadata ends in advance.

fileID = fopen(filename,'rt');
d = string;
n = 0;

while ~feof(fileID)
    n = n + 1;
    d(n) = fgetl(fileID);
end

% CLOSE FILE (again)
fclose(fileID);

%==========================================================================
%% Get headers and metadata
% first find where the meta data ends by finding the start of the column
% headers for the tabulated data:
for i = 1:length(d)
    if any(regexpi(d(i),'Date')) && any(regexpi(d(i),'Time')) && any(regexpi(d(i),'File'))
        lengthofmeta = i-1;
        break
    end
end

% Split Meta and Data:
meta = d(1:lengthofmeta);
d = d(lengthofmeta+2:end);

% SORT OUT METADATA =======================================================
Meta = struct;

for i = 1:length(meta)
    
    metaline = strsplit(meta(i),' '); % split each line at spaces
    
    % get Variable name:
    varname = lower(char(metaline{1}));
    
    % capitalise metadata name for neatness :)
    varname(1) = upper(varname(1));
    capitals = regexp(varname,'[^a-z ^A-Z]')+1;
    varname(capitals(capitals < length(varname))) = upper(varname(capitals(capitals < length(varname))));
    varname = varname(regexp(varname,'[a-z A-Z]'));
    
    % get metadata value
    val = metaline(2:end);
    
    % assign
    switch lower(varname)
        case 'location' % specific for location
            Meta.(varname) = str2double(strsplit(val,','));
        otherwise % all other variables:
            
            % try to convert to numerics:
            if all(~isnan(str2double(val)))
                Meta.(varname) = str2double(val);
            else % if not, leave as cell object:
                if length(val) == 1 % must just be one string:
                    Meta.(varname) = char(val);
                else % must be mixture of numbers and strings:
                    Meta.(varname) = cellstr(val);
                end
            end
    end
    
end

% IMPORTANT: The location field seems to be wrong, in that the longitude
% doesn't seem to be negative (actually this may only be the South Georgia
% radar at KEP), so fix:
if strcmpi(Meta.Sitename,'kep')
    Meta.Location(2) = -Meta.Location(2);
end
% could just use the GPS status field, but not sure if all files will
% always have that?


%==========================================================================
%% NOW THE REST OF THE TABULATED DATA:
% extract the data matrix (ignore column titles after metadata):

try
    d = cellstr(squeeze(split(d)));
catch
    % fix the weird problem where the range column in the mpd file is over 4
    % digits, which means there aren't the same number of columns.
    dc = cell(length(d),18);
    goodinds = ones(length(d),1);
    for i = 1:length(d)
        strcell = strsplit(d{i},' ');
        if length(strcell) == 18
            dc(i,:) = strcell;
            %             numgood = numgood+1;
            %             dc = [dc ; strcell];
            %             if i == 5
            %                 disp('')
            %             end
        else % and remove the end row
            goodinds(i) = 0;
            %             dc = dc(1:end-1,:);
        end
    end
    dc = dc(goodinds == 1,:);
    d = dc;
end

% first column should be empty i think:
if isempty(d{1,1})
    d = d(:,2:end);
end


%% Assign Variables =======================================================

% Site
MPD.Site = lower(Meta.Sitename);

% Date
MPD.Date      = datestr(d(1,1),'dd-mmm-yyyy');

% Assign MetaData
MPD.Meta = Meta;

% Initialise numeric matrix:
MPD.Data = struct;

% Basic File Variables:
MPD.Data.Time               = datenum(d(:,1),'yyyy/mm/dd') + rem(datenum(d(:,2),'HH:MM:SS.FFF'),1);% Datenum
MPD.Data.Range              = str2double(d(:,4)); % Range (km) LOS distance
MPD.Data.RadialVelocity     = str2double(d(:,6)); % Vrad (m/s). +ve is away from, -ve is towards.
MPD.Data.ZenithAngle        = str2double(d(:,8)); % Zenith angle from vertical (deg)
MPD.Data.AngleNorthFromEast = str2double(d(:,9)); % NORTHWARDS FROM EAST!!!! (deg) It's in the mpd files, useless as it is...
MPD.Data.Azimuth            = wrapTo360(360-(MPD.Data.AngleNorthFromEast-90)); % Actual Azimuth
MPD.Data.MaxAmplitude       = str2double(d(:,14)); % amax Amplitude of meteor echo (digiser units)
MPD.Data.Ambiguous          = str2double(d(:,10)); % Ambiguous
MPD.Data.Tau                = str2double(d(:,15)); % Tau (decay time to half of peak amplitude)
MPD.Data.MaxPhaseError      = str2double(d(:,11)); % max phase error (?!)

% Lat Lon Alt:
% NOTE: I've checked, and it seems the height in the MPD file already has a curved
% earth correction in it. So it's the local Altitude.
% method: use geometry to get the minor arc length to the point on the
% ground directly beneath the meteor trail, then use reckon:
Re = 6371; % spherical earth
arc_angle = asind((MPD.Data.Range./(str2double(d(:,5))+Re)) .* sind(180-MPD.Data.ZenithAngle));
arclen = km2deg(Re * 2*pi*(arc_angle/360));
[MPD.Data.Lat,MPD.Data.Lon] = reckon(MPD.Meta.Location(1),MPD.Meta.Location(2),arclen,MPD.Data.Azimuth);
MPD.Data.Alt                = sqrt(MPD.Data.Range.^2 + Re.^2 - (2.*Re.*MPD.Data.Range.*cosd(180-MPD.Data.ZenithAngle))) - Re;

ground_range = deg2km(arclen);
MPD.Data.x = ground_range .* sind(MPD.Data.Azimuth);
MPD.Data.y = ground_range .* cosd(MPD.Data.Azimuth);

% Local Wind:
% need to include a curved earth correction. This gives the local
% perpendicular wind speed (u) at the point where the meteor trail was seen.
% method: find the horizonatal wind component wrt the radar then, using
% geometry, find the projection of the local horz wind component wrt where
% the meteor was seen. Difficult to explain without a diagram, but it
% collapses quite nicely:
% EDIT: 20180502 - Seems like I might have been projecting to the wrong
% component. The HorzVelocity should always be bigger than the
% RadialVelocity, so you need to do a /cos():
MPD.Data.HorzVelocity = MPD.Data.RadialVelocity ./ cosd(90 - MPD.Data.ZenithAngle + arc_angle);
% MPD.Data.HorzVelocity = MPD.Data.RadialVelocity .* (sind(MPD.Data.ZenithAngle).^2);
% MPD.Data.VertVelocity = MPD.Data.RadialVelocity .* sind(MPD.Data.ZenithAngle) .* cosd(MPD.Data.ZenithAngle);



%==========================================================================
%% Remove outliers and other corrections:
% okay I admit this would have been cleaner with the previous big data
% array method in mossy/robin/dave's code. But only here! Not for working
% with later!

% removals = 'yes';
% 
% switch removals
%     case 'yes'
        
        totalmeteors = length(MPD.Data.Time);
        inds2keep = 1:totalmeteors;
        
        % Take only unambiguous meteors:
        ambvalue = 1;
        ambinds = find(MPD.Data.Ambiguous == ambvalue);
        inds2keep = intersect(inds2keep,ambinds);
        
        % Take only reasonable zenith limits:
%         zenithlimits = [15 65]; % or try [15 60], mossy used [15 75]
        zenithlimits = [0 90];
        zeninds = find(MPD.Data.ZenithAngle >= zenithlimits(1) & MPD.Data.ZenithAngle <= zenithlimits(2));
        inds2keep = intersect(inds2keep,zeninds);
        
        % Exclude Tau noise spike:
        taulimits = [0.015 1];
        tauinds = find(MPD.Data.Tau >= taulimits(1) & MPD.Data.Tau < taulimits(2));
        inds2keep = intersect(inds2keep,tauinds);
        
        % Choose only velocities < 200ms-1
        vellimits = [0 200];
        velinds = find(abs(MPD.Data.RadialVelocity) >= vellimits(1) & abs(MPD.Data.RadialVelocity) <= vellimits(2));
        inds2keep = intersect(inds2keep,velinds);
        
        
        % % % % Take only sequential meteors that occured at different times:
        % % % inds2keep = intersect(inds2keep,find(diff(MPD.Data.Time) == 0));
        % % % % note: it seems sometimes you get two meteors occurring at exactly the
        % % % % same time with the same MaxAmplitude, RadialVelocity etc...
        % % %
        % % % % Take only sequential meteors with different MaxAmplitudes:
        % % % inds2keep = intersect(inds2keep,find(diff(MPD.Data.MaxAmplitude) == 0));
        
        % number of meteors:
        MPD.Meteors = length(inds2keep);
        
        % Removal Criteria:
        MPD.RemovalCriteria = struct;
        MPD.RemovalCriteria.Ambiguous = ambvalue;
        MPD.RemovalCriteria.ZenithLimits = zenithlimits;
        MPD.RemovalCriteria.TauLimits = taulimits;
        MPD.RemovalCriteria.RadialVelocityLimits = vellimits;
        
        % Count up total meteors remaining:
        MPD.RemovedMeteors = struct;
        MPD.RemovedMeteors.Total = totalmeteors - length(inds2keep);
        MPD.RemovedMeteors.Ambiguous = totalmeteors - length(ambinds);
        MPD.RemovedMeteors.ZenithLimits = totalmeteors - length(zeninds);
        MPD.RemovedMeteors.TauLimits = totalmeteors - length(tauinds);
        MPD.RemovedMeteors.RadialVelocityLimits = totalmeteors - length(velinds);
        
        % compile:
        inds2keep = unique(inds2keep);
        f = fieldnames(MPD.Data);
        for n = 1:length(f)
            MPD.Data.(f{n}) = MPD.Data.(f{n})(inds2keep);
        end
        
%     otherwise
%         
%         MPD.Meteors = length(MPD.Data.Time);
%         MPD.RemovedMeteors = 0;
%         
% end

% Date converted:
MPD.CreatedOn = datestr(today);


% Remove unwanted fields:
MPD.Data = rmfield(MPD.Data,'Ambiguous');
% MPD.Data = rmfield(MPD.Data,'AngleNorthFromEast');








