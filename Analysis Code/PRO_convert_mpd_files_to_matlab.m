
%% CONVERT MPD FILES TO MATLAB

% note to self: stop applying the zenith limits in the MPD conversion, and
% intead apply them in the HWD conversion I think.


%% SITE AND TIME ==========================================================

site = 'ascension';
timerange =  datenum('01-Jan-2001'):datenum('31-Dec-2012');

%% GET DIRECS =============================================================

localpath = '/Users/neil/data/MeteorRadar/'; % either mac or scratch
% localpath = '/Volumes/NJMitchell-Scratch/Data/MeteorRadar/MPD/'; % either mac or scratch
% localpath = '/Volumes/SDRed/data/MeteorRadar/';

direc = [localpath site '/mpd/'];
savedirec = [localpath site '/matlab/'];

if ~exist(savedirec,'dir') % isfolder
    mkdir(savedirec);
end

%% DAILY MPD FILES 2 MATLAB ===============================================

tic

nmets = nan(size(timerange));

for i = 1:length(timerange)
    
    timestr = datestr(timerange(i),'yyyymmdd');
    
    filestr = [direc 'mp' timestr '.' site '.mpd'];
    
    disp(['Converting ' filestr ' to matlab format...'])
    
    try
        MPD = nph_mpd2mat(filestr);
        save([savedirec datestr(MPD.Date,'yyyymmdd') '_' site '_mpd.mat'],'MPD')
    catch err
        disp(['Can''t convert ' timestr ': ' err.identifier])
        continue
    end
    
    timesofar = toc;
    totaltime = timesofar / (i/length(timerange));
    eta = (1 - i/length(timerange)) .* totaltime;
    
    disp([num2str(round(eta/60)) 'mins remaining...'])
    
end

disp('Done!')




return

























