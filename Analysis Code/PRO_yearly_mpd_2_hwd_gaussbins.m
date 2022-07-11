

% new version 20180318 - compute HWD only from a year's worth of MPD files.
% key aspects:
% - sine fit of HORIZONTAL component of radial winds for each height gate
% - compute u,v wind hourly using 2hrs worth of data sliding window
% - added daily solar tidal components.
% - compute GW variance each day from u,v wind interpolants.

% - NEW METHOD FOR BINNING!!!!! Use a Gaussian weighting function to focus
% in on specific heights and times from all available meteors. This
% massively increases the number of meteors available for fits, but doesn't
% give undue emphasis on ones far away.

% EDIT: Early 2019 - Winds are now subtracted from the MPD files to yield
% GW residuals in the MPD files themselves, I was tired of doing this every
% time! The MPD are then saved with the GW Residuals in, ready for
% variances.

% EDIT: The tidal removal has been long overdue a refit. The change here is
% that tides are fitted simultaneously (although it is mathematically
% equivalent to fit one then subtract if using a matrix inversion method)
% and then they are removed from the resolved wind field. Other
% improvements include now we use Gaussian weightings to fit the tides in
% both time and height.

% EDIT: May 2019 - It's time for an overhaul. When I was fitting and
% removing tides, I was noticing discontinuities between 0 and 24 hours
% from one day to the next. I think this is because I only fit the winds
% each day, so the edges are not well contrained by the next. To fix this,
% I'm going to load all MPD files for a year and fit a true sliding
% gaussian weighting function through it, using a composite day approach.
% For each day, a composite time-height bin is made using winds from the
% same time from next and previous days, and the tides are fitted to this.

% EDIT: June 2019 - I now have separated Tides and Mean Winds (anything
% greater than 36hours I think), so when you subtract these you get
% resolved gravity waves. You can then do resolved GW variance, which is
% nice!
%
% EDIT: Jan 2020 - Shaun/Nick rightly pointed out that when you use the
% Gaussian weightings, the wind measurement you get isn't exactly where you
% centred the Gaussian. So, take a weighted mean of the altitudes of all
% meteors used and subscribe this to a matrix, then interpolate onto a
% regular altitude grid afterwards. The mean height thing is done slightly
% differently for the tidal fit - means are taken of each hourly weighted
% mean height for the composite day window. So you're taking the mean of
% some weighted-mean heights. It's basically the same :)
%
% EDIT: Jan 2020 Fancy adding monthly composite winds and tides? Yeah why
% not!
%
% EDIT: Jan 2020 again... Need to think carefully if we want to continue
% doing the interpolation onto a regular height grid in-house or if we want
% to just export the weighted altitude matrix in the HWD files and leave
% it for plotting. Currently, the u and v fields have their own weighted altitudes and
% also the tidal components and GWs have them, so wouldn't be any extra
% effort. The main thing is that updating this code is becoming EXTREMELY
% complicated due to the need to reinterpolate every little thing back onto
% a regular height grid in-house.
%
% EDIT: Included meteors from next and previous years to reduce tidal
% fitting artefacts at the changes in the years.
%
% EDIT: June 2021 - Added a monthly composite record of meteor distribution
% peak height and FWHM, and also did it hourly so we can fit tides to it to
% determine tidal amplitudes :) - note: I don't think this will work
% actually, due to the diurnal cycle of temperatures being much larger than
% any tidal apmlitudes, which may be falling below the noise. Still, gonna
% keep the stuff here anyway in case it's useful.


%% SITE AND TIME ==========================================================

% sites = {'rothera-sk','ascension','riogrande','esrange','kep','bearlake'};
% sites = {'rothera-sk','kep'};
sites = {'kep'};

% years = [2001:2007 2009:2012];
years = 2016:2020;


for s = 1:length(sites)
    
    site = sites{s};
    
    
    %% DIRECTORIES ============================================================
    
    mpddirec = ['/Users/neil/data/MeteorRadar/' site '/matlab/'];
%     mpddirec = ['/Volumes/SDRed/data/MeteorRadar/' site '/matlab/'];
    
    for yr = years
        
        yrstr = num2str(yr);
        
        %% INITIALISE HWD ==========================================================
        HWD = struct;
        HWD.Site = site;
        HWD.Year = yr;
        HWD.Meta = struct;
        HWD.Data = struct;
        
        hourrange = datenum(yr,01,01):(1/24):datenum(yr,12,31)+(23/24);
        
        %% LOAD MPD FILES =========================================================
        % LOAD ALL MPD FILES FOR THIS YEAR
        
        disp(['Loading all MPDs for ' site ' ' yrstr '...'])
        
        mpdfiles = dir([mpddirec yrstr '*' site '_mpd.mat' ]);
        
        allMPD = struct;
        
        dayrange = datenum(yr,01,01):datenum(yr,12,31);
        
        %%%% EXTRA METEORS IN ADJACENT YEARS NEEDED FOR SLIDING TIDAL FIT
        % also load extra meteors from previous and next years if you can:
        extradayrange = datenum(yr-1,12,25):datenum(yr+1,01,05);
        
        for d = 1:length(extradayrange)
            
            yyyymmdd = datestr(extradayrange(d),'yyyymmdd');
            
            % pre assign length
            allMPD(d).Time = [];
            
            try
                load([mpddirec yyyymmdd '_' site '_mpd.mat'])
                
                % if it exists, subscribe data:
                % Subscribe:
                allMPD(d).Time                  = MPD.Data.Time;
                allMPD(d).Alt                   = MPD.Data.Alt;
                allMPD(d).RadialVelocity        = MPD.Data.RadialVelocity;
                allMPD(d).ZenithAngle           = MPD.Data.ZenithAngle;
                allMPD(d).Azimuth               = MPD.Data.Azimuth;
                %allMPD(d).AngleNorthFromEast    = wrapTo360(90 - MPD.Data.Azimuth);
                
            catch err
                disp(['Couldn''t load ' yyyymmdd '... ' err.identifier])
                continue
            end
        end
        
        % Concatenate useful quantities:
        zen     = cat(1,allMPD(:).ZenithAngle);
        az      = cat(1,allMPD(:).Azimuth);
        tim     = cat(1,allMPD(:).Time);
        alt     = cat(1,allMPD(:).Alt);
        
        vrad    = cat(1,allMPD(:).RadialVelocity);
        vhorz   = vrad ./ cosd(zen);
        
        
        % Apply Zenith Angle Limits:
        HWD.Thresholds.ZenithLimits = [15 65]; % use [15 65] for computing mean winds, impose stricter for GWs.
        inds = zen >= HWD.Thresholds.ZenithLimits(1) & zen <= HWD.Thresholds.ZenithLimits(2);
        
        % avoid silly horz velocities:
        inds = inds & vhorz < 200;
        
        az = az(inds);
        tim = tim(inds);
        alt = alt(inds);
        vrad = vrad(inds);
        vhorz = vhorz(inds);
        zen = zen(inds);
        
        % Add MPD meta data!
        HWD.Meta = MPD.Meta;
        
        % record a datevec too for months:
        dv = datevec(tim);
        % and hours:
        hr = 24 .* mod(tim,1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% WIND FITTING SPECIFICATIONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % ORIGINAL HEIGHT GATES:
        % [78,83 ; 83,86 ; 86,89 ; 89,92 ; 92,95 ; 95,100];
        % EDIT: estimated height gate center points as given in Mitchell et al. (2002) for
        % Esrange: [81.05 84.58 87.51 90.45 93.38 96.93];
        
        % HEIGHT AND TIME WINDOW STDs:
        std_z = 1.275; % gives approx 3km FWHM
        std_time = 0.85/24; % DAYS, gives approx 2hrs FWHM
        HWD.Thresholds.TimeWindowSTD = std_time;
        HWD.Thresholds.HeightWindowSTD = std_z;
        
        zrange = 76:1:105; % final height range for output matrix
        %         zrange = 72:1:108; % final height range for output matrix
        
        zlen = length(zrange);
        
        % NUMBER OF METEORS THRESHOLD FOR FIT:
        nmet_threshold = 20;
        HWD.Thresholds.NMeteors = nmet_threshold;
        
        HWD.Data.Time = reshape(hourrange,[24 length(dayrange)]);
        HWD.Data.Alt = zrange;
        HWD.Data.OldHeightGates = [78,83 ; 83,86 ; 86,89 ; 89,92 ; 92,95 ; 95,100];
        HWD.Data.Hour = 0:23;
        HWD.Data.Day = dayrange;
        HWD.Data.u = nan(zlen,24,length(dayrange));
        HWD.Data.v = nan(zlen,24,length(dayrange));
        % save the composite tide winds too :)
        HWD.Data.Comp.u = nan(zlen,24,length(dayrange));
        HWD.Data.Comp.v = nan(zlen,24,length(dayrange));
        HWD.Data.walt   = nan(zlen,24,length(dayrange));% weighted mean ALTITUDE of the meteors in each gaussian window.
        HWD.Data.wtime  = nan(zlen,24,length(dayrange));% weighted mean TIME of the meteors in gaussian window.
        HWD.Data.MeteorsUsedForFit = nan(zlen,24,length(dayrange));
        HWD.Data.TotalUsableMeteors = nan(1,length(dayrange));
        
        HWD.Data.VertMetDist = struct;
        HWD.Data.VertMetDist.Day    = nan(1,length(dayrange));
        HWD.Data.VertMetDist.Center = nan(1,length(dayrange));
        HWD.Data.VertMetDist.FWHM   = nan(1,length(dayrange));
        
        % Composite day winds and tides for each month:
        HWD.Data.MonthlyComp.u                      = nan(zlen,24,12);
        HWD.Data.MonthlyComp.v                      = nan(zlen,24,12);
        HWD.Data.MonthlyComp.walt                   = nan(zlen,24,12);
        HWD.Data.MonthlyComp.wtime                  = nan(zlen,24,12);
        HWD.Data.MonthlyComp.nmets                  = nan(zlen,24,12);    
        % monthly comp "mean winds" will just be winds - all tides, since
        % GWs should average out as noise.
        HWD.Data.MonthlyComp.SubVolumeGWVariance         = nan(zlen,12);
        HWD.Data.MonthlyComp.SubVolumeGWVarianceZonal    = nan(zlen,12);
        HWD.Data.MonthlyComp.SubVolumeGWVarianceMerid    = nan(zlen,12);
        % and the vertical distribution of meteors:
        hourvec = 0:0.5:23.5;
        HWD.Data.MonthlyComp.VertMetDist.Center          = nan(1,12);
        HWD.Data.MonthlyComp.VertMetDist.FWHM            = nan(1,12);       
        HWD.Data.MonthlyComp.VertMetDist.Hour            = hourvec;
        HWD.Data.MonthlyComp.VertMetDist.HourlyCenter    = nan(length(hourvec),12);
        HWD.Data.MonthlyComp.VertMetDist.HourlyFWHM      = nan(length(hourvec),12);
        
        
        HWD.CreatedOn = datestr(now);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% TIDES AND MEAN WINDS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tidalnames = {'Diurnal','Semi','Ter','Quar'};
        tidalperiods = [24 12 8 6] ./ 24;
        
%         tidal_std_z = 1.275; % gives approx 3km FWHM
        tidalnmets = 20; % number of meteors required for tidal fitting
        
        tidalcompositedaywindow = 4; % days, for the new composite day method
        %         tidal_std_time = 1; % 1 gives FWHM approx 2.355 days
        
        HWD.Data.Tides = struct;
        HWD.Data.Tides.AllTides.u = nan(length(zrange),24,length(dayrange));
        HWD.Data.Tides.AllTides.v = nan(length(zrange),24,length(dayrange));
        HWD.Data.Tides.OneDay.AllTides.u = nan(length(zrange),24,length(dayrange));
        HWD.Data.Tides.OneDay.AllTides.v = nan(length(zrange),24,length(dayrange));
        
        HWD.Data.Tides.TidalComponents = struct;
        for abc = 1:3
            HWD.Data.Tides.TidalComponents.u.(upper(alphabet(abc))) = nan(length(zrange),length(dayrange),length(tidalperiods));
            HWD.Data.Tides.TidalComponents.v.(upper(alphabet(abc))) = nan(length(zrange),length(dayrange),length(tidalperiods));
            HWD.Data.Tides.OneDay.TidalComponents.u.(upper(alphabet(abc))) = nan(length(zrange),length(dayrange),length(tidalperiods));
            HWD.Data.Tides.OneDay.TidalComponents.v.(upper(alphabet(abc))) = nan(length(zrange),length(dayrange),length(tidalperiods));
        end
        
        % True weighted average height:
        HWD.Data.Tides.walt = nan(length(zrange),length(dayrange));
        HWD.Data.Tides.OneDay.walt = nan(length(zrange),length(dayrange));
        
        HWD.Data.Tides.Meta.Periods = tidalperiods * 24;
        %HWD.Data.Tides.Meta.StdOfFitWindowInDays = tidal_std_time;
        HWD.Data.Tides.Meta.TidalCompositeDayWindow = tidalcompositedaywindow;
        %HWD.Data.Tides.Meta.NumberOfMeteorsForFit = tidalnmets;
        
        % MEAN WIND SMOOTHING
        % try smoothing up to the GW inertial period?
        ip = inertialperiod(MPD.Meta.Location(1));
        if ip < 24 % less than a day
            meanwind_std = ((ip/24)./2.355); % days
        else
            meanwind_std = (1./2.355); % FWHM = 1 day.
        end
        
        
        HWD.Data.MeanWinds.u = nan(length(zrange),24,length(dayrange));
        HWD.Data.MeanWinds.v = nan(length(zrange),24,length(dayrange));
        
        
        % MONTHLY COMPOSITE DAYS
        HWD.Data.MonthlyComp.Tides.AllTides.u           = nan(zlen,24,12);
        HWD.Data.MonthlyComp.Tides.AllTides.v           = nan(zlen,24,12);
        HWD.Data.MonthlyComp.Tides.TidalComponents.u.A  = nan(zlen,12,length(tidalperiods));
        HWD.Data.MonthlyComp.Tides.TidalComponents.u.B  = nan(zlen,12,length(tidalperiods));
        HWD.Data.MonthlyComp.Tides.TidalComponents.u.C  = nan(zlen,12,length(tidalperiods));
        HWD.Data.MonthlyComp.Tides.TidalComponents.v.A  = nan(zlen,12,length(tidalperiods));
        HWD.Data.MonthlyComp.Tides.TidalComponents.v.B  = nan(zlen,12,length(tidalperiods));
        HWD.Data.MonthlyComp.Tides.TidalComponents.v.C  = nan(zlen,12,length(tidalperiods));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% GWS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        HWD.Data.GWs = struct;
        
        HWD.Data.SubVolumeGWVariance = nan(zlen,length(dayrange));
        HWD.Data.ResolvedGWVariance  = nan(zlen,length(dayrange));
%         HWD.Data.ResolvedGWVariance.Meta = 'Add these components linearly (not in quadrature) to get total GW variance. They''re already projected into the line of sight!';
        % Note that you can't split these perturbations into zonal and
        % meridional if they are already projected into the meteor's line
        % of sight.
        % You can however compute the sliding variance of the residual wind
        % fields in the zonal and meridional direction (raw winds - tides -
        % mean winds = residual winds). Remember though this isn't
        % comparable to the sub volume GW variance as used by Mitchell and
        % Beldon (2009) and others. But could still be interesting.
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% CALCULATE GAUSSIAN WEIGHTINGS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        disp('Preparing Gaussian windows...')
        
        Inds                =  cell(size(HWD.Data.u));
        Weights_z           =  cell(size(HWD.Data.u));
        Weights_time        =  cell(size(HWD.Data.u));
        
        % To save a bit of time, pre-compute the gaussian weighting
        % functions for the height bins, they'll be the same each time:
        Gauss_z = struct;
        for z = 1:length(zrange)
            Gauss_z(z).vec = exp(- ((alt - zrange(z)).^2) ./ (2 * std_z^2));
        end
        
        
        % Run through meteors and fit winds at hourly intervals...
        
        %% %% FOR EACH DAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        prev_month = [];

        for d = 1:length(dayrange)
                        
            if ~strcmpi(datestr(dayrange(d),'mmm'),prev_month)
                prev_month = datestr(dayrange(d),'mmm');
                disp([prev_month ' ' yrstr])
            end
            
            % First, flatten everything
            center_time = dayrange(d); % MIDNIGHT on the day in question
            
            % Extract a reasonable time window - include a few hours in the previous and next days:
            % This makes it quicker to compute the time gaussians below, since
            % there are fewer meteors.
            dayinds = find(tim > (center_time-0.1) & tim < (center_time+1.1));
            % ^ 10% into the prev and next days
            
            tim_day = tim(dayinds);
            alt_day = alt(dayinds);
            
            % for each height step
            for z = 1:length(zrange)
                
                % Gaussian HEIGHT weighting
                gauss_z = Gauss_z(z).vec(dayinds);
                
                % for each hour
                for h = 0:23
                    
                    % compute Gaussian TIME weighting
                    hour_time = center_time+(h/24);
                    gauss_time = exp(- ((tim_day - hour_time).^2) ./ (2 * std_time^2));
                    
                    % combine the two weightings:
                    % w = gauss_time .* gauss_z;
                    % (no longer combining, keep separate for doing composite days etc without hourly breakdowns)
                    
                    % only choose meteors within 2STDs:
                    rng = (gauss_time .* gauss_z) > 0.05;
                    
                    % Subscribe to cell storage
                    Inds{z,h+1,d}               = dayinds(rng);
                    Weights_z{z,h+1,d}          = single(gauss_z(rng));
                    Weights_time{z,h+1,d}       = single(gauss_time(rng));
                    % singles for storage, but remember we need to double
                    % them for any weighted computation, especially datenums.
                    
                end % next hour
                
            end % next height
            
        end % next day
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% %% NOW CALCULATE HOURLY WINDS AND COMPOSITE DAILY TIDES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        disp('Calculating Hourly Winds and Fitting Tides...')
        
        % Assemble composite day temp variable
        CompDay = struct;
        
        prev_month = [];
        
        % for each day
        for d = 1:length(dayrange)
%         for d = 15

            
            % Display what month it is...
            if ~strcmpi(datestr(dayrange(d),'mmm'),prev_month)
                prev_month = datestr(dayrange(d),'mmm');
                disp([prev_month ' ' yrstr])
            end
            
            %--------------------------------------------------------------
            % Define a tidal range of days with which to build our
            % composite day for tidal fitting:
            tw = floor( tidalcompositedaywindow / 2 );
            
            if tidalcompositedaywindow ~= 1
                tidaldayrange = (d-tw):(d+tw-1);
            else
                tidaldayrange = d;
            end
            
            % Reset storage:
            CompDay.u       = nan(length(zrange),24);
            CompDay.v       = nan(length(zrange),24);
            CompDay.walt    = nan(length(zrange),24);
            CompDay.wtime   = nan(length(zrange),24);
            %--------------------------------------------------------------
            
            % for each height step
            for z = 1:length(zrange)
                
                % for each hour:
                for h = 0:23
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% HOURLY WIND FIT
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    inds_now    = double(Inds{z,h+1,d});
                    weights_now = double(Weights_z{z,h+1,d} .* Weights_time{z,h+1,d});
                    
                    az_now      = az(inds_now);
                    vhorz_now   = vhorz(inds_now);
                    alts_now    = alt(inds_now);
                    tim_now     = tim(inds_now);
                    
                    % find number of useful meteors within 1 FWHM:
                    nmets = length(inds_now);
                    HWD.Data.MeteorsUsedForFit(z,h+1,d) = nmets;
                    
                    % Choose whether to do the sine fit:
                    if nmets >= nmet_threshold
                        
                        [yfit,F] = nph_sinefit(az_now,vhorz_now,360,'weights',weights_now);
                        
                        %                         % Also Hourly STD from the wind fit:
                        %                         HWD.Data.SubVolumeGWVariance(z,h+1,d) = nanvar(vhorz_now - yfit,weights_now);
                        
                        % Subscribe the irregular weighted mean altitude location for this
                        % bunch of meteors. We can then interpolate onto a
                        % regular altitude grid later.
                        HWD.Data.walt(z,h+1,d) = wmean(alts_now,weights_now);
                        
                        % does the time weighting centre make much difference?
                        % Let's find out:
                        HWD.Data.wtime(z,h+1,d) = wmean(tim_now,weights_now);
                        
                        
                    else
                        F = [NaN NaN NaN];
                    end
                    
                    % F(1) = cos part, F(2) = sin part, F(3) = mean.
                    % so, if my definition of azimuth is correct, F(2) the sine
                    % part should be u and F(1) cos part should be v.
                    
                    % Check for silly numbers:
                    F(abs(F) > 200) = NaN;
                    
                    % Subscribe!
                    HWD.Data.u(z,h+1,d) = F(2);
                    HWD.Data.v(z,h+1,d) = F(1);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% COMPOSITE DAY TIDAL ASSEMBLY
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % EDIT: Also, now that we've arranged the tide fitting
                    % like this I think I can reliably add a gaussian
                    % weighting to the tidal window, such that we use e.g.
                    % 4 days of data but with a FWHM of 2 days.
                    
                    %%%%% DEALING WITH FIRST/LAST DAYS OF THE YEAR.
                    % we need to deal with this, it's come back to bite me
                    % in the arse. if the tidal comp day window stretches
                    % into another year, we need to include meteors from
                    % that year if we've got them. I've already prepared
                    % the big arrays of meteors above to contain 5 days of
                    % overlap each way, now we just need to include those
                    % meteors in this fit.
                    if any(~inrange(tidaldayrange,[1 length(dayrange)]))
                        
                        if tidalcompositedaywindow ~= 1
                            overlapinds = find(inrange(tim,dayrange(d)-tw,dayrange(d)+tw));
                        else
                            overlapinds = find(inrange(tim,dayrange(d),dayrange(d)+1));
                        end
                        % Now need to recalculate weightings for this
                        % height and time.
                        
                        % First, make them look like they're all on the
                        % same day:
                        tim_day = tim(overlapinds);
                        %                         tim_day = dayrange(d) + mod(tim_day,1);
                        
                        gauss_time = zeros(size(tim_day));
                        for dd = (tidaldayrange-1)
                            % compute extra composite time weightings
                            hour_time = (dayrange(d)+dd)+(h/24);
                            gauss_time = cat(2,gauss_time,exp(- ((tim_day - hour_time).^2) ./ (2 * std_time^2)));
                        end
                        gauss_time = nanmax(gauss_time,[],2);
                        
                        % compute extra composite altitude weightings
                        alt_day = alt(overlapinds);
                        gauss_z = exp(- ((alt_day - zrange(z)).^2) ./ (2 * std_z^2));
                        
                        %                         % NEW: compute daily weightings based on the
                        %                         % specified window width and std:
                        %                         tim_day = tim(overlapinds);
                        %                         gauss_day = exp(- ((tim_day - dayrange(d)).^2) ./ (2 * tidal_std_time^2));
                        
                        % combine
                        w = gauss_time .* gauss_z;% .* gauss_day;
                        
                        % again, only choose meteors within 2STDs:
                        rng = w > 0.05;
                        
                        % and use these as your composite
                        % indeces for this sliding window:
                        comp_inds       = double(overlapinds(rng));
                        comp_weights    = double(w(rng));
                        % this is INSTEAD of the indeces and weightings you
                        % computed earlier, which will only be for the days
                        % of the selected year and won't have any overlap.
                        
                        %                         if z == 17
                        %                             return
                        %                         end
                        
                    else
                        % Gather all the meteors at this hour and height in the
                        % time window specified:
                        comp_inds               = double(cat(1,Inds{z,h+1,tidaldayrange}));
                        comp_weights_z          = double(cat(1,Weights_z{z,h+1,tidaldayrange}));
                        comp_weights_time       = double(cat(1,Weights_time{z,h+1,tidaldayrange}));
                        comp_weights            = double(comp_weights_z .* comp_weights_time);
                        %                         comp_inds       = cat(1,Inds{z,h+1,[d-1 d+1]});
                        %                         comp_weights    = cat(1,Weights{z,h+1,[d-1 d+1]});
                        %                         if d == 31 && z == 17
                        %                             return
                        %                         end
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Now that we have selected all the indeces of all the
                    % meteors in our sliding window, arrange them and fit
                    % winds to this composite day:
                    
                    % Choose whether to do the sine fit:
                    if length(comp_inds) >= nmet_threshold
                        
                        [yfit,F] = nph_sinefit(az(comp_inds),vhorz(comp_inds),360,'weights',comp_weights);
                        
                        % get the composite altitude to which this tidal
                        % fit corresponded to:
                        CompDay.walt(z,h+1) = wmean(alt(comp_inds),double(comp_weights));
                        
                        % and also the weighted time IN HOURS:
                        % CompDay.wtime(z,h+1) = wmean(mod(tim(comp_inds),1).*24,comp_weights);
                        % You didn't think this through... What if you have
                        % HOURS of 23 etc from the previous day and 00, 01
                        % from the current day? Huh? You think of that? No.
                        % So now I'm gonna have to fix it.
                        % We need a weighted time method that can cope with not
                        % only the wide day windows, but also wraparound
                        % between 23/01 hours.
                        
                        % new method: calculate difference from nearest
                        % "THIS HOUR" of each day. This *should* always result in a nice
                        % neat gaussian around the selected hour.
                        diffs = tim(comp_inds)-round(tim(comp_inds)+(h/24));
                        % take the weighted mean (the mod is there to cope
                        % with the unusual scenario where you get an
                        % average of more than a day since the selected
                        % time?!?)
                        wmean_of_diffs = mod(wmean(diffs,comp_weights),1);
                        % if it's less than zero, circulate it:
                        if wmean_of_diffs < 0
                            wmean_of_diffs = 1 + wmean_of_diffs;
                        end
                        % and subscribe!
                        CompDay.wtime(z,h+1) = dayrange(d) + wmean_of_diffs;
                        
%                         if d == 3 && z == 5 && h == 2
%                             datestr(CompDay.wtime(z,h+1))
%                             return
%                         end
                        
                    else
                        F = [NaN NaN NaN];
                    end
                    
                    % Check for silly numbers for tides:
                    F(abs(F) > 200) = NaN;
                    
                    % And subscribe to today's temp composite day storage:
                    CompDay.u(z,h+1) = F(2);
                    CompDay.v(z,h+1) = F(1);
                    
                end % next hour
                
                
                % Also save the composite jus for funzies
                HWD.Data.Comp.u(:,:,d) = CompDay.u;
                HWD.Data.Comp.v(:,:,d) = CompDay.v;
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% COMPOSITE DAY DAILY TIDAL FITTING
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % We can do this in the same loop, since we have all time
                % info for a given height.
                
                % and extract a list of u and v components and their
                % times:
                uu = CompDay.u(z,:);
                vv = CompDay.v(z,:);
%                 tt = (0:23)/24;
                tt = CompDay.wtime(z,:);
                
%                 hold on;
%                 cla;
%                 plot((0:23)/24,tt,'.k');
%                 xlim([0 1])
%                 ylim([0 1])
%                 grid on;
%                 
%                 title(['d = ' num2str(d) ', z = ' num2str(z)])
%                 
%                 pause
%                 
%                 continue
                
                
                % and fit sine waves to them!
                if sum(~isnan(uu)) >= 18
                    % % % % %                     [~,zii] = min(abs(HWD.Data.Alt - 89));
                    % % % % %                     if sum(isnan(CompDay.u(zii,:))) >= 18
                    % % % % %                         % require at least 18 hours of data at z=89km
                    [~,Fu] = nph_sinefit(tt,uu,tidalperiods);
                    [~,Fv] = nph_sinefit(tt,vv,tidalperiods);
                else
                    Fu      = nan(length(tidalperiods),3);
                    Fv      = nan(length(tidalperiods),3);
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Check for silly values in the tidal fit
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % First check for large values:
                Fu(abs(Fu) > 100) = NaN;
                Fv(abs(Fv) > 100) = NaN;
                Fu(quadadd(Fu(:,1),Fu(:,2)) > 100,1:2) = NaN;
                Fv(quadadd(Fv(:,1),Fv(:,2)) > 100,1:2) = NaN;
                
                % If the radar was offline today (i.e. don't use
                % adjacent days in the sliding fit if the radar was off
                % for the day in question)
                if all(isnan(linearise(HWD.Data.u(:,:,d))))
                    Fu = nan(length(tidalperiods),3);
                    Fv = nan(length(tidalperiods),3);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % % % % %                     if any(abs(Fu(:)) > 100) || any(abs(Fv(:)) > 100)
                % % % % %                         yfitu   = nan(size(tt));
                % % % % %                         yfitv   = nan(size(tt));
                % % % % %                         Fu      = nan(length(tidalperiods),3);
                % % % % %                         Fv      = nan(length(tidalperiods),3);
                % % % % %                     end
                % % % % %                     % Now check to see if the radar is off today - to do
                % % % % %                     % this, check that we have a full 24hrs of data at
                % % % % %                     % z=89km:
                % % % % %                     if all(isnan(linearise(HWD.Data.u((zii-1):(zii+1),:,d))))
                % % % % %                         yfitu   = nan(size(tt));
                % % % % %                         yfitv   = nan(size(tt));
                % % % % %                         Fu      = nan(length(tidalperiods),3);
                % % % % %                         Fv      = nan(length(tidalperiods),3);
                % % % % %                     end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Reconstruct JUST the tides - NOT THE MEAN WINDS TOO
                ttt = repmat(tt,length(tidalperiods),1);
                tps = repmat(tidalperiods',1,24);
                
                uamp = repmat(quadadd(Fu(:,1),Fu(:,2)),1,24);
                vamp = repmat(quadadd(Fv(:,1),Fv(:,2)),1,24);
                
                uphase = repmat(atan2(Fu(:,1),Fu(:,2)),1,24);
                vphase = repmat(atan2(Fv(:,1),Fv(:,2)),1,24);
                
                reconu = nansum(uamp .* sin(((2.*pi.*ttt) ./ tps) + uphase),1);
                reconv = nansum(vamp .* sin(((2.*pi.*ttt) ./ tps) + vphase),1);
                
                % deal with all nans:
                if all(reconu == 0), reconu = nan(size(reconu)); end
                if all(reconv == 0), reconv = nan(size(reconv)); end
                
                HWD.Data.Tides.AllTides.u(z,:,d) = reconu;
                HWD.Data.Tides.AllTides.v(z,:,d) = reconv;
                
                % Apply some sensible limits of size of tides:
                HWD.Data.Tides.AllTides.u(abs(HWD.Data.Tides.AllTides.u) > 120) = NaN;
                HWD.Data.Tides.AllTides.v(abs(HWD.Data.Tides.AllTides.v) > 120) = NaN;
                
                % Store Tidal Components:
                for td = 1:length(tidalperiods)
                    for abc = 1:3
                        HWD.Data.Tides.TidalComponents.u.(upper(alphabet(abc)))(z,d,td) = Fu(td,abc);
                        HWD.Data.Tides.TidalComponents.v.(upper(alphabet(abc)))(z,d,td) = Fv(td,abc);
                    end
                end
                %
                %                     if any(abs(Fu(:)) > 100) || any(abs(Fv(:)) > 100)
                %                         disp('beep')
                %                     end
                
                % Finally, record the average weighted-mean height for
                % all of the hours in this composite day window:
                HWD.Data.Tides.walt(z,d) = nanmean(CompDay.walt(z,:));
                
            end % next height
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% VERTICAL DISTRIBUTION FWHM FIT
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Doing this daily for now:
            % find the fwhm of the roughly gaussian distribution of meteors
            % with height:
            
%             % select altitude bins:
%             bins = 70:110; % 1km bins
%             
%             % make histogram:
%             [a,b] = hist(alt(inrange(tim,dayrange(d),dayrange(d)+1)),bins);
%             
%             % fit with a 1-term Gaussian, bit slow sorry :(
%             try
%                 f = fit(b(:),a(:),'gauss1');
%             catch err
%                 % found a weird error with this during Jan 2010 for rothera.
%                 % The fit yielded a NaN despite no nans going in. Weird.
%                 f = struct;
%                 f.a1 = NaN; f.b1 = NaN; f.c1 = NaN;
%                 % oh well, just ignore it if so!
%             end
%             
%             % compute FWHM:
%             fwhm = 2 * sqrt(ln(2)) * f.c1; % in km
%             
%             % subscribe:
%             HWD.Data.VertMetDist.Bins = bins(:);
%             HWD.Data.VertMetDist.Amplitude(d) = f.a1; % number of meteors (1km bin)
%             HWD.Data.VertMetDist.Center(d) = f.b1; % which 1km bin
%             HWD.Data.VertMetDist.FWHM(d) = fwhm; % fwhm width in km
            
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           %%%% New way of doing this - you don't need to fit a gaussian,
           %%%% just take the mean and the standard deviation. These
           %%%% quantities assume a gaussian distribution anyway, so they
           %%%% should be the same (ish)
           
           % centre on midnight:
           alt_day = alt(inrange(tim,dayrange(d) + pm(0.5)));
           
           % trim to sensible altitudes:
           alt_day(~inrange(alt_day,[70 110])) = NaN;
           
           HWD.Data.VertMetDist.Day(d)      = dayrange(d);
           HWD.Data.VertMetDist.Center(d)   = nanmean(alt_day);
           HWD.Data.VertMetDist.FWHM(d)     = 2.355.*nanstd(alt_day);
             
            
        end % next day
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ONE DAY DAILY TIDAL FITTING
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        disp('Additional one-day tidal fit...')
        
        % Also so a one day tidal sine fit for comparison.
        
        for d = 1:length(dayrange)
            
            % for each height step
            for z = 1:length(zrange)
                
                % Just select a day of raw winds and fit tides to them...
                uu = HWD.Data.u(z,:,d);
                vv = HWD.Data.v(z,:,d);
                
                % exactly as before
                %             uu = CompDay.u(z,:);
                %             vv = CompDay.v(z,:);
%                 tt = (0:23)/24;
                tt = mod(HWD.Data.wtime(z,:,d),1);
                
                % and fit sine waves to them!
                if sum(~isnan(uu)) >= 18
                    [~,Fu] = nph_sinefit(tt,uu,tidalperiods);
                    [~,Fv] = nph_sinefit(tt,vv,tidalperiods);
                else
                    Fu      = nan(length(tidalperiods),3);
                    Fv      = nan(length(tidalperiods),3);
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Check for silly values in the tidal fit
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % First check for large values:
                Fu(abs(Fu) > 100) = NaN;
                Fv(abs(Fv) > 100) = NaN;
                Fu(quadadd(Fu(:,1),Fu(:,2)) > 100,1:2) = NaN;
                Fv(quadadd(Fv(:,1),Fv(:,2)) > 100,1:2) = NaN;
                
                % If the radar was offline today (i.e. don't use
                % adjacent days in the sliding fit if the radar was off
                % for the day in question)
                if all(isnan(linearise(HWD.Data.u(:,:,d))))
                    Fu = nan(length(tidalperiods),3);
                    Fv = nan(length(tidalperiods),3);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Reconstruct JUST the tides - NOT THE MEAN WINDS TOO
                ttt = repmat(tt,length(tidalperiods),1);
                tps = repmat(tidalperiods',1,24);
                
                uamp = repmat(quadadd(Fu(:,1),Fu(:,2)),1,24);
                vamp = repmat(quadadd(Fv(:,1),Fv(:,2)),1,24);
                
                uphase = repmat(atan2(Fu(:,1),Fu(:,2)),1,24);
                vphase = repmat(atan2(Fv(:,1),Fv(:,2)),1,24);
                
                reconu = nansum(uamp .* sin(((2.*pi.*ttt) ./ tps) + uphase),1);
                reconv = nansum(vamp .* sin(((2.*pi.*ttt) ./ tps) + vphase),1);
                
                % deal with all nans:
                if all(reconu == 0), reconu = nan(size(reconu)); end
                if all(reconv == 0), reconv = nan(size(reconv)); end
                
                HWD.Data.Tides.OneDay.AllTides.u(z,:,d) = reconu;
                HWD.Data.Tides.OneDay.AllTides.v(z,:,d) = reconv;
                
                % Apply some sensible limits of size of tides:
                HWD.Data.Tides.OneDay.AllTides.u(abs(HWD.Data.Tides.OneDay.AllTides.u) > 120) = NaN;
                HWD.Data.Tides.OneDay.AllTides.v(abs(HWD.Data.Tides.OneDay.AllTides.v) > 120) = NaN;
                
                % Store Tidal Components:
                for td = 1:length(tidalperiods)
                    for abc = 1:3
                        HWD.Data.Tides.OneDay.TidalComponents.u.(upper(alphabet(abc)))(z,d,td) = Fu(td,abc);
                        HWD.Data.Tides.OneDay.TidalComponents.v.(upper(alphabet(abc)))(z,d,td) = Fv(td,abc);
                    end
                end
                
            end % next height
            
            % Finally, record the average weighted-mean height for
            % all of the hours in this composite day window:
            HWD.Data.Tides.OneDay.walt(:,d) = nanmean(HWD.Data.walt(:,:,d),2);
            
        end % next day
        
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% %%% COMPOSITE MONTHLY WINDS AND TIDES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % first, compute a composite day for a whole month's data of wind
        % against height.
        % then, fit tidal periods to this for each height.
        
        disp('Computing monthly composite winds and tides...')
        
        % for each month
        for m = 1:12
            
            dayrange = daynumber(datenum(yr,m,01):datenum(yr,m+1,0));
            
            % find middle of the month:
            middle_day_of_month = floor(nanmean(tim(month(tim) == m)));
            
            % for each height step
            for z = 1:length(zrange)
                
                % for each hour:
                for h = 0:23
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%% HOURLY WIND FIT
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    inds_now = []; weights_now = [];
                    for d = dayrange
                        inds_now        = double(cat(1,inds_now,Inds{z,h+1,d}));
                        weights_now     = double(cat(1,weights_now,Weights_z{z,h+1,d} .* Weights_time{z,h+1,d}));
                    end
                    
                    az_now      = az(inds_now);
                    vhorz_now   = vhorz(inds_now);
                    alts_now    = alt(inds_now);
                    tim_now     = tim(inds_now);
                    
                    % find number of useful meteors within 1 FWHM:
                    nmets = length(inds_now);
                    
                    % Choose whether to do the sine fit:
                    if nmets >= nmet_threshold
                        
                        [yfit,F] = nph_sinefit(az_now,vhorz_now,360,'weights',weights_now);
                        
                        HWD.Data.MonthlyComp.nmets(z,h+1,m) = length(notnan(vhorz_now(weights_now > 0.5)));
                        
                        % Subscribe the irregular weighted mean altitude location for this
                        % bunch of meteors. We can then interpolate onto a
                        % regular altitude grid later.
                        HWD.Data.MonthlyComp.walt(z,h+1,m) = wmean(alts_now,weights_now);
                        
%                         % does the time weighting centre make much difference?
%                         % Let's find out:
%                         HWD.Data.MonthlyComp.wtime(z,h+1,m) = ...
%                             middle_day_of_month + wmean(mod(tim_now,1),weights_now);
%                         
                        % and also the weighted time IN HOURS:
                        % new method: calculate difference from nearest
                        % "THIS HOUR" of each day. This *should* always result in a nice
                        % neat gaussian around the selected hour.
                        diffs = tim_now-round(tim_now+(h/24));
                        % take the weighted mean (the mod is there to cope
                        % with the unusual scenario where you get an
                        % average of more than a day since the selected
                        % time?!?)
                        wmean_of_diffs = mod(wmean(diffs,weights_now),1);
                        % if it's less than zero, circulate it:
                        if wmean_of_diffs < 0
                            wmean_of_diffs = 1 + wmean_of_diffs;
                        end
                        
                        % and subscribe!                        
                        HWD.Data.MonthlyComp.wtime(z,h+1,m) = ...
                            middle_day_of_month + wmean_of_diffs;
                        
                        
                    else
                        F = [NaN NaN NaN];
                    end
                    
                    % F(1) = cos part, F(2) = sin part, F(3) = mean.
                    % so, if my definition of azimuth is correct, F(2) the sine
                    % part should be u and F(1) cos part should be v.
                    
                    % Check for silly numbers:
                    F(abs(F) > 200) = NaN;
                    
                    % Subscribe!
                    HWD.Data.MonthlyComp.u(z,h+1,m) = F(2);
                    HWD.Data.MonthlyComp.v(z,h+1,m) = F(1);
                    
                end % next hour
                
                % =================================================================
                %%%% NEXT, FIT TIDAL COMPONENTS TO THESE MONTHLY COMP WINDS:
                % =================================================================
                
                % We can do this in the same loop, since we have all time
                % info for a given height.
                
                uu = HWD.Data.MonthlyComp.u(z,:,m);
                vv = HWD.Data.MonthlyComp.v(z,:,m);
%                 tt = (0:23)/24;
                tt = HWD.Data.MonthlyComp.wtime(z,:,m);
                
                % and fit sine waves to them!
                if sum(~isnan(uu)) >= 18
                    [~,Fu] = nph_sinefit(tt,uu,tidalperiods);
                    [~,Fv] = nph_sinefit(tt,vv,tidalperiods);
                else
                    Fu      = nan(length(tidalperiods),3);
                    Fv      = nan(length(tidalperiods),3);
                end
                
                % First check for large values:
                Fu(abs(Fu) > 150) = NaN;
                Fv(abs(Fv) > 150) = NaN;
                Fu(quadadd(Fu(:,1),Fu(:,2)) > 150,1:2) = NaN;
                Fv(quadadd(Fv(:,1),Fv(:,2)) > 150,1:2) = NaN;
                
                % NOW RECONSTRUCT THE TIDES...
                ttt = repmat(tt,length(tidalperiods),1);
                tps = repmat(tidalperiods',1,24);
                
                uamp = repmat(quadadd(Fu(:,1),Fu(:,2)),1,24);
                vamp = repmat(quadadd(Fv(:,1),Fv(:,2)),1,24);
                
                uphase = repmat(atan2(Fu(:,1),Fu(:,2)),1,24);
                vphase = repmat(atan2(Fv(:,1),Fv(:,2)),1,24);
                
                reconu = nansum(uamp .* sin(((2.*pi.*ttt) ./ tps) + uphase),1);
                reconv = nansum(vamp .* sin(((2.*pi.*ttt) ./ tps) + vphase),1);
                
                % deal with all nans:
                if all(reconu == 0), reconu = nan(size(reconu)); end
                if all(reconv == 0), reconv = nan(size(reconv)); end
                
                HWD.Data.MonthlyComp.Tides.AllTides.u(z,:,m) = reconu;
                HWD.Data.MonthlyComp.Tides.AllTides.v(z,:,m) = reconv;
                
                % Apply some sensible limits of size of tides:
                HWD.Data.MonthlyComp.Tides.AllTides.u(abs(HWD.Data.MonthlyComp.Tides.AllTides.u) > 150) = NaN;
                HWD.Data.MonthlyComp.Tides.AllTides.v(abs(HWD.Data.MonthlyComp.Tides.AllTides.v) > 150) = NaN;
                
                % Store Tidal Components:
                for td = 1:length(tidalperiods)
                    for abc = 1:3
                        HWD.Data.MonthlyComp.Tides.TidalComponents.u.(upper(alphabet(abc)))(z,m,td) = Fu(td,abc);
                        HWD.Data.MonthlyComp.Tides.TidalComponents.v.(upper(alphabet(abc)))(z,m,td) = Fv(td,abc);
                    end
                end
                
                
            end % next height
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% COMPOSITE MONTHLY METEOR DISTRIBUTIONS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % record monthly composite meteor distributions so that we can
        % determine tidal temperature amplitudes etc.
        
        monthinds       = dv(:,1) == yr & dv(:,2) == m;
        alt_month       = alt(monthinds);
        hour_month      = hr(monthinds);
        
        % trim to sensible altitudes:
        alt_month(~inrange(alt_month,[70 110])) = NaN;
        
        % take the avg center and FWHM for the whole month:
        HWD.Data.MonthlyComp.VertMetDist.Center(m) = mean(alt_month(:),'omitnan');
        HWD.Data.MonthlyComp.VertMetDist.FWHM(m)   = 2.355.*std(alt_month(:),'omitnan');
        
        % now for each hour, let's take a weighted mean and std of the
        % meteor distribution:
        Gauss_h = struct;
        for h = 1:length(hourvec)
            w_hour = exp(- ((hour_month - hourvec(h)).^2) ./ (2 * std_time^2));
            inds = w_hour > 0.05;
            HWD.Data.MonthlyComp.VertMetDist.HourlyCenter(h,m) = wmean(alt_month(inds),w_hour(inds));
            HWD.Data.MonthlyComp.VertMetDist.HourlyFWHM(h,m)   = 2.355.*std(linearise(alt_month(inds)),linearise(w_hour(inds)),'omitnan');
        end
            
        end % next month
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% % RESHAPE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        sz = size(HWD.Data.u);
        
        HWD.Data.u      = reshape(HWD.Data.u,[sz(1) sz(2)*sz(3)]);
        HWD.Data.v      = reshape(HWD.Data.v,[sz(1) sz(2)*sz(3)]);
        HWD.Data.walt   = reshape(HWD.Data.walt,[sz(1) sz(2)*sz(3)]);
        HWD.Data.wtime  = reshape(HWD.Data.wtime,[sz(1) sz(2)*sz(3)]);
        HWD.Data.Time   = reshape(HWD.Data.Time,[1 sz(2)*sz(3)]);
        
        HWD.Data.Tides.AllTides.u   = reshape(HWD.Data.Tides.AllTides.u,[sz(1) sz(2)*sz(3)]);
        HWD.Data.Tides.AllTides.v   = reshape(HWD.Data.Tides.AllTides.v,[sz(1) sz(2)*sz(3)]);
        HWD.Data.Tides.OneDay.AllTides.u   = reshape(HWD.Data.Tides.OneDay.AllTides.u,[sz(1) sz(2)*sz(3)]);
        HWD.Data.Tides.OneDay.AllTides.v   = reshape(HWD.Data.Tides.OneDay.AllTides.v,[sz(1) sz(2)*sz(3)]);
        
        HWD.Data.MeteorsUsedForFit = reshape(HWD.Data.MeteorsUsedForFit,[sz(1) sz(2)*sz(3)]);
        
        % % % % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % % % %         %% %%% INTERPOLATE WEIGHTED HEIGHTS ONTO REGULAR HEIGHT GRID
        % % % % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % % % %
        % % % % %         %%%%% DO THIS *BEFORE* MEAN WINDS, PWS AND GWS!!!!!!!!!
        % % % % %         % you need the regridded heights to be subtracted correctly or
        % % % % %         % you'll introduce signals due to the match up not being right.
        % % % % %
        % % % % %         % as we know, our weighted Gaussian method doesn't always return
        % % % % %         % the wind at the height/time at the centre of the Gaussian, so we
        % % % % %         % need to interpolate back from the irregular weighted-mean height
        % % % % %         % grid to the regular height grid that we specified at the start.
        % % % % %
        % % % % %         disp('Regridding gaussian-weighted heights and times to regular grids...')
        % % % % %
        % % % % %         %%%%% WINDS
        % % % % %
        % % % % %         % linearise vectors for the scattered interpolant:
        % % % % %         tlin = HWD.Data.wtime(:);
        % % % % %         zlin = HWD.Data.walt(:);
        % % % % %         nanlin = double(linearise(isnan(HWD.Data.u) | isnan(HWD.Data.v)));
        % % % % %
        % % % % %         % compute lengths for vectors
        % % % % %         timelen         = size(HWD.Data.u,2);
        % % % % %         zlen            = length(HWD.Data.Alt);
        % % % % %         indmap_z        = repmat((1:zlen)',1,timelen);
        % % % % %         indmap_time     = repmat(1:timelen,zlen,1);
        % % % % %
        % % % % %         % create a gridded interpolant using the indeces of the irregular
        % % % % %         % data, pretending it was gridded:
        % % % % %         Fz = scatteredInterpolant(tlin(nanlin == 0),zlin(nanlin == 0),indmap_z(nanlin == 0),'linear','none');
        % % % % %         indmapi_z = Fz({HWD.Data.Time,HWD.Data.Alt})';
        % % % % %         Ftime = scatteredInterpolant(tlin(nanlin == 0),zlin(nanlin == 0),indmap_time(nanlin == 0),'linear','none');
        % % % % %         indmapi_time = Ftime({HWD.Data.Time,HWD.Data.Alt})';
        % % % % %         % this way you only have to do ONE scattered interpolant.
        % % % % %
        % % % % %         % now sample the original data irregularly using the above grid:
        % % % % %         Fu = griddedInterpolant({1:timelen,1:zlen},HWD.Data.u','linear','none');
        % % % % %         ui = Fu(indmapi_time,indmapi_z);
        % % % % %         Fv = griddedInterpolant({1:timelen,1:zlen},HWD.Data.v','linear','none');
        % % % % %         vi = Fv(indmapi_time,indmapi_z);
        % % % % %
        % % % % %         % Reassign the interpolated winds:
        % % % % %         HWD.Data.u = ui;
        % % % % %         HWD.Data.v = vi;
        % % % % %         HWD.Data = rmfield(HWD.Data,{'walt','wtime'});
        % % % % %
        % % % % %         % and also the reconstructed tides
        % % % % %         Fu = griddedInterpolant({1:timelen,1:zlen},HWD.Data.Tides.AllTides.u','linear','none');
        % % % % %         HWD.Data.Tides.AllTides.u = Fu(indmapi_time,indmapi_z);
        % % % % %         Fv = griddedInterpolant({1:timelen,1:zlen},HWD.Data.Tides.AllTides.v','linear','none');
        % % % % %         HWD.Data.Tides.AllTides.v = Fv(indmapi_time,indmapi_z);
        % % % % %
        % % % % %         %%% AND THE COMPOSITE MONTHLY WINDS AND TIDAL COMPONENTS...
        % % % % %         % just use basic interps for now, it's only small.
        % % % % %         for m = 1:12 % months
        % % % % %             for h = 0:23 % hours of composite day
        % % % % %
        % % % % %                 % WINDS
        % % % % %                 uu = HWD.Data.MonthlyComp.u(:,h+1,m);
        % % % % %                 vv = HWD.Data.MonthlyComp.v(:,h+1,m);
        % % % % %                 zz = HWD.Data.MonthlyComp.walt(:,h+1,m);
        % % % % %
        % % % % %                 nanlocs = isnan(uu) | isnan(vv) | isnan(zz);
        % % % % %
        % % % % %                 uui = interp1(zz(~nanlocs),uu(~nanlocs),HWD.Data.Alt);
        % % % % %                 vvi = interp1(zz(~nanlocs),vv(~nanlocs),HWD.Data.Alt);
        % % % % %
        % % % % %                 HWD.Data.MonthlyComp.u(:,h+1,m) = uui;
        % % % % %                 HWD.Data.MonthlyComp.v(:,h+1,m) = vvi;
        % % % % %
        % % % % %                 % and all tides:
        % % % % %                 nanlocs = isnan(zz);
        % % % % %                 uu = HWD.Data.MonthlyComp.Tides.AllTides.u(:,h+1,m);
        % % % % %                 vv = HWD.Data.MonthlyComp.Tides.AllTides.v(:,h+1,m);
        % % % % %                 HWD.Data.MonthlyComp.Tides.AllTides.u(:,h+1,m) = interp1(zz(~nanlocs),uu(~nanlocs),HWD.Data.Alt);
        % % % % %                 HWD.Data.MonthlyComp.Tides.AllTides.v(:,h+1,m) = interp1(zz(~nanlocs),vv(~nanlocs),HWD.Data.Alt);
        % % % % %
        % % % % %             end % next hour
        % % % % %
        % % % % %             % once all hours are done for this height,
        % % % % %             % do TIDES...
        % % % % %             uv = {'u','v'};
        % % % % %             nanlocs = isnan(zz);
        % % % % %             for tide = 1:length(tidalperiods) % for each tide
        % % % % %                 for w = 1:2 % for u and v
        % % % % %                     A = HWD.Data.MonthlyComp.Tides.TidalComponents.(uv{w}).A(:,m,tide);
        % % % % %                     B = HWD.Data.MonthlyComp.Tides.TidalComponents.(uv{w}).B(:,m,tide);
        % % % % %                     C = HWD.Data.MonthlyComp.Tides.TidalComponents.(uv{w}).C(:,m,tide);
        % % % % %
        % % % % %                     HWD.Data.MonthlyComp.Tides.TidalComponents.(uv{w}).A(:,m,tide) = interp1(zz(~nanlocs),A(~nanlocs),HWD.Data.Alt);
        % % % % %                     HWD.Data.MonthlyComp.Tides.TidalComponents.(uv{w}).B(:,m,tide) = interp1(zz(~nanlocs),B(~nanlocs),HWD.Data.Alt);
        % % % % %                     HWD.Data.MonthlyComp.Tides.TidalComponents.(uv{w}).C(:,m,tide) = interp1(zz(~nanlocs),B(~nanlocs),HWD.Data.Alt);
        % % % % %                 end
        % % % % %             end
        % % % % %
        % % % % %         end
        % % % % %
        % % % % %         %% %%% AND TIDAL COMPONENTS...
        % % % % %
        % % % % %         % linearise vectors for the scattered interpolant:
        % % % % %         tlin = linearise(repmat(HWD.Data.Day,zlen,1));
        % % % % %         zlin = HWD.Data.Tides.walt(:);
        % % % % %         nanlin = double(isnan(HWD.Data.Tides.walt(:)));
        % % % % %
        % % % % %         % compute lengths for vectors
        % % % % %         timelen         = length(HWD.Data.Day);
        % % % % %         zlen            = length(HWD.Data.Alt);
        % % % % %         indmap_z        = repmat((1:zlen)',1,timelen);
        % % % % %         indmap_time     = repmat(1:timelen,zlen,1);
        % % % % %
        % % % % %         % create a gridded interpolant using the indeces of the irregular
        % % % % %         % data, pretending it was gridded:
        % % % % %         Fz = scatteredInterpolant(tlin(nanlin == 0),zlin(nanlin == 0),indmap_z(nanlin == 0),'linear','none');
        % % % % %         indmapi_z = Fz({HWD.Data.Day,HWD.Data.Alt})';
        % % % % %         Ftime = scatteredInterpolant(tlin(nanlin == 0),zlin(nanlin == 0),indmap_time(nanlin == 0),'linear','none');
        % % % % %         indmapi_time = Ftime({HWD.Data.Day,HWD.Data.Alt})';
        % % % % %         % this way you only have to do ONE scattered interpolant.
        % % % % %
        % % % % %         % now sample the original data irregularly using the above grid:
        % % % % %         % now for each TIDAL COMPONENT:
        % % % % %         parts = {'A','B','C'};
        % % % % %         ui = nan(zlen,timelen,length(tidalperiods));
        % % % % %         vi = nan(zlen,timelen,length(tidalperiods));
        % % % % %         for p = 1:3
        % % % % %             for tid = 1:length(tidalperiods)
        % % % % %
        % % % % %                 % extract the tidal component:
        % % % % %                 FARTu.(parts{p}) = HWD.Data.Tides.TidalComponents.u.(parts{p})(:,:,tid);
        % % % % %                 FARTv.(parts{p}) = HWD.Data.Tides.TidalComponents.v.(parts{p})(:,:,tid);
        % % % % %
        % % % % %                 Fu = griddedInterpolant({1:timelen,1:zlen},FARTu.(parts{p})','linear','none');
        % % % % %                 ui(:,:,tid) = Fu(indmapi_time,indmapi_z);
        % % % % %
        % % % % %                 Fv = griddedInterpolant({1:timelen,1:zlen},FARTv.(parts{p})','linear','none');
        % % % % %                 vi(:,:,tid) = Fv(indmapi_time,indmapi_z);
        % % % % %
        % % % % %             end
        % % % % %
        % % % % %             % and assign:
        % % % % %             HWD.Data.Tides.TidalComponents.u.(parts{p}) = ui;
        % % % % %             HWD.Data.Tides.TidalComponents.v.(parts{p}) = vi;
        % % % % %
        % % % % %         end
        % % % % %
        % % % % %         HWD.Data.Tides = rmfield(HWD.Data.Tides,'WeightedAltitude');
        % % % % %
        
        
        
        
        % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % % % %
        % % % % % %     new approach
        % % % % % % we have a certain number of grids in the HWD file. Some are 30x8760
        % % % % % % grids, some are 30x365 grids etc.
        % % % % % % create a scattered inteprolant for each, then change the points to be all
        % % % % % % the data you need.
        % % % % %
        % % % % % return
        % % % % %
        % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % % % % % FIRST GRID: 30x8760 or similar (hourly data sampling)
        % % % % % G1.walt = HWD.Data.walt;
        % % % % % G1.wtime = HWD.Data.wtime;
        % % % % % % fix the end of time:
        % % % % % G1.wtime(:,end) = HWD.Data.Time(:,end);
        % % % % % % find any nans
        % % % % % G1.nanlocs = isnan(G1.walt) | isnan(G1.wtime);
        % % % % % % % % % G1.dummy_walt = repmat(nanmean(G1.walt,2),1,length(HWD.Data.Time));
        % % % % % % % % % G1.walt(isnan(G1.walt)) = G1.dummy_walt(isnan(G1.walt));
        % % % % % % % % % G1.dummy_wtime = repmat(HWD.Data.Time,zlen,1);
        % % % % % % % % % G1.wtime(isnan(G1.wtime)) = G1.dummy_wtime(isnan(G1.wtime));
        % % % % % % create interpolant
        % % % % % dummy_in = ones(size(G1.nanlocs));
        % % % % % F1 = scatteredInterpolant(linearise(G1.wtime(~G1.nanlocs)),linearise(G1.walt(~G1.nanlocs)),linearise(dummy_in(~G1.nanlocs)),'linear','none');
        % % % % %
        % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % % % % % SECOND GRID: 30x365 or similar (daily data sampling)
        % % % % % G2.walt = HWD.Data.Tides.walt;
        % % % % % G2.wtime = repmat(HWD.Data.Day,zlen,1);
        % % % % % % find any nans
        % % % % % G2.nanlocs = isnan(G2.walt) | isnan(G2.wtime);
        % % % % % % create interpolant
        % % % % % dummy_in = ones(size(G2.nanlocs));
        % % % % % F2 = scatteredInterpolant(linearise(G2.wtime(~G2.nanlocs)),linearise(G2.walt(~G2.nanlocs)),linearise(dummy_in(~G2.nanlocs)),'linear','none');
        % % % % %
        % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % % % % % NOW GO THROUGH AND REGRID EVERYTHING BASED WITH THE RIGHT GRID:
        % % % % %
        % % % % % uv = {'u','v'};
        % % % % % ABC = {'A','B','C'};
        % % % % %
        % % % % % % raw winds:
        % % % % % for w = 1:2
        % % % % %     nanlocs = isnan(HWD.Data.(uv{w})) | G1.nanlocs;
        % % % % %     F1.Values = linearise(HWD.Data.(uv{w})(~nanlocs));
        % % % % %     HWD.Data.(uv{w}) = F1({HWD.Data.Time,HWD.Data.Alt})';
        % % % % % end
        % % % % %
        % % % % % % tides:
        % % % % % for w = 1:2
        % % % % %     % use grid 1 for the "all tides"
        % % % % %     nanlocs = isnan(HWD.Data.Tides.AllTides.(uv{w})) | G1.nanlocs;
        % % % % %     F1.Values = HWD.Data.Tides.AllTides.(uv{w})(~nanlocs);
        % % % % %     HWD.Data.Tides.AllTides.(uv{w}) = F1({HWD.Data.Time,HWD.Data.Alt})';
        % % % % %
        % % % % %     % and grid 2 for the tidal components
        % % % % %     for tide = 1:4
        % % % % %         for abc = 1:3
        % % % % %             nanlocs = isnan(HWD.Data.Tides.TidalComponents.(uv{w}).(ABC{abc})) | G2.nanlocs;
        % % % % %             F1.Values = HWD.Data.Tides.TidalComponents.(uv{w}).(ABC{abc})(~nanlocs);
        % % % % %         end
        % % % % %     end
        % % % % %
        % % % % % end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% COMPUTE MEAN WINDS (which will include PWs) AND RESOLVED GWS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Great! We've now got hourly raw winds and fitted tidal
        % components. Now let's compute the mean winds, such that we can
        % determine only GW perturbations, both resolved and sub volume.
        
        % This is so that we can extract the tides, PWs and mean winds from
        % the observed radial velocities to leave all resolved and
        % sub-volume GWs.
        
        smoothing_stdvs = [0.05 24*meanwind_std];
        HWD.Data.MeanWinds.Meta.SmoothingSTDs = smoothing_stdvs;
        
        % First remove the tides:
%         smoo.u = (HWD.Data.u - HWD.Data.Tides.AllTides.u);
%         smoo.v = (HWD.Data.v - HWD.Data.Tides.AllTides.v);
        smoo.u = (HWD.Data.u - HWD.Data.Tides.OneDay.AllTides.u);
        smoo.v = (HWD.Data.v - HWD.Data.Tides.OneDay.AllTides.v);
        
        % now deal with NaNs using super smoothed data:
        nanlocs = isnan(smoo.u) | isnan(smoo.v);
        
        supersmoothu = movmean(smoo.u,(24*30)+1,2,'omitnan');
        supersmoothv = movmean(smoo.v,(24*30)+1,2,'omitnan');
        
        smoo.u(nanlocs) = supersmoothu(nanlocs);
        smoo.v(nanlocs) = supersmoothv(nanlocs);
        
        % fix any stragglers
        smoo.u(isnan(smoo.u)) = nanmean(smoo.u(:));
        smoo.v(isnan(smoo.v)) = nanmean(smoo.v(:));
        
        % apply filter to separate mean winds
        smoo.u = imgaussfilt(smoo.u,smoothing_stdvs,'padding','replicate');
        smoo.v = imgaussfilt(smoo.v,smoothing_stdvs,'padding','replicate');
        
        % and reapply NaNs:
        smoo.u(nanlocs) = NaN;
        smoo.v(nanlocs) = NaN;
        
        HWD.Data.MeanWinds.u = smoo.u;
        HWD.Data.MeanWinds.v = smoo.v;
        
        %%% AND FINALLY, RESOLVED GRAVITY WAVES IN u AND v
%         HWD.Data.GWs.u = (HWD.Data.u - HWD.Data.Tides.AllTides.u) - HWD.Data.MeanWinds.u;
%         HWD.Data.GWs.v = (HWD.Data.v - HWD.Data.Tides.AllTides.v) - HWD.Data.MeanWinds.v;
        HWD.Data.GWs.u = (HWD.Data.u - HWD.Data.Tides.OneDay.AllTides.u) - HWD.Data.MeanWinds.u;
        HWD.Data.GWs.v = (HWD.Data.v - HWD.Data.Tides.OneDay.AllTides.v) - HWD.Data.MeanWinds.v;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% SAVE GRAVITY WAVE RESIDUALS IN MPD FILE FOR EACH METEOR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Build an interpolant for the complete winds, and also for the
        % resolved gravity wave perturbations only.
        % To get the sub-volume perts, subtract all of the resolved wind
        % projected on to the time and height of each meteor.
        % For the resolved GW perts, project the resolved GW winds onto the
        % time height of each metoer.
        % Then you can get the variance of both in the same way.
        
        disp('Computing daily (and monthly composite daily) GW variances...')
        
        F = struct;
        
        F.u = griddedInterpolant({HWD.Data.Alt,HWD.Data.Time},HWD.Data.u,'linear','none');
        F.v = griddedInterpolant({HWD.Data.Alt,HWD.Data.Time},HWD.Data.v,'linear','none');
        
        F.ugw = griddedInterpolant({HWD.Data.Alt,HWD.Data.Time},HWD.Data.GWs.u,'linear','none');
        F.vgw = griddedInterpolant({HWD.Data.Alt,HWD.Data.Time},HWD.Data.GWs.v,'linear','none');
        
        % Now compute Sub-Volume and Resolved GW residuals all at once for
        % the whole year:
        uflow = F.u(alt,tim);
        vflow = F.v(alt,tim);
        flow_radu = uflow .* sind(az) .* sind(zen);
        flow_radv = vflow .* cosd(az) .* sind(zen);
        svgw = vrad - (flow_radu + flow_radv);
        
        ugw = F.ugw(alt,tim);
        vgw = F.vgw(alt,tim);
        gw_radu = ugw .* sind(az) .* sind(zen);
        gw_radv = vgw .* cosd(az) .* sind(zen);
        rgw = gw_radu + gw_radv;
%         rgwu = gw_radu;
%         rgwv = gw_radv;
        
        % Now for each day get daily variance for each height:
        dayrange = datenum(yr,01,01):datenum(yr,12,31);
        
        for d = 1:length(dayrange)
            
            % now for each height step
            for z = 1:length(zrange)
                
                % Use the weightings and indeces from above, otherwise it'll
                % take forever to run again.
                inds_now            = double(cat(1,Inds{z,:,d}));
                weights_now_z       = double(cat(1,Weights_z{z,:,d}));
                weights_now_time    = double(cat(1,Weights_time{z,:,d}));
                weights_now         = double(weights_now_z .* weights_now_time);
                
                rr = rgw(inds_now);
                sv = svgw(inds_now);
                
                % set radial velocity limits:
                rr(rr > 150) = NaN;
                sv(sv > 150) = NaN;
                
                %                 % and percentile limits for outliers
                %                 rr(abs(rr) > prctile(abs(rr),85)) = NaN;
                %                 sv(abs(sv) > prctile(abs(sv),85)) = NaN;
                
                % compute weighted variance
                rgw_var = nanvar(rr,weights_now);
                svgw_var = nanvar(sv,weights_now);
                
                % Resolved GW wind variance of 625 ~ a std of 25m/s.
                if length(rr) >= 20 && rgw_var < 1000
                    HWD.Data.ResolvedGWVariance(z,d)    = rgw_var;
                end
                
                % Sub-Volume GW variance of 4000 ~ a std of 65m/s.
                if length(sv) >= 20 && svgw_var < 4000
                    HWD.Data.SubVolumeGWVariance(z,d)   = svgw_var;
                end
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % you know, it'd be super useful to do the above but for composite
        % monthly sub-volume GW variance at each height. As in, a composite
        % day of each month. I think this'd be pretty robust.
        
        for mn = 1:12
            
            % now for each height step
            for z = 1:length(zrange)
                
                % choose ONLY the height weightings:
                inds_now        = double(cat(1,Inds{z,:,month(dayrange) == mn}));
                weights_now     = double(cat(1,Weights_z{z,:,month(dayrange) == mn}));
                
                % compute weighted variance:
                HWD.Data.MonthlyComp.SubVolumeGWVariance(z,mn) = ...
                    var(svgw(inds_now),weights_now,'omitnan');
                
%                 % also try for some directional variance:
%                 HWD.Data.MonthlyComp.SubVolumeGWVarianceZonal(z,mn) = ...
%                     var(svgw(inds_now),weights_now.*abs(sind(az(inds_now))),'omitnan');
%                 HWD.Data.MonthlyComp.SubVolumeGWVarianceMerid(z,mn) = ...
%                     var(svgw(inds_now),weights_now.*abs(cosd(az(inds_now))),'omitnan');
                
                % NEW: try sin^2 instead of abs(sin) for directional variance:
                HWD.Data.MonthlyComp.SubVolumeGWVarianceZonal(z,mn) = ...
                    var(svgw(inds_now),weights_now.*sind(az(inds_now)).^2,'omitnan');
                HWD.Data.MonthlyComp.SubVolumeGWVarianceMerid(z,mn) = ...
                    var(svgw(inds_now),weights_now.*cosd(az(inds_now)).^2,'omitnan');
                
                
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Project sub-volume GW perts into radial LOS direction for each
        % meteor echo, interpolated to the time and height of each:
        
        %% NOW LOAD EACH MPD FILE AND ADD THE GW RESIDUALS
        % It's annoying you have to essentially do this twice, but I want
        % the GW residuals in the MPD files to exactly correpond to all the
        % other fields, i.e. with all the meteors included not just those
        % within zenith limits etc that we use above.
        
        disp(['Computing GW residuals for all MPD .mat files for ' num2str(yr) '...'])
        
        dayrange = datenum(yr,01,01):datenum(yr,12,31);
        
        for d = 1:length(dayrange)
            
            yyyymmdd = datestr(dayrange(d),'yyyymmdd');
            
            try
                load([mpddirec yyyymmdd '_' site '_mpd.mat'])
            catch err
                %                 disp(['Couldn''t load: ' err.message])
                continue
            end
            
            % I'm pretty sure now that these should NOT be added in
            % quadrature. Once projected into the radial line of sight, the
            % vectors simply get added, right?
            
            % Project and assign the perturbations for every meteor:
            uflow = F.u(MPD.Data.Alt,MPD.Data.Time);
            vflow = F.v(MPD.Data.Alt,MPD.Data.Time);
            flow_radu = uflow .* sind(MPD.Data.Azimuth) .* sind(MPD.Data.ZenithAngle);
            flow_radv = vflow .* cosd(MPD.Data.Azimuth) .* sind(MPD.Data.ZenithAngle);
            
            MPD.Data.GWResiduals = MPD.Data.RadialVelocity - (flow_radu + flow_radv);
            
            % and do the same for the resolved GW perturbations:
            ugw = F.ugw(MPD.Data.Alt,MPD.Data.Time);
            vgw = F.vgw(MPD.Data.Alt,MPD.Data.Time);
            gw_radu = ugw .* sind(MPD.Data.Azimuth) .* sind(MPD.Data.ZenithAngle);
            gw_radv = vgw .* cosd(MPD.Data.Azimuth) .* sind(MPD.Data.ZenithAngle);
            
%             MPD.Data.ResolvedGWs = struct;
            MPD.Data.ResolvedGWs = gw_radu + gw_radv;
%             MPD.Data.ResolvedGWs.v = gw_radv;
%             MPD.Data.ResolvedGWs.Meta = 'Add these components linearly (not in quadrature) to get total GW variance. They''re already projected into the line of sight!';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
%             MPD.GWResidualMethod = 'gaussbins';
            
            % SAVE!
            save([mpddirec yyyymmdd '_' site '_mpd.mat'],'MPD')
            
        end
        
        
        %% SAVE ===============================================================
        
        savedirec = [mpddirec 'hwd/'];
        
        if ~exist(savedirec,'dir')
            mkdir(savedirec)
        end
        
        savename = [savedirec num2str(yr) '_' site '_hwd_gaussbins.mat'];
        save(savename,'HWD')
        
        disp(['Saving to: ' savename '...'])
        
        
    end % next year
    
end % next site




























