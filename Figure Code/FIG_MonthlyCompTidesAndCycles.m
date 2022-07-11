

%%%% EDIT: ALSO PLOTS ANNUAL AND SEMIANNUAL CYCLES OF TIDAL COMPONENTS

%%%% Plot long time period TIDAL AMPLITUDES and PHASES from MWR.

% Plot the monthly composite tidal amplitudes against height to show the
% large scale variation during the year. Do this for 24, 12, 8 and 6h
% tides.

runagain = 1;

comp30flag = 1;

if runagain
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SITE AND TIME
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    site = 'kep';
    
    years = 2016:2020;
    
    direc = '/Users/neil/data/MeteorRadar/';
    % direc = '/Volumes/SDRed/data/MeteorRadar/';
    
    % First, stick all the monthly composite mean winds together for the
    % all-time plot:
    
    Comp.u = [];
    Comp.v = [];
    Comp.gwvar = [];
    Comp.walt = [];
    Comp.wtime = [];
    
    % also tides?
    Comp.Tides.u.A = [];
    Comp.Tides.u.B = [];
    Comp.Tides.v.A = [];
    Comp.Tides.v.B = [];
    
    Comp.Monthly.Tides.u.A = [];
    Comp.Monthly.Tides.u.B = [];
    Comp.Monthly.Tides.v.A = [];
    Comp.Monthly.Tides.v.B = [];
    
    Comp.Tides.OneDay.u.A = [];
    Comp.Tides.OneDay.u.B = [];
    Comp.Tides.OneDay.u.C = [];
    Comp.Tides.OneDay.v.A = [];
    Comp.Tides.OneDay.v.B = [];
    Comp.Tides.OneDay.v.C = [];
    
    % also raw hourly winds?
    Comp.uu = [];
    Comp.vv = [];
    Comp.tt = [];
    
    % and some mean winds?
    Comp.MeanWinds.u = [];
    Comp.MeanWinds.v = [];
    
    for yr = years
        
        yrstr = num2str(yr);
        
        load([direc site '/matlab/hwd/' yrstr '_' site '_hwd_gaussbins.mat'])
        
        Comp.u      = cat(2,Comp.u,sq(nanmean(HWD.Data.MonthlyComp.u,2)));
        Comp.v      = cat(2,Comp.v,sq(nanmean(HWD.Data.MonthlyComp.v,2)));
        
        Comp.Tides.u.A = cat(2,Comp.Tides.u.A,HWD.Data.Tides.TidalComponents.u.A);
        Comp.Tides.u.B = cat(2,Comp.Tides.u.B,HWD.Data.Tides.TidalComponents.u.B);
        Comp.Tides.v.A = cat(2,Comp.Tides.v.A,HWD.Data.Tides.TidalComponents.v.A);
        Comp.Tides.v.B = cat(2,Comp.Tides.v.B,HWD.Data.Tides.TidalComponents.v.B);
        
        Comp.Monthly.Tides.u.A = cat(2,Comp.Monthly.Tides.u.A,HWD.Data.MonthlyComp.Tides.TidalComponents.u.A);
        Comp.Monthly.Tides.u.B = cat(2,Comp.Monthly.Tides.u.B,HWD.Data.MonthlyComp.Tides.TidalComponents.u.B);
        Comp.Monthly.Tides.v.A = cat(2,Comp.Monthly.Tides.v.A,HWD.Data.MonthlyComp.Tides.TidalComponents.v.A);
        Comp.Monthly.Tides.v.B = cat(2,Comp.Monthly.Tides.v.B,HWD.Data.MonthlyComp.Tides.TidalComponents.v.B);
        
        Comp.Tides.OneDay.u.A = cat(2,Comp.Tides.OneDay.u.A,HWD.Data.Tides.OneDay.TidalComponents.u.A);
        Comp.Tides.OneDay.u.B = cat(2,Comp.Tides.OneDay.u.B,HWD.Data.Tides.OneDay.TidalComponents.u.B);
        Comp.Tides.OneDay.u.C = cat(2,Comp.Tides.OneDay.u.C,HWD.Data.Tides.OneDay.TidalComponents.u.C);
        Comp.Tides.OneDay.v.A = cat(2,Comp.Tides.OneDay.v.A,HWD.Data.Tides.OneDay.TidalComponents.v.A);
        Comp.Tides.OneDay.v.B = cat(2,Comp.Tides.OneDay.v.B,HWD.Data.Tides.OneDay.TidalComponents.v.B);
        Comp.Tides.OneDay.v.C = cat(2,Comp.Tides.OneDay.v.C,HWD.Data.Tides.OneDay.TidalComponents.v.C);
        
        
        Comp.uu      = cat(2,Comp.uu,HWD.Data.u);
        Comp.vv      = cat(2,Comp.vv,HWD.Data.v);
        Comp.tt      = cat(2,Comp.tt,HWD.Data.Time);
        
        Comp.MeanWinds.u = cat(2,Comp.MeanWinds.u,HWD.Data.u-HWD.Data.Tides.OneDay.AllTides.u);
        Comp.MeanWinds.v = cat(2,Comp.MeanWinds.v,HWD.Data.v-HWD.Data.Tides.OneDay.AllTides.v);
        
        Comp.walt = cat(2,Comp.walt,sq(nanmean(HWD.Data.MonthlyComp.walt,2)));
        Comp.wtime = cat(2,Comp.wtime,sq(nanmean(HWD.Data.MonthlyComp.wtime,2)));
        
    end
    
    %%%%
    % next, load all the mpd files for the year range and cat them into a
    % MASSIVE array. This will be HUGE I hope it fits in memory.
    
    disp(['Loading all MPDs...'])
    mpddirec = ['/Users/neil/data/MeteorRadar/' site '/matlab/'];
    
    allMPD = struct;
    dayrange = datenum(years(1),01,01):datenum(years(end),12,31);
    
    for d = 1:length(dayrange)
        
        yyyymmdd = datestr(dayrange(d),'yyyymmdd');
        
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
    MWD.Thresholds.ZenithLimits = [15 65]; % use [15 65] for computing mean winds, impose stricter for GWs.
    inds = zen >= MWD.Thresholds.ZenithLimits(1) & zen <= MWD.Thresholds.ZenithLimits(2);
    
    % avoid silly horz velocities:
    inds = inds & vhorz < 200;
    
    az = az(inds);
    tim = tim(inds);
    alt = alt(inds);
    vrad = vrad(inds);
    vhorz = vhorz(inds);
    zen = zen(inds);
    
    hourofday = mod(tim,1) * 24;
    dn = daynumber(tim);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% get rolling 30 day composite tidal fit...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if comp30flag
        
    disp('Assembling a rolling 30 day tidal fit...')
    
    days = datenum(min(years),01,01):10:datenum(max(years)+1,01,01);
    
    daywindow = 30;
    
    std_z = 1.275; % FWHM for height gaussian
    zrange = 76:1:105;
    
    Comp30 = struct;
    Comp30.u = nan(length(zrange),length(days));
    Comp30.v = nan(length(zrange),length(days));
    Comp30.walt = nan(length(zrange),length(days));
    
    Comp30.Tides.u.A = nan(length(zrange),length(days),4);
    Comp30.Tides.u.B = nan(length(zrange),length(days),4);
    Comp30.Tides.v.A = nan(length(zrange),length(days),4);
    Comp30.Tides.v.B = nan(length(zrange),length(days),4);
    
    std_z = 1.275; % gives approx 3km FWHM
    std_time = 0.85/24; % DAYS, gives approx 2hrs FWHM
    nmet_threshold = 20;
    
    tidalperiods = [24 12 8 6] ./ 24;
    tidalnmets = 20; % number of meteors required for tidal fitting
    
    for d = 1:length(days)
        
        dd = days(d);
        
        disp(datestr(dd,'yyyymmdd'))
        
        % extract day range
        drange = (dd-daywindow/2):(dd+daywindow/2);
        
        %     % deal with wraparound (stolen from wrapTo360):
        %     positiveInput = (drange > 0);
        %     drange = mod(drange, 365);
        %     drange((drange == 0) & positiveInput) = 365;
        
        % find the indeces of this day range:
        dayinds = ismember(floor(tim),drange);
        % thats a lot of meteors...
        
        altday       = alt(dayinds);
        timday       = tim(dayinds);
        azday       = az(dayinds);
        vhorzday    = vhorz(dayinds);
        
        %%%% Now let's do our wind and tidal fits:
        
        % To save a bit of time, pre-compute the gaussian weighting
        % functions for the height bins, they'll be the same each time:
        Gauss_z = struct;
        for z = 1:length(zrange)
            Gauss_z(z).vec = exp(- ((altday - zrange(z)).^2) ./ (2 * std_z^2));
        end
        Gauss_h = struct;
        for h = 1:24
            Gauss_h(h).vec = exp(- ((mod(timday,1) - ((h-1)/24)).^2) ./ (2 * std_time^2));
        end
        
        % temp, gets rewritten every time
        CompDay.u       = nan(length(zrange),24);
        CompDay.v       = nan(length(zrange),24);
        CompDay.walt    = nan(length(zrange),24);
        CompDay.wtime    = nan(length(zrange),24);
        
        % for each height step
        for z = 1:length(zrange)
            
            % for each hour
            for h = 0:23
                
                w = Gauss_z(z).vec .* Gauss_h(h+1).vec;
                
                % only choose meteors within 2STDs:
                rng = w > 0.05;
                
                % Choose whether to do the sine fit:
                if sum(rng) >= nmet_threshold
                    
                    [yfit,F] = nph_sinefit(azday(rng),vhorzday(rng),360,'weights',w(rng));
                    
                    % get the composite altitude to which this tidal
                    % fit corresponded to:
                    CompDay.walt(z,h+1) = wmean(altday(rng),w(rng));
%                     CompDay.wtime(z,h+1) = wmean(timday(rng),w(rng));
                    CompDay.wtime(z,h+1) = wmean(mod(timday(rng),1).*24,w(rng));
                    
                else
                    F = [NaN NaN NaN];
                end
                
                % Check for silly numbers for tides:
                F(abs(F) > 200) = NaN;
                
                % And subscribe to today's temp composite day storage:
                CompDay.u(z,h+1) = F(2);
                CompDay.v(z,h+1) = F(1);
                
            end % next hour
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Great, we should have a line of winds for this height. Let's
            %%%% fit some tides:
            
            % and extract a list of u and v components and their
            % times:
            uu = CompDay.u(z,:);
            vv = CompDay.v(z,:);
%             tt = (0:23)/24;
            tt = CompDay.wtime(z,:);
%             if all(~isnan(tt))
%                 return
%             end

            
            % and fit sine waves to them!
            if sum(~isnan(uu)) >= 18
                % require at least 18 hours of data at z=89km
                [~,Fu] = nph_sinefit(tt,uu,tidalperiods);
                [~,Fv] = nph_sinefit(tt,vv,tidalperiods);
            else
                Fu      = nan(length(tidalperiods),3);
                Fv      = nan(length(tidalperiods),3);
            end
            
            % Check for silly values in the tidal fit
            % First check for large values:
            Fu(abs(Fu) > 100) = NaN;
            Fv(abs(Fv) > 100) = NaN;
            Fu(quadadd(Fu(:,1),Fu(:,2)) > 100,1:2) = NaN;
            Fv(quadadd(Fv(:,1),Fv(:,2)) > 100,1:2) = NaN;
            
            %%%% SUBSCRIBE!
            %         Comp.Tides.u.A(z,d,:) = Fu;
            %         Comp.Tides.u.A(z,d,:) = Fu;
            for td = 1:length(tidalperiods)
                for abc = 1:3
                    Comp30.Tides.u.(upper(alphabet(abc)))(z,d,td) = Fu(td,abc);
                    Comp30.Tides.v.(upper(alphabet(abc)))(z,d,td) = Fv(td,abc);
                end
            end
            
            
            
            %%%% Finally, use the fitted tidal daily mean to subscribe to the
            %%%% zona, and merid wind:
            Comp30.u(z,d) = Fu(length(tidalperiods),3);
            Comp30.v(z,d) = Fv(length(tidalperiods),3);
            
            % get the composite altitude to which this tidal
            % fit corresponded to:
            Comp30.walt(z,d)    = nanmean(linearise(CompDay.walt(z,:)));
            Comp30.wtime(z,d)   = days(d) + nanmean(linearise(CompDay.wtime(z,:)));
            
        end % next height
        
    end
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% get an all-years composite for the specified years...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('Assembling an all-years composite for the specified years...')
    
    daynumbers = ceil(linspace(1,365,60));
    % days = sort([365 daynumber(sort([datenum(min(years),1:12,01) datenum(min(years),1:12,11) datenum(min(years),1:12,21)]))]);
    
    daywindow = 30;
    
    std_z = 1.275; % FWHM for height gaussian
    zrange = 76:1:105;
    
    CompYear = struct;
    CompYear.u = nan(length(zrange),length(daynumbers));
    CompYear.v = nan(length(zrange),length(daynumbers));
    CompYear.walt = nan(length(zrange),length(daynumbers));
    
    CompYear.Tides.u.A = nan(length(zrange),length(daynumbers),4);
    CompYear.Tides.u.B = nan(length(zrange),length(daynumbers),4);
    CompYear.Tides.v.A = nan(length(zrange),length(daynumbers),4);
    CompYear.Tides.v.B = nan(length(zrange),length(daynumbers),4);
    
    std_z = 1.275; % gives approx 3km FWHM
    std_time = 0.85/24; % DAYS, gives approx 2hrs FWHM
    nmet_threshold = 20;
    
    tidalperiods = [24 12 8 6] ./ 24;
    tidalnmets = 20; % number of meteors required for tidal fitting
    
    
    for d = 1:length(daynumbers)
        
        dd = daynumbers(d);
        
        disp(num2str(dd))
        
        % extract day range
        drange = (dd-daywindow/2):(dd+daywindow/2);
        
        % deal with wraparound (stolen from wrapTo360):
        positiveInput = (drange > 0);
        drange = mod(drange, 365);
        drange((drange == 0) & positiveInput) = 365;
        
        % find the indeces of this day range:
        dayinds = ismember(dn,drange);
        % thats a lot of meteors...
        
        
        altday       = alt(dayinds);
        timday       = tim(dayinds);
        azday       = az(dayinds);
        vhorzday    = vhorz(dayinds);
        
        %%%% Now let's do our wind and tidal fits:
        
        % To save a bit of time, pre-compute the gaussian weighting
        % functions for the height bins, they'll be the same each time:
        Gauss_z = struct;
        for z = 1:length(zrange)
            Gauss_z(z).vec = exp(- ((altday - zrange(z)).^2) ./ (2 * std_z^2));
        end
        Gauss_h = struct;
        for h = 1:24
            Gauss_h(h).vec = exp(- ((mod(timday,1) - ((h-1)/24)).^2) ./ (2 * std_time^2));
        end
        
        % temp, gets rewritten every time
        CompDay.u       = nan(length(zrange),24);
        CompDay.v       = nan(length(zrange),24);
        CompDay.walt    = nan(length(zrange),24);
        CompDay.wtime   = nan(length(zrange),24);
        
        % for each height step
        for z = 1:length(zrange)
            
            % for each hour
            for h = 0:23
                
                w = Gauss_z(z).vec .* Gauss_h(h+1).vec;
                
                % only choose meteors within 2STDs:
                rng = w > 0.05;
                
                % Choose whether to do the sine fit:
                if sum(rng) >= nmet_threshold
                    
                    [yfit,F] = nph_sinefit(azday(rng),vhorzday(rng),360,'weights',w(rng));
                    
                    %                 % get the composite altitude to which this tidal
                    %                 % fit corresponded to:
                    %                 CompDay.walt(z,h+1) = wmean(altday(rng),w(rng));
                    CompDay.walt(z,h+1)  = wmean(altday(rng),w(rng));
                    CompDay.wtime(z,h+1) = wmean(mod(timday(rng),1).*24,w(rng));
                    
                    
                else
                    F = [NaN NaN NaN];
                end
                
                % Check for silly numbers for tides:
                F(abs(F) > 200) = NaN;
                
                % And subscribe to today's temp composite day storage:
                CompDay.u(z,h+1) = F(2);
                CompDay.v(z,h+1) = F(1);
                
            end % next hour
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Great, we should have a line of winds for this height. Let's
            %%%% fit some tides:
            
            % and extract a list of u and v components and their
            % times:
            uu = CompDay.u(z,:);
            vv = CompDay.v(z,:);
%             tt = (0:23)/24;
            tt = CompDay.wtime(z,:);
%             if all(~isnan(tt))
%                 return
%             end
            
            % and fit sine waves to them!
            if sum(~isnan(uu)) >= 18
                % require at least 18 hours of data at z=89km
                [~,Fu] = nph_sinefit(tt,uu,tidalperiods);
                [~,Fv] = nph_sinefit(tt,vv,tidalperiods);
            else
                Fu      = nan(length(tidalperiods),3);
                Fv      = nan(length(tidalperiods),3);
            end
            
            % Check for silly values in the tidal fit
            % First check for large values:
            Fu(abs(Fu) > 100) = NaN;
            Fv(abs(Fv) > 100) = NaN;
            Fu(quadadd(Fu(:,1),Fu(:,2)) > 100,1:2) = NaN;
            Fv(quadadd(Fv(:,1),Fv(:,2)) > 100,1:2) = NaN;
            
            %%%% SUBSCRIBE!
            %         Comp.Tides.u.A(z,d,:) = Fu;
            %         Comp.Tides.u.A(z,d,:) = Fu;
            for td = 1:length(tidalperiods)
                for abc = 1:3
                    CompYear.Tides.u.(upper(alphabet(abc)))(z,d,td) = Fu(td,abc);
                    CompYear.Tides.v.(upper(alphabet(abc)))(z,d,td) = Fv(td,abc);
                end
            end
            
            %%%% Finally, use the fitted tidal daily mean to subscribe to the
            %%%% zona, and merid wind:
            CompYear.u(z,d) = Fu(length(tidalperiods),3);
            CompYear.v(z,d) = Fv(length(tidalperiods),3);
            
            %                 % get the composite altitude to which this tidal
            %                 % fit corresponded to:
            %                 CompDay.walt(z,h+1) = wmean(altday(rng),w(rng));
            
            
        end % next height
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% also fit some annual and semiannual cycles to the data:
    
    cycleperiods = 365 .* [1 0.5];
    tidalperiods = [24 12 8 6] ./ 24;
    
    Cycles.u = nan(length(zrange),length(cycleperiods)+1);
    Cycles.v = nan(length(zrange),length(cycleperiods)+1);
    Cycles.Tides.u = nan(length(zrange),length(cycleperiods),length(tidalperiods));
    Cycles.Tides.v = nan(length(zrange),length(cycleperiods),length(tidalperiods));
    Cycles.Tides.uph = nan(length(zrange),length(cycleperiods),length(tidalperiods));
    Cycles.Tides.vph = nan(length(zrange),length(cycleperiods),length(tidalperiods));
    
    
    %%%% fit annual and semiannual cycles to the composite year
    for z = 1:length(zrange)
        [~,Fu] = nph_sinefit(daynumbers,CompYear.u(z,:),cycleperiods);
        [~,Fv] = nph_sinefit(daynumbers,CompYear.v(z,:),cycleperiods);
        % Annual and Semiannual cycles in wind...
        Cycles.u(z,1) = Fu(2,3);
        Cycles.v(z,1) = Fv(2,3);
        Cycles.u(z,2:3) = quadadd(Fu(1:2,1),Fu(1:2,2))';
        Cycles.v(z,2:3) = quadadd(Fv(1:2,1),Fv(1:2,2))';
        % do the same for tidal periods...
        for td = 1:length(tidalperiods)
            [~,Fu] = nph_sinefit(daynumbers,quadadd(CompYear.Tides.u.A(z,:,td),CompYear.Tides.u.B(z,:,td)),cycleperiods);
            [~,Fv] = nph_sinefit(daynumbers,quadadd(CompYear.Tides.v.A(z,:,td),CompYear.Tides.v.B(z,:,td)),cycleperiods);
            Cycles.Tides.u(z,1,td) = Fu(2,3);
            Cycles.Tides.v(z,1,td) = Fv(2,3);
            Cycles.Tides.u(z,2:3,td) = quadadd(Fu(1:2,1),Fu(1:2,2))';
            Cycles.Tides.v(z,2:3,td) = quadadd(Fv(1:2,1),Fv(1:2,2))';
            Cycles.Tides.uph(z,2:3,td) = atan2(Fu([2 1],1),Fu([2 1],2))';
            Cycles.Tides.vph(z,2:3,td) = atan2(Fv([2 1],1),Fv([2 1],2))';
        end
    end
    
    %%%% now to work out standard deviations, do the same but for all 5 years
    %%%% we have:
    days = datenum(min(years),01,01):10:datenum(max(years)+1,01,01);
    
    Cycles.Years.u = nan(length(zrange),length(cycleperiods)+1,length(years));
    Cycles.Years.v = nan(length(zrange),length(cycleperiods)+1,length(years));
    Cycles.Years.Tides.u = nan(length(zrange),length(cycleperiods),length(tidalperiods),length(years));
    Cycles.Years.Tides.v = nan(length(zrange),length(cycleperiods),length(tidalperiods),length(years));
    Cycles.Years.Tides.uph = nan(length(zrange),length(cycleperiods),length(tidalperiods),length(years));
    Cycles.Years.Tides.vph = nan(length(zrange),length(cycleperiods),length(tidalperiods),length(years));
    
    for y = 1:length(years)
        yr = years(y);
        % select this year:
        dayinds = inrange(days,[datenum(yr,01,01) datenum(yr+1,01,01)]);
        for z = 1:length(zrange)
            [~,Fu] = nph_sinefit(days(dayinds),Comp30.u(z,dayinds),cycleperiods);
            [~,Fv] = nph_sinefit(days(dayinds),Comp30.v(z,dayinds),cycleperiods);
            % Annual and Semiannual cycles in wind...
            Cycles.Years.u(z,1,y) = Fu(2,3);
            Cycles.Years.v(z,1,y) = Fv(2,3);
            Cycles.Years.u(z,2:3,y) = quadadd(Fu(1:2,1),Fu(1:2,2))';
            Cycles.Years.v(z,2:3,y) = quadadd(Fv(1:2,1),Fv(1:2,2))';
            % do the same for tidal periods...
            for td = 1:length(tidalperiods)
                [~,Fu] = nph_sinefit(days(dayinds),quadadd(Comp30.Tides.u.A(z,dayinds,td),Comp30.Tides.u.B(z,dayinds,td)),cycleperiods);
                [~,Fv] = nph_sinefit(days(dayinds),quadadd(Comp30.Tides.v.A(z,dayinds,td),Comp30.Tides.v.B(z,dayinds,td)),cycleperiods);
                Cycles.Years.Tides.u(z,1,td,y) = Fu(2,3);
                Cycles.Years.Tides.v(z,1,td,y) = Fv(2,3);
                Cycles.Years.Tides.u(z,2:3,td,y) = quadadd(Fu(1:2,1),Fu(1:2,2))';
                Cycles.Years.Tides.v(z,2:3,td,y) = quadadd(Fv(1:2,1),Fv(1:2,2))';
                Cycles.Years.Tides.uph(z,2:3,td,y) = atan2(Fu([2 1],1),Fu([2 1],2))';
                Cycles.Years.Tides.vph(z,2:3,td,y) = atan2(Fv([2 1],1),Fv([2 1],2))';
            end
        end
    end
    
    %%%% finally compute those STDs...
    Cycles.STDs.u = std(Cycles.Years.u,[],3,'omitnan');
    Cycles.STDs.v = std(Cycles.Years.v,[],3,'omitnan');
    for td = 1:length(tidalperiods)
        Cycles.STDs.Tides.u = std(Cycles.Years.Tides.u,[],4,'omitnan');
        Cycles.STDs.Tides.v = std(Cycles.Years.Tides.v,[],4,'omitnan');
        Cycles.STDs.Tides.uph = circ_std(Cycles.Years.Tides.uph,[],[],4);
        Cycles.STDs.Tides.vph = circ_std(Cycles.Years.Tides.vph,[],[],4);
    end
    
    %% or just load a previous run?
    % save('/Users/neil/Drive/MATLAB/MonthlyCompTidesData.mat','Comp','Comp30','CompDay','CompYear','Cycles','site','years','direc','HWD','daynumbers','zrange');
    
    
else
    disp('loading...')
    load('/Users/neil/Drive/MATLAB/MonthlyCompTidesData.mat')
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PLOT OF TIDAL AMPLITUDES/PHASES AGAINST ALL TIME

for td = 4
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure; hold all; whitefig; figpos([1 0.475])
    
    %-------------------------------------------------------
    vert_gap = 0.08;        horz_gap = 0.03;
    lower_marg = 0.075;     upper_marg = 0.025;
    left_marg = 0.065;      right_marg = 0.075;
    
    rows = 2; cols = 8;
    
    subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);
    
    %--------------------------------------------------------
    
    fs = 20;
    
    %--------------------------------------------------------
    
    % select tide:
    % td = 1;
    
    % bump for the axes labels
    bump = (td-1) * 4;
    
    % select amplitude or phase:
    ap = 'p'; % a | p
    
    %%%% First, assemble and over-interpolate the tidal amplitudes and phases:
    tdh = 24/td;
    
    TD.uA = Comp30.Tides.u.A(:,:,td);
    TD.uB = Comp30.Tides.u.B(:,:,td);
    TD.vA = Comp30.Tides.v.A(:,:,td);
    TD.vB = Comp30.Tides.v.B(:,:,td);
    
    %%%% Over-Interpolation?
    odays = datenum(2016,01,01):datenum(2021,01,01);
    ozrange = 76:0.25:105;
    for uv = 'uv'
        for ab = 'AB'
            nanlocs = isnan(Comp30.wtime) | isnan(Comp30.walt) | isnan(TD.([uv ab]));
            %         F = griddedInterpolant({zrange,days},TD.([uv ab]),'linear','none');
            F = scatteredInterpolant(linearise(Comp30.walt(~nanlocs)),linearise(Comp30.wtime(~nanlocs)),linearise(TD.([uv ab])(~nanlocs)),'linear','none');
            TD.([uv ab]) = F({ozrange,odays});
        end
    end
    
    %%%% first, compute tidal amps and phases:
    umag = quadadd(TD.uA,TD.uB);
    vmag = quadadd(TD.vA,TD.vB);
    uph  = atan2(TD.uB,TD.uA);
    vph  = atan2(TD.vB,TD.vA);
    
    % convert to local time:
    uph = wrapToPi(uph + ((HWD.Meta.TimeZone/24) * (2*pi)));
    vph = wrapToPi(vph + ((HWD.Meta.TimeZone/24) * (2*pi)));
    % convert to hours...
    uph = (uph./(2*pi)) .* tdh;
    vph = (vph./(2*pi)) .* tdh;
    % wrap:
    uph(uph < 0) = tdh + uph(uph < 0);
    vph(vph < 0) = tdh + vph(vph < 0);
    
    
    %%%% Now assemble our all-years composite...
    odaynumbers = 0:0.5:365; % and over-interpolate (for phase)
    TD.uAy = interp1(daynumbers,interp1(zrange,CompYear.Tides.u.A(:,:,td),ozrange)',odaynumbers)';
    TD.uBy = interp1(daynumbers,interp1(zrange,CompYear.Tides.u.B(:,:,td),ozrange)',odaynumbers)';
    TD.vAy = interp1(daynumbers,interp1(zrange,CompYear.Tides.v.A(:,:,td),ozrange)',odaynumbers)';
    TD.vBy = interp1(daynumbers,interp1(zrange,CompYear.Tides.v.B(:,:,td),ozrange)',odaynumbers)';
    % umagy = quadadd(CompYear.Tides.u.A(:,:,td),CompYear.Tides.u.B(:,:,td));
    % vmagy = quadadd(CompYear.Tides.v.A(:,:,td),CompYear.Tides.v.B(:,:,td));
    % uphy  = atan2(CompYear.Tides.u.B(:,:,td),CompYear.Tides.u.A(:,:,td));
    % vphy  = atan2(CompYear.Tides.v.B(:,:,td),CompYear.Tides.v.A(:,:,td));
    umagy = quadadd(TD.uAy,TD.uBy);
    vmagy = quadadd(TD.vAy,TD.vBy);
    uphy  = atan2(TD.uBy,TD.uAy);
    vphy  = atan2(TD.vBy,TD.vAy);
    % convert to local time:
    uphy = wrapToPi(uphy + ((HWD.Meta.TimeZone/24) * (2*pi)));
    vphy = wrapToPi(vphy + ((HWD.Meta.TimeZone/24) * (2*pi)));
    % convert to hours...
    uphy = (uphy./(2*pi)) .* tdh;
    vphy = (vphy./(2*pi)) .* tdh;
    % wrap:
    uphy(uphy < 0) = tdh + uphy(uphy < 0);
    vphy(vphy < 0) = tdh + vphy(vphy < 0);
    
    
    % =========================================================================
    
    T = tiledlayout(2,6,'tilespacing','compact');
    
    spans = {...
        [1 5],...
        [1 5],...
        [1 1],...
        [1 1]};
    
    % different limits for different tides...
    switch td
        case 1 % 24h tide
            ampclims = [0 20];
            ampclevs = 0:50;
            ampclines = 0:5:50;
            ampcbarticks = 0:5:60;
            ampcbarminorticks = 0:2.5:60;
            phaseminorclines = 0:3:tdh;
            darkclines = [10 15];
        case 2 % 12h tide
            ampclims = [0 70];
            ampclevs = 0:5:100;
            ampclines = 0:20:100;
            ampcbarticks = 0:20:100;
            ampcbarminorticks = 0:10:100;
            phaseminorclines = 0:1:tdh;
            darkclines = [40 40];
        case 3 % 8h tide
            ampclims = [0 20];
            ampclevs = 0:50;
            ampclines = 0:5:50;
            ampcbarticks = 0:5:60;
            ampcbarminorticks = 0:2.5:60;
            phaseminorclines = 0:1:tdh;
            darkclines = [5 10];
        case 4 % 6h tide
            ampclims = [0 6];
            ampclevs = 0:0.5:20;
            ampclines = 0:2:10;
            ampcbarticks = 0:2:10;
            ampcbarminorticks = ampclevs;
            phaseminorclines = 0:1:tdh;
            darkclines = [4 4];
    end
    phaseclims = [0 tdh];
    phaseclevs = 0:1:tdh;
    phaseclines = tdh .* [0.25 0.5 0.75];
    phasecbarticks = 0:(tdh/4):(tdh);
    phasecbarminorticks = 0:(tdh/8):(tdh);
    
    
    
    switch ap
        case 'a'
            TBP = {...
                umag,...
                vmag,...
                umagy,...
                vmagy};
            clims = ampclims;
            clevs = ampclevs;
            clines = ampclines;
            cbarticks = ampcbarticks;
            cbarminorticks = ampcbarminorticks;
            clabs = ['   \bf{ ms^{-1} }'];
            cmap = nph_saturate(cbrew('nph_modspectral',18),1.2);
            ampphase = 'amp';
            textcolor = 'w';
            textbordercolor = 'k';
        case 'p'
            TBP = {...
                uph,...
                vph,...
                uphy,...
                vphy};
            clims = phaseclims;
            clevs = phaseclevs;
            clines = phaseclines;
            minorclines = phaseminorclines;
            cbarticks = phasecbarticks;
            cbarminorticks = phasecbarminorticks;
            clabs = '\bf{Hours (Local Time)}';
            cmap = nph_saturate(cbrew('nph_CyclicRainbow',24),1);
            ampphase = 'phase';
            textcolor = 'w';
            textbordercolor = 'k';
    end
    
    % clabs = {...
    %     'ms^{-1}','ms^{-1}','Hours (LT)','Hours (LT)'};
    
    % tits = {...
    %     [num2str(tdh) 'h TIDE, ZONAL COMPONENT'],...
    %     [num2str(tdh) 'h TIDE, MERIDIONAL COMPONENT'],...
    %     ' ',...
    %     ' '};
    tits = {...
        'ZONAL COMPONENT',...
        'MERIDIONAL COMPONENT',...
        'ZONAL',...
        'MERID.'};
    
    letternum = [1 3 2 4];
    
    for ax = 1:4
        
        %     subplot(rows,cols,axs{ax})
        nexttile(spans{ax})
        
        %     if ax == 3, continue; end
        %     if ax == 4, continue; end
        
        axx = gca;
        
        % choose wot 2 plot n wot not
        switch ax
            case {1,2}
                %             X = Comp.wtime; Y = Comp.walt;
                [X,Y] = meshgrid(odays,ozrange);
                smoo1 = [0.25 0.75]; smoo2 = [0.25 0.75];
            case {3,4}
                [X,Y] = meshgrid(datenum(2019,00,odaynumbers),ozrange);
                smoo1 = [0.25 0.75]; smoo2 = [0.25 0.75];
        end
        
        % allll the contours...
        tbp = TBP{ax};
        nanlocs = isnan(tbp); smoo = movmean(tbp,30,2,'omitnan');
        tbp(nanlocs) = smoo(nanlocs);
        switch ap
            case 'a'
                tbp = imgaussfilt(imgaussfilt(tbp,smoo1),smoo2);
                if td == 4
                    tbp = imgaussfilt(imgaussfilt(tbp,smoo1),smoo2);
                end
                % filled contours:
                hold on; contourf(X,Y,tbp,clevs,'edgecolor','none');
                % contour lines:
                hold on; [C,h] = contour(X,Y,tbp,clines,'edgecolor',rgbtrip(.25),'linewi',1.25);
                clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
                % and bold zero line
                hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','w','linewi',4);
                clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000,'color','w','fontweight','bold');
                hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linewi',2);
                clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000);
                
                % and dark contour lines:
                % contour lines:
                hold on; [C,h] = contour(X,Y,tbp,darkclines,'edgecolor',rgbtrip(.25),'linewi',1.25);
                clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.25),'fontweight','bold');
                
                
            case 'p'
                % pcolor instead for phase?
                hold on; pcolor(X,Y,tbp); shat;
                
                %%%% shall we try some contours for phase too?
                %%%% just try a few...
                
                % and some minor clines:
                minorclines = minorclines(~ismember(minorclines,clines));
                hold on; [C,h] = contour(X,Y,tbp,minorclines,'edgecolor',rgbtrip(.5),'linewi',0.5);
                
                % first zero phase:
                tbp2 = tbp; tbp2(tbp > tdh/2) = -tdh + tbp2(tbp > tdh/2);
                nanmask = abs(diff(tbp2([1 1:end],:),[],1)) > tdh/2 | abs(diff(tbp2(:,[1 1:end]),[],2)) > tdh/2;
                tbp2(nanmask) = NaN;
                hold on; [C,h] = contour(X,Y,tbp2,[0 0],'edgecolor',rgbtrip(.25),'linewi',1.25);
                clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
                
                % before you do the others, exclude any wraparound bits:
                nanmask = abs(diff(tbp([1 1:end],:),[],1)) > tdh/2 | abs(diff(tbp(:,[1 1:end]),[],2)) > tdh/2;
                tbp(nanmask) = NaN;
                hold on; [C,h] = contour(X,Y,tbp,clines,'edgecolor',rgbtrip(.25),'linewi',1.25);
                clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
                
        end
        
        
        % LIMITS
        ylim([77 103])
        axx.YTick = 70:5:110;
        axx.YMinorTick = 'on';
        %     switch ax
        %         case {1,2}
        axx.YAxis.MinorTickValues = 70:1:110;
        %         case {3,4}
        %             axx.YAxis.MinorTickValues = 70:110;
        %     end
        
        
        switch ax
            case {1,2}
                xlim([datenum(2016,02,15) datenum(2020,11,15)])
                % Just show the Januarys as major ticks
                axx.XTick = datenum(years,01,01);
                axx.XMinorTick = 'on';
                axx.XAxis.MinorTickValues = datenum(2016,1:(length(years)*12),01);
                datetick('x',' ','keepticks','keeplimits')
                %             axx.XTickLabel = {};
                % other months as ticks using text:
                xtix = datenum(years(1),1:(length(years)*12),15);
                for xt = xtix
                    if inrange(xt,xlim) && ~any(xt == axx.XTick)
                        mn = monthname(month(xt));
                        text(xt,75.5,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
                    end
                end
                if ax == 2
                    % Years as lower bold numbers:
                    xtix = datenum(years,07,01);
                    for xt = xtix
                        if inrange(xt,xlim) && ~any(xt == axx.XTick)
                            yrstr = datestr(xt,'yyyy');
                            text(xt,72.5,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
                        end
                    end
                end
                
            case {3,4}
                xlim(datenum(2019,[1 13],01))
                % Just show the Januarys as major ticks
                axx.XTick = datenum(2019,1:13,01);
                axx.XMinorTick = 'on';
                axx.XAxis.MinorTickValues = datenum(2016,1:(length(years)*12),01);
                datetick('x',' ','keepticks','keeplimits')
                % other months as ticks using text:
                xtix = datenum(years(1),1:(length(years)*12),16);
                for xt = xtix
                    if inrange(xt,xlim) && ~any(xt == axx.XTick)
                        mn = monthname(month(xt));
                        text(xt,76,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
                    end
                end
                %             axx.XTick = datenum(2019,1:13,01);
                %             datetick('x','m','keepticks','keeplimits')
                %             xlim(datenum(2019,[1 13],01))
                if ax == 4
                    xlabel('Average Year','fontsize',0.8*fs,'fontweight','bold')
                end
        end
        
        
        switch ax
            case {1,2}
                ylabel('Altitude (km)')
        end
        
        % line around the axis
        hold on; plot(axx.XLim([1 2 2]),axx.YLim([2 2 1]),'color','k','linewi',1.5)
        
        drawnow;
        
        % COLORI
        %     cmap = cbrew('RdBu',length(clevs{ax}));
        %     cmap = cbrew('nph_RdYlBu',length(clevs{ax}));
        %     cmap = cbrew('nph_BuOr',length(clevs{ax}));
        %     cmap = cbrew('BrBG',length(clevs{ax}));
        %     cmap = flipud(cbrew('PiYG',length(clevs{ax})));
        %     cmap = cbrew('nph_Spectral',length(clevs{ax}));
        %     cmap = cbrew('nph_modspectral',length(clevs{ax}));
        %     cmap = nph_saturate(cbrew('nph_modspectral',18),1.2);
        %     cmap = cbrew('nph_RdBuPastel',length(clevs{ax}));
        %     colormap(gca,cmap)
        
        % sort out double white in the middle:
        %     cmap = [cmap(1:floor(length(clevs{ax})/2),:) ; [1 1 1] ; cmap(floor(length(clevs{ax})/2)+1:end,:)];
        
        % saturate?
        %     cmap_hsv = rgb2hsv(cmap);
        %     cmap_hsv(:,2) = 1.2.*cmap_hsv(:,2); cmap_hsv(cmap_hsv > 1) = 1;
        %     cmap_sat = hsv2rgb(cmap_hsv);
        
        colormap(gca,cmap)
        
        clim(clims)
        
        % COLORBAR
        switch ax
            %         case {1,2}
            %             cbar = nph_colorbar;
            %             cbar.Ticks = cbarticks{ax};
            %             cbar.Position = cbar.Position .* [0.99 1 0.5 1];
            %             cbar.Label.String = ['\bf{' clabs{ax} '}'];
            %             cbar.Label.Rotation = 0;
            %             cbar.Label.VerticalAlignment = 'middle';
            % %             cbar.Label.HorizontalAlignment = 'left';
            %             cbar.Title.String = '   ms^{-1}';
            case {3,4}
                cbar = nph_colorbar;
                cbar.Ticks = cbarticks;
                cbar.Position = cbar.Position .* [1 1 3 1];
                cbar.Position = cbar.Position +  [cbar.Position(3)/3 0*cbar.Position(4) 0 0];
                %             cbar.Label.String = ['\bf{' clabs '}'];
                %             cbar.Label.Rotation = 0;
                %             cbar.Label.VerticalAlignment = 'middle';
                %             cbar.Label.HorizontalAlignment = 'left';
                drawnow;
                if strcmpi(ap,'p')
                    cbar.Label.String = clabs;
                    cbar.Label.HorizontalAlignment = 'center';
                    %                 cbar.Label.Position(1) = 0.1*cbar.Title.Position(1);
                    cbar.Label.Rotation = 90;
                    cbar.Label.FontSize = 0.8*fs;
                else
                    cbar.Title.String = clabs;
                end
        end
        
        
        
        % AXES LETTER
        switch ax
            case {1,2}
                hold on; nph_text([0.005 0.85],['(' alphabet(letternum(ax)+bump) ')'],'fontsize',           1.5*fs,'color',textcolor,'textborder',textbordercolor,'horizontalalignment','left','verticalalignment','middle');
                hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','color',textcolor,'textborder',textbordercolor,'horizontalalignment','left','verticalalignment','middle');
            case 3 % no idea why these are different
                hold on; nph_text([0.075 0.85],['(' alphabet(letternum(ax)+bump) ')'],'fontsize',           1.5*fs,'color',textcolor,'textborder',textbordercolor,'horizontalalignment','left','verticalalignment','middle');
                %             hold on; nph_text([0.075 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','color','k','textborder','w','horizontalalignment','left','verticalalignment','middle');
            case 4 % no idea why these are different
                hold on; nph_text([0.025 0.85],['(' alphabet(letternum(ax)+bump) ')'],'fontsize',           1.5*fs,'color',textcolor,'textborder',textbordercolor,'horizontalalignment','left','verticalalignment','middle');
                %             hold on; nph_text([0.025 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','color','k','textborder','w','horizontalalignment','left','verticalalignment','middle');
        end
        
        
        %     switch ax
        %         case {1,2}
        % %             text(datenum(2015,10,01),90,tits{ax},'fontsize',1*fs,'fontweight','bold','HorizontalAlignment','center','rotation',90)
        %             text(datenum(2016,4,15),107.25,tits{ax},'fontsize',1*fs,'fontweight','bold','HorizontalAlignment','left','rotation',0)
        %
        %         case {3,4}
        %             text(datenum(2019,02,01),104.75,tits{ax},'fontsize',1*fs,'fontweight','bold','HorizontalAlignment','left')
        %     end
        
        
        setfont(fs)
        
        set(gca,'linewi',1.5,'layer','top','tickdir','out')
        
        
    end
    
    
%     return
    
    % %% EXPORT??? ==============================================================
    
    drawnow;
    
    % savename = ['~/Desktop/' upper(site) '_TidesAllYearsComp_' ampphase];
    savename = ['~/Desktop/' upper(site) '_' num2str(tdh) 'hTideComp_' ampphase];
    
    disp(['Saving as "' savename '"...'])
    
    nph_saveas(gcf,savename,'png')
    
    disp('Saved.')
    
    
    
    close all hidden
    
end

return


%% JUST DO A COMP YEAR OF TIDAL AMPS AND PHASES AGAINST HEIGHT AND TIME

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig;

figpos([1 0.7])
% figpos([0.9 1])

%-------------------------------------------------------
vert_gap = 0.08;        horz_gap = 0.03;
lower_marg = 0.075;     upper_marg = 0.025;
left_marg = 0.065;      right_marg = 0.075;

rows = 3; cols = 9;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 20;

%--------------------------------------------------------

T = tiledlayout(rows,cols,'tilespacing','compact');

spannnyspan = [1 2];

axcounter = 0;

for td = 1:3
    
    % select tide:
    % td = 1;
    
    %     % bump for the axes labels
    %     bump = (td-1) * 4;
    
    %     % select amplitude or phase:
    %     ap = 'a'; % a | p
    
    %%%% First, assemble and over-interpolate the tidal amplitudes and phases:
    tdh = 24/td;
    
    %     TD.uA = Comp30.Tides.u.A(:,:,td);
    %     TD.uB = Comp30.Tides.u.B(:,:,td);
    %     TD.vA = Comp30.Tides.v.A(:,:,td);
    %     TD.vB = Comp30.Tides.v.B(:,:,td);
    
    %     %%%% Over-Interpolation?
%         odays = datenum(2016,01,01):datenum(2021,01,01);
        ozrange = 76:0.25:105;
    %     for uv = 'uv'
    %         for ab = 'AB'
    %             nanlocs = isnan(Comp30.wtime) | isnan(Comp30.walt) | isnan(TD.([uv ab]));
    %             %         F = griddedInterpolant({zrange,days},TD.([uv ab]),'linear','none');
    %             F = scatteredInterpolant(linearise(Comp30.walt(~nanlocs)),linearise(Comp30.wtime(~nanlocs)),linearise(TD.([uv ab])(~nanlocs)),'linear','none');
    %             TD.([uv ab]) = F({ozrange,odays});
    %         end
    %     end
    
    %     %%%% first, compute tidal amps and phases:
    %     umag = quadadd(TD.uA,TD.uB);
    %     vmag = quadadd(TD.vA,TD.vB);
    %     uph  = atan2(TD.uB,TD.uA);
    %     vph  = atan2(TD.vB,TD.vA);
    %
    %     % convert to local time:
    %     uph = wrapToPi(uph + ((HWD.Meta.TimeZone/24) * (2*pi)));
    %     vph = wrapToPi(vph + ((HWD.Meta.TimeZone/24) * (2*pi)));
    %     % convert to hours...
    %     uph = (uph./(2*pi)) .* tdh;
    %     vph = (vph./(2*pi)) .* tdh;
    %     % wrap:
    %     uph(uph < 0) = tdh + uph(uph < 0);
    %     vph(vph < 0) = tdh + vph(vph < 0);
    
    %%%% Now assemble our all-years composite...
    odaynumbers = 0:0.5:365; % and over-interpolate (for phase)
    TD.uAy = interp1(daynumbers,interp1(zrange,CompYear.Tides.u.A(:,:,td),ozrange)',odaynumbers)';
    TD.uBy = interp1(daynumbers,interp1(zrange,CompYear.Tides.u.B(:,:,td),ozrange)',odaynumbers)';
    TD.vAy = interp1(daynumbers,interp1(zrange,CompYear.Tides.v.A(:,:,td),ozrange)',odaynumbers)';
    TD.vBy = interp1(daynumbers,interp1(zrange,CompYear.Tides.v.B(:,:,td),ozrange)',odaynumbers)';
    % umagy = quadadd(CompYear.Tides.u.A(:,:,td),CompYear.Tides.u.B(:,:,td));
    % vmagy = quadadd(CompYear.Tides.v.A(:,:,td),CompYear.Tides.v.B(:,:,td));
    % uphy  = atan2(CompYear.Tides.u.B(:,:,td),CompYear.Tides.u.A(:,:,td));
    % vphy  = atan2(CompYear.Tides.v.B(:,:,td),CompYear.Tides.v.A(:,:,td));
    umagy = quadadd(TD.uAy,TD.uBy);
    vmagy = quadadd(TD.vAy,TD.vBy);
    uphy  = atan2(TD.uBy,TD.uAy);
    vphy  = atan2(TD.vBy,TD.vAy);
    % convert to local time:
    uphy = wrapToPi(uphy + ((HWD.Meta.TimeZone/24) * (2*pi)));
    vphy = wrapToPi(vphy + ((HWD.Meta.TimeZone/24) * (2*pi)));
    % convert to hours...
    uphy = (uphy./(2*pi)) .* tdh;
    vphy = (vphy./(2*pi)) .* tdh;
    % wrap:
    uphy(uphy < 0) = tdh + uphy(uphy < 0);
    vphy(vphy < 0) = tdh + vphy(vphy < 0);
    
    % =========================================================================
    
    % different limits for different tides...
    switch td
        case 1 % 24h tide
            ampclims = [0 20];
            ampclevs = 0:1:20;
            ampclines = 0:5:50;
            ampcbarticks = 0:5:60;
            ampcbarminorticks = ampclevs;
            phaseminorclines = 0:3:tdh;
            darkclines = [10 10];
        case 2 % 12h tide
            ampclims = [0 70];
            ampclevs = 0:5:100;
            ampclines = 0:20:100;
            ampcbarticks = 0:20:100;
            ampcbarminorticks = ampclevs;
            phaseminorclines = 0:1:tdh;
            darkclines = [40 40];
        case 3 % 8h tide
            ampclims = [0 15];
            ampclevs = 0:1:20;
            ampclines = 0:5:50;
            ampcbarticks = 0:5:60;
            ampcbarminorticks = ampclevs;
            phaseminorclines = 0:1:tdh;
            darkclines = [10 10];
        case 4 % 6h tide
            ampclims = [0 5];
            ampclevs = 0:0.5:20;
            ampclines = 0:2.5:10;
            ampcbarticks = 0:1:10;
            ampcbarminorticks = ampclevs;
            phaseminorclines = 0:1:tdh;
            darkclines = [2.5 2.5];
    end
    
%     switch td
%         case {1,2,3}
            phaseclims = [0 tdh];
            phaseclevs = 0:1:tdh;
            phaseclines = tdh .* [0.25 0.5 0.75];
            phasecbarticks = 0:(tdh/4):(tdh);
            phasecbarminorticks = 0:(tdh/8):(tdh);
%         case 4
%             phaseclims = [0 tdh];
%             phaseclevs = 0:1:tdh;
%             phaseclines = tdh .* [0.5 0.5];
%             phasecbarticks = 0:(tdh/2):(tdh);
%             phasecbarminorticks = 0:(tdh/4):(tdh);
%     end
    
    AP = {'a','a','p','p'};
    
    spans = {spannnyspan,spannnyspan,spannnyspan,spannnyspan};
    
    TBP = {umagy,vmagy,uphy,vphy};
    
    tits = {...
        'ZONAL COMPONENT',...
        'MERIDIONAL COMPONENT',...
        'ZONAL',...
        'MERID.'};
    
%     letternum = [1 3 2 4];
    
%     axs = {1:4, 5:8, 9:12};
    Axes = struct;
    
    for ax = 1:4
        
        axcounter = axcounter + 1;
        
        % gap in the middle of the columns
        switch ax
            case 3
                nexttile([1 1])
                set(gca,'color','none','xcolor','none','ycolor','none');
        end
        
        %     subplot(rows,cols,axs{ax})
        nexttile(spans{ax})
        
        axx = gca;
        
        Axes(ax).axx = gca;
        
        %         % choose wot 2 plot n wot not
        %         switch ax
        %             case {1,2}
        %                 %             X = Comp.wtime; Y = Comp.walt;
        %                 [X,Y] = meshgrid(odays,ozrange);
        %                 smoo1 = [0.25 0.75]; smoo2 = [0.25 0.75];
        %             case {3,4}
        [X,Y] = meshgrid(datenum(2019,00,odaynumbers),ozrange);
        smoo1 = [0.25 0.75]; smoo2 = [0.25 0.75];
        %         end
        
        ap = AP{ax};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SKIPPY SKIP FOR TESTING
%         if strcmp(ap,'p')
%             continue
%         end
%         if td ~= 1
%             continue
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        switch ap
            case 'a'
                clims = ampclims;
                clevs = ampclevs;
                clines = ampclines;
                cbarticks = ampcbarticks;
                cbarminorticks = ampcbarminorticks;
                clabs = ['   \bf{ ms^{-1} }'];
                cmap = nph_saturate(cbrew('nph_modspectral',32),1.2);
                ampphase = 'amp';
                textcolor = 'w';
                textbordercolor = 'k';
            case 'p'
                clims = phaseclims;
                clevs = phaseclevs;
                clines = phaseclines;
                minorclines = phaseminorclines;
                cbarticks = phasecbarticks;
                cbarminorticks = phasecbarminorticks;
                clabs = '\bf{Hours}';
                cmap = nph_saturate(cbrew('nph_CyclicRainbow',24),1);
                ampphase = 'phase';
                textcolor = 'w';
                textbordercolor = 'k';
        end
        
        
        
        % allll the contours...
        tbp = TBP{ax};
        nanlocs = isnan(tbp); smoo = movmean(tbp,30,2,'omitnan');
        tbp(nanlocs) = smoo(nanlocs);
        switch ap
            case 'a'
                tbp = imgaussfilt(imgaussfilt(tbp,smoo1),smoo2);
                if td == 4
                    tbp = imgaussfilt(imgaussfilt(tbp,smoo1),smoo2);
                end
                % filled contours:
                hold on; contourf(X,Y,tbp,clevs,'edgecolor','none');
                % contour lines:
                hold on; [C,h] = contour(X,Y,tbp,clines,'edgecolor',rgbtrip(.25),'linewi',1.25);
                clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
                % and bold zero line
                hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','w','linewi',4);
                clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000,'color','w','fontweight','bold');
                hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linewi',2);
                clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000);
                
                % and dark contour lines:
                % contour lines:
                hold on; [C,h] = contour(X,Y,tbp,darkclines,'edgecolor',rgbtrip(.25),'linewi',1.25);
                clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.25),'fontweight','bold');
                
                
            case 'p'
                % pcolor instead for phase?
                hold on; pcolor(X,Y,tbp); shat;
                
                %%%% shall we try some contours for phase too?
                %%%% just try a few...
%                 if td == 4
%                     cospart = imgaussfilt(cos(twopi.*(tbp/tdh)),smoo1);
%                     sinpart = imgaussfilt(sin(twopi.*(tbp/tdh)),smoo1);
%                     tbp     = atan2(cospart,sinpart);
%                     % convert to hours...
%                     tbp = (tbp./(2*pi)) .* tdh;
%                     % wrap:
%                     tbp(tbp < 0) = tdh + tbp(tbp < 0);
%                 end
                
                % and some minor clines:
                minorclines = minorclines(~ismember(minorclines,clines));
                hold on; [C,h] = contour(X,Y,tbp,minorclines,'edgecolor',rgbtrip(.5),'linewi',0.5);
                
                % first zero phase:
                tbp2 = tbp; tbp2(tbp > tdh/2) = -tdh + tbp2(tbp > tdh/2);
                nanmask = abs(diff(tbp2([1 1:end],:),[],1)) > tdh/2 | abs(diff(tbp2(:,[1 1:end]),[],2)) > tdh/2;
                tbp2(nanmask) = NaN;
                hold on; [C,h] = contour(X,Y,tbp2,[0 0],'edgecolor',rgbtrip(.25),'linewi',1.25);
                clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
                
                % before you do the others, exclude any wraparound bits:
                nanmask = abs(diff(tbp([1 1:end],:),[],1)) > tdh/2 | abs(diff(tbp(:,[1 1:end]),[],2)) > tdh/2;
                tbp(nanmask) = NaN;
                hold on; [C,h] = contour(X,Y,tbp,clines,'edgecolor',rgbtrip(.25),'linewi',1.25);
                clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
                
        end
        
        
        % LIMITS
        ylim([77 103])
        axx.YTick = 70:5:110;
        axx.YMinorTick = 'on';
        %     switch ax
        %         case {1,2}
        axx.YAxis.MinorTickValues = 70:1:110;
        %         case {3,4}
        %             axx.YAxis.MinorTickValues = 70:110;
        %     end
        
        
        %         switch ax
        %             case {1,2}
        %                 xlim([datenum(2016,02,15) datenum(2020,11,15)])
        %                 % Just show the Januarys as major ticks
        %                 axx.XTick = datenum(years,01,01);
        %                 axx.XMinorTick = 'on';
        %                 axx.XAxis.MinorTickValues = datenum(2016,1:(length(years)*12),01);
        %                 datetick('x',' ','keepticks','keeplimits')
        %                 %             axx.XTickLabel = {};
        %                 % other months as ticks using text:
        %                 xtix = datenum(years(1),1:(length(years)*12),15);
        %                 for xt = xtix
        %                     if inrange(xt,xlim) && ~any(xt == axx.XTick)
        %                         mn = monthname(month(xt));
        %                         text(xt,75.5,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
        %                     end
        %                 end
        %                 if ax == 2
        %                     % Years as lower bold numbers:
        %                     xtix = datenum(years,07,01);
        %                     for xt = xtix
        %                         if inrange(xt,xlim) && ~any(xt == axx.XTick)
        %                             yrstr = datestr(xt,'yyyy');
        %                             text(xt,72.5,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
        %                         end
        %                     end
        %                 end
        %
        %             case {3,4}
        xlim(datenum(2019,[1 13],01))
        % Just show the Januarys as major ticks
        axx.XTick = datenum(2019,1:13,01);
        axx.XMinorTick = 'on';
        axx.XAxis.MinorTickValues = datenum(2016,1:(length(years)*12),01);
        datetick('x',' ','keepticks','keeplimits')
        % other months as ticks using text:
        xtix = datenum(years(1),1:(length(years)*12),16);
        for xt = xtix
            if inrange(xt,xlim) && ~any(xt == axx.XTick)
                mn = monthname(month(xt));
                text(xt,76,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
            end
        end
        %             axx.XTick = datenum(2019,1:13,01);
        %             datetick('x','m','keepticks','keeplimits')
        %             xlim(datenum(2019,[1 13],01))
%         if ax == 4
%             xlabel('Average Year','fontsize',0.8*fs,'fontweight','bold')
%         end
        %         end
        
        
        switch ax
            case 1
                ylabel('Altitude (km)')
        end
        
        % line around the axis
        hold on; plot(axx.XLim([1 2 2]),axx.YLim([2 2 1]),'color','k','linewi',1.5)
        
        drawnow;
        
        % COLORI
        %     cmap = cbrew('RdBu',length(clevs{ax}));
        %     cmap = cbrew('nph_RdYlBu',length(clevs{ax}));
        %     cmap = cbrew('nph_BuOr',length(clevs{ax}));
        %     cmap = cbrew('BrBG',length(clevs{ax}));
        %     cmap = flipud(cbrew('PiYG',length(clevs{ax})));
        %     cmap = cbrew('nph_Spectral',length(clevs{ax}));
        %     cmap = cbrew('nph_modspectral',length(clevs{ax}));
        %     cmap = nph_saturate(cbrew('nph_modspectral',18),1.2);
        %     cmap = cbrew('nph_RdBuPastel',length(clevs{ax}));
        %     colormap(gca,cmap)
        
        % sort out double white in the middle:
        %     cmap = [cmap(1:floor(length(clevs{ax})/2),:) ; [1 1 1] ; cmap(floor(length(clevs{ax})/2)+1:end,:)];
        
        % saturate?
        %     cmap_hsv = rgb2hsv(cmap);
        %     cmap_hsv(:,2) = 1.2.*cmap_hsv(:,2); cmap_hsv(cmap_hsv > 1) = 1;
        %     cmap_sat = hsv2rgb(cmap_hsv);
        
        colormap(gca,cmap)
        
        clim(clims)
        
        drawnow;
        
        % COLORBAR
        switch ax
            %         case {1,2}
            %             cbar = nph_colorbar;
            %             cbar.Ticks = cbarticks{ax};
            %             cbar.Position = cbar.Position .* [0.99 1 0.5 1];
            %             cbar.Label.String = ['\bf{' clabs{ax} '}'];
            %             cbar.Label.Rotation = 0;
            %             cbar.Label.VerticalAlignment = 'middle';
            % %             cbar.Label.HorizontalAlignment = 'left';
            %             cbar.Title.String = '   ms^{-1}';
            case {2,4}
                cbar = nph_colorbar;
                cbar.Ticks = cbarticks;
                cbar.Position = cbar.Position .* [1 1 3 1];
                if ax == 4
                    cbar.Position = cbar.Position +  [cbar.Position(3)/3 0*cbar.Position(4) 0 0];
                end
                if ax == 2 && td > 1
                    cbar.Position([1 3]) = prevcbar.Position([1 3]);
                end
                cbar.LineWidth = 1.5;
                cbar.TickDirection = 'out';
                cbar.Ruler.MinorTick = 'on';
                switch ap
                    case 'a'
                        cbar.Ruler.MinorTickValues = ampcbarminorticks;                        
                    case 'p'
                        cbar.Ruler.MinorTickValues = sort([cbar.Ticks cbar.Ticks+mean(diff(cbar.Ticks)/2)]);
                end
                %             cbar.Label.String = ['\bf{' clabs '}'];
                %             cbar.Label.Rotation = 0;
                %             cbar.Label.VerticalAlignment = 'middle';
                %             cbar.Label.HorizontalAlignment = 'left';
                drawnow;
%                 if strcmpi(ap,'p')
%                     cbar.Label.String = clabs;
                    cbar.Label.HorizontalAlignment = 'center';
                    %                 cbar.Label.Position(1) = 0.1*cbar.Title.Position(1);
                    cbar.Label.Rotation = 90;
                    cbar.Label.FontSize = 0.8*fs;
%                 else
                    cbar.Title.String = clabs;
%                 end
                if ax == 2
                    prevcbar = cbar;
                end
                
        end
        
        drawnow;
        
        %         % AXES LETTER
        %         switch ax
        %             case {1,2}
%         hold on; nph_text([0.005 0.85],['(' alphabet(axcounter) ')'],'fontsize',           1.5*fs,'color',textcolor,'textborder',textbordercolor,'horizontalalignment','left','verticalalignment','middle');
        %                 hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','color',textcolor,'textborder',textbordercolor,'horizontalalignment','left','verticalalignment','middle');
        %             case 3 % no idea why these are different
        %                 hold on; nph_text([0.075 0.85],['(' alphabet(letternum(ax)+bump) ')'],'fontsize',           1.5*fs,'color',textcolor,'textborder',textbordercolor,'horizontalalignment','left','verticalalignment','middle');
        %                 %             hold on; nph_text([0.075 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','color','k','textborder','w','horizontalalignment','left','verticalalignment','middle');
        %             case 4 % no idea why these are different
        %                 hold on; nph_text([0.025 0.85],['(' alphabet(letternum(ax)+bump) ')'],'fontsize',           1.5*fs,'color',textcolor,'textborder',textbordercolor,'horizontalalignment','left','verticalalignment','middle');
        %                 %             hold on; nph_text([0.025 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','color','k','textborder','w','horizontalalignment','left','verticalalignment','middle');
        %         end
        
        
        %     switch ax
        %         case {1,2}
        % %             text(datenum(2015,10,01),90,tits{ax},'fontsize',1*fs,'fontweight','bold','HorizontalAlignment','center','rotation',90)
        %             text(datenum(2016,4,15),107.25,tits{ax},'fontsize',1*fs,'fontweight','bold','HorizontalAlignment','left','rotation',0)
        %
        %         case {3,4}
        %             text(datenum(2019,02,01),104.75,tits{ax},'fontsize',1*fs,'fontweight','bold','HorizontalAlignment','left')
        %     end
        
        
        setfont(fs)
        
        set(gca,'linewi',1.5,'layer','top','tickdir','out')
        
        
        
% % % %         % gap at the end of the columns
% % % %         switch ax
% % % %             case 4
% % % %                 nexttile([1 1])
% % % %                 set(gca,'color','none','xcolor','none','ycolor','none');
% % % %         end
        
        
    end
    
    % LABELS AFTER
    for ax = 1:4
        axes(Axes(ax).axx);
        bump = (td-1) * 4;
%         hold on; nph_text([0.005 0.85],['(' alphabet(ax+bump) ')'],'fontsize',           1.5*fs,'color',textcolor,'textborder',textbordercolor,'horizontalalignment','left','verticalalignment','middle');
        hold on; nph_text([0.005 0.85],['(' alphabet(ax+bump) ')'],'fontsize',           1.5*fs,'color','k','textborder','w','horizontalalignment','left','verticalalignment','middle');
    end
    
    


    
end


return


%% EXPORT??? ==============================================================

% drawnow;

% savename = ['~/Desktop/' upper(site) '_TidesAllYearsComp_' ampphase];
savename = ['~/Desktop/' upper(site) '_TidesCompYearAmpPhase'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')



% close all hidden








return







%% PLOT OF ANNUAL AND SEMIANNUAL CYCLES AND PHASES


figure; hold all; whitefig; figpos([0.75 0.75])

%-------------------------------------------------------
vert_gap = 0.08;        horz_gap = 0.03;
lower_marg = 0.075;     upper_marg = 0.025;
left_marg = 0.065;      right_marg = 0.075;

rows = 2; cols = 8;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 20;

%--------------------------------------------------------

part = 2;

switch part
    case 1
        axrange = 1:3;
    case 2
        axrange = 4:6;
end

T = tiledlayout(2,3,'tilespacing','compact');

spans = {...
    [1 1],...
    [1 1],...
    [1 1],...
    [1 1],...
    [1 1],...
    [1 1]};

xlims = {...
    [0 40],...
    [0 20],...
    [0 20],...
    datenum([2019 2020],01,01),...
    datenum([2019 2020],01,01),...
    datenum([2019 2020],01,01)};

xticks = {...
    [0:10:40],...
    [0:5:20],...
    [0:5:20],...
    datenum(2019,1:13,01),...
    datenum(2019,1:13,01),...
    datenum(2019,1:13,01)};

xminorticks = {...
    [0:2.5:40],...
    [0:1:20],...
    [0:1:20],...
    [],...
    [],...
    []};

cmap = nph_saturate(cbrew('nph_modspectral',22),2);

tidecolors = {...
    cmap(3,:),...
    cmap(16,:),...
    mcolor(3),...
    mcolor(4)};
% tidecolors = {...
%     mcolor(1),...
%     mcolor(2),...
%     mcolor(3),...
%     mcolor(4)};

tits = {...
    'ANNUAL MEAN',...
    'ANNUAL CYCLE',...
    'SEMIANNUAL CYCLE',...
    ' ',...
    'ANNUAL MAXIMA',...
    'SEMIANNUAL MAXIMA'};


for ax = axrange
    
    nexttile(spans{ax})
    axx = gca;
    
    %%%% TICKS AND LIMITS
    xlim(xlims{ax})
    xtick(xticks{ax})
    xminortick(xminorticks{ax})
    
    ylim([76 103])
    ytick(70:5:110);
    yminortick(70:110);
        
    switch ax
        case 1
            ylabel('Altitude (km)')
            xl = xlabel('Wind Speed (ms^-^1)');
            xl.FontSize = 0.9*fs;
%             xl.Position(2) = xl.Position(2) - 1;
        case {1,2,3}
            xl = xlabel('Wind Speed (ms^-^1)');
%             xl.Position(2) = xl.Position(2) - 1;
        case {4,5,6}
            xl = xlabel('Month');
%             if ax == 5, ylabel('Altitude (km)'); end
    end
    
    if part == 2 && ax == 4
        ylabel('Altitude (km)')
    end
    
    switch ax
        case {4,5,6}
            datetick('x','      m','keepticks','keeplimits')
            axx.XTickLabel(13,:) = repmat(' ',1,length(axx.XTickLabel(13,:)));
    end
    
%     %%%% SKIP AXIS 4
%     if ax == 4
%         axx.Color = 'none';
%         axx.XColor = 'none';
%         axx.YColor = 'none';
%         title(' ')
%         continue
%     end
    
    %%%% determine what we're gonna plot:
    linestyle.u = {'linewi',3,'linest','-'};
    linestyle.v = {'linewi',3,'linest','--'};
    patchlinestyle.u = 'none';
    patchlinestyle.v = ':';
    
    switch ax
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% TIDAL CYCLE AMPLITUDES
        case {1,2,3}
            %             %%%% PLOT WIND
            %             hold on; plot(Cycles.u(:,ax),zrange,zonallinestyle{:}); grid on;
            %             hold on; plot(Cycles.v(:,ax),zrange,meridlinestyle{:}); grid on;
            
            %%%% PLOT FILLED TIDAL STDs:
            for td = [3 2 1]
                for uv = 'uv'
                    hold on; F = fill(...
                        [Cycles.Tides.(uv)(:,ax,td)+Cycles.STDs.Tides.(uv)(:,ax,td) ; ...
                        reverse(Cycles.Tides.(uv)(:,ax,td)-Cycles.STDs.Tides.(uv)(:,ax,td))],...
                        [zrange reverse(zrange)]','k'); grid on;
                    F.LineStyle = patchlinestyle.(uv);
                    F.LineWidth = 1;
                    F.EdgeColor = tidecolors{td};
                    F.FaceColor = tidecolors{td};
                    F.FaceAlpha = 0.15;
                end
            end
            %%%% THEN PLOT LINES ON TOP
            for td = [3 2 1]
                for uv = 'uv'
%                     hold on; plot(Cycles.Tides.(uv)(:,ax,td),zrange,linestyle.(uv){:},'color','w','linewi',3.5); grid on;
                    hold on; plot(Cycles.Tides.(uv)(:,ax,td),zrange,linestyle.(uv){:},'color',tidecolors{td}); grid on;
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% TIDAL CYCLE PHASES
        case {4,5,6}
            %%%% need to do +/-pi duplicates for the (semi)annual cycle(s)
            if ax == 5, shifts = [0 2*pi]; end
            if ax == 6, shifts = [0 pi 2*pi]; end
            
            % Filled STDs, including the shifted ones:
            for td = [1 3 2]
                for uv = 'uv'
                    % extract phase and STDs
                    ph = unwrap(Cycles.Tides.([uv 'ph'])(:,ax-3,td));
                    ph_std = Cycles.STDs.Tides.([uv 'ph'])(:,ax-3,td);
                    % convert to time of year in 2019, why not:
                    ph = datenum(2019,01,01) + ((ph./(2*pi)) * 365);
                    % convert to decimal days:
                    ph_std = ((ph_std./(2*pi)) * 365);
                    for s = 1:length(shifts)
                        shift = ((shifts(s)/(2*pi)) * 365);
                        hold on; F = fill(shift+[ph+ph_std ; reverse(ph-ph_std)],[zrange reverse(zrange)]','k'); grid on;
                        F.LineStyle = patchlinestyle.(uv);
                        F.LineWidth = 1;
                        F.EdgeColor = tidecolors{td};
                        F.FaceColor = tidecolors{td};
                        F.FaceAlpha = 0.15;
                    end
                end
            end
            % and now the lines:
            for td = [1 3 2]
                for uv = 'uv'
                    % extract phase and STDs
                    ph = unwrap(Cycles.Tides.([uv 'ph'])(:,ax-3,td));
                    ph_std = Cycles.STDs.Tides.([uv 'ph'])(:,ax-3,td);
                    % convert to time of year in 2019, why not:
                    ph = datenum(2019,01,01) + ((ph./(2*pi)) * 365);
                    % convert to decimal days:
                    ph_std = ((ph_std./(2*pi)) * 365);
                    for s = 1:length(shifts)
                        shift = ((shifts(s)/(2*pi)) * 365);
%                         hold on; plot(ph+shift,zrange,linestyle.(uv){:},'color','w','linewi',3.5); grid on;
                        hold on; plot(ph+shift,zrange,linestyle.(uv){:},'color',tidecolors{td}); grid on;
                    end
                end
            end
            
    end
    
    
    %%%% line around the axes
    hold on; plot(axx.XLim([1 2 2 1 1]),axx.YLim([1 1 2 2 1]),'color','k','linewi',1.5)
    
    %%%% TITLES
    tit = title(tits{ax});
    tit.Position(2) = 1.005 * tit.Position(2);
    tit.FontSize = 0.9*fs;
    
    %%%% labels
    switch ax
        case {1,2,3}
            hold on; nph_text([0 1.03],['(' alphabet(ax) ')'],'textborder','w','fontsize',1.75*fs,'fontweight','normal','horizontalalignment','left');
        case {4,5,6}
            hold on; nph_text([0 1.03],['(' alphabet(ax-1) ')'],'textborder','w','fontsize',1.75*fs,'fontweight','normal','horizontalalignment','left');
    end
    
    set(gca,'linewi',1.5,'tickdir','out');
    
    setfont(fs)
    
end

% return


%%%% EXPORT??? ==============================================================

% savename = ['~/Desktop/' upper(site) '_TidesAllYearsComp_' ampphase];
savename = ['~/Desktop/' upper(site) '_TidalAnnualCycles_part' num2str(part)];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')

close all hidden

% end % next part




return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% DISPLAY AVERAGE AMPS AND PHASES OF TIDES

% zlims = [80 85 90 95 ; 85 90 95 100];
zlims = [90 ; 100];

cycleperiods = 365 .* [1 0.5];
tidalperiods = [24 12 8 6] ./ 24;

% do the same for tidal periods...
for td = 1:length(tidalperiods(1:3))
% for td = 2
    disp(['~~~~~' num2str(tidalperiods(td)*24) 'h tide ~~~~~'])
    
    disp({[num2str(zlims(1,a)) '-' num2str(zlims(2,a)) 'km']})
    disp({'         ','Annual Mean','Annual','Semiannual'})
    disp({'         ','Zon, Mer','Zon, Mer','Zon, Mer'})
    
    for a = 1:size(zlims,2)
        zinds = inrange(zrange,zlims([1 2],a));
        [~,Fu] = nph_sinefit(daynumbers,quadadd(nanmean(CompYear.Tides.u.A(zinds,:,td),1),nanmean(CompYear.Tides.u.B(zinds,:,td),1)),cycleperiods);
        [~,Fv] = nph_sinefit(daynumbers,quadadd(nanmean(CompYear.Tides.v.A(zinds,:,td),1),nanmean(CompYear.Tides.v.B(zinds,:,td),1)),cycleperiods);
        disp({...
            'Amplitude', ...
            [num2str(round(Fu(2,3),1)) ', ' num2str(round(Fv(2,3),1))], ...
            [num2str(round(quadadd(Fu(1,1),Fu(1,2)),1)) ', ' num2str(round(quadadd(Fv(1,1),Fv(1,2)),1))], ...
            [num2str(round(quadadd(Fu(2,1),Fu(2,2)),1)) ', ' num2str(round(quadadd(Fv(2,1),Fv(2,2)),1))], ...
            })
        
        uph = unwrap(atan2(Fu(:,1),Fu(:,2)));
        uph = datenum(2019,01,01) + ((uph./(2*pi)) * 365);
        
        vph = unwrap(atan2(Fv(:,1),Fv(:,2)));
        vph = datenum(2019,01,01) + ((vph./(2*pi)) * 365);
        
%         phdate1_u = datestr(datenum(2019,00,(atan2(Fu(1,1),Fu(1,2))/twopi)*365),'dd mmm');
%         phdate1_v = datestr(datenum(2019,00,(atan2(Fv(1,1),Fv(1,2))/twopi)*365),'dd mmm');
%         phdate2_u = datestr(datenum(2019,00,pi/2+(atan2(Fu(2,1),Fu(2,2))/twopi)*365),'dd mmm');
%         phdate2_v = datestr(datenum(2019,00,pi/2+(atan2(Fv(2,1),Fv(2,2))/twopi)*365),'dd mmm');
%         
        disp({...
            'Phase', ...
            ['    ,    '], ...
            [datestr(uph(1),'dd mmm') ', ' datestr(vph(1),'dd mmm')], ...
            [datestr(uph(2),'dd mmm') ', ' datestr(vph(2),'dd mmm')], ...
            })
%         return
        %             Cycles.Tides.u(z,1,td) = Fu(2,3);
        %             Cycles.Tides.v(z,1,td) = Fv(2,3);
        %             Cycles.Tides.u(z,2:3,td) = quadadd(Fu(1:2,1),Fu(1:2,2))';
        %             Cycles.Tides.v(z,2:3,td) = quadadd(Fv(1:2,1),Fv(1:2,2))';
        %             Cycles.Tides.uph(z,2:3,td) = atan2(Fu([2 1],1),Fu([2 1],2))';
        %             Cycles.Tides.vph(z,2:3,td) = atan2(Fv([2 1],1),Fv([2 1],2))';
    end
end

    






%% OLD - PLOT TIDAL AMPLITUDES AGAINST TIME AND HEIGHT
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure; hold all; whitefig; figpos([1 0.475])
%
% %-------------------------------------------------------
% vert_gap = 0.08;        horz_gap = 0.03;
% lower_marg = 0.075;     upper_marg = 0.025;
% left_marg = 0.065;      right_marg = 0.075;
%
% rows = 2; cols = 8;
%
% subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);
%
% %--------------------------------------------------------
%
% fs = 18;
%
% td = 2;
%
% alldays = datenum(min(years),01,01):1:datenum(max(years)+1,01,01);
% % uA = interp1(days,Comp30.Tides.u.A(:,:,td)',alldays)';
% % uB = interp1(days,Comp30.Tides.u.B(:,:,td)',alldays)';
% % vA = interp1(days,Comp30.Tides.v.A(:,:,td)',alldays)';
% % vB = interp1(days,Comp30.Tides.v.B(:,:,td)',alldays)';
% % uamp = quadadd(uA,uB);
% % vamp = quadadd(vA,vB);
% uamp = quadadd(Comp30.Tides.u.A(:,:,td),Comp30.Tides.u.A(:,:,td));
% vamp = quadadd(Comp30.Tides.v.B(:,:,td),Comp30.Tides.v.B(:,:,td));
%
%
% T = tiledlayout(2,6,'tilespacing','compact');
%
% spans = {...
%     [1 5],...
%     [1 5],...
%     [1 1],...
%     [1 1]};
%
% % axs = {...
% %     [1  6],...
% %     [9 14],...
% %     [7 8],...
% %     [15 16]};
%
% TBP = {...
%     uamp,...
%     vamp,...
%     CompYear.u,...
%     CompYear.v};
%
% clims = {...
%     [0 50],...
%     [0 20],...
%     [-50 50],...
%     [-20 20]};
%
% clevs = {...
%     [-50:5:50],...
%     [-20:2.5:20],...
%     [-50:5:50],...
%     [-20:2.5:20]};
%
% clabs = {...
%     'u','v','u','v'};
%
% cbarticks = {...
%     -60:20:60,...
%     -20:10:20,...
%     -60:20:60,...
%     -20:10:20};
%
% tits = {...
%     'ZONAL WIND',...
%     'MERIDIONAL WIND',...
%     'ZONAL WIND',...
%     'MERIDIONAL WIND'};
%
%
% % also contour lines:
% clines = {...
%     [-40 -20 20 40],...
%     [-20 -10 10 20],...
%     [-40 -20 20 40],...
%     [-20 -10 10 20]};
%
%
% for ax = [1 3 2 4]
%
% %     subplot(rows,cols,axs{ax})
%     nexttile(spans{ax})
%
%     if ax == 3, continue; end
%     if ax == 4, continue; end
%
%     axx = gca;
%
%     % choose wot 2 plot n wot not
%     switch ax
%         case {1,2}
%             X = Comp30.wtime; Y = Comp30.walt;
%             smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
%         case {3,4}
%             [X,Y] = meshgrid(datenum(2019,00,days),nanmean(CompYear.walt,2));
%             smoo1 = [0.25 0.75]; smoo2 = [0.25 0.75];
%     end
%
%     % allll the contours...
%     tbp = TBP{ax}; nanlocs = isnan(tbp); tbp(nanlocs) = nanmean(tbp(:));
%     tbp = imgaussfilt(imgaussfilt(tbp,smoo1),smoo2);
%     tbp(nanlocs) = NaN;
%     % filled contours:
%     hold on; contourf(X,Y,tbp,clevs{ax},'edgecolor','none');
%     % contour lines:
%     hold on; [C,h] = contour(X,Y,tbp,clines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
%     clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
%     % and bold zero line
%     hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','w','linewi',4);
%     clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000,'color','w','fontweight','bold');
%     hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linewi',2);
%     clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000);
%
%
%     % LIMITS
%     ylim([77 103])
%     axx.YTick = 70:5:110;
%     axx.YMinorTick = 'on';
% %     switch ax
% %         case {1,2}
%             axx.YAxis.MinorTickValues = 70:1:110;
% %         case {3,4}
% %             axx.YAxis.MinorTickValues = 70:110;
% %     end
%
%
%     switch ax
%         case {1,2}
%
%             xlim([datenum(2016,02,15) datenum(2020,11,15)])
%             % Just show the Januarys as major ticks
%             axx.XTick = datenum(years,01,01);
%             axx.XMinorTick = 'on';
%             axx.XAxis.MinorTickValues = datenum(2016,1:(length(years)*12),01);
%             datetick('x','m','keepticks','keeplimits')
%             % other months as ticks using text:
%             xtix = datenum(years(1),1:(length(years)*12),01);
%             for xt = xtix
%                if inrange(xt,xlim) && ~any(xt == axx.XTick)
%                    mn = monthname(month(xt));
%                    text(xt,75.5,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
%                end
%             end
%             if ax == 2
%             % Years as lower bold numbers:
%             xtix = datenum(years,07,01);
%             for xt = xtix
%                if inrange(xt,xlim) && ~any(xt == axx.XTick)
%                    yrstr = datestr(xt,'yyyy');
%                    text(xt,72.5,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
%                end
%             end
%             end
%         case {3,4}
%             xlim(datenum(2019,[1 13],01))
%             % Just show the Januarys as major ticks
%             axx.XTick = datenum(2019,1:13,01);
%             axx.XMinorTick = 'on';
%             axx.XAxis.MinorTickValues = datenum(2016,1:(length(years)*12),01);
%             datetick('x','m','keepticks','keeplimits')
%             % other months as ticks using text:
%             xtix = datenum(years(1),1:(length(years)*12),01);
%             for xt = xtix
%                if inrange(xt,xlim) && ~any(xt == axx.XTick)
%                    mn = monthname(month(xt));
%                    text(xt,75.5,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
%                end
%             end
% %             axx.XTick = datenum(2019,1:13,01);
% %             datetick('x','m','keepticks','keeplimits')
% %             xlim(datenum(2019,[1 13],01))
%             if ax == 4
%             xlabel('Composite Year','fontsize',0.8*fs,'fontweight','bold')
%             end
%     end
%
%
%     switch ax
%         case {1,2}
%             ylabel('Altitude (km)')
%     end
%
%     % line around the axis
%     hold on; plot(axx.XLim([1 2 2]),axx.YLim([2 2 1]),'color','k','linewi',1.5)
%
%     drawnow;
%
%     % COLORI
% %     cmap = cbrew('RdBu',length(clevs{ax}));
% %     cmap = cbrew('nph_RdYlBu',length(clevs{ax}));
% %     cmap = cbrew('nph_BuOr',length(clevs{ax}));
% %     cmap = cbrew('BrBG',length(clevs{ax}));
% %     cmap = flipud(cbrew('PiYG',length(clevs{ax})));
% %     cmap = cbrew('nph_Spectral',length(clevs{ax}));
% %     cmap = cbrew('nph_modspectral',length(clevs{ax}));
%     cmap = nph_saturate(cbrew('nph_modspectral',18),1.2);
% %     cmap = cbrew('nph_RdBuPastel',length(clevs{ax}));
% %     colormap(gca,cmap)
%
%     % sort out double white in the middle:
% %     cmap = [cmap(1:floor(length(clevs{ax})/2),:) ; [1 1 1] ; cmap(floor(length(clevs{ax})/2)+1:end,:)];
%
%     % saturate?
% %     cmap_hsv = rgb2hsv(cmap);
% %     cmap_hsv(:,2) = 1.2.*cmap_hsv(:,2); cmap_hsv(cmap_hsv > 1) = 1;
% %     cmap_sat = hsv2rgb(cmap_hsv);
%
%     colormap(gca,cmap)
%
%     clim(clims{ax})
%
%     % COLORBAR
%     switch ax
% %         case {1,2}
% %             cbar = nph_colorbar;
% %             cbar.Ticks = cbarticks{ax};
% %             cbar.Position = cbar.Position .* [0.99 1 0.5 1];
% %             cbar.Label.String = ['\bf{' clabs{ax} '}'];
% %             cbar.Label.Rotation = 0;
% %             cbar.Label.VerticalAlignment = 'middle';
% % %             cbar.Label.HorizontalAlignment = 'left';
% %             cbar.Title.String = '   ms^{-1}';
%        case {3,4}
%             cbar = nph_colorbar;
%             cbar.Ticks = cbarticks{ax};
%             cbar.Position = cbar.Position .* [1 1 3 1];
%             cbar.Position = cbar.Position +  [cbar.Position(3)/3 0*cbar.Position(4) 0 0];
%             cbar.Label.String = ['\bf{' clabs{ax} '}'];
%             cbar.Label.Rotation = 0;
%             cbar.Label.VerticalAlignment = 'middle';
% %             cbar.Label.HorizontalAlignment = 'left';
%             cbar.Title.String = '   ms^{-1}';
%     end
%
%
%
%     % AXES LETTER
%     switch ax
%         case {1,2}
%             hold on; nph_text([0.005 0.85],['(' alphabet(ax) ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
%             hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
%         case 3 % no idea why these are different
%             hold on; nph_text([0.075 0.85],['(' alphabet(ax) ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
%         case 4 % no idea why these are different
%             hold on; nph_text([0.025 0.85],['(' alphabet(ax) ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
%
%             %             hold on; nph_text([2*0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
%
%     end
%
%
% %     switch ax
% %         case {1,2}
% % %             text(datenum(2015,10,01),90,tits{ax},'fontsize',1*fs,'fontweight','bold','HorizontalAlignment','center','rotation',90)
% %             text(datenum(2016,4,15),107.25,tits{ax},'fontsize',1*fs,'fontweight','bold','HorizontalAlignment','left','rotation',0)
% %
% %         case {3,4}
% %             text(datenum(2019,02,01),104.75,tits{ax},'fontsize',1*fs,'fontweight','bold','HorizontalAlignment','left')
% %     end
%
%
%     setfont(fs)
%
%     set(gca,'linewi',1.5,'layer','top','tickdir','out')
%
%
% end
%
%
%
% return
%
% %% PLOT OF TIDAL AMPLITUDES/PHASES FOR THE ALL-YEARS COMPOSITE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % think I need to split amp/phase into two figs and then glue in ppt
%
% figure; hold all; whitefig; figpos([0.75 0.8])
%
% %-------------------------------------------------------
% vert_gap = 0.04;        horz_gap = 0.015;
% lower_marg = 0.075;     upper_marg = 0.04;
% left_marg = 0.065;      right_marg = 0.05;
%
% rows = 4; cols = 3;
%
% subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);
%
% %--------------------------------------------------------
%
% fs = 17;
%
% %--------------------------------------------------------
%
% days(end) = 366;
%
% ap = 'p'; % amplitude or phase?
%
% switch lower(ap)
%     case {'amp','a'}
%         axr = 1:2;
%         ampphase = 'amp';
%     case {'ph','phase','p'}
%         axr = 3:4;
%         ampphase = 'phase';
% end
%
% tidenames = {'Diurnal','Semidiurnal','Terdiurnal','Quardiurnal'};
%
% % [COLS,ROWS] = meshgrid(1:cols,1:rows);
% % ROWS = ROWS(:);
% % COLS = linearise(COLS');
%
% %%%%
% %%%% FOR EACH TIDE...
% for td = 1:3
%
%     axrange = (1:cols:(rows*cols)) + (td-1);
%
%     % select tide:
%     tdh = 24/td;
%
%     TD.uA = CompYear.Tides.u.A(:,:,td);
%     TD.uB = CompYear.Tides.u.B(:,:,td);
%     TD.vA = CompYear.Tides.v.A(:,:,td);
%     TD.vB = CompYear.Tides.v.B(:,:,td);
%
%     %%%% Over-Interpolation?
%     odays = 1:0.25:366;
%     ozrange = 76:0.25:105;
%     for uv = 'uv'
%         for ab = 'AB'
%             F = griddedInterpolant({zrange,days},TD.([uv ab]),'linear','none');
%             TD.([uv ab]) = F({ozrange,odays});
% %             smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
% %             tbp = TD.([uv ab]); nanlocs = isnan(tbp); tbp(nanlocs) = nanmean(tbp(:));
% %             tbp = imgaussfilt(imgaussfilt(tbp,smoo1),smoo2);
% %             tbp(nanlocs) = NaN;
% %             TD.([uv ab]) = tbp;
%         end
%     end
%
%     %%%% first, compute tidal amps and phases:
%     umag = quadadd(TD.uA,TD.uB);
%     vmag = quadadd(TD.vA,TD.vB);
%
%     uph  = atan2(TD.uB,TD.uA);
%     vph  = atan2(TD.vB,TD.vA);
%
%     % convert to local time:
%     uph = wrapToPi(uph + ((HWD.Meta.TimeZone/24) * (2*pi)));
%     vph = wrapToPi(vph + ((HWD.Meta.TimeZone/24) * (2*pi)));
%
%     % % wrap
%     % uph(uph < 0) = (2*pi) + uph(uph < 0);
%     % vph(vph < 0) = (2*pi) + vph(vph < 0);
%
%     % convert to hours...
%     uph = (uph./(2*pi)) .* tdh;
%     vph = (vph./(2*pi)) .* tdh;
%
%
%     uph(uph < 0) = tdh + uph(uph < 0);
%     vph(vph < 0) = tdh + vph(vph < 0);
%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%% TBP and settings
%
%     TBP = {...
%         umag,...
%         vmag,...
%         uph,...
%         vph};
%
%     % different limits for different tides...
%     phaseclims = [0 tdh];
%     phaseclevs = 0:1:tdh;
%     phaseclines = 0:(tdh/4):(tdh);
%     phasecbarticks = 0:(tdh/4):(tdh);
%     phasecbarminorticks = 0:(tdh/8):(tdh);
%     switch td
%         case 1 % 24h tide
%             clims = {[0 20],[0 20],phaseclims,phaseclims};
%             clevs = {0:1:50,0:1:50,phaseclevs,phaseclevs};
%             clines = {0:5:50,0:5:50,phaseclines,phaseclines};
%             cbarticks = {0:5:60,0:5:60,phasecbarticks,phasecbarticks};
%             cbarminorticks = {0:2.5:60,0:2.5:60,phasecbarminorticks,phasecbarminorticks};
%         case 2 % 12h tide
%             clims = {[0 70],[0 70],phaseclims,phaseclims};
%             clevs = {0:5:100,0:5:100,phaseclevs,phaseclevs};
%             clines = {0:20:100,0:20:100,phaseclines,phaseclines};
%             cbarticks = {0:20:100,0:20:100,phasecbarticks,phasecbarticks};
%             cbarminorticks = {0:10:100,0:10:100,phasecbarminorticks,phasecbarminorticks};
%         case 3 %  8h tide
%             clims = {[0 20],[0 20],phaseclims,phaseclims};
%             clevs = {0:1:50,0:1:50,phaseclevs,phaseclevs};
%             clines = {0:5:50,0:5:50,phaseclines,phaseclines};
%             cbarticks = {0:5:60,0:5:60,phasecbarticks,phasecbarticks};
%             cbarminorticks = {0:2.5:60,0:2.5:60,phasecbarminorticks,phasecbarminorticks};
%         case 4 %  6h tide
%             clims = {[0 10],[0 10],phaseclims,phaseclims};
%             clevs = {0:0.5:20,0:0.5:20,phaseclevs,phaseclevs};
%             clines = {0:1:20,0:1:20,phaseclines,phaseclines};
%             cbarticks = {0:5:60,0:5:60,phasecbarticks,phasecbarticks};
%             cbarminorticks = {0:1:60,0:1:60,phasecbarminorticks,phasecbarminorticks};
%     end
%
%
%     clabs = {...
%         'ms^{-1}','ms^{-1}','Hours (LT)','Hours (LT)'};
%
%     tits = {...
%         'ZONAL',...
%         'MERIDIONAL',...
%         'ZONAL',...
%         'MERIDIONAL'};
%
%
%     % and colour maps
%     cmaps = {...
%         nph_saturate(cbrew('nph_modspectral',20),1.2),...
%         nph_saturate(cbrew('nph_modspectral',20),1.2),...
%         nph_saturate(cbrew('nph_CyclicRainbow',24),1),...
%         nph_saturate(cbrew('nph_CyclicRainbow',24),1)};
%     % cmaps = {...
%     %     nph_saturate(cbrew('nph_modspectral',20),1.2),...
%     %     nph_saturate(cbrew('nph_modspectral',20),1.2),...
%     %     [flipud(cbrew('Blues',9)) ; cbrew('Blues',9)],...
%     %     nph_saturate(cbrew('nph_CyclicRainbow',18),1)};
%
%
%     for ax = axr
%
%         subplot(rows,cols,axrange(ax))
%
%         axx = gca;
%
%         % choose wot 2 plot n wot not
%         switch ax
%             case {1,2}
% %                 [X,Y] = meshgrid(datenum(min(years),00,days),zrange);
%                 [X,Y] = meshgrid(datenum(min(years),00,odays),ozrange);
%                 smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
%                 tbp = TBP{ax}; nanlocs = isnan(tbp); tbp(nanlocs) = nanmean(tbp(:));
%                 tbp = imgaussfilt(imgaussfilt(tbp,smoo1),smoo2);
%                 tbp(nanlocs) = NaN;
%             case {3,4} % no smoothing for phase
% %                 [X,Y] = meshgrid(datenum(min(years),00,days),zrange);
%                 [X,Y] = meshgrid(datenum(min(years),00,odays),ozrange);
%                 tbp = TBP{ax};
%         end
%
%
%
%         % allll the contours...
%         switch ax
%             case {1,2}
%                 %%%% filled contours:
%                 hold on; [C1,h1] = contourf(X,Y,tbp,clevs{ax},'edgecolor','none');
%                 %%%% contour lines:
%                 hold on; [C,h] = contour(X,Y,tbp,clines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
%                 clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
%                 %%%% and bold zero line
%                 hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','w','linewi',4);
%                 clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000,'color','w','fontweight','bold');
%                 hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linewi',2);
%                 clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000);
%             case {3,4}
%
%                 % pcolor instead for phase?
%                 hold on; pcolor(X,Y,tbp); shat;
%
%                 %%%% shall we try some contours for phase too?
%                 %%%% just try a few...
%
%                 % first zero phase:
%                 tbp2 = tbp; tbp2(tbp > tdh/2) = -tdh + tbp2(tbp > tdh/2);
%                 nanmask = abs(diff(tbp2([1 1:end],:),[],1)) > 6 | abs(diff(tbp2(:,[1 1:end]),[],2)) > 6;
%                 tbp2(nanmask) = NaN;
%                 hold on; [C,h] = contour(X,Y,tbp2,[0 0],'edgecolor',rgbtrip(.25),'linewi',1.25);
%                 clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
%
%                 % before you do the others, exclude any wraparound bits:
%                 nanmask = abs(diff(tbp([1 1:end],:),[],1)) > 6 | abs(diff(tbp(:,[1 1:end]),[],2)) > 6;
%                 tbp(nanmask) = NaN;
%
%                 hold on; [C,h] = contour(X,Y,tbp,tdh .* [0.25 0.5 0.75],'edgecolor',rgbtrip(.25),'linewi',1.25);
%                 clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
%
%         end
%
%
%
%
%         %%%% also draw on some contour labels in grey:
%         switch td
%             case 1
%                 extraclines = [10 10];
%             case 2
%                 extraclines = [40 40];
%             case 3
%                 extraclines = [50 50];
%             otherwise
%                 extraclines = [100 100];
%         end
%         switch ax
%             case {1,2}
%                 % contour lines:
%                 hold on; [C,h] = contour(X,Y,tbp,extraclines,'edgecolor',rgbtrip(.25),'linewi',1.25);
%                 clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.8),'fontweight','bold');
%         end
%
%
%         % LIMITS
%         ylim([77 103])
%         axx.YTick = 70:5:110;
%         axx.YMinorTick = 'on';
%         axx.YAxis.MinorTickValues = 70:110;
%
%         xlim([datenum(min(years),01,01) datenum(min(years),13,01)])
%
%         % Just show the Januarys as major ticks
%         axx.XTick = datenum(min(years),1:13,01);
%         datetick('x','m','keepticks','keeplimits')
%
%         switch axrange(ax)
%             case {1,2,3,7,8,9}
%                 % remove end Jan and nudge along to middle of ticks
%                 axx.XTickLabel(end) = ' ';
%                 axx.XTickLabel = cellstr(axx.XTickLabel);
%                 for lab = 1:length(axx.XTick)
%                     axx.XTickLabel{lab} = ['      ' axx.XTickLabel{lab}];
%                 end
%             otherwise
%                 axx.XTickLabel = {};
%         end
%
%         % only label for LHS
%         switch td
%             case 1
%                 ylabel('Altitude (km)')
%             otherwise
%                 axx.YTickLabel = {};
%         end
%
%         % line around the axis
%         hold on; plot(axx.XLim([1 2 2]),axx.YLim([2 2 1]),'color','k','linewi',1.5)
%
%         drawnow;
%
%         %     % AXES LETTER
%         hold on; nph_text([0.005 0.85],['(' alphabet(axrange(ax)) ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
%         hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','color','k','textborder','w','horizontalalignment','left','verticalalignment','middle');
%
%
%         % I COLORI
%         cmap = cmaps{ax};
%         colormap(gca,cmap)
%
%         clim(clims{ax})
%
%         % COLORBAR
%         switch axrange(ax)
%             case {4,5,6,10,11,12}
%                 cbar = nph_colorbar;
%                 cbar.Location = 'southoutside';
%                 cbar.Ticks = cbarticks{ax};
%                 cbar.Position = cbar.Position .* [1 1 1 1]; drawnow;
%                 cbar.Position = [...
%                     cbar.Position(1)+0.2*cbar.Position(3) ...
%                     cbar.Position(2)-4.5*cbar.Position(4) ...
%                     0.6*cbar.Position(3) ...
%                     1.6*cbar.Position(4)];
%
%                 % minor ticks
%                 cbar.Ruler.MinorTick = 'on';
%                 cbar.Ruler.MinorTickValues = cbarminorticks{ax};
%
%                 % put Label at the end of colorbar:
%                 cbar.Label.String = ['\bf{' clabs{ax} '}'];
%                 switch ax
%                     case {1,2}
%                         cbar.Label.Position = cbar.Label.Position .* [2.3 -1.5 1];
%                     case {3,4}
%                         cbar.Label.Position(1) = tdh/2;
%                         cbar.Label.Position = cbar.Label.Position .* [2.4 -1.25 0];
%                 end
%                 %             cbar.Label.Rotation = 0;
%                 %             cbar.Label.VerticalAlignment = 'middle';
%                 %             cbar.Label.HorizontalAlignment = 'left';
%                 cbar.TickDirection = 'out';
%                 cbar.LineWidth = 1.5;
%         end
%
%         switch ax
%             case 1
% %                 title([num2str(tdh) 'h Tide'])
%                 title(tidenames{td})
%         end
%
%         setfont(fs)
%
%         set(gca,'linewi',1.5,'layer','top','tickdir','out')
%
%     end
%
% end
%
%
%
% return
%
%
% %% EXPORT??? ==============================================================
%
% % savename = ['~/Desktop/' upper(site) '_TidesAllYearsComp_' ampphase];
% savename = ['~/Desktop/' upper(site) '_' num2str(tdh) 'hTideComp_' ampphase];
%
% disp(['Saving as "' savename '"...'])
%
% nph_saveas(gcf,savename,'png')
%
% disp('Saved.')
%
%
% % end
%
%
%
%
% return
%
%
%
%
%
%
%































