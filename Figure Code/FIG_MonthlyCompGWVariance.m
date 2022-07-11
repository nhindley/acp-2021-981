
%%%% NOTE: Also contains pseudo-wavefield plots to show the difference
%%%% between zonal/meridional GW variance


% Plot GW variance (monthly composites) for the whole KEP mission so far.

% Then, underneath, plot composite years of GW variance and direction

% EDIT Feb 2020: So thanks to the foresight of past me, there are now
% monthly composite days in the HWD files. Fancy that!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SITE AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


site = 'kep';

years = 2016:2020;

% First, stick all the monthly composite mean winds together for the
% all-time plot:

Comp.u = [];
Comp.v = [];
Comp.gwvar = [];
Comp.gwvar_zonal = [];
Comp.gwvar_merid = [];
Comp.walt = [];
Comp.wtime = [];
Comp.time = [];

Comp.gwres_zonal = [];
Comp.gwres_merid = [];
Comp.gwres_walt = [];
Comp.gwres_wtime = [];

for AB = 'ABC'
    for uv = 'uv'
        Comp.Tides.OneDay.(uv).(AB) = [];
        Comp.Tides.(uv).(AB) = [];
    end
end
    
for yr = years
    
    yrstr = num2str(yr);
    
    load(['/Users/neil/data/MeteorRadar/' site '/matlab/hwd/' yrstr '_' site '_hwd_gaussbins.mat'])
    
%     disp(yrstr);
%     disp(length(HWD.Data.Alt));
%     continue

    Comp.time = cat(2,Comp.time,repmat(datenum(yr,1:12,15),30,1));
    
    Comp.u              = cat(2,Comp.u,sq(nanmean(HWD.Data.MonthlyComp.u,2)));
    Comp.v              = cat(2,Comp.v,sq(nanmean(HWD.Data.MonthlyComp.v,2)));
    Comp.gwvar          = cat(2,Comp.gwvar,HWD.Data.MonthlyComp.SubVolumeGWVariance);
    Comp.gwvar_zonal    = cat(2,Comp.gwvar_zonal,HWD.Data.MonthlyComp.SubVolumeGWVarianceZonal);
    Comp.gwvar_merid    = cat(2,Comp.gwvar_merid,HWD.Data.MonthlyComp.SubVolumeGWVarianceMerid);
    
    Comp.walt = cat(2,Comp.walt,sq(nanmean(HWD.Data.MonthlyComp.walt,2)));
    Comp.wtime = cat(2,Comp.wtime,sq(nanmean(HWD.Data.MonthlyComp.wtime,2)));
    
%     Comp.uu             = cat(2,Comp.uu,HWD.Data.u);
%     Comp.vv             = cat(2,Comp.vv,HWD.Data.v);
    
    for AB = 'ABC'
        for uv = 'uv'
            Comp.Tides.OneDay.(uv).(AB) = cat(2,Comp.Tides.OneDay.(uv).(AB),HWD.Data.Tides.OneDay.TidalComponents.(uv).(AB));
            Comp.Tides.(uv).(AB) = cat(2,Comp.Tides.(uv).(AB),HWD.Data.Tides.TidalComponents.(uv).(AB));
        end
    end
    
    alltidesu = HWD.Data.Tides.AllTides.u;
    alltidesv = HWD.Data.Tides.AllTides.v;
    
    Comp.gwres_zonal    = cat(2,Comp.gwres_zonal,HWD.Data.u-alltidesu);
    Comp.gwres_merid    = cat(2,Comp.gwres_merid,HWD.Data.v-alltidesv);
    
    Comp.gwres_walt = cat(2,Comp.gwres_walt,HWD.Data.walt);
    Comp.gwres_wtime = cat(2,Comp.gwres_wtime,HWD.Data.wtime);
    
end

% %%%% fill in the missing one day tidal components with the 4 day fits
% for AB = 'ABC'
%     for uv = 'uv'
%         nanlocs = isnan(Comp.Tides.OneDay.(uv).(AB));
%         Comp.Tides.OneDay.(uv).(AB)(nanlocs) = Comp.Tides.(uv).(AB)(nanlocs);
%         %%%% and the moving average:
%         nanlocs = isnan(Comp.Tides.OneDay.(uv).(AB));
%         smoosmoo = movmean(Comp.Tides.OneDay.(uv).(AB),30*24,2,'omitnan');
%         Comp.Tides.OneDay.(uv).(AB)(nanlocs) = smoosmoo(nanlocs);
%     end
% end

% return


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
        allMPD(d).GWResiduals           = MPD.Data.GWResiduals;
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
gwres   = cat(1,allMPD(:).GWResiduals);

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
gwres = gwres(inds);

hourofday = mod(tim,1) * 24;
dn = daynumber(tim);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% now make a composite year... :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('Assembling an all-years composite for the specified years...')

days = ceil(linspace(1,365,120));

daywindow = 30;

std_z = 1.275; % FWHM for height gaussian
% zrange = 76:1:105;
zrange = HWD.Data.Alt;

CompYear                = struct;
CompYear.gwvar          = nan(length(zrange),length(days));
CompYear.gwvar_zonal    = nan(length(zrange),length(days));
CompYear.gwvar_merid    = nan(length(zrange),length(days));
CompYear.walt           = nan(length(zrange),length(days));

for d = 1:length(days)
    
    dd = days(d);
    
    disp(num2str(dd))
    
    drange = (dd-daywindow/2):(dd+daywindow/2);
    
    % fix wraparounds:
    drange(drange < 0)      = 365 + drange(drange < 0);
    drange(drange > 365)    = drange(drange > 365) - 365;
    
    % select indeces of this day window
    dayinds = ismember(dn,drange);

    
    % select todays composite meteors:
    az_day      = az(dayinds);
    gwres_day   = gwres(dayinds);
    alt_day     = alt(dayinds);
    tim_day     = tim(dayinds);
    dn_day      = dn(dayinds);
    
    % To save a bit of time, pre-compute the gaussian weighting
    % functions for the height bins, they'll be the same each time:
    Gauss_z = struct;
    for z = 1:length(zrange)
        Gauss_z(z).vec = exp(- ((alt_day - zrange(z)).^2) ./ (2 * std_z^2));
    end
    
    
    % fit!
    for z = 1:length(zrange)
        
        % altitude weightings:
        w_z = Gauss_z(z).vec;
        
        % zonal/merid weightings
%         w_zonal = w_z .* abs(sind(az_day));
%         w_merid = w_z .* abs(cosd(az_day));
        w_zonal = w_z .* sind(az_day).^2;
        w_merid = w_z .* cosd(az_day).^2;
        
        % only select weights above a threshold:
        inds_z = w_z > 0.05;
        inds_zonal = w_zonal > 0.05;
        inds_merid = w_merid > 0.05;

        % Subscribe!
        CompYear.gwvar(z,d) = var(gwres_day(inds_z),w_z(inds_z),'omitnan');
        CompYear.gwvar_zonal(z,d) = var(gwres_day(inds_zonal),w_zonal(inds_zonal),'omitnan');
        CompYear.gwvar_merid(z,d) = var(gwres_day(inds_merid),w_merid(inds_merid),'omitnan');

        % Subscribe the irregular weighted mean altitude location for this
        % bunch of meteors. We can then interpolate onto a
        % regular altitude grid later.
        CompYear.walt(z,d) = wmean(alt_day(inds_z),w_z(inds_z));
        
    end

end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return


%% PLOTTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Zonal and merid GW Variance
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

T = tiledlayout(2,6,'tilespacing','compact');

spans = {...
    [1 5],...
    [1 5],...
    [1 1],...
    [1 1]};

% axs = {...
%     [1  6],...
%     [9 14],...
%     [7 8],...
%     [15 16]};

TBP = {...
    log10(Comp.gwvar_zonal),...
    log10(Comp.gwvar_merid),...
    log10(CompYear.gwvar_zonal),...
    log10(CompYear.gwvar_merid)};

clims = {...
    [2.2 3.08],...
    [2.2 3.08],...
    [2.2 3.08],...
    [2.2 3.08]};

% powerspace(2,3.2,50,2)
clevs = {...
    2:0.05:3,...
    2:0.05:3,...
    2:0.05:3,...
    2:0.05:3};

clabs = {...
    'Variance (m^{2}s^{-2})','Variance (m^{2}s^{-2})','Variance (m^{2}s^{-2})','Variance (m^{2}s^{-2})'};

cbarticks = {...
    log10([0 200 400 800 1200]),...
    log10([0 200 400 800 1200]),...
    log10([0 200 400 800 1200]),...
    log10([0 200 400 800 1200])};

cbarlabs = {};
for i = 1:length(cbarticks)
    cbarlabs{i} = strtrim(cellstr(num2str(10.^cbarticks{i}')));
end

tits = {...
    'ZONAL GW VARIANCE',...
    'MERIDIONAL GW VARIANCE',...
    'ZONAL GW VARIANCE',...
    'MERIDIONAL GW VARIANCE'};


% also contour lines:
clines = {...
    [0 200 400 600 800 1200 1600],...
    [0 200 400 600 800 1200 1600],...
    [0 200 400 600 800 1200 1600],...
    [0 200 400 600 800 1200 1600]};

darkclines = {...
    [400 400],...
    [400 400],...
    [400 400],...
    [400 400]};

letters = {'a','b','c','d'};

ynudge = 0;


for ax = 1:4
    
%     subplot(rows,cols,axs{ax})
    nexttile(spans{ax})
    
    axx = gca;
    
    % choose wot 2 plot n wot not
    switch ax
        case {1,2}
            X = Comp.wtime; Y = Comp.walt;
            smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
        case {3,4}
            [X,Y] = meshgrid(datenum(2019,00,days),nanmean(CompYear.walt,2));
            smoo1 = [0.25 0.75]; smoo2 = [0.25 0.75];
    end
    
    % allll the contours...
    tbp = TBP{ax}; nanlocs = isnan(tbp); tbp(nanlocs) = nanmean(tbp(:));
    tbp = imgaussfilt(imgaussfilt(tbp,smoo1),smoo2);
    tbp(nanlocs) = NaN;
    
%     %%%%%%%%%%%%%%%%%
%     %%%% possibly apply scale factor for height?
%     sf = exp( - (Y-90) ./ (2*7));
%     tbp = tbp .* sf;
%     %%%%%%%%%%%%%%%%%
    
    % filled contours:
    hold on; contourf(X,Y,tbp,clevs{ax},'edgecolor','none');
    % contour lines:
    hold on; [C,h] = contour(X,Y,10.^tbp,clines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',1000,'color',rgbtrip(.95),'fontweight','bold');
    % and bold zero line
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','w','linewi',4);
    clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000,'color','w','fontweight','bold');
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linewi',2);
    clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000);
    
    % and dark contour lines:
    % contour lines:
    hold on; [C,h] = contour(X,Y,10.^tbp,darkclines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',1000,'color',rgbtrip(.25),'fontweight','bold');
    
    
    
    
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
            switch site
                case 'kep'
                    xlim([datenum(2016,02,15) datenum(2020,11,15)])
            end
            % Just show the Januarys as major ticks
            axx.XTick = datenum(years,01,01);
            axx.XMinorTick = 'on';
            axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
            datetick('x',' ','keepticks','keeplimits')
%             axx.XTickLabel = {};
            % other months as ticks using text:
            xtix = datenum(years(1),1:(length(years)*12),15);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,75.5+ynudge,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
               end
            end
            if ax == 2
            % Years as lower bold numbers:
            xtix = datenum(years,07,01);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   yrstr = datestr(xt,'yyyy');
                   text(xt,72.5+ynudge-0,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
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
               if inrange(xt,axx.XLim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,76+ynudge,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
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
    
    
    % I COLORI
    cmap = nph_saturate(cbrew('nph_modspectral',18),1.2);
    colormap(gca,cmap)
    clim(clims{ax})
    
    
    % COLORBAR
    switch ax
        
       case {3,4}
            cbar = nph_colorbar;
            cbar.Ticks = cbarticks{ax};
            cbar.Position = cbar.Position .* [1 1 3 1];
            cbar.Position = cbar.Position +  [cbar.Position(3)/3 0*cbar.Position(4) 0 0];
            
%             cbar.Title.String = {['\bf{' clabs{ax} '}'],''};
%             cbar.Title.String = {['\bf{GW Variance}'],'(m^{2}s^{-2})'};
            
            cbar.TickDirection = 'out';
            
            cbar.TickLabels = cbarlabs{ax};
            cbar.Ruler.MinorTick = 'on';
            cbar.Ruler.MinorTickValues = log10(0:100:2000);
%             cbar.Label.Rotation = 0;
%             cbar.Label.VerticalAlignment = 'middle';
%             cbar.Label.HorizontalAlignment = 'left';
            cbar.Title.String = '   m^{2}s^{-2}';     
    end
    
    
    
    % AXES LETTER
    %letter = alphabet(ax);
    letter = letters{ax};
    switch ax
        case {1,2}
            hold on; nph_text([0.005 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
            hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        case 3 % no idea why these are different
            hold on; nph_text([0.075 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
        case 4 % no idea why these are different
            hold on; nph_text([0.025 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');

            %             hold on; nph_text([2*0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        
    end

    
    
    setfont(fs)
    
    set(gca,'linewi',1.5,'layer','top','tickdir','out')
    
    
    
end


return

%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_GWVarianceZonalMerid_thin'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return
























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% RESOLVED GW Variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% create some monthly resolved GW variance too...
gwres = quadadd(Comp.gwres_zonal,Comp.gwres_merid);
% nanlocs = isnan(gwres); smoo = movmean(gwres,6,2,'omitnan');
% gwres(nanlocs) = smoo(nanlocs);
% gwres(isnan(gwres)) = 0;
gwres_filt = gwres - movmean(gwres,6,2,'omitnan');
% gwres_filt = gwres - imgaussfilt(gwres,[0.1 6./2.355]);
% gwres_var = nan(length(HWD.Data.Alt),length(years)*12);
% for yr = 1:length(years)
%     for mn = 1:12
%         data = gwres_filt(:,inrange(nanmean(Comp.gwres_wtime,1),datenum(years(yr),[mn mn+1],1)));
%         varrr = std(data,[],2,'omitnan');
%         if ~isempty(varrr)
%             gwres_var(:,mn+((yr-1)*12)) = varrr;
%         end
%     end
% end

gwres_var = movstd(gwres_filt,30*24,[],2,'omitnan').^2;

monthly_timevec = datenum(min(years),1:(12*length(years)),01);
[XX,YY] = meshgrid(monthly_timevec,HWD.Data.Alt);

Comp.gwres_var = nph_bin2mat(Comp.gwres_wtime,Comp.gwres_walt,gwres_var,monthly_timevec,HWD.Data.Alt)';
CompYear.gwres_var = nph_bin2mat(daynumber(Comp.gwres_wtime),Comp.gwres_walt,gwres_var,days,HWD.Data.Alt)';



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

T = tiledlayout(2,6,'tilespacing','compact');

spans = {...
    [1 5],...
    [1 5],...
    [1 1],...
    [1 1]};

% axs = {...
%     [1  6],...
%     [9 14],...
%     [7 8],...
%     [15 16]};

TBP = {...
    log10(Comp.gwres_var),...
    log10(Comp.gwres_var),...
    log10(CompYear.gwres_var),...
    log10(CompYear.gwres_var)};

clims = {...
    [1.25 2.45],...
    [1.25 2.45],...
    [1.25 2.45],...
    [1.25 2.45]};

clevs = {...
    log10([10:2:50 50:5:100 100:25:750]),...
    log10([10:2:50 50:5:100 100:25:750]),...
    log10([10:2:50 50:5:100 100:25:750]),...
    log10([10:2:50 50:5:100 100:25:750])};

clabs = {...
    'Variance (m^{2}s^{-2})','Variance (m^{2}s^{-2})','Variance (m^{2}s^{-2})','Variance (m^{2}s^{-2})'};

cbarticks = {...
    log10([0 25 50 100 200 400]),...
    log10([0 25 50 100 200 400]),...
    log10([0 25 50 100 200 400]),...
    log10([0 25 50 100 200 400])};

cbarlabs = {};
for i = 1:length(cbarticks)
    cbarlabs{i} = strtrim(cellstr(num2str(10.^cbarticks{i}')));
end

tits = {...
    'RESOLVED GW VARIANCE',...
    'GW VARIANCE (NORMALISED)',...
    'RESOLVED GW VARIANCE',...
    'GW VARIANCE (NORMALISED)'};

% also contour lines:
clines = {...
    [0 50 100:100:500],...
    [0 50 100:100:500],...
    [0 50 100:100:500],...
    [0 50 100:100:500]};

darkclines = {...
    [50 100],...
    [50 100],...
    [50 100],...
    [50 100]};

letters = {'a','c','b','d'};


for ax = [1 3]
    
%     subplot(rows,cols,axs{ax})
    nexttile(spans{ax})
    
    axx = gca;
    
    % choose wot 2 plot n wot not
    switch ax
        case {1,2}
%             X = Comp.time; Y = Comp.walt;
            X = XX; Y = Comp.walt;
            smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
        case {3,4}
            [X,Y] = meshgrid(datenum(2019,00,days),nanmean(CompYear.walt,2));
            smoo1 = [0.25 0.75]; smoo2 = [0.25 0.75];
    end
    
    % allll the contours...
    tbp = TBP{ax};
    nanlocs = isnan(tbp);
    switch ax
        case {1,2}
            smoosmoo = movmean(tbp,3,2,'omitnan');
            tbp(nanlocs) = smoosmoo(nanlocs);
    end
    nanlocs = isnan(tbp);
    tbp(nanlocs) = nanmean(tbp(:));
    tbp = imgaussfilt(imgaussfilt(tbp,smoo1),smoo2);
    tbp(nanlocs) = NaN;
    
    %%%%%%%%%%%%%%%%%
    %%%% possibly apply scale factor for height?
    switch ax
        case {2,4}
            sf = exp( - (Y-90) ./ (2*7));
            tbp = log10((10.^tbp) .* sf);
    end
    %%%%%%%%%%%%%%%%%
    
    % filled contours:
    hold on; contourf(X,Y,tbp,clevs{ax},'edgecolor','none');
    % contour lines:
    hold on; [C,h] = contour(X,Y,10.^tbp,clines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',1000,'color',rgbtrip(.95),'fontweight','bold');
    % and bold zero line
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','w','linewi',4);
    clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000,'color','w','fontweight','bold');
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linewi',2);
    clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000);
    
    % and dark contour lines:
    % contour lines:
    hold on; [C,h] = contour(X,Y,10.^tbp,darkclines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',1000,'color',rgbtrip(.25),'fontweight','bold');
    
    
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
            switch site
                case 'kep'
                    xlim([datenum(2016,02,15) datenum(2020,11,15)])
            end
%             xlim([datenum(2016,02,15) datenum(2020,11,15)])
            % Just show the Januarys as major ticks
            axx.XTick = datenum(years,01,01);
            axx.XMinorTick = 'on';
            axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
            datetick('x',' ','keepticks','keeplimits')
%             axx.XTickLabel = {};
            % other months as ticks using text:
            xtix = datenum(years(1),1:(length(years)*12),15);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,75.5+ynudge,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
               end
            end
            if ax == 2
            % Years as lower bold numbers:
            xtix = datenum(years,07,01);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   yrstr = datestr(xt,'yyyy');
                   text(xt,72.5+ynudge-0,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
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
               if inrange(xt,axx.XLim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,76+ynudge,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
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
    
    
    % I COLORI
    cmap = nph_saturate(cbrew('nph_modspectral',18),1.2);
    colormap(gca,cmap)
    clim(clims{ax})
    
    
    % COLORBAR
    switch ax
        
       case {3,4}
            cbar = nph_colorbar;
            cbar.Ticks = cbarticks{ax};
            cbar.Position = cbar.Position .* [1 1 3 1];
            cbar.Position = cbar.Position +  [cbar.Position(3)/3 0*cbar.Position(4) 0 0];
            
%             cbar.Title.String = {['\bf{' clabs{ax} '}'],''};
%             cbar.Title.String = {['\bf{GW Variance}'],'(m^{2}s^{-2})'};
            
            cbar.TickDirection = 'out';
            
            cbar.TickLabels = cbarlabs{ax};
            cbar.Ruler.MinorTick = 'on';
            cbar.Ruler.MinorTickValues = unique(log10([25:25:500]));
%             cbar.Label.Rotation = 0;
%             cbar.Label.VerticalAlignment = 'middle';
%             cbar.Label.HorizontalAlignment = 'left';
            cbar.Title.String = '   m^{2}s^{-2}';     
    end
    
    
%     
%     % AXES LETTER
%     %letter = alphabet(ax);
%     letter = letters{ax};
%     switch ax
%         case {1,2}
%             hold on; nph_text([0.005 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
%             hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
%         case 3 % no idea why these are different
%             hold on; nph_text([0.075 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
%         case 4 % no idea why these are different
%             hold on; nph_text([0.025 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
% 
%             %             hold on; nph_text([2*0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
%         
%     end
    % AXES LETTER
    %letter = alphabet(ax);
    letter = letters{ax};
    switch ax
        case {1,2}
            hold on; nph_text([0.005 0.85],['(' letter ')'],'fontsize',1.5*fs,  'color','w','textborder','k','horizontalalignment','left','verticalalignment','middle');
            hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'color','w','textborder','k','fontweight','bold','horizontalalignment','left','verticalalignment','middle');
        case 3 % no idea why these are different
            hold on; nph_text([0.075 0.85],['(' letter ')'],'fontsize',1.5*fs,  'color','w','textborder','k','horizontalalignment','left','verticalalignment','middle');
        case 4 % no idea why these are different
            hold on; nph_text([0.025 0.85],['(' letter ')'],'fontsize',1.5*fs,  'color','w','textborder','k','horizontalalignment','left','verticalalignment','middle');

            %             hold on; nph_text([2*0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        
    end
    
    
    setfont(fs)
    
    set(gca,'linewi',1.5,'layer','top','tickdir','out')
    
    
    
end


return

%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_ResolvedGWVariance'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Directional Tendencies of GW Variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% use a new approach to derive the directional tendencies of GW
%%%% variance. By writing the equation below, we can easily relate the
%%%% zonal and meridional GW variances are fractions of each other, which
%%%% makes it much easier to describe in text (e.g. there is a 10%
%%%% tenedency to the zonal etc)
% zonal = merid + (E / 100 * merid);
% zonal = (1 + E/100) * merid;
% E = 100 .* (Comp.gwvar_zonal./Comp.gwvar_merid - 1);

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

T = tiledlayout(2,6,'tilespacing','compact');

spans = {...
    [1 5],...
    [1 5],...
    [1 1],...
    [1 1]};

% axs = {...
%     [1  6],...
%     [9 14],...
%     [7 8],...
%     [15 16]};

TBP = {...
    100 .* (Comp.gwvar_zonal./Comp.gwvar_merid - 1),...
    100 .* (Comp.gwvar_zonal./Comp.gwvar_merid - 1),...
    100 .* (CompYear.gwvar_zonal./CompYear.gwvar_merid - 1),...
    100 .* (CompYear.gwvar_zonal./CompYear.gwvar_merid - 1)};

clims = {...
    [-25 25],...
    [-25 25],...
    [-25 25],...
    [-25 25]};

cbarlims = {...
    [-25 25],...
    [-25 25],...
    [-25 25],...
    [-25 25]};

clevs = {...
    -30:30,...
    -30:30,...
    -30:30,...
    -30:30};

clabs = {...
    'Fraction','Fraction','Fraction','Fraction'};

cbarticks = {...
    -30:10:30,...
    -30:10:30,...
    -30:10:30,...
    -30:10:30};

cbarminorticks = {...
    -30:5:30,...
    -30:5:30,...
    -30:5:30,...
    -30:5:30};

tits = {...
    'DIRECTIONAL TENDENCY',...
    'DIRECTIONAL TENDENCY',...
    'DIRECTIONAL TENDENCY',...
    'DIRECTIONAL TENDENCY'};

% also contour lines:
clines = {...
    [-30:10:30],...
    [-30:10:30],...
    [-30:10:30],...
    [-30:10:30]};

darkclines = {...
    [0 0],...
    [0 0],...
    [0 0],...
    [0 0]};

letters = {'e','e','f','f'};

ynudge = 0;
% ylims = [77 103];
ylims = [77+ynudge 103];

for ax = 1:4
    
%     subplot(rows,cols,axs{ax})
    nexttile(spans{ax})
    
    axx = gca;
    
    % choose wot 2 plot n wot not
    switch ax
        case {1,2}
            X = Comp.wtime; Y = Comp.walt;
            smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
        case {3,4}
            [X,Y] = meshgrid(datenum(2019,00,days),nanmean(CompYear.walt,2));
            smoo1 = [0.25 0.75]; smoo2 = [0.25 0.75];
    end
    
    % allll the contours...
    tbp = TBP{ax}; nanlocs = isnan(tbp); tbp(nanlocs) = nanmean(tbp(:));
    tbp = imgaussfilt(imgaussfilt(tbp,smoo1),smoo2);
    tbp(nanlocs) = NaN;
    
    % filled contours:
    hold on; contourf(X,Y,tbp,clevs{ax},'edgecolor','none');
    % contour lines:
    hold on; [C,h] = contour(X,Y,tbp,clines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',1000,'color',rgbtrip(.95),'fontweight','bold');
%     % and bold zero line
%     hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','w','linewi',4);
%     clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000,'color','w','fontweight','bold');
%     hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linewi',2);
%     clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000);
%     
    % and dark contour lines:
    % contour lines:
    hold on; [C,h] = contour(X,Y,tbp,darkclines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',1000,'color',rgbtrip(.25),'fontweight','bold');
    
    
    
    
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
        switch site
            case 'kep'
                xlim([datenum(2016,02,15) datenum(2020,11,15)])
        end            % Just show the Januarys as major ticks
        axx.XTick = datenum(years,01,01);
        axx.XMinorTick = 'on';
        axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
        datetick('x',' ','keepticks','keeplimits')
        %             axx.XTickLabel = {};
        % other months as ticks using text:
        xtix = datenum(years(1),1:(length(years)*12),15);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,75.5+ynudge,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
               end
            end
            if ax == 2
            % Years as lower bold numbers:
            xtix = datenum(years,07,01);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   yrstr = datestr(xt,'yyyy');
                   text(xt,72.5+ynudge-0,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
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
               if inrange(xt,axx.XLim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,76+ynudge,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
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
    
    
    % I COLORI
    cmap = nph_saturate(cbrew('nph_modspectral',18),1.2);
    colormap(gca,cmap)
    clim(clims{ax})
    
    
    % COLORBAR
    switch ax
        
       case {3,4}
            cbar = nph_colorbar;
            cbar.Ticks = cbarticks{ax};
            cbar.Position = cbar.Position .* [1 1 3 1];
            cbar.Position = cbar.Position +  [cbar.Position(3)/3 0*cbar.Position(4) 0 0];
            
%             cbar.Title.String = {['\bf{' clabs{ax} '}'],''};
%             cbar.Title.String = {['\bf{GW Variance}'],'(m^{2}s^{-2})'};
%             cbar.Title.String = '     \bf{^{Fraction}}';  
            
            cbar.TickDirection = 'out';
            cbar.Ruler.MinorTick = 'on';
            cbar.Ruler.MinorTickValues = cbarminorticks{ax};
            
%             cbar.Limits = cbarlims{ax};
            
%             %%%% convert cbarticks:
%             cbarlabs = {};
%             for i = 1:length(cbarticks)
%                 cbarlabs{i} = strtrim(cellstr(num2str(10.^cbarticks{i}')));
%             end
            
            cbar.Label.String = '\bf{%}';
            cbar.Label.Rotation = 0;
            cbar.Label.VerticalAlignment = 'middle';
            cbar.Label.HorizontalAlignment = 'left';
               
    end
    
    
    
    % AXES LETTER
    %letter = alphabet(ax);
    letter = letters{ax};
    switch ax
        case {1,2}
            hold on; nph_text([0.005 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
            hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        case 3 % no idea why these are different
            hold on; nph_text([0.075 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
        case 4 % no idea why these are different
            hold on; nph_text([0.025 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');

            %             hold on; nph_text([2*0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        
    end

    
    
    setfont(fs)
    
    set(gca,'linewi',1.5,'layer','top','tickdir','out')
    
    
    
end


return

%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_GWVarDirectionalTendencies'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')


return





%% OLD - Fraction of GW variance, zonal/merid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zonal = merid + (E / 100 * merid);
% zonal = (1 + E/100) * merid;
% E = 100 .* (Comp.gwvar_zonal./Comp.gwvar_merid - 1);

figure; hold all; whitefig; figpos([1 0.475])

%-------------------------------------------------------
vert_gap = 0.08;        horz_gap = 0.03;
lower_marg = 0.075;     upper_marg = 0.025;
left_marg = 0.065;      right_marg = 0.075;

rows = 2; cols = 8;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 18;

%--------------------------------------------------------

T = tiledlayout(2,6,'tilespacing','compact');

spans = {...
    [1 5],...
    [1 5],...
    [1 1],...
    [1 1]};

% axs = {...
%     [1  6],...
%     [9 14],...
%     [7 8],...
%     [15 16]};

TBP = {...
    log10(Comp.gwvar_zonal ./ Comp.gwvar_merid),...
    log10(Comp.gwvar_zonal ./ Comp.gwvar_merid),...
    log10(CompYear.gwvar_zonal ./ CompYear.gwvar_merid),...
    log10(CompYear.gwvar_zonal ./ CompYear.gwvar_merid)};

clims = {...
    [-0.09 0.09],...
    [-0.09 0.09],...
    [-0.09 0.09],...
    [-0.09 0.09]};

cbarlims = {...
    [-0.0975 0.0975],...
    [-0.0975 0.0975],...
    [-0.0975 0.0975],...
    [-0.0975 0.0975]};

% powerspace(2,3.2,50,2)
clevs = {...
    -0.5:0.001:0.5,...
    -0.5:0.001:0.5,...
    -0.5:0.001:0.5,...
    -0.5:0.001:0.5};

clabs = {...
    'Fraction','Fraction','Fraction','Fraction'};

cbarticks = {...
    log10([0.8 0.9 1 1.11 1.25]),...
    log10([0.8 0.9 1 1.11 1.25]),...
    log10([0.8 0.9 1 1.11 1.25]),...
    log10([0.8 0.9 1 1.11 1.25])};

cbarminorticks = {...
    -0.1:0.005:0.1,...
    -0.1:0.005:0.1,...
    -0.1:0.005:0.1,...
    -0.1:0.005:0.1};

tits = {...
    'ZONAL VARIANCE / MERIDIONAL VARIANCE',...
    'ZONAL VARIANCE / MERIDIONAL VARIANCE',...
    'ZONAL VARIANCE / MERIDIONAL VARIANCE',...
    'ZONAL VARIANCE / MERIDIONAL VARIANCE'};


% also contour lines:
clines = {...
    [0.8 0.85 0.9 1 1.1 1.15 1.2],...
    [0.8 0.85 0.9 1 1.1 1.15 1.2],...
    [0.8 0.85 0.9 1 1.1 1.15 1.2],...
    [0.8 0.85 0.9 1 1.1 1.15 1.2]};

darkclines = {...
    [1 1],...
    [1 1],...
    [1 1],...
    [1 1]};

letters = {'a','a','b','b'};

ynudge = 0;
% ylims = [77 103];
ylims = [77+ynudge 103];

for ax = 1:4
    
%     subplot(rows,cols,axs{ax})
    nexttile(spans{ax})
    
    axx = gca;
    
    % choose wot 2 plot n wot not
    switch ax
        case {1,2}
            X = Comp.wtime; Y = Comp.walt;
            smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
        case {3,4}
            [X,Y] = meshgrid(datenum(2019,00,days),nanmean(CompYear.walt,2));
            smoo1 = [0.25 0.75]; smoo2 = [0.25 0.75];
    end
    
    % allll the contours...
    tbp = TBP{ax}; nanlocs = isnan(tbp); tbp(nanlocs) = nanmean(tbp(:));
    tbp = imgaussfilt(imgaussfilt(tbp,smoo1),smoo2);
    tbp(nanlocs) = NaN;
    
    % filled contours:
    hold on; contourf(X,Y,tbp,clevs{ax},'edgecolor','none');
    % contour lines:
    hold on; [C,h] = contour(X,Y,10.^tbp,clines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',1000,'color',rgbtrip(.95),'fontweight','bold');
%     % and bold zero line
%     hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','w','linewi',4);
%     clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000,'color','w','fontweight','bold');
%     hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linewi',2);
%     clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000);
    
    % and dark contour lines:
    % contour lines:
    hold on; [C,h] = contour(X,Y,10.^tbp,darkclines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',1000,'color',rgbtrip(.25),'fontweight','bold');
    
    
    
    
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
            axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
            datetick('x',' ','keepticks','keeplimits')
%             axx.XTickLabel = {};
            % other months as ticks using text:
            xtix = datenum(years(1),1:(length(years)*12),15);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,75.5+ynudge,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
               end
            end
            if ax == 2
            % Years as lower bold numbers:
            xtix = datenum(years,07,01);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   yrstr = datestr(xt,'yyyy');
                   text(xt,72.5+ynudge-0,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
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
               if inrange(xt,axx.XLim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,76+ynudge,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
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
    
    
    % I COLORI
    cmap = nph_saturate(cbrew('nph_modspectral',18),1.2);
    colormap(gca,cmap)
    clim(clims{ax})
    
    
    % COLORBAR
    switch ax
        
       case {3,4}
            cbar = nph_colorbar;
            cbar.Ticks = cbarticks{ax};
            cbar.Position = cbar.Position .* [1 1 3 1];
            cbar.Position = cbar.Position +  [cbar.Position(3)/3 0*cbar.Position(4) 0 0];
            
%             cbar.Title.String = {['\bf{' clabs{ax} '}'],''};
%             cbar.Title.String = {['\bf{GW Variance}'],'(m^{2}s^{-2})'};
            cbar.Title.String = '     \bf{^{Fraction}}';  
            
            cbar.TickDirection = 'out';
            cbar.Ruler.MinorTick = 'on';
            cbar.Ruler.MinorTickValues = cbarminorticks{ax};
            
            cbar.Limits = cbarlims{ax};
            
            %%%% convert cbarticks:
            cbarlabs = {};
            for i = 1:length(cbarticks)
                cbarlabs{i} = strtrim(cellstr(num2str(10.^cbarticks{i}')));
            end
            
            cbar.TickLabels = cbarlabs{ax};
%             cbar.Label.Rotation = 0;
%             cbar.Label.VerticalAlignment = 'middle';
%             cbar.Label.HorizontalAlignment = 'left';
               
    end
    
    
    
    % AXES LETTER
    %letter = alphabet(ax);
    letter = letters{ax};
    switch ax
        case {1,2}
            hold on; nph_text([0.005 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
            hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        case 3 % no idea why these are different
            hold on; nph_text([0.075 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
        case 4 % no idea why these are different
            hold on; nph_text([0.025 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');

            %             hold on; nph_text([2*0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        
    end

    
    
    setfont(fs)
    
    set(gca,'linewi',1.5,'layer','top','tickdir','out')
    
    
    
end


return

%% OLD - EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_GWVarianceFraction_thin'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return


%% GW Variance, standard and scaled for altitude
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

T = tiledlayout(2,6,'tilespacing','compact');

spans = {...
    [1 5],...
    [1 5],...
    [1 1],...
    [1 1]};

% axs = {...
%     [1  6],...
%     [9 14],...
%     [7 8],...
%     [15 16]};

TBP = {...
    log10(Comp.gwvar),...
    log10(Comp.gwvar),...
    log10(CompYear.gwvar),...
    log10(CompYear.gwvar)};

clims = {...
    [2.2 3.08],...
    [2.2 3.08],...
    [2.2 3.08],...
    [2.2 3.08]};

% powerspace(2,3.2,50,2)
clevs = {...
    2:0.05:3,...
    2:0.05:3,...
    2:0.05:3,...
    2:0.05:3};

clabs = {...
    'Variance (m^{2}s^{-2})','Variance (m^{2}s^{-2})','Variance (m^{2}s^{-2})','Variance (m^{2}s^{-2})'};

cbarticks = {...
    log10([0 200 400 800 1200]),...
    log10([0 200 400 800 1200]),...
    log10([0 200 400 800 1200]),...
    log10([0 200 400 800 1200])};

cbarlabs = {};
for i = 1:length(cbarticks)
    cbarlabs{i} = strtrim(cellstr(num2str(10.^cbarticks{i}')));
end

tits = {...
    'SUB-VOLUME GW VARIANCE',...
    'SUB-VOLUME GW VARIANCE (NORMALISED)',...
    'SUB-VOLUME GW VARIANCE',...
    'SUB-VOLUME GW VARIANCE (NORMALISED)'};


% also contour lines:
clines = {...
    [0 200 400 600 800 1200 1600],...
    [0 200 400 600 800 1200 1600],...
    [0 200 400 600 800 1200 1600],...
    [0 200 400 600 800 1200 1600]};

darkclines = {...
    [400 400],...
    [400 400],...
    [400 400],...
    [400 400]};

letters = {'c','e','d','f'};


for ax = 1:4
    
%     subplot(rows,cols,axs{ax})
    nexttile(spans{ax})
    
    axx = gca;
    
    % choose wot 2 plot n wot not
    switch ax
        case {1,2}
            X = Comp.wtime; Y = Comp.walt;
            smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
        case {3,4}
            [X,Y] = meshgrid(datenum(2019,00,days),nanmean(CompYear.walt,2));
            smoo1 = [0.25 0.75]; smoo2 = [0.25 0.75];
    end
    
    % allll the contours...
    tbp = TBP{ax}; nanlocs = isnan(tbp); tbp(nanlocs) = nanmean(tbp(:));
    tbp = imgaussfilt(imgaussfilt(tbp,smoo1),smoo2);
    tbp(nanlocs) = NaN;
    
    %%%%%%%%%%%%%%%%%
    %%%% possibly apply scale factor for height?
    switch ax
        case {2,4}
            sf = exp( - (Y-90) ./ (2*7));
            tbp = log10((10.^tbp) .* sf);
    end
    %%%%%%%%%%%%%%%%%
    
    % filled contours:
    hold on; contourf(X,Y,tbp,clevs{ax},'edgecolor','none');
    % contour lines:
    hold on; [C,h] = contour(X,Y,10.^tbp,clines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',1000,'color',rgbtrip(.95),'fontweight','bold');
    % and bold zero line
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','w','linewi',4);
    clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000,'color','w','fontweight','bold');
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linewi',2);
    clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000);
    
    % and dark contour lines:
    % contour lines:
    hold on; [C,h] = contour(X,Y,10.^tbp,darkclines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',1000,'color',rgbtrip(.25),'fontweight','bold');
    
    
    
    
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
            axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
            datetick('x',' ','keepticks','keeplimits')
%             axx.XTickLabel = {};
            % other months as ticks using text:
            xtix = datenum(years(1),1:(length(years)*12),15);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,75.5+ynudge,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
               end
            end
            if ax == 2
            % Years as lower bold numbers:
            xtix = datenum(years,07,01);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   yrstr = datestr(xt,'yyyy');
                   text(xt,72.5+ynudge-0,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
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
               if inrange(xt,axx.XLim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,76+ynudge,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
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
    
    
    % I COLORI
    cmap = nph_saturate(cbrew('nph_modspectral',18),1.2);
    colormap(gca,cmap)
    clim(clims{ax})
    
    
    % COLORBAR
    switch ax
        
       case {3,4}
            cbar = nph_colorbar;
            cbar.Ticks = cbarticks{ax};
            cbar.Position = cbar.Position .* [1 1 3 1];
            cbar.Position = cbar.Position +  [cbar.Position(3)/3 0*cbar.Position(4) 0 0];
            
%             cbar.Title.String = {['\bf{' clabs{ax} '}'],''};
%             cbar.Title.String = {['\bf{GW Variance}'],'(m^{2}s^{-2})'};
            
            cbar.TickDirection = 'out';
            
            cbar.TickLabels = cbarlabs{ax};
            cbar.Ruler.MinorTick = 'on';
            cbar.Ruler.MinorTickValues = log10(0:100:2000);
%             cbar.Label.Rotation = 0;
%             cbar.Label.VerticalAlignment = 'middle';
%             cbar.Label.HorizontalAlignment = 'left';
            cbar.Title.String = '   m^{2}s^{-2}';     
    end
    
    
    
    % AXES LETTER
    %letter = alphabet(ax);
    letter = letters{ax};
    switch ax
        case {1,2}
            hold on; nph_text([0.005 0.85],['(' letter ')'],'fontsize',1.5*fs,  'color','w','textborder','k','horizontalalignment','left','verticalalignment','middle');
            hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'color','w','textborder','k','fontweight','bold','horizontalalignment','left','verticalalignment','middle');
        case 3 % no idea why these are different
            hold on; nph_text([0.075 0.85],['(' letter ')'],'fontsize',1.5*fs,  'color','w','textborder','k','horizontalalignment','left','verticalalignment','middle');
        case 4 % no idea why these are different
            hold on; nph_text([0.025 0.85],['(' letter ')'],'fontsize',1.5*fs,  'color','w','textborder','k','horizontalalignment','left','verticalalignment','middle');

            %             hold on; nph_text([2*0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        
    end

    
    
    setfont(fs)
    
    set(gca,'linewi',1.5,'layer','top','tickdir','out')
    
    
    
end


return

%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_GWVarianceAndScaled_thin'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%
%% PSUEDO-WAVEFIELD PLOTS 

% red blue wavefield plots to show the differen between zonal/merid
% variance:

[X,Y] = meshgrid(-500:500,-500:500);

kx = twopi./150;
ky = twopi./20000;

Z = 0.9*sin(kx.*X + ky.*Y);

clev = -1:0.05:1;

figure; whitefig; figpos([0.5 0.5])

hold on; pcolor(X,Y,Z); shat;
hold on; contourf(X,Y,Z,clev,'edgecolor','none');


clim([-1 1]);
colormap(gca,nph_saturate(cbrew('RdBu',2*length(clev)),0.8))

axis square;

set(gca,'color','none','xcolor','none','ycolor','none')

drawnow;


savename = ['~/Desktop/pseudo_wavefield'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')






























