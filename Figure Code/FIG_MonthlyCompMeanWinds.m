

%%%% EDIT: trying to standardise figure style for the KEP paper...

% Plot mean winds (monthly composites) for the whole KEP mission so far.

% Then, underneath, plot composite years of zonal and meridional winds
% from all years.

% EDIT Feb 2020: So thanks to the foresight of past me, there are now
% monthly composite days in the HWD files. Fancy that!

% ALSO PLOTS MLS TEMPERATURE IF SPECIFIED :D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SITE AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

site = 'kep';

years = 2016:2020;

% First, get some composite of HWD stuff anyway...
Comp.u = [];
Comp.v = [];
Comp.uu = [];
Comp.vv = [];
Comp.day = [];
Comp.monthlyfwhm = [];
Comp.monthlycenter = [];
Comp.month = [];
Comp.gwvar = [];
Comp.walt = [];
Comp.wtime = [];

% Comp.alt = [];
% Comp.time = [];
% Comp.waltwalt = [];
% Comp.wtimewtime = [];

for yr = years
    yrstr = num2str(yr);
    load(['/Users/neil/data/MeteorRadar/' site '/matlab/hwd/' yrstr '_' site '_hwd_gaussbins.mat'])
    
    HWD.Data.MonthlyComp.Alt = repmat(HWD.Data.Alt',24,12);

    Comp.u      = cat(2,Comp.u,sq(nanmean(HWD.Data.MonthlyComp.u,2)));
    Comp.v      = cat(2,Comp.v,sq(nanmean(HWD.Data.MonthlyComp.v,2)));
%     Comp.u      = cat(2,Comp.u,sq(HWD.Data.MonthlyComp.Tides.TidalComponents.u.C(:,:,4)));
%     Comp.v      = cat(2,Comp.v,sq(HWD.Data.MonthlyComp.Tides.TidalComponents.v.C(:,:,4)));
%     
%     Comp.alt   = cat(2,Comp.alt,sq(nanmean(HWD.Data.MonthlyComp.alt,2)));
%     Comp.time  = cat(2,Comp.time,sq(nanmean(HWD.Data.MonthlyComp.time,2)));
    
    Comp.walt   = cat(2,Comp.walt,sq(nanmean(HWD.Data.MonthlyComp.walt,2)));
    Comp.wtime  = cat(2,Comp.wtime,sq(nanmean(HWD.Data.MonthlyComp.wtime,2)));
    Comp.uu  = cat(2,Comp.uu,HWD.Data.u);
    Comp.vv  = cat(2,Comp.vv,HWD.Data.v);
    
%     Comp.alt   = cat(2,Comp.alt,repmat(HWD.Data.Alt',1,length(HWD.Data.Time)));
%     Comp.time  = cat(2,Comp.time,repmat(HWD.Data.Time,length(HWD.Data.Alt),1));
%     Comp.waltwalt  = cat(2,Comp.waltwalt,HWD.Data.walt);
%     Comp.wtimewtime  = cat(2,Comp.wtimewtime,HWD.Data.wtime);
    
    %     Comp.fwhm   = cat(2,Comp.fwhm,HWD.Data.VertMetDist.FWHM);
%     Comp.center = cat(2,Comp.center,HWD.Data.VertMetDist.Center);
    Comp.day    = cat(2,Comp.day,HWD.Data.Day);
    Comp.monthlyfwhm   = cat(2,Comp.monthlyfwhm,HWD.Data.MonthlyComp.VertMetDist.FWHM);
    Comp.monthlycenter = cat(2,Comp.monthlycenter,HWD.Data.MonthlyComp.VertMetDist.Center);
    Comp.month    = cat(2,Comp.month,datenum(yr,1:12,01));
end

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


%% %% LOAD MLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



mlsyears = 2016:2020;
mlsT = [];
mlsTdiff = [];
mlslat = [];
mlstime = [];

Comp.T          = nan(55,12*length(years));
Comp.Tdiff      = nan(55,12*length(years));
Comp.Tmonth     = nan(1,12*length(years));

for yr = mlsyears
    
    yrstr = num2str(yr);
    
    disp(['Loading MLS for ' yrstr '...'])
    
    load(['/Users/neil/data/MLS/' yrstr '_MLS_T.mat']);
    
    %%%% select region:
%     inds = ...
%         inrange(MLS.Data.Longitude,[-45 -30]) & ...
%         inrange(MLS.Data.Latitude ,[-60 -50]);
    inds = ...
        inrange(MLS.Data.Longitude,[-90 -40]) & ...
        inrange(MLS.Data.Latitude ,[-75 -55]);
    
    TT = MLS.Data.T(:,inds);
    latlat = MLS.Data.Latitude(inds);
    timtim = MLS.Data.Time(inds);
    mlsz = MLS.Data.Altitude;
    
    clear MLS
    
%     % sort out temperature gradient with latitude:
%     TTdiff = diff(TT,[],2); TTdiff = cat(2,TTdiff,TTdiff(:,end));
%     latdiff = diff(latlat); latdiff = cat(1,latdiff,latdiff(end));
%     TTdiff(latdiff < 0) = -TTdiff(latdiff < 0);
    
    %%%% stick it together:
    mlsT = cat(2,mlsT,TT);
%     mlsTdiff = cat(2,mlsTdiff,TTdiff);
    mlslat = cat(1,mlslat,latlat);
    mlstime = cat(1,mlstime,timtim);
    
    %%%% get monthly mean...
    for mn = 1:12
        minds = inrange(mlstime,datenum(yr,[mn mn+1],01));
        subind = (12*(yr  - mlsyears(1))) + mn;
        Comp.T(:,subind) = nanmean(mlsT(:,minds),2);
        Comp.Tmonth(subind) = datenum(yr,mn,15);
        
        % latitudinal gradient?
        tt = mlsT(:,minds);
        ll = mlslat(minds);
        for z = 1:55
            P = polyfit(ll,tt(z,:),1);
            Comp.Tdiff(z,subind) = P(1);
        end
        
    end
    
end

%%

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recalculate a rolling 30-day composite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Assembling an all-years composite for the specified years...')

days = ceil(linspace(1,365,60));

daywindow = 30;

std_z = 1.275; % FWHM for height gaussian
zrange = 76:1:105;
% zrange = 65:1:106;

CompYear = struct;
CompYear.u = nan(length(zrange),length(days));
CompYear.v = nan(length(zrange),length(days));
CompYear.walt = nan(length(zrange),length(days));

dn_mls          = daynumber(mlstime);
dn_mls_month    = daynumber(Comp.Tmonth);
CompYear.T      = nan(55,length(days));
CompYear.Tdiff  = nan(55,length(days));

% return

for d = 1:length(days)
    
    dd = days(d);
    
    disp(num2str(dd))
    
    drange = (dd-daywindow/2):(dd+daywindow/2);
    
    % deal with wraparound for years
    drange_neg = drange(drange <= 0);
    drange = drange(drange > 0);
    drange_pos = drange(drange >= 366);
    drange = drange(drange < 366);
    
    % if there a low wraparoud:
    if ~isempty(drange_neg)
        drange_neg = 365+drange_neg;
        timerange = inrange(dn,min(drange_neg),366) | inrange(dn,minmax(drange));
    else
        timerange = inrange(dn,minmax(drange));
    end
    
    % if there a low wraparoud:
    if ~isempty(drange_pos)
        drange_pos = drange_pos-365;
        timerange = inrange(dn,minmax(drange_pos)) | inrange(dn,minmax(drange));
    else
        timerange = inrange(dn,minmax(drange));
    end
        
    
    timerange_mls = ismember(dn_mls,drange);
    CompYear.T(:,d) = nanmean(mlsT(:,timerange_mls),2);
    
%     timerange_mls_month = ismember(dn_mls_month,drange);
%     CompYear.Tdiff(:,d) = nanmean(Comp.Tdiff(:,timerange_mls_month),2);
    % or try and redo the diff?
    tt = mlsT(:,timerange_mls);
    ll = mlslat(timerange_mls);
    for z = 1:55
        P = polyfit(ll,tt(z,:),1);
        CompYear.Tdiff(z,d) = P(1);
    end
    

    % select todays composite meteors:
    az_day = az(timerange);
    vhorz_day = vhorz(timerange);
    alt_day = alt(timerange);
    tim_day = tim(timerange);
    dn_day = dn(timerange);
    
%     %%%% as a sanity check for a reviewer concerned about the direction of
%     %%%% winds above 100km, remove any meteors below 100km and see what
%     %%%% happens:
%     vhorz_day(alt_day < 100) = NaN;
    
    % To save a bit of time, pre-compute the gaussian weighting
    % functions for the height bins, they'll be the same each time:
    Gauss_z = struct;
    for z = 1:length(zrange)
        Gauss_z(z).vec = exp(- ((alt_day - zrange(z)).^2) ./ (2 * std_z^2));
    end
    
    % we only care about the mean winds here, so just fit the vhorz/az sine
    % wave as we would usually and get the cos and sine terms which give us
    % u and v. no tides, no sliding hour window.
    
    % fit!
    for z = 1:length(zrange)
        
        % only select weights above a threshold:
        w = Gauss_z(z).vec;
        inds = w > 0.05;
        
        [yfit,F] = nph_sinefit(az_day(inds),vhorz_day(inds),360,'weights',w(inds));
        
        % Subscribe the irregular weighted mean altitude location for this
        % bunch of meteors. We can then interpolate onto a
        % regular altitude grid later.
        CompYear.walt(z,d) = wmean(alt_day(inds),w(inds));
        
        % F(1) = cos part, F(2) = sin part, F(3) = mean.
        % so, if my definition of azimuth is correct, F(2) the sine
        % part should be u and F(1) cos part should be v.
        
        % Check for silly numbers:
        F(abs(F) > 200) = NaN;
        
        % Subscribe!
        CompYear.u(z,d) = F(2);
        CompYear.v(z,d) = F(1);
        
    end
    

    
end







return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% ZONAL AND MERID WINDS
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
    Comp.u,...
    Comp.v,...
    CompYear.u,...
    CompYear.v};

clims = {...
    [-50 50],...
    [-20 20],...
    [-50 50],...
    [-20 20]};

clevs = {...
    [-80:5:80],...
    [-30:2.5:30],...
    [-80:5:80],...
    [-30:2.5:30]};

clabs = {...
    ' \bf{u}',' \bf{v}',' \bf{u}',' \bf{v}'};

cbarticks = {...
    -60:20:60,...
    -20:10:20,...
    -60:20:60,...
    -20:10:20};

tits = {...
    'ZONAL WIND',...
    'MERIDIONAL WIND',...
    'ZONAL WIND',...
    'MERIDIONAL WIND'};


% also contour lines:
clines = {...
    [-60 -40 -20 20 40 60],...
    [-20 -10 10 20],...
    [-60 -40 -20 20 40 60],...
    [-20 -10 10 20]};


letters = {'a','c','b','d'};
% letters = {'a','c','a','c'};

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
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
    % and bold zero line
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','w','linewi',4);
    clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000,'color','w','fontweight','bold');
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linewi',2);
    clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000);
    
    
    % LIMITS
%     ylim([77 103])
    ylim(ylims)
    axx.YTick = 60:5:110;
    axx.YMinorTick = 'on';
%     switch ax
%         case {1,2}
            axx.YAxis.MinorTickValues = 60:1:110;
%         case {3,4}
%             axx.YAxis.MinorTickValues = 70:110;
%     end
    
    monthlabels = 0;
    nudgenudge = 0;
    
    switch ax
        case {1,2}
            switch site
                case 'kep'
                    xlim([datenum(2016,02,15) datenum(2020,11,15)])
                case 'sodankyla'
                    xlim([datenum(2008,11,01) datenum(2015,05,01)])    
                otherwise
                    xlim([datenum(min(years),01,01) datenum(max(years)+1,01,01)])
            end
            
            % Just show the Januarys as major ticks
            axx.XTick = datenum(years,01,01);
            axx.XMinorTick = 'on';
            axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
            datetick('x',' ','keepticks','keeplimits')
            
            if monthlabels
                % other months as ticks using text:
                xtix = datenum(years(1),1:(length(years)*12),15);
                for xt = xtix
                    if inrange(xt,xlim) && ~any(xt == axx.XTick)
                        mn = monthname(month(xt));
                        text(xt,75.5+ynudge,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
                    end
                end
            else
                nudgenudge = 3;
            end
            
            axrng = abs([2 1-monthlabels]);
            axrng = sort(axrng(axrng ~= 0));
            if any(ax == axrng)
            % Years as lower bold numbers:
            xtix = datenum(years,07,01);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   yrstr = datestr(xt,'yyyy');
                   text(xt,72.5+nudgenudge-0,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
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
    
    % COLORI
%     cmap = cbrew('RdBu',length(clevs{ax}));
%     cmap = cbrew('nph_RdYlBu',length(clevs{ax}));
%     cmap = cbrew('nph_BuOr',length(clevs{ax}));
%     cmap = cbrew('BrBG',length(clevs{ax}));
%     cmap = flipud(cbrew('PiYG',length(clevs{ax})));
%     cmap = cbrew('nph_Spectral',length(clevs{ax}));
%     cmap = cbrew('nph_modspectral',length(clevs{ax}));
    cmap = nph_saturate(cbrew('nph_modspectral',18),1.2);
%     cmap = cbrew('nph_RdBuPastel',length(clevs{ax}));
%     colormap(gca,cmap)

    % sort out double white in the middle:
%     cmap = [cmap(1:floor(length(clevs{ax})/2),:) ; [1 1 1] ; cmap(floor(length(clevs{ax})/2)+1:end,:)];
    
    % saturate?
%     cmap_hsv = rgb2hsv(cmap);
%     cmap_hsv(:,2) = 1.2.*cmap_hsv(:,2); cmap_hsv(cmap_hsv > 1) = 1;
%     cmap_sat = hsv2rgb(cmap_hsv);
    
    colormap(gca,cmap)
    
    clim(clims{ax})
    
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
            cbar.Ticks = cbarticks{ax};
            cbar.Position = cbar.Position .* [1 1 3 1];
            cbar.Position = cbar.Position +  [cbar.Position(3)/3 0*cbar.Position(4) 0 0];
%             cbar.Label.String = ['\bf{' clabs{ax} '}'];
            cbar.Label.Rotation = 0;
            cbar.Label.VerticalAlignment = 'middle';
%             cbar.Label.HorizontalAlignment = 'left';
%             cbar.Title.String = '   ms^{-1}';
            cbar.Title.String = clabs{ax};
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
%             hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        case 4 % no idea why these are different
            hold on; nph_text([0.025 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
%             hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');

            %             hold on; nph_text([2*0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        
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



return


%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_CompMeanWinds' num2str(min(years)) '-' num2str(max(years)) '_thin'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return



%% MLS TEMPERATURES
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

T = tiledlayout(2,6,'tilespacing','compact');

spans = {...
    [1 5],...
    [1 5],...
    [1 1],...
    [1 1]};

% clims.T = [145 255];
% clevs.T = [100:5:300];
% clines.T = [100:20:300];
% darkclines.T = [200 200];
% clabs.T = '    K';
% cbarticks.T = 100:20:300;
% cbarlabel.T = '\bf{T}';
% tits.T = {'TEMPERATURE','TEMPERATURE','TEMPERATURE','TEMPERATURE'};

% axs = {...
%     [1  6],...
%     [9 14],...
%     [7 8],...
%     [15 16]};

TBP = {...
    Comp.T,...
    Comp.T,...
    CompYear.T,...
    CompYear.T};

clims = {...
    [140 215],...
    [140 215],...
    [140 215],...
    [140 215]};

clevs = {...
    [100:2:300],...
    [100:2:300],...
    [100:2:300],...
    [100:2:300]};

clabs = {...
    '   ','   ','   ','   '};

cbarticks = {...
    110:20:300,...
    110:20:300,...
    110:20:300,...
    110:20:300};

tits = {...
    'TEMPERATURE',...
    'TEMPERATURE',...
    'TEMPERATURE',...
    'TEMPERATURE'};


% also contour lines:
clines = {...
    [100:20:300],...
    [100:20:300],...
    [100:20:300],...
    [100:20:300]};

darkclines = {...
    [180 180],...
    [180 180],...
    [180 180],...
    [180 180]};


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
            X = repmat(Comp.wtime(30,:),55,1); Y = repmat(mlsz,1,12*length(years));
            smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
        case {3,4}
            [X,Y] = meshgrid(datenum(2019,00,days),mlsz);
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
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',300,'color',rgbtrip(.95),'fontweight','bold');
    
    % dark contour lines:
    hold on; [C,h] = contour(X,Y,tbp,darkclines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',300,'color',rgbtrip(.25),'fontweight','bold');
    
    
    % LIMITS
%     ylim([77 103])
    ylim(ylims)
    axx.YTick = 60:5:110;
    axx.YMinorTick = 'on';
%     switch ax
%         case {1,2}
            axx.YAxis.MinorTickValues = 60:1:110;
%         case {3,4}
%             axx.YAxis.MinorTickValues = 70:110;
%     end
    
    
    switch ax
        case {1,2}
            switch site
                case 'kep'
                    xlim([datenum(2016,02,15) datenum(2020,11,15)])
                case 'sodankyla'
                    xlim([datenum(2008,11,01) datenum(2015,05,01)])    
                otherwise
                    xlim([datenum(min(years),01,01) datenum(max(years)+1,01,01)])
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
    
    % COLORI
    cmap = nph_saturate(cbrew('RdBu',18),1);
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
    
    clim(clims{ax})
    
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
            cbar.Ticks = cbarticks{ax};
            cbar.Position = cbar.Position .* [1 1 3 1];
            cbar.Position = cbar.Position +  [cbar.Position(3)/3 0*cbar.Position(4) 0 0];
            cbar.Label.String = [ clabs{ax} ];
            cbar.Label.Rotation = 0;
            cbar.Label.VerticalAlignment = 'middle';
%             cbar.Label.HorizontalAlignment = 'left';
            cbar.Title.String = ' \bf{T}';     
            cbar.TickDirection = 'out';
            
            cbar.Ruler.MinorTick = 'on';
            cbar.Ruler.MinorTickValues = 100:10:300;
            
%             cbar.Limits = [140 220];
            
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
%             hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        case 4 % no idea why these are different
            hold on; nph_text([0.025 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
%             hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');

            %             hold on; nph_text([2*0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        
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



return


%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_MLSTemperature' num2str(min(years)) '-' num2str(max(years)) '_thin'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return



%% MLS TEMPERATURE GRADIENTS WITH LATITUDE
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

T = tiledlayout(2,6,'tilespacing','compact');

spans = {...
    [1 5],...
    [1 5],...
    [1 1],...
    [1 1]};

% clims.T = [145 255];
% clevs.T = [100:5:300];
% clines.T = [100:20:300];
% darkclines.T = [200 200];
% clabs.T = '    K';
% cbarticks.T = 100:20:300;
% cbarlabel.T = '\bf{T}';
% tits.T = {'TEMPERATURE','TEMPERATURE','TEMPERATURE','TEMPERATURE'};

% axs = {...
%     [1  6],...
%     [9 14],...
%     [7 8],...
%     [15 16]};

TBP = {...
    Comp.Tdiff,...
    Comp.Tdiff,...
    CompYear.Tdiff,...
    CompYear.Tdiff};

clims = {...
    [-1.2 1.2],...
    [-1.2 1.2],...
    [-1.2 1.2],...
    [-1.2 1.2]};

clevs = {...
    [-5:0.1:5],...
    [-5:0.1:5],...
    [-5:0.1:5],...
    [-5:0.1:5]};

clabs = {...
    '   ','   ','   ','   '};

cbarticks = {...
    [-5:0.5:5],...
    [-5:0.5:5],...
    [-5:0.5:5],...
    [-5:0.5:5]};

tits = {...
    'MERIDIONAL TEMPERATURE GRADIENT',...
    'MERIDIONAL TEMPERATURE GRADIENT',...
    'MERIDIONAL TEMPERATURE GRADIENT',...
    'MERIDIONAL TEMPERATURE GRADIENT'};


% also contour lines:
clines = {...
    [-5:0.5:5],...
    [-5:0.5:5],...
    [-5:0.5:5],...
    [-5:0.5:5]};

darkclines = {...
    [0 0],...
    [0 0],...
    [0 0],...
    [0 0]};


letters = {'g','g','h','h'};

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
            X = repmat(Comp.wtime(30,:),55,1); Y = repmat(mlsz,1,12*length(years));
            smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
        case {3,4}
            [X,Y] = meshgrid(datenum(2019,00,days),mlsz);
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
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',300,'color',rgbtrip(.95),'fontweight','bold');
    
    % dark contour lines:
    hold on; [C,h] = contour(X,Y,tbp,darkclines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',300,'color',rgbtrip(.25),'fontweight','bold');
    
    
    % LIMITS
%     ylim([77 103])
    ylim(ylims)
    axx.YTick = 60:5:110;
    axx.YMinorTick = 'on';
%     switch ax
%         case {1,2}
            axx.YAxis.MinorTickValues = 60:1:110;
%         case {3,4}
%             axx.YAxis.MinorTickValues = 70:110;
%     end
    
    
    switch ax
        case {1,2}
            switch site
                case 'kep'
                    xlim([datenum(2016,02,15) datenum(2020,11,15)])
                case 'sodankyla'
                    xlim([datenum(2008,11,01) datenum(2015,05,01)])    
                otherwise
                    xlim([datenum(min(years),01,01) datenum(max(years)+1,01,01)])
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
    
    % COLORI
    cmap = nph_saturate(cbrew('RdBu',18),1);
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
    
    clim(clims{ax})
    
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
            cbar.Ticks = cbarticks{ax};
            cbar.Position = cbar.Position .* [1 1 3 1];
            cbar.Position = cbar.Position +  [cbar.Position(3)/3 0*cbar.Position(4) 0 0];
            cbar.Label.String = [ clabs{ax} ];
            cbar.Label.Rotation = 0;
            cbar.Label.VerticalAlignment = 'middle';
%             cbar.Label.HorizontalAlignment = 'left';
            cbar.Title.String = '  \bf{dT/d\phi}';     
            cbar.TickDirection = 'out';
            
            cbar.Ruler.MinorTick = 'on';
            cbar.Ruler.MinorTickValues = 100:10:300;
            
%             cbar.Limits = [140 220];
            
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
%             hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        case 4 % no idea why these are different
            hold on; nph_text([0.025 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
%             hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');

            %             hold on; nph_text([2*0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        
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



return


%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_MLSTemperatureGradient' num2str(min(years)) '-' num2str(max(years)) '_thin'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return






































