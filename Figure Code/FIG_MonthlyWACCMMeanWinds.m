
%%

%%%% EDIT: trying to standardise figure style for the KEP paper...

% Plot mean winds from WACCM for KEP.
% Currently using WACCM specified dynamics winds from Anne Smith. Cheers
% Anne!

%%%% ALSO compare a composite year of winds from KEP to an average year
%%%% from WACCM.

%%%% ALSO LOAD AND COMPARE MLS TEMPERATURES?




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SITE AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

site = 'kep';

years = 2000:2021;

%%%% LOAD WACCM:
disp(['Loading WACCM...'])

nc = getnet(['/Users/neil/data/WACCM/' site '/WACCM_monthly_2000-2014_324E_54S.nc']);
Waccm = struct;
Waccm.u = sq(nanmean(nc.Data.u,2));
Waccm.v = sq(nanmean(nc.Data.v,2));
Waccm.T = sq(nanmean(nc.Data.T,2));
Waccm.wtime = datenum(num2str(nc.Data.month),'yyyymm') + 15; % centre on the 15th of each month
Waccm.walt = nc.Data.altitude;

%%%% also create an average year:
Waccm.CompYear.u        = nanmean(reshape(Waccm.u,[51 12 15]),3);
Waccm.CompYear.v        = nanmean(reshape(Waccm.v,[51 12 15]),3);
Waccm.CompYear.T        = nanmean(reshape(Waccm.T,[51 12 15]),3);
Waccm.CompYear.wtime    = datenum(2019,1:12,15);
Waccm.CompYear.walt     = Waccm.walt;

% %
% % % First, get some composite of HWD stuff anyway...
% % Waccm.u = [];
% % Waccm.v = [];
% % Waccm.gwvar = [];
% % Waccm.walt = [];
% % Waccm.wtime = [];
% % for yr = years
% %     yrstr = num2str(yr);
% %     load(['/Users/neil/data/MeteorRadar/' site '/matlab/hwd/' yrstr '_' site '_hwd_gaussbins.mat'])
% %     Waccm.u      = cat(2,Waccm.u,sq(nanmean(HWD.Data.MonthlyComp.u,2)));
% %     Waccm.v      = cat(2,Waccm.v,sq(nanmean(HWD.Data.MonthlyComp.v,2)));
% %     Waccm.walt = cat(2,Waccm.walt,sq(nanmean(HWD.Data.MonthlyComp.walt,2)));
% %     Waccm.wtime = cat(2,Waccm.wtime,sq(nanmean(HWD.Data.MonthlyComp.wtime,2)));
% % end
% %

%%%% LOAD RADAR:

% load all the mpd files for the year range and cat them into a massive
% array.
disp(['Loading RADAR...'])
mpddirec = ['/Users/neil/data/MeteorRadar/' site '/matlab/'];

allMPD = struct;
dayrange = datenum(2016,01,01):datenum(2020,12,31);

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
%%%% LOAD MLS

mlsyears = 2016:2020;
mlsT = [];
mlstime = [];

for yr = mlsyears
    
    yrstr = num2str(yr);
    
    disp(['Loading MLS for ' yrstr '...'])
    
    load(['/Users/neil/data/MLS/' yrstr '_MLS_T.mat']);
    
    %%%% select region:
    inds = ...
        inrange(MLS.Data.Longitude,[-45 -30]) & ...
        inrange(MLS.Data.Latitude ,[-60 -50]);
    
    TT = MLS.Data.T(:,inds);
    timtim = MLS.Data.Time(inds);
    mlsz = MLS.Data.Altitude;
    
    clear MLS
    
    %%%% stick it together:
    mlsT = cat(2,mlsT,TT);
    mlstime = cat(1,mlstime,timtim);
    
end

% %%%% now assemble it into an average year:
% 
% %%%%for each month...
% CompYear.T          = nan(55,12);
% % CompYear.mlstime    = datenum(min(mlsyears),01,15):datenum(max(mlsyears),12,15);
% CompYear.mlsalt     = mlsz;
% for mn = 1:12
%     minds = month(mlstime) == mn;
%     CompYear.T(:,mn) = nanmean(mlsT(:,minds),2);
% end


% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Recalculate a rolling 30-day composite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Assembling an all-years composite for the specified years...')

method = 'short'; % long | short
% i only find the tiniest of difference between the results for each

days = ceil(linspace(1,365,60));

daywindow = 30;
std_z = 1.275;      % FWHM for height gaussian
std_time = 0.85;    % FWHM for time gaussian
zrange = 70:1:105;
% zrange = 70:1:110;

CompYear = struct;
CompYear.u = nan(length(zrange),length(days));
CompYear.v = nan(length(zrange),length(days));
CompYear.walt = nan(length(zrange),length(days));

dn_mls = daynumber(mlstime);
CompYear.T = nan(55,length(days));

for d = 1:length(days)
    
    dd = days(d);
    disp(num2str(dd))
    drange = (dd-daywindow/2):(dd+daywindow/2);
    
    % try dealing with wraparound in a different way:
    drange(drange < 0) = 365+drange(drange < 0);
    drange(drange > 365) = drange(drange > 365) - 365;
    timerange = ismember(dn,drange);
    
    timerange_mls = ismember(dn_mls,drange);
    
    CompYear.T(:,d) = nanmean(mlsT(:,timerange_mls),2);
    
    %     % deal with wraparound for years
    %     drange_neg = drange(drange <= 0);
    %     drange = drange(drange > 0);
    %     drange_pos = drange(drange >= 366);
    %     drange = drange(drange < 366);
    %
    %     % if there a low wraparoud:
    %     if ~isempty(drange_neg)
    %         drange_neg = 365+drange_neg;
    %         timerange = inrange(dn,min(drange_neg),366) | inrange(dn,minmax(drange));
    %     else
    %         timerange = inrange(dn,minmax(drange));
    %     end
    %     % if there a low wraparoud:
    %     if ~isempty(drange_pos)
    %         drange_pos = drange_pos-365;
    %         timerange = inrange(dn,minmax(drange_pos)) | inrange(dn,minmax(drange));
    %     else
    %         timerange = inrange(dn,minmax(drange));
    %     end
    
    % select todays composite meteors:
    az_day = az(timerange);
    vhorz_day = vhorz(timerange);
    alt_day = alt(timerange);
    tim_day = tim(timerange);
    dn_day = dn(timerange);
    
    hour_day = mod(tim_day,1);
    
    % To save a bit of time, pre-compute the gaussian weighting
    % functions for the height bins, they'll be the same each time:
    Gauss_z = struct;
    for z = 1:length(zrange)
        Gauss_z(z).vec = exp(- ((alt_day - zrange(z)).^2) ./ (2 * std_z^2));
    end
    
    
    % we only care about the mean winds here, so just fit the vhorz/az sine
    % wave as we would usually and get the cos and sine terms which give us
    % u and v. no tides, no sliding hour window.
    
    %     if d == 1
    %         return
    %     end
    
    switch method
        case 'long'
            % We should also try to introduce something to deal with the fact
            % that so many more meteors come from one time of day that others. This
            % means that, becuase there are so many more meteors from e.g. the
            % morning, then the winds will look more like the morning winds rather
            % than the daily average.
            Gauss_h = struct;
            for h = 0:23
                Gauss_h(h+1).vec = exp(- ((hour_day - (h/24)).^2) ./ (2 * std_time^2));
            end
            %%%% try binning into a composite day of u and v?
            CompDay.u       = nan(length(zrange),24);
            CompDay.v       = nan(length(zrange),24);
            hbins = (0:23)./24;
            for z = 1:length(zrange)
                for h = 0:23
                    w = Gauss_z(z).vec .* Gauss_h(h+1).vec;
                    inds = w > 0.05;
                    [yfit,F] = nph_sinefit(az_day(inds),vhorz_day(inds),360,'weights',w(inds));
                    F(abs(F) > 200) = NaN;
                    CompDay.u(z,h+1) = F(2);
                    CompDay.v(z,h+1) = F(1);
                end
                CompYear.walt(z,d) = wmean(...
                    alt_day(Gauss_z(z).vec > 0.05),...
                    Gauss_z(z).vec(Gauss_z(z).vec > 0.05));
            end
            
            %%%% Now fit sine waves to this to get the average:
            for z = 1:length(zrange)
                [~,F] = nph_sinefit(0:23,CompDay.u(z,:),24);
                F(abs(F) > 200) = NaN;
                CompYear.u(z,d) = F(3);
                [~,F] = nph_sinefit(0:23,CompDay.v(z,:),24);
                F(abs(F) > 200) = NaN;
                CompYear.v(z,d) = F(3);
            end
            
            
        case 'short'
            
            % just take all meteors detected on a day and fit them. Don't
            % worry about the fact that there will be more meteors at some
            % times of day than others.
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
                % Check for silly numbers:
                F(abs(F) > 200) = NaN;
                % Subscribe!
                CompYear.u(z,d) = F(2);
                CompYear.v(z,d) = F(1);
            end
            
    end
    
end % next day


%%%% WACCM vertical pressure levels:

Waccm.pres = reverse(linearise([5.9603000e-06   9.8269002e-06   1.6201851e-05   2.6712251e-05   4.4041000e-05
   7.2612748e-05   0.00011971900   0.00019738000   0.00032542250   0.00053653251
   0.00088460249    0.0014584575    0.0024045750    0.0039782499    0.0065568256
     0.010813826     0.017898000     0.029557750     0.048730750     0.079910751
      0.12827325      0.19812000      0.29202500      0.41016750      0.55346999
      0.73047998      0.95594750       1.2447950       1.6128500       2.0793250
       2.6674250       3.4048751       4.3245750       5.4653999       6.8728500
       8.5997251       10.707050       13.264750       16.351751       20.056751
       24.479000       29.727999       35.923250       43.193750       51.677499
       61.520498       73.750958       87.821230       103.31713       121.54724
       142.99404       168.22508       197.90809       232.82862       273.91082
       322.24190       379.10090       445.99257       524.68717       609.77869
       691.38943       763.40448       820.85837       859.53477       887.02025
       912.64455       936.19840       957.48548       976.32541       992.55610]'));
   




return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% PLOTTING - WACCM WINDS TIME SERIES AGAINST ALTITUDE
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

T = tiledlayout(2,32,'tilespacing','compact');

uv = 'T';

spans = {...
    [1 32],...
    [1 24],...
    [1 5]};

TBP = {...
    Waccm.(uv),...
    Waccm.(uv),...
    repmat(Waccm.CompYear.(uv),1,3)};

switch uv
    case 'u'
        clims = [-50 50];
        clevs = [-60:5:60];
        clabs = 'u';
        cbarticks = -60:20:60;
        tits = {'WACCM, ZONAL WIND',' ',' '};
        clines = [-80 -60 -40 -20 20 40 60 80];
        blackclines = [-100 100];
        clabelspacing = [1000 1000 200];
        letters = {'a',' ','b'};
        textcolor = 'k';
        textbordercolor = 'w';
    case 'v'
        clims = [-20 20];
        clevs = [-20:2.5:20];
        clabs = 'v';
        cbarticks = -40:10:40;
        tits = {'WACCM, MERIDIONAL WIND',' ',' '};
        clines = [-30 -20 -10 10 20 30 40];
        blackclines = [-100 100];
        clabelspacing = [1000 1000 200];
        letters = {'c',' ','d'};
        textcolor = 'k';
        textbordercolor = 'w';
    case 'T'
        clims = [160 240];
        clevs = [100:5:280];
        clabs = '  T';
        cbarticks = 100:20:280;
        tits = {'WACCM, TEMPERATURE',' ',' '};
        clines = [100:20:280];
        blackclines = [200 200];
        clabelspacing = [500 500 200];
        letters = {'e',' ','f'};
        textcolor = 'k';
        textbordercolor = 'w';
end

% xlims = {...
%     [datenum(2000,01,15) datenum(2008,02,15)],...
%     [datenum(2008,02,15) datenum(2014,12,15)],...
%     [datenum(2019,01,01) datenum(2020,01,01)]};
xlims = {...
    [datenum(2001,01,15) datenum(2009,01,01)],...
    [datenum(2009,01,15) datenum(2014,12,15)],...
    [datenum(2019,01,01) datenum(2020,01,01)]};

ynudge = -12;
ylims = [77+ynudge 110];


for ax = [1 2 3]
    
    %     subplot(rows,cols,axs{ax})
    
    % blank space?
    if ax == 3
        nexttile([1 1])
        set(gca,'color','none','xcolor','none','ycolor','none');
    end
    
    nexttile(spans{ax})
    
    axx = gca;
    
    % choose wot 2 plot n wot not
    switch ax
        case {1,2}
            X = Waccm.wtime; Y = Waccm.walt;
            smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
        case {3,4}
            %             [X,Y] = meshgrid(datenum(2019,00,days),nanmean(CompYear.walt,2));
            [X,Y] = meshgrid(datenum(2018,1:36,15),Waccm.CompYear.walt);
            smoo1 = [0.25 0.75]; smoo2 = [0.25 0.75];
    end
    
    
    
    
    % allll the contours...
    tbp = TBP{ax}; %nanlocs = isnan(tbp); tbp(nanlocs) = nanmean(tbp(:));
    %     tbp = imgaussfilt(imgaussfilt(tbp,smoo1),smoo2);
    %     tbp(nanlocs) = NaN;
    
    % got a problem where the clabels are drawing outside the axes limits,
    % then there are not very many labels in the axes limits. Going to set
    % some lower altitudes to NaN and see what happens
    switch ax
        case {1,2}
            tbp(~inrange(nc.Data.altitude,ylims+pm(1)),:) = NaN;
        case 3
            tbp(~inrange(nc.Data.altitude,ylims+pm(1)),:) = NaN;
            tbp(~inrange(X,[datenum(2018,12,01) datenum(2020,02,01)])) = NaN;
    end
    
    % filled contours:
    hold on; contourf(X,Y,tbp,clevs,'edgecolor','none');
    % contour lines:
    hold on; [C,h] = contour(X,Y,tbp,clines,'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',clabelspacing(ax),'color',rgbtrip(.95),'fontweight','bold');
    %     % minor contour lines:
    %     hold on; [C,h] = contour(X,Y,tbp,minorclines,'edgecolor',rgbtrip(.25),'linewi',1.25);
    % and bold zero line
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','w','linewi',4);
    clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000,'color','w','fontweight','bold');
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linewi',2);
    clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000);
    
    % and extra contour lines:
    hold on; [C,h] = contour(X,Y,tbp,blackclines,'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',clabelspacing(ax),'color',rgbtrip(.25),'fontweight','bold');
    
    
    % LIMITS
    ylim(ylims)
    %     ylim([77 103])
    axx.YTick = 60:10:110;
    axx.YMinorTick = 'on';
    %     switch ax
    %         case {1,2}
    axx.YAxis.MinorTickValues = 60:2:110;
    %         case {3,4}
    %             axx.YAxis.MinorTickValues = 70:110;
    %     end
    
    
    switch ax
        case {1,2}
            %             xlim([datenum(2016,02,15) datenum(2020,11,15)])
            xlim(xlims{ax})
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
            %             if ax == 2
            % Years as lower bold numbers:
            xtix = datenum(years,07,01);
            for xt = xtix
                if inrange(xt,xlim) && ~any(xt == axx.XTick)
                    yrstr = datestr(xt,'yyyy');
                    text(xt,72.5+ynudge-1,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
                end
            end
            %             end
            
        case 3
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
                    text(xt,76+ynudge,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
                end
            end
            %             axx.XTick = datenum(2019,1:13,01);
            %             datetick('x','m','keepticks','keeplimits')
            %             xlim(datenum(2019,[1 13],01))
            if ax == 3
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
            cbar.Label.String = ['\bf{' clabs '}'];
            cbar.Label.Rotation = 0;
            cbar.Label.VerticalAlignment = 'middle';
            %             cbar.Label.HorizontalAlignment = 'left';
            switch uv
                case 'T'
                    cbar.Title.String = 'K';
                otherwise
                    cbar.Title.String = '   ms^{-1}';
            end
    end
    
    
    
    % AXES LETTER
    %letter = alphabet(ax);
    letter = letters{ax};
    switch ax
        case 1
            hold on; nph_text([0.005 0.85],['(' letter ')'],'fontsize',1.5*fs,'color',textcolor,'textborder',textbordercolor,'horizontalalignment','left','verticalalignment','middle');
            hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','color',textcolor,'textborder',textbordercolor,'horizontalalignment','left','verticalalignment','middle');
        case 2
%             hold on; nph_text([0.005 0.83],['(' letter ')'],'fontsize',1.5*fs,'color',textcolor,'textborder',textbordercolor,'horizontalalignment','left','verticalalignment','middle');
%             hold on; nph_text([0.005 0.82],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','color',textcolor,'textborder',textbordercolor,'horizontalalignment','left','verticalalignment','middle');
        case 3 % no idea why these are different
            hold on; nph_text([0.085 0.83],['(' letter ')'],'fontsize',1.5*fs,'color',textcolor,'textborder',textbordercolor,'horizontalalignment','left','verticalalignment','middle');
            
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
    
    
    % blank space?
    if ax == 3
        nexttile([1 1])
        set(gca,'color','none','xcolor','none','ycolor','none');
    end
    
    
    
    
end



return


%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_WACCMMeanWinds_' uv '_thin'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% WACCM/RADAR COMPOSITE YEAR COMPARISON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([1 0.7125])

%-------------------------------------------------------
vert_gap = 0.08;        horz_gap = 0.03;
lower_marg = 0.075;     upper_marg = 0.025;
left_marg = 0.065;      right_marg = 0.075;

rows = 2; cols = 8;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 24;

T = tiledlayout(2,23,'tilespacing','compact');

spans = {...
    [1 10],...
    [1 10],...
    [1 10],...
    [1 10]};

uv = {'u','u','v','v'};
% uv = {'T','T','T','T'};

TBP = {...
    CompYear.(uv{1}),...
    repmat(Waccm.CompYear.(uv{2}),1,3),...
    CompYear.(uv{3}),...
    repmat(Waccm.CompYear.(uv{4}),1,3)};

clims.u = [-50 50];
clevs.u = [-60:5:60];
clines.u = [-80 -60 -40 -20 20 40 60 80];
darkclines.u = [-100 100];
clabs.u = '    ms^{-1}';
cbarticks.u = -60:20:60;
cbarminorticks.u = -60:20:60;
cbarlabel.u = '\bf{u}';
tits.u = {'ZONAL','ZONAL','MERIDIONAL','MERIDIONAL'};

clims.v = [-20 20];
clevs.v = [-20:2.5:20];
clines.v = [-30 -20 -10 10 20 30 40];
darkclines.v = [-100 100];
clabs.v = '    ms^{-1}';
cbarticks.v = -40:10:40;
cbarminorticks.v = -40:10:40;

cbarlabel.v = '\bf{v}';
tits.v = {'ZONAL','ZONAL','MERIDIONAL','MERIDIONAL'};

clims.T = [145 255];
clevs.T = [100:5:300];
clines.T = [100:20:300];
darkclines.T = [200 200];
clabs.T = '    K';
cbarticks.T = 100:20:300;
cbarminorticks.T = 100:10:300;
cbarlabel.T = '\bf{T}';
tits.T = {'TEMPERATURE','TEMPERATURE','TEMPERATURE','TEMPERATURE'};


% xlims = {...
%     [datenum(2000,01,15) datenum(2008,02,15)],...
%     [datenum(2008,02,15) datenum(2014,12,15)],...
%     [datenum(2019,01,01) datenum(2020,01,01)]};
xlims = {...
    [datenum(2001,01,15) datenum(2009,01,01)],...
    [datenum(2009,01,15) datenum(2014,12,15)],...
    [datenum(2019,01,01) datenum(2020,01,01)]};

ynudge = -17;
ylims = [77+ynudge 107.5];

if strcmpi(uv{1},'T')
    letters = {'e','f','e','f'};
else
    letters = {'a','b','c','d'};
end

Axes = struct;

for ax = 1:4
    
    %     subplot(rows,cols,axs{ax})
    
%     % blank space?
%     switch ax
%         case {1,3}
%         nexttile([1 1])
%         set(gca,'color','none','xcolor','none','ycolor','none');
%     end
    
    nexttile(spans{ax})
    
    axx = gca;
    Axes(ax).axx = axx;
    
    % choose wot 2 plot n wot not
    switch ax
        case {2,4}
            %             X = Waccm.CompYear.wtime; Y = Waccm.CompYear.walt;
            [X,Y] = meshgrid(datenum(2018,1:36,15),Waccm.CompYear.walt);
            smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
        case {1,3}
            [X,Y] = meshgrid(datenum(2019,00,days),nanmean(CompYear.walt,2));
            %             [X,Y] = meshgrid(datenum(2018,1:36,15),Waccm.CompYear.walt);
            smoo1 = [0.25 0.75]; smoo2 = [0.25 0.75];
            % overwrite this if T...
            switch uv{ax}
                case 'T'
                    [X,Y] = meshgrid(datenum(2019,00,days),mlsz);
                    smoo1 = [0.25 0.75]; smoo2 = [0.25 0.75];
            end
    end
    
    % allll the contours...
    tbp = TBP{ax}; %nanlocs = isnan(tbp); tbp(nanlocs) = nanmean(tbp(:));
    %     tbp = imgaussfilt(imgaussfilt(tbp,smoo1),smoo2);
    %     tbp(nanlocs) = NaN;
    
    % got a problem where the clabels are drawing outside the axes limits,
    % then there are not very many labels in the axes limits. Going to set
    % some lower altitudes to NaN and see what happens
    if ~strcmpi(uv{ax},'T')
        switch ax
            case {2,4}
                tbp(~inrange(nc.Data.altitude,ylims+pm(1)),:) = NaN;
                tbp(~inrange(X,[datenum(2018,12,01) datenum(2020,02,01)])) = NaN;
            case {1,3}
                % tbp(~inrange(zrange,[73 103]),:) = NaN;
        end
    end
    
    %%%% ticks first...
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
            text(xt,76+ynudge,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
        end
    end
    
    % LIMITS
    ylim(ylims)
    %     ylim([77 103])
    axx.YTick = 60:5:110;
    axx.YMinorTick = 'on';
    axx.YAxis.MinorTickValues = 60:1:110;
    
    %%%% grid lines under the contours...
    for i = 1:length(axx.XTick)
        hold on; plot(axx.XTick([i i]),axx.YLim,'color',[0 0 0 0.15],'linewi',1);
    end
    
    % filled contours:
    hold on; contourf(X,Y,tbp,clevs.(uv{ax}),'edgecolor','none');
    % contour lines:
    hold on; [C,h] = contour(X,Y,tbp,clines.(uv{ax}),'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',500,'color',rgbtrip(.95),'fontweight','bold');
    %     % minor contour lines:
    %     hold on; [C,h] = contour(X,Y,tbp,minorclines,'edgecolor',rgbtrip(.25),'linewi',1.25);
    % and bold zero line
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','w','linewi',4);
    clabel(C,h,'fontsize',0.8*fs,'labelspacing',500,'color','w','fontweight','bold');
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linewi',2);
    clabel(C,h,'fontsize',0.8*fs,'labelspacing',500);
    
    % darkcontour lines:
    hold on; [C,h] = contour(X,Y,tbp,darkclines.(uv{ax}),'edgecolor',rgbtrip(.25),'linewi',1.25);
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',500,'color',rgbtrip(.25),'fontweight','bold');
    
    
 
%     % LIMITS
%     ylim(ylims)
%     %     ylim([77 103])
%     axx.YTick = 60:5:110;
%     axx.YMinorTick = 'on';
%     %     switch ax
%     %         case {1,2}
%     axx.YAxis.MinorTickValues = 60:1:110;
%     %         case {3,4}
%     %             axx.YAxis.MinorTickValues = 70:110;
%     %     end
    
    
    %     switch ax
    %         case {1,2}
    % %             xlim([datenum(2016,02,15) datenum(2020,11,15)])
    %             xlim(xlims{ax})
    %             % Just show the Januarys as major ticks
    %             axx.XTick = datenum(years,01,01);
    %             axx.XMinorTick = 'on';
    %             axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
    %             datetick('x',' ','keepticks','keeplimits')
    % %             axx.XTickLabel = {};
    %             % other months as ticks using text:
    %             xtix = datenum(years(1),1:(length(years)*12),15);
    %             for xt = xtix
    %                if inrange(xt,xlim) && ~any(xt == axx.XTick)
    %                    mn = monthname(month(xt));
    %                    text(xt,75.5+ynudge,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
    %                end
    %             end
    % %             if ax == 2
    %             % Years as lower bold numbers:
    %             xtix = datenum(years,07,01);
    %             for xt = xtix
    %                if inrange(xt,xlim) && ~any(xt == axx.XTick)
    %                    yrstr = datestr(xt,'yyyy');
    %                    text(xt,72.5+ynudge-1,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
    %                end
    %             end
    % %             end
    %
    %         case 3
    %             axx.XTick = datenum(2019,1:13,01);
    %             datetick('x','m','keepticks','keeplimits')
    %             xlim(datenum(2019,[1 13],01))
    
    %     end
    
    switch ax
        case 3
            xlabel('Average (Composite) Year','fontsize',0.8*fs,'fontweight','bold')
        case 4
            xlabel('Average Year','fontsize',0.8*fs,'fontweight','bold')
    end
    
    %%%% overwrite if T and do labels for all...
    if strcmpi(uv{ax},'T')
        xlabel('Average Year','fontsize',0.8*fs,'fontweight','bold')
    end
    
    

%     switch ax
%         case {1,3}
            ylabel('Altitude (km)')
%     end
    switch ax
        case 2
            title('WACCM','fontsize',fs);
        case 1
            if strcmpi(uv{ax},'T')
                title('MLS','fontsize',fs)  
            else
                title('RADAR','fontsize',fs)
            end
    end
    
    drawnow;
    
    % line around the axis
    hold on; plot(axx.XLim([1 2 2]),axx.YLim([2 2 1]),'color','k','linewi',1.5)

    
    % COLORI
    %     cmap = cbrew('RdBu',length(clevs{ax}));
    %     cmap = cbrew('nph_RdYlBu',length(clevs{ax}));
    %     cmap = cbrew('nph_BuOr',length(clevs{ax}));
    %     cmap = cbrew('BrBG',length(clevs{ax}));
    %     cmap = flipud(cbrew('PiYG',length(clevs{ax})));
    %     cmap = cbrew('nph_Spectral',length(clevs{ax}));
    %     cmap = cbrew('nph_modspectral',length(clevs{ax}));
    if strcmpi(uv{ax},'T')
        cmap = cbrew('RdBu',length(clines.(uv{ax})));
    else
        cmap = nph_saturate(cbrew('nph_modspectral',18),1.2);
    end
    %     cmap = cbrew('nph_RdBuPastel',length(clevs{ax}));
    %     colormap(gca,cmap)
    
    % sort out double white in the middle:
    %     cmap = [cmap(1:floor(length(clevs{ax})/2),:) ; [1 1 1] ; cmap(floor(length(clevs{ax})/2)+1:end,:)];
    
    % saturate?
    %     cmap_hsv = rgb2hsv(cmap);
    %     cmap_hsv(:,2) = 1.2.*cmap_hsv(:,2); cmap_hsv(cmap_hsv > 1) = 1;
    %     cmap_sat = hsv2rgb(cmap_hsv);
    
    colormap(gca,cmap)
    
    clim(clims.(uv{ax}))
    
    % COLORBAR
%     switch ax
        %         case {1,2}
        %             cbar = nph_colorbar;
        %             cbar.Ticks = cbarticks{ax};
        %             cbar.Position = cbar.Position .* [0.99 1 0.5 1];
        %             cbar.Label.String = ['\bf{' clabs{ax} '}'];
        %             cbar.Label.Rotation = 0;
        %             cbar.Label.VerticalAlignment = 'middle';
        % %             cbar.Label.HorizontalAlignment = 'left';
        %             cbar.Title.String = '   ms^{-1}';
%         case {3,4}
%             cbar = nph_colorbar;
%             cbar.Ticks = cbarticks.(uv{ax});
%             cbar.Position = cbar.Position .* [1 1 1.5 1];
% %             cbar.Position = cbar.Position +  [cbar.Position(3)/3 0*cbar.Position(4) 0 0];
%             cbar.Label.String = ['\bf{' clabs.(uv{ax}) '}'];
%             cbar.Label.Rotation = 0;
%             cbar.Label.VerticalAlignment = 'middle';
%             %             cbar.Label.HorizontalAlignment = 'left';
%             cbar.Title.String = '   ms^{-1}';
%     end
    
    
    
    % AXES LETTER
    %letter = alphabet(ax);
    letter = letters{ax};
%     switch ax
%         case 1
%             hold on; nph_text([0.005 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
%             hold on; nph_text([0.005 0.84],['       ' tits.(uv{ax})],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
%         case 2
%             hold on; nph_text([0.005 0.83],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
%             hold on; nph_text([0.005 0.82],['       ' tits.(uv{ax})],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
%         case 3 % no idea why these are different
%             hold on; nph_text([0.085 0.83],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
%             
%     end
    
    
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
%     set(gca,'linewi',1.5,'tickdir','out')
    
%     %%%%% blank space?
    switch ax
        case {1,3}
        nexttile([1 2])
        set(gca,'color','none','xcolor','none','ycolor','none');
        case {2,4}
        nexttile([1 1])
        set(gca,'color','none','xcolor','none','ycolor','none');
    end
    
    
    
    
end

% do the colorbars outside the looping as they move the axes around
for ax = 1:length(Axes)
    axes(Axes(ax).axx);
    
    cbar = nph_colorbar;
    cbar.Ticks = cbarticks.(uv{ax});
    cbar.Position = cbar.Position .* [1 1 1.5 1];
    cbar.Position = cbar.Position +  [0.25*cbar.Position(3) 0 0 0];
    cbar.Label.String = [ clabs.(uv{ax}) ];
    cbar.Label.Rotation = 0;
    cbar.LineWidth = 1.5;
    cbar.Label.VerticalAlignment = 'middle';
    cbar.Title.String = ['   ' cbarlabel.(uv{ax})];
    cbar.FontSize = fs;
    cbar.TickDirection = 'out';
    
    cbar.Ruler.MinorTick = 'on';
    cbar.Ruler.MinorTickValues = cbarminorticks.(uv{ax});
    
    % axes letters too
    letter = letters{ax};
    hold on; nph_text([0.005 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
%     hold on; nph_text([0.005 0.84],['       ' tits.(uv{ax}){ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');

end

for ax = 1:length(Axes)
    %%%% add waccm pressure level tick marks:
    switch ax
        case {2,4}
            axx2 = axes('position',Axes(ax).axx.Position,...
                'color','none','xcolor','none',...
                'xlim',Axes(ax).axx.XLim,'ylim',Axes(ax).axx.YLim,...
                'yaxislocation','right','linewi',1.5);
            axx2.YTick      = sort(round(pres2alt(Waccm.pres),1));
            axx2.YTickLabel = {};
            axx2.TickDir = 'out';
    end
end





return


%% EXPORT??? ==============================================================

if strcmpi(uv{ax},'T')
    savename = ['~/Desktop/' upper(site) '_WACCM_Comparison_T'];    
else
    savename = ['~/Desktop/' upper(site) '_WACCM_Comparison_uv'];
end

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')


%%


return











































