

% HISTOGRAMS



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SITE AND TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

site = 'kep';

years = 2016:2020;

% load all the mpd files for the year range and cat them into a
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
        allMPD(d).x                     = MPD.Data.x;
        allMPD(d).y                     = MPD.Data.y;
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
xdist   = cat(1,allMPD(:).x);
ydist   = cat(1,allMPD(:).y);

vrad    = cat(1,allMPD(:).RadialVelocity);
vhorz   = vrad ./ cosd(zen);

% % % % Apply Zenith Angle Limits:
% % % MWD.Thresholds.ZenithLimits = [15 65]; % use [15 65] for computing mean winds, impose stricter for GWs.
% % % inds = zen >= MWD.Thresholds.ZenithLimits(1) & zen <= MWD.Thresholds.ZenithLimits(2);

% avoid silly horz velocities:
% inds = inds & vhorz < 200;
inds = vhorz < 200;

az = single(az(inds));
tim = double(tim(inds));
alt = single(alt(inds));
vrad = single(vrad(inds));
vhorz = single(vhorz(inds));
zen = single(zen(inds));
gwres = single(gwres(inds));
xdist = single(xdist(inds));
ydist = single(ydist(inds));

hourofday = mod(tim,1) * 24;
dn = daynumber(tim);

ndays = length(unique(floor(tim)));

% %%%% bin for range plot:
% xdist = rng.*sind(az);
% ydist = rng.*cosd(az);


%%%% compute avg nmets per day each month

mrange = datenum(2016,1:(12*length(years)+1),01);

metsperday = nan(1,length(mrange)-1);

for m = 1:length(mrange)-1
    inds = inrange(tim,[mrange(m) mrange(m+1)]);
    lendays = length(mrange(m):mrange(m+1));
    metsperday(m) = floor(nansum(inds) ./ lendays);
end





return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTY PLOT!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([0.75 0.8])

%-------------------------------------------------------
vert_gap = 0.075;       horz_gap = 0.075;
lower_marg = 0.075;     upper_marg = 0.05;
left_marg = 0.06;       right_marg = 0.04;

rows = 3; cols = 10;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 16;


barcolor = mcolor(1);
clen = 12;
cmap = nph_saturate(cbrew('nph_CyclicRainbow',clen),0.8);


axs = {...
    [1:4 11:14],...
    5:7,...
    8:10,...
    15:17,...
    18:20,...
    [21:30]};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Meteor map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(rows,cols,axs{1})

switch site
    case 'kep'
        clat = -54.283715; clon = -36.495695;
        selected_day = datenum('21-Jun-2018'); % choose an equinox day
    otherwise
%         clat = MPD.Meta.Location(1); clon = MPD.Meta.Location(2);
        clat = -53.786049; clon = -67.722423;
%         selected_day = datenum('21-Jun-2010'); % choose an equinox day
        selected_day = datenum('21-Jun-2020'); % choose an equinox day
end


inds = floor(tim) == selected_day;

% get local time:
localtime = mod(tim + (MPD.Meta.TimeZone/24),1) * 24;
localtime = localtime(inds);
xd = xdist(inds);
yd = ydist(inds);

xd(quadadd(xd,yd) > 500) = NaN;
yd(quadadd(xd,yd) > 500) = NaN;

%%%% plot by time of day:
csteps = 0:(24/clen):(24-(24/clen));

for i = 1:length(csteps)
    cinds = inrange(localtime,[csteps(i) csteps(i)+(24/clen)]);
    markerspec = {'marker','.','markersize',14,'color',cmap(i,:)};
    hold on; plot3(xd(cinds),yd(cinds),rand(size(xd(cinds))),'linest','none',markerspec{:})
end

axis square
% grid on;
% box on;

axx = gca;
axx.XColor = 'none';
axx.YColor = 'none';

xlim(pm(375))
ylim(pm(375))

%%%% COASTLINE:
[LatBox,~] = reckon(clat,clon,km2deg(abs(axx.YLim)),[0 180]);
[~,LonBox] = reckon(clat,clon,km2deg(abs(axx.XLim)),[90 270]);
LonBox = sort(LonBox); LatBox = sort(LatBox);
% LonBox = clon + pm(30); LatBox = clat + pm(20);

C = nph_draw_coastline([LonBox(1) LatBox(1) ; LonBox(2) LatBox(2)],0,0,'noplot','color','k');
for q = 1:length(C)
    if length(C(q).Lon) > 1
        [d,az] = distance(clat,clon,C(q).Lat,C(q).Lon);
        hold on; plot3(deg2km(d).*sind(az),deg2km(d).*cosd(az),10*ones(size(d)),'color',[1 1 1 1],'linewi',3);
        hold on; plot3(deg2km(d).*sind(az),deg2km(d).*cosd(az),10*ones(size(d)),'color',[0 0 0 1],'linewi',2);
    end
end



%%%% let's draw some range rings and lines:
gridlinespec = {'linewi',2,'color',[0.15 0.15 0.15 0.15]};
hold on; plot(1.1*axx.XLim,[0 0],gridlinespec{:});
hold on; plot([0 0],1.1*axx.YLim,gridlinespec{:});

ranges = [200 300 400];
azz = 0:0.01:(2*pi);
for r = ranges
    hold on; plot3(r*sin(azz),r*cos(azz),10*ones(size(azz)),gridlinespec{:});
    % label the range lines
    textaz = 300;
    hold on; text(r*sind(textaz),r*cosd(textaz),10,[num2str(r) 'km'],'fontsize',1*fs,'horizontalalignment','center','VerticalAlignment','middle')
end

%%%% labels
hold on; text(0,1.1*axx.YLim(2),'N','fontsize',fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','bottom');
hold on; text(1.1*axx.XLim(2),0,'E','fontsize',fs,'fontweight','bold','horizontalalignment','left','VerticalAlignment','middle');
hold on; text(0,1.1*axx.YLim(1),'S','fontsize',fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
hold on; text(1.1*axx.XLim(1),0,'W','fontsize',fs,'fontweight','bold','horizontalalignment','right','VerticalAlignment','middle');

%%%% KEP label
for n = (1*[-1 +1])
    hold on; text(n,-n,10,{[' ' upper(site)],'/'},'color','w','fontsize',fs,'fontweight','bold','horizontalalignment','left','VerticalAlignment','bottom');
    hold on; text(n,+n,10,{[' ' upper(site)],'/'},'color','w','fontsize',fs,'fontweight','bold','horizontalalignment','left','VerticalAlignment','bottom');
end
hold on; text(0,0,10,{[' ' upper(site)],'/'},'fontsize',fs,'fontweight','bold','horizontalalignment','left','VerticalAlignment','bottom');


set(gca,'layer','top','clipping','off')

% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Meteors against Altitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(rows,cols,axs{2})

hold on; H = histogram(alt,65:115,'visible','off');
hold on; barh(H.BinEdges(1:end-1)+0.5,H.Values./ndays,'barwidth',1,'facecolor',barcolor);
grid on;

switch site
    case 'kep'
        xlim([0 400])
        ylim([60 120])
        ytick(60:15:120)
        yminortick(60:5:120)
    otherwise
        xlim([0 1200])
        ylim([60 120])
        ytick(60:15:120)
        yminortick(60:5:120)
end

%%%% LABELS
xlabel('N_{met}')
ylabel('Altitude (km)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Meteors against Zenith Angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(rows,cols,axs{3})


hold on; H = histogram(zen,0:2.5:90,'visible','off');
hold on; bar(H.BinEdges(1:end-1),H.Values./ndays,'barwidth',1,'facecolor',barcolor);
grid on;

switch site
    case 'kep'
        xlim([0 90])
        ylim([0 400])
        xtick(0:15:90)
        ytick(0:100:500)
        yminortick(0:50:500)
        %%%% add zenith angle limits
        axx = gca;
        hold on; plot([15 15],axx.YLim,'linest','--','color',[0 0 0],'linewi',1.5);
        hold on; plot([65 65],axx.YLim,'linest','--','color',[0 0 0],'linewi',1.5);
        
    otherwise
        xlim([0 90])
        ylim([0 1400])
        xtick(0:15:90)
        ytick(0:200:2000)
        yminortick(0:100:2000)
end


%%%% LABELS
ylabel('N_{met}')
xlabel('Zenith Angle (deg)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Meteors against Ground Range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(rows,cols,axs{4})

hold on; H = histogram(quadadd(xdist,ydist),0:10:500,'visible','off');
hold on; bar(H.BinEdges(1:end-1)+0.5,H.Values./ndays,'barwidth',1,'facecolor',barcolor);
grid on;

switch site
    case 'kep'
        xlim([0 500])
        ylim([0 400])
        ytick(0:100:500)
        yminortick(0:50:500)
    otherwise
        xlim([0 500])
        ylim([0 2000])
        ytick(0:400:2000)
        yminortick(0:200:2000)
end


%%%% LABELS
ylabel('N_{met}')
xlabel('Ground Range (km)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Meteors against Time of day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(rows,cols,axs{5})

localtime = mod(tim + (MPD.Meta.TimeZone/24),1) * 24;

hold on; H = histogram(localtime,0:24,'visible','off');
binedges = H.BinEdges(1:end-1);
% hold on; bar(H.BinEdges(1:end-1)+0.5,H.Values./ndays,'barwidth',1,'facecolor',mcolor(6));

% colour by time:
for i = 1:length(csteps)
    cinds = inrange(binedges,[csteps(i) csteps(i)+(24/clen)]);
    hold on; bar(binedges(cinds)+0.5,H.Values(cinds)./ndays,'barwidth',1,'facecolor',cmap(i,:));
end

grid on;

switch site
    case 'kep'
        xlim([0 24])
        ylim([0 400])
        xtick(0:6:24)
        xminortick(0:3:24)
        ytick(0:100:500)
        yminortick(0:50:500)
    otherwise
        xlim([0 24])
        ylim([0 1000])
        xtick(0:6:24)
        xminortick(0:3:24)
        ytick(0:100:1000)
        yminortick(0:50:1000)
end


%%%% LABELS
ylabel('N_{met}')
xlabel('Local Time (hours)')


% % % % 
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %% %% Meteors against Azimuth
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % 
% % % % subplot(rows,cols,axs{3})
% % % % 
% % % % hold on; H1 = histogram(az,0:10:360,'visible','off');
% % % % hold on; H2 = histogram(az(~inrange(zen,[15 65])),0:10:360,'visible','off');
% % % % 
% % % % hold on; bar(H1.BinEdges(1:end-1)+5,H1.Values./ndays,'barwidth',1,'facecolor',rgbtrip(0.85));
% % % % hold on; bar(H1.BinEdges(1:end-1)+5,(H1.Values-H2.Values)./ndays,'barwidth',1,'facecolor',barcolor);
% % % % 
% % % % 
% % % % grid on;
% % % % box on;
% % % % 
% % % % 
% % % % axx = gca;
% % % % 
% % % % xlim([0 360])
% % % % ylim([0 250])
% % % % 
% % % % xtick(0:90:360)
% % % % xminortick(0:30:360)
% % % % 
% % % % setfont(fs);
% % % % set(gca,'linewi',1.5,'tickdir','out');
% % % % 
% % % % ylabel('N_{met}')
% % % % xlabel({'Azimuth','(degrees clockwise from North)'})
% % % % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Meteors against time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(rows,cols,axs{6})

hold on; H = histogram(tim,unique(floor(tim)),'visible','off');
hold on; bar(H.BinEdges(1:end-1)+0.5,(H.Values./1000),'barwidth',1,'facecolor',barcolor);
% hold on; stairs(H.BinEdges(1:end-1)+0.5,H.Values./ndays,'color','k');

grid on;

switch site
    case 'kep'
        ylim([0 15000]./1000)
        ytick([0 5000 10000 15000]./1000)
        yminortick((0:2500:15000)./1000)
    otherwise
        ylim([0 45000]./1000)
        ytick([0:10000:40000]./1000)
        yminortick((0:5000:50000)./1000)
end

% xlim([datenum(2016,01,01) datenum(2020,12,01)])

%%%% complicated xticks as month letters...
% Just show the Januarys as major ticks
axx = gca;
axx.XTick = datenum([years years(end)+1],01,01);
axx.XMinorTick = 'on';
axx.XAxis.MinorTickValues = datenum(min(years),1:((length(years)+1)*12),01);
datetick('x','m','keepticks','keeplimits')
axx.XTickLabel = {};

% % other months as ticks using text:
% xtix = datenum(years(1),1:(length(years)*12),01);
% for xt = xtix
%     if inrange(xt,xlim) && ~any(xt == axx.XTick)
%         mn = monthname(month(xt));
%         text(xt,-500./1000,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
%     end
% end
%
% % Years as lower bold numbers:
% xtix = datenum(years,07,01);
% for xt = xtix
%     if inrange(xt,xlim) && ~any(xt == axx.XTick)
%         yrstr = datestr(xt,'yyyy');
%         text(xt,-2000./1000,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
%     end
% end

xlim([datenum(min(years),01,01) datenum(max(years)+1,01,01)])
% % Just show the Januarys as major ticks
% axx.XTick = datenum([years years(end)+1],01,01);
% axx.XMinorTick = 'on';
% axx.XAxis.MinorTickValues = datenum(2016,1:((length(years)+1)*12),01);
% datetick('x','m','keepticks','keeplimits')
%             axx.XTickLabel = {};
switch site
    case 'kep'
        % other months as ticks using text:
        xtix = datenum(years(1),1:(length(years)*12),15);
        for xt = xtix
            if inrange(xt,xlim) && ~any(xt == axx.XTick)
                mn = monthname(month(xt));
                text(xt,-500./1000,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
            end
        end
end

% Years as lower bold numbers:
xtix = datenum(years,07,01);
for xt = xtix
    if inrange(xt,xlim) && ~any(xt == axx.XTick)
        yrstr = datestr(xt,'yyyy');
        text(xt,-2000./1000,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
    end
end


%%%% LABELS
ylabel('N_{met}  x10^{3}')
% ylabel('Meteor Counts  x10^{3}')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GLOBAL FORMATTING FOR ALL AXES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ax = 1:6
    
    subplot(rows,cols,axs{ax});
    
    axx = gca;
    
    setfont(fs);
    set(gca,'linewi',1.5,'tickdir','out');
    
    switch ax
        case 1
            hold on; nph_text([-0.025 0.915],['(' alphabet(ax) ')'],'fontsize',1.75*fs,'textborder','w');
        case {2 3 4 5}
            hold on; nph_text([ 0.025 0.865],['(' alphabet(ax) ')'],'fontsize',1.75*fs,'textborder','w');
        case 6
            hold on; nph_text([-0.03 0.865],['(' alphabet(ax) ')'],'fontsize',1.75*fs,'textborder','w');
    end
    
end








%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_HistogramsAndMap_2020-2022'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return





































return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(rows,cols,3)

hold on; contourf(X,Y,Z./ndays,31,'edgecolor','none'); shat;
colormap(gca,cbrew('nph_RdYlBuGrey',31));

axis square
grid on;
box on;
set(gca,'linewi',1.5,'layer','top')

xlim(minmax(X))
ylim(minmax(Y))

%%%% COASTLINE
switch site
    case 'kep'
        LonBox = [-47.5 -26.5];
        LatBox = [-59 -49];
    case 'rothera-sk'
        LonBox = [-80 -50];
        LatBox = [-70 50];
end
% and actually plot the coastline
C = nph_draw_coastline([LonBox(1) LatBox(1) ; LonBox(2) LatBox(2)],0,'noplot');
for i = 1:length(C)
    %     if length(C(i).Lon) > 10
    [dd,azz] = distance(MPD.Meta.Location(1),MPD.Meta.Location(2),C(i).Lat,C(i).Lon);
    xdd = deg2km(dd) .* cosd(azz);
    ydd = deg2km(dd) .* sind(azz);
    hold on; plot(xdd,ydd,'color','w','linewi',2.5);
    hold on; plot(xdd,ydd,'color','k','linewi',1.5);
    %     end
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subplot(rows,cols,[4 5 6])

axx = gca;

hold on; H = histogram(floor(tim),datenum(2016,01,01):datenum(2020,12,31),'visible','off');
hold on; bar(H.BinEdges(1:end-1)+0.5,H.Values,'barwidth',1,'facecolor',mcolor(6));
grid on;
box on;

xlim([datenum(2016,01,01) datenum(2020,02,01)])

yyaxis('right')
axx.YColor = 'k';

[yy,mm,~,~,~,~] = datevec(tim);
for yr = years
    inds = yy == yr;
    hold on; H = histogram(mm(yy),1:13,'visible','off');
    %     hold on; bar(H.BinEdges(1:end-1)+15,H.Values,'barwidth',1,'facecolor',mcolor(6));
    hold on; bar(datenum(yr,1:12,01),H.Values,'barwidth',1,'facecolor',mcolor(6));
end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now scroll though and get a rolling 30-day composite:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('Assembling an all-years composite for the specified years...')

days = ceil(linspace(1,365,120));

daywindow = 30;

std_z = 1.275; % FWHM for height gaussian
zrange = 76:1:105;

CompYear = struct;
CompYear.u = nan(length(zrange),length(days));
CompYear.v = nan(length(zrange),length(days));
CompYear.walt = nan(length(zrange),length(days));

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
    
    % select todays composite meteors:
    az_day = az(timerange);
    vhorz_day = vhorz(timerange);
    alt_day = alt(timerange);
    tim_day = tim(timerange);
    dn_day = dn(timerange);
    
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











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTY PLOT!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([0.75 0.8])

%-------------------------------------------------------
vert_gap = 0.085;   horz_gap = 0.08;
lower_marg = 0.075;   upper_marg = 0.05;
left_marg = 0.065;   right_marg = 0.075;

rows = 4; cols = 2;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 18;

axs = {...
    [1 2],...
    [3 4],...
    [5 7],...
    [6 8]};

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
    [-50:5:50],...
    [-20:2.5:20],...
    [-50:5:50],...
    [-20:2.5:20]};

clabs = {...
    'u','v','u','v'};

cbarticks = {...
    -50:25:50,...
    -20:10:20,...
    -50:25:50,...
    -20:10:20};

tits = {...
    'ZONAL',...
    'MERIDIONAL',...
    'ZONAL',...
    'MERIDIONAL'};


% also contour lines:
clines = {...
    [-50 -25 25 50],...
    [-20 -10 10 20],...
    [-50 -25 25 50],...
    [-20 -10 10 20]};


for ax = 1:4
    
    subplot(rows,cols,axs{ax})
    
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
    hold on; [C,h] = contour(X,Y,tbp,clines{ax},'edgecolor',rgbtrip(.25));
    clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
    % and bold zero line
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','w','linewi',4);
    clabel(C,h,'fontsize',0.8*fs,'labelspacing',200,'color','w','fontweight','bold');
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linewi',2);
    clabel(C,h,'fontsize',0.8*fs,'labelspacing',1000);
    
    
    switch ax
        case {1,2}
            %             axx.XTick = datenum(2016,1:(length(years)*12),01);
            axx.XTick = [datenum(2017,[1],01) datenum(2018,[1],01) datenum(2019,[1],01) datenum(2020,[1],01)];
            axx.XMinorTick = 'on';
            axx.XAxis.MinorTickValues = datenum(2016,1:(length(years)*12),01);
            datetick('x','mmm yyyy','keepticks','keeplimits')
            xlim([datenum(2016,02,15) datenum(2020,02,01)])
        case {3,4}
            axx.XTick = datenum(2019,1:13,01);
            datetick('x','m','keepticks','keeplimits')
            xlim(datenum(2019,[1 13],01))
    end
    
    % LIMITS
    ylim([77 103])
    axx.YTick = 70:5:110;
    axx.YMinorTick = 'on';
    switch ax
        case {1,2}
            axx.YAxis.MinorTickValues = 70:2.5:110;
        case {3,4}
            axx.YAxis.MinorTickValues = 70:110;
    end
    
    switch ax
        case {1,2,3}
            ylabel('Altitude (km)')
    end
    
    % line around the axis
    hold on; plot(axx.XLim([1 2 2]),axx.YLim([2 2 1]),'color','k','linewi',1.5)
    
    drawnow;
    
    % AXES LETTER
    switch ax
        case {1,2}
            hold on; nph_text([-0.03 1.115],['(' alphabet(ax) ')'],'fontsize',1.5*fs,'textborder','w');
        case {3,4}
            hold on; nph_text([-0.01 1.02],['(' alphabet(ax) ')'],'fontsize',1.5*fs,'textborder','w');
    end
    
    % COLORI
    %     cmap = cbrew('RdBu',length(clevs{ax}));
    %     cmap = cbrew('nph_RdYlBu',length(clevs{ax}));
    %     cmap = cbrew('nph_BuOr',length(clevs{ax}));
    %     cmap = cbrew('BrBG',length(clevs{ax}));
    %     cmap = flipud(cbrew('PiYG',length(clevs{ax})));
    %     cmap = cbrew('nph_Spectral',length(clevs{ax}));
    %     cmap = cbrew('nph_modspectral',length(clevs{ax}));
    cmap = cbrew('nph_modspectral',16);
    %     cmap = cbrew('nph_RdBuPastel',length(clevs{ax}));
    %     colormap(gca,cmap)
    
    % sort out double white in the middle:
    %     cmap = [cmap(1:floor(length(clevs{ax})/2),:) ; [1 1 1] ; cmap(floor(length(clevs{ax})/2)+1:end,:)];
    
    % saturate?
    cmap_hsv = rgb2hsv(cmap);
    cmap_hsv(:,2) = 1.2.*cmap_hsv(:,2); cmap_hsv(cmap_hsv > 1) = 1;
    cmap_sat = hsv2rgb(cmap_hsv);
    colormap(gca,cmap_sat)
    
    clim(clims{ax})
    
    % COLORBAR
    switch ax
        case {1,2}
            cbar = nph_colorbar;
            cbar.Ticks = cbarticks{ax};
            cbar.Position = cbar.Position .* [0.99 1 0.5 1];
            cbar.Label.String = ['\bf{' clabs{ax} '}'];
            cbar.Label.Rotation = 0;
            cbar.Label.VerticalAlignment = 'middle';
            %             cbar.Label.HorizontalAlignment = 'left';
            cbar.Title.String = '   ms^{-1}';
        case {3,4}
            cbar = nph_colorbar;
            cbar.Ticks = cbarticks{ax};
            cbar.Position = cbar.Position .* [1 1 1 0.6];
            cbar.Position = cbar.Position +  [0 0.2*cbar.Position(4) 0 0];
            cbar.Label.String = ['\bf{' clabs{ax} '}'];
            cbar.Label.Rotation = 0;
            cbar.Label.VerticalAlignment = 'middle';
            %             cbar.Label.HorizontalAlignment = 'left';
            cbar.Title.String = '   ms^{-1}';
            
            
    end
    
    switch ax
        case {1,2}
            %             text(datenum(2015,10,01),90,tits{ax},'fontsize',1*fs,'fontweight','bold','HorizontalAlignment','center','rotation',90)
            text(datenum(2016,4,15),107.25,tits{ax},'fontsize',1*fs,'fontweight','bold','HorizontalAlignment','left','rotation',0)
            
        case {3,4}
            text(datenum(2019,02,01),104.75,tits{ax},'fontsize',1*fs,'fontweight','bold','HorizontalAlignment','left')
    end
    
    
    setfont(fs)
    
    set(gca,'linewi',1.5,'layer','top','tickdir','out')
    
    
end



return




%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_CompYearMeanWinds'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return











































% new approach for this - use a monthly value based on using a COMPOSITE
% day, or rather composite time and height bins and then take the mean
% component.
%
%
% plots of MONTHLY MEAN COMPOSITE wind speed for time against height
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% SITE AND TIME ==========================================================

site = 'kep';

years = 2016:2019;


%% DIRECTORIES ============================================================

mpddirec = ['/Users/neil/data/MeteorRadar/' site '/matlab/'];
%     mpddirec = ['/Volumes/SDRed/data/MeteorRadar/' site '/matlab/'];

YWD = struct;


for yr = years
    
    yrstr = num2str(yr);
    
    hourrange = datenum(yr,01,01):(1/24):datenum(yr,12,31)+(23/24);
    
    %% LOAD MPD FILES =========================================================
    % LOAD ALL MPD FILES FOR THIS YEAR
    
    disp(['Loading all MPDs for ' yrstr '...'])
    
    mpdfiles = dir([mpddirec yrstr '*' site '_mpd.mat' ]);
    
    allMPD = struct;
    
    dayrange = datenum(yr,01,01):datenum(yr,12,31);
    
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
    [~,mmm,~,~,~,~] = datevec(tim);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% WIND FITTING SPECIFICATIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % HEIGHT AND TIME WINDOW STDs:
    std_z = 1.275; % gives approx 3km FWHM
    %     std_time = 0.85/24; % DAYS, gives approx 2hrs FWHM
    std_time = 0.85; % HOURS, gives approx 2hrs FWHM
    
    %     zrange = 76:1:105; % final height range for output matrix
    zrange = 72:1:108; % final height range for output matrix
    zlen = length(zrange);
    
    % NUMBER OF METEORS THRESHOLD FOR FIT:
    nmet_threshold = 20;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CALCULATE GAUSSIAN WEIGHTINGS AND WINDS IN EACH COMPOSITE TIME HEIGHT BIN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % To save a bit of time, pre-compute the gaussian weighting
    % functions for the height bins, they'll be the same each time:
    Gauss_z = struct;
    
    for z = 1:length(zrange)
        Gauss_z(z).vec = exp(- ((alt - zrange(z)).^2) ./ (2 * std_z^2));
    end
    
    
    % FOR EACH MONTH:
    for m = 1:12
        
        MWD.Data(m).u = nan(length(zrange),24);
        MWD.Data(m).v = nan(length(zrange),24);
        MWD.Data(m).nmets = nan(length(zrange),24);
        
        % display current month we're working on...
        disp(datestr(datenum(yr,m,01),'mmm yyyy'))
        
        % find inds for this month:
        monthinds = mmm == m;
        
        hourofdaymonth = hourofday(monthinds);
        altmonth = alt(monthinds);
        azmonth = az(monthinds);
        vhorzmonth = vhorz(monthinds);
        
        % for each height step
        for z = 1:length(zrange)
            
            % Gaussian HEIGHT weighting
            gauss_z = Gauss_z(z).vec(monthinds);
            %             gauss_z = exp(- ((alt - zrange(z)).^2) ./ (2 * std_z^2));
            
            % for each hour
            for h = 0:23
                
                % compute composite Gaussian TIME weighting
                gauss_time = exp(- ((h - hourofdaymonth).^2) ./ (2 * std_time^2));
                
                % combine the two weightings:
                w = gauss_time .* gauss_z;
                
                % only choose meteors within 2STDs:
                rng = w > 0.05;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% HOURLY WIND FIT
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                nmets = sum(rng);
                
                % Choose whether to do the sine fit:
                if nmets > nmet_threshold
                    [yfit,F] = nph_sinefit(azmonth(rng),vhorzmonth(rng),360,'weights',w(rng));
                else
                    F = [NaN NaN NaN];
                end
                
                % F(1) = cos part, F(2) = sin part, F(3) = mean.
                % so, if my definition of azimuth is correct, F(2) the sine
                % part should be u and F(1) cos part should be v.
                
                % Check for silly numbers:
                F(abs(F) > 200) = NaN;
                
                % Subscribe!
                MWD.Data(m).u(z,h+1) = F(2);
                MWD.Data(m).v(z,h+1) = F(1);
                
                MWD.Data(m).nmets(z,h+1) = nmets;
                
                
            end % next hour
            
        end % next height
        
        
        % Nice. Now to plot a year's mean winds, we need to assemble and take the
        % mean of each month.
        
        
        YWD.(['yr' yrstr]).u(:,m) = nanmean(MWD.Data(m).u,2);
        YWD.(['yr' yrstr]).v(:,m) = nanmean(MWD.Data(m).v,2);
        
        
    end % next MONTH
    
end % next year




return







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTY PLOT!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([0.6 1])

%-------------------------------------------------------
vert_gap = 0.07;   horz_gap = 0.025;
lower_marg = 0.075;   upper_marg = 0.04;
left_marg = 0.125;   right_marg = 0.075;

rows = 4; cols = 2;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 16;

% assemble all years together:

% winds = 'zonal'; % zonal | meridional
%
% tbp = [];
% switch winds
%     case 'zonal'
%         for ax = 1:length(years)
%             tbp = cat(2,tbp,YWD.(['yr' num2str(years(ax))]).u);
%         end
%     case 'meridional'
%         for ax = 1:length(years)
%             tbp = cat(2,tbp,YWD.(['yr' num2str(years(ax))]).v);
%         end
% end
% [X,Y] = meshgrid(datenum(years(1),1:(12*length(years)),01),zrange);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZONAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

winds = 'zonal';

tbp = [];
for ax = 1:length(years)
    tbp = cat(2,tbp,YWD.(['yr' num2str(years(ax))]).u);
end

[X,Y] = meshgrid(datenum(years(1),1:(12*length(years)),01),zrange);



for ax = 1:4
    
    subplot(rows,cols,ax)
    
    yr = years(ax);
    
    clims = [-45 45];
    clev = -45:5:45;
    
    tbp(tbp < -40) = -40;
    %     tbp = nph_fix2clims(tbp,clims);
    
    hold on; contourf(X,Y,tbp,clev,'edgecolor','none'); grid on;
    
    hold on; [C,h] = contour(X,Y,tbp,[-40:10:-10 10:10:40],'edgecolor',[0.25 0.25 0.25],'linest','-','linewi',1);
    
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linest','-','linewi',1.5);
    hold on; clabel(C,h,'fontsize',16);
    
    colormap(gca,cbrew('nph_BuOr',length(clev)+3))
    clim(clims)
    
    xlim(datenum(yr,[1 13],01))
    ylim([73 107])
    ytick(70:5:110)
    
    xtick(datenum(yr,1:13,01))
    datetick('x','mmm','keepticks','keeplimits')
    xlabel(num2str(yr),'fontweight','bold')
    
    set(gca,'tickdir','out','linewi',1.5)
    
    setfont(fs)
    
    % overlay box
    xl = get(gca,'xlim'); yl = get(gca,'ylim');
    hold on; plot(xl([1 2 2 1 1]),yl([1 1 2 2 1]),'color','k','linewi',1.5)
    
    switch ax
        case {2,4}
            cbar = nph_colorbar;
            cbar.Label.String = 'ms^-^1';
            set(gca,'yticklabel',{})
        otherwise
            ylabel('Altitude (km)')
    end
    
    if ax == 1, text(xl(1)-87.5,67.5,'Zonal Wind','fontsize',1.5*fs,'fontweight','bold','rotation',90,'HorizontalAlignment','center'); end
    
    nph_text([0.04 1.125],['(' alphabet(ax) ')'],'fontsize',1.5*fs,'textborder','w');
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MERIDIONAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

winds = 'meridional';

tbp = [];
for ax = 1:length(years)
    tbp = cat(2,tbp,YWD.(['yr' num2str(years(ax))]).v);
end

[X,Y] = meshgrid(datenum(years(1),1:(12*length(years)),01),zrange);



for ax = 5:8
    
    subplot(rows,cols,ax)
    
    yr = years(ax-4);
    
    clims = [-20 20];
    clev = -20:2.5:20;
    
    tbp(tbp < -20) = -20;
    %     tbp = nph_fix2clims(tbp,clims);
    
    hold on; contourf(X,Y,tbp,clev,'edgecolor','none'); grid on;
    
    hold on; [C,h] = contour(X,Y,tbp,[-20:5:-5 5:5:20],'edgecolor',[0.25 0.25 0.25],'linest','-','linewi',1);
    
    hold on; [C,h] = contour(X,Y,tbp,[0 0],'edgecolor','k','linest','-','linewi',1.5);
    hold on; clabel(C,h,'fontsize',16);
    
    colormap(gca,cbrew('nph_BuOr',length(clev)+3))
    clim(clims)
    
    xlim(datenum(yr,[1 13],01))
    ylim([73 107])
    ytick(70:5:110)
    
    xtick(datenum(yr,1:13,01))
    datetick('x','mmm','keepticks','keeplimits')
    xlabel(num2str(yr),'fontweight','bold')
    
    set(gca,'tickdir','out','linewi',1.5)
    
    setfont(fs)
    
    % overlay box
    xl = get(gca,'xlim'); yl = get(gca,'ylim');
    hold on; plot(xl([1 2 2 1 1]),yl([1 1 2 2 1]),'color','k','linewi',1.5)
    
    switch ax
        case {6,8}
            cbar = nph_colorbar;
            cbar.Label.String = 'ms^-^1';
            cbar.Ticks = -20:10:20;
            set(gca,'yticklabel',{})
        otherwise
            ylabel('Altitude (km)')
    end
    
    if ax == 5, text(xl(1)-87.5,67.5,'Meridional Wind','fontsize',1.5*fs,'fontweight','bold','rotation',90,'HorizontalAlignment','center'); end
    
    nph_text([0.04 1.125],['(' alphabet(ax) ')'],'fontsize',1.5*fs,'textborder','w');
    
    
end


return

%% EXPORT? ================================================================

savename = ['~/Desktop/' upper(site) '_MonthlyCompWinds'];
disp(['Exporting to "' savename '"...'])
nph_saveas(gcf,savename,'png')

disp('Done.')


return
































%% LOAD HWD ===============================================================

disp('Loading HWD...')

uu = []; vv = []; tt = [];

for yr = years
    
    yrstr = num2str(yr);
    
    % Draw number of Meteors versus time of day for KEP
    
    load(['/Users/neil/data/MeteorRadar/' site '/matlab/hwd/' yrstr '_' site '_hwd.mat'])
    %     load(['/Volumes/SDCard/data/MeteorRadar/' site '/matlab/hwd/' yrstr '_' site '_hwd.mat'])
    
    uu = cat(3,uu,MWD.Data.u);
    vv = cat(3,vv,MWD.Data.v);
    tt = cat(2,tt,MWD.Data.Day);
    
    %     nan_inds = all(isnan(uu),1);
    %
    %     uu_sm = movmean(uu,24*30,2,'omitnan');
    %
    %     uu_sm_nan = uu_sm;
    %     uu_sm_nan(:,nan_inds) = NaN;
    
end

% Instead of reshape, simply take the daily mean wind:
uu = squeeze(nanmean(uu,2));
vv = squeeze(nanmean(vv,2));

% nans...
nandays = isnan(movmean(uu(4,:),5,2,'omitnan'));

% smooth...
uu = movmean(uu,31,2,'omitnan');
vv = movmean(vv,31,2,'omitnan');

% replace nans for spline interpolation:
uu(isnan(uu)) = 0;
vv(isnan(vv)) = 0;

% interpolate...
z =  [78 MWD.Data.Alt 100];
zi = 78:1:100;
uui = interp1(z,uu([1 1:6 6],:),zi,'spline');
vvi = interp1(z,vv([1 1:6 6],:),zi,'spline');

% replace nans...
uu(:,nandays) = NaN;
vv(:,nandays) = NaN;

% replace nans again to make sure...
uui(:,nandays) = NaN;
vvi(:,nandays) = NaN;

% and now trim to timerange...
uui(:,~inrange(tt,min(timerange),max(timerange))) = NaN;
vvi(:,~inrange(tt,min(timerange),max(timerange))) = NaN;





%% LOAD MLS ===============================================================

disp('Loading MLS...')

mlsdirec = ['/Users/neil/data/MLS/matlab/'];

years = unique(year(timerange));

mlslat = []; mlslon = []; mlstime = []; mlstemp = [];

for yr = years
    
    yrstr = num2str(yr);
    
    % LOAD MLS:
    load([mlsdirec yrstr '_mls_t.mat'])
    
    % First, regionally subset:
    r = 500; % km from radar
    [d,~] = distance(MWD.Meta.Location(1),MWD.Meta.Location(2),MLS.Data.Latitude,MLS.Data.Longitude);
    
    distinds = find(deg2km(d) < r);
    
    % now temporally subset:
    timeinds = find(inrange(MLS.Data.Time,timerange(1),timerange(end)));
    
    % combine:
    inds = intersect(distinds,timeinds);
    
    mlslat = cat(1,mlslat,MLS.Data.Latitude(inds));
    mlslon = cat(1,mlslon,MLS.Data.Longitude(inds));
    mlstime = cat(1,mlstime,MLS.Data.Time(inds));
    mlstemp = cat(2,mlstemp,MLS.Data.Temperature(:,inds));
    
    mlsalt = MLS.Data.Altitude;
    mlspres = MLS.Data.Pressure;
    
    
    % %     %% build an interpolant for evaulating radar height:
    % %     % EDIT: Don't worry about this any more. The centre height doesn't seem
    % %     % to vary much around 89km, and when it does it's less that 2.5km,
    % %     % which when MLS's vertical resolution here is around 5km won't matter
    % %     % much at all.
    % %     Fcen = griddedInterpolant(HWD.Data.Day,HWD.Data.VertMetDist.Center,'linear','none');
    % %     ceni = Fcen(time);
    % %
    % %     Ffwhm = griddedInterpolant(HWD.Data.Day,HWD.Data.VertMetDist.FWHM,'linear','none');
    % %     fwhmi = Ffwhm(time);
    % %
    % %     alti = ceni; % 89.5km
    % %     Fmls = griddedInterpolant({alt,time},double(temp),'linear','none');
    % %     mlstempi = Fmls(ceni,time);
    % %
    % %     Store.(['Time' yrstr]) = time;
    % %     Store.(['Center' yrstr]) = ceni;
    % %     Store.(['FWHM' yrstr]) = fwhmi;
    % %     Store.(['MLSTemp' yrstr]) = mlstempi;
    
end

%% Assemble MLS temp into a matrix:

F = griddedInterpolant({mlsalt,mlstime},mlstemp,'linear','none');

mlstempi = F({zi,tt});

mlstempi = movmean(mlstempi,31,2,'omitnan');


return




%% plotting ===============================================================

figure; hold all; whitefig; figpos([1 1])

%-------------------------------------------------------
vert_gap = 0.075;       horz_gap = 0.01;
lower_marg = 0.075;   upper_marg = 0.04;
left_marg = 0.075;   right_marg = 0.1;

rows = 3; cols = 1;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 16;
nclevels = 16;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tbp = {};
tbp{1} = uui;
tbp{2} = vvi;
tbp{3} = mlstempi;

clims{1} = [-50 50];
clims{2} = [-20 20];
clims{3} = [150 210];

clevels{1} = [-50:5:50];
clevels{2} = [-20:2.5:20];
clevels{3} = [140:5:210];

clines{1} = [-40 -20 20 40];
clines{2} = [-20 -10 10 20];
clines{3} = [140:10:210];

clines_lighter{1} = [-30 -10 10 30];
clines_lighter{2} = [-15 -5 5 15];
clines_lighter{3} = [];

cmaps{1} = cbrew('nph_BuOr',length(clevels{1}));
cmaps{2} = cbrew('nph_BuOr',length(clevels{2}));
cmaps{3} = cbrew('RdBu',length(clevels{3}));

labels = {
    '\bf{u} \rm{(ms^{-1})}'
    '\bf{v} \rm{(ms^{-1})}'
    '\bf{T} \rm{(K)}'};

for ax = 1:3
    
    Axes(ax).ax = subplot(rows,cols,ax);
    axx = gca;
    
    % set limits:
    tbp{ax}(tbp{ax} > clims{ax}(2)) = clims{ax}(2);
    tbp{ax}(tbp{ax} < clims{ax}(1)) = clims{ax}(1);
    
    % a quick smooth:
    tbp{ax} = movmean(tbp{ax},5,2,'omitnan');
    
    % contourf
    hold on; contourf(tt,zi,tbp{ax},clevels{ax},'edgecolor','none'); grid on;
    
    % some selected contour lines
    hold on; [C,h] = contour(tt,zi,tbp{ax},clines{ax},'color',[.25 .25 .25],'linewi',2); grid on;
    clabel(C,h,'fontsize',fs,'labelspacing',2000,'color',[.25 .25 .25])
    hold on; contour(tt,zi,tbp{ax},clines_lighter{ax},'color',[.25 .25 .25],'linewi',1); grid on;
    
    % zero wind line:
    hold on; [C,h] = contour(tt,zi,tbp{ax},[0 0],'color','k','linewi',4); grid on;
    clabel(C,h,'fontsize',fs,'labelspacing',2000)
    
    % ticks and limits
    ylim([79 99])
    ytick(80:2:100)
    axx.YAxis.MinorTick = 'on';
    axx.YAxis.MinorTickValues = 75:100;
    xlim([datenum(min(years),month(min(timerange)),01) datenum(max(years),month(max(timerange))+1,01)])
    xtick(datenum(min(years),1:(12*(diff(minmax(years))+1)+1),01))
    datetick('x','m','keepticks','keeplimits')
    set(gca,'tickdir','out','linewi',1.5)
    
    % colors
    colormap(gca,cmaps{ax})
    clim(clims{ax})
    cbar = nph_colorbar;
    cbar.Position = cbar.Position .* [0.99 1 0.75 1];
    %     cbar.Label.String = ['\bf{' labels{ax} '} \rm{(ms^{-1})}'];
    cbar.Label.String = labels{ax};
    cbar.Ticks = clims{ax}(1):10:clims{ax}(2);
    %     cbar.Title.String = ['\bf{' labels{ax} '}'];
    %     cbar.Title.String = {'Zonal','Wind'};
    
    % font and labels
    setfont(fs)
    ylabel('Altitude (km)')
    
    
    % YEAR LABELS
    ypos = 75;
    for yr = years
        %     text(datenum(2016,06,15),ypos,'2016','fontsize',1*fs,'HorizontalAlignment','center','VerticalAlignment','middle')
        %     text(datenum(2017,06,15),ypos,'2017','fontsize',1*fs,'HorizontalAlignment','center','VerticalAlignment','middle')
        %     text(datenum(2018,06,15),ypos,'2018','fontsize',1*fs,'HorizontalAlignment','center','VerticalAlignment','middle')
        text(datenum(yr,12,31),ypos,'|','fontsize',0.75*fs,'HorizontalAlignment','center','VerticalAlignment','middle')
        text(datenum(yr,06,15),ypos,num2str(yr),'fontsize',fs,'HorizontalAlignment','center','VerticalAlignment','middle')
    end
    if ax == 3
        text(datenum(years(1),01,01),ypos,'|','fontsize',0.75*fs,'HorizontalAlignment','center','VerticalAlignment','middle')
    end
    
    
    % AXES LETTER
    text(min(timerange),99.5,['(' alphabet(ax) ')'],'fontsize',1.75*fs,'VerticalAlignment','bottom','HorizontalAlignment','left');
    
    % TITLE
    if ax == 1
        title(upper(site))
    end
    
end




return




%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_MeanWinds_' num2str(years(1)) '-' num2str(years(end)) '_plus_mls'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return



























































