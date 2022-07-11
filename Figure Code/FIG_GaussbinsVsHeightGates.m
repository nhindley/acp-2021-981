

% Comparing the new Gaussian weighting Method and the old height gates method


%% LOAD HWD

site = 'kep';

yr = 2017; yrstr = num2str(yr);

dayrange = datenum(yr,04,01):datenum(yr,06,01);

% Load HWD
gb = load(['/Users/neil/data/MeteorRadar/' site '/matlab/hwd/' yrstr '_' site '_hwd_gaussbins.mat']);
hg = load(['/Users/neil/data/MeteorRadar/' site '/matlab/hwd/' yrstr '_' site '_hwd.mat']);


% yyyymmdd = datestr(dy,'yyyymmdd');
% direc = '/Users/neil/data/MeteorRadar/kep/matlab/hwd';
% load([direc yyyymmdd '_' site '_mpd.mat'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; figpos([0.8 0.55]); whitefig;

%-------------------------------------------------------
vert_gap = 0.085;           horz_gap = 0.05;
lower_marg = 0.125;         upper_marg = 0.075;
left_marg = 0.08;           right_marg = 0.08;

rows = 2; cols = 2;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);
%--------------------------------------------------------

fs = 18;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORMATTING

uv = {'u','u','v','v'};

for ax = 1:4
    
    subplot(rows,cols,ax)
%     cmap = cbrew('RdBu',31);
    cmap = nph_saturate(cbrew('nph_modspectral',31),1.2);
    colormap(gca,cmap);
%     cmap = colormap(gca,[flipud(cbrew('PuBu',11)) ; [1 1 1] ; cbrew('YlOrRd',11)]);
    
    xlim(datenum(yr,05,[0.75 6.25])+2)
%     xlim(datenum(yr,05,[1 12]))
    ylim([75 105])
    
    axx = gca;
    axx.YMinorTick = 'on';
    axx.XMinorTick = 'on';
    axx.XAxis.MinorTickValues = dayrange(1):0.1:dayrange(end);
    
    clim([-100 100])
    
    set(gca,'xtick',dayrange(1):1:dayrange(end),'ytick',70:5:110)
    datetick('x','dd mmm','keepticks','keeplimits')
    
    set(gca,'linewi',1.5,'tickdir','out','layer','top','color',rgbtrip(1))
    
    % outer box
    xl = get(gca,'xlim'); yl = get(gca,'ylim'); 
    hold on; plot(xl([1 2 2 1 1]),yl([1 1 2 2 1]),'linewi',1,'color','k');
%     hold on; plot(xl([2 2 1]),yl([1 2 2]),'linewi',1.5,'color','k');

    % hatched area below:
    h = patch(axx.XLim([1 2 2 1]),axx.YLim([1 1 2 2]),'red');
    hh = hatchfill(h, 'cross', 45, 8);
    set(hh,'color',[0 0 0 0.15],'linewi',1)


%     for t = dayrange(1):0.1:dayrange(end)
%         hold on; plot([t t],axx.YLim,'color',[0 0 0 0.25],'linewi',1);
%     end
    
    switch ax
        case {2,4}
            cbar = nph_colorbar;
            cbar.Label.String = ['\bf{' uv{ax} '} \rm{(ms^-^1)}'];
            cbar.FontSize = fs;
            % cbar.Ticks = -80:20:80;
            cbar.Ticks = -90:30:90;
            cbar.Units = 'pixels';
            cbar.Position(1) = 1 * cbar.Position(1);
        case {1,3}
            ylabel('Altitude (km)')
    end
    
    if ax == 1, t1 = title('HEIGHT GATES METHOD','VerticalAlignment','bottom','fontweight','bold'); end
    if ax == 2, t2 = title('GAUSSIAN-WEIGHTING METHOD','VerticalAlignment','bottom','fontweight','bold'); end
    if ax == 1, text(axx.XLim(1)-0.85,90,'Zonal Wind','fontsize',1.25*fs,'fontweight','bold','Rotation',90,'HorizontalAlignment','center'); end
    if ax == 3, text(axx.XLim(1)-0.85,90,'Meridional Wind','fontsize',1.25*fs,'fontweight','bold','Rotation',90,'HorizontalAlignment','center'); end
    
    
    % label
    nph_text([-0.005  0.86],['(' alphabet(ax) ')'],'fontsize',1.75*fs,'textborder','w');
    
    setfont(fs)
    
% sort out local time
switch site
    case 'kep'
        localtimeshift = (-2 / 24);
%         localtimeshift = 0;
        if ax == 3 || ax == 4, xlabel('Local Time, 2017','fontweight','bold'); end
    otherwise
        localtimeshift = 0;
        if ax == 3 || ax == 4, xlabel('Local Time, 2017','fontweight','bold'); end
end

end

% return


%%% HEIGHT GATES METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(rows,cols,1)

hg.timeinds = inrange(hg.HWD.Data.Time,min(dayrange),max(dayrange));

for i = 1:6
    nanmask = ~isnan(hg.HWD.Data.u([i i],hg.timeinds));
    [X,Y] = meshgrid(hg.HWD.Data.Time(hg.timeinds)+localtimeshift,hg.HWD.Data.HeightGates(i,1:2));
    hold on; pcolor(X,Y,hg.HWD.Data.u([i i],hg.timeinds)); shat;
end

% Height gate lines:
yyaxis right
set(gca,'ycolor','k','linewi',1.5,'color','none','ylim',axx.YLim);
set(gca,'ytick',unique(hg.HWD.Data.HeightGates(:)),'yticklabel',{})
% for z = unique(hg.HWD.Data.HeightGates(:))
%     hold on; plot(dayrange,z.*ones(size(dayrange)),'color',[0 0 0 0.5],'linewi',1,'linest','-');
% end

%--------------------------------------------------------------------------

subplot(rows,cols,3)

hg.timeinds = inrange(hg.HWD.Data.Time,min(dayrange),max(dayrange));

for i = 1:6
    nanmask = ~isnan(hg.HWD.Data.u([i i],hg.timeinds));
    [X,Y] = meshgrid(hg.HWD.Data.Time(hg.timeinds)+localtimeshift,hg.HWD.Data.HeightGates(i,1:2));
    hold on; pcolor(X,Y,hg.HWD.Data.v([i i],hg.timeinds)); shat;
end

% Height gate lines:
yyaxis right
set(gca,'ycolor','k','linewi',1.5,'color','none','ylim',axx.YLim);
set(gca,'ytick',unique(hg.HWD.Data.HeightGates(:)),'yticklabel',{})
% for z = unique(hg.HWD.Data.HeightGates(:))
%     hold on; plot(dayrange,z.*ones(size(dayrange)),'color',[0 0 0 0.5],'linewi',1,'linest','-');
% end


%%% GAUSSIAN WEIGHTING METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(rows,cols,2)

gb.timeinds = inrange(gb.HWD.Data.Time,min(dayrange),max(dayrange));

nanmask = ~isnan(gb.HWD.Data.u(:,gb.timeinds));
% hold on; imagesc(gb.HWD.Data.Time(gb.timeinds),gb.HWD.Data.Alt,gb.HWD.Data.u(:,gb.timeinds),'AlphaData',nanmask); ydir;
hold on; pcolor(gb.HWD.Data.wtime(:,gb.timeinds)+localtimeshift,gb.HWD.Data.walt(:,gb.timeinds),gb.HWD.Data.u(:,gb.timeinds)); shat;


%--------------------------------------------------------------------------

subplot(rows,cols,4)

gb.timeinds = inrange(gb.HWD.Data.Time,min(dayrange),max(dayrange));

nanmask = ~isnan(gb.HWD.Data.u(:,gb.timeinds));
% hold on; imagesc(gb.HWD.Data.Time(gb.timeinds),gb.HWD.Data.Alt,gb.HWD.Data.v(:,gb.timeinds),'AlphaData',nanmask); ydir;
hold on; pcolor(gb.HWD.Data.wtime(:,gb.timeinds)+localtimeshift,gb.HWD.Data.walt(:,gb.timeinds),gb.HWD.Data.v(:,gb.timeinds)); shat;




% outer box
for ax = 1:4
    subplot(rows,cols,ax)
    xl = get(gca,'xlim'); yl = get(gca,'ylim'); 
    hold on; plot(xl([2 2 1]),yl([1 2 2]),'linewi',1.5,'color','k');
end



return

%% EXPORT? ================================================================

savename = ['~/Desktop/' upper(site) '_GaussbinsVsHeightGates'];
disp(['Exporting to "' savename '"...'])
nph_saveas(gcf,savename,'png')

disp('Done.')


return



%% FIGURE SHOWING THE ALT/TIME DIFFERENCE OF USING THE GAUSSIAN WEIGHTING METHOD:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; figpos([0.8 0.55]); whitefig;

%-------------------------------------------------------
vert_gap = 0.085;           horz_gap = 0.05;
lower_marg = 0.125;         upper_marg = 0.075;
left_marg = 0.08;           right_marg = 0.08;

rows = 2; cols = 2;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);
%--------------------------------------------------------


fs = 18;

[Z,T] = ndgrid(HWD.Data.Alt,(HWD.Data.Time));

TBP = {...
    24*60*(HWD.Data.wtime-T),...
    HWD.Data.walt-Z,...
    HWD.Data.walt-Z,...
    24*60*(HWD.Data.wtime-T)};

clims = {...
    [-30 30],...
    [-1.5 1.5],...
    [-1.5 1.5],...
    [-30 30]};

cbarticks = {...
    [-30:10:30],...
    [-1.5:0.5:1.5],...
    [-1.5:0.5:1.5],...
    [-30:10:30]};

clabels = {...
    'minutes',...
    'km',...
    'km',...
    'minutes'};


for ax = [1 2 4]
    
    subplot(rows,cols,ax)
%     cmap = cbrew('RdBu',31);
    cmap = nph_saturate(cbrew('nph_modspectral',31),1.2);
    colormap(gca,cmap);
%     cmap = colormap(gca,[flipud(cbrew('PuBu',11)) ; [1 1 1] ; cbrew('YlOrRd',11)]);
    
    xlim(datenum(yr,05,[0.75 6.25])+2)
%     xlim(datenum(yr,05,[1 12]))
    ylim([75 105])
    
    axx = gca;
    axx.YMinorTick = 'on';
    axx.XMinorTick = 'on';
    axx.XAxis.MinorTickValues = dayrange(1):0.1:dayrange(end);
    
    clim(clims{ax})
    
    % colorbar
    switch ax
        case 2
            cbar = nph_colorbar('westoutside');
            cbar.Label.String = clabels{ax};
            cbar.FontSize = fs;
            cbar.TickDirection = 'out';
            cbar.LineWidth = 1.5;
            % cbar.Ticks = -80:20:80;
            cbar.Ticks = cbarticks{ax};
            %             cbar.Units = 'pixels';
            %             cbar.Position(1) = 1 * cbar.Position(1);
            cbar1 = cbar;
            drawnow;
        case 4
            
            cbar = colorbar('westoutside');
            cbar.Label.String = clabels{ax};
            cbar.FontSize = fs;
            cbar.TickDirection = 'out';
            cbar.LineWidth = 1.5;
            cbar.Ticks = cbarticks{ax};
            
            cbar.Position(3:4) = cbar1.Position(3:4);
            
            set(gca,'xcolor','none','ycolor','none')
            continue
            
        case {1,3}
            ylabel('Altitude (km)')
    end
    
    
    set(gca,'xtick',dayrange(1):1:dayrange(end),'ytick',70:5:110)
    datetick('x','dd mmm','keepticks','keeplimits')
    
    set(gca,'linewi',1.5,'tickdir','out','layer','top','color',rgbtrip(1))
    
    % hatched area below:
    h = patch(axx.XLim([1 2 2 1]),axx.YLim([1 1 2 2]),'red');
    hh = hatchfill(h, 'cross', 45, 8);
    set(hh,'color',[0 0 0 0.15],'linewi',1)

    %%%% PLOT DATA
    hold on; pcolor(HWD.Data.wtime,HWD.Data.walt,TBP{ax}); shat
    
    % outer box
    xl = get(gca,'xlim'); yl = get(gca,'ylim');
    hold on; plot(xl([1 2 2 1 1]),yl([1 1 2 2 1]),'linewi',1.5,'color','k');
    %     hold on; plot(xl([2 2 1]),yl([1 2 2]),'linewi',1.5,'color','k');

    
%     for t = dayrange(1):0.1:dayrange(end)
%         hold on; plot([t t],axx.YLim,'color',[0 0 0 0.25],'linewi',1);
%     end
    
    
    if ax == 1, t1 = title('TIME DIFFERENCE','VerticalAlignment','bottom','fontweight','bold'); end
    if ax == 2, t2 = title('HEIGHT DIFFERENCE','VerticalAlignment','bottom','fontweight','bold'); end
%     if ax == 1, text(axx.XLim(1)-0.85,90,'Zonal Wind','fontsize',1.25*fs,'fontweight','bold','Rotation',90,'HorizontalAlignment','center'); end
%     if ax == 3, text(axx.XLim(1)-0.85,90,'Meridional Wind','fontsize',1.25*fs,'fontweight','bold','Rotation',90,'HorizontalAlignment','center'); end
    
    
    % label
    nph_text([-0.005  0.86],['(' alphabet(ax+4) ')'],'fontsize',1.75*fs,'textborder','w');
    
    setfont(fs)
    
% sort out local time
switch site
    case 'kep'
        localtimeshift = (-2 / 24);
%         localtimeshift = 0;
        if ax == 3 || ax == 4, xlabel('Local Time, 2017','fontweight','bold'); end
    otherwise
        localtimeshift = 0;
        if ax == 3 || ax == 4, xlabel('Local Time, 2017','fontweight','bold'); end
end


end

return


%% EXPORT? ================================================================

savename = ['~/Desktop/' upper(site) '_GaussbinsVsHeightGates_Difference'];
disp(['Exporting to "' savename '"...'])
nph_saveas(gcf,savename,'png')

disp('Done.')


return
















%% METEOR DOT MAP

subplot(rows,cols,[2:6 8:12])

for i = 1:1:nmets
    hold on; plot(tim1(i),alt1(i),'markerfacecolor',cmapi1({ww1(i),1:3}),'marker','o','linest','none','markersize',meteormarkersize,'markeredgecolor',[0.75 0.75 0.75]); grid on;
    hold on; plot(tim2(i),alt2(i),'markerfacecolor',cmapi2({ww2(i),1:3}),'marker','o','linest','none','markersize',meteormarkersize,'markeredgecolor',[0.75 0.75 0.75]); grid on;
end

ylim([73 107])
xlim([dy dy+1])

xl = get(gca,'xlim');
yl = get(gca,'ylim');
hold on; plot(xl([1 1 2 2 1]),yl([1 2 2 1 1]),'linewi',1.5,'color','k');

set(gca,'xtick',dy + (0:2/24:1),'ytick',70:5:120,'layer','top','linewi',1.5)
datetick('x','HH:MM','keepticks','keeplimits')

axx = gca;
axx.XMinorTick = 'on';
axx.XAxis.MinorTickValues = dy + (0:1/24:1);
axx.YMinorTick = 'on';
axx.YAxis.MinorTickValues = 70:2.5:120;

xlabel('Time Of Day')

% over plot bin edges:
hold on; plot((dy+(hr1/24)) + ([-1 +1 +1 -1 -1]./24),z1 + [-1.5 -1.5 +1.5 +1.5 -1.5],dashedlinespec1{:})
hold on; plot((dy+(hr2/24)) + ([-1 +1 +1 -1 -1]./24),z2 + [-1.5 -1.5 +1.5 +1.5 -1.5],dashedlinespec2{:})

% label
nph_text([0.04 0.92],'(a)','fontsize',1.8*fs,'textborder','w');

% legend
ii = fix(nmets ./ 10);
hold on; p = plot(tim2(ii),alt2(ii),'markerfacecolor','w','marker','o','linest','none','markersize',10,'markeredgecolor',[0.75 0.75 0.75]); grid on;
L = legend(p,'Meteor Detection','fontsize',fs);

% colorbar
clim([0 1])
colormap(gca,cmap1);
cbar1 = nph_colorbar;
cbar1.Label.String = 'Weighting';
cbar1.FontSize = fs;
cbar1.Ticks = 0:0.2:1;
cbar1.Units = 'pixels';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HEIGHT GAUSSIAN

subplot(rows,cols,[1 7])

hold on; plot(gauss_z1,MPD.Data.Alt,'linest','none','marker','.','markersize',5,'color',mcolor('blue')); grid on;
hold on; plot(gauss_z2,MPD.Data.Alt,'linest','none','marker','.','markersize',5,'color',mcolor('red')); grid on;

ylim(yl)
xlim([0 1.1])

set(gca,'xtick',0:0.2:1,'ytick',70:5:120,'layer','top','linewi',1.5)

% overlay bin edges:
hold on; plot([0 1.1],z1 + [-1.5 -1.5],dashedlinespec1{:},'linewi',1);
hold on; plot([0 1.1],z1 + [+1.5 +1.5],dashedlinespec1{:},'linewi',1);
hold on; plot([0 1.1],z2 + [-1.5 -1.5],dashedlinespec2{:},'linewi',1);
hold on; plot([0 1.1],z2 + [+1.5 +1.5],dashedlinespec2{:},'linewi',1);

xlabel('Weighting')
ylabel({'Altitude','(km)'})

% label
nph_text([0.8 0.1],'(b)','fontsize',1.8*fs,'textborder','w');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TIME GAUSSIAN

subplot(rows,cols,[14:18])

hold on; plot(MPD.Data.Time,gauss_time1,'linest','none','marker','.','markersize',5,'color',mcolor('blue')); grid on;
hold on; plot(MPD.Data.Time,gauss_time2,'linest','none','marker','.','markersize',5,'color',mcolor('red')); grid on;

ylim([0 1.1])
xlim(xl)

set(gca,'ytick',0:0.2:1,'xtick',dy + (0:2/24:1),'layer','top','linewi',1.5)
datetick('x','HH:MM','keepticks','keeplimits')

% overlay bin edges:
hold on; plot(dy+(hr1/24) + ([-1 -1]./24),[0 1.1],dashedlinespec1{:},'linewi',1);
hold on; plot(dy+(hr1/24) + ([+1 +1]./24),[0 1.1],dashedlinespec1{:},'linewi',1);
hold on; plot(dy+(hr2/24) + ([-1 -1]./24),[0 1.1],dashedlinespec2{:},'linewi',1);
hold on; plot(dy+(hr2/24) + ([+1 +1]./24),[0 1.1],dashedlinespec2{:},'linewi',1);

ylabel({'Weighting',' '})
xlabel('Time Of Day')

% label
nph_text([0.04 1.15],'(c)','fontsize',1.8*fs,'textborder','w');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SINE FITTING

subplot(rows,cols,[20:24 26:30])

for i = 1:1:nmets
    hold on; plot(az1(i),vhorz1(i),'markerfacecolor',cmapi1({ww1(i),1:3}),'marker','o','linest','none','markersize',meteormarkersize,'markeredgecolor',[0.75 0.75 0.75]); grid on;
    hold on; plot(az2(i),vhorz2(i),'markerfacecolor',cmapi2({ww2(i),1:3}),'marker','o','linest','none','markersize',meteormarkersize,'markeredgecolor',[0.75 0.75 0.75]); grid on;
end

ylim([-200 200])
xlim([0 360])

xl = get(gca,'xlim');
yl = get(gca,'ylim');
hold on; plot(xl([1 1 2 2 1]),yl([1 2 2 1 1]),'linewi',1.5,'color','k');

set(gca,'xtick',0:60:360,'ytick',-200:50:200,'layer','top','linewi',1.5)


% Over plot fitted sine waves
inds1 = w1 > 0.1;
inds2 = w2 > 0.1;

[~,F1] = nph_sinefit(az1(inds1),vhorz1(inds1),360,'weights',w1(inds1));
[~,F2] = nph_sinefit(az2(inds2),vhorz2(inds2),360,'weights',w2(inds2));

x = 0:360;
y1 = F1(1).*cosd(x) + F1(2).*sind(x) + F1(3);
y2 = F2(1).*cosd(x) + F2(2).*sind(x) + F2(3);

hold on; plot(x,y1,dashedlinespec1{:});
hold on; plot(x,y2,dashedlinespec2{:});

ylabel({'Horizontal Velocity','(m/s)'})
xlabel('Azimuth')

% label
nph_text([0.04 0.92],'(d)','fontsize',1.8*fs,'textborder','w');

% legend
hold on; p = plot(az2(ii),vhorz2(ii),'markerfacecolor','w','marker','o','linest','none','markersize',9,'markeredgecolor',[0.75 0.75 0.75]); grid on;
L = legend(p,'Meteor Detection','fontsize',fs);

% colorbar (move it to panel (a))
clim([0 1])
colormap(gca,cmap2);
cbar2 = nph_colorbar;
cbar2.Units = 'pixels';
cbar2.Position = get(cbar1,'position');
cbar2.Position = cbar2.Position .* [1 1 0.5 1];
cbar2.Ticks = [];


setfont(fs)



return

%% EXPORT? ================================================================

savename = ['~/Desktop/' yyyymmdd '_' upper(site) '_GaussianBinningMethod'];
disp(['Exporting to "' savename '"...'])
nph_saveas(gcf,savename,'png')

disp('Done.')


return
















































