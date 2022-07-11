
%%%% 4x1 plots of vertical wind profiles for the 4 seasons, featuring
%%%% ERA5, WACCM and the KEP radar over South Georgia.


site = 'kep';

years = 2016:2020;

seasons = {'djf','mam','jja','son'};
season_numbers = {[12 1 2],[3 4 5],[6 7 8],[9 10 11]};

%%%% LOAD ERA5
disp(['Loading ERA5...'])
load(['/Users/neil/Drive/MATLAB/ERA5/era5_kep_daily_uvt.mat'])

%%%% LOAD WACCM
disp(['Loading WACCM...'])
nc = getnet(['/Users/neil/data/WACCM/' site '/WACCM_monthly_2000-2014_324E_54S.nc']);
Waccm = struct;
Waccm.u = sq(nanmean(nc.Data.u,2));
Waccm.v = sq(nanmean(nc.Data.v,2));
Waccm.T = sq(nanmean(nc.Data.T,2));
Waccm.wtime = datenum(num2str(nc.Data.month),'yyyymm') + 15; % centre on the 15th of each month
Waccm.walt = nc.Data.altitude;

%%%% LOAD KEP RADAR
disp(['Loading KEP radar...'])
Comp.u = [];
Comp.v = [];
Comp.walt = [];
Comp.wtime = [];
for yr = years
    yrstr = num2str(yr);
    load(['/Users/neil/data/MeteorRadar/' site '/matlab/hwd/' yrstr '_' site '_hwd_gaussbins.mat'])
    Comp.u      = cat(2,Comp.u,sq(nanmean(HWD.Data.MonthlyComp.u,2)));
    Comp.v      = cat(2,Comp.v,sq(nanmean(HWD.Data.MonthlyComp.v,2)));
    Comp.walt   = cat(2,Comp.walt,sq(nanmean(HWD.Data.MonthlyComp.walt,2)));
    Comp.wtime  = cat(2,Comp.wtime,sq(nanmean(HWD.Data.MonthlyComp.wtime,2)));
end

inds = struct;

era5 = struct;
era5.alt = pres2alt(ecmwfpressurelevels);

waccm = struct;
waccm.alt = Waccm.walt;

kep = struct;
kep.alt = nanmean(Comp.walt,2);

for s = 1:4
    % era5...
    inds.era5.(seasons{s}) = ismember(month(ERA5.dayrange),season_numbers{s});
    era5.(seasons{s}).u = nanmean(ERA5.u(:,inds.era5.(seasons{s})),2);
    era5.(seasons{s}).v = nanmean(ERA5.v(:,inds.era5.(seasons{s})),2);
    
    % waccm...
    inds.waccm.(seasons{s}) = ismember(month(Waccm.wtime),season_numbers{s});
    waccm.(seasons{s}).u = nanmean(Waccm.u(:,inds.waccm.(seasons{s})),2);
    waccm.(seasons{s}).v = nanmean(Waccm.v(:,inds.waccm.(seasons{s})),2);
    
    % radar...
    inds.kep.(seasons{s}) = ismember(month(nanmean(Comp.wtime,1)),season_numbers{s});
    kep.(seasons{s}).u = nanmean(Comp.u(:,inds.kep.(seasons{s})),2);
    kep.(seasons{s}).v = nanmean(Comp.v(:,inds.kep.(seasons{s})),2);
    
    %%%% trim to heights:
    era5.(seasons{s}).u(~inrange(era5.alt,[0 75])) = NaN;
    era5.(seasons{s}).v(~inrange(era5.alt,[0 75])) = NaN;
    
    %%%% trim to heights:
    era5.(seasons{s}).u(~inrange(era5.alt,[0 75])) = NaN;
    era5.(seasons{s}).v(~inrange(era5.alt,[0 75])) = NaN;
    
end





%% PLOTTINGS ==============================================================

figure; hold all; whitefig; figpos([1 0.45])

%-------------------------------------------------------
vert_gap = 0.05;        horz_gap = 0.075;
lower_marg = 0.15;     upper_marg = 0.15;
left_marg = 0.05;      right_marg = 0.05;

rows = 1; cols = 8;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 22;

% T = tiledlayout(1,4);

% spans = {...
%     [1 1],...
%     [1 1],...
%     [1 1],...
%     [1 1]};

axs = {...
    [1 2],...
    [3 4],...
    [5 6],...
    [7 8]};

TBP = {...
    'djf',...
    'mam'...
    'jja',...
    'son'};

% colors = {'k',mcolor(2),mcolor(6)};
% colors = {mcolor(1),mcolor(1),mcolor(1)};

Lines = struct;

meridcolor = mcolor(2);

for ax = 1:4
    
    subplot(rows,cols,axs{ax})
    %     nexttile(spans{ax})
    
    axx = gca;
    
    grid on;
    
    ylim([0 110])
    ytick(0:20:120);
    yminortick(0:5:120);
    
    %%%% axes label
    %     text(-20,135,['(' alphabet(ax) ')'],'fontsize',1.75*fs,'fontweight','normal','HorizontalAlignment','left');
    %     text(0,135,upper(seasons{ax}),'fontsize',1.5*fs,'fontweight','bold','HorizontalAlignment','center');
    %     text(-20,100,['(' alphabet(ax) ')'],'fontsize',1.75*fs,'fontweight','normal','HorizontalAlignment','left');
    %     text(-15,100,upper(seasons{ax}),'fontsize',1.5*fs,'fontweight','bold','HorizontalAlignment','center');
    nph_text([0.025 0.875],['(' alphabet(ax+4) ')'],'fontsize',1.35*fs,'fontweight','normal','HorizontalAlignment','left','verticalalignment','middle');
    
    switch ax % deal with f being thin
        case 2
            nph_text([0.025 0.868],['     ' upper(seasons{ax})],'fontsize',1.1*fs,'fontweight','bold','HorizontalAlignment','left','verticalalignment','middle');
        otherwise
            nph_text([0.025 0.868],['      ' upper(seasons{ax})],'fontsize',1.1*fs,'fontweight','bold','HorizontalAlignment','left','verticalalignment','middle');
    end
    
    % zero line
    hold on; plot([0 0],axx.YLim,'linewi',2,'linest','--','color',rgbtrip(0.5))
    
%     switch ax
%         case 1
            ylabel('Altitude (km)')
%     end
    
    zonallinespec = {'linewi',3};
    meridlinespec = {'linewi',3};
    
    setfont(fs)
    set(gca,'linewi',1.5,'tickdir','out')
    
    %%%% merid wind...
    axx.XAxisLocation = 'top';
    axx.TickDir = 'out';
    axx.XColor = meridcolor;
    axx.GridColor = [0.15 0.15 0.15];
    xlim([-25 25])
    xtick(-40:10:40)
    xminortick(-40:5:40)
    xlabel('\bf{\it{v}} \rm{(ms^{-1})}','fontsize',fs)
%     %%%% waccm
%     hold on; plot(axx,waccm.(seasons{ax}).v,waccm.alt,meridlinespec{:},'color',mcolor(2),'linest','--');
    %%%% era5
    hold on; Lines(1).line = plot(axx,era5.(seasons{ax}).v,era5.alt,meridlinespec{:},'color',meridcolor,'linest','-');
    %%%% kep
    hold on; Lines(2).line = plot(axx,kep.(seasons{ax}).v,kep.alt,meridlinespec{:},'color',meridcolor,'linest',':');
    
%     % line around the axis
%     hold on; plot(axx.XLim([1 2 2]),axx.YLim([2 2 1]),'color','k','linewi',1.5)
%     
    
    %%%% zonal wind...
    axx2 = nph_copyaxes(axx);
    axx2.XColor = 'k';
    axx2.GridColor = [0.15 0.15 0.15];
    axx2.YTick = [];
    xlim([-100 100])
    xtick(-200:40:200)
    xminortick(-200:10:200)
    xlabel('\bf{\it{u}} \rm{(ms^{-1})}','fontsize',fs)
%     %%%% waccm
%     hold on; plot(axx2,waccm.(seasons{ax}).u,waccm.alt,zonallinespec{:},'color',mcolor(1),'linest','--');    
    %%%% era5
    hold on; Lines(3).line = plot(axx2,era5.(seasons{ax}).u,era5.alt,zonallinespec{:},'color','k','linest','-');    
    %%%% kep
    hold on; Lines(4).line = plot(axx2,kep.(seasons{ax}).u,kep.alt,zonallinespec{:},'color','k','linest',':');    
    
    % line around the axis
    hold on; plot(axx2.XLim([2 2]),axx2.YLim([1 2]),'color','k','linewi',1.5)
    
    
    switch ax
        case 4
            % legend?
            L = legend([Lines([4 2 3 1]).line],{' \bf\it{u}\rm{, Radar} ',' \bf\it{v}\rm{, Radar} ',' \bf\it{u}\rm{, ERA5}   ',' \bf\it{v}\rm{, ERA5}  '});
    end
    
%     
%     % fix the xlims for merid...
%     axx.XLim = [-20 20];
%     axx.XTick = -40:10:40;
    
    drawnow;
    
end


return


%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_seasonal_wind_profiles'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')
















