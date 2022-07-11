

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
Comp.tt = [];
Comp.day = [];
Comp.monthlyfwhm = [];
Comp.monthlycenter = [];
Comp.month = [];
Comp.gwvar = [];
Comp.walt = [];
Comp.wtime = [];

for yr = years
    yrstr = num2str(yr);
    load(['/Users/neil/data/MeteorRadar/' site '/matlab/hwd/' yrstr '_' site '_hwd_gaussbins.mat'])
    Comp.u      = cat(2,Comp.u,sq(nanmean(HWD.Data.MonthlyComp.u,2)));
    Comp.v      = cat(2,Comp.v,sq(nanmean(HWD.Data.MonthlyComp.v,2)));
    Comp.walt   = cat(2,Comp.walt,sq(nanmean(HWD.Data.MonthlyComp.walt,2)));
    Comp.wtime  = cat(2,Comp.wtime,sq(nanmean(HWD.Data.MonthlyComp.wtime,2)));
    Comp.uu  = cat(2,Comp.uu,HWD.Data.u);
    Comp.vv  = cat(2,Comp.vv,HWD.Data.v);
    Comp.tt  = cat(2,Comp.tt,HWD.Data.Time);
    %     Comp.fwhm   = cat(2,Comp.fwhm,HWD.Data.VertMetDist.FWHM);
%     Comp.center = cat(2,Comp.center,HWD.Data.VertMetDist.Center);
    Comp.day    = cat(2,Comp.day,HWD.Data.Day);
    Comp.monthlyfwhm   = cat(2,Comp.monthlyfwhm,HWD.Data.MonthlyComp.VertMetDist.FWHM);
    Comp.monthlycenter = cat(2,Comp.monthlycenter,HWD.Data.MonthlyComp.VertMetDist.Center);
    Comp.month    = cat(2,Comp.month,datenum(yr,1:12,01));
end

% return


% % % % % next, load all the mpd files for the year range and cat them into a
% % % % % MASSIVE array. This will be HUGE I hope it fits in memory.
% % % % disp(['Loading all MPDs...'])
% % % % mpddirec = ['/Users/neil/data/MeteorRadar/' site '/matlab/'];
% % % % 
% % % % allMPD = struct;
% % % % dayrange = datenum(years(1),01,01):datenum(years(end),12,31);
% % % % 
% % % % for d = 1:length(dayrange)
% % % %     
% % % %     yyyymmdd = datestr(dayrange(d),'yyyymmdd');
% % % %     
% % % %     % pre assign length
% % % %     allMPD(d).Time = [];
% % % %     
% % % %     try
% % % %         load([mpddirec yyyymmdd '_' site '_mpd.mat'])
% % % %         
% % % %         % if it exists, subscribe data:
% % % %         % Subscribe:
% % % %         allMPD(d).Time                  = MPD.Data.Time;
% % % %         allMPD(d).Alt                   = MPD.Data.Alt;
% % % %         allMPD(d).RadialVelocity        = MPD.Data.RadialVelocity;
% % % %         allMPD(d).ZenithAngle           = MPD.Data.ZenithAngle;
% % % %         allMPD(d).Azimuth               = MPD.Data.Azimuth;
% % % %         %allMPD(d).AngleNorthFromEast    = wrapTo360(90 - MPD.Data.Azimuth);
% % % %         
% % % %     catch err
% % % %         disp(['Couldn''t load ' yyyymmdd '... ' err.identifier])
% % % %         continue
% % % %     end
% % % %     
% % % %     
% % % % end
% % % % 
% % % % % Concatenate useful quantities:
% % % % zen     = cat(1,allMPD(:).ZenithAngle);
% % % % az      = cat(1,allMPD(:).Azimuth);
% % % % tim     = cat(1,allMPD(:).Time);
% % % % alt     = cat(1,allMPD(:).Alt);
% % % % 
% % % % vrad    = cat(1,allMPD(:).RadialVelocity);
% % % % vhorz   = vrad ./ cosd(zen);
% % % % 
% % % % % Apply Zenith Angle Limits:
% % % % MWD.Thresholds.ZenithLimits = [15 65]; % use [15 65] for computing mean winds, impose stricter for GWs.
% % % % inds = zen >= MWD.Thresholds.ZenithLimits(1) & zen <= MWD.Thresholds.ZenithLimits(2);
% % % % 
% % % % % avoid silly horz velocities:
% % % % inds = inds & vhorz < 200;
% % % % 
% % % % az = az(inds);
% % % % tim = tim(inds);
% % % % alt = alt(inds);
% % % % vrad = vrad(inds);
% % % % vhorz = vhorz(inds);
% % % % zen = zen(inds);
% % % % 
% % % % hourofday = mod(tim,1) * 24;
% % % % dn = daynumber(tim);
% % % % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOMB-SCARGLE

disp('Lomb-Scargle...')

zind = 15;

uu = Comp.uu(zind,:); vv = Comp.vv(zind,:); tt = Comp.tt;
nanlocs = isnan(uu) | isnan(vv) | isnan(tt);

freqs = [0.02:0.0025:10];
ftype = 'normalized';
[pxxu,f] = plomb(uu(~nanlocs),tt(~nanlocs),freqs,ftype);
[pxxv,f] = plomb(vv(~nanlocs),tt(~nanlocs),freqs,ftype);

% sz = size(uu);
% uu = nanmean(reshape(uu,[2 sz(2)/2]),1);
% vv = nanmean(reshape(vv,[2 sz(2)/2]),1);

% % %%

fsamp = 400;

% try and do a confidence interval?
Pfa = [50 10 1 0.01]/100;
Pd = 1-Pfa;
[~,~,pthu] = plomb(uu,fsamp,ftype,'Pd',Pd);
[~,~,pthv] = plomb(vv,fsamp,ftype,'Pd',Pd);


return



%% ZONAL AND MERID WINDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([1 0.475])

%-------------------------------------------------------
vert_gap = 0.05;        horz_gap = 0.05;
lower_marg = 0.15;       upper_marg = 0.125;
left_marg = 0.05;      right_marg = 0.05;

rows = 1; cols = 2;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 20;

% T = tiledlayout(1,2,'tilespacing','compact');

spans = {...
    [1 1],...
    [1 1],...
    [1 1],...
    [1 1]};

% axs = {...
%     [1  6],...
%     [9 14],...
%     [7 8],...
%     [15 16]};

TBP = {...
    sqrt(pxxu),...
    sqrt(pxxv),...
    Comp.u,...
    Comp.v};

linecolors = {...
    mcolor(1),...
    (mcolor(3)+mcolor(7))./2};

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
    'ZONAL',...
    'MERIDIONAL',...
    'ZONAL WIND',...
    'MERIDIONAL WIND'};


% also contour lines:
clines = {...
    [-60 -40 -20 20 40 60],...
    [-20 -10 10 20],...
    [-60 -40 -20 20 40 60],...
    [-20 -10 10 20]};


letters = {'a','c','b','d'};

ynudge = 0;
% ylims = [77 103];
ylims = [77+ynudge 103];


for ax = 1:2
    
%     nexttile(spans{ax})
    subplot(rows,cols,ax)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Cycles per day axes:
    axx = gca;
    tbp = TBP{ax};
    
    grid on;
    
    xlim([0.01 20])
    axx.XTickLabel = strtrim(cellstr(num2str(axx.XTick')));
    xlabel('Frequency [cycles per day]')
    
    ylim([1E-2 1E2])
    axx.YMinorGrid = 'off';
    axx.YTick = 10.^(-8:8);
    axx.YTickLabel = strtrim(cellstr(num2str(log10(axx.YTick)')));
    ylabel('Relative Power (10^x)')
    
    % first, shaded regions for PWs!
    pws = 1./[1.67 2.5 ; 4 6 ; 8.5 11.5 ; 13 18];
    for p = 1:size(pws,1)
        hold on; P = patch(pws(p,[1 2 2 1 1]),axx.YLim([1 1 2 2 1]),'k');
        P.FaceAlpha = 0.05;
        P.EdgeColor = 'none';
    end
    
    % approx 90th confidence interval
    switch ax
        case 1
            hold on; plot(axx.XLim,0.6*sqrt(pthu([2 2])),'linewi',1,'color',[rgbtrip(0) 0.5],'linest','--');
        case 2
            hold on; plot(axx.XLim,0.6*sqrt(pthv([2 2])),'linewi',1,'color',[rgbtrip(0) 0.5],'linest','--');
    end
    
    % intertial freq
    hold on; plot([24 24]./inertialperiod(-54),axx.YLim,'linewi',1,'color',[0 0 0 1],'linest',':');
    
    % now the PSD data...
    hold on; plot(f,movmean(tbp,3,1),'color',linecolors{ax},'linewi',2);
    
    %%%% add a -5/3 line:
%     x = linspace(axx.XLim(1),axx.XLim(2),100);
%     x = linspace(axx.XLim(1),axx.XLim(2),100);
    x = [6 24./inertialperiod(-54)];
    y = x.^(-5/3);
    hold on; plot(x,4.8.*y,'r--');
    
%     %%%% fit a straight line to the GW portion:
%     intp = 24./inertialperiod(-54);
%     inds = f >= intp;
%     x = f(inds); y = tbp(inds);
%     %     P = polyfit(log10(x),log10(y),1);
%     switch ax
%         case 1
%             hold on; plot(x,(x.^(-5/3))*3.5,'linewi',2,'color',[mcolor(7) 1],'linest','--');
%         case 2
%             hold on; plot(x,(x.^(-5/3))*3.75,'linewi',2,'color',[mcolor(7) 1],'linest','--');
%             
%     end

    setfont(fs);
    
    % line around the axis
    hold on; plot(axx.XLim([1 2 2]),axx.YLim([2 2 1]),'color','k','linewi',1.5)
    
    % axes letter
    hold on; nph_text([0.02 0.875],['(' alphabet(ax) ')'],'fontsize',1.75*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
    hold on; nph_text([0.08 0.87],tits{ax},'fontsize',1.25*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
    
    set(gca,'linewi',1.5,'tickdir','out')
    
    logscale('x','y')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Periods axes:
    axx2 = nph_copyaxes(axx);
    axx2.YColor = 'none';
    axx2.XAxisLocation = 'top';
    
    axx2.XScale = 'log';
    
    perioddays = [50 20 10 5 2];
    periodhours = [1 2 3 4 6];
    axx2.XTick = [1./perioddays periodhours];
    axx2.XTickLabels = cat(1,...
        cellstr([num2str(perioddays')  ('ddddd')']),...
        cellstr([num2str(24./periodhours') ('hhhhh')']));
    xlabel('Period [days, hours]')
    
    set(gca,'linewi',1.5,'tickdir','out')
    
    
end


return


%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_log_power_spectrum_with_line'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return

















