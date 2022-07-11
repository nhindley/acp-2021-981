

%%%% description: non-linear tidal planteray wave coupling

%%%% Firstly, makes a figure of mean winds u,v from KEP, then below that an
%%%% S-transform spectrum of both u,v to reveal PW activity.


runagain = 1;

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
    
    % raw hourly winds
    Comp.uu     = [];
    Comp.vv     = [];
    Comp.tt     = [];
    Comp.walt   = [];
    Comp.wtime  = [];
    Comp.dd     = [];
    
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
    
%     Comp.Tides.OneDay.AllTides.u = [];
%     Comp.Tides.OneDay.AllTides.v = [];
    
    % and some mean winds?
    Comp.MeanWinds.u = [];
    Comp.MeanWinds.v = [];
    Comp.MeanWindsWithDailyMean.u = [];
    Comp.MeanWindsWithDailyMean.v = [];
    
    disp('Loading HWD(s)...')
    
    for yr = years
        
        yrstr = num2str(yr);
        
        load([direc site '/matlab/hwd/' yrstr '_' site '_hwd_gaussbins.mat'])
        
        % raw wind
        Comp.uu      = cat(2,Comp.uu,HWD.Data.u);
        Comp.vv      = cat(2,Comp.vv,HWD.Data.v);
        Comp.tt      = cat(2,Comp.tt,HWD.Data.Time);
        Comp.dd      = cat(2,Comp.dd,HWD.Data.Day);
        
        Comp.walt = cat(2,Comp.walt,HWD.Data.walt);
        Comp.wtime = cat(2,Comp.wtime,HWD.Data.wtime);
        
%         % monthly mean winds
%         Comp.u      = cat(2,Comp.u,sq(nanmean(HWD.Data.MonthlyComp.u,2)));
%         Comp.v      = cat(2,Comp.v,sq(nanmean(HWD.Data.MonthlyComp.v,2)));
        
        % 4-day Tides
        Comp.Tides.u.A = cat(2,Comp.Tides.u.A,HWD.Data.Tides.TidalComponents.u.A);
        Comp.Tides.u.B = cat(2,Comp.Tides.u.B,HWD.Data.Tides.TidalComponents.u.B);
        Comp.Tides.v.A = cat(2,Comp.Tides.v.A,HWD.Data.Tides.TidalComponents.v.A);
        Comp.Tides.v.B = cat(2,Comp.Tides.v.B,HWD.Data.Tides.TidalComponents.v.B);
        
        % 30-day tides
        Comp.Monthly.Tides.u.A = cat(2,Comp.Monthly.Tides.u.A,HWD.Data.MonthlyComp.Tides.TidalComponents.u.A);
        Comp.Monthly.Tides.u.B = cat(2,Comp.Monthly.Tides.u.B,HWD.Data.MonthlyComp.Tides.TidalComponents.u.B);
        Comp.Monthly.Tides.v.A = cat(2,Comp.Monthly.Tides.v.A,HWD.Data.MonthlyComp.Tides.TidalComponents.v.A);
        Comp.Monthly.Tides.v.B = cat(2,Comp.Monthly.Tides.v.B,HWD.Data.MonthlyComp.Tides.TidalComponents.v.B);
        
        % 1-day tides
        Comp.Tides.OneDay.u.A = cat(2,Comp.Tides.OneDay.u.A,HWD.Data.Tides.OneDay.TidalComponents.u.A);
        Comp.Tides.OneDay.u.B = cat(2,Comp.Tides.OneDay.u.B,HWD.Data.Tides.OneDay.TidalComponents.u.B);
        Comp.Tides.OneDay.u.C = cat(2,Comp.Tides.OneDay.u.C,HWD.Data.Tides.OneDay.TidalComponents.u.C);
        Comp.Tides.OneDay.v.A = cat(2,Comp.Tides.OneDay.v.A,HWD.Data.Tides.OneDay.TidalComponents.v.A);
        Comp.Tides.OneDay.v.B = cat(2,Comp.Tides.OneDay.v.B,HWD.Data.Tides.OneDay.TidalComponents.v.B);
        Comp.Tides.OneDay.v.C = cat(2,Comp.Tides.OneDay.v.C,HWD.Data.Tides.OneDay.TidalComponents.v.C);
        
        % residual winds after tides removed
        Comp.MeanWinds.u = cat(2,Comp.MeanWinds.u,HWD.Data.u-HWD.Data.Tides.OneDay.AllTides.u);
        Comp.MeanWinds.v = cat(2,Comp.MeanWinds.v,HWD.Data.v-HWD.Data.Tides.OneDay.AllTides.v);
        
        % or, remove tides but NOT the daily mean:
        dailymean_u = interp1(HWD.Data.Day,HWD.Data.Tides.OneDay.TidalComponents.u.C(:,:,4)',HWD.Data.Time)';
        dailymean_v = interp1(HWD.Data.Day,HWD.Data.Tides.OneDay.TidalComponents.v.C(:,:,4)',HWD.Data.Time)';
        Comp.MeanWindsWithDailyMean.u = cat(2,Comp.MeanWindsWithDailyMean.u,HWD.Data.u-HWD.Data.Tides.OneDay.AllTides.u+dailymean_u);
        Comp.MeanWindsWithDailyMean.v = cat(2,Comp.MeanWindsWithDailyMean.v,HWD.Data.v-HWD.Data.Tides.OneDay.AllTides.v+dailymean_v);
        
        
    end
    
    if strcmpi(site,'rothera-sk')
        site = 'rothera';
    end
    eval([upper(site) '.Comp = Comp;'])
%     disp('Loading riogrande too...')
%     load('/Users/neil/Drive/MATLAB/riograndePWs.mat')
    
    save(['/Users/neil/Drive/MATLAB/' site 'PWs.mat'],upper(site));
    
    
else
%     disp('loading...')
%     load('/Users/neil/Drive/MATLAB/MonthlyCompTidesData.mat')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOW S-TRANSFORM BOTH THE WINDS AND THE TIDAL AMPLITUDES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('NDST...')

zlev = 20; % 20 = 95km, 15 = 90km

meanwinds_u = movmean(Comp.MeanWinds.u(zlev,:),6,2,'omitnan');
meanwinds_v = movmean(Comp.MeanWinds.v(zlev,:),6,2,'omitnan');

Comp.dd6h = datenum(min(years),01,01):(6/24):datenum(max(years)+1,01,01);
meanwinds_u = interp1(Comp.tt,meanwinds_u,Comp.dd6h);
meanwinds_v = interp1(Comp.tt,meanwinds_v,Comp.dd6h);

% meanwinds_u = Comp.Tides.OneDay.u.C(zlev,:,4);
% meanwinds_v = Comp.Tides.OneDay.v.C(zlev,:,4);

nanlocs = isnan(meanwinds_u) | isnan(meanwinds_u);
meanwinds_u(nanlocs) = nanmean(meanwinds_u(:));
meanwinds_v(nanlocs) = nanmean(meanwinds_v(:));

%%%% ST the raw winds
c = 1;
point_spacing = 6/24;
scales = unique(floor(length(meanwinds_u)./((1:0.05:100)./point_spacing)));

% %%%% take away a 30 day smooth...
% meanwinds_u = meanwinds_u - movmean(meanwinds_u,30/point_spacing,2,'omitnan');
% meanwinds_v = meanwinds_v - movmean(meanwinds_v,30/point_spacing,2,'omitnan');

Comp.MeanWinds.STu = nph_ndst(meanwinds_u,scales,point_spacing,c);
Comp.MeanWinds.STv = nph_ndst(meanwinds_v,scales,point_spacing,c);

% replace nans:
Comp.MeanWinds.STu.ST(:,nanlocs) = NaN;
Comp.MeanWinds.STv.ST(:,nanlocs) = NaN;
meanwinds_u(nanlocs) = NaN;
meanwinds_v(nanlocs) = NaN;


%%%% ST the daily tidal amplitudes
% scales = unique(floor(length(meanwinds_u)./(1:0.05:100)));
c = 1;
point_spacing = 6/24;
scales = unique(floor(length(meanwinds_u)./((1:0.05:100)./point_spacing)));


td = 2; % 12h tide
tidalamp_u = quadadd(Comp.Tides.OneDay.u.A(zlev,:,td),Comp.Tides.OneDay.u.B(zlev,:,td));
tidalamp_v = quadadd(Comp.Tides.OneDay.v.A(zlev,:,td),Comp.Tides.OneDay.v.B(zlev,:,td));

nanlocs = isnan(tidalamp_u) | isnan(tidalamp_v);
tidalamp_u(nanlocs) = nanmean(tidalamp_u(:));
tidalamp_v(nanlocs) = nanmean(tidalamp_v(:));

Comp.Tides.OneDay.STu = nph_ndst(tidalamp_u,scales,point_spacing,c);
Comp.Tides.OneDay.STv = nph_ndst(tidalamp_v,scales,point_spacing,c);

% replace nans:
Comp.Tides.OneDay.STu.ST(:,nanlocs) = NaN;
Comp.Tides.OneDay.STv.ST(:,nanlocs) = NaN;
tidalamp_u(nanlocs) = NaN;
tidalamp_v(nanlocs) = NaN;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PW AMPLITUDES AND S-TRANSFORM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([1 0.8])

%-------------------------------------------------------
vert_gap = 0.08;        horz_gap = 0.03;
lower_marg = 0.075;     upper_marg = 0.025;
left_marg = 0.065;      right_marg = 0.075;

rows = 2; cols = 8;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 20;

T = tiledlayout(3,6,'tilespacing','compact');

spans = {...
    [1 6],...
    [1 6],...
    [1 6],...
    [1 6]};

% axs = {...
%     [1  6],...
%     [9 14],...
%     [7 8],...
%     [15 16]};

TBP = {...
    [meanwinds_u ; meanwinds_v],...
    abs(Comp.MeanWinds.STu.ST),...
    abs(Comp.MeanWinds.STv.ST),...
    []};
% TBP = {...
%     [meanwinds_u ; meanwinds_v],...
%     sqrt(abs(Cu)),...
%     sqrt(abs(Cv)),...
%     []};

clims = {...
    [-20 80],...
    [-1 20],...
    [-1 20],...
    [-1 20]};

clevs = {...
    [0:20],...
    [0:20],...
    [0:20],...
    [-20:2.5:20]};

clabs = {...
    'u','v','u','v'};

cbarticks = {...
    0:5:20,...
    0:5:20,...
    0:5:20,...
    0:5:20};

tits = {...
    'WIND (TIDES REMOVED)',...
    'ZONAL AMPLITUDE',...
    'MERIDIONAL AMPLITUDE',...
    ' '};

cmaps = {...
    cbrew('nph_OSOrRdPu'),...
    cbrew('nph_OSBlues'),...
    cbrew('nph_OSOrRdPu'),...
    cbrew('nph_OSOrRdPu')};


% also contour lines:
blackclines = {...
    [0 0],...
    [0 10 20],...
    [0 10 20],...
    [0 0]};
whiteclines = {...
    [0 0],...
    [20 30],...
    [20 30],...
    [0 0]};

letters = {'a','b','d','d'};

coiflag = 1;

for ax = 1:3
    
    nexttile(spans{ax})
    
    axx = gca;
    
    xlim([datenum(2016,02,15) datenum(2020,11,15)])
    
    switch ax
        case 1

            tbp = TBP{ax};
            hold on; plots.v = plot(Comp.dd6h,tbp(2,:),'color',mcolor(2),'linewi',1.5);
            hold on; plots.u = plot(Comp.dd6h,tbp(1,:),'color',mcolor(1),'linewi',1.5);
            
            grid on;
            
            ylim([-50 90])
            ytick(-100:20:100);
            yminortick(-100:10:100);
            
            ylabel('ms^-^1')
            
            L = legend([plots.u plots.v],{'Zonal Wind','Meridional Wind'},'fontsize',0.8*fs,'AutoUpdate','off');
            
            %%%% MONTHS
            axx.XTick = datenum(years,01,01);
            axx.XMinorTick = 'on';
            axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
            datetick('x',' ','keepticks','keeplimits')
            % other months as ticks using text:
            xtix = datenum(years(1),1:(length(years)*12),15);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,axx.YLim(1)-5,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
               end
            end
%             % Years as lower bold numbers:
%             xtix = datenum(years,07,01);
%             for xt = xtix
%                if inrange(xt,xlim) && ~any(xt == axx.XTick)
%                    yrstr = datestr(xt,'yyyy');
%                    text(xt,axx.YLim(1)-25,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
%                end
%             end
        
        % allll the contours...
        case {2,3}
            
            tbp = TBP{ax}; nanlocs = isnan(tbp);
            % smooth nans...
            smoo = movmean(tbp,7*4,2,'omitnan');
            tbp(nanlocs) = smoo(nanlocs);
            % put a few nans back in...
            smoosmoo = movmean(tbp,7,2);
            tbp(isnan(smoosmoo)) = NaN;
            
            % grids...
            [X,Y] = meshgrid(Comp.dd6h,1./Comp.MeanWinds.STu.freqs);
            smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
            
            ylim([1.4 35])
            axx.YTick = [0 2 5 10 15 20 30 40 50 100];
            axx.YAxis.MinorTickValues = [0:1:20 25:5:100];
            
            % hatched area below:
%             axx2 = axes;
%             axes(axx2);
%             set(axx2,'color','none','xcolor','none','ycolor','none',...
%                 'position',axx.Position,'xlim',axx.XLim);
            h = patch(axx.XLim([1 2 2 1]),axx.YLim([1 1 2 2]),'red');
            hh = hatchfill(h, 'cross', 45, 8);
            set(hh,'color',[0 0 0 0.15],'linewi',1)
%             axes(axx);
%             drawnow;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% CONE OF INFLUENCE!
            if coiflag
            % scroll through, and every time you hit a NaN, put a cone of
            % influence mask there:
            startmarkers = datenum(2016,01,01);
            endmarkers = [];
            for i = 2:size(tbp,2)
               if ~all(isnan(tbp(:,i-1))) && all(isnan(tbp(:,i)))
                   startmarkers(end+1) = X(1,i);
               end
               if all(isnan(tbp(:,i-1))) && ~all(isnan(tbp(:,i)))
                   endmarkers(end+1) = X(1,i);
               end
            end
            endmarkers(end+1) = datenum(2020,12,31);
            % now set these limits to NaN:
            for i = 1:length(startmarkers)
%             for i = 2
                coi = X > (startmarkers(i) - Y) & X < (endmarkers(i) + Y);
%                 coi = X > (startmarkers(i) + Y) & X < (endmarkers(i) - Y);
                tbp(coi) = NaN;
            end
            end
%             return
%             coi = X > (axx.XLim(1) + Y) & X < (axx.XLim(2) - Y);
%             tbp(~coi) = NaN;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % filled contours:
            hold on; contourf(X,Y,tbp,clevs{ax},'edgecolor','none');
            
            % contour lines:
            hold on; [C,h] = contour(X,Y,tbp,blackclines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
            clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.25),'fontweight','bold');
            
            hold on; [C,h] = contour(X,Y,tbp,[20 20],'edgecolor',rgbtrip(.95),'linewi',1.25);
%             clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
            
            hold on; [C,h] = contour(X,Y,tbp,[30 30],'edgecolor',rgbtrip(.95),'linewi',1.25);
            clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
            
            
            logscale('y')
            ylabel('Period (days)')
            
            axx.YGrid = 'on';
            axx.YMinorGrid = 'off';
            
            clim(clims{ax})
            
            %%%% COLORS
            cmap = [1 1 1 ; cmaps{ax}];
            nclev = 11;
            cmapi = interp1(1:size(cmap,1),cmap,linspace(1,size(cmap,1),nclev),'linear');
            colormap(gca,cmapi)
            
            %%%% COLORBAR
            cbar = nph_colorbar;
            cbar.Ticks = cbarticks{ax};
            cbar.TickDirection = 'out';
            cbar.Position = cbar.Position .* [1 1 0.4 1];
            cbar.Position = cbar.Position + [-1*cbar.Position(3) 0 0 0];
%             cbar.Label.String = ['\bf{' clabs{ax} '}'];
            cbar.Label.Rotation = 0;
            cbar.Label.VerticalAlignment = 'middle';
%             cbar.Label.HorizontalAlignment = 'left';
            cbar.Limits(1) = 0;
            cbar.Title.String = '   ms^{-1}';
            
            cbar.Ruler.MinorTick = 'on';
            cbar.Ruler.MinorTickValues = 0:2.5:30;
            
            %%%% MONTHS
            axx.XTick = datenum(years,01,01);
            axx.XMinorTick = 'on';
            axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
            datetick('x',' ','keepticks','keeplimits')
            % other months as ticks using text:
            xtix = datenum(years(1),1:(length(years)*12),15);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,1.25,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
               end
            end
            % Years as lower bold numbers:
            if ax == 3
            xtix = datenum(years,07,01);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   yrstr = datestr(xt,'yyyy');
                   text(xt,0.75,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
               end
            end
            end
            
    end
    
    
    % AXES LETTER
    %letter = alphabet(ax);
    letter = letters{ax};
    switch ax
        case 1
            hold on; nph_text([0.005 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
            hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        case {2,3}
            hold on; nph_text([0.005 0.81],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
            hold on; nph_text([0.005 0.80],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
    end

    setfont(fs)
    
    set(gca,'linewi',1.5,'layer','top','tickdir','out')
    
    %%%% bold line around:
    hold on; plot(axx.XLim([1 2 2 1 1]),axx.YLim([1 1 2 2 1]),'linewi',1.5,'color','k');
    
    
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_PWs'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% PW PHASES AND E/W MODES USING RIOGRANDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('reloading data...')

load('/Users/neil/Drive/MATLAB/kepPWs.mat')
load('/Users/neil/Drive/MATLAB/riograndePWs.mat')
load('/Users/neil/Drive/MATLAB/rotheraPWs.mat')

disp('NDST...')

zlev = 20; % 20 = 95km, 15 = 90km

for i = 1:3
    
switch i
    case 1
        data = KEP.Comp;
    case 2
        data = RIOGRANDE.Comp;
    case 3
        data = ROTHERA.Comp;    
end
        
% data.meanwinds_u = movmean(data.MeanWinds.u(zlev,:),6,2,'omitnan');
% data.meanwinds_v = movmean(data.MeanWinds.v(zlev,:),6,2,'omitnan');
data.meanwinds_u = movmean(data.MeanWindsWithDailyMean.u(zlev,:),6,2,'omitnan');
data.meanwinds_v = movmean(data.MeanWindsWithDailyMean.v(zlev,:),6,2,'omitnan');


data.dd6h = datenum(min(years),01,01):(6/24):datenum(max(years)+1,01,01);
data.meanwinds_u = interp1(data.tt,data.meanwinds_u,data.dd6h);
data.meanwinds_v = interp1(data.tt,data.meanwinds_v,data.dd6h);

nanlocs = isnan(data.meanwinds_u) | isnan(data.meanwinds_u);
data.meanwinds_u(nanlocs) = nanmean(data.meanwinds_u(:));
data.meanwinds_v(nanlocs) = nanmean(data.meanwinds_v(:));

%%%% ST the raw winds
c = 1;
point_spacing = 6/24;
scales = unique(floor(length(data.meanwinds_u)./((1:0.05:100)./point_spacing)));

% %%%% take away a 30 day smooth...
% meanwinds_u = meanwinds_u - movmean(meanwinds_u,30/point_spacing,2,'omitnan');
% meanwinds_v = meanwinds_v - movmean(meanwinds_v,30/point_spacing,2,'omitnan');

data.MeanWinds.STu = nph_ndst(data.meanwinds_u,scales,point_spacing,c);
data.MeanWinds.STv = nph_ndst(data.meanwinds_v,scales,point_spacing,c);

% replace nans:
data.MeanWinds.STu.ST(:,nanlocs) = NaN;
data.MeanWinds.STv.ST(:,nanlocs) = NaN;
data.meanwinds_u(nanlocs) = NaN;
data.meanwinds_v(nanlocs) = NaN;

switch i
    case 1
        KEP.Comp = data;
    case 2
        RIOGRANDE.Comp = data;
    case 3
        ROTHERA.Comp = data;    
end

end

Cu = RIOGRANDE.Comp.MeanWinds.STu.ST.*conj(KEP.Comp.MeanWinds.STu.ST);
Cv = RIOGRANDE.Comp.MeanWinds.STv.ST.*conj(KEP.Comp.MeanWinds.STv.ST);

% Cu = ROTHERA.Comp.MeanWinds.STu.ST.*conj(KEP.Comp.MeanWinds.STu.ST);
% Cv = ROTHERA.Comp.MeanWinds.STv.ST.*conj(KEP.Comp.MeanWinds.STv.ST);

%%

%-------------------------------------------------------


figure; hold all; whitefig; figpos([1 0.8])

%-------------------------------------------------------
vert_gap = 0.08;        horz_gap = 0.03;
lower_marg = 0.075;     upper_marg = 0.025;
left_marg = 0.065;      right_marg = 0.075;

rows = 2; cols = 8;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 20;

T = tiledlayout(3,6,'tilespacing','compact');

spans = {...
    [1 6],...
    [1 6],...
    [1 6],...
    [1 6]};

TBPamp = {...
    [meanwinds_u ; meanwinds_v],...
    abs(Comp.MeanWinds.STu.ST),...
    abs(Comp.MeanWinds.STv.ST),...
    []};
TBP = {...
    [meanwinds_u ; meanwinds_v],...
    Cu,...
    Cv,...
    []};

clims = {...
    [-20 80],...
    [-pi +pi],...
    [-pi +pi],...
    [-1 20]};

clevs = {...
    [0:20],...
    [-pi:pi/6:pi],...
    [-pi:pi/6:pi],...
    [-20:2.5:20]};

clabs = {...
    'u','v','u','v'};

cbarticks = {...
    0:5:20,...
    -pi:pi/6:pi,...
    -pi:pi/6:pi,...
    0:5:20};

cbarticklabels = {...
    {},...
    {'-\pi','-2\pi/3','-\pi/3','0','+\pi/3','+2\pi/3','+\pi'},...
    {'-\pi','-2\pi/3','-\pi/3','0','+\pi/3','+2\pi/3','+\pi'},...
    {}};

tits = {...
    ' ',...
    'ZONAL PHASE rel. RIO GRANDE',...
    'MERIDIONAL PHASE rel. RIO GRANDE',...
    ' '};

cmaps = {...
    cbrew('nph_OSOrRdPu'),...
    [flipud(cbrew('BuPu')) ; cbrew('RdPu')],...
    [flipud(cbrew('BuPu')) ; cbrew('RdPu')],...
    cbrew('nph_OSOrRdPu')};

% cmaps{2} = [flipud(cbrew('BuPu')) ; cbrew('RdPu')];

% also contour lines:
blackclines = {...
    [0 0],...
    [0 10],...
    [0 10],...
    [0 0]};
whiteclines = {...
    [0 0],...
    [0 0],...
    [0 0],...
    [0 0]};

letters = {'-','c','e','-'};

coiflag = 1;

for ax = 1:3
    
    nexttile(spans{ax})
    
    axx = gca;
    
    xlim([datenum(2016,02,15) datenum(2020,11,15)])
    
    switch ax
        case 1

            tbp = TBP{ax};
            hold on; plots.v = plot(Comp.dd6h,tbp(2,:),'color',mcolor(2),'linewi',1.5);
            hold on; plots.u = plot(Comp.dd6h,tbp(1,:),'color',mcolor(1),'linewi',1.5);
            
            grid on;
            
            ylim([-50 90])
            ytick(-100:20:100);
            yminortick(-100:10:100);
            
            ylabel('ms^-^1')
            
            L = legend([plots.u plots.v],{'Zonal Wind','Meridional Wind'},'fontsize',0.8*fs,'AutoUpdate','off');
            
            %%%% MONTHS
            axx.XTick = datenum(years,01,01);
            axx.XMinorTick = 'on';
            axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
            datetick('x',' ','keepticks','keeplimits')
            % other months as ticks using text:
            xtix = datenum(years(1),1:(length(years)*12),15);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,axx.YLim(1)-5,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
               end
            end
%             % Years as lower bold numbers:
%             xtix = datenum(years,07,01);
%             for xt = xtix
%                if inrange(xt,xlim) && ~any(xt == axx.XTick)
%                    yrstr = datestr(xt,'yyyy');
%                    text(xt,axx.YLim(1)-25,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
%                end
%             end
        
        % allll the contours...
        case {2,3}
            
            tbp = TBP{ax}; nanlocs = isnan(tbp);
            % smooth nans...
            smoo = movmean(tbp,7*4,2,'omitnan');
            tbp(nanlocs) = smoo(nanlocs);
            % put a few nans back in...
            smoosmoo = movmean(tbp,7,2);
            tbp(isnan(smoosmoo)) = NaN;
            
            %%%% convert to angle:
            tbp = angle(tbp);
            
            %%%% fix 2 clim:
            tbp = nph_fix2clims(tbp,clims{ax});
            
            % grids...
            [X,Y] = meshgrid(Comp.dd6h,1./Comp.MeanWinds.STu.freqs);
            smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
            
            ylim([1.4 35])
            axx.YTick = [0 2 5 10 15 20 30 40 50 100];
            axx.YAxis.MinorTickValues = [0:1:20 25:5:100];
            
            % hatched area below:
%             axx2 = axes;
%             axes(axx2);
%             set(axx2,'color','none','xcolor','none','ycolor','none',...
%                 'position',axx.Position,'xlim',axx.XLim);
            h = patch(axx.XLim([1 2 2 1]),axx.YLim([1 1 2 2]),'red');
            hh = hatchfill(h, 'cross', 45, 8);
            set(hh,'color',[0 0 0 0.15],'linewi',1)
%             axes(axx);
%             drawnow;
            
            % AMPLITUDE MASK
            % the contour LINES should be the covarying amplitude, or at
            % least the amplitude of the PWs at KEP:
            tbpamp = TBPamp{ax}; nanlocs = isnan(tbpamp);
            % smooth nans...
            smoo = movmean(tbpamp,7*4,2,'omitnan');
            tbpamp(nanlocs) = smoo(nanlocs);
            % put a few nans back in...
            smoosmoo = movmean(tbpamp,7,2);
            tbpamp(isnan(smoosmoo)) = NaN;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% CONE OF INFLUENCE!
            if coiflag
            % scroll through, and every time you hit a NaN, put a cone of
            % influence mask there:
            startmarkers = datenum(2016,01,01);
            endmarkers = [];
            for i = 2:size(tbpamp,2)
               if ~all(isnan(tbpamp(:,i-1))) && all(isnan(tbpamp(:,i)))
                   startmarkers(end+1) = X(1,i);
               end
               if all(isnan(tbpamp(:,i-1))) && ~all(isnan(tbpamp(:,i)))
                   endmarkers(end+1) = X(1,i);
               end
            end
            endmarkers(end+1) = datenum(2020,12,31);
            % now set these limits to NaN:
            for i = 1:length(startmarkers)
                %             for i = 2
                coi = X > (startmarkers(i) - Y) & X < (endmarkers(i) + Y);
                %                 coi = X > (startmarkers(i) + Y) & X < (endmarkers(i) - Y);
                tbp(coi) = NaN;
                tbpamp(coi) = NaN;
            end
            end
            %             return
            %             coi = X > (axx.XLim(1) + Y) & X < (axx.XLim(2) - Y);
            %             tbp(~coi) = NaN;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % filled contours:
            % need to set regions where the amplitude is low to be white:
            tbp(tbpamp < 6) = 0;
            
            % plot the negs first, then pos to avoid weird contouring colour problems
            % without white in the middle.
            tbppos = tbp; tbpneg = tbp;
            tbppos(tbppos <= 0) = NaN; tbpneg(tbpneg >= 0) = NaN;
            
            hold on; P = pcolor(X,Y,tbp);
            P.EdgeColor = 'none';
            P.FaceColor = 'w';
            P.AlphaData = double(~isnan(tbp));
            P.FaceAlpha = 'flat';
            
            hold on; contourf(X,Y,tbpneg,clevs{ax},'edgecolor','none');
            hold on; contourf(X,Y,tbppos,clevs{ax},'edgecolor','none');
            
            % contour lines:
            hold on; [C,h] = contour(X,Y,tbpamp,blackclines{ax},'edgecolor','w','linewi',2.5);
            hold on; [C,h] = contour(X,Y,tbpamp,blackclines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
            
            %             clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.25),'fontweight','bold');
            %             hold on; [C,h] = contour(X,Y,tbp,whiteclines{ax},'edgecolor',rgbtrip(.95),'linewi',1.25);
            %             clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
            
            logscale('y')
            ylabel('Period (days)')
            
            axx.YGrid = 'on';
            axx.YMinorGrid = 'off';
            
            clim(clims{ax})
            
            %%%% COLORS
            cmap = cmaps{ax};
            colormap(gca,cmap);
%             cmap = [1 1 1 ; cmaps{ax}];
            nclev = 13;
            cmapi = interp1(1:size(cmap,1),cmap,linspace(1,size(cmap,1),nclev),'linear');
            cmapi(ceil(nclev/2),:) = [];
            colormap(gca,cmapi)
            
            %%%% COLORBAR
            cbar = nph_colorbar;
            cbar.Ticks = cbarticks{ax};
            cbar.TickDirection = 'out';
            cbar.Position = cbar.Position .* [1 1 0.4 1];
            cbar.Position = cbar.Position + [-1*cbar.Position(3) 0 0 0];
%             cbar.Label.String = ['\bf{' clabs{ax} '}'];
            cbar.Label.Rotation = 0;
            cbar.Label.VerticalAlignment = 'middle';
%             cbar.Label.HorizontalAlignment = 'left';
%             cbar.Limits(1) = 0;
            cbar.Title.String = '   rad';
            
%             cbar.TickLabelInterpreter = 'latex';
%             cbar.TickLabels = cbarticklabels{ax};     
            cbar.TickLabels = {};     
            
            cbar.Ruler.MinorTick = 'on';
            cbar.Ruler.MinorTickValues = -pi:pi/6:pi;
            
            %%%% MONTHS
            axx.XTick = datenum(years,01,01);
            axx.XMinorTick = 'on';
            axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
            datetick('x',' ','keepticks','keeplimits')
            % other months as ticks using text:
            xtix = datenum(years(1),1:(length(years)*12),15);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,1.25,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
               end
            end
            % Years as lower bold numbers:
            if ax == 3
            xtix = datenum(years,07,01);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   yrstr = datestr(xt,'yyyy');
                   text(xt,0.75,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
               end
            end
            end
            
    end
    
    
    % AXES LETTER
    %letter = alphabet(ax);
    letter = letters{ax};
    switch ax
        case 1
            hold on; nph_text([0.005 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
        case {2,3}
            hold on; nph_text([0.005 0.81],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
            hold on; nph_text([0.005 0.80],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
    end

    setfont(fs)
    
    set(gca,'linewi',1.5,'layer','top','tickdir','out')
    
    %%%% bold line around:
    hold on; plot(axx.XLim([1 2 2 1 1]),axx.YLim([1 1 2 2 1]),'linewi',1.5,'color','k');
    
    
end

return

%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_PW_EWphases'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% TIDAL AMPLITUDES AND S-TRANSFORM COVARIANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([1 0.8])

%-------------------------------------------------------
vert_gap = 0.08;        horz_gap = 0.03;
lower_marg = 0.075;     upper_marg = 0.025;
left_marg = 0.065;      right_marg = 0.075;

rows = 2; cols = 8;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 20;

T = tiledlayout(3,6,'tilespacing','compact');

spans = {...
    [1 6],...
    [1 6],...
    [1 6],...
    [1 6]};

% axs = {...
%     [1  6],...
%     [9 14],...
%     [7 8],...
%     [15 16]};

TBP = {...
    [tidalamp_u ; tidalamp_v],...
    sqrt(abs(Comp.MeanWinds.STu.ST .* conj(Comp.Tides.OneDay.STu.ST))),...
    sqrt(abs(Comp.MeanWinds.STv.ST .* conj(Comp.Tides.OneDay.STv.ST))),...
    []};

clims = {...
    [],...
    [-1 15],...
    [-1 15],...
    [-1 10]};

clevs = {...
    [0:20],...
    [0:20],...
    [0:20],...
    [-20:2.5:20]};

clabs = {...
    'u','v','u','v'};

cbarticks = {...
    0:5:20,...
    0:5:20,...
    0:5:20,...
    0:5:20};

tits = {...
    'TIDAL AMPLITUDE',...
    'CO-VARIANCE, ZONAL',...
    'CO-VARIANCE, MERIDIONAL',...
    ' '};

cmaps = {...
    cbrew('nph_OSOrRdPu'),...
    cbrew('nph_OSBlues'),...
    cbrew('nph_OSOrRdPu'),...
    cbrew('nph_OSOrRdPu')};


% also contour lines:
blackclines = {...
    [0 0],...
    [0 10 20],...
    [0 10 20],...
    [0 0]};
whiteclines = {...
    [0 0],...
    [0 15],...
    [0 15],...
    [0 0]};

letters = {'a','b','c','d'};


for ax = 1:3
    
    nexttile(spans{ax})
    
    axx = gca;
    
    xlim([datenum(2016,02,15) datenum(2020,11,15)])
    
    switch ax
        case 1

            tbp = TBP{ax};
            hold on; plots.v = plot(Comp.dd,tbp(2,:),'color',mcolor(2),'linewi',1.5);
            hold on; plots.u = plot(Comp.dd,tbp(1,:),'color',mcolor(1),'linewi',1.5);
            
            grid on;
            
            ylim([0 100])
            ytick(-100:20:100);
            yminortick(-100:10:100);
            
            ylabel('ms^-^1')
            
            L = legend([plots.u plots.v],{'Zonal Component','Meridional Component'},'fontsize',0.8*fs,'AutoUpdate','off');
            
            %%%% MONTHS
            axx.XTick = datenum(years,01,01);
            axx.XMinorTick = 'on';
            axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
            datetick('x',' ','keepticks','keeplimits')
            % other months as ticks using text:
            xtix = datenum(years(1),1:(length(years)*12),15);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,axx.YLim(1)-5,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
               end
            end
%             % Years as lower bold numbers:
%             xtix = datenum(years,07,01);
%             for xt = xtix
%                if inrange(xt,xlim) && ~any(xt == axx.XTick)
%                    yrstr = datestr(xt,'yyyy');
%                    text(xt,axx.YLim(1)-25,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
%                end
%             end
        
        % allll the contours...
        case {2,3}
            
            tbp = TBP{ax}; nanlocs = isnan(tbp);
            % smooth nans...
            smoo = movmean(tbp,7,2,'omitnan');
            tbp(nanlocs) = smoo(nanlocs);
            % put a few nans back in...
            smoosmoo = movmean(tbp,7,2);
            tbp(isnan(smoosmoo)) = NaN;
            
            % grids...
            [X,Y] = meshgrid(Comp.dd,1./Comp.MeanWinds.STu.freqs);
            smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
            
            ylim([1.4 35])
            axx.YTick = [0 2 5 10 20 30 40 50 100];
            axx.YAxis.MinorTickValues = [0:1:20 25:5:100];
            
            % hatched area below:
%             axx2 = axes;
%             axes(axx2);
%             set(axx2,'color','none','xcolor','none','ycolor','none',...
%                 'position',axx.Position,'xlim',axx.XLim);
            h = patch(axx.XLim([1 2 2 1]),axx.YLim([1 1 2 2]),'red');
            hh = hatchfill(h, 'cross', 45, 8);
            set(hh,'color',[0 0 0 0.15],'linewi',1)
%             axes(axx);
%             drawnow;
            
            %%%% CONE OF INFLUENCE!
            coi = X > (axx.XLim(1) + Y) & X < (axx.XLim(2) - Y);
            tbp(~coi) = NaN;
            
            % filled contours:
            hold on; contourf(X,Y,tbp,clevs{ax},'edgecolor','none');
            
            % contour lines:
            hold on; [C,h] = contour(X,Y,tbp,blackclines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
            clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.25),'fontweight','bold');
            hold on; [C,h] = contour(X,Y,tbp,whiteclines{ax},'edgecolor',rgbtrip(.95),'linewi',1.25);
            clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
            
            logscale('y')
            ylabel('Period (days)')
            
            axx.YGrid = 'on';
            axx.YMinorGrid = 'off';
            
            clim(clims{ax})
            
            %%%% COLORS
            cmap = [1 1 1 ; cmaps{ax}];
            nclev = 18;
            cmapi = interp1(1:size(cmap,1),cmap,linspace(1,size(cmap,1),nclev),'linear');
            colormap(gca,cmapi)
            
            %%%% COLORBAR
            cbar = nph_colorbar;
            cbar.Ticks = cbarticks{ax};
            cbar.Position = cbar.Position .* [1 1 0.4 1];
            cbar.Position = cbar.Position + [-1*cbar.Position(3) 0 0 0];
%             cbar.Label.String = ['\bf{' clabs{ax} '}'];
            cbar.Label.Rotation = 0;
            cbar.Label.VerticalAlignment = 'middle';
%             cbar.Label.HorizontalAlignment = 'left';
            cbar.Limits(1) = 0;
            cbar.Title.String = '   ms^{-1}'; 
            
            %%%% MONTHS
            axx.XTick = datenum(years,01,01);
            axx.XMinorTick = 'on';
            axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
            datetick('x',' ','keepticks','keeplimits')
            % other months as ticks using text:
            xtix = datenum(years(1),1:(length(years)*12),15);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,1.25,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
               end
            end
            % Years as lower bold numbers:
            if ax == 3
            xtix = datenum(years,07,01);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   yrstr = datestr(xt,'yyyy');
                   text(xt,0.75,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
               end
            end
            end
            
    end
    
    
    % AXES LETTER
    %letter = alphabet(ax);
    letter = letters{ax};
    switch ax
        case 1
            hold on; nph_text([0.005 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
            hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        case {2,3}
            hold on; nph_text([0.005 0.81],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
            hold on; nph_text([0.005 0.80],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
    end
    
    setfont(fs)
    
    set(gca,'linewi',1.5,'layer','top','tickdir','out')
    
    
end

return

%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_TidalAmplitudes'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return


















return



%% CASE STUDY TIMESERIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold all; whitefig; figpos([1 0.8])

%-------------------------------------------------------
vert_gap = 0.08;        horz_gap = 0.03;
lower_marg = 0.075;     upper_marg = 0.025;
left_marg = 0.065;      right_marg = 0.075;

rows = 2; cols = 8;

subplot = @(rows,cols,p) subtightplot (rows,cols,p,[vert_gap horz_gap],[lower_marg upper_marg],[left_marg right_marg]);

%--------------------------------------------------------

fs = 20;

T = tiledlayout(3,2,'tilespacing','compact');

spans = {...
    [1 1],...
    [1 1],...
    [1 1],...
    [1 1]};

xlims = {...
    [datenum(2017,05,01) datenum(2017,10,05)],...
    [datenum(2018,12,25) datenum(2019,02,05)],...
    [datenum(2016,02,15) datenum(2020,11,15)],...
    [datenum(2016,02,15) datenum(2020,11,15)]};
% axs = {...
%     [1  6],...
%     [9 14],...
%     [7 8],...
%     [15 16]};

% TBP = {...
%     [tidalamp_u ; tidalamp_v],...
%     sqrt(abs(Comp.MeanWinds.STu.ST .* Comp.Tides.OneDay.STu.ST)),...
%     sqrt(abs(Comp.MeanWinds.STv.ST .* Comp.Tides.OneDay.STv.ST)),...
%     []};

clims = {...
    [],...
    [-1 15],...
    [-1 15],...
    [-1 10]};

clevs = {...
    [0:20],...
    [0:20],...
    [0:20],...
    [-20:2.5:20]};

clabs = {...
    'u','v','u','v'};

cbarticks = {...
    0:5:20,...
    0:5:20,...
    0:5:20,...
    0:5:20};

tits = {...
    'TIDAL AMPLITUDE',...
    'CO-VARIANCE, ZONAL',...
    'CO-VARIANCE, MERIDIONAL',...
    ' '};

cmaps = {...
    cbrew('nph_OSOrRdPu'),...
    cbrew('nph_OSBlues'),...
    cbrew('nph_OSOrRdPu'),...
    cbrew('nph_OSOrRdPu')};


% also contour lines:
blackclines = {...
    [0 0],...
    [0 10 20],...
    [0 10 20],...
    [0 0]};
whiteclines = {...
    [0 0],...
    [0 15],...
    [0 15],...
    [0 0]};

letters = {'a','b','c','d'};


for ax = 1:2
    
    nexttile(spans{ax})
    
    axx = gca;
    
%     switch ax
%         case 1
% 
%             tbp = TBP{ax};

            hold on; plots.tidev = plot(Comp.dd,tidalamp_v-movmean(tidalamp_v,31,2,'omitnan'),'color',mcolor(2),'linewi',1.5,'linest','-');
%             hold on; plots.tideu = plot(Comp.dd,tidalamp_u-movmean(tidalamp_u,31,2,'omitnan'),'color',mcolor(1),'linewi',1.5,'linest','-');
            
            hold on; plots.windv = plot(Comp.dd,meanwinds_v-movmean(meanwinds_v,31,2,'omitnan'),'color',mcolor(2),'linewi',2.5,'linest',':');
%             hold on; plots.windu = plot(Comp.dd,meanwinds_u-movmean(meanwinds_u,31,2,'omitnan'),'color',mcolor(1),'linewi',2.5,'linest',':');
            
            
            xlim(xlims{ax})
            
            grid on;
            
            ylim([0 100])
            ytick(-100:20:100);
            yminortick(-100:10:100);
            
            ylabel('ms^-^1')
            
%             L = legend([plots.u plots.v],{'Zonal Component','Meridional Component'},'fontsize',0.8*fs,'AutoUpdate','off');
            
            %%%% MONTHS
            axx.XTick = datenum(years,01,01);
            axx.XMinorTick = 'on';
            axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
            datetick('x',' ','keepticks','keeplimits')
            % other months as ticks using text:
            xtix = datenum(years(1),1:(length(years)*12),15);
            for xt = xtix
               if inrange(xt,xlim) && ~any(xt == axx.XTick)
                   mn = monthname(month(xt));
                   text(xt,axx.YLim(1)-5,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
               end
            end
%             % Years as lower bold numbers:
%             xtix = datenum(years,07,01);
%             for xt = xtix
%                if inrange(xt,xlim) && ~any(xt == axx.XTick)
%                    yrstr = datestr(xt,'yyyy');
%                    text(xt,axx.YLim(1)-25,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
%                end
%             end
        
%         % allll the contours...
%         case {2,3}
%             
%             tbp = TBP{ax}; nanlocs = isnan(tbp);
%             % smooth nans...
%             smoo = movmean(tbp,7,2,'omitnan');
%             tbp(nanlocs) = smoo(nanlocs);
%             % put a few nans back in...
%             smoosmoo = movmean(tbp,7,2);
%             tbp(isnan(smoosmoo)) = NaN;
%             
%             % grids...
%             [X,Y] = meshgrid(Comp.dd,1./Comp.MeanWinds.STu.freqs);
%             smoo1 = [0.25 0.25]; smoo2 = [0.25 0.25];
%             
%             ylim([1.4 35])
%             axx.YTick = [0 2 5 10 20 30 40 50 100];
%             axx.YAxis.MinorTickValues = [0:1:20 25:5:100];
%             
%             % hatched area below:
% %             axx2 = axes;
% %             axes(axx2);
% %             set(axx2,'color','none','xcolor','none','ycolor','none',...
% %                 'position',axx.Position,'xlim',axx.XLim);
%             h = patch(axx.XLim([1 2 2 1]),axx.YLim([1 1 2 2]),'red');
%             hh = hatchfill(h, 'cross', 45, 8);
%             set(hh,'color',[0 0 0 0.15],'linewi',1)
% %             axes(axx);
% %             drawnow;
%             
%             %%%% CONE OF INFLUENCE!
%             coi = X > (axx.XLim(1) + Y) & X < (axx.XLim(2) - Y);
%             tbp(~coi) = NaN;
%             
%             % filled contours:
%             hold on; contourf(X,Y,tbp,clevs{ax},'edgecolor','none');
%             
%             % contour lines:
%             hold on; [C,h] = contour(X,Y,tbp,blackclines{ax},'edgecolor',rgbtrip(.25),'linewi',1.25);
%             clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.25),'fontweight','bold');
%             hold on; [C,h] = contour(X,Y,tbp,whiteclines{ax},'edgecolor',rgbtrip(.95),'linewi',1.25);
%             clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
%             
%             logscale('y')
%             ylabel('Period (days)')
%             
%             axx.YGrid = 'on';
%             axx.YMinorGrid = 'off';
%             
%             clim(clims{ax})
%             
%             %%%% COLORS
%             cmap = [1 1 1 ; cmaps{ax}];
%             nclev = 18;
%             cmapi = interp1(1:size(cmap,1),cmap,linspace(1,size(cmap,1),nclev),'linear');
%             colormap(gca,cmapi)
%             
%             %%%% COLORBAR
%             cbar = nph_colorbar;
%             cbar.Ticks = cbarticks{ax};
%             cbar.Position = cbar.Position .* [1 1 0.4 1];
%             cbar.Position = cbar.Position + [-1*cbar.Position(3) 0 0 0];
% %             cbar.Label.String = ['\bf{' clabs{ax} '}'];
%             cbar.Label.Rotation = 0;
%             cbar.Label.VerticalAlignment = 'middle';
% %             cbar.Label.HorizontalAlignment = 'left';
%             cbar.Limits(1) = 0;
%             cbar.Title.String = '   ms^{-1}'; 
%             
%             %%%% MONTHS
%             axx.XTick = datenum(years,01,01);
%             axx.XMinorTick = 'on';
%             axx.XAxis.MinorTickValues = datenum(min(years),1:(length(years)*12),01);
%             datetick('x',' ','keepticks','keeplimits')
%             % other months as ticks using text:
%             xtix = datenum(years(1),1:(length(years)*12),15);
%             for xt = xtix
%                if inrange(xt,xlim) && ~any(xt == axx.XTick)
%                    mn = monthname(month(xt));
%                    text(xt,1.25,mn(1),'fontsize',0.9*fs,'horizontalalignment','center','VerticalAlignment','top');
%                end
%             end
%             % Years as lower bold numbers:
%             if ax == 3
%             xtix = datenum(years,07,01);
%             for xt = xtix
%                if inrange(xt,xlim) && ~any(xt == axx.XTick)
%                    yrstr = datestr(xt,'yyyy');
%                    text(xt,0.75,yrstr,'fontsize',1*fs,'fontweight','bold','horizontalalignment','center','VerticalAlignment','top');
%                end
%             end
%             end
%             
%     end
    
    
    % AXES LETTER
    %letter = alphabet(ax);
    letter = letters{ax};
    switch ax
        case 1
            hold on; nph_text([0.005 0.85],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
            hold on; nph_text([0.005 0.84],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
        case {2,3}
            hold on; nph_text([0.005 0.81],['(' letter ')'],'fontsize',1.5*fs,'textborder','w','horizontalalignment','left','verticalalignment','middle');
            hold on; nph_text([0.005 0.80],['       ' tits{ax}],'fontsize',1*fs,'fontweight','bold','textborder','w','horizontalalignment','left','verticalalignment','middle');
    end
    
    setfont(fs)
    
    set(gca,'linewi',1.5,'layer','top','tickdir','out')
    
    
end

return

%% EXPORT??? ==============================================================

savename = ['~/Desktop/' upper(site) '_TidalAmplitudes'];

disp(['Saving as "' savename '"...'])

nph_saveas(gcf,savename,'png')

disp('Saved.')







return


















return




























%% PLOT OF TIDAL AMPLITUDES/PHASES AGAINST ALL TIME

for td = 1:4
    
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
    ap = 'a'; % a | p
    
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
            ampclims = [0 5];
            ampclevs = 0:0.5:20;
            ampclines = 0:10;
            ampcbarticks = 0:1:10;
            ampcbarminorticks = 0:10;
            phaseminorclines = 0:1:tdh;
            darkclines = [4:10];
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
                nanmask = abs(diff(tbp2([1 1:end],:),[],1)) > 6 | abs(diff(tbp2(:,[1 1:end]),[],2)) > 6;
                tbp2(nanmask) = NaN;
                hold on; [C,h] = contour(X,Y,tbp2,[0 0],'edgecolor',rgbtrip(.25),'linewi',1.25);
                clabel(C,h,'fontsize',0.6*fs,'labelspacing',2000,'color',rgbtrip(.95),'fontweight','bold');
                
                % before you do the others, exclude any wraparound bits:
                nanmask = abs(diff(tbp([1 1:end],:),[],1)) > 6 | abs(diff(tbp(:,[1 1:end]),[],2)) > 6;
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
    
    
    % return
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































