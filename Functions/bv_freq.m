

% Estimate the Brunt-Vaisala Frequency at a given altitide and temperature.
% Using the wiki page for standard atmosphere.

% Input:  Altitude in KM.
% Output: BV freq in RAD/S.

% EDIT Jan 2019: Updated to be array safe!!

% EDIT EDIT: Okay, ffs, now gonna use the lapse rate formulation to compute
% BV freq, cos it's much easier.


% z = 0:100;

function N = bv_freq(z,varargin)

ref_z = [...
    0
    11
    20
    32
    47
    51
    71];

ref_temp = [...
    288.15
    216.65
    216.65
    228.65
    270.65
    270.65
    214.65];

ref_pres = [...
    1013.25
    226.321
    54.7489
    8.68019
    1.10906
    0.669389
    0.0395642];

% First, get a spline-interpolated curve for the US standard atmosphere.
% Because the temp doesn't change between a couple of the layers, the lapse
% rate is zero and the BV period is infinite. Not a problem, but you don't
% want zeros creeping in.
F_temp = griddedInterpolant(ref_z,ref_temp,'spline');

% Now compute the lapse rate for this curve at a high resolution:
zi = linspace(0,100,1000); % some height range between 0 and 100km.
gamma = - diff(F_temp(zi));
gamma = gamma([1 1:end]); % add the first

% Use this an an interpolant for our BV computation:
F_gamma = griddedInterpolant(zi,gamma,'spline');

% And compute BV freq:
g = 9.8; g_km = g ./ 1000;
gamma_d = 9.8; % K/km

N = sqrt((g_km./F_temp(z)) .* (gamma_d - F_gamma(z)));


if ~isempty(varargin)
    if any(strcmp(varargin,'plot'))
        zz = 0:100;
        NN = sqrt((g_km./F_temp(zz)) .* (gamma_d - F_gamma(zz)));
        
        figure; hold all;
        
        subplot(1,2,1);
        hold on; plot(NN,zz); grid on;
        hold on; plot(N,z,'marker','+','linest','none'); grid on;
        xlabel('BV Frequency (rad/s)')
        ylabel('Altitude (km)')
        legend('Standard Atmosphere','Input Altitudes')
        set(gca,'fontsize',12)
        
        subplot(1,2,2);
        hold on; plot(1./(NN./(2*pi)),zz); grid on;
        hold on; plot(1./(N./(2*pi)),z,'marker','+','linest','none'); grid on;
        xlabel('BV Period (seconds)')
        ylabel('Altitude (km)')
        legend('Standard Atmosphere','Input Altitudes')
        set(gca,'fontsize',12)
        
    end
end

return

% % % % % 
% % % % % % ALSO: Turns out that potential temperature increases with height quite
% % % % % % smoothly, even in the standard atmosphere. This means that it's very good
% % % % % % for creating an interpolant using a spline, which is much better for
% % % % % % getting smooth profiles of BV freq. You can then take nice smooth values
% % % % % % for rate of change of potential temp with height.
% % % % % 
% % % % % % Only works between 11km and 71km, currently. You get some weird numbers
% % % % % % low down in the bottom layer due to the spline interpolant.
% % % % % 
% % % % % % find reference potential temperature: (wrt to ground)
% % % % % R = 8.31; Cp = 29.07;
% % % % % ref_ptemp = ref_temp .* (ref_pres(1) ./ ref_pres) .^ (R./Cp);
% % % % % 
% % % % % % get interpolants for potential temperature and rate of change of
% % % % % % potential temperature:
% % % % % F_ptemp = griddedInterpolant(ref_z,ref_ptemp,'linear');
% % % % % zz = linspace(ref_z(1),ref_z(end),100);
% % % % % ref_ptempi = F_ptemp(zz);
% % % % % 
% % % % % dptempdz = abs(diff(ref_ptempi) ./ diff(zz));
% % % % % dptempdz = dptempdz([1 1:end]);
% % % % % F_dptempdz = griddedInterpolant(zz,dptempdz,'linear','none');
% % % % % 
% % % % % % finally compute N by evaluating the potential and rate of change of
% % % % % % potential temperatures, and computing:
% % % % % g = 9.8;
% % % % % 
% % % % % N = sqrt( (g./F_ptemp(z)) .* (F_dptempdz(z)) );
% % % % % 
% % % % % return


% % % % % 
% % % % % 
% % % % % % make an over-interpolated version for the rate of change, it's just
% % % % % % weirdly cleaner I think.
% % % % % zz = linspace(ref_z(1),ref_z(end),100);
% % % % % ref_ptempi = interp1(ref_z,ref_ptemp,zz,'linear');
% % % % % dptempdz = diff(ref_ptempi) ./ diff(zz);
% % % % % dptempdz = dptempdz([1 1:end]);
% % % % % 
% % % % % % then use these to build the interpolant (reduces artifacts)
% % % % % % midpoints = ref_z(1:6) + (diff(ref_z)./2); % use nearest neighbour from midpoints
% % % % % % dptempdz = diff(ref_ptemp)./diff(ref_z);
% % % % % % dptempdz = [dptempdz(1) ; dptempdz ; dptempdz(end)]; % expand to get rate of change for end points
% % % % % % F_dptempdz = griddedInterpolant([ref_z(1) ; midpoints ; ref_z(end)],dptempdz,'linear','none');
% % % % % F_dptempdz = griddedInterpolant(zz,dptempdz,'linear','none');
% % % % % 
% % % % % 
% % % % % % finally compute N by evaluating the potential and rate of change of
% % % % % % potential temperatures, and computing:
% % % % % g = 9.8;
% % % % % 
% % % % % N = sqrt( (g./F_ptemp(z)) .* (F_dptempdz(z)) );


% % % % % 
% % % % % return
% % % % % % find rate of change of density with altitude:
% % % % % 
% % % % % diff_z = diff(ref_z);
% % % % % diff_rho = diff(ref_rho);
% % % % % 
% % % % % drho_dz = diff_rho ./ diff_z;
% % % % % drho_dz(end+1) = drho_dz(end); % extend last layer to be same as previous above 71km.
% % % % % 
% % % % % % Interpolate density values and evaulate at input altitude:
% % % % % F_rho = griddedInterpolant(ref_z,ref_rho,'linear','none');
% % % % % rho_z = F_rho(z);
% % % % % 
% % % % % midpoints = ref_z(1:6) + (diff_z./2);
% % % % % F_drho_dz = griddedInterpolant(midpoints,diff_rho,'nearest');
% % % % % drho_z = F_drho_dz(z);
% % % % % 
% % % % % % Assume rate of change of density is constant for each layer,
% % % % % layers = [
% % % % %     0  11 20 32 47 51 71 ;
% % % % %     11 20 32 47 51 71 Inf] .* 1000;
% % % % % 
% % % % % layer = z >= layers(1,:) & z < layers(2,:);
% % % % % 
% % % % % g = 9.81;
% % % % % 
% % % % % 
% % % % % 
% % % % % N = sqrt( - (g./rho_z) .* drho_dz(layer));
% % % % % 
% % % % % z
% % % % % 
% % % % % rho_z
% % % % % 
% % % % % drho_dz(layer)
% % % % % 
% % % % % N
% % % % % 
% % % % % 
% % % % % return
% % % % % 
% % % % % 
% % % % % 
% % % % % % now find rho for the input altitude, then you will have all you need to
% % % % % % compute local Brunt-Vaisala frequency:
% % % % % g = 9.81;
% % % % % 
% % % % % N = NaN;
% % % % % 
% % % % % if z >= ref_z(1) && z < ref_z(2)
% % % % %     layer = 1;
% % % % %     rho_at_z = ref_rho(layer) + drho_dz(layer) .* (z - ref_z(layer));
% % % % %     N = sqrt( - (g./rho_at_z) .* drho_dz(layer));
% % % % % end
% % % % % if z >= ref_z(2) && z < ref_z(3)
% % % % %     layer = 2;
% % % % %     rho_at_z = ref_rho(layer) + drho_dz(layer) .* (z - ref_z(layer));
% % % % %     N = sqrt( - (g./rho_at_z) .* drho_dz(layer));
% % % % % end
% % % % % if z >= ref_z(3) && z < ref_z(4)
% % % % %     layer = 3;
% % % % %     rho_at_z = ref_rho(layer) + drho_dz(layer) .* (z - ref_z(layer));
% % % % %     N = sqrt( - (g./rho_at_z) .* drho_dz(layer));
% % % % % end
% % % % % if z >= ref_z(4) && z < ref_z(5)
% % % % %     layer = 4;
% % % % %     rho_at_z = ref_rho(layer) + drho_dz(layer) .* (z - ref_z(layer));
% % % % %     N = sqrt( - (g./rho_at_z) .* drho_dz(layer));
% % % % % end
% % % % % if z >= ref_z(5) && z < ref_z(6)
% % % % %     layer = 5;
% % % % %     rho_at_z = ref_rho(layer) + drho_dz(layer) .* (z - ref_z(layer));
% % % % %     N = sqrt( - (g./rho_at_z) .* drho_dz(layer));
% % % % % end
% % % % % if z >= ref_z(6) && z < ref_z(7)
% % % % %     layer = 6;
% % % % %     rho_at_z = ref_rho(layer) + drho_dz(layer) .* (z - ref_z(layer));
% % % % %     N = sqrt( - (g./rho_at_z) .* drho_dz(layer));
% % % % % end
% % % % % if z >= ref_z(7) && z < ref_z(8)
% % % % %     layer = 7;
% % % % %     rho_at_z = ref_rho(layer) + drho_dz(layer) .* (z - ref_z(layer));
% % % % %     N = sqrt( - (g./rho_at_z) .* drho_dz(layer));
% % % % % end
% % % % % if z >= ref_z(7) % set hard limit to density beyond 100km.
% % % % % %     layer = 7;
% % % % % %     rho_at_z = ref_rho(layer) + drho_dz(layer) .* (z - ref_z(layer));
% % % % % %     if rho_at_z < 0, rho_at_z = ref_rho(layer); end
% % % % % %     N = sqrt( - (g./rho_at_z) .* drho_dz(layer));
% % % % % N = NaN;
% % % % % end
% % % % % 
% % % % % 
% % % % % return
% % % % % 
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % 
% % % % % 

% % % % % 
% % % % % 
% % % % % 
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % OLD METHOD, USES POTENTIAL TEMPERATURE:
% % % % % 
% % % % % % zero altitudes or z > 71km don't work yet...
% % % % % 
% % % % % z(z == 0) = 1;
% % % % % z(z > 71) = 71;
% % % % % 
% % % % % z = round(z);
% % % % % 
% % % % % 
% % % % % % Make an approximate, smoothed temperature profile for the atmosphere:
% % % % % 
% % % % % R = 8.31; Cp = 29.07;
% % % % % 
% % % % % % Use standard atmosphere, spline interpolated:
% % % % % 
% % % % % %t = [288.15 216.65 216.65 228.65 270.65 270.65 214.65];
% % % % % %t = interp1([0 11 20 32 47 51 71],t,0:71,'spline');
% % % % % 
% % % % % t = [288.15000000000,275.970181330425,265.156519742410,255.632284120071, ...
% % % % %     247.320743347529,240.145166308902,234.028821888309,228.894978969868, ...
% % % % %     224.666906437697,221.267873175917,218.621148068645,216.650000000000, ...
% % % % %     215.277697854101,214.427510515066,214.022706867015,213.986555794065, ...
% % % % %     214.242326180336,214.713286909947,215.322706867015,215.993854935660, ...
% % % % %     216.650000000000,217.233639312062,217.764183595503,218.280271941887, ...
% % % % %     218.820543442780,219.423637189746,220.128192274350,220.972847788157, ...
% % % % %     221.996242822732,223.237016469638,224.733807820442,226.525255966708, ...
% % % % %     228.650000000000,231.130044567674,233.920856540249,236.961268344032, ...
% % % % %     240.190112405332,243.546221150458,246.968427005717,250.395562397420, ...
% % % % %     253.766459751873,257.019951495387,260.094870054269,262.930047854827, ...
% % % % %     265.464317323371,267.636510886209,269.385460969649,270.650000000000, ...
% % % % %     271.384730807997,271.607337842084,271.351275955128,270.650000000000, ...
% % % % %     269.536964829569,268.045625296704,266.209436254274,264.061852555150, ...
% % % % %     261.636329052199,258.966320598292,256.085282046298,253.026668249087, ...
% % % % %     249.823934059527,246.510534330487,243.119923914838,239.685557665449, ...
% % % % %     236.240890435188,232.819377076926,229.454472443532,226.179631387874, ...
% % % % %     223.028308762823,220.033959421247,217.230038216016,214.650000000000];
% % % % % 
% % % % % pt = t .* (1013.25./alt2pres(0:71)) .^(R/Cp);
% % % % % 
% % % % % dpt = abs(diff(pt)); dpt(end+1) = dpt(end);
% % % % % 
% % % % % %===================================
% % % % % 
% % % % % % N = zeros(size(z));
% % % % % 
% % % % % g = 9.8 ./ 1000; %ms-2 Could reduce g with height to get a slightly more accurate answer...
% % % % % 
% % % % % 
% % % % % N = (((g./pt(z)).*dpt(z)).^0.5)';
% % % % % 
% % % % % 
% % % % % % 
% % % % % % for i = 1:length(z),
% % % % % % 
% % % % % %    ind = (roundto(z(i),0)+1);
% % % % % %     
% % % % % %     
% % % % % % 
% % % % % % 
% % % % % 
% % % % % % end % end cycling through all of z's element(s)...
% % % % % 
% % % % % 
% % % % % 

















