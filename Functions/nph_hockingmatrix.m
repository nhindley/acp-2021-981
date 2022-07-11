
% Application of Hocking 2005.
%
% Compute u'^2, v'^2, w'^2, u'v', u'w', v'w' from radial velocity (minus
% mean winds), zenith angle (theta) and angle anticlockwise from the east (phi)
%
% THETA     = ZENITH ANGLE (deg)
% PHI       = ANGLE ANITCLOCKWISE FROM EAST (deg).    
%
% inputs in M/S and DEGREES.
%
% EDIT: Need to be able to cope with NaNs in the input. Obviously not all
% NaNs, but some.
%
% EDIT: re-typing of robins code. Let's be more careful, and arrange the
% matrices EXACTLY as they are given in Hocking 2005, except obv matlab
% counts from the bottom left:
% 
% % | a(6,1)                             ...    |  | b(6,1) |     | c(6,1) |
% % | a(5,1)                      ...           |  | b(5,1) |     | c(5,1) |
% % | a(4,1)               ...                  |  | b(4,1) |     | c(4,1) |
% % | a(3,1)        ...                         |  | b(3,1) |  =  | c(3,1) |
% % | a(2,1) ...                                |  | b(2,1) |     | c(2,1) |
% % | a(1,1) a(2,1) a(3,1) a(4,1) a(5,1) a(6,1) |  | b(1,1) |     | c(1,1) |
% %
% %


function OUT = nph_hockingmatrix(vradprime,theta,phi)

% A * B = C, where B = [up_sq, vp_sq, wp_sq, upvp, upwp, vpwp]

% the matrix from Hocking 05
vp2 = vradprime.^2;

st = sind(theta);   s2t = sind(theta).^2;   s3t = sind(theta).^3;   s4t = sind(theta).^4;
ct = cosd(theta);   c2t = cosd(theta).^2;   c3t = cosd(theta).^3;   c4t = cosd(theta).^4;

sp = sind(phi);     s2p = sind(phi).^2;     s3p = sind(phi).^3;     s4p = sind(phi).^4;
cp = cosd(phi);     c2p = cosd(phi).^2;     c3p = cosd(phi).^3;     c4p = cosd(phi).^4;

A = [...
    nansum(s4t.*c4p)            nansum(s4t.*c2p.*s2p)       nansum(s2t.*c2t.*c2p)       nansum(2*s4t.*c3p.*sp)      nansum(2*s3t.*ct.*c3p)      nansum(2*s3t.*ct.*c2p.*sp)  ; ...
    nansum(s4t.*c2p.*s2p)       nansum(s4t.*s4p)            nansum(s2t.*c2t.*s2p)       nansum(2*s4t.*s3p.*cp)      nansum(2*s3t.*ct.*s2p.*cp)  nansum(2*s3t.*ct.*s3p)      ; ...
    nansum(s2t.*c2t.*c2p)       nansum(s2t.*c2t.*s2p)       nansum(c4t)                 nansum(2*s2t.*c2t.*cp.*sp)  nansum(2*c3t.*st.*cp)       nansum(2*c3t.*st.*sp)       ; ...
    nansum(2*s4t.*c3p.*sp)      nansum(2*s4t.*s3p.*cp)      nansum(2*s2t.*c2t.*cp.*sp)  nansum(4*s4t.*c2p.*s2p)     nansum(4*s3t.*ct.*c2p.*sp)  nansum(4*s3t.*ct.*s2p.*cp)  ; ...
    nansum(2*s3t.*ct.*c3p)      nansum(2*s3t.*ct.*s2p.*cp)  nansum(2*c3t.*st.*cp)       nansum(4*s3t.*ct.*c2p.*sp)  nansum(4*s2t.*c2t.*c2p)     nansum(4*s2t.*c2t.*cp.*sp)  ; ...
    nansum(2*s3t.*ct.*c2p.*sp)  nansum(2*s3t.*ct.*s3p)      nansum(2*c3t.*st.*sp)       nansum(4*s3t.*ct.*s2p.*cp)  nansum(4*s2t.*c2t.*cp.*sp)  nansum(4*s2t.*c2t.*s2p)     ];

C = [...
    nansum(vp2.*s2t.*c2p)        ; ...
    nansum(vp2.*s2t.*s2p)        ; ...
    nansum(vp2.*c2t)             ; ...
    nansum(2*vp2.*s2t.*cp.*sp)   ; ...
    nansum(2*vp2.*st.*ct.*cp)    ; ...
    nansum(2*vp2.*st.*ct.*sp)    ];

% % % % % % top left quarter:
% % % % % A = [...
% % % % %     nansum(s4t.*c4p)            nansum(s4t.*c2p.*s2p)       nansum(s2t.*c2t.*c2p)  ; ...
% % % % %     nansum(s4t.*c2p.*s2p)       nansum(s4t.*s4p)            nansum(s2t.*c2t.*s2p) ; ...
% % % % %     nansum(s2t.*c2t.*c2p)       nansum(s2t.*c2t.*s2p)       nansum(c4t)];
% % % % % 
% % % % % C = [...
% % % % %     nansum(2*vp2 .*s2t .*c2p)   ; ...
% % % % %     nansum(2*vp2 .*s2t .*s2p)    ; ...
% % % % %     nansum(2*vp2 .*c2t      )    ];


OUT = A \ C;

% OUT = lsqlin(A,C,[],[],[],[],[0 0 0 0 0 0]);










return


% ROBIN'S TYPING UP OF HOCKING MATRIX: (thanks robin! saved me some time!)

% the matrix from Hocking 05
% st = sind(theta);
% ct = cosd(theta);
% sp = sind(phi);
% cp = cosd(phi);
% s2t = (sind(theta)).^2;
% c2t = (cosd(theta)).^2;
% s2p = (sind(phi)).^2;
% c2p = (cosd(phi)).^2;
% s3t = (sind(theta)).^3;
% c3t = (cosd(theta)).^3;
% s3p = (sind(phi)).^3;
% c3p = (cosd(phi)).^3;
% s4t = (sind(theta)).^4;
% c4t = (cosd(theta)).^4;
% s4p = (sind(phi)).^4;
% c4p = (cosd(phi)).^4;
% v2  = vradprime.^2;
% L1 = [nansum(s4t.*c4p) nansum(s4t.*c2p.*s2p) nansum(s2t.*c2t.*c2p) nansum(2*s4t.*c3p.*sp) nansum(2*s3t.*ct.*c3p) nansum(2*s3t.*ct.*c2p.*sp)];
% L2 = [nansum(s4t.*c2p.*s2p) nansum(s4t.*s4p) nansum(s2t.*c2t.*s2p) nansum(2*s4t.*s3p.*cp) nansum(2*s3t.*ct.*s2p.*cp) nansum(2*s3t.*ct.*s3p)];
% L3 = [nansum(s2t.*c2t.*c2p) nansum(s2t.*c2t.*s2p) nansum(c4t) nansum(2*s2t.*c2t.*cp.*sp) nansum(2*c3t.*st.*cp) nansum(2*c3t.*st.*sp)];
% L4 = [nansum(2*s4t.*c3p.*sp) nansum(2*s4t.*s3p.*cp) nansum(2*s2t.*c2t.*cp.*sp) nansum(4*s4t.*c2p.*s2p) nansum(4*s3t.*ct.*c2p.*sp) nansum(4*s3t.*ct.*s2p.*cp)];
% L5 = [nansum(2*s3t.*ct.*c3p) nansum(2*s3t.*ct.*s2p.*cp) nansum(2*c3t.*st.*cp) nansum(4*s3t.*ct.*c2p.*sp) nansum(4*s2t.*c2t.*c2p) nansum(4*s2t.*c2t.*cp.*sp)];
% L6 = [nansum(2*s3t.*ct.*c2p.*sp) nansum(2*s3t.*ct.*s3p) nansum(2*c3t.*st.*sp) nansum(4*s3t.*ct.*s2p.*cp) nansum(4*s2t.*c2t.*cp.*sp) nansum(4*s2t.*c2t.*s2p)];
% 
% A = [L1;L2;L3;L4;L5;L6];
% 
% b = [...
%     nansum(v2.*s2t.*c2p); ...
%     nansum(v2.*s2t.*s2p); ...
%     nansum(v2.*c2t); ...
%     nansum(2*v2.*s2t.*cp.*sp); ...
%     nansum(2*v2.*st.*ct.*cp); ...
%     nansum(2*v2.*st.*ct.*sp)];
% 
% % OUT = [up_sq,vp_sq,wp_sq,upvp,upwp,vpwp];
% OUT = A\b;
              
end





























