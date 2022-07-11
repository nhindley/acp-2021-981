

% Designed to convert MLS Temperature data from a .he5 file to a .mat file in
% a nice structure.
%
% Note that the documnetation for how the he5 file is laid out is sketchy,
% so I had to find each variable myself - so I'm not sure how robust this
% will be in the future.
%
%


function MLS = nph_mls2mat(filename)

MLS = struct;

% Read attributes for filename and date:

fn = strsplit(filename,'/');

dy = h5readatt(filename,'/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/','GranuleDay');
mth = h5readatt(filename,'/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/','GranuleMonth');
yr = h5readatt(filename,'/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/','GranuleYear');

MLS.Date = datestr(datenum(double(yr),double(mth),double(dy)));
MLS.Filename = fn{end};
MLS.Meta = struct;
MLS.Data = struct;

% Read data fields:
MLS.Data.Latitude   = double(h5read(filename,'/HDFEOS/SWATHS/Temperature/Geolocation Fields/Latitude'));
MLS.Data.Longitude  = double(h5read(filename,'/HDFEOS/SWATHS/Temperature/Geolocation Fields/Longitude'));
MLS.Data.Altitude   = double(pres2alt(h5read(filename,'/HDFEOS/SWATHS/Temperature/Geolocation Fields/Pressure')));
MLS.Data.Pressure   = double(h5read(filename,'/HDFEOS/SWATHS/Temperature/Geolocation Fields/Pressure'));
MLS.Data.Time       = double(h5read(filename,'/HDFEOS/SWATHS/Temperature/Geolocation Fields/Time'));

MLS.Data.Temperature = double(h5read(filename,'/HDFEOS/SWATHS/Temperature/Data Fields/L2gpValue'));

% CONVERT TIME:
% Think it's seconds since midnight on 1st Jan 1993:
MLS.Data.Time = datenum(1993,01,01) + (MLS.Data.Time ./ 60 ./ 60 ./ 24);

% Created:
MLS.CreatedOn = datestr(today);




end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% output of h5disp on an .he5 Aura MLS file %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HDF5 MLS-Aura_L2GP-Temperature_v04-23-c01_2017d224.he5 
% Group '/' 
%     Group '/HDFEOS' 
%         Group '/HDFEOS/ADDITIONAL' 
%             Group '/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES' 
%                 Attributes:
%                     'OrbitNumber':  69546 69547 69548 69549 69550 69551 69552 69553 69554 69555 69556 69557 69558 69559 69560 -1 
%                     'OrbitPeriod':  5933.000398 5932.965705 5932.967781 5933.020213 5933.100966 5933.112973 5933.085400 5933.006408 5932.956255 5932.959607 5932.978691 5933.050554 5933.127371 5933.072164 5933.072164 0.000000 
%                     'InstrumentName':  'MLS Aura'
%                     'HostName':  'jackal.jpl.nasa.gov'
%                     'ProcessLevel':  'L2'
%                     'PGEVersion':  'V04-23'
%                     'StartUTC':  '2017-08-12T00:00:00.000000Z'
%                     'EndUTC':  '2017-08-12T23:59:59.999999Z'
%                     'GranuleMonth':  8 
%                     'GranuleDay':  12 
%                     'GranuleDayOfYear':  224 
%                     'GranuleYear':  2017 
%                     'TAI93At0zOfGranule':  776649610.000000 
%                     'FirstMAF':  6886041 
%                     'LastMAF':  6889548 
%                     'A Priori l2gp':  ' '
%                     'A Priori l2aux':  ' '
%                     'A Priori ncep':  ' '
%                     'A Priori gmao':  ' '
%                     'A Priori geos5':  '/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0000.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0000.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0000.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0300.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0300.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0600.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0600.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0900.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0900.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_1200.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_1200.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_1500.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_1500.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_1800.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_1800.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_2100.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_2100.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170813_0000.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170813_0000.V01.nc4'
%                     'geos5 type':  'geos5_7,geos5_7,geos5_7,geos5_7,geos5_7,geos5_7,geos5_7,geos5_7,geos5_7,geos5_7,geos5_7,geos5_7,geos5_7,geos5_7,geos5_7,geos5_7,geos5_7,geos5_7,geos5_7'
%                     'MiscNotes':  'apriori(geos5),/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0000.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0000.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0000.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0300.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0300.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0600.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0600.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0900.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_0900.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_1200.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_1200.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_1500.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_1500.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_1800.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_1800.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_2100.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170812_2100.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170813_0000.V01.nc4,/workops/jobs/science/1502642672.0000001/DFPITI3NVASMMLS.fpit.asm.inst3_3d_asm_Nv.GEOS5124.20170813_0000.V01.nc4'
%                     'PCF1':  '
% 
%                     'identifier_product_DOI':  '10.5067/AURA/MLS/DATA2001'
%                     'ProductionLocation':  'jackal.jpl.nasa.gov'
% 
% 
% 
%         Group '/HDFEOS/SWATHS' 
%             Group '/HDFEOS/SWATHS/Temperature' 
%                 Attributes:
%                     'Pressure':  1000.000000 825.404175 681.292053 562.341309 464.158875 383.118683 316.227753 261.015717 215.443466 177.827942 146.779922 121.152763 100.000000 82.540421 68.129204 56.234131 46.415890 38.311867 31.622776 26.101572 21.544348 17.782795 14.677993 12.115276 10.000000 8.254042 6.812921 5.623413 4.641589 3.831187 3.162278 2.610157 2.154435 1.778279 1.467799 1.211528 1.000000 0.681292 0.464159 0.316228 0.215443 0.146780 0.100000 0.046416 0.021544 0.010000 0.004642 0.002154 0.001000 0.000464 0.000215 0.000100 0.000046 0.000022 0.000010 
%                     'VerticalCoordinate':  'Pressure'
%                 Dataset 'nLevels' 
%                     Size:  55
%                     MaxSize:  55
%                     Datatype:   H5T_IEEE_F32LE (single)
%                     ChunkSize:  []
%                     Filters:  none
%                     FillValue:  0.000000
%                     Attributes:
%                         'CLASS':  'DIMENSION_SCALE'
%                         'NAME':  'nLevels'
%                         'REFERENCE_LIST':  H5T_COMPOUND
%                 Dataset 'nTimes' 
%                     Size:  3496
%                     MaxSize:  3496
%                     Datatype:   H5T_IEEE_F32LE (single)
%                     ChunkSize:  []
%                     Filters:  none
%                     FillValue:  0.000000
%                     Attributes:
%                         'CLASS':  'DIMENSION_SCALE'
%                         'NAME':  'nTimes'
%                         'REFERENCE_LIST':  H5T_COMPOUND
%                 Dataset 'nTimesTotal' 
%                     Size:  3496
%                     MaxSize:  3496
%                     Datatype:   H5T_IEEE_F32LE (single)
%                     ChunkSize:  []
%                     Filters:  none
%                     FillValue:  0.000000
%                     Attributes:
%                         'CLASS':  'DIMENSION_SCALE'
%                         'NAME':  'nTimesTotal'
%                 Group '/HDFEOS/SWATHS/Temperature/Data Fields' 
%                     Link:  'Temperature'
%                         Type:  'soft link'
%                         Target:  'L2gpValue'
%                     Link:  'TemperaturePrecision'
%                         Type:  'soft link'
%                         Target:  'L2gpPrecision'
%                     Dataset 'AscDescMode' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_STD_I32LE (int32)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999
%                         Attributes:
%                             'Units':  'NoUnits'
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  -999 
%                             'Title':  'TemperatureAscDescMode'
%                             'MissingValue':  0 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                     Dataset 'Convergence' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'TemperatureConvergence'
%                             'Units':  'NoUnits'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'L2gpPrecision' 
%                         Size:  55x3496
%                         MaxSize:  55x3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  55x120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'TemperaturePrecision'
%                             'Units':  'K'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'L2gpValue' 
%                         Size:  55x3496
%                         MaxSize:  55x3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  55x120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'Temperature'
%                             'Units':  'K'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'Quality' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  -999.989990 
%                             'Title':  'TemperatureQuality'
%                             'Units':  'NoUnits'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                     Dataset 'Status' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_STD_I32LE (int32)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  513
%                         Attributes:
%                             'Units':  'NoUnits'
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  513 
%                             'Title':  'TemperatureStatus'
%                             'MissingValue':  513 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                 Group '/HDFEOS/SWATHS/Temperature/Geolocation Fields' 
%                     Dataset 'ChunkNumber' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_STD_I32LE (int32)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999
%                         Attributes:
%                             'Units':  'NoUnits'
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  -999 
%                             'Title':  'ChunkNumber'
%                             'MissingValue':  -999 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                     Dataset 'Latitude' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'Latitude'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'LineOfSightAngle' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'LineOfSightAngle'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'LocalSolarTime' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             'Units':  'h'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                             '_FillValue':  -999.989990 
%                             'Title':  'LocalSolarTime'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'Longitude' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'Longitude'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'OrbitGeodeticAngle' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'OrbitGeodeticAngle'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'Pressure' 
%                         Size:  55
%                         MaxSize:  55
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  55
%                         Filters:  none
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'Pressure'
%                             'Units':  'hPa'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'Aura-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'SolarZenithAngle' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'SolarZenithAngle'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'Time' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F64LE (double)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             'Title':  'Time'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'Aura-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  -999.989990 
%                             'Units':  's'
%             Group '/HDFEOS/SWATHS/Temperature-APriori' 
%                 Attributes:
%                     'Pressure':  1000.000000 825.404175 681.292053 562.341309 464.158875 383.118683 316.227753 261.015717 215.443466 177.827942 146.779922 121.152763 100.000000 82.540421 68.129204 56.234131 46.415890 38.311867 31.622776 26.101572 21.544348 17.782795 14.677993 12.115276 10.000000 8.254042 6.812921 5.623413 4.641589 3.831187 3.162278 2.610157 2.154435 1.778279 1.467799 1.211528 1.000000 0.681292 0.464159 0.316228 0.215443 0.146780 0.100000 0.046416 0.021544 0.010000 0.004642 0.002154 0.001000 0.000464 0.000215 0.000100 0.000046 0.000022 0.000010 
%                     'VerticalCoordinate':  'Pressure'
%                 Dataset 'nLevels' 
%                     Size:  55
%                     MaxSize:  55
%                     Datatype:   H5T_IEEE_F32LE (single)
%                     ChunkSize:  []
%                     Filters:  none
%                     FillValue:  0.000000
%                     Attributes:
%                         'CLASS':  'DIMENSION_SCALE'
%                         'REFERENCE_LIST':  H5T_COMPOUND
%                         'NAME':  'nLevels'
%                 Dataset 'nTimes' 
%                     Size:  3496
%                     MaxSize:  3496
%                     Datatype:   H5T_IEEE_F32LE (single)
%                     ChunkSize:  []
%                     Filters:  none
%                     FillValue:  0.000000
%                     Attributes:
%                         'CLASS':  'DIMENSION_SCALE'
%                         'NAME':  'nTimes'
%                         'REFERENCE_LIST':  H5T_COMPOUND
%                 Dataset 'nTimesTotal' 
%                     Size:  3496
%                     MaxSize:  3496
%                     Datatype:   H5T_IEEE_F32LE (single)
%                     ChunkSize:  []
%                     Filters:  none
%                     FillValue:  0.000000
%                     Attributes:
%                         'CLASS':  'DIMENSION_SCALE'
%                         'NAME':  'nTimesTotal'
%                 Group '/HDFEOS/SWATHS/Temperature-APriori/Data Fields' 
%                     Link:  'Temperature-APriori'
%                         Type:  'soft link'
%                         Target:  'L2gpValue'
%                     Link:  'Temperature-APrioriPrecision'
%                         Type:  'soft link'
%                         Target:  'L2gpPrecision'
%                     Dataset 'AscDescMode' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_STD_I32LE (int32)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999
%                         Attributes:
%                             'Units':  'NoUnits'
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  -999 
%                             'Title':  'Temperature-APrioriAscDescMode'
%                             'MissingValue':  0 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                     Dataset 'Convergence' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'Temperature-APrioriConvergence'
%                             'Units':  'NoUnits'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'L2gpPrecision' 
%                         Size:  55x3496
%                         MaxSize:  55x3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  55x120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'Temperature-APrioriPrecision'
%                             'Units':  'K'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'L2gpValue' 
%                         Size:  55x3496
%                         MaxSize:  55x3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  55x120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'Temperature-APriori'
%                             'Units':  'K'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'Quality' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  -999.989990 
%                             'Title':  'Temperature-APrioriQuality'
%                             'Units':  'NoUnits'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                     Dataset 'Status' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_STD_I32LE (int32)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  513
%                         Attributes:
%                             'Units':  'NoUnits'
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  513 
%                             'Title':  'Temperature-APrioriStatus'
%                             'MissingValue':  513 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                 Group '/HDFEOS/SWATHS/Temperature-APriori/Geolocation Fields' 
%                     Dataset 'ChunkNumber' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_STD_I32LE (int32)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999
%                         Attributes:
%                             'Units':  'NoUnits'
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  -999 
%                             'Title':  'ChunkNumber'
%                             'MissingValue':  -999 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                     Dataset 'Latitude' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'Latitude'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'LineOfSightAngle' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'LineOfSightAngle'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'LocalSolarTime' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  -999.989990 
%                             'Title':  'LocalSolarTime'
%                             'Units':  'h'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                     Dataset 'Longitude' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'Longitude'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'OrbitGeodeticAngle' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'OrbitGeodeticAngle'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'Pressure' 
%                         Size:  55
%                         MaxSize:  55
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  55
%                         Filters:  none
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'Pressure'
%                             'Units':  'hPa'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'Aura-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'SolarZenithAngle' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'SolarZenithAngle'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'Time' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F64LE (double)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             'Title':  'Time'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'Aura-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  -999.989990 
%                             'Units':  's'
%             Group '/HDFEOS/SWATHS/WMOTPPressure' 
%                 Attributes:
%                     'Pressure':  -999.989990 
%                     'VerticalCoordinate':  'Pressure'
%                 Dataset 'nTimes' 
%                     Size:  3496
%                     MaxSize:  3496
%                     Datatype:   H5T_IEEE_F32LE (single)
%                     ChunkSize:  []
%                     Filters:  none
%                     FillValue:  0.000000
%                     Attributes:
%                         'CLASS':  'DIMENSION_SCALE'
%                         'NAME':  'nTimes'
%                         'REFERENCE_LIST':  H5T_COMPOUND
%                 Dataset 'nTimesTotal' 
%                     Size:  3496
%                     MaxSize:  3496
%                     Datatype:   H5T_IEEE_F32LE (single)
%                     ChunkSize:  []
%                     Filters:  none
%                     FillValue:  0.000000
%                     Attributes:
%                         'CLASS':  'DIMENSION_SCALE'
%                         'NAME':  'nTimesTotal'
%                 Group '/HDFEOS/SWATHS/WMOTPPressure/Data Fields' 
%                     Link:  'WMOTPPressure'
%                         Type:  'soft link'
%                         Target:  'L2gpValue'
%                     Link:  'WMOTPPressurePrecision'
%                         Type:  'soft link'
%                         Target:  'L2gpPrecision'
%                     Dataset 'AscDescMode' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_STD_I32LE (int32)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999
%                         Attributes:
%                             'Units':  'NoUnits'
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  -999 
%                             'Title':  'WMOTPPressureAscDescMode'
%                             'MissingValue':  0 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                     Dataset 'Convergence' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'WMOTPPressureConvergence'
%                             'Units':  'NoUnits'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'L2gpPrecision' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'WMOTPPressurePrecision'
%                             'Units':  'hPa'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'L2gpValue' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'WMOTPPressure'
%                             'Units':  'hPa'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'Quality' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  -999.989990 
%                             'Title':  'WMOTPPressureQuality'
%                             'Units':  'NoUnits'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                     Dataset 'Status' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_STD_I32LE (int32)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  513
%                         Attributes:
%                             'Units':  'NoUnits'
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  513 
%                             'Title':  'WMOTPPressureStatus'
%                             'MissingValue':  513 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                 Group '/HDFEOS/SWATHS/WMOTPPressure/Geolocation Fields' 
%                     Dataset 'ChunkNumber' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_STD_I32LE (int32)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999
%                         Attributes:
%                             'Units':  'NoUnits'
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  -999 
%                             'Title':  'ChunkNumber'
%                             'MissingValue':  -999 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                     Dataset 'Latitude' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'Latitude'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'LineOfSightAngle' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'LineOfSightAngle'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'LocalSolarTime' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  -999.989990 
%                             'Title':  'LocalSolarTime'
%                             'Units':  'h'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                     Dataset 'Longitude' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'Longitude'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'OrbitGeodeticAngle' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'OrbitGeodeticAngle'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'MLS-Specific'
%                             'DIMENSION_LIST':  H5T_VLEN
%                     Dataset 'SolarZenithAngle' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F32LE (single)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             'DIMENSION_LIST':  H5T_VLEN
%                             '_FillValue':  -999.989990 
%                             'Title':  'SolarZenithAngle'
%                             'Units':  'deg'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'HIRDLS-MLS-TES-Shared'
%                     Dataset 'Time' 
%                         Size:  3496
%                         MaxSize:  3496
%                         Datatype:   H5T_IEEE_F64LE (double)
%                         ChunkSize:  120
%                         Filters:  deflate(1)
%                         FillValue:  -999.989990
%                         Attributes:
%                             '_FillValue':  -999.989990 
%                             'Title':  'Time'
%                             'Units':  's'
%                             'MissingValue':  -999.989990 
%                             'UniqueFieldDefinition':  'Aura-Shared'
%                             'DIMENSION_LIST':  H5T_VLEN
%     Group '/HDFEOS INFORMATION' 
%         Attributes:
%             'HDFEOSVersion':  'HDFEOS_5.1.14'
%         Dataset 'StructMetadata.0' 
%             Size:  scalar
%             Datatype:   H5T_STRING
%                 String Length: 32000
%                 Padding: H5T_STR_NULLTERM
%                 Character Set: H5T_CSET_ASCII
%                 Character Type: H5T_C_S1
%             ChunkSize:  []
%             Filters:  none
%             FillValue:  '                                            
% 
% 
% 
% 
% 
% 
% 
% 
% 


