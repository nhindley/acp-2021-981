



% Read a netcdf file's variables into a structure


% filepath = '~/ionPrf_C001.2008.183.00.13.G10_2013.3520_nc';





function nc = getnet_octave(filepath,varnums,attnums)

nci = ncinfo(filepath);

if nargin == 1, varnums = 0:length(nci.Variables)-1; end;
% varnums = 0:length(nci.Variables)-1;
if nargin <= 2, attnums = 0:length(nci.Attributes)-1; end;
% attnums = 0:length(nci.Attributes)-1;

varnames = {nci.Variables.Name};

% Cope with dashes in variable names:
varnames = strrep(varnames,'-','_');
% varnames = regexprep(varnames,'-','_');
        

if ~isempty(nci.Attributes),
    attnames = {nci.Attributes.Name};
end

ncid = netcdf_open(filepath,'nc_nowrite');

for v = varnums,
    
    nc.(varnames{v+1}) = double(netcdf_getVar(ncid,v));

end

for a = attnums,
    
    nc.Attributes.(attnames{a+1}) = double(netcdf_getAtt(ncid,netcdf_getConstant('NC_GLOBAL'),attnames{a+1}));
    
end


netcdf_close(ncid); clear ncid;

end


