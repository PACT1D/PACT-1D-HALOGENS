function [  ] = open_emission_file( model_path, Times, nlev, ntim, DateStrLen, emissions, BOX_WALL, BOXCH, output_file_comment, output_file_created_by)
%open a file for writing out the reaction rates

%global fill value number set in the initialization routines
global fill_value_netcdf

mech_Parameters;

%----------------------- write chemical rates to netcdf file --------------------
disp(['------>Creating netcdf output files - emissions.nc']);
%write netcdf file with initial concentration values
ncid = netcdf.create([model_path '/output/emissions.nc'], 'CLOBBER');

% define dimensions
z_dimid = netcdf.defDim(ncid,'z_levels',nlev);
dateStrLen_dimid = netcdf.defDim(ncid,'DateStrLen',DateStrLen);
time_dimid = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));

%define global attributes
varid = netcdf.getConstant('GLOBAL');
netcdf.putAtt(ncid,varid,'creation_date',datestr(now));
netcdf.putAtt(ncid,varid,'created_by', output_file_created_by);
netcdf.putAtt(ncid,varid,'model_start',output_file_comment);
netcdf.putAtt(ncid,varid,'note','for 1D model on model levels');
clear varid;

% define variables
varid_timeStr = netcdf.defVar(ncid,'Times','char',[dateStrLen_dimid,time_dimid]);
netcdf.putAtt(ncid,varid_timeStr,'units','Time and date (local)')
netcdf.putAtt(ncid,varid_timeStr,'long_name','Character string - current time')

varid_box = netcdf.defVar(ncid,'BOX_WALL','double',[z_dimid]);
netcdf.putAtt(ncid,varid_box,'units','m')
netcdf.putAtt(ncid,varid_box,'long_name','Vertical grid spacing, upper boundary of grid box in meters')

varid_boxch = netcdf.defVar(ncid,'BOXCH','double',[z_dimid]);
netcdf.putAtt(ncid,varid_boxch,'units','m')
netcdf.putAtt(ncid,varid_boxch,'long_name','Vertical grid spacing, box center point hight in meters')

%get the number of emissions
[NVAR,nlev,ntp] = size(emissions);

emi_name={''};
spec_names = get_spec_names( model_path );
%write reaction names and numbers to file, make space for rates
for j = 1:NVAR
  emi_name{j}=strrep(spec_names{j},spec_names{j},['E_' spec_names{j}] );
  varid(j)  = netcdf.defVar(ncid,emi_name{j},'double',[z_dimid,time_dimid]); 
  netcdf.putAtt(ncid,varid(j),'units','mol m-2 s-1')
  netcdf.putAtt(ncid,varid(j),'long_name',['emission rates - ' emi_name{j}])
  netcdf.putAtt(ncid,varid(j),'_FillValue',fill_value_netcdf)
end

netcdf.endDef(ncid)
%done defining file format

% ----------WRITE DATA TO THE FILE - BOX HEIGHTS, Times, species concentraitons

% box center points
netcdf.putVar(ncid,varid_box,[0],[nlev],BOX_WALL);
netcdf.putVar(ncid,varid_boxch,[0],[nlev],BOXCH);

%Time string for each value
netcdf.putVar(ncid,varid_timeStr,[0, 0],[DateStrLen, 1],Times(1,:)');

%add reaction rates
for j=1:NVAR
    tmp = squeeze(emissions(j,:,1));
    netcdf.putVar(ncid,varid(j),[0,0],[nlev,1],tmp);
    clearvars tmp;
end

clear varid;

netcdf.close(ncid)

end

