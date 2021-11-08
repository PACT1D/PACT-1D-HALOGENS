function [  ] = open_surf_source_file( model_path, Times, nlev, ntim, DateStrLen, surface_source_emissions, add_surface_source_HONO,...
    add_surface_source_Cl_Br, BOX_WALL, BOXCH, output_file_comment, output_file_created_by)
%open a file for writing out the reaction rates

%global fill value number set in the initialization routines
global fill_value_netcdf

mech_Parameters;

%----------------------- write chemical rates to netcdf file --------------------
disp(['------>Creating netcdf output files - surf_source.nc']);
%write netcdf file with initial concentration values
ncid = netcdf.create([model_path '/output/surf_source.nc'], 'CLOBBER');

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

if (add_surface_source_HONO == 1)
    %get the number of emissions
    %[NVAR,nlev,ntp] = size(surface_source_emissions.HONO_total_surf_source);
    
    emi_name={''};
    spec_names = get_spec_names( model_path );
    j=get_ind('HONO');
    
    %write reaction names and numbers to file, make space for rates
    emi_name{j}=strrep(spec_names{j},spec_names{j},['total_source_' spec_names{j}] );
    varid(j)  = netcdf.defVar(ncid,emi_name{j},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(j),'units','number cm-3 s-1')
    netcdf.putAtt(ncid,varid(j),'long_name',['total surface source term for - ' emi_name{j}])
    netcdf.putAtt(ncid,varid(j),'_FillValue',fill_value_netcdf)

    %write reaction names and numbers to file, make space for rates
    emi_name{j}=strrep(spec_names{j},spec_names{j},['dark_source_' spec_names{j}] );
    varid(j)  = netcdf.defVar(ncid,emi_name{j},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(j),'units','number cm-3 s-1')
    netcdf.putAtt(ncid,varid(j),'long_name',['dark mechanism surface source term for - ' emi_name{j}])
    netcdf.putAtt(ncid,varid(j),'_FillValue',fill_value_netcdf)

    %write reaction names and numbers to file, make space for rates
    emi_name{j}=strrep(spec_names{j},spec_names{j},['photEnh_source_' spec_names{j}] );
    varid(j)  = netcdf.defVar(ncid,emi_name{j},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(j),'units','number cm-3 s-1')
    netcdf.putAtt(ncid,varid(j),'long_name',['photo-enhanced mechanism surface source term for - ' emi_name{j}])
    netcdf.putAtt(ncid,varid(j),'_FillValue',fill_value_netcdf)
    
    %write reaction names and numbers to file, make space for rates
    emi_name{j}=strrep(spec_names{j},spec_names{j},['HNO3phot_source_' spec_names{j}] );
    varid(j)  = netcdf.defVar(ncid,emi_name{j},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(j),'units','number cm-3 s-1')
    netcdf.putAtt(ncid,varid(j),'long_name',['HNO3 photolysis mechanism surface source term for - ' emi_name{j}])
    netcdf.putAtt(ncid,varid(j),'_FillValue',fill_value_netcdf)
    
    %write reaction names and numbers to file, make space for rates
    emi_name{j}=strrep(spec_names{j},spec_names{j},['soil_source_' spec_names{j}] );
    varid(j)  = netcdf.defVar(ncid,emi_name{j},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(j),'units','number cm-3 s-1')
    netcdf.putAtt(ncid,varid(j),'long_name',['soil emission mechanism surface source term for - ' emi_name{j}])
    netcdf.putAtt(ncid,varid(j),'_FillValue',fill_value_netcdf)
    
    %write reaction names and numbers to file, make space for rates
    emi_name{j}=strrep(spec_names{j},spec_names{j},['acidDisp_source_' spec_names{j}] );
    varid(j)  = netcdf.defVar(ncid,emi_name{j},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(j),'units','number cm-3 s-1')
    netcdf.putAtt(ncid,varid(j),'long_name',['surface source term for - ' emi_name{j}])
    netcdf.putAtt(ncid,varid(j),'_FillValue',fill_value_netcdf)
end

if (add_surface_source_Cl_Br == 1)
    %get the number of emissions
    %[NVAR,nlev,ntp] = size(surface_source_emissions.HONO_total_surf_source);
    
    emi_name={''};
    spec_names = get_spec_names( model_path );
    j=get_ind('Cl2');
    n=get_ind('Br2');
    
    %write reaction names and numbers to file, make space for rates
    emi_name{j}=strrep(spec_names{j},spec_names{j},['total_source_' spec_names{j}] );
    varid(j)  = netcdf.defVar(ncid,emi_name{j},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(j),'units','number cm-3 s-1')
    netcdf.putAtt(ncid,varid(j),'long_name',['total surface source term for - ' emi_name{j}])
    netcdf.putAtt(ncid,varid(j),'_FillValue',fill_value_netcdf)
       
    %write reaction names and numbers to file, make space for rates
    emi_name{j}=strrep(spec_names{j},spec_names{j},['photo_surf_source_' spec_names{j}] );
    varid(j)  = netcdf.defVar(ncid,emi_name{j},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(j),'units','number cm-3 s-1')
    netcdf.putAtt(ncid,varid(j),'long_name',['surface source term based on J value - ' emi_name{j}])
    netcdf.putAtt(ncid,varid(j),'_FillValue',fill_value_netcdf)
    
    %write reaction names and numbers to file, make space for rates
    emi_name{j}=strrep(spec_names{j},spec_names{j},['ClONO2_recycling_' spec_names{j}] );
    varid(j)  = netcdf.defVar(ncid,emi_name{j},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(j),'units','number cm-3 s-1')
    netcdf.putAtt(ncid,varid(j),'long_name',['surface conversion of ClONO2 to - ' emi_name{j}])
    netcdf.putAtt(ncid,varid(j),'_FillValue',fill_value_netcdf)
    
    %write reaction names and numbers to file, make space for rates
    emi_name{j}=strrep(spec_names{j},spec_names{j},['HOCl_recycling_' spec_names{j}] );
    varid(j)  = netcdf.defVar(ncid,emi_name{j},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(j),'units','number cm-3 s-1')
    netcdf.putAtt(ncid,varid(j),'long_name',['surface conversion of HOCl to - ' emi_name{j}])
    netcdf.putAtt(ncid,varid(j),'_FillValue',fill_value_netcdf)
    
    %write reaction names and numbers to file, make space for rates
    emi_name{n}=strrep(spec_names{n},spec_names{n},['total_source_' spec_names{n}] );
    varid(n)  = netcdf.defVar(ncid,emi_name{n},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(n),'units','number cm-3 s-1')
    netcdf.putAtt(ncid,varid(n),'long_name',['total surface source term for - ' emi_name{n}])
    netcdf.putAtt(ncid,varid(n),'_FillValue',fill_value_netcdf)
    
    %write reaction names and numbers to file, make space for rates
    emi_name{n}=strrep(spec_names{n},spec_names{n},['photo_surf_source_' spec_names{n}] );
    varid(n)  = netcdf.defVar(ncid,emi_name{n},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(n),'units','number cm-3 s-1')
    netcdf.putAtt(ncid,varid(n),'long_name',['surface source term based on J value - ' emi_name{n}])
    netcdf.putAtt(ncid,varid(n),'_FillValue',fill_value_netcdf)
    
    %write reaction names and numbers to file, make space for rates
    emi_name{n}=strrep(spec_names{n},spec_names{n},['BrONO2_recycling_' spec_names{n}] );
    varid(n)  = netcdf.defVar(ncid,emi_name{n},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(n),'units','number cm-3 s-1')
    netcdf.putAtt(ncid,varid(n),'long_name',['surface conversion of BrONO2 to - ' emi_name{n}])
    netcdf.putAtt(ncid,varid(n),'_FillValue',fill_value_netcdf)
    
    %write reaction names and numbers to file, make space for rates
    emi_name{n}=strrep(spec_names{n},spec_names{n},['HOBr_recycling_' spec_names{n}] );
    varid(n)  = netcdf.defVar(ncid,emi_name{n},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(n),'units','number cm-3 s-1')
    netcdf.putAtt(ncid,varid(n),'long_name',['surface conversion of HOBr to - ' emi_name{n}])
    netcdf.putAtt(ncid,varid(n),'_FillValue',fill_value_netcdf)
end

emi_name={''};
spec_names = get_spec_names( model_path );
%write reaction names and numbers to file, make space for rates
for j = 1:NVAR
    emi_name{j}=strrep(spec_names{j},spec_names{j},['loss_' spec_names{j}] );
    varid(j)  = netcdf.defVar(ncid,emi_name{j},'double',[time_dimid]); 
    netcdf.putAtt(ncid,varid(j),'units','molec cm-2')
    netcdf.putAtt(ncid,varid(j),'long_name',['total loss to ground for - ' emi_name{j}])
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
if (add_surface_source_HONO == 1)
    j=get_ind('HONO');
    tmp = squeeze(surface_source_emissions.HONO_total_surf_source(1));
    netcdf.putVar(ncid,varid(j),[0],[1],tmp);
    clearvars tmp;
end

if (add_surface_source_Cl_Br == 1)
    j=get_ind('Cl2');
    tmp = squeeze(surface_source_emissions.Cl2_total_surf_source(1));
    netcdf.putVar(ncid,varid(j),[0],[1],tmp);
    clearvars tmp;
end

clear varid;

netcdf.close(ncid)

end

