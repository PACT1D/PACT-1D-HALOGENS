%initialize model - read all input files & give errors if this doesn't work 

%  global variable for the list of the photolysis rates 
   global jrates_list
   
%  global fill value, for missing values in netcdf files and initial values
%  for all arrays
   global fill_value_netcdf
   fill_value_netcdf = 9.9e+36;

% read in master time and level file
   time_level_file   = [model_path '/input/master_time_lev.nc'];
   [BOX_WALL, BOXCH, Times]    = read_times_levels(time_level_file);       %BOX_WALL - box upper boundary, BOXCH box centerpoint

 % for defining the size of arrays 
   [NTIM_RUN, timeStrLen] = size(Times);                                   %number of times for the run (including initial and final time)
   NTIM_INTEGRATION = NTIM_RUN-1;                                          %number of time steps to perform integration
   [years, months, days, hours, minutes, seconds, datenum_input] ...
       = get_times(Times);                                                 %get the date/time for each timestep, datenum_input is in whole fractional day
   dt_input            ...
       = double(int16((datenum_input(2)-datenum_input(1))*24*3600));       %in seconds, assume equally spaced timesteps in the input file, can be checked later
   [NLEV, ~]           = size(BOX_WALL);                                   %number of model levels
   %these are defined in mechanism/mech_Parameters.m
   NSPEC;                                                                  %number of species - defined by kpp , defined in kpp mech_Parameters  - variable and fixed concentration
   NVAR;                                                                   %number of variable concentraiton species - defined by kpp, defined in kpp mech_Parameters 
   NFIX;                                                                   %number of fixed concentraiton species - defined by kpp, defined in kpp mech_Parameters 
   %NJRATES
   NREACT;                                                                 %number of reactions - defined by kpp, defined in kpp mech_Parameters     
 
 % initialize model time steps
   NTIM_CHEM    = NTIM_INTEGRATION*(dt_input/dt_chem)+1;                   %get the number of timesteps in terms of the chemistry timestep, one extra for the inital conc
   err_msg      = check_timesteps(dt_input, dt_chem, dt_kpp);
   if err_msg
       clean_exit('initialize_model.m',err_msg)
   end
  
   %calculate the time in elapsed seconds the model run ends
   tfinal       = (datenum_input(NTIM_RUN)-datenum_input(1))*3600*24;      %model run time in elapsed seconds
   model_times  = linspace( 0, tfinal, NTIM_CHEM );                        %model time in seconds
   datenum_chem = datenum_input(1)+model_times/24./3600.;                  %the datenum (fractional days) for each timestep where we will call chemistry
   
  % initialize the temperature and pressure, currently just numbers, fix
   pressure          = zeros(NLEV,NTIM_CHEM);
   pressure(:,:)     = fill_value_netcdf;
   temperature       = zeros(NLEV,NTIM_CHEM);
   temperature(:,:)  = fill_value_netcdf;
   relative_humidity = zeros(NLEV,NTIM_CHEM);
   relative_humidity(:,:)  = fill_value_netcdf;
   atmfile = [model_path '/input/atm_properties.nc'];
   [pressure_input,temperature_input,relative_humidity_input,err_msg] = read_atm(atmfile,BOX_WALL,Times);
   if err_msg
       clean_exit('initialize_model.m',err_msg)
   end
  
   %interpolate T, P, RH in time ot the chemistry integration times
   for l=1:NLEV
       pressure(l,:)          = interp1(datenum_input,pressure_input(l,:),datenum_chem,'linear');
       relative_humidity(l,:) = interp1(datenum_input,relative_humidity_input(l,:),datenum_chem,'linear');
       temperature(l,:)       = interp1(datenum_input,temperature_input(l,:),datenum_chem,'linear');
   end
       
   % create an empty array for the species, only store the variable concentration species (NVAR)
   % these are the arrays for output

   spec                     = zeros(NVAR,NLEV,NTIM_CHEM);                  %variable concentration species
   spec(:,:,:)              = fill_value_netcdf;
   spec_fixed               = zeros(NFIX,NLEV,NTIM_CHEM);                  %fixed concentraiton species (assume same concentration in all levels at all times
   spec_fixed(:,:,:)        = fill_value_netcdf;
   VT                       = zeros(NVAR,NLEV,NTIM_CHEM);                  %vertical transport rates
   VT(:,:,1)                = 0.;
   depo                     = zeros(NVAR,NTIM_CHEM);                       %deposition rate to surface
   depo(:,1)                = 0.;
   total_loss_to_ground     = zeros(NVAR,1);                               %number of molecules lost to the surface - cumulative over the model run
   total_loss_to_ground(1)  = 0.;
   rates                    = zeros(NREACT,NLEV,NTIM_CHEM);                %rates from solving the chemistry rates
   rates(:,:,:)             = fill_value_netcdf;
   rates(:,:,1)             = 0.;
   rate_constants           = zeros(NREACT,NLEV,NTIM_CHEM);                %rate constants used online
   rate_constants(:,:,:)    = fill_value_netcdf;
   rate_constants(:,:,1)    = 0.;
   emissions                = zeros(NVAR,NLEV,NTIM_CHEM);                  %emission rates 
   emissions(:,:,:)         = fill_value_netcdf;
   emissions(:,:,1)         = 0.;
   soil_emi_NO              = zeros(1,NTIM_CHEM);                          %emission rates
   surface_source_emissions = struct;                                      %create empty structure for surface source emissions
   
   %Create empty arrays for HONO surface source terms in
   %surface_source_emissions structure if option is selected
   if (add_surface_source_HONO == 1)
       %Define array for HONO at 1st model level for all timesteps
       surface_source_emissions.HONO_total_surf_source             = zeros(1,1,NTIM_CHEM);             %surface source - calculated online of HONO and other species 
       surface_source_emissions.HONO_total_surf_source(:,1,1)      = 0.;
       surface_source_emissions.HONO_dark_surf_source              = zeros(1,1,NTIM_CHEM);             %surface source - calculated online of HONO and other species 
       surface_source_emissions.HONO_dark_surf_source(:,1,1)       = 0.;
       surface_source_emissions.HONO_photEnh_surf_source           = zeros(1,1,NTIM_CHEM);             %surface source - calculated online of HONO and other species 
       surface_source_emissions.HONO_photEnh_surf_source(:,1,1)    = 0.;
       surface_source_emissions.HONO_HNO3phot_surf_source          = zeros(1,1,NTIM_CHEM);             %surface source - calculated online of HONO and other species 
       surface_source_emissions.HONO_HNO3phot_surf_source(:,1,1)   = 0.;
       surface_source_emissions.HONO_soil_surf_source              = zeros(1,1,NTIM_CHEM);             %surface source - calculated online of HONO and other species 
       surface_source_emissions.HONO_soil_surf_source(:,1,1)       = 0.;
       surface_source_emissions.HONO_acidDisp_surf_source          = zeros(1,1,NTIM_CHEM);             %surface source - calculated online of HONO and other species 
       surface_source_emissions.HONO_acidDisp_surf_source(:,1,1)   = 0.;
   end
   
   %Create empty arrays for Cl and Br surface source terms in
   %surface_source_emissions structure if option is selected
   if (add_surface_source_Cl_Br == 1)
       %Define array for Cl and Br at 1st model level for all timesteps
       surface_source_emissions.Cl2_total_surf_source             = zeros(1,1,NTIM_CHEM);             %surface source - Summed total for chlorine surface emissions 
       surface_source_emissions.Cl2_total_surf_source(:,1,1)      = 0.;
       surface_source_emissions.Cl2_photo_surf_source             = zeros(1,1,NTIM_CHEM);             %surface source - Total for photochemical emission of chlorine 
       surface_source_emissions.Cl2_photo_surf_source(:,1,1)      = 0.;
       surface_source_emissions.Cl2_clono2_surf_source            = zeros(1,1,NTIM_CHEM);             %surface source - Surface recycling conversion of ClONO2 to Cl2 
       surface_source_emissions.Cl2_clono2_surf_source(:,1,1)     = 0.;
       surface_source_emissions.Cl2_hocl_surf_source              = zeros(1,1,NTIM_CHEM);             %surface source - Surface recycling conversion of HOCl to Cl2  
       surface_source_emissions.Cl2_hocl_surf_source(:,1,1)       = 0.;
       surface_source_emissions.Br2_total_surf_source             = zeros(1,1,NTIM_CHEM);             %surface source - Summed total for bromine surface emissions 
       surface_source_emissions.Br2_total_surf_source(:,1,1)      = 0.;
       surface_source_emissions.Br2_photo_surf_source             = zeros(1,1,NTIM_CHEM);             %surface source - Total for photochemical emission of bromine 
       surface_source_emissions.Br2_photo_surf_source(:,1,1)      = 0.;
       surface_source_emissions.Br2_brono2_surf_source            = zeros(1,1,NTIM_CHEM);             %surface source - Surface recycling conversion of BrONO2 to Br2 
       surface_source_emissions.Br2_brono2_surf_source(:,1,1)     = 0.;
       surface_source_emissions.Br2_hobr_surf_source              = zeros(1,1,NTIM_CHEM);             %surface source - Surface recycling conversion of HOBr to Br2  
       surface_source_emissions.Br2_hobr_surf_source(:,1,1)       = 0.;
   end
   
   %Space to add other surface source options here (as above)
   %if (add_surface_source_Iodine == 1)
   %
   %end
   
   % read Kz values from netcdf file
   kz_file                              = [model_path '/input/Kz.nc'];
   [Kz_input,err_msg]                   = read_kz(kz_file,BOX_WALL,Times);
   if err_msg
       clean_exit('initialize_model.m',err_msg)
   end
   
   %make space for the Kz values at the chemistry model times
   Kz = zeros(NLEV,NTIM_CHEM);
   Kz(:,:,:) = fill_value_netcdf;
   %interpolate in time to the chemistry times
   for l=1:NLEV
       Kz(l,:) = interp1(datenum_input,Kz_input(l,:),datenum_chem,'linear');
   end
   
 % read jrate file
   jrate_file                           = [model_path '/input/J_values.nc'];
   [jrates_in, jrates_list, err_msg]    = read_jrates(jrate_file,BOX_WALL,Times);
   if err_msg
       clean_exit('initialize_model.m',err_msg)
   end
   
   %set up space for the photolysis rates
   [NJRATES, ~, ~] = size(jrates_in);
   jrates = zeros(NJRATES,NLEV,NTIM_CHEM);
   jrates(:,:,:) = fill_value_netcdf;
   %interplate jrates in time, to the chemistry timesteps
   for i=1:NJRATES
       for l=1:NLEV   
           jrates(i,l,:)    = interp1(squeeze(datenum_input),squeeze(jrates_in(i,l,:)),squeeze(datenum_chem),'linear');
       end
   end

  %read aerosol data 
  aerosol_file                              = [model_path '/input/aerosol.nc'];
  [aerosol_radius_in, aerosol_number_in, err_msg] = read_aerosol(aerosol_file,BOX_WALL,Times);
  %if there is an error, display it and exit
  if err_msg
      clean_exit('initialize_model.m',err_msg)
  end
  
  %make space for the aerosol radius and aerosol number concentration values
  aerosol_radius = zeros(NLEV,NTIM_CHEM);
  aerosol_radius(:,:,:) = fill_value_netcdf;
  aerosol_number = zeros(NLEV,NTIM_CHEM);
  aerosol_number(:,:,:) = fill_value_netcdf;
  %interplate the aerosol information in time, to the chemistry timesteps
  for l=1:NLEV   
     aerosol_radius(l,:)    = interp1(squeeze(datenum_input),squeeze(aerosol_radius_in(l,:)),squeeze(datenum_chem),'linear');
     aerosol_number(l,:)    = interp1(squeeze(datenum_input),squeeze(aerosol_number_in(l,:)),squeeze(datenum_chem),'linear'); 
  end
 
  %get the list of all reactants, this is used to calculate the reaciton rates
  global reaction_comment
  reaction_comment = get_rxn_desc( model_path );
  global reactant_list;
  reactant_list = get_list_of_reactants( );
  
  % initialize the name of species and inital species concentrations for everything that has a concentration that varies
  spec_name_file                       = [model_path '/mechanism/mech.spc'];
  initial_spec_file                    = [model_path '/input/initial_spec_profiles.nc'];
  temperature_t0                       = temperature(:,1);
  pressure_t0                          = pressure(:,1);
  relative_humidity_t0                 = relative_humidity(:,1);
  [spec,spec_fixed,err_msg]            = init_spec(initial_spec_file, NLEV, spec, spec_fixed, Times, BOX_WALL, ...
                                                    temperature_t0, pressure_t0, relative_humidity_t0);    %set species names and spec cocnetration at NTIM=1
  % if there was a problem, display the error message 
  if err_msg
      clean_exit('initialize_model.m',err_msg)
  end
  
  %only read the emissions file if we will add emisisons, otherwise the
  %emissions array will contain missing values
  if (addemissions==1)
      % initialize the emissions
      emission_rates_file                  = [model_path '/input/Emission_rates.nc'];
      [emissions_input, err_msg] = add_emissions(emission_rates_file,NLEV, Times, BOX_WALL);
      if err_msg
          clean_exit('initialize_model.m',err_msg)
      end
      %[NEMISSIONS, ~, ~] = size(emissions_input);
   
      %interpolate emissions in time ot the chemistry integration times
      for i=1:NVAR
          for l=1:NLEV
              emissions(i,l,:)  = interp1(squeeze(datenum_input),squeeze(emissions_input(i,l,:)),squeeze(datenum_chem),'linear');
          end
      end
      
      %Add NO soil emissions if option is selected
      if (add_NO_soil_emissions == 1)
          soil_emi_file                              = [model_path '/input/NO_soil_emi.nc'];
          [soil_emi_NO_input,err_msg]                   = read_soil_emi(soil_emi_file,Times,BOX_WALL);
          if err_msg
              clean_exit('initialize_model.m',err_msg)
          end
          
          soil_emi_NO(1,:) = interp1(datenum_input,soil_emi_NO_input,datenum_chem,'linear');
      end
  end
   
  % set diffusion constants directly with interpolated values of pressure
  % and temperature
  global MW;   %global array with molecular weights
  diffusion_constant = set_diffusion_constant(temperature,pressure,model_path);  %diffusion constant
  %surface reaction probilibites & index for heterogenous reactions, read in from the mechanism
  [aerosol_rxn_probabilities, reactant_name] = get_aerosol_reaction_probabilities(model_path);
  %calcualte surface reaction rate constants in 1/s (K_het)
  K_het = set_heterogenous_reactions(temperature, diffusion_constant, aerosol_radius, aerosol_number, aerosol_rxn_probabilities, reactant_name);

  
  % set surface uptake rates (deposition), to be used in the diffusion routines
  % surface reaction probilibites at the ground 
  [ground_uptake_probability] = get_surface_uptake_probabilities(model_path);

  % calculate surface uptake, to be mulitpled by the concentrations in the
  % diffusion routine
  spec_names = get_spec_names(model_path);
  Eff_dep_surf = set_surface_velocity(temperature,  ground_uptake_probability, spec_names);
  
  % initialize the surface storage term if option is selected
  if (initialize_surface_concentrations == 1)
      spec_name_file                       = [model_path '/mechanism/mech.spc'];
      initial_surf_stor_file               = [model_path '/input/initial_surf_storage.nc'];
      [init_total_loss_to_ground,err_msg]  = init_surf_stor(initial_surf_stor_file , total_loss_to_ground(1), Times);   %set species names and spec cocnetration at NTIM=1
      % if there was a problem, display the error message 
      if err_msg
          clean_exit('initialize_model.m',err_msg)
      end
      total_loss_to_ground  = init_total_loss_to_ground ;
  end

  % Read in observed ozone concentrations
  ozoneobs = fopen('input/o3obs1min.txt');
  output = textscan(ozoneobs, '%s');
  fclose(ozoneobs);
  obs_o3 = str2double(output{1}(1:1:end));
  
  % Interpolate observed ozone onto model time grid
  obs_o3_interp = zeros(NTIM_CHEM);
  t1 = datenum(2009,3,18,0,0,0);
  t2 = datenum(2009,3,20,0,0,0);
  t = t1:(1/1440):t2;
  obs_o3_interp = interp1(t,obs_o3,datenum_chem,'linear');  
