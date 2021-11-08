function  [spec,surface_source] = surface_source_cl_br(spec, dt, BOX_WALL_cm, depo_rate,o3_measured)

% 18/05/2021, Shaddy Ahmed

% Routine to calculate a snowpack surface emission and recycling source
% for chlorine and bromine based on 2 separate mechanisms:
% (1) Photochemical emission source scaled to surface ozone observations 
% (2) Snowpack surface recycling scaled to XONO2 and HOX (X = Cl, Br) deposition  

%global variable
global GStruct TIME

%set the new concentrations equal to the old concentrations
spec_before_surfSource = spec;

% set up array to save surface source for this time step
total_surf_source_t=zeros(2,1);         % temp array for total surface source, this is summed (Chlorine saved in cell 1, bromine saved in cell 2)
photo_surf_source_t=zeros(2,1);         % temp array for photochemcial emission (chlorine is saved in cell 1, bromine saved in cell 2)
clono2_surf_reactions_t=zeros(1,1);     % temp array for surface recycling of clono2
hocl_surf_reactions_t=zeros(1,1);       % temp array for surface recycling of hocl
brono2_surf_reactions_t=zeros(1,1);     % temp array for surface recycling of brono2
hobr_photo_source_t=zeros(1,1);         % temp array for surface recycling of hobr

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chlorine and bromine surface sources
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%source setup, emit into level 1
Cl2_ind   = get_ind('Cl2');
ClONO2_ind   = get_ind('ClONO2');
HOCl_ind   = get_ind('HOCl');
Br2_ind   = get_ind('Br2');
BrONO2_ind   = get_ind('BrONO2');
HOBr_ind   = get_ind('HOBr');

J_Cl2    = jrates_name('J_Cl2', TIME,GStruct);
J_Br2    = jrates_name('J_Br2', TIME,GStruct);

level1 = 1;      %Model level to put a portion of the emissions

%define dT/dZ scaling factor (s/cm) - used to distribute emissions within
%specified level for this timestep
facLev1=double(dt)/BOX_WALL_cm(level1);


%-------------------- Primary photochemical emission --------------------
                %----------- Chlorine -----------
% Set scaling factor for chlorine primary emission
cl2_emis = 20;

% Read in measured ozone concentrations at 1.5 m (molec cm-3)
o3val = o3_measured;

photo_surf_source_t(1,level1) = cl2_emis*(J_Cl2^0.5)*o3val*double(dt);
total_surf_source_t(1,level1) = total_surf_source_t(1,level1)+photo_surf_source_t(1,level1);

                %----------- Bromine -----------
% Release of Br2 from photochemical source
photo_surf_source_t(2,level1) = J_Br2*o3val*double(dt);
total_surf_source_t(2,level1) = total_surf_source_t(2,level1)+photo_surf_source_t(2,level1);


%-------------------- Surface recycling reactions --------------------
                %----------- Chlorine -----------
% Set conversion efficiency of ClONO2 and HOCl to Cl2
ConvEff_Cl = 0.1;   

% Release of Cl2 from surface conversion of ClONO2 and HOCl to Cl2
clono2_surf_reactions_t(1,level1) = depo_rate(ClONO2_ind)*facLev1*ConvEff_Cl;
hocl_surf_reactions_t(1,level1) = depo_rate(HOCl_ind)*facLev1*ConvEff_Cl;
total_surf_source_t(1,level1) = total_surf_source_t(1,level1)+clono2_surf_reactions_t(1,level1)+hocl_surf_reactions_t(1,level1);

                %----------- Bromine -----------
% Set conversion efficiency of BrONO2 and HOBr to Br2
ConvEff_Br = 0.6;

% Release of Br2 from surface conversion of HOBr and BrONO2 to Br2
brono2_surf_reactions_t(1,level1) = depo_rate(BrONO2_ind)*facLev1*ConvEff_Br;
hobr_surf_reactions_t(1,level1) = depo_rate(HOBr_ind)*facLev1*ConvEff_Br;
total_surf_source_t(2,level1) = total_surf_source_t(2,level1)+brono2_surf_reactions_t(1,level1)+hobr_surf_reactions_t(1,level1);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      END of Cl and Br surface sources
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Update spec_new %
spec(Cl2_ind,level1) = spec_before_surfSource(Cl2_ind,level1)+total_surf_source_t(1,level1);
spec(Br2_ind,level1) = spec_before_surfSource(Br2_ind,level1)+total_surf_source_t(2,level1);

% calculate the surface source rate - in molec cm-3 s-1
surface_source.total_surf_source_t = total_surf_source_t/double(dt);
surface_source.photo_surf_source_t = photo_surf_source_t/double(dt);
surface_source.clono2_surf_reactions_t = clono2_surf_reactions_t/double(dt);
surface_source.hocl_surf_reactions_t = hocl_surf_reactions_t/double(dt);
surface_source.brono2_surf_reactions_t = brono2_surf_reactions_t/double(dt);
surface_source.hobr_surf_reactions_t = hobr_surf_reactions_t/double(dt);

end