% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                                                                  
% The Reaction Rates File                                          
%                                                                  
% Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor  
%       (http://www.cs.vt.edu/~asandu/Software/KPP)                
% KPP is distributed under GPL, the general public licence         
%       (http://www.gnu.org/copyleft/gpl.html)                     
% (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa           
% (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech            
%     With important contributions from:                           
%        M. Damian, Villanova University, USA                      
%        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
%                                                                  
% File                 : mech_Rates.m                              
% Time                 : Mon Jan 24 15:35:18 2022                  
% Working directory    : /home/sahmed/mechanism                    
% Equation file        : mech.kpp                                  
% Output root filename : mech                                      
%                                                                  
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




% Begin Rate Law Functions from KPP_HOME/util/UserRateLaws         

%  User-defined Rate Law functions
%  Note: insert this file at the end of Update_RCONST

%---  Arrhenius
   function [rate] =  ARR( A0,B0,C0 )
      global TEMP CFACTOR
      rate =  (A0) * exp(-(B0)/TEMP) * (TEMP/300.0)^(C0) ;            
   return %  ARR        

%--- Simplified Arrhenius, with two arguments
%--- Note: The argument B0 has a changed sign when compared to ARR
   function [rate] =  ARR2( A0,B0 )
      global TEMP CFACTOR
      rate =  (A0) * exp( (B0)/TEMP ) ;             
   return %  ARR2          

   function [rate] =  EP2(A0,C0,A2,C2,A3,C3)
      global TEMP CFACTOR                       
      K0 = (A0) * exp(-C0/TEMP);
      K2 = (A2) * exp(-C2/TEMP);
      K3 = (A3) * exp(-C3/TEMP);
      K3 = K3*CFACTOR*1.0e+6;
      rate = K0 + K3/(1.0+K3/K2) ;       
   return %  EP2

   function [rate] =  EP3(A1,C1,A2,C2) 
      global TEMP CFACTOR               
      K1 = (A1) * exp(-(C1)/TEMP);
      K2 = (A2) * exp(-(C2)/TEMP);
      rate = K1 + K2*(1.0e+6*CFACTOR);      
   return %  EP3 

   function [rate] =  FALL ( A0,B0,C0,A1,B1,C1,CF)
      global TEMP CFACTOR                      
      K0 = A0 * exp(-B0/TEMP)* (TEMP/300.0)^(C0);
      K1 = A1 * exp(-B1/TEMP)* (TEMP/300.0)^(C1);
      K0 = K0*CFACTOR*1.0e+6;
      K1 = K0/K1;
      rate = (K0/(1.0+K1))*(CF)^(1.0/(1.0+(log(K1))^2));        
   return %  FALL


% End Rate Law Functions from KPP_HOME/util/UserRateLaws           


% Begin INLINED Rate Law Functions                                 


%//////////////////////////  RACM_THERMAL ///////////////////////////
function [rate_constant] =  RACM_THERMAL(temperature,A0,B0)
%////   RACM2 reaction rates have the form K = A * exp(-B / T)
%// for example - for the reaction: HO+O3 -> O2 + O2 
%// rate_constant = 1.70 x 10-12 exp(-940/T)
%// A0 = 1.70 x 10-12
%// B0 = 940
   rate_constant =  (A0) * exp(-(B0)/temperature);
return

%//////////////////////////  RACM_THERMAL_T2 ///////////////////////////
function [rate_constant] =  RACM_THERMAL_T2(temperature,A0,B0)
%///    NEEDS TO BE UPDATED
%////   RACM2 reaction rates have the form K = A * exp(-B / T)
%// for example - for the reaction: HO+O3 -> O2 + O2 
%// rate_constant = 1.70 x 10-12 exp(-940/T)
%// A0 = 1.70 x 10-12
%// B0 = 940
   rate_constant =  temperature^2. * (A0) * exp(-(B0)/temperature);
return

%//////////////////////////  RACM_TROE ///////////////////////////
function [rate_constant] =  RACM_TROE(temperature,M_conc,K0,N,KINF,M)
%//Following the definition in the NASA JPL data eval 2011, section 2.1
%//The rate constant is given by  k_f([M],T)
%//k_0_T = K0 * (T/300.)^(-N)
%//k_inf_T = KINF * (T/300.)^(-M)
%//FACTOR = [1+(log10(k_0_T AIR_CONC/k_inf_T))^2]^-1
%//rate_constant = (k_0_T AIR_CONC )/ (1 + (k_0_T AIR_CONC/k_inf_T) ) * 0.6^(FAC1.0
  k_0_T   = K0 * (temperature/300.)^(-N);
  k_inf_T = KINF * (temperature/300.)^(-M);
%//M_conc = get_M_CONC();
  FACTOR = ( 1  +  (  log10(  k_0_T*M_conc/k_inf_T   )  )^2. )^(-1.);
  rate_constant = (k_0_T*M_conc)/(1+(k_0_T*M_conc/k_inf_T)) * 0.6^(FACTOR); 
return

%//////////////////////////  RACM_TROE_EQUIL ///////////////////////////
function [rate_constant] =  RACM_TROE_EQUIL(temperature,M_conc,K0,N,KINF,M,A0,B01.0
%//Following the definition in the NASA JPL data eval 2011
  K_TROE = RACM_TROE(temperature,M_conc,K0,N,KINF,M); 
  Keq = A0*exp(B0/temperature) ;
  rate_constant = K_TROE/Keq; 
return

%//////////////////////////  RACM_reaction51 ///////////////////////////
function [rate_constant] =  RACM_reaction51(temperature)
%//table 2f from Goloffi electronic supplement - for reaction 51, pressure in Pa1.0
%//units of pressure - PRES_TIME_LEV is Pa
  global PRES_TIME_LEV ;
  TMP3 = 3.45E-12*exp(270./temperature);
  TMP4 = (530./temperature)+(4.8E-6)*PRES_TIME_LEV-1.73;
  rate_constant = TMP3*TMP4/100;
return

%//////////////////////////  RACM_reaction57 ///////////////////////////
function [rate_constant] =  RACM_reaction57(temperature, CM)
%//table 2f from Goloffi electronic supplement - for reaction 57
%//typo in the electornic supplement table 2f, listed as reaction 56
%//  CM = get_M_CONC();
  TMP0 = 2.4E-14*exp(460/temperature);
  TMP2 = 2.7E-17*exp(2199/temperature);
  TMP3 = 6.5E-34*exp(1335/temperature)*CM;
  rate_constant = TMP0 + TMP3/(1+TMP3/TMP2);
return


% End INLINED Rate Law Functions                                   

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%                                                                  
% Update_SUN - update SUN light using TIME                         
%   Arguments :                                                    
%                                                                  
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function Update_SUN( )

global TIME SUN
 
   SunRise = 4.5;
   SunSet  = 19.5;
   Thour = TIME/3600.;
   Tlocal = Thour - floor(Thour/24)*24;

   if ( (Tlocal>=SunRise) & (Tlocal<=SunSet) ) 
     Ttmp = (2.0*Tlocal-SunRise-SunSet)/(SunSet-SunRise);
     if (Ttmp>0) 
       Ttmp =  Ttmp*Ttmp;
     else
       Ttmp = -Ttmp*Ttmp;
     end 
     SUN = ( 1.0 + cos(pi*Ttmp) )/2.0 ;
   else
     SUN = 0.0;
   end 

return % Update_SUN

% End of Update_SUN function                                       
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


