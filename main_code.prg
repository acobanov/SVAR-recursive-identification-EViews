'************************************************************************************************************************************
'  Toolbox for estimating SVAR model with recursive identification method 
'  Main code
'************************************************************************************************************************************

' Author:  ANTE COBANOV (ante.cobanov@yahoo.com)
' GitHub: acobanov

'************************************************************************************************************************************
' List of functions:
'    {create_system}          -> creates and estimates EViews system of equations
'    {var_representation}  -> calculates all VAR representations (MA, structural, reduced-form,VAR(1))
'    {hd}                                -> calculates hystorical decompositions of all endogenous variables
'    {irf}                                 -> calculates impulse responses and accumulated impulse responses
'    {spillover_index}        ->  calculates spillover index and creates spillover table (Diebold -Yilmaz  measure)

'***********************************************************************************************************************************
' Define your working directory
cd "..."
' Add useful functions (SVAR representation)
include "...\create_system"
include "...\var_representation"
include "...\irf"
include "...\hd"
include "...\fevd"
include "...\spillover_index"

' Define workfile name
%wfname="Spillover_index"       
' Read excel data
%path="...\data\dy_ej2009.xls" 
' Define excel sheet name
%xls_sheet="weekly_realreturns"                           
' Define page names
%pgname1="input"                   
%pgname2="test_1"
%pgname3="test_2"
' Define time span (chosen from excel file)
%obs_beg_all="2000M01"
%obs_end_all=@datestr(@dateadd(@dateval(%obs_beg_all),828, "mm"),"yyyy[M]mm")               
' Variable names
%variables_all="rrdjia rrftse rrfra rrger rrhkg	rrjpn rraus rridn rrkor rrmys rrphl  rrsgp rrtai rrtha rrarg rrbra rrchl rrmex	rrtur"    
' Create working file
wfcreate(wf=%wfname) m {%obs_beg_all} {%obs_end_all}        
' Create page     
pagecreate(page=%pgname1) m {%obs_beg_all} {%obs_end_all}    
' Read excel data
read(c3,s={%xls_sheet}) %path {%variables_all}   
  
'**************************************************************************************************************************************
'**************************************************************************************************************************************
'  TEST 1 
'*************************************************************************************************************************************
'  To test toolbox functions, code below replicates spillover index table (Table 3)  from the following working paper :
'  { Source:    Economic Journal
'    Date:        January 1, 2009
'    Title:         "Measuring Financial asset return and volatility spillovers, with application to global equity markets"
'    Authors:    Francis X. Diebold and Kamil Yilmaz }

%variables=%variables_all
%egzo_variables=" "                                               ' no exogenous parameters
scalar nr_var=@wcount(%variables)                  ' number of endogenous variables      
scalar nr_egzo=@wcount(%egzo_variables)   '  number of exogenous variables
scalar nr_lag=2                                                       ' number of lags 
scalar nr_egzo_res=0                                           ' number of restrictions
scalar nr_obs=@obs(rrdjia)                                '  number of observations

' IRF (impulse response function) information
' !hzr                     - >   impulse response horizont
' !unit_shocks    ->   for unitshock put one, for shock of one standard deviation put zero
!hzr=10
!unit_shocks=0
' Creating system of equations:
string variables=%variables
string egzo_variables=%egzo_variables
' Creating group of endogenous variables
group  grupa_var
grupa_var.add {%variables}
' Creating group of exogenous variables
group  grupa_egzo
grupa_egzo.add {%egzo_variables}
' Initialization of system parameters (will be estimated later)
coef(5000) c_param=0

%nnaammee1="struct_system"
%nnaammee2="system_res"
system {%nnaammee1}

' egzo_m  ->  matrix of ones if all exogenous variables appear in all equations
'                 ->  set matrix elements to zero to exclude exogenous variables in certain equation
if nr_egzo<>0 then
    matrix(nr_var,nr_egzo) egzo_m=1 
else
    matrix(nr_var,1) egzo_m=0  ' initialization if there is no exogenous variables
endif

' Function create_system creates system of equations
call create_system(grupa_var,grupa_egzo,nr_lag,nr_egzo_res,egzo_m,c_param,{%nnaammee1},%nnaammee2)
' Function var_representation calculates structural, reduced-form and MA representation 
call var_representation(grupa_var,grupa_egzo,nr_lag,nr_egzo_res,egzo_m,!hzr,!unit_shocks,c_param,{%nnaammee2})
' Function spillover_index calculates spillover index and creates spillover table
call spillover_index(%variables,"psi_tilda","spillover_table","spillover_index")

' Move selected objects to output page of working file
pagecreate(page=%pgname2) m {%obs_beg_all} {%obs_end_all}    
%objects_to_move="spillover_table spillover_index"
for %s {%objects_to_move}
	copy {%pgname1}\{%s}  {%pgname2}\{%s} 
next

'**************************************************************************************************************************************************
'**************************************************************************************************************************************************                 
' TEST 2
'*************************************************************************************************************************************************
' This test example has no economic meaning (to illustrate how to include exogenous variables)

' Create new page and copy selected endogenous and exogenous variables 
pagecreate(page=%pgname3) m {%obs_beg_all} {%obs_end_all}   
pageselect { %pgname3}
%variables="rrdjia rrftse rrfra"
%egzo_variables="rrchl rrmex rrtur"
for %s {%variables}
      copy {%pgname1}\{%s}  {%pgname3}\{%s} 
next
for %s {%egzo_variables}
      copy {%pgname1}\{%s}  {%pgname3}\{%s} 
next

scalar nr_var=@wcount(%variables)                  ' number of endogenous variables      
scalar nr_egzo=@wcount(%egzo_variables)   '  number of exogenous variables
scalar nr_lag=2                                                       ' number of lags 
scalar nr_egzo_res=0                                           ' number of restrictions
scalar nr_obs=@obs(rrdjia)                                '  number of observations

' IRF (impulse response function) information
' !hzr                     - >  impulse response horizont
' !unit_shocks    ->   for unitshock put one, for shock of one standard deviation put zero
!hzr=10
!unit_shocks=0

' egzo_m  ->  matrix of ones if all exogenous variables appear in all equations
'                 ->  set matrix elements to zero to exclude exogenous variables in certain equation
if nr_egzo<>0 then
    matrix(nr_var,nr_egzo) egzo_m=1 
else
    matrix(nr_var,1) egzo_m=0  ' initialization if there is no exogenous variables
endif
egzo_m(1,2)=0
egzo_m(3,2)=0
egzo_m(1,3)=0
egzo_m(2,3)=0

' Creating system of equations:
string variables=%variables
string egzo_variables=%egzo_variables
' Creating group of endogenous variables
group  grupa1
grupa1.add {%variables}
' Creating group of exogenous variables
group  grupa2
grupa2.add {%egzo_variables}
' Initialization of system parameters (will be estimated)
coef(5000) c_param=0

%nnaammee1="struct_system"
%nnaammee2="system_res"
system {%nnaammee1}

' Function create_system creates system of equations
call create_system(grupa1,grupa2,nr_lag,nr_egzo_res,egzo_m,c_param,{%nnaammee1},%nnaammee2)
' Function var_representation calculates structural, reduced-form and MA representation 
call var_representation(grupa1,grupa2,nr_lag,nr_egzo_res,egzo_m,!hzr,!unit_shocks,c_param,{%nnaammee2})
' Function hd calculates hystorical decomposition
call hd("psi","other_objects",%variables,%egzo_variables,{%nnaammee2},nr_lag)
' Function irf calculates impulse responses and accumulated impulse responses
call irf("psi_tilda")
' Function spillover_index calculates spillover index and creates spillover table
call spillover_index(%variables,"psi_tilda","spillover_table","spillover_index")

pageselect {%pgname2}

' Stop execution here
stop

' *********************************************************************************************************************************************

'

