subroutine create_system(group grupa_var,group grupa_egzo,scalar num_lag, scalar num_egzo_res,matrix egzo, coef c_param,system struct_system,string system_res)
'*****************************************************************************************************************************************************************************************************
' Create and estimate system of equations
'***************************************************************************************************************************************************************************************************

' Author:    Ante Cobanov (ante.cobanov@yahoo.com)
' GitHub:   acobanov

'*******************************************************************************************************************************************************************
	' INPUT:
     ' grupa_var             ->   group of endogenous variables
     ' grupa_egzo          ->   group of exogenous variables
     ' num_lag               ->   number of lags
     ' num_egzo_res   ->    number of block restrictions 
	' egzo                      ->   matrix to put zero restriction for exogenous variables
	' length                   ->   number of periods (for computing psi matrices for MA representation)
	' c_param              ->   initialize series of estimated parameters  
     ' struct_system     ->   initilazite system of equations
	' sistem_res         ->   name of group of structural form residuals

	' OUTPUT:
     ' c_param              ->   final series of estimated parameters  
	' sistem_res         ->   group of estimated structural form residuals

'*************************************************************************************************************************************************************
string data_var=grupa_var.@members
string data_egzo=grupa_egzo.@members

'SVAR in form suitable for EVIews equation by equation estimation

'	number of endogenous variables in SVAR model (num_var)
'	number of exogenous variables in SVAR model (num_egzo)
!num_var=@wcount(data_var)
!num_egzo=@wcount(data_egzo)

'	vector of estimated parameters (c_param)
'coef(5000) c_param=0

'	start value for number of coefficients:
!nr_coeff = 0

' creating string for equations ( str_eq_1, str_eq_2, itd ) and put left-hend side of equations in string
for !iii=1 to !num_var
	%str_name = "str_eq_" + @str(!iii)
	string {%str_name}	
	{%str_name} = @word(data_var, !iii)+"="
next

' contemporaneous effects
for !iii=2 to !num_var
	%str_name = "str_eq_" + @str(!iii)
	for !jjj=1 to !iii-1
		!nr_coeff = !nr_coeff + 1
		{%str_name} = {%str_name}  + "c_param(" + @str(!nr_coeff) + ")*" + @word(data_var, !jjj) + "+"
	next
next

' adding constant
for !iii=1 to !num_var
	%str_name = "str_eq_" + @str(!iii)
	!nr_coeff = !nr_coeff + 1
	{%str_name} = {%str_name}  + "c_param(" + @str(!nr_coeff) + ")+"
next

' adding lags:
for !kkk=1 to num_lag
	for !jjj=1 to num_egzo_res 'for first num_egzo_res variables (block with zeros)
		%str_name = "str_eq_" + @str(!jjj)		
		for !iii=1 to num_egzo_res 'only first num_egzo_res variables have non-zero coefficients
			!nr_coeff = !nr_coeff + 1	
			{%str_name} = {%str_name}  +  "c_param(" + @str(!nr_coeff) + ")* " + @word(data_var, !iii) + "(-" + @str(!kkk) + ")" + "+"
		next
	next

	for !jjj=num_egzo_res+1 to !num_var
		%str_name = "str_eq_" + @str(!jjj)	
		for !iii=1 to !num_var
			!nr_coeff = !nr_coeff + 1			
			{%str_name} = {%str_name}  +  "c_param(" + @str(!nr_coeff) + ")* " + @word(data_var, !iii) + "(-" + @str(!kkk) + ")" + "+"
		next
	next
next

'adding exogenous variables (egzo matrix defines equations for which exogenous variables appears)
if  @isempty(data_egzo)=0 then 
	for !kkk=1 to !num_egzo
		for !iii=1 to !num_var
			%str_name = "str_eq_" + @str(!iii)	
			if egzo(!iii,!kkk)<>0 then
				!nr_coeff = !nr_coeff + 1	
				{%str_name} = {%str_name}  +  "c_param(" + @str(!nr_coeff) + ")* " + @word(data_egzo, !kkk)+ "+"
			endif
		next
	next
endif		

' removing + from the end
for !iii=1 to !num_var
	%str_name =  "str_eq_" + @str(!iii)
	!str_len= @length({%str_name})
	{%str_name} = @left({%str_name},!str_len-1)
next

'********************************************************************************************************************************************
''Creating system of equations:
'%nnaammee1="struct_system"
'system {%nnaammee1}
' system creation
for !iii=1 to !num_var
	%str_name =  "str_eq_" + @str(!iii)
	struct_system.append {{%str_name}}
	struct_system.append "   "
next

'Estimate system of equations:
struct_system.ls 

%ime=system_res
'make residuals from equations:
'struct_system.makeresids(n={%ime})
struct_system.makeresids(n={%ime})
'system_res={%ime}


'delete equations and system:
for !iii=1 to !num_var
	%str_eq = "str_eq_" + @str(!iii)
	delete {%str_eq}	
next
'   ' delete {%nnaammee1}

endsub
'*******************************************************************************************************************************************


