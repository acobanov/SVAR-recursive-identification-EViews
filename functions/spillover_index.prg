subroutine  spillover_index(string %list,string %object,string %table,string %index)
'***********************************************************************************************************************************************************************************
' Spillover index and spillover table
'***********************************************************************************************************************************************************************************

' Author:    Ante Cobanov (ante.cobanov@yahoo.com)
' GitHub:   acobanov

'**************************************************************************************************************************************************************************************************
	'INPUT:
	'string %list          ->  string list of variables in SVAR model (be careful with ordering!!! ) (we need list of variables to put names of series and shocks into final spillover table)
	'string  %object   ->  name of user defined object  (psi_tilda matrices from MA representation of SVAR model -- subroutine "var_representation" creates  user object with psi_tilda matrices)
	'string %table      ->  string name for output spillover table 
	'string %index     ->  string name for scalar number that represents overall spillover index
		
	'OUTPUT:
	' table {%table}      ->  spillover table (table contains all contributions from and to others)
	' scalar {%index}   ->  overall spillover index
	
'*************************************************************************************************************************************************************************************
	'string list that contains all names of matrices collected in object {%object}
	%object_members={%object}.@members
	'number of variables in model:
	!num_var=@wcount(%list)

	'	Caluclating spillover index:
	matrix(!num_var,!num_var) spillover=0
	for !iii=1 to !num_var
		!suma2=0

		for %k {%object_members} 
			{%object}.extract {%k} pomocna
			for !jjj=1 to !num_var
				!suma2=!suma2+pomocna(!iii,!jjj)^2
			next
			delete pomocna
		next
         

		for !jjj=1 to !num_var
			!suma1=0
			for %k {%object_members} 
				{%object}.extract {%k} pomocna
				!suma1=!suma1+pomocna(!iii,!jjj)^2
				delete pomocna
			next	
			spillover(!iii,!jjj)=!suma1/!suma2*100
		next
	next

	'Creating spillover table:
       %table_name=%table
       table {%table_name}
       '    Uèitavanje prvog polja tablice:
             {%table_name}(1,1)=  %table_name
        ' Uèitavanje prvog stupca tablice:
           for %s {%list}
		!iii=@wfind(%list,%s)
           {%table_name}(!iii+1,1)=%s
           next
        '   Uèitavanje prvog reda tablice:
         for %s {%list}
		!jjj=@wfind(%list,%s)
        {%table_name}(1,!jjj+1)=%s
        next
	'	Uèitavanje matrice procjena  u tablicu:
		for !iii=1 to !num_var
			for !jjj=1 to !num_var
				{%table_name}(!iii+1,!jjj+1)=spillover(!iii,!jjj)
			next
		next
'************************************************************************************************
	  {%table_name}(1,!num_var+2)="contribution from others"
		!suma2=0
     for !iii=1 to !num_var
		!suma1=0
		for !jjj=1 to !num_var
			if !iii<>!jjj then
				!suma1=!suma1+@val({%table_name}(!iii+1,!jjj+1))
			endif
            next
		  {%table_name}(!iii+1,!num_var+2)=!suma1
		  !suma2=!suma2+!suma1
	next

'****************************************************************************************
	{%table_name}(!num_var+2,!num_var+1)="spillover index: "
	{%table_name}(!num_var+2,!num_var+2)=!suma2/!num_var

	scalar {%index}
	{%index}=!suma2/!num_var
'*********************************************************************************************

 {%table_name}(!num_var+2,1)="contribution to others"
	for !jjj=1 to !num_var
		!suma1=0
		for !iii=1 to !num_var
			if !iii<>!jjj then
				!suma1=!suma1+@val({%table_name}(!iii+1,!jjj+1))
			endif
            next
		 {%table_name}(!num_var+2,!jjj+1)=!suma1
	next

	delete spillover 

endsub

'*************************************************************************************************************************************************************************************************************************************************************************************************************************
'*************************************************************************************************************************************************************************************************************************************************************************************************************************

'


