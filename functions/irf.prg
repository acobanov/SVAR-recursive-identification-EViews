subroutine irf(string %ime1)
'*************************************************************************************************************************************
' IMPULSE RESPONSE FUNCTION 
'************************************************************************************************************************************

' Author:    Ante Cobanov (ante.cobanov@yahoo.com)
' GitHub:   acobanov

'*************************************************************************************************************************************

' INPUT:      %ime1                        ->     name of user-define object which contains "psi_tilda" matrices from MA representation

' OUTPUT:   irf_matrix                   ->     matrix of all impulse responses of all variables
			   '  irf_matrix_accum    ->      matrix of all accumulated impulse responses of all variables

'*************************************************************************************************************************************
%members={%ime1}.@members
!length=@wcount(%members)
%first_member=@word(%members,1)
{%ime1}.extract {%first_member} pomocna
!num_var=@rows(pomocna)
delete pomocna

%imeime3="irf_matrix"
%imeime4="irf_matrix_accum"
%imeime5="irf"

 'Impulse responses: 
     matrix(!length, !num_var*!num_var)  {%imeime3}=na
     for !kkk=0 to !length-1 ' periods
			%kth_member_name=@word(%members,!kkk+1)
			{%ime1}.extract {%kth_member_name} pomocna
            	for !iii = 1 to !num_var' variables
              		  for !jjj=1 to !num_var'variables
               	    		  {%imeime3}(!kkk+1, ((!iii-1)*!num_var)+!jjj) =pomocna(!iii, !jjj) ' response{!k}(!i, !j)= response of i-th variable to shock in j-th variable after !k periods
                next
           next
		delete pomocna
     next

' Accumulated impulse responses:
	matrix(!length, !num_var*!num_var)  {%imeime4}={%imeime3}
	 for !kkk = 2 to !length
                for !jjj = 1 to (!num_var^2)
                     {%imeime4}(!kkk,!jjj) = {%imeime4}(!kkk-1,!jjj) + {%imeime4}(!kkk,!jjj) 
                next
       next

Userobj {%imeime5}
for !kkk=1 to !num_var
	!jjj1=1+(!kkk-1)*!num_var
	!iii2=!length
	!jjj2=!jjj1+!num_var-1
	matrix(!length,!num_var) irf_{!kkk}=@subextract({%imeime3},1,!jjj1,!iii2,!jjj2)
	matrix(!length,!num_var) irf_accum_{!kkk}=@subextract({%imeime4},1,!jjj1,!iii2,!jjj2)
next

for !kkk=1 to !num_var
	{%imeime5}.add(d)  irf_{!kkk}
	{%imeime5}.add(d)  irf_accum_{!kkk}
next 

endsub

'*************************************************************************************************************************************************************************************************************************************************************************************************************************
'*************************************************************************************************************************************************************************************************************************************************************************************************************************


