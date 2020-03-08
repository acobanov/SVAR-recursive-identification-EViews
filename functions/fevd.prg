subroutine fevd(string %psi)
'************************************************************************************************************************************************
' FORECAST ERROR VARIANCE DECOMPOSITION
'************************************************************************************************************************************************

' Author:    Ante Cobanov (ante.cobanov@yahoo.com)
' GitHub:   acobanov

'*************************************************************************************************************************************************

	'INPUT:
		'%psi --- name of  user object "psi" (output from "var_representation" subroutine)

 	'OUTPUT
		' user object fevd
		' contains n matrices, fevd_1,...,fevd_n
		' fevd_i matrix contains forecast error variance decomposition for i-th variable
		' h-row of matrix is fevd for h-horizont (sum to 1)
		' j-th column contains contributions from j-th shock

'****************************************************************************************************************************************************
%imeime1="fevd_"
%imeime2="fevd_norm_"
%imeime3="fevd"

%members={%psi}.@members
!length=@wcount(%members)
%first_member=@word(%members,1)
{%psi}.extract {%first_member} pomocna
!num_var=@rows(pomocna)
delete pomocna

Userobj  {%imeime3}

' calculate fevd's

for !kkk=1 to !num_var
	matrix(!length-1,!num_var) {%imeime1}{!kkk}
next
for !kkk=1 to !num_var
	matrix(!length-1,!num_var) {%imeime2}{!kkk}
next

for !iii=1 to !num_var
	for !kkk=1 to !length-1
		for !jjj=1 to !num_var
			!sum=0
			for !sss=0 to !kkk-1
				%sssth_member=@word(%members,!sss+1)
				{%psi}.extract {%sssth_member} pomocna
				!pomocni=pomocna(!iii,!jjj)*pomocna(!iii,!jjj)
				!sum=!sum+!pomocni
				delete pomocna
			next
			{%imeime1}{!iii}(!kkk,!jjj)=!sum
		next
		vector pomocni=@rowextract({%imeime1}{!iii},!kkk)
		!sumsum=@sum(pomocni)
		delete pomocni
		for !jjj=1 to !num_var
			{%imeime2}{!iii}(!kkk,!jjj)={%imeime1}{!iii}(!kkk,!jjj)/!sumsum
		next
		
	next
next

for !kkk=1 to !num_var
	{%imeime3}.add(d) {%imeime1}{!kkk}
	{%imeime3}.add(d) {%imeime2}{!kkk}
next

endsub
''********************************************************************************************************************************************************************


