subroutine hd(string %psi,string %other,string %series,string %egzo,group gr_error_str,scalar num_lag)
'*****************************************************************************************************************************************************************************************************'*************************************************

' Author:    Ante Cobanov (ante.cobanov@yahoo.com)
' GitHub:   acobanov

'*******************************************************************************************************************************************************************************************************************************************************

'  Modife info!!! hd_z contribution from exogenous variables added!!!

' 	Historical decomposition calculation requires call of "var_representation" subroutine to calculate MA representation (psi matrices), J1 matrix and G matrix
'	Subroutine calculates hd from time t=1 to time T, t=1 is first time point with non missing observations in group of structural residuals ant T is defined by horizon to calculate psi matrices (change horizont to calculate hd for whole sample)

'    NOTATION:
' 	y(t) = J1 * G^(t) * z(0)  +  sum(i=1 to t-1) ( psi(i)*theta(t-i) )
'	y-part:  sum(i=1 to t-1) (psi(i)*theta(t-i) ) is contribution from shocks 
'	x-part:   J1 * G^(t) * z(0)  is contribution from VAR structure
' 	theta(t):  structural form residuals 
'	y(t) contains !num_var variables, subroutine calculates hd for each variable for each contribution x and y
'	hd_x_i  and  hd_y_i are x and y contributions for ith variable, hd_x_i is matrix with 1 column and hd_y_i is matrix with numbers of columns equal to number of shocks
	
	' INPUT:
		' %psi                 ->     name of user object "psi" from var_representation subroutine
           ' %other             ->     name of user object "other_objects" from var_representation subroutine
		' %series           ->    string list containing endogenous variables in SVAR model in correct order
          ' %egzo              ->     string list containing exogenous variables in SVAR model
		' gr_error_str    ->    group of structural form residuals
		' num_lag         ->     number of lags in SVAR model
		
 	  'OUTPUT
		' "hd" ->  name for user object hd (contains matrices hd_x and matrices hd_y for each varibale)
		' contains matrices: hd_y_1,...,hd_y_n ( y is contribution from structural shocks)
		' contains matrices: hd_x_1,...,hd_x_n ( x is contribution from  VAR structure )
		' hd_x_i  and hd_y_i matrices contains historical decomposition for i-th variable
		' in h-row of matrix there is hd for h time point
		' for matrix hd_y_i    j-th coloumn contains contributions from j-th shock to variable i

'*************************************************************************************************************************************************************************************************************************************************************
%imeime1="hd"
%imeime2="hd_x_"
%imeime3="hd_y_"
%imeime4="hd_z_"
'
{%other}.extract G G
{%other}.extract J1 J1
{%other}.extract H H
'{%other}.extract ex_red ex_red

%members={%psi}.@members
!length=@wcount(%members)
%first_member=@word(%members,1)
{%psi}.extract {%first_member} pomocna
!num_var=@rows(pomocna)
delete pomocna
!num_egzo=@wcount(%egzo)
'calculate hd_y:
for !iii=1 to !num_var
	matrix(!length-1,!num_var) {%imeime3}{!iii}
next
matrix pomocna2
stom(gr_error_str,pomocna2)
for !kkk=1 to !num_var ' for each variable
	for !iii=1 to !length-1 ' for given time point
		for !jjj=1 to !num_var ' for given shock
			!sum=0
			for !sss=0 to !iii-1
				%sssth_member=@word(%members,!sss+1)
				{%psi}.extract {%sssth_member} pomocna1
				!rrr=!iii-!sss
				!pomocni=pomocna1(!kkk,!jjj)*pomocna2(!rrr,!jjj)
				!sum=!sum+!pomocni
				delete pomocna1
			next
			{%imeime3}{!kkk}(!iii,!jjj)=!sum
		next
	next
next

delete pomocna2 

'calculate hd_Z:
for !iii=1 to !num_var
	matrix(!length-1,!num_egzo) {%imeime4}{!iii}
next
matrix pomocna2
group pomocna3
pomocna3.add {%egzo}
stom(pomocna3,pomocna2)
for !kkk=1 to !num_var ' for each variable
	for !iii=1 to !length-1 ' for given time point
		for !jjj=1 to !num_egzo ' for given shock
			!sum=0
			for !sss=0 to !iii-1
				%sssth_member=@word(%members,!sss+1)
				{%psi}.extract {%sssth_member} pomocna1
				!rrr=!iii-!sss
				pomocna1=pomocna1*H
				!pomocni=pomocna1(!kkk,!jjj)*pomocna2(!rrr,!jjj)
				!sum=!sum+!pomocni
				delete pomocna1
			next
			{%imeime4}{!kkk}(!iii,!jjj)=!sum
		next
	next
next

delete pomocna2 pomocna3
'***********************************************************************************************************************************************************************************************************************************************************************************************************************
'calculate hd_x contributions:

for !iii=1 to !num_var
	matrix(!length-1,1) {%imeime2}{!iii}
next

'create group of series in SVAR model and transform group into matrix object
group pomocna1.add {%series}
matrix pomocna2
stom(pomocna1,pomocna2)


'calculate z(0)
matrix(num_lag*!num_var+1,1) z0
z0(1,1)=1
!brojac=-1
for !kkk=num_lag to 1 step -1
	!brojac=!brojac+1
	matrix v=@rowextract(pomocna2,!kkk)
	v=@transpose(v)
	!iii=2+!num_var*!brojac
	matplace(z0,v,!iii,1)
next
delete pomocna1 pomocna2

'calculate contributions:
matrix pomocna1=@identity(!num_var*num_lag+1,!num_var*num_lag+1)
for !iii=1 to !length-1
	pomocna1=pomocna1*G
	matrix pomocna2=J1*pomocna1*z0
	for !kkk=1 to !num_var 
		{%imeime2}{!kkk}(!iii,1)=pomocna2(!kkk,1)
	next
next
'*****************************************************************************************************************************************************************************************************************************************************************************************************************************
Userobj {%imeime1}

for !iii=1 to !num_var
	hd.add(d) {%imeime2}{!iii}
	hd.add(d) {%imeime3}{!iii}
	hd.add(d) {%imeime4}{!iii}
next

endsub
' ****************************************************************************************************************************************************************************************************************************************************************************************************************************************
' *****************************************************************************************************************************************************************************************************************************************************************************************************************************************

'


