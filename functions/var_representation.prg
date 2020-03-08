subroutine var_representation(group grupa_var,group grupa_egzo,scalar num_lag, scalar num_egzo_res,matrix egzo,scalar length,scalar shocks,coef c_param,group sistem_res)
'*****************************************************************************************************************************************************************************************************

' Author:    Ante Cobanov (ante.cobanov@yahoo.com)
' GitHub:   acobanov

'*****************************************************************************************************************************************************************************************************
	' NOTATION: 
	' Structural representation :   D0 * y(t) = const_structural + D1 * y(t-1) + D2 * y(t-2) +  ...  + DP * y(t-p) +  E1 * x1(t)  +  ... + Em *x m(t)  + error_structural  
	'                                                  (D0 is unit lower triangular matrix, unit elements on the main diagonal;  x1(t) to xm(t) are m exogenous variables)

	' Reduced representation :   y(t) = D0^(-1) * const_structural + D0^(-1) * D1 * y(t-1) + D0^(-1) *D2 * y(t-2) +  ...  + D0^(-1) * DP * y(t-p) + D0^(-1) * E1 * x1(t)  +  ... + D0^(-1) * E2 * xm(t) +  D0^(-1) * error_structural 
	'                                                  ( i.e. y(t) = const_reduced + A1 * y(t-1)+A2*y(t-2)+...+AP*y(t-p)+E1_red*x1(t) + ...+Em_red * xm(t)  + error_reduced

	' Notation:  D_matrica = D0^(-1);  error_reduced = D_matrica * error_structural ; cov (error_structural) = gamma_matrix  (gamma_matrix is diagonal  with variances of residuals on the main diagonal)
	'                   cov(error_structural)  = B * cov(I) * B'  ("I" is identity matrix) so that B needs to be diagonal with st.deviations on the main diagonal
	'                   If  B is equal to identity matrix, it will be propagation of unit shock

	' VAR(1) representation:   Y(t) = const + A * Y(t-1) + error or z(t) = G * z(t-1) + error

	' MA representation:  y(t) = const + sum(i=0 to infinity) ( phi(i)*u(t-i) )  = const + sum(i=0 to infinity) ( psi(i)*v(t-i) ) = const + sum(i=0 to infinity) ( psi_tilda(i)v'(t-i) ) 

	' u(t) - >   reduced form shocks
     '  v(t)   ->  structural shocks (covariance matrix diagonal)
     '  v'(t)  ->  structural shocks (covariance matrix identity matrix)

'************************************************************************************************************************************************************************************
	' INPUT:
     ' grupa_var              ->   group of endogenous variables
     ' grupa_egzo           ->   group of exogenous variables
     ' num_lag                ->   number of lags
     ' num_egzo_res    ->    number of block restrictions 
	' egzo                       ->   matrix to put zero restriction for exogenous variables
	' length                    ->   number of periods (for computing psi matrices for MA representation)
	' c_param               ->   series of estimated parameters 
     ' shocks                  ->   put "0" for shocks of one standard deviation or "1" for unity shocks
	' sistem_res          ->   group of structural form residuals

	' OUTPUT:
	' output is organized with user defined objects:
	' user object "red" for reduced - form representation     (contains matrices:  const_reduced    and   A1, A2, ..., AP)
	' user object "struct" for structural representation           (contains matrices:  D0,  const_structural   and   D0, D1, D2, ..., DP )
	' user object "phi"  for phi matrices                                    (contains matrices: phi0, phi1, .. ., phi%length-1%)
	' user object "psi"  for psi matrices                                    (contains matrices: psi0, psi1, .. ., psi%length-1%)
	' user object "psi_tilda"  for psi_tilda matrices                (contains matrices: psi_tilda0, psi_tilda1, .. ., psi_tilda%length-1%)
	' user object "other objects"                                                (contains matrices for VAR(1) representation: A, G  matrices: B, J, J1, D_matrica (inverse of matrix D0)

'******************************************************************************************************************************************************************************
%phi="phi"
%psi="psi"
%psi_tilda="psi_tilda"
%reduced="red"
%structural="struct"
%other_objects="other_objects"

'********************************************************************************************************************************************************************************
' SVAR model in form suitable for EVIews estimation equation by equation 

string data_var=grupa_var.@members
string data_egzo=grupa_egzo.@members

' number of endogenous variables in SVAR model (!num_var)
!num_var=@wcount(data_var)
' number of exogenous variables in SVAR model (!num_egzo)
!num_egzo=@wcount(data_egzo)

'*************************************************************************************************************************************************************************************************************************************************************************
' we need c_param (vector of estimated coefficients to fill matrices) and group of estimated residuals (contained in object {%nnaammee2} )

'	initialize matrix B
' 	if we want unit shocks put diag(B) to identity, if we want std dev shocks put diag(B) to std devs
if (shocks = 0) then
	matrix(!num_var, !num_var) B=0
	for !iii=1 to !num_var
		'	Calculating standard deviations of structural form residual series:
		%name_of_series=sistem_res.@seriesname(!iii)
		B(!iii,!iii)=@stdev({%name_of_series})
	next
else '		shocks=1
	matrix B=@identity(!num_var)
endif
'*******************************************************************************************************************************************************************************
' Structural representation ( D0*y(t)=const_structural+D1*y(t-1)+D2*y(t-2)+ ... +DP*y(t-p)+ex_structural* egzo(t) + error_structural)
' Matrices: D0, D1, D2, ..., DP , const_structural, ex_structural
	
matrix(!num_var, !num_var) D0
matrix(!num_var, 1) const_struct
for !kkk=1 to num_lag
	matrix(!num_var, !num_var) D{!kkk}
next
if !num_egzo<>0 then
	matrix(!num_var,!num_egzo) ex_struct
endif

' Coefficients counter:
!br=0
' define matrix D0 in terms of estimated "C" parameters
for !kkk=1 to !num_var
	D0(!kkk,!kkk) = 1
next
for !iii=2 to !num_var
	for !jjj=1 to !iii-1
		!br = !br+1
		D0(!iii,!jjj) = -c_param(!br)
	next
next

' define vector const_structural in terms of estimated  "C" parameters
for !iii=1 to !num_var
	!br = !br +1
	const_struct(!iii,1) = c_param(!br)
next

' define matrices D1, D2, ..., Dp in terms of estimated  "C" parameters
for !kkk=1 to num_lag
	%matrix_name = "D" + @str(!kkk)	
	for !iii=1 to num_egzo_res
		for !jjj=1 to num_egzo_res
			!br =!br +1
			{%matrix_name}(!iii,!jjj) = c_param(!br)
		next
	next
	for !iii=num_egzo_res+1 to !num_var
		for !jjj=1 to !num_var
			!br = !br + 1
			{%matrix_name}(!iii, !jjj) = c_param(!br)
		next	
	next
next

' define matrix ex_structural in terms of estimated  "C" parameters
for !kkk=1 to !num_egzo
	for !iii=1 to !num_var
		if egzo(!iii,!kkk)=0 then
				ex_struct(!iii,!kkk)=0
		else
			!br=!br+1
			ex_struct(!iii,!kkk)=c_param(!br)
		endif
	next
next
	
'*******************************************************************************************************************************************************************************
''    USES  VAR(1) REPRESENTATION OF VAR(P) MODEL to calculate MA representation 

'	Calculating inverse of D0:
matrix D_matrica = @inverse(D0)

'	Calculating matrices, A1=D0^(-1)*D1,  A2=D0^(-1)*D2,  ...  ,  AP=D0^(-1)*DP
for !kkk= 1 to num_lag
	matrix A{!kkk} = D_matrica*D{!kkk}
next
	
'	Calculating matrix A, from matrices A1, A2, ... , AP ( VAR(1) representation )
matrix(!num_var*num_lag,!num_var*num_lag) A=0

for !kkk=1 to num_lag
	!jjj=1+!num_var*(!kkk-1)
	'insert matrix Ai into matrix A at row 1 and column !j
	matplace(A,A{!kkk},1,!jjj)
next

if num_lag>1 then
	!kkk=!num_var*(num_lag-1) 'dimension of big identity matrix to insert into matrix A if number of lags is greater than one
	!rrr=!num_var*(num_lag-1) 'number of rows of zero matrix to insert into matrix A if number of lags is greater than one
	!ccc=!num_var  'number of columns of zero matrix to insert into matrix A if number of lags is greater than one
	matrix insert1=@identity(!kkk)
	matrix(!rrr,!ccc) insert2=0
	!iii=!num_var+1
	matplace(A,insert1,!iii,1)
	!iii=!num_var+1
	!jjj=(num_lag-1)*!num_var+1
	matplace(A,insert2,!iii,!jjj)
	delete insert1 insert2
endif
'********************************************************************************************************************************************************************************************************************************************************************************************************************************************************
'Calculate matrices J,J1 and G
matrix(!num_var,!num_var*num_lag) J=0
matrix(!num_var,!num_var*num_lag+1) J1=0
matrix insert=@identity(!num_var)
matplace(J,insert,1,1)
matplace(J1,insert,1,2)
delete insert

'Calculate matrices J1_tilda
if !num_egzo<>0 then
	matrix(!num_egzo,!num_egzo*num_lag) J_tilda=0
	matrix(!num_egzo,!num_egzo*num_lag+1) J1_tilda=0
	matrix insert=@identity(!num_egzo)
	matplace(J_tilda,insert,1,1)
	matplace(J1_tilda,insert,1,2)
	delete insert
endif
matrix(!num_var*num_lag+1,!num_var*num_lag+1) G=0
G(1,1)=1
matrix(!num_var, 1) const_red=D_matrica*const_struct
matplace(G,const_red,2,1)
matplace(G,A,2,2)
'*******************************************************************************************************************************************************************************
if !num_egzo<>0 then
	matrix(!num_var,!num_egzo) ex_red=D_matrica*ex_struct
	matrix(!num_var*num_lag,!num_egzo) ex_red_big=0
	matplace(ex_red_big,ex_red,1,1)
	matrix(!num_var*num_lag+1,!num_egzo) ex_red_big1=0
	matplace(ex_red_big1,ex_red,2,1)
endif

if !num_egzo<>0 then
	matrix J1_tilda_tr=@transpose(J1_tilda)
	matrix  H=D0*ex_red
endif
'Calculate  phi, psi, psi_tilda matrices:
matrix phi0=@identity(!num_var)
matrix pomocna=@identity(!num_var*num_lag)
matrix J_transpose=@transpose(J)
for !iii=1 to length-1
	pomocna=pomocna*A
	matrix  phi{!iii} = J*pomocna*J_transpose
next
delete pomocna J_transpose
for !iii=0 to length-1
	matrix  psi{!iii} =phi{!iii}*D_matrica
next
for !iii=0 to length-1
	matrix  psi_tilda{!iii} =psi{!iii}*B
next
'*************************************************************************************************************************************************************************************************************************************************************************************************	
'Creating user object "psi" for matrices psi0,psi1,...,psi%length-1%
Userobj  {%phi}
for !iii=0 to length-1
	{%phi}.add(d) phi{!iii}
next

Userobj  {%psi}
for !iii=0 to length-1
	 {%psi}.add(d) psi{!iii}
next
Userobj {%psi_tilda}
for !iii=0 to length-1
	 {%psi_tilda}.add(d) psi_tilda{!iii}
next

'Creating user object for reduced form representation
Userobj {%reduced}
{%reduced}.add(d) const_red
for !iii=1 to num_lag
	{%reduced}.add(d) A{!iii}
next 

if !num_egzo<>0 then
	'{%reduced}.add(d) ex_red
endif

'Creating user object for structural form representation
Userobj {%structural}
{%structural}.add(d) D0
{%structural}.add(d) const_struct
for !iii=1 to num_lag 
	{%structural}.add(d) D{!iii}
next 

if !num_egzo<>0 then
	{%structural}.add(d) ex_struct
endif

if !num_egzo<>0 then
	Userobj {%other_objects}
	%pomocni_string="B D_matrica A J J1 G H ex_red_big ex_red_big1 "
	for %s {%pomocni_string}
 	{%other_objects}.add(d) {%s}
	next
else 
    Userobj {%other_objects}
	%pomocni_string="B D_matrica A J J1 G"
	for %s {%pomocni_string}
 	{%other_objects}.add(d) {%s}
     next
endif


endsub
'*************************************************************************************************************************************************************************************************************************************************************************************************************************

