program Aging


!implicit none
Integer, Parameter      :: nkmax = 2**12, ntc=10, ntmax=1000, nt2max=1000000, nrmax=250
Real(8), Parameter      ::    R=1.d0, alfa=1.305d0 
Real(8), Parameter      ::  D=1.d0
Real(8)			::pi, fact, a, b, etha, rho, k,gama, kc, suma, tol, tr1, ethavw,lamdak, rhovw, du, tn, dk, dr, jj, u, u0
Real(8)			:: ethamax, errorF, tolF,  h, t, suma1, suma2, Dl,large, error, ethamin, Friccion, ymax, ymin, suma7,umax,umin
Real(8), Dimension(nkmax) :: Sb,landa, landaI, SbI, kap, kaps, gamalist
Real(8), Dimension(ntmax) :: Dz, bdt, wm2
Real(8), Dimension(nrmax) ::rd,gdr
Real(8), Dimension(nt2max)	:: t2, Fs1, bdt2
Real(8)			:: Dt, sumacri,sumagr
Real(8), Dimension(nkmax,ntmax)	:: Fn, Fs
Complex(8), Parameter	:: I = Cmplx(0.d0,1.d0)
Integer		        :: flag, nk, it, nt, l, decim, ml, tim1, mlm, npi, my1, my2, my3, my4, my5, my6, my7,ii
Real(8), Dimension(nkmax) :: rk,Ski, Skf, Sk0, alpha, kd, rk2
 character*4 sub


!npi=32


pi=4.d0*atan(1.)

fact =1.d0/(6.d0*(pi**2))



etha=0.63d0! + (lx-1)*dx0

rho = (6.0d0*etha)/(pi*R**3)


!stop

!open(5,file='Gdraging'//subx//'.dat',status="Unknown")
open(9,file='Skaging.dat',status="Unknown")
open(10,file="S(k)_ini.dat")			!%%%% phi0=0.30, T0=10.d0 %%%%
open(11,file='S(k)_fin.dat')

do nk=1,nkmax
read(10,*) rk(nk),Ski(nk)
read(11,*) rk2(nk),Skf(nk)

kd(nk)=rk(nk)


alpha(nk)=2.d0*(kd(nk)**2)*D/Skf(nk)					!definición de la alfa(k)


enddo
close(10)
close(11)


dk = kd(2)-kd(1)


!stop


! !!!!!!!!!
 Umin=0.0d0
 umax=1.d0
 u= (umin +umax)/2.d0
!!!!!!!!!1

 tolF=1.d-6
 errorF=1.d0
  
 do while(errorF >= tolF)
 
 

do nk=1,nkmax
					
Sk0(nk)=Ski(nk)*exp(-alpha(nk)*u) + Skf(nk)*(1.d0 - exp(-alpha(nk)*u))			!definicion del factor de estructura dependiente del tiempo

 Sb(nk)=Sk0(nk)			!factores de estructura que le entran al programa dinámico
 SbI(nk)=1.d0/Sb(nk)

 call first_minimum( ymax, ymin) 

enddo


kc=2.d0*pi*alfa


!   definicion de la funcion interpoladora lamda %%%%%%
do nk=1,nkmax

landa(nk)= 1.d0/(1.d0 + (kd(nk)/kc)**2)
landaI(nk)= 1.d0 + (kd(nk)/kc)**2

enddo

call evalua_flag

	print*, "flag=",flag, u



 if( flag==-1) then
 
 umax = u
 
   else
    umin = u
    
  end if

	errorF= abs((((umin + umax)/2.d0) - u)/u)

	u= (umin +umax)/2.d0
enddo

write(*,*)'la u de arresto es',sngl(u), "gamma=",sngl(gama), "flag=",flag

do nk=1,nkmax
					
Sk0(nk)=Ski(nk)*exp(-alpha(nk)*u) + Skf(nk)*(1.d0 - exp(-alpha(nk)*u))			!definicion del factor de estructura dependiente del tiempo

Sb(nk)=Sk0(nk)			!factores de estructura que le entran al programa dinámico

write(9,*) kd(nk), Sb(nk)
enddo




contains




Subroutine Indice(ml,sub)      

     Implicit Integer (i-n)
     Implicit Real (A-H,O-Z)
     Character*10 Doschr
     Character*4 sub
     Character*1 Str1
     Character*2 Str2
     Character*3 Str3
     Character*4 Str4
     Character*5 Str5
     
     If (ml.Lt.10) Then
        Write(Str1,933)ml
        sub=Str1
     Elseif (ml.Lt.100) Then  
        Write(Str2,935)ml
        sub=Str2
     Elseif (ml.Lt.1000) Then
        Write(Str3,937)ml
        sub=Str3
     Elseif (ml.Lt.10000) Then
        Write(Str4,939)ml
        Sub=Str4
     Elseif (ml.Lt.100000) Then
        Write(Str5,940)ml
        sub=Str5
     Endif

933  Format(I1)
935  Format(I2)
937  Format(I3)
939  Format(I4)
940  Format(I5)
     Return
End subroutine




! Subroutine Indice(mlx,subx)      
! 
!      Implicit Integer (i-n)
!      Implicit Real (A-H,O-Z)
!      Character*10 Doschr
!      Character*4 subx
!      Character*1 Str1
!      Character*2 Str2
!      Character*3 Str3
!      Character*4 Str4
!      Character*5 Str5
!      
!      If (mlx.Lt.10) Then
!         Write(Str1,933)mlx
!         subx=Str1
!      Elseif (mlx.Lt.100) Then  
!         Write(Str2,935)mlx
!         subx=Str2
!      Elseif (mlx.Lt.1000) Then
!         Write(Str3,937)mlx
!         subx=Str3
!      Elseif (mlx.Lt.10000) Then
!         Write(Str4,939)mlx
!         Subx=Str4
!      Elseif (mlx.Lt.100000) Then
!         Write(Str5,940)mlx
!         subx=Str5
!      Endif
! 
! 933  Format(I1)
! 935  Format(I2)
! 937  Format(I3)
! 939  Format(I4)
! 940  Format(I5)
!      Return
! End subroutine





subroutine first_minimum(ymax, ymin)		!calcula el maximo y el primer minimo despues del maximo del factor de estructura	lyr
	Implicit none
	
	Double precision  ymax, ymin
	integer nk, imax, imin
	Double precision Smax, Smin
        
       
	Smax= sb(1)						!aqui se calcula el maximo del factor de estructura
	do nk=1, nkmax
	If (sb(nk).ge.Smax) then 
	Smax=Sb(nk)
	imax=nk

	end if
	end do
   	ymax= kd(imax)

	Smin= Sb(imax)
	do nk=imax, nkmax					!aqui calculamos el minimo del factor de estructura :p
	If (Sb(nk).le.Smin) then 
	Smin=Sb(nk)
	imin=nk
	end if
	end do
   	ymin= kd(imin)

	return 
end  subroutine first_minimum





subroutine evalua_flag				
integer	:: nk, it
real(8)	:: k, suma

tol=1.d-6
large=1.d+6
gama=1.d-8
it=0
error=1.d0

do while ( error > tol  .and. gama < large)
	it= it + 1

	sumacri= 0.0d0

	do nk=1, nkmax

				
		sumacri= sumacri + ((kd(nk))**4)*(((Sb(nk) - 1.d0)**2)*((landa(nk))**2))/(((landa(nk))*(Sb(nk)) + &
		((kd(nk))**2)*gama)*(landa(nk) + ((kd(nk))**2)*gama))

	end do

	sumacri=(fact/rho)*sumacri*dk

	error= abs(((1.d0/sumacri) - gama)/gama)

	gama=1.d0/sumacri

	gamalist(it)= gama

!print*, "it, gama =", it, gama

end do


	flag =sign(1.d0, gamalist(it) - 2.d0*gamalist(it - 1) + gamalist(it - 2))

end subroutine evalua_flag


end program Aging
