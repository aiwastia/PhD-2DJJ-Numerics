program PD

! calculate spectrum using Friedel sum rule (see Eq. (4) in PRL 107, 196804 (2011))

use io
use selectedsubr
use scans

implicit none
!Parameters (input)
real*8		:: Deltar,Deltai,dL0,Deltaphase,bmag,btheta,mu0,SCmass
real*8		:: Deltar1,Deltai1,dL1,Deltaphase1,bmag1,btheta1,mu01
real*8		:: gammatot,bmagx,kx,kxmax,kxmin,kFmax,dkx,gap,phi
integer		:: Opnbr,Opnbr1,Nseg,Nseg1,nkxmax,Ndata
character(len=200)	:: varphase

real		::tarray(2), result, tottime
integer		:: n, i, jj
real*8,dimension(5)	:: vardata
character(len=3) 	:: tempstr

character(len=4),dimension(6) :: vartitle
integer :: Ai
real*8,allocatable,dimension(:) :: Aphiarray,ABarray,AmuSCarray,Amuarray,ASCmarray

integer :: g,gg,ggg,gggg,ggggg, g2,lgth
real*8 :: Aloopparam(6)
logical :: testprint



tottime=0

!!!!!!!!!! PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!! system !!!!!!!!!!!!!
call get_infile_parameter

Deltaphase=0d0 !SC phase (units of 2 pi)
Deltaphase1=0d0
btheta=0d0 !magnetic phase (units of 2 pi)
btheta1=0d0
gammatot = 0 !disorder amplitude
varphase='scLR' !specify the system (SC on left and right)
Deltar = Deltamag * cos(Deltaphase*pi) !SC parameters
Deltai = Deltamag * sin(Deltaphase*pi)
Deltar1 = Deltamag1 * cos(Deltaphase1*pi)
Deltai1 = Deltamag1 * sin(Deltaphase1*pi)
!Space discretization
 if (Lmax.ge.dL) then !two leads
     Opnbr=min(30,int(log(Lmax/dL)/log(2d0))+1)
     Nseg = 2**Opnbr
     dL0 = Lmax / dble(Nseg)
 else
     dL0=dL
     Nseg=0
 endif 
 
 if (Lmax1.gt.dLJ) then !junction
     Opnbr1=min(30,int(log(Lmax1/dLJ)/log(2d0))+1)
     Nseg1 = 2**Opnbr1
     dL1 = Lmax1 / dble(Nseg1)
 else
     dL1=dLJ
     Nseg1=0
 endif


!!!!!!!!!! loops !!!!!!!!!!!!!
call get_infile_newparameters

open(unit=30,file=dataname,form='formatted',status='new')

vartitle=(/ '#mu ','SCm ','muSC','B   ','phi ','gap ' /)
write(30,"(6(a,a,a))") (vartitle(i),achar(9),achar(9),i=1,6)

call makelist(trim(Aphilist),Aphiarray,'phi')
Anstepphi=size(Aphiarray)
call makelist(trim(ABlist),ABarray,'B')
AnstepB=size(ABarray)
call makelist(trim(AmuSClist),AmuSCarray,'muSC')
AnstepmuSC=size(AmuSCarray)
call makelist(trim(Amulist),Amuarray,'mu')
Anstepmu=size(Amuarray)
call makelist(trim(ASCmlist),ASCmarray,'SCmass')
AnstepSCm=size(ASCmarray)


!!!!!!!!!! LOOPS TO COMPUTE THE GAP USING THE SCATTERING CODE !!!!!!!!!!!!!!!!!!!!

testprint=.True.
!if change parameter order here, change in Aloopparam, vartitle and write(30 '') at the end
do gggg=1,Anstepmu
	Aloopparam(1)=Amuarray(gggg)
	mu01=Amuarray(gggg)
do ggggg=1,AnstepSCm
	Aloopparam(2)=ASCmarray(ggggg)
	SCmass=ASCmarray(ggggg)
do ggg=1,AnstepmuSC
	Aloopparam(3)=AmuSCarray(ggg)
	mu0=AmuSCarray(ggg)
do gg=1,AnstepB
	Aloopparam(4)=ABarray(gg)
	bmag=ABarray(gg)
do g=1,Anstepphi
	Aloopparam(5)=Aphiarray(g)
	phi=Aphiarray(g)


bmag1=bmag !set same magnetic field in lead and junction
bmag=0. !set to 0 in leads
if(testprint) then
if (bmag.eq.bmag1) then
	print*, 'same magnetic field in the SC leads and in the junction'
else
	print*, 'no magnetic field in the SC leads'
endif
endif

kFmax=mass*alpha1+sqrt(2d0*mass*(mu01+bmag1)+mass**2*alpha1**2) !in junction
if(testprint) print*, 'kFmax= ', kFmax
kxmax=1.1*kFmax
kxmin=0.
dkx=kFmax/(kFmax*Lmax1)**2 / dkfrac
nkxmax=int((kxmax-kxmin)/dkx)


!!!!!!!!!!!!!!!!!!! VERTICAL SECANTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!call vertical_secante(Ndata,epsmin,epsmax,&
!			gammatot,varphase,potshift,phi,&
!			mu0,Deltar,Deltai,alpha,bmag,btheta,SCmass,Opnbr,dL0,&
!			mu01,Deltar1,Deltai1,alpha1,bmag1,btheta1,mass,Opnbr1,dL1,&
!			gap,kxminbuff,kxmaxbuff,secant)


!!!!!!!!!!!!!!!!!!! VERTICAL DYNAMICAL SECANTE !!!!!!!!!!!!!!!!!!!!!!
call dynamical_secante(nkxmax,kxmin,kxmax,&
            		gammatot,varphase,potshift,phi,&
			mu0,Deltar,Deltai,alpha,bmag,btheta,SCmass,Opnbr,dL0,&
			mu01,Deltar1,Deltai1,alpha1,bmag1,btheta1,mass,Opnbr1,dL1,&
            		gap,epsmin,epsmax,secant,5) !!can ask for the maximal depth, will stop at the last relevant


!sum up the running time of each loop
call dtime(tarray, result)
tottime=tottime+result

	Aloopparam(6)=gap
	write(30,'(F8.4,a,F8.4,a,F8.3,a,F8.4,a,F8.6,a,F15.12,a)') (Aloopparam(g2),achar(9),g2=1,6) !fill in each column of data to plot
	flush(30)
testprint=.False.

end do !phi
end do !B
end do !muSC
end do !SCmass
end do !mu
close(30)

open(unit=31,file='infile_JJ_2',status='old') !overwrite the number of steps in each loop (to be used in plots)
do i=1,9 !newlines are counted in
	read(31,*)
enddo
write(tempstr,'(I3)') Anstepphi
write(31,'(a,a,a,a)') 'Anstepphi',achar(9),'=',tempstr
write(tempstr,'(I3)') AnstepB
write(31,'(a,a,a,a,a)') 'AnstepB',achar(9),achar(9),'=',tempstr
write(tempstr,'(I3)') AnstepmuSC
write(31,'(a,a,a,a)') 'AnstepmuSC',achar(9),'=',tempstr
write(tempstr,'(I3)') Anstepmu
write(31,'(a,a,a,a)') 'Anstepmu',achar(9),'=',tempstr
write(tempstr,'(I3)') AnstepSCm
write(31,'(a,a,a,a)') 'AnstepSCm',achar(9),'=',tempstr
write(31,'(a)') '/'
close(31)

call dtime(tarray, result)
tottime=tottime+result
print*,'#total number of loops for Phi, B, muSC, mu, SCmass!',Anstepphi,AnstepB,AnstepmuSC,Anstepmu,AnstepSCm
print*,'#total time(min) = ', tottime/60.0

end program
