program PD

! calculate spectrum using Friedel sum rule (see Eq. (4) in PRL 107, 196804 (2011))

use io
use selectedsubr
use scans

implicit none

real*8              :: Deltar,Deltai,Deltar1,Deltai1
real*8              :: gammatot, bmagx, kx, mubuff, mubuff1, kxmax, kxmin, phi, normres, eps,&
			            Ekx, Ekx1, Ekx2, kFmax, dkx, gap
real                ::tarray(2), result, tottime
integer             :: n, i, jj, Opnbr, Opnbr1,nkxmax

real*8,dimension(5)		:: vardata
character(len=3) 	:: tempstr

character(len=4),dimension(6) :: vartitle
integer :: Ai
real*8,allocatable,dimension(:) :: Aphiarray,ABarray,AmuSCarray,Amuarray,ASCmarray

integer :: g,gg,ggg,gggg,ggggg, g2,lgth
real*8 :: Aloopparam(6)
character(len=100) :: Alecture
integer :: Apointer


!!!!!!! start
tottime=0


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


!! LOOPS TO COMPUTE THE GAP USING THE SCATTERING CODE

!if change order here, change in Aloopparam , write(22), vartitle and write(30 '') at the end
do gggg=1,Anstepmu
	Aloopparam(1)=Amuarray(gggg)
do ggggg=1,AnstepSCm
	Aloopparam(2)=ASCmarray(ggggg)
do ggg=1,AnstepmuSC
	Aloopparam(3)=AmuSCarray(ggg)
do gg=1,AnstepB
	Aloopparam(4)=ABarray(gg)
do g=1,Anstepphi
	Aloopparam(5)=Aphiarray(g)
	
	!! MODIFYING THE INPUT FILE
	Apointer=0
	open(unit=22,file='infile_JJ_bands_1',form='formatted',status='old')
	do while (Apointer.eq.0)
		read(22,'(a4)') Alecture
		if (Alecture.eq.'mass') then
			Apointer=1
			write(22,*) 'Deltaphase2',achar(9),'=',Aloopparam(5)
			write(22,*) 'bmag',achar(9),achar(9),'=',Aloopparam(4)
			write(22,*) 'mu0',achar(9),achar(9),'=',Aloopparam(3)
			write(22,*) 'mu01',achar(9),achar(9),'=',Aloopparam(1)
			write(22,*) 'SCmass',achar(9),achar(9),'=',Aloopparam(2)
			write(22,*) '/'
		end if
	end do
	close(22)

!!!!!!!!!!!!!!! COMPUTE THE GAP WITH THE NEW INPUT => will be stored in Aloopparam(5) for printing by line

call get_infile_parameter

!!!!!!!!!! PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Deltaphase=0d0 !common shift of the phases
Deltaphase1=0d0 ! in junction
btheta=0d0
btheta1=0d0
bmag1=bmag ! same magnetic field in lead and junction
bmag=0. !set to 0 in leads
alpha1=alpha !same SOC in the lead and in the junction
!alpha=0. !set to 0 in leads

gammatot = 0 !disorder amplitude
varphase='scLR' !specify the system (SC on left and right)

! SC phase difference in units of 2*pi
phi = Deltaphase2

Deltar = Deltamag * cos(Deltaphase*pi)
Deltai = Deltamag * sin(Deltaphase*pi)
Deltar1 = Deltamag1 * cos(Deltaphase1*pi)
Deltai1 = Deltamag1 * sin(Deltaphase1*pi)

kFmax=mass*alpha1+sqrt(2d0*mass*(mu01+bmag1)+mass**2*alpha1**2) !in junction
kxmax=1.1*kFmax
kxmin=0. !1.0*kFmax-2*kFmax/(kFmax*Lmax1)**2
dkx=kFmax/(kFmax*Lmax1)**2 / dkfrac
nkxmax=int((kxmax-kxmin)/dkx)

!!!!!!!!!!! SPACE DISCRETIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !two leads
 if (Lmax.ge.dL) then
     Opnbr=min(30,int(log(Lmax/dL)/log(2d0))+1)
     Nseg = 2**Opnbr
     dL0 = Lmax / dble(Nseg)
 else
     dL0=dL
     Nseg=0
 endif 

 !junction
 if (Lmax1.gt.dL) then
     Opnbr1=min(30,int(log(Lmax1/dL)/log(2d0))+1)
     Nseg1 = 2**Opnbr1
     dL1 = Lmax1 / dble(Nseg1)
 else
     dL1=dL
     Nseg1=0
 endif


!!!!!!!!!!!!!!!!!!! VERTICAL SECANTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!call vertical_secante(nkxmax,kxmin,kxmax,dL0,dL1,&
!			gammatot,mu0,mu01,Deltamag,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
!			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
!			Opnbr,Opnbr1,pot,Lj,phi,varphase,gap,&
!			epsmin,epsmax,secant)


!!!!!!!!!!!!!!!!!!! VERTICAL DYNAMICAL SECANTE !!!!!!!!!!!!!!!!!!!!!!
call dynamical_secante(nkxmax,kxmin,kxmax,dL0,dL1,&
            gammatot,mu0,mu01,Deltamag,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
            bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
            Opnbr,Opnbr1,pot,Lj,phi,varphase,gap,&
            epsmin,epsmax,secant,5)  !!can ask for the maximal depth, will stop at the last relevant


!sum up the running time of each loop
call dtime(tarray, result)
tottime=tottime+result

	Aloopparam(6)=gap
	write(30,'(F8.4,a,F8.4,a,F8.3,a,F8.4,a,F8.6,a,F15.12,a)') (Aloopparam(g2),achar(9),g2=1,6) !fill in each column of data to plot
    flush(30)

end do !SCmass
end do !mu
end do !muSC
end do !B
end do !phi

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
write(*,*) "#total time(min) = ", tottime/60.0

end program
