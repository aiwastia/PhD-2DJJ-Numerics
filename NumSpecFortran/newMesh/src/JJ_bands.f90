program PD
! calculate spectrum using Friedel sum rule (see Eq. (4) in PRL 107, 196804 (2011))
use io
use selectedsubr
use scans

implicit none

complex*16          :: dS(8,8), S1(8,8), S2(8,8), S3(8,8), Sbar(8,8), unitmat(8,8), Sb(8,8)
complex*16          :: detR, detRold
real*8              :: Deltar,Deltai,Deltar1,Deltai1
real*8              :: gammatot, bmagx, kx, mubuff, mubuff1, kxmax, kxmin, phi, normres, eps,&
			            Ekx, Ekx1, Ekx2, kFmax, Bigdkx, dkx, Vslope
real                ::tarray(2), result, tottime
integer             :: n, i, j, jj, keps, kroot, nkx, Opnbr, Opnbr1, nl

real*8				:: gap, newlevel
real*8				:: vardata(5)
real*8,allocatable,dimension(:) :: fullBZ,newlevelLIST
real*8,allocatable,dimension(:,:) :: fulllevels
character(len=3) 	:: tempstr

character(len=4),dimension(5) :: vartitle
integer :: Ai
real*8,allocatable,dimension(:) :: Aphiarray,ABarray,AmuSCarray,Amuarray

integer :: g,gg,ggg,gggg, g2,lgth
real*8 :: Aloopparam(5)
character(len=100) :: Alecture,tempname
integer :: Apointer




call get_infile_newparameters

open(unit=30,file=dataname,form='formatted',status='new')

vartitle=(/ '#phi','B   ','muSC','mu  ','gap ' /)
write(30,"(5(a,a))") (vartitle(i),achar(9),i=1,5)

call makelist(trim(Aphilist),Aphiarray,'phi')
Anstepphi=size(Aphiarray)
call makelist(trim(ABlist),ABarray,'B')
AnstepB=size(ABarray)
call makelist(trim(AmuSClist),AmuSCarray,'muSC')
AnstepmuSC=size(AmuSCarray)
call makelist(trim(Amulist),Amuarray,'mu')
Anstepmu=size(Amuarray)

tottime=0

!! LOOPS TO COMPUTE THE GAP USING THE SCATTERING CODE

do g=1,Anstepphi
	Aloopparam(1)=Aphiarray(g)
do gg=1,AnstepB
	Aloopparam(2)=ABarray(gg)
do ggg=1,AnstepmuSC
	Aloopparam(3)=AmuSCarray(ggg)
do gggg=1,Anstepmu
	Aloopparam(4)=Amuarray(gggg)
	
	!! MODIFYING THE INPUT FILE
	Apointer=0
	open(unit=22,file='infile_JJ_bands_1',form='formatted',status='old')
	do while (Apointer.eq.0)
		read(22,'(a4)') Alecture
		if (Alecture.eq.'mass') then
			Apointer=1
			write(22,*) 'Deltaphase2',achar(9),'=',Aloopparam(1)
			write(22,*) 'bmag',achar(9),achar(9),'=',Aloopparam(2)
			write(22,*) 'mu0',achar(9),achar(9),'=',Aloopparam(3)
			write(22,*) 'mu01',achar(9),achar(9),'=',Aloopparam(4)
			write(22,*) '/'
		end if
	end do
	close(22)

!!!!!!!!!!!!!!! COMPUTE THE GAP WITH THE NEW INPUT => will be stored in Aloopparam(5) for printing by line

call get_infile_parameter

! change name in infile_1 to name in infile_2
!lgth=len_trim(dataname)
!tempname=trim(dataname(1:lgth-4))
!outfname=trim(tempname) // '.dat'

write(outfname,'(a,I0,a,I0,a,I0,a,I0,a)') 'test',g,'_',gg,'_',ggg,'_',gggg,'.dat'

open(unit=9,file=outfname,status='unknown')

write(9,*)
write(9,*) "# by wire-concat.f90"
write(9,*) "#",outfname
write(9,*) "#mu=", mu0
write(9,*) "#|Delta|(1,2)=", deltamag,deltamag1
write(9,*) "#phase difference=", Deltaphase2
write(9,*) "#bmag=", bmag
write(9,*) "#alpha=", alpha
write(9,*) "#Lmax(1,2,3)=", Lmax,Lmax1,Lmax
write(9,*) "#epsmin=", epsmin
write(9,*) "#epsmax=", epsmax
write(9,*) "#Ndata=", Ndata
write(9,*) "#potshift=", potshift
write(9,*) "#pot=", pot
write(9,*) "#dL=", dL
write(9,*) "#Lj=", Lj
write(9,*) "#SCmass=", SCmass
write(9,*) "#mass=", mass
write(9,*) "#nkxmax=", nkxmax


!!!!!!!!!! PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Deltaphase=0d0
Deltaphase1=0d0
btheta=0d0
btheta1=0d0
bmag1=bmag ! same magnetic field in lead and junction, but bmag=>bmag*0 later
alpha1=alpha !same SOC in the lead and in the junction

gammatot = 0 !disorder amplitude
varphase='scLR' !specify the system (SC on left and right)

! SC phase difference in units of 2*pi
phi = Deltaphase2

 Deltar = Deltamag * cos(Deltaphase*pi)
 Deltai = Deltamag * sin(Deltaphase*pi)
 Deltar1 = Deltamag1 * cos(Deltaphase1*pi)
 Deltai1 = Deltamag1 * sin(Deltaphase1*pi)

 kFmax=mass*alpha+sqrt(2d0*mass*(mu01+bmag1)+mass**2*alpha**2) !keep junction's mass
 kxmax=1.1*kFmax
 kxmin=0. !1.0*kFmax-2*kFmax/(kFmax*Lmax1)**2
 Bigdkx=0.1!2*Lmax1*min((kxmax-kxmin)/dble(nkxmax),kFmax/(kFmax*Lmax1)**2)
print*,'Bigdkx=',Bigdkx
 dkx=min(kFmax/(kFmax*Lmax1)**2, Bigdkx/2.)

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


write(9,*) "#dL=", dL0,dL1,dL0
write(9,*) "#Nseg(1,2,3)=", Nseg,Nseg1,Nseg


!!!!!!!!!!!!!!!! REGULAR SCANNING along kx !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(9,'(A)') '# kx, energy'
!allocate(newlevelLIST(5))
!call horizontal_scanning(Ndata,epsmin,epsmax,dL0,dL1,&
!			gammatot,mu0,mu01,Deltamag,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
!			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
!			Opnbr,Opnbr1,pot,Lj,phi,varphase,kroot,newlevelLIST,nl,gap,&
!			kxmin,kxmax,nkxmax)
!deallocate(newlevelLIST)


!!!!!!!!!!!!!!!! BACKWARD SCANNING with first slope !!!!!!!!!!!!!!!!!!!!!
!write(9,'(A)') '# kx, energy'
!allocate(newlevelLIST(5))
!call backward_scanning_slope(Ndata,epsmin,epsmax,dL0,dL1,&
!			gammatot,mu0,mu01,Deltamag,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
!			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
!			Opnbr,Opnbr1,pot,Lj,phi,varphase,kroot,newlevelLIST,nl,gap,&
!			kxmin,kxmax,dkx,nkxmax)
!deallocate(newlevelLIST)


!!!!!!!!!!!!!!!!dynamical discretization with flat !!!!!!!!!!!!!!!!!!!!!!
!write(9,'(A)') '# kx, energy'
!allocate(fullBZ(flatness*nkxmax))
!allocate(fulllevels(flatness*nkxmax,5))
!
!call dyn_scan_flat(Ndata,epsmin,epsmax,dL0,dL1,&
!			gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
!			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
!			Opnbr,Opnbr1,pot,Lj,phi,varphase,gap,&
!			dkx,Bigdkx,fulllevels,fullBZ,kxmax,kxmin,flatness,nkxmax)
!deallocate(fullBZ)
!deallocate(fulllevels)


!!!!!!!!!!!!!!!!dynamical discretization with curvature (_curv) or only min/max focus (_min) !!!!!!!!!!!!!!!!!
write(9,'(A)') '# kx, energy'
allocate(fullBZ(1000000))
allocate(fulllevels(1000000,5))

call dyn_scan_min(Ndata,epsmin,epsmax,dL0,dL1,&
            gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
            bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
            Opnbr,Opnbr1,pot,Lj,phi,varphase,gap,&
            Bigdkx,fulllevels,fullBZ,kxmax,kxmin,nkxmax)
deallocate(fullBZ)
deallocate(fulllevels)


!!!!!!!!!!!!!!!! VERTICAL SCANNING fix grid !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(9,'(A)') '# energy, kx'
!allocate(newlevelLIST(50))
!call vertical_scanning(nkxmax,kxmin,kxmax,dL0,dL1,&
!			gammatot,mu0,mu01,Deltamag,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
!			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
!			Opnbr,Opnbr1,pot,Lj,phi,varphase,kroot,newlevelLIST,nl,gap,&
!			epsmin,epsmax,Ndata)
!deallocate(newlevelLIST)


!!!!!!!!!!!!!!!!!!! VERTICAL SECANTE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(9,'(A)') '# no data about kx and E, gap only with secante method'
!allocate(newlevelLIST(20))
!call vertical_secante(nkxmax,kxmin,kxmax,dL0,dL1,&
!			gammatot,mu0,mu01,Deltamag,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
!			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
!			Opnbr,Opnbr1,pot,Lj,phi,varphase,kroot,newlevelLIST,nl,gap,&
!			epsmin,epsmax,Ndata)
!deallocate(newlevelLIST)


call flush(9)

write(9,'(1X,A,F15.13)') "#gap =", gap

!sum up the running time of each loop
call dtime(tarray, result)
tottime=tottime+result
write(9,*) "#time(min) = ", result/60.0
close(9)


	Aloopparam(5)=gap
	write(30,'(4(F6.3,a),F15.12,a)') (Aloopparam(g2),achar(9),g2=1,5) !fill in each column of data to plot

end do !mu
end do !muSC
end do !B
end do !phi

close(30)
open(unit=31,file='infile_JJ_2',status='old') !overwrite the number of steps in each loop (to be used in plots)
do i=1,8 !newlines are counted in
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
write(31,'(a)') '/'
close(31)

print*,'#total number of loops for Phi, B, muSC, mu!',Anstepphi,AnstepB,AnstepmuSC,Anstepmu
write(*,*) "#total time(min) = ", tottime/60.0

end program
