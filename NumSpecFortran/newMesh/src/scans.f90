module scans
use selectedsubr

implicit none


contains ! the different methods to create the spectrum and compute the gap depending on the required precision or time


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! COMPUTE ENERGY !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	part of the initial program turned into subroutine for several calls
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
			gammatot,mubuff,mubuff1,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)

integer, intent(in)	:: Ndata,Opnbr,Opnbr1
real*8,intent(in) 	:: epsmin,epsmax,dL0,dL1,&
			Deltar,Deltai,Deltar1,Deltai1,alpha,alpha1,mubuff,mubuff1,SCmass,mass,&
			gammatot,potshift,pot,Lj,bmag,bmag1,btheta,btheta1,phi,kx
character(len=200),intent(in) :: varphase
real*8,intent(out)	:: newlevelLIST(5)
integer, intent(out) :: kroot
real*8,intent(inout) :: gap
integer,intent(inout) :: nl
complex*16,intent(inout) :: detR

integer :: keps,n,i,j
real*8 :: eps,newlevel,bmagx,bmagx1,mu0,mu01
complex*16,dimension(8,8) :: dS,S1,S2,S3,Sbar,unitmat,Sb,massScattL1J,massScattJL2
complex*16 :: detRold



nl=0
kroot=0
newlevelLIST=epsmax

bmagx =-alpha*kx
bmagx1=-alpha1*kx
mu0 = mubuff-kx**2/2d0/SCmass-SCmass*alpha**2/2d0
mu01 = mubuff1-kx**2/2d0/mass-mass*alpha1**2/2d0
!if (mu0+potshift<0) write(*,*) 'mu0+potshift negative: ', mu0, '+',potshift,'<0'
!if (mu01+potshift<0) write(*,*) 'mu0+potshift negative: ', mu01, '+',potshift,'<0'

do keps=0,Ndata !energy loop

!add DGAF for gap
   if (keps.gt.(Ndata/10)) then
	gap=min(gap,epsmax/10.)
	exit
   endif
!!
   eps=epsmin + (epsmax - epsmin)/dble(Ndata) * dble(keps) !scan in energy
   
   ! initialize matrix as unit matrix
   unitmat=dcmplx(0,0)
   forall (i=1:8) unitmat(i,i)=dcmplx(1,0)
   S1=unitmat
   S2=unitmat
   S3=unitmat
   Sbar=0
   do i=1,4
         Sbar(i,i+4) = dcmplx(-1.d0,0.d0)
         Sbar(i+4,i) = dcmplx(-1.d0,0.d0)
   enddo


   !lead 1
   call deltaSwire_1ch(dS,1,eps,dL0,gammatot,mu0,Deltar,Deltai,alpha,bmag,bmagx,btheta,potshift*SCmass/mass,SCmass)
   call concat(dS,S1,S1,4)
   do n=1,Opnbr
      call concat(S1,S1,S1,4)
   enddo

   !junction
   call deltaSwire_1ch(dS,1,eps,dL1,gammatot,mu01,Deltar1,Deltai1,alpha1,bmag1,bmagx1,btheta1,potshift,mass)
   call concat(dS,S2,S2,4)
   do n=1,Opnbr1
      call concat(S2,S2,S2,4)
   enddo

   !lead 2
   S3 = S1

   !add potential barrier at both interfaces to include some normal reflection
   !call Sbarrier(Sb,mu0, eps, pot, Lj, potshift, mass, SCmass)
   !call concat(Sb,S2,S2,4)
   !call concat(S2,Sb,S2,4)

   !create the scattering between regions of different masses (step)
   !call massStep(massScattL1J,SCmass,mass)
   !call concat(S1,massScattL1J,S1,4)
   !call massStep(massScattJL2,mass,SCmass)
   !call concat(massScattJL2,S3,S3,4)

   !compute determinant
   detRold = detR
   detR = detringconcat(S1,S2,S3,Sbar,phi/2.,varphase) 

   ! detR = det(1-S), this changes sign when there is a state
   ! the change in sign for both real and imaginary parts of the determinant is a sufficient condition for non-degenerate states
   if ((keps.gt.0).and.(dreal(detR)*dreal(detRold).lt.0d0).and.(dimag(detR)*dimag(detRold).lt.0d0) ) then 
     kroot=1
     newlevel=eps - (epsmax - epsmin)/dble(Ndata)/2d0 !backward extrapolation of a half step for better accuracy
     nl=nl+1

	!crucial condition because of memory overwriting without segmentation fault !!!!
     if (nl.le.size(newlevelLIST)) then
		newlevelLIST(nl)=newlevel
	 else
		nl=size(newlevelLIST)
	 endif
     !defines the gap as the lowest energy (all bands included)
     if (newlevel.lt.gap) gap=newlevel

   endif
enddo


end subroutine computeEkx




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! HORIZONTAL SCANNING of BZ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine horizontal_scanning(Ndata,epsmin,epsmax,dL0,dL1,&
			gammatot,mu0,mu01,Deltamag,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kroot,newlevelLIST,nl,gap,&
			kxminbuff,kxmaxbuff,nkxmax)

integer, intent(in)	:: Ndata,Opnbr,Opnbr1,nkxmax
real*8,intent(in) 	:: epsmin,epsmax,dL0,dL1,&
			Deltamag,Deltar,Deltai,Deltar1,Deltai1,alpha,alpha1,mu0,mu01,SCmass,mass,&
			gammatot,potshift,pot,Lj,bmag,bmag1,btheta,btheta1,phi,&
			kxminbuff,kxmaxbuff
character(len=200),intent(in) :: varphase
real*8,intent(out)	:: newlevelLIST(5)
integer, intent(out) :: kroot
real*8,intent(inout) :: gap
integer,intent(inout) :: nl

integer :: nkx,jj
real*8 :: kx,kxmin,kxmax
complex*16 :: detR


kxmin=kxminbuff
kxmax=kxmaxbuff

gap=epsmax

do nkx=0,nkxmax !scan of right-most part of the spectrum
  kx = kxmin + dble(nkx)/dble(nkxmax)*(kxmax-kxmin)
  call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
			gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,&
			detR)

  if (kroot.gt.0) then !start a new line
    write(9,'(F20.13,1x,$)') kx !! $ is used to keep writing on the same line
    do jj=1,nl
	write(9,'(F20.13,1x,$)') newlevelLIST(jj)
    enddo
    write(9,*) !! used to finish the current line
  endif
  call flush(9)
enddo ! end of loop over transverse momenta kx
print*,'gap',gap

end subroutine horizontal_scanning




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! BACKWARD SCANNING of BZ with first slope !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine backward_scanning_slope(Ndata,epsmin,epsmax,dL0,dL1,&
			gammatot,mu0,mu01,Deltamag,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kroot,newlevelLIST,nl,gap,&
			kxminbuff,kxmaxbuff,dkx,nkxmax)

integer, intent(in)	:: Ndata,Opnbr,Opnbr1,nkxmax
real*8,intent(in) 	:: epsmin,epsmax,dL0,dL1,&
			Deltamag,Deltar,Deltai,Deltar1,Deltai1,alpha,alpha1,mu0,mu01,SCmass,mass,&
			gammatot,potshift,pot,Lj,bmag,bmag1,btheta,btheta1,phi,&
			kxminbuff,kxmaxbuff,dkx
character(len=200),intent(in) :: varphase
real*8,intent(out)	:: newlevelLIST(5)
integer, intent(out) :: kroot
real*8,intent(inout) :: gap
integer,intent(inout) :: nl

integer :: keps,n,i,j,jj,nkx
real*8 :: kx,kxmin,kxmax,eps,newlevel,Ekx1,Ekx2,Vslope
complex*16 :: detR


kxmin=kxminbuff
kxmax=kxmaxbuff

!first scan to find the highest subgap energy level from the right
Ekx1=2*Deltamag
kx=kxmax
do while (Ekx1.ge.Deltamag)
  kx = kx - dkx/10.
  call computeEkx(Ndata/10,epsmin,epsmax,dL0,dL1,&
			gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,&
			detR)
  Ekx1=minval(newlevelLIST(1:nl))
enddo

kxmax=kx

!compute the slope on the right of the spectrum
  kx = kxmax-dkx/100.
  call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
			gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,&
			detR)
  Ekx2=minval(newlevelLIST(1:nl))
  Vslope=(Ekx1-Ekx2)/(dkx/100.)
  kxmin=kxmax-Ekx1/Vslope

 gap=epsmax

!uses the new mesh around the gap
do nkx=0,nkxmax !scan of right-most part of the spectrum
  kx = kxmin + dble(nkx)/dble(nkxmax)*(kxmax-kxmin)
  call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
			gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,&
			detR)

  if (kroot.gt.0) then !start a new line
    write(9,'(F20.13,1x,$)') kx !! $ is used to keep writing on the same line
    do jj=1,nl
	write(9,'(F20.13,1x,$)') newlevelLIST(jj)
    enddo
    write(9,*) !! used to finish the current line
  endif
  call flush(9)
enddo ! end of loop over transverse momenta kx


end subroutine backward_scanning_slope




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! DYNAMICAL discretization with FLAT !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! flatness is an integer giving the number of levels in a plateau considered as a flat region
!!!!!!!! if flatness=1, all the BZ is discretized with Bigdkx
!!!!!!!! if flatness=nkxmax, " with dkx

subroutine dyn_scan_flat(Ndata,epsmin,epsmax,dL0,dL1,&
			gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,gap,&
			dkx,Bigdkx,fulllevels,fullBZ,kxmax,kxmin,flatness,nkxmax)

integer, intent(in)	:: Ndata,Opnbr,Opnbr1,flatness,nkxmax
real*8,intent(in) 	:: epsmin,epsmax,dL0,dL1,&
			Deltar,Deltai,Deltar1,Deltai1,alpha,alpha1,mu0,mu01,SCmass,mass,&
			gammatot,potshift,pot,Lj,bmag,bmag1,btheta,btheta1,phi,&
			dkx,Bigdkx,kxmax,kxmin
character(len=200),intent(in) :: varphase
real*8,intent(inout) :: gap
real*8,intent(out) :: fulllevels(flatness*nkxmax,5)
real*8,intent(out) :: fullBZ(flatness*nkxmax)

integer :: ff,fff,dd,ddd,nl,kroot,i,jj
complex*16 :: detR
real*8,dimension(flatness) :: Eplateau
real*8,dimension(5) :: newlevelLIST
real*8 :: kx
logical :: test


fullBZ=kxmax
fulllevels=epsmax
kroot=0

kx=kxmin
fullBZ(1)=kx
fff=1
call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
			gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)
Eplateau(1)=minval(newlevelLIST(1:nl))
do ddd=1,5
	fulllevels(1,ddd)=newlevelLIST(ddd)
enddo
ff=1

if (flatness.ge.2) then
	do while (kx.le.(kxmax)) !(Eplateau(1).le.epsmax)
		do dd=1,flatness-1
			kx=kx+dkx
			fff=fff+1
			fullBZ(fff)=kx
			call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
					gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
					bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
					Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)
			if (kroot.gt.0) then !start a new line
    			write(9,'(F20.13,1x,$)') kx !! $ is used to keep writing on the same line
			    do jj=1,nl
					write(9,'(F20.13,1x,$)') newlevelLIST(jj)
			    enddo
			    write(9,*) !! used to finish the current line
			endif
			Eplateau(dd+1)=minval(newlevelLIST(1:nl))
			ff=ff+1
			do ddd=1,5
				fulllevels(ff,ddd)=newlevelLIST(ddd)
			enddo
		enddo

		test=.true.
		do i=1,flatness-1
			test=Eplateau(i).eq.Eplateau(i+1)
			if (.not.test) exit
		enddo
		if(test) then
			kx=kx+Bigdkx
			fff=fff+1
			fullBZ(fff)=kx
			call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
					gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
					bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
					Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)
			if (kroot.gt.0) then !start a new line
    			write(9,'(F20.13,1x,$)') kx !! $ is used to keep writing on the same line
			    do jj=1,nl
					write(9,'(F20.13,1x,$)') newlevelLIST(jj)
			    enddo
			    write(9,*) !! used to finish the current line
			endif
			Eplateau(1)=minval(newlevelLIST(1:nl))
			ff=ff+1
			do ddd=1,5
				fulllevels(ff,ddd)=newlevelLIST(ddd)
			enddo
		else
			Eplateau(1)=Eplateau(flatness)
		endif
	call flush(9)
	enddo
else if (flatness.eq.1) then
	do while (kx.le.(kxmax)) !(Eplateau(1).le.epsmax)
		kx=kx+Bigdkx
		fff=fff+1
		fullBZ(fff)=kx
		call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
				gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
				bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
				Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)
		if (kroot.gt.0) then !start a new line
			write(9,'(F20.13,1x,$)') kx !! $ is used to keep writing on the same line
		    do jj=1,nl
				write(9,'(F20.13,1x,$)') newlevelLIST(jj)
		    enddo
		    write(9,*) !! used to finish the current line
		endif
		Eplateau(1)=minval(newlevelLIST(1:nl))
		ff=ff+1
		do ddd=1,5
			fulllevels(ff,ddd)=newlevelLIST(ddd)
		enddo
	call flush(9)
	enddo
else
	print*, "Enter a plateau length, flatness > 0."
endif

gap=minval(fulllevels)

end subroutine dyn_scan_flat



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! DYNAMICAL discretization with CURV !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!second derivative!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dyn_scan_curv(Ndata,epsmin,epsmax,dL0,dL1,&
        gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
        bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
        Opnbr,Opnbr1,pot,Lj,phi,varphase,gap,&
        Bigdkx,fulllevels,fullBZ,kxmax,kxmin,nkxmax)

integer, intent(in)    :: Ndata,Opnbr,Opnbr1,nkxmax
real*8,intent(in)     :: epsmin,epsmax,dL0,dL1,&
                Deltar,Deltai,Deltar1,Deltai1,alpha,alpha1,mu0,mu01,SCmass,mass,&
                gammatot,potshift,pot,Lj,bmag,bmag1,btheta,btheta1,phi,&
                Bigdkx,kxmax,kxmin
character(len=200),intent(in) :: varphase
real*8,intent(inout) :: gap
real*8,intent(out) :: fulllevels(1000000,5)
real*8,intent(out) :: fullBZ(1000000)

integer :: ff,fff,dd,ddd,nl,nlprev,kroot,i,jj,pp
complex*16 :: detR
real*8,dimension(3) :: Eplateau
real*8,dimension(5) :: newlevelLIST,prevlevel
real*8 :: dkx,kx,kx3,curvmin,curvmax, curvature,Ept,prevcurv
logical :: startwhite,endwhite,testhole,checkplat


fullBZ=kxmax
fulllevels=epsmax
kroot=0
prevlevel=epsmax
nlprev=1
prevcurv=0.

!defines the level of curvature to be resolved
curvmin=-0.2!-1d10  !don't care about resolving local maxima
curvmax=0.2

startwhite=.false.
endwhite=.true.
testhole=.false.

kx=kxmin
fff=0
ff=0

do while (kx.le.(kxmax))

	if(endwhite) then !initialize Eplateau after an empty part of spectrum
		!! INITIAL 1
		fff=fff+1
		fullBZ(fff)=kx
		call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
				gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
				bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
				Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)
		if (kroot.gt.0) then
			write(9,'(F20.13,1x,$)') kx
			do jj=1,nl
				write(9,'(F20.13,1x,$)') newlevelLIST(jj)
			enddo
			write(9,*)
		endif
		Eplateau(1)=minval(newlevelLIST(1:nl))
		ff=ff+1
		do ddd=1,5
			fulllevels(ff,ddd)=newlevelLIST(ddd)
		enddo

		!! INITIAL 2
		kx=kx+Bigdkx
		fff=fff+1
		fullBZ(fff)=kx
		call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
				gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
				bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
				Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)
		if (kroot.gt.0) then
			write(9,'(F20.13,1x,$)') kx
			do jj=1,nl
				write(9,'(F20.13,1x,$)') newlevelLIST(jj)
			enddo
			write(9,*)
		endif
		Eplateau(2)=minval(newlevelLIST(1:nl))
		ff=ff+1
		do ddd=1,5
			fulllevels(ff,ddd)=newlevelLIST(ddd)
		enddo

		startwhite=.false.
		endwhite=.false.

	endif

	!! LOOP over 3
	kx=kx+Bigdkx
	kx3=kx
	call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
		gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
		bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
		Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)
	Eplateau(3)=minval(newlevelLIST(1:nl))

	!compute local curvature (discrete second derivative)
	!print*,Eplateau
	curvature=(Eplateau(3)-2*Eplateau(2)+Eplateau(1))/Bigdkx**2
	print*, curvature

	!when needed, see if lowest level is higher than ground state (missed point)
	if((prevcurv.le.curvmax).and.(prevcurv.ge.curvmin).and.((curvature.gt.curvmax).or.(curvature.lt.curvmin)).and.(nlprev.ge.2)) then
		do jj=2,nlprev
			testhole=testhole.or.(abs(Eplateau(3)-prevlevel(jj)).le.10**(-12))
		enddo
	endif

	if(testhole) then !insert the missed point by prolongation
		do i=1,nl
			newlevelLIST(nl+2-i)=newlevelLIST(nl+1-i)
		enddo
		newlevelLIST(1)=minval(prevlevel(1:nlprev))
		nl=nl+1
		Eplateau(3)=newlevelLIST(1)
		curvature=(Eplateau(3)-2*Eplateau(2)+Eplateau(1))/Bigdkx**2
		testhole=.false.
	endif


	!keep or start over  the previous segment E2-E3

	if((curvature.le.curvmax).and.(curvature.ge.curvmin)) then
		if(startwhite) then
			endwhite=.true.
		else
			nlprev=nl
			prevlevel=newlevelLIST
			!!!!! print if previous energy accepted
			fff=fff+1
			fullBZ(fff)=kx
			if (kroot.gt.0) then
				write(9,'(F20.13,1x,$)') kx
				do jj=1,nl
		            write(9,'(F20.13,1x,$)') newlevelLIST(jj)
		        enddo
		        write(9,*)
		    endif
		    ff=ff+1
		    do ddd=1,5
		        fulllevels(ff,ddd)=newlevelLIST(ddd)
		    enddo

		    Eplateau(1)=Eplateau(2)
		    Eplateau(2)=Eplateau(3)
		endif

	else if (.not.(curvature.eq.(curvature-1))) then !test for .not.Infinite curvature
		if(startwhite) then
			endwhite=.true.
		else
			!!!!! re-descretize the previous segment
			kx=fullBZ(fff)
			dkx=Bigdkx/(1000.*log(1+abs(curvature)/100.))
			do while (kx.le.kx3)
				kx=kx+dkx
				fff=fff+1
				fullBZ(fff)=kx
				call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
					gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
					bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
					Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)
				Ept=minval(newlevelLIST(1:nl))
				!again check the holes in the subdivision of the spectrum
				if(nlprev.ge.2) then
					do jj=2,nlprev
						testhole=testhole.or.(abs(Ept-prevlevel(jj)).le.10**(-12))
					enddo
				endif
				if(testhole) then
					do i=1,nl
						newlevelLIST(nl+2-i)=newlevelLIST(nl+1-i)
					enddo
					newlevelLIST(1)=minval(prevlevel(1:nlprev))
					nl=nl+1
					Eplateau(3)=newlevelLIST(1)
					testhole=.false.
				endif
				nlprev=nl
				prevlevel=newlevelLIST

				if (kroot.gt.0) then
					write(9,'(F20.13,1x,$)') kx
					do jj=1,nl
						write(9,'(F20.13,1x,$)') newlevelLIST(jj)
					enddo
					write(9,*)
				endif
				ff=ff+1
				do ddd=1,5
					fulllevels(ff,ddd)=newlevelLIST(ddd)
				enddo

				checkplat=.true.
				if (ff.gt.20) then
				do pp=ff-20,ff-1
					checkplat=checkplat.and.(fulllevels(pp,1).eq.fulllevels(pp+1,1))
				enddo
				if (checkplat) then
					dkx=5*dkx
				else
					dkx=Bigdkx/(1000.*log(1+abs(curvature)/100.))
				endif
				endif

			enddo
			Eplateau(1)=Eplateau(2)
			Eplateau(2)=Eplateau(3)
		endif

	else
		startwhite=.true.
	endif

	prevcurv=curvature

	!if (ff>1108) stop
	call flush(9)
	!print*,ff

enddo


gap=minval(fulllevels)

end subroutine dyn_scan_curv







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! DYNAMICAL discretization with MIN !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!second derivative!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dyn_scan_min(Ndata,epsmin,epsmax,dL0,dL1,&
        gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
        bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
        Opnbr,Opnbr1,pot,Lj,phi,varphase,gap,&
        Bigdkx,fulllevels,fullBZ,kxmax,kxmin,nkxmax)

integer, intent(in)    :: Ndata,Opnbr,Opnbr1,nkxmax
real*8,intent(in)     :: epsmin,epsmax,dL0,dL1,&
                Deltar,Deltai,Deltar1,Deltai1,alpha,alpha1,mu0,mu01,SCmass,mass,&
                gammatot,potshift,pot,Lj,bmag,bmag1,btheta,btheta1,phi,&
                Bigdkx,kxmax,kxmin
character(len=200),intent(in) :: varphase
real*8,intent(inout) :: gap
real*8,intent(out) :: fulllevels(1000000,5)
real*8,intent(out) :: fullBZ(1000000)

integer :: ff,fff,dd,ddd,nl,nlprev,kroot,i,jj,selbot,stspec,igap,platcoord
complex*16 :: detR
real*8,dimension(3) :: Eplateau
real*8,dimension(5) :: newlevelLIST,prevlevel
real*8,dimension(50) :: spec,diffgap
real*8 :: dkx,kx,kx3,kx2,kx1,tempdkx,intergap,intergap1,intergap2
real*8 :: curvmin,curvmax,curvature,Ept,prevcurv,slope1,slope2,prevslope,BBigdkx,slopemin,prevbdkx
logical :: startwhite,endwhite,testhole,checkdown,checkgap,checkslope2,getout


BBigdkx=Bigdkx
fullBZ=kxmax
fulllevels=epsmax
spec=epsmax
diffgap=epsmax
kroot=0
prevlevel=epsmax
nlprev=1
prevcurv=0.

!defines the level of curvature to be resolved
curvmin=-10d13 !-1d10  !don't care about resolving local maxima
curvmax=10.
slopemin=10.*(epsmax-epsmin)/(kxmax-kxmin)

startwhite=.false.
endwhite=.true.
testhole=.false.
getout=.false.

kx=kxmin-Bigdkx
fff=0
ff=0
print*,'kxmax',kxmax
do while ((kx.le.kxmax).and..not.getout)

	do while(endwhite) !initialize Eplateau after an empty part of spectrum
		!! INITIAL 1
		kx=kx+Bigdkx
		fff=fff+1
		fullBZ(fff)=kx
print*,'E1k',kx
		call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
				gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
				bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
				Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)
		if (kroot.gt.0) then !implies nl>1
			write(9,'(F20.13,1x,$)') kx
print*,"E1",kx
			do jj=1,nl
				write(9,'(F20.13,1x,$)') newlevelLIST(jj)
			enddo
			write(9,*)
			Eplateau=minval(newlevelLIST(1:nl))
			ff=ff+1
			do ddd=1,5
				fulllevels(ff,ddd)=newlevelLIST(ddd)
			enddo
			startwhite=.false.
			endwhite=.false.
		endif

		!! INITIAL 2
		kx=kx+BBigdkx
		fff=fff+1
		fullBZ(fff)=kx
print*,'E2k',kx
		call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
				gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
				bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
				Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)
		if (kroot.gt.0) then
			write(9,'(F20.13,1x,$)') kx
print*,"E2",kx
			do jj=1,nl
				write(9,'(F20.13,1x,$)') newlevelLIST(jj)
			enddo
			write(9,*)
			if(endwhite) then
				 Eplateau=minval(newlevelLIST(1:nl))
			else
				Eplateau(2)=minval(newlevelLIST(1:nl))
			endif
			ff=ff+1
			do ddd=1,5
				fulllevels(ff,ddd)=newlevelLIST(ddd)
			enddo
			startwhite=.false.
			endwhite=.false.
		endif
	slope1=(Eplateau(2)-Eplateau(1))/BBigdkx

	if (kx.gt.kxmax) then
print*,"get out"
		getout=.true.
		exit
	endif
	enddo

	!slope1 is copied from slope2
	!print*, 's',slope1
	prevbdkx=BBigdkx

	!if ((abs(slope1).lt.slopemin).or.(abs(prevcurv).lt.0.01)) then
	!	BBigdkx=Bigdkx*min(slopemin/abs(slope1),50.)
	!else if (abs(slope1).eq.(abs(slope1)-1)) then
	!	BBigdkx=10.*Bigdkx
	!else
	!	BBigdkx=Bigdkx
	!endif
	!print*,BBigdkx


	!store the location of kx2
	platcoord=fff

	!! LOOP over 3
	kx=kx+BBigdkx
print*,'k',kx
	kx3=kx
	call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
		gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
		bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
		Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)
	if(nl.ge.1) Eplateau(3)=minval(newlevelLIST(1:nl))

	slope2=(Eplateau(3)-Eplateau(2))/BBigdkx

	!compute local curvature (discrete second derivative)
!print*,Eplateau
	curvature=(slope2-slope1)/(BBigdkx+prevbdkx)*2d0
!print*, 'c',curvature


	!reporting slope2 to avoid plateaux
	if(slope2.eq.0d0) then
print*,"use slope2"

		!store and print kx3,E3
				fff=fff+1 
				fullBZ(fff)=kx
				if (kroot.gt.0) then
					write(9,'(F20.13,1x,$)') kx
print*, "BZ kx3 ",kx
					do jj=1,nl
						write(9,'(F20.13,1x,$)') newlevelLIST(jj)
					enddo
					write(9,*)
				endif
				ff=ff+1
				do ddd=1,5
					fulllevels(ff,ddd)=newlevelLIST(ddd)
				enddo

		!loop to find the end of the plateau + storage
		do while (slope2.eq.0d0)
			kx=kx+BBigdkx
			if (kx.gt.kxmax) then
print*,"get out"
				getout=.true.
				exit
			endif
					fff=fff+1
					fullBZ(fff)=kx
print*,"BZ slope2 ",kx
					call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
						gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
						bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
						Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)
					if(nl.ge.1) Eplateau(3)=minval(newlevelLIST(1:nl))
					if (kroot.gt.0) then
						write(9,'(F20.13,1x,$)') kx
						do jj=1,nl
							write(9,'(F20.13,1x,$)') newlevelLIST(jj)
						enddo
						write(9,*)
					endif
					ff=ff+1
					do ddd=1,5
						fulllevels(ff,ddd)=newlevelLIST(ddd)
					enddo

			slope2=(Eplateau(3)-Eplateau(2))/BBigdkx !only E3 changes in the loop
		enddo

		!check if it's a plateau or a local max/min
		if(slope2*slope1.ge.0d0) then
			checkslope2=.true.
			curvature=(slope2-slope1)/((fff-platcoord)*Bigdkx) !average curvature
			Eplateau(1)=Eplateau(2)
			Eplateau(2)=Eplateau(3)
		else
			!prepare to discretize the min
			kx3=kx
		endif
	endif
!print*,"slope1=",slope1
!print*,"slope2=",slope2

	!When needed, see if lowest level is higher than ground state (missed point)
	if((.not.checkslope2).and.(.not.(curvature.eq.(curvature-1)))&
.and.((slope1*slope2.le.0).and.(slope1.le.0.))) then !ie, if "Discretizing"
!(((prevcurv.le.curvmax).and.(prevcurv.ge.curvmin)).or.(abs(prevcurv).eq.(abs(prevcurv)-1)))&
!((curvature.gt.curvmax).or.(curvature.lt.curvmin))&
!.and.(nlprev.ge.2)&
!.and..not.(slope2.eq.0d0)) then
		!discretize a bit more
print*,"passe"
		kx2=fullBZ(fff)
		tempdkx=(kx3-kx2)/10.
		stspec=0
		checkgap=.false.
		do i=1,10
			kx=kx2+i*tempdkx
			call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
				gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
				bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
				Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)
			if(nl.gt.1) nl=1 !! needed to avoid taking into account a higher gap (łe2> - łe1>)
! and relying on luck to find łe1> when not łg>.
			if(nl.ge.1) spec(stspec+1:stspec+nl)=newlevelLIST(1:nl)
			stspec=stspec+nl
		enddo
		kx=kx3
		!sort all the values
		call quicksort(spec,1,stspec)
		!define the existence of a gap and if the previous Energy3 was above it
		forall(i=1:stspec-1) diffgap(i)=spec(i+1)-spec(i)
		igap=maxloc(diffgap(1:stspec-1),1)
		intergap=diffgap(igap)
		intergap1=maxval(diffgap(1:igap-1))
		intergap2=maxval(diffgap(igap+1:stspec-1))
		if ((intergap1.gt.0).and.(intergap2.gt.0)&
.and.(intergap.gt.(3*intergap1)).and.(intergap.gt.(3*intergap2))) checkgap=.true.
		if ((checkgap).and.(Eplateau(3).gt.spec(igap-1))) then
print*,"Changed E3"
			Eplateau(3)=Infinity
			slope2=Infinity
			curvature=Infinity
		endif
	endif


	!keep or start over the 2 previous segments E1-E2 & E2-E3
	
	if (checkslope2) then
		checkslope2=.false.
	else if (.not.(curvature.eq.(curvature-1))) then !test for .not.Infinite curvature
		if((slope1*slope2.gt.0).or.(slope1.gt.0.)) then
			if(startwhite) then
				endwhite=.true.
			else
				nlprev=nl
				prevlevel=newlevelLIST
				!!!!! write if previous energy accepted
				fff=fff+1
				fullBZ(fff)=kx
				if (kroot.gt.0) then
					write(9,'(F20.13,1x,$)') kx
print*,"normal E",kx
					do jj=1,nl
						write(9,'(F20.13,1x,$)') newlevelLIST(jj)
					enddo
					write(9,*)
				endif
				ff=ff+1
				do ddd=1,5
				    fulllevels(ff,ddd)=newlevelLIST(ddd)
				enddo

				Eplateau(1)=Eplateau(2)
				Eplateau(2)=Eplateau(3)
			endif
		else
print*,"Discretizing"
			if(startwhite) then
				endwhite=.true.
			else
				selbot=0
				checkdown=.false.
				!!!!! re-descretize the 2 previous segments
				kx=fullBZ(platcoord-1)
print*,curvature
				dkx=min(BBigdkx/1.5,BBigdkx/(100.*log(1+abs(curvature)/100.)))  !!!!!!! here /curvmax
!print*,dkx
				do while (kx.lt.kx3)
					kx=kx+dkx
					fff=fff+1
					fullBZ(fff)=kx
					call computeEkx(Ndata,epsmin,epsmax,dL0,dL1,&
						gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
						bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
						Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gap,detR)
					if(nl.ge.1) Ept=minval(newlevelLIST(1:nl))

					!remove possible local max
					if (Ept.lt.minval(prevlevel(1:nlprev))) checkdown=.true.
					if (checkdown.and.(Ept.gt.minval(prevlevel(1:nlprev)))) selbot=selbot+1
				
					!prevslope=(Ept-minval(prevlevel(1:nlprev)))/dkx
					nlprev=nl
					prevlevel=newlevelLIST

					if (kroot.gt.0) then
						write(9,'(F20.13,1x,$)') kx
print*,"disc E",kx
						do jj=1,nl
							write(9,'(F20.13,1x,$)') newlevelLIST(jj)
						enddo
						write(9,*)
					endif
					ff=ff+1
					do ddd=1,5
						fulllevels(ff,ddd)=newlevelLIST(ddd)
					enddo
				
					if (selbot.eq.5) then !! avoid going upwards after discretizing a local min
						kx=kx3
						exit
					endif
				enddo
				Eplateau(1)=Eplateau(2)
				Eplateau(2)=Eplateau(3)
print*,"end"
			endif
		endif
	else
		startwhite=.true.
	endif

	prevcurv=curvature
	!prevslope=slope1
	slope1=slope2

	call flush(9)

enddo


gap=minval(fulllevels)

end subroutine dyn_scan_min






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! ENERGY SCANNING kx(E) at fixed spacing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine computeKxe(Ndata,epsmin,epsmax,& !!meaning nkxmax kxmin kxmax
			dL0,dL1,&
			gammatot,mubuff,mubuff1,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gaptest,detR) !!kx meaning energy level En

integer, intent(in)	:: Ndata,Opnbr,Opnbr1
real*8,intent(in) 	:: epsmin,epsmax,dL0,dL1,&
			Deltar,Deltai,Deltar1,Deltai1,alpha,alpha1,mubuff,mubuff1,SCmass,mass,&
			gammatot,potshift,pot,Lj,bmag,bmag1,btheta,btheta1,phi,kx
character(len=200),intent(in) :: varphase
real*8,intent(out)	:: newlevelLIST(20)
integer, intent(out) :: kroot
integer,intent(inout) :: gaptest
integer,intent(inout) :: nl
complex*16,intent(inout) :: detR

integer :: keps,n,i,j
real*8 :: eps,newlevel,bmagx,mu0,mu01
complex*16 :: dS(8,8),S1(8,8),S2(8,8),S3(8,8),Sbar(8,8),unitmat(8,8),Sb(8,8),massScattL1J(8,8),massScattJL2(8,8)
complex*16 :: detRold

!! comment : (*) means kx <-> eps

nl=0
kroot=0
newlevelLIST=epsmax

do keps=0,Ndata !energy loop

   eps=epsmin + (epsmax - epsmin)/dble(Ndata) * dble(keps) !scan in energy
	bmagx =-alpha*eps !(*)
	mu0 = mubuff-eps**2/2d0/SCmass !(*)
	mu01 = mubuff1-eps**2/2d0/mass !(*)
	if (mu0+potshift<0) write(*,*) 'mu0+potshift negative: ', mu0, '+',potshift,'<0'
	if (mu01+potshift<0) write(*,*) 'mu0+potshift negative: ', mu01, '+',potshift,'<0'
   
   ! initialize matrix as unit matrix
   unitmat=dcmplx(0,0)
   forall (i=1:8) unitmat(i,i)=dcmplx(1,0)
   S1=unitmat
   S2=unitmat
   S3=unitmat
   Sbar=0
   do i=1,4
         Sbar(i,i+4) = dcmplx(-1.d0,0.d0)
         Sbar(i+4,i) = dcmplx(-1.d0,0.d0)
   enddo


   !lead 1
   call deltaSwire_1ch(dS,1,kx,dL0,gammatot,mu0,Deltar,Deltai,alpha,0d0*bmag,bmagx,btheta,potshift,SCmass) !(*)
   call concat(dS,S1,S1,4)
   do n=1,Opnbr
      call concat(S1,S1,S1,4)
   enddo

   !junction
   call deltaSwire_1ch(dS,1,kx,dL1,gammatot,mu01,Deltar1,Deltai1,alpha1,bmag1,bmagx,btheta1,potshift,mass) !(*)
   call concat(dS,S2,S2,4)
   do n=1,Opnbr1
      call concat(S2,S2,S2,4)
   enddo

   !lead 2
   S3 = S1

   !add potential barrier at both interfaces to include some normal reflection
   !call Sbarrier(Sb,mu0, kx, pot, Lj, potshift, mass, SCmass)
   !call concat(Sb,S2,S2,4)
   !call concat(S2,Sb,S2,4)

   !create the scattering between regions of different masses (step)
   call massStep(massScattL1J,SCmass,mass)
   call concat(S1,massScattL1J,S1,4)
   call massStep(massScattJL2,mass,SCmass)
   call concat(massScattJL2,S3,S3,4)

   !compute determinant
   detRold = detR
   detR = detringconcat(S1,S2,S3,Sbar,phi/2.,varphase) 

   ! detR = det(1-S), this changes sign when there is a state
   ! the change in sign for both real and imaginary parts of the determinant is a sufficient condition for non-degenerate states
   if ((keps.gt.0).and.(dreal(detR)*dreal(detRold).lt.0d0).and.(dimag(detR)*dimag(detRold).lt.0d0) ) then 
     kroot=1
     newlevel=eps - (epsmax - epsmin)/dble(Ndata)/2d0 !backward extrapolation of a half step for better accuracy
     nl=nl+1
     !crucial condition because of memory overwriting without segmentation fault !!!!
     if (nl.le.size(newlevelLIST)) then
		newlevelLIST(nl)=newlevel
	 else
		nl=size(newlevelLIST)
	 endif
   endif
enddo

if(kroot.eq.1) gaptest=gaptest+1


end subroutine computeKxe




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! VERTICAL scanning of kx(E) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vertical_scanning(Ndata,epsmin,epsmax,& !! meaning nkxmax,kxmin,kxmax
			dL0,dL1,&
			gammatot,mu0,mu01,Deltamag,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kroot,newlevelLIST,nl,gap,&
			kxminbuff,kxmaxbuff,nkxmax) !! meaning epsmin,epsmax,Ndata

integer, intent(in)	:: Ndata,Opnbr,Opnbr1,nkxmax
real*8,intent(in) 	:: epsmin,epsmax,dL0,dL1,&
			Deltamag,Deltar,Deltai,Deltar1,Deltai1,alpha,alpha1,mu0,mu01,SCmass,mass,&
			gammatot,potshift,pot,Lj,bmag,bmag1,btheta,btheta1,phi,&
			kxminbuff,kxmaxbuff
character(len=200),intent(in) :: varphase
real*8,intent(out)	:: newlevelLIST(50)
integer, intent(out) :: kroot
real*8,intent(out) :: gap
integer,intent(inout) :: nl

integer :: jj,nkx,gaptest
real*8 :: dkx,kx,kxmin,kxmax
complex*16 :: detR


kxmin=kxminbuff
kxmax=kxmaxbuff
dkx=(kxmax-kxmin)/dble(nkxmax)
gaptest=0

do nkx=0,nkxmax
  kx = kxmin + dble(nkx)*dkx
  call computeKxe(Ndata,epsmin,epsmax,dL0,dL1,&
			gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gaptest,detR)
  if (nl.gt.size(newlevelLIST)) print*, "More roots kx(E) than allowed by the length of newlevelLIST."
  if (kroot.gt.0) then !start a new line
    write(9,'(F20.13,1x,$)') kx !! $ is used to keep writing on the same line
    do jj=1,min(nl,size(newlevelLIST))
	write(9,'(F20.13,1x,$)') newlevelLIST(jj)
    enddo
    write(9,*) !! used to finish the current line
  endif

!  if (gaptest.eq.1) then
!	gap=kx
!  	exit
!  endif
  call flush(9)
enddo ! end of loop over ENERGY kx

end subroutine vertical_scanning




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! VERTICAL secante along E !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vertical_secante(Ndata,epsmin,epsmax,& !! meaning nkxmax,kxmin,kxmax
			dL0,dL1,&
			gammatot,mu0,mu01,Deltamag,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kroot,newlevelLIST,nl,gap,&
			kxminbuff,kxmaxbuff,nkxmax) !! meaning epsmin,epsmax,Ndata

integer, intent(in)	:: Ndata,Opnbr,Opnbr1,nkxmax
real*8,intent(in) 	:: epsmin,epsmax,dL0,dL1,&
			Deltamag,Deltar,Deltai,Deltar1,Deltai1,alpha,alpha1,mu0,mu01,SCmass,mass,&
			gammatot,potshift,pot,Lj,bmag,bmag1,btheta,btheta1,phi,&
			kxminbuff,kxmaxbuff
character(len=200),intent(in) :: varphase
real*8,intent(out)	:: newlevelLIST(20)
integer, intent(out) :: kroot
real*8,intent(out) :: gap
integer,intent(inout) :: nl

integer :: jj,nkx,gaptest,sec,secant
real*8 :: kx,kxmin,kxmax
complex*16 :: detR


secant=12

kxmin=kxminbuff
kxmax=kxmaxbuff
gaptest=0

do sec=1,secant
	kx = (kxmax+kxmin)/2.
	call computeKxe(Ndata,epsmin,epsmax,dL0,dL1,&
			gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gaptest,detR)
	
	if(gaptest.eq.1) then
		kxmax=kx
		gaptest=0
	else
		kxmin=kx
		gaptest=0
	endif

enddo

gap=(kxmax+kxmin)/2.

end subroutine vertical_secante




end module
