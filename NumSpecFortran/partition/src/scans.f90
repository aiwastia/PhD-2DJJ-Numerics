module scans
use selectedsubr

implicit none


contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! ENERGY SCANNING kx(E) at fixed spacing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine computeKxe(Ndata,epsmin,epsmax,& !!meaning nkxmax kxmin kxmax
			dL0,dL1,&
			gammatot,mubuff,mubuff1,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gaptest,&
            detR,active,insite) !!kx meaning energy level En

integer, intent(in)	:: Ndata,Opnbr,Opnbr1,active
real*8,intent(in) 	:: epsmin,epsmax,dL0,dL1,&
			Deltar,Deltai,Deltar1,Deltai1,alpha,alpha1,mubuff,mubuff1,SCmass,mass,&
			gammatot,potshift,pot,Lj,bmag,bmag1,btheta,btheta1,phi,kx
character(len=200),intent(in) :: varphase
real*8,intent(out)	:: newlevelLIST(20)
integer, intent(out) :: insite,nl
integer,intent(out) :: gaptest
complex*16,intent(inout) :: detR

integer :: keps,n,i,kroot
real*8 :: eps,newlevel,bmagx,bmagx1,mu0,mu01
complex*16 :: dS(8,8),S1(8,8),S2(8,8),S3(8,8),Sbar(8,8),unitmat(8,8),Sb(8,8),massScattL1J(8,8),massScattJL2(8,8)
complex*16 :: detRold
logical :: intest

!! comment : (*) means kx <-> eps

nl=0
kroot=0
newlevelLIST=epsmax
intest=.false.
insite=0

gaptest=0

do keps=0,Ndata !energy loop

	!eps=epsmin + (epsmax - epsmin)/dble(Ndata) * dble(keps) !scan in energy
	eps=epsmax - (epsmax - epsmin)/dble(Ndata) * dble(keps) !scan in energy
	bmagx =-alpha*eps !(*)
	bmagx1=-alpha1*eps !(*)
	!ma^2/2 because affects velocity which is the only physical dependency of PD
	mu0 = mubuff-eps**2/2d0/SCmass-SCmass*alpha**2/2d0 !(*)
	mu01 = mubuff1-eps**2/2d0/mass-mass*alpha1**2/2d0 !(*)
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
	call deltaSwire_1ch(dS,1,kx,dL0,gammatot,mu0,Deltar,Deltai,alpha,bmag,bmagx,btheta,potshift*SCmass/mass,SCmass) !(*)
	call concat(dS,S1,S1,4)
	do n=1,Opnbr
		call concat(S1,S1,S1,4)
	enddo

	!junction
	call deltaSwire_1ch(dS,1,kx,dL1,gammatot,mu01,Deltar1,Deltai1,alpha1,bmag1,bmagx1,btheta1,potshift,mass) !(*)
	call concat(dS,S2,S2,4)
	do n=1,Opnbr1
		call concat(S2,S2,S2,4)
	enddo

	!lead 2
	S3 = S1

	!compute determinant
	detRold = detR
	detR = detringconcat(S1,S2,S3,Sbar,phi/2.,varphase)

	if ((active.eq.1).and.(intest)) insite=insite+1

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
		intest=.not.intest
	endif
   
	if(insite.ge.3) then
		!print*,"exit after ", keps, " steps"
		exit
	endif
enddo

if(kroot.eq.1) gaptest=1

end subroutine computeKxe



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! VERTICAL secante along E !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vertical_secante(Ndata,epsmin,epsmax,& !! meaning nkxmax,kxmin,kxmax
			dL0,dL1,&
			gammatot,mu0,mu01,Deltamag,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,gap,&
			kxminbuff,kxmaxbuff,secant) !! meaning epsmin,epsmax,Ndata

integer, intent(in)	:: Ndata,Opnbr,Opnbr1,secant
real*8,intent(in) 	:: epsmin,epsmax,dL0,dL1,&
			Deltamag,Deltar,Deltai,Deltar1,Deltai1,alpha,alpha1,mu0,mu01,SCmass,mass,&
			gammatot,potshift,pot,Lj,bmag,bmag1,btheta,btheta1,phi,&
			kxminbuff,kxmaxbuff
character(len=200),intent(in) :: varphase
real*8,intent(out) :: gap

integer :: jj,nkx,gaptest,sec,insite,kroot,nl
real*8 :: kx,kxmin,kxmax,deps
real*8 :: newlevelLIST(20)
complex*16 :: detR


kxmin=kxminbuff
kxmax=kxmaxbuff
gaptest=0

do sec=1,secant
	kx = (kxmax+kxmin)/2d0
	call computeKxe(Ndata,epsmin,epsmax,dL0,dL1,&
			gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
			bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
			Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gaptest,&
            detR,0,insite)
	
    if(gaptest.eq.1) then
        kxmax=kx
        gaptest=0
    else
        kxmin=kx
        gaptest=0
    endif

enddo

gap=(kxmax+kxmin)/2d0

end subroutine vertical_secante




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! VERTICAL secante & DYNAMICAL kx !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dynamical_secante(Ndatabuff,epsminbuff,epsmaxbuff,& !! meaning nkxmax,kxmin,kxmax
            dL0,dL1,&
            gammatot,mu0,mu01,Deltamag,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
            bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
            Opnbr,Opnbr1,pot,Lj,phi,varphase,gap,&
            kxminbuff,kxmaxbuff,secant,depth) !! meaning epsmin,epsmax,Ndata

integer, intent(in)    :: Ndatabuff,Opnbr,Opnbr1,secant,depth
real*8,intent(in)     :: epsminbuff,epsmaxbuff,dL0,dL1,&
            Deltamag,Deltar,Deltai,Deltar1,Deltai1,alpha,alpha1,mu0,mu01,SCmass,mass,&
            gammatot,potshift,pot,Lj,bmag,bmag1,btheta,btheta1,phi,&
            kxminbuff,kxmaxbuff
character(len=200),intent(in) :: varphase
real*8,intent(out) :: gap

integer :: jj,nkx,gaptest,sec,insite,Ndata,depthc,kroot,nl
real*8 :: kx,kxmin,kxmax,deps,epsmin,epsmax,slope
real*8 :: newlevelLIST(20)
complex*16 :: detR
!logical :: convtest


kxmin=kxminbuff
kxmax=kxmaxbuff
epsmin=epsminbuff
epsmax=epsmaxbuff
Ndata=Ndatabuff
!depthc=depth

!convtest=.false.

do sec=1,secant
	kx = (kxmax+kxmin)/2d0
	!print*,"eps=",kx
	call computeKxe(Ndata,epsmin,epsmax,dL0,dL1,&
	            gammatot,mu0,mu01,Deltar,Deltar1,Deltai,Deltai1,alpha,alpha1,&
                bmag,bmag1,btheta,btheta1,potshift,SCmass,mass,&
                Opnbr,Opnbr1,pot,Lj,phi,varphase,kx,kroot,newlevelLIST,nl,gaptest,&
                detR,1,insite)
    !print*,"insite=", insite
    !print*,"kxmin,max=", epsmin,epsmax
	if(gaptest.eq.1) then !ie if kx cut the spectrum
		if((insite.eq.1).and.(nl.eq.2)) then !ie if there is a single point of the grid at the bottom of the minimum
			gap=(kxmax+kxmin)/2d0
			!depthc=depthc-1 !light control to avoid loosing yourself far in the depth
			!if (depthc.lt.0) then
			!	!convtest=.true.
			!	exit
			!endif
			deps=(epsmax-epsmin)/dble(Ndata)
			epsmin=newlevelLIST(1)-deps*0.5d0
			epsmax=newlevelLIST(1)+deps*1.5d0
			!print*,"new kxmin,max=",epsmin,epsmax
			Ndata=50
			kxmin=kxminbuff !can't extrapolate a closer min
			kxmax=kx
		else  !!ie insite > 1 because insite=0 is excluded by gaptest=1
			kxmax=kx
			gap=(kxmax+kxmin)/2d0
		endif

	else !ie if kx doesn't cut the spectrum
		if(abs(kxmin-kxmax).le.1d-13) then !doesn't mean the precision on the gap
			!print*,'check'
			gap=kxmin
			exit
		endif
		kxmin=kx
	endif

	!print*,"epsmin,max= ", kxmin,kxmax
enddo

!if (convtest) then
!    gap=(kxmax+kxmin)/2.
!else
!    gap=kxmaxbuff
!    write(30,*) "Gap did not converge. secant is too low for this depth."
!endif

end subroutine dynamical_secante




end module



