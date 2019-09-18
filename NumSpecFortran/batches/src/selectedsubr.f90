module selectedsubr

!!!! all the other subroutines for scanning the spectrum are in the module scans.f90

implicit none
real*8,parameter	:: pi=3.141592653589793d0

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! LINALG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     inverts matrix H of size nn x nn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine invert(H,n)

    implicit none
    integer,intent(in)      ::n
    complex*16,intent(in)   ::H(n,n)
    integer                 :: Info,Ipiv(n)
    complex*16              :: work(n)

    call zgetrf(n,n,H,n,Ipiv,Info)
    call zgetri(n,H,n,Ipiv,work,n,Info)

end subroutine invert


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!routinte to find the determinant of a nxn square matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function finddetN(S,n)

  implicit none
  integer, intent(in)		:: n
  integer	:: ipiv(n), info, i
  complex*16, intent(in)	:: S(n,n)
  complex*16	:: V(n,n)
  complex*16	:: finddetN

  V = S
  call ZGETRF( n, n, V, n, IPIV, INFO )
  
  if (info.NE.0) then
    write(*,*) 'error in subroutine finddetN, info= ',info
  endif
  finddetN = 1d0
  
  do i=1,n
    finddetN = finddetN * V(i,i)
    if (ipiv(i).NE.i) finddetN = - finddetN
  enddo

end function finddetN



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! SCATTERING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     returns concatenation S12 = S1 x S2 of S-matrices S1 and S2
!     S matrices are of the form 
!     S = {t , rp}
!         {r , tp}
!     p stands for primed matrices
!
!     S1 describes a scatterer on the left
!     S2    "                         right
!
!     S -matrices have size 2*nn
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine concat(S1,S2,S12,nn)
    implicit none
    integer, intent(in)     :: nn  ! size of t matrix
    complex*16,intent(in)   :: S1(2*nn,2*nn), S2(2*nn,2*nn)
    complex*16,intent(out)  :: S12(2*nn,2*nn)
    complex*16              :: aux1(nn,nn), aux2(nn,nn), den(nn,nn)
    complex*16              :: r1(nn,nn), rp1(nn,nn), t1(nn,nn), tp1(nn,nn)
    complex*16              :: r2(nn,nn), rp2(nn,nn), t2(nn,nn), tp2(nn,nn)
    complex*16              :: r12(nn,nn),rp12(nn,nn),t12(nn,nn),tp12(nn,nn)
    integer                 :: i, j


    ! take reflection and transmission matrices from S1 and S2
    do i=1,nn
     do j=1,nn
        r1(i,j) = S1(nn+i,j)
        r2(i,j) = S2(nn+i,j)
        t1(i,j) = S1(i,j)
        t2(i,j) = S2(i,j)
        tp1(i,j) = S1(nn+i,nn+j)
        tp2(i,j) = S2(nn+i,nn+j)
        rp1(i,j) = S1(i,nn+j)
        rp2(i,j) = S2(i,nn+j)
     end do
    end do

    den=-dcmplx(1d0,0d0) *matmul(r2,rp1)
!    call minmult(r2,rp1,den,nn)  ! den = - r2 . rp1
    do i=1,nn
     den(i,i) = den(i,i) + dcmplx(1.d0,0.d0)
    end do                      ! den = 1 - r2 . rp1
    call invert(den,nn)          ! den = [1- r2 . rp1]^-1
    aux1=matmul(den,tp2)          ! aux1= [1- r2 . rp1]^-1 .tp2
    tp12=matmul(tp1,aux1)         ! tp12= tp1 .[1- r2 . rp1]^-1 .tp2
    aux2=matmul(rp1,aux1)         ! aux2= rp1. [1- r2 . rp1]^-1 .tp2
    aux1=matmul(t2,aux2)          ! aux1 = t2 . rp1. [1- r2 . rp1]^-1 .tp2
    rp12=rp2+aux1                 ! rp12 = rp2 +  t2 . rp1. [1- r2 . rp1]^-1 .tp2

    den=-dcmplx(1d0,0d0) *matmul(rp1,r2)
!    call minmult(rp1,r2,den,nn)  ! den = - rp1 . r2
    do i=1,nn
    den(i,i) = den(i,i) + dcmplx(1.d0,0.d0)   ! den = 1- rp1 . r2
    end do
    call invert(den,nn)          ! den = [1- rp1 . r2]^-1
    aux1=matmul(den,t1)           ! aux1= [1- rp1 . r2]^-1 . t1
    t12=matmul(t2,aux1)           ! t12= t2 .[1- rp1 . r2]^-1 . t1
    aux2=matmul(r2,aux1)          ! aux2 = r2 .[1- rp1 . r2]^-1 . t1
    aux1=matmul(tp1,aux2)         ! aux1 = tp1 . r2 .[1- rp1 . r2]^-1 . t1
    r12=r1+aux1                   ! r12 = r1 +  tp1 . r2 .[1- rp1 . r2]^-1 . t1

    ! reassemble concatenated S-matrix
    do i=1,nn
        do j=1,nn
            S12(i,j) = t12(i,j)
            S12(i+nn,j) = r12(i,j)
            S12(i,nn+j) = rp12(i,j)
            S12(i+nn,nn+j) = tp12(i,j)
        end do
    end do

end subroutine concat




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! SCATTERING_WIRE !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     returns in dS the S-matrix of a segment with length dL=L/N 
!     with singlet sc order parameter Delta
!     and spin-orbit interaction alpha
!     The matrix has size 2x 2x 2 x nc
!     for in- out- going states, spin, eh-chacter resp.
!     nc is the number of transverse channels
!     
!     the form is:
!     S = {t1, rp }
!         {r , t1p}
!     off-diagonal elements of r originate from Andreev reflection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine deltaSwire_1ch(dS2,nc,eps,dL,gamma1,mu,Deltar, Deltai,alpha,bmag, bmagx, btheta, potshift,mass)
  implicit none
    real*8,intent(in)       ::Deltar, Deltai, alpha, bmag, bmagx,dL, eps, gamma1,btheta
    real*8,intent(in)       ::potshift, mass
    integer, intent(in)     :: nc ! number of transverse channels
    complex*16,intent(out)  ::dS2(8*nc,8*nc)
    real*8                  ::x, y, vF, kFe, kFh, mu, vpot,eso, ww ! comps of the disorder potential and kF
    complex*16              :: iz 
    complex*16              ::Sn(8*nc,8*nc),Sd(8*nc,8*nc),V(8*nc,8*nc),W(8*nc,8*nc)
    integer                 :: i,j
    complex*16              ::tr1pp,tr1pm,tr1mm,tr1mp

    ! write(*,*) "call deltaSAR4"
    ! write(*,*) "nc=", nc
    ! write(*,*) "eps=", eps
    ! write(*,*) "dL=", dL
    ! write(*,*) "st=", st
    ! write(*,*) "|Delta|=", sqrt(Deltar**2 + Deltai**2)
    
    iz = dcmplx(0.d0, 1.d0)

    if (nc .gt. 1) then
       write(*,*) "Only for nc=1."
       stop
    end if

    !vF = sqrt(2.0 * (sqrt(bmag**2 + alpha**4) + alpha**2)) &
    !     - sqrt(2.0 * (sqrt(bmag**2 + alpha**4) - alpha**2)) * alpha**2/bmag
    ! write (*,*) "vF = " , vF

    kFe = sqrt(2d0 * mass * (mu + potshift + eps) )
    kFh = sqrt(2d0 * mass * (mu + potshift - eps) )
    


    
    ! calculate correction to unit matrix
    do i=1,8*nc   ! start with zero
       do j=1,8*nc
          V(i,j) = dcmplx(0.d0,0.d0) ! s-matrix corr. due to Andreev refl
       end do
    end do


    !phases due to finite energies

    ! transmission from the left
    !sigma,s = ++ (e up left moving)
    V(1,1) = kFe * dL 
    !sigma,s = +- (e down left moving)
    V(2,2) = kFe * dL 
    !sigma,s = -+ (h up left moving)
    V(3,3) = - kFh * dL 
    !sigma,s = -- (h down left moving)
    V(4,4) = - kFh * dL 

    ! transmission from the right
    !sigma,s = ++ (e up right moving)
    V(5,5) = kFe * dL 
    !sigma,s = +- (e down right moving)
    V(6,6) = kFe * dL 
    !sigma,s = -+ (h up right moving)
    V(7,7) = - kFh * dL 
    !sigma,s = -- (h down right moving)
    V(8,8) = - kFh * dL 


    !include constant potential to subtract chemical potential
    !mu = kF**2/2.0_8


    call normalen(1,ww)
    !write(*,*) "before ww=", ww, "st", st, "dL", dL, "vF", vF
    !ww = (st/dL * kF) * ww
    ww = sqrt(gamma1*dL)  * ww

    vpot =  potshift

    !transmission from the left
    V(1,1) = V(1,1)  - ( vpot * mass / kFe * dL  + ww * mass / kFe )
    V(2,2) = V(2,2)  - ( vpot * mass / kFe * dL  + ww * mass / kFe )
    V(3,3) = V(3,3)  + ( vpot * mass / kFh * dL  + ww * mass / kFh )
    V(4,4) = V(4,4)  + ( vpot * mass / kFh * dL  + ww * mass / kFh )
    
    !transmission from the right
    V(5,5) = V(5,5)  - ( vpot * mass / kFe * dL  + ww * mass / kFe )
    V(6,6) = V(6,6)  - ( vpot * mass / kFe * dL  + ww * mass / kFe )
    V(7,7) = V(7,7)  + ( vpot * mass / kFh * dL  + ww * mass / kFh )
    V(8,8) = V(8,8)  + ( vpot * mass / kFh * dL  + ww * mass / kFh )

    !reflection from the left
    V(5,1) = V(5,1)  - ( vpot * mass / kFe * dL  + ww * mass / kFe )
    V(6,2) = V(6,2)  - ( vpot * mass / kFe * dL  + ww * mass / kFe )
    V(7,3) = V(7,3)  + ( vpot * mass / kFh * dL  + ww * mass / kFh )
    V(8,4) = V(8,4)  + ( vpot * mass / kFh * dL  + ww * mass / kFh )

    !reflection from the right
    V(1,5) = V(1,5)  - ( vpot * mass / kFe * dL  + ww * mass / kFe )
    V(2,6) = V(2,6)  - ( vpot * mass / kFe * dL  + ww * mass / kFe )
    V(3,7) = V(3,7)  + ( vpot * mass / kFh * dL  + ww * mass / kFh )
    V(4,8) = V(4,8)  + ( vpot * mass / kFh * dL  + ww * mass / kFh )


    !call norm(V,8,ww)
    !write(*,*) "norm2 V=", ww


    ! include singlet pairing correllations

    ! reflection from the left
    V(5,4) = - dcmplx(Deltar, Deltai) / sqrt( kFe * kFh ) * mass * dL
    V(6,3) = dcmplx(Deltar, Deltai) / sqrt( kFe * kFh ) * mass * dL
    V(7,2) = dcmplx(Deltar, -Deltai) / sqrt( kFe * kFh ) * mass * dL
    V(8,1) = - dcmplx(Deltar, -Deltai) / sqrt( kFe * kFh ) * mass * dL  

    ! reflection from the right
    V(1,8) = - dcmplx(Deltar, Deltai) / sqrt( kFe * kFh ) * mass * dL 
    V(2,7) = dcmplx(Deltar, Deltai) / sqrt( kFe * kFh ) * mass * dL
    V(3,6) = dcmplx(Deltar, -Deltai) / sqrt( kFe * kFh ) * mass * dL
    V(4,5) = - dcmplx(Deltar, -Deltai) / sqrt( kFe * kFh ) * mass * dL 

    ! transmission from the left
    V(1,4) = V(1,4) - dcmplx(Deltar, Deltai) / sqrt( kFe * kFh ) * mass * dL 
    V(2,3) = V(2,3) + dcmplx(Deltar, Deltai) / sqrt( kFe * kFh ) * mass * dL
    V(3,2) = V(3,2) + dcmplx(Deltar, -Deltai) / sqrt( kFe * kFh ) * mass * dL
    V(4,1) = V(4,1) - dcmplx(Deltar, -Deltai) / sqrt( kFe * kFh ) * mass * dL

    ! transmission from the right
    V(5,8) = V(5,8) - dcmplx(Deltar, Deltai) / sqrt( kFe * kFh ) * mass * dL 
    V(6,7) = V(6,7) + dcmplx(Deltar, Deltai) / sqrt( kFe * kFh ) * mass * dL
    V(7,6) = V(7,6) + dcmplx(Deltar, -Deltai) / sqrt( kFe * kFh ) * mass * dL
    V(8,5) = V(8,5) - dcmplx(Deltar, -Deltai) / sqrt( kFe * kFh ) * mass * dL



!     ! include magnetic field in the z-direction: B sigma_z tau_z
      ! This is important: Magnetic field along spin-orbit direction
! 
    ! reflection from the left
    V(5,1) = V(5,1) - bmag/kFe * dL*mass
    V(6,2) = V(6,2) + bmag/kFe * dL*mass
    V(7,3) = V(7,3) + bmag/kFh * dL*mass
    V(8,4) = V(8,4) - bmag/kFh * dL*mass

    ! reflection from the right
    V(1,5) = V(1,5) - bmag/kFe * dL*mass
    V(2,6) = V(2,6) + bmag/kFe * dL*mass
    V(3,7) = V(3,7) + bmag/kFh * dL*mass
    V(4,8) = V(4,8) - bmag/kFh * dL*mass


    ! transmission from the left
    V(1,1) = V(1,1) - bmag/kFe * dL*mass
    V(2,2) = V(2,2) + bmag/kFe * dL*mass
    V(3,3) = V(3,3) + bmag/kFh * dL*mass
    V(4,4) = V(4,4) - bmag/kFh * dL*mass

    ! transmission from the right
    V(5,5) = V(5,5) - bmag/kFe * dL*mass
    V(6,6) = V(6,6) + bmag/kFe * dL*mass
    V(7,7) = V(7,7) + bmag/kFh * dL*mass
    V(8,8) = V(8,8) - bmag/kFh * dL*mass


    ! include effective (constant) spin-orbit field in the x-direction 

    ! reflection from the left
    V(5,2) = V(5,2) - bmagx/kFe * dL*mass
    V(6,1) = V(6,1) - bmagx/kFe * dL*mass
    V(7,4) = V(7,4) - bmagx/kFh * dL*mass
    V(8,3) = V(8,3) - bmagx/kFh * dL*mass

    ! reflection from the right
    V(1,6) = V(1,6) - bmagx/kFe * dL*mass
    V(2,5) = V(2,5) - bmagx/kFe * dL*mass
    V(3,8) = V(3,8) - bmagx/kFh * dL*mass
    V(4,7) = V(4,7) - bmagx/kFh * dL*mass


    ! transmission from the left
    V(1,2) = V(1,2) - bmagx/kFe * dL*mass
    V(2,1) = V(2,1) - bmagx/kFe * dL*mass
    V(3,4) = V(3,4) - bmagx/kFh * dL*mass
    V(4,3) = V(4,3) - bmagx/kFh * dL*mass

    ! transmission from the right
    V(5,6) = V(5,6) - bmagx/kFe * dL*mass
    V(6,5) = V(6,5) - bmagx/kFe * dL*mass
    V(7,8) = V(7,8) - bmagx/kFh * dL*mass
    V(8,7) = V(8,7) - bmagx/kFh * dL*mass

!     ! include magnetic field in the x-direction with magnetic phase
! 
!     ! reflection from the left
!     V(5,2) = V(5,2) - exp( dcmplx(0d0,btheta*pi) ) * bmag/kFe * mass * dL
!     V(6,1) = V(6,1) - exp(-dcmplx(0d0,btheta*pi) ) * bmag/kFe * mass * dL
!     V(7,4) = V(7,4) + exp(-dcmplx(0d0,btheta*pi) ) * bmag/kFh * mass * dL
!     V(8,3) = V(8,3) + exp( dcmplx(0d0,btheta*pi) ) * bmag/kFh * mass * dL
! 
!     ! reflection from the right
!     V(1,6) = V(1,6) - exp( dcmplx(0d0,btheta*pi) ) * bmag/kFe * mass * dL
!     V(2,5) = V(2,5) - exp(-dcmplx(0d0,btheta*pi) ) * bmag/kFe * mass * dL
!     V(3,8) = V(3,8) + exp(-dcmplx(0d0,btheta*pi) ) * bmag/kFh * mass * dL
!     V(4,7) = V(4,7) + exp( dcmplx(0d0,btheta*pi) ) * bmag/kFh * mass * dL
! 
! 
!     ! transmission from the left
!     V(1,2) = V(1,2) - exp( dcmplx(0d0,btheta*pi) ) * bmag/kFe * mass * dL
!     V(2,1) = V(2,1) - exp(-dcmplx(0d0,btheta*pi) ) * bmag/kFe * mass * dL
!     V(3,4) = V(3,4) + exp(-dcmplx(0d0,btheta*pi) ) * bmag/kFh * mass * dL
!     V(4,3) = V(4,3) + exp( dcmplx(0d0,btheta*pi) ) * bmag/kFh * mass * dL
! 
!     ! transmission from the right
!     V(5,6) = V(5,6) - exp( dcmplx(0d0,btheta*pi) ) * bmag/kFe * mass * dL
!     V(6,5) = V(6,5) - exp(-dcmplx(0d0,btheta*pi) ) * bmag/kFe * mass * dL
!     V(7,8) = V(7,8) + exp(-dcmplx(0d0,btheta*pi) ) * bmag/kFh * mass * dL
!     V(8,7) = V(8,7) + exp( dcmplx(0d0,btheta*pi) ) * bmag/kFh * mass * dL


    ! include spin-orbit interaction
    ! include BC for the velocity mismatch O1(dL)
    ! h_SO = alpha kF sigma_z
    eso = mass * alpha**2/2.0_8
    ! transmission from the left (positive k for electrons, negative for holes)
    V(1,1) = V(1,1) + ( -alpha + eso / kFe ) * mass * dL
    V(2,2) = V(2,2) + (  alpha + eso / kFe ) * mass * dL
    V(3,3) = V(3,3) + (  alpha - eso / kFh ) * mass * dL
    V(4,4) = V(4,4) + ( -alpha - eso / kFh ) * mass * dL
    
    ! transmission from the left
    V(5,5) = V(5,5) + (  alpha + eso / kFe ) * mass * dL
    V(6,6) = V(6,6) + ( -alpha + eso / kFe ) * mass * dL
    V(7,7) = V(7,7) + ( -alpha - eso / kFh ) * mass * dL
    V(8,8) = V(8,8) + (  alpha - eso / kFh ) * mass * dL

    ! reflection from the left
    V(5,1) = V(5,1) + ( eso / kFe ) * mass * dL
    V(6,2) = V(6,2) + ( eso / kFe ) * mass * dL
    V(7,3) = V(7,3) - ( eso / kFh ) * mass * dL
    V(8,4) = V(8,4) - ( eso / kFh ) * mass * dL

    ! reflection from the right
    V(1,5) = V(1,5) + ( eso / kFe ) * mass * dL
    V(2,6) = V(2,6) + ( eso / kFe ) * mass * dL
    V(3,7) = V(3,7) - ( eso / kFh ) * mass * dL
    V(4,8) = V(4,8) - ( eso / kFh ) * mass * dL


 
    ! generate random Hermitean backscattering matrix
    !call randomV4(W, nc, st)


    ! trick: return unitary matrix
    do i=1,8*nc
        do j=1,8*nc
           Sn(i,j) = iz * 0.5d0 * (V(i,j) ) !+    W(i,j)
           Sd(i,j) = -iz *0.5d0 * (V(i,j) ) !+    W(i,j)
           if (i .EQ. j) then
              Sn(i,j) = Sn(i,j) + dcmplx(1.d0, 0.d0)
              Sd(i,j) = Sd(i,j) + dcmplx(1.d0, 0.d0)
           end if
        end do
    end do
            
    call invert(Sd,8*nc)
    
    dS2=matmul(Sn,Sd)

end subroutine deltaSwire_1ch



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!    returns scattering matrix S of a potential barrier
!    potential in junction and mu are shifted by potshift
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
subroutine Sbarrier(Sb,mu, eps, pot, Lj, potshift, mass, SCmass)
implicit none
  complex*16, intent(out)   :: Sb(8,8)
  real*8                    :: eps, pot, Lj, potshift, mu, mass, SCmass
  integer                   :: i,j
  complex*16                :: ke, qe, den, ree, tee, thh, rhh, denh, keh, qeh

!initialize Sb
  do i=1,8 
     do j=1,8
        Sb(i,j) = dcmplx(0.0_8, 0.0_8)
     end do
  end do
  
  !electrons
  ke = sqrt(2d0*mass* (mu+potshift+eps))
  qe = sqrt(2d0*mass* (mu+potshift+eps-(pot+potshift)))

  if (pot .gt. mu + eps) then
     qe = dcmplx(0.0d0, sqrt(2d0*mass*(-(mu+potshift)-eps+(pot+potshift))))
  end if

  !holes
  keh = sqrt(2d0*mass* (mu+potshift-eps))
  qeh = sqrt(2d0*mass* (mu+potshift-eps-(pot+potshift)))

  if (pot .gt. mu - eps) then
     qeh = dcmplx(0.0d0, sqrt(2d0*mass*(-(mu+potshift)+eps+(pot+potshift))))
  end if


  den = 2.0d0 * dcmplx(0.0d0, 1.0d0) * ke * qe * cos(qe * Lj) &
       & +  (ke**2 + qe**2) * sin(qe * Lj)
  
  denh = 2.0d0 * dcmplx(0.0d0, 1.0d0) * keh * qeh * cos(qeh * Lj) &
       & +  (keh**2 + qeh**2) * sin(qeh * Lj) 

  tee =       2.0d0 * dcmplx(0.0d0, 1.0d0) * ke  * qe / den
  thh = conjg(2.0d0 * dcmplx(0.0d0, 1.0d0) * keh * qeh/ denh )
  ree =       (ke**2  - qe**2 ) * sin(qe  *  Lj) /den
  rhh = conjg((keh**2 - qeh**2) * sin(qeh *  Lj) /denh )

! hard wall reflection matrix
! tee=0d0
! thh=0d0
! ree=-1d0
! rhh=-1d0

! ! transparent wall transmission matrix
! tee=1d0
! thh=1d0
! ree=0d0
! rhh=0d0


  ! add backscattering terms
  Sb(1,1) = Sb(1,1) + tee
  Sb(2,2) = Sb(2,2) + tee
  Sb(3,3) = Sb(3,3) + thh
  Sb(4,4) = Sb(4,4) + thh
  Sb(5,5) = Sb(5,5) + tee
  Sb(6,6) = Sb(6,6) + tee
  Sb(7,7) = Sb(7,7) + thh
  Sb(8,8) = Sb(8,8) + thh

  Sb(5,1) = Sb(5,1) + ree
  Sb(6,2) = Sb(6,2) + ree
  Sb(7,3) = Sb(7,3) + rhh
  Sb(8,4) = Sb(8,4) + rhh
  Sb(1,5) = Sb(1,5) + ree
  Sb(2,6) = Sb(2,6) + ree
  Sb(3,7) = Sb(3,7) + rhh
  Sb(4,8) = Sb(4,8) + rhh

end subroutine Sbarrier



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!    function: returns det[ 1 - S1 x S2 x S3 x Sb x BC]
!    Si are the 8x8 scattering matrices of the semiconductor wire segments
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11

complex*16 function detringconcat(S1,S2,S3,Sb,phi,varphase)
implicit none
complex*16, intent(in)	:: S1(8,8),S2(8,8),S3(8,8),Sb(8,8)
complex*16		:: Smag(8,8),Ssc(8,8),Stot(8,8)
real*8, intent(in)	:: phi
integer			:: i,j
character(len=200)  	::varphase

  call concat_wire(Stot,S1,S2,S3,Sb,phi,varphase)   
  Stot=-Stot
  do i=1,8
     Stot(i,i) = Stot(i,i) + dcmplx(1.0_8, 0.0_8)
  end do

  detringconcat = finddetN(Stot,8)
end function detringconcat


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	concatenate the matrices for the function detringconcat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine concat_wire(Stot,S1,S2,S3,Sb,phi,varphase)
complex*16, intent(in)	:: S1(8,8),S2(8,8),S3(8,8),Sb(8,8)
complex*16		:: Smag(8,8),Ssc(8,8)
complex*16, intent(out)	:: Stot(8,8)
complex*16		:: SphaseL(8,8),SphaseR(8,8)
real*8, intent(in)	:: phi
integer			:: i,j
character(len=200), intent (in)  	::varphase


!initialization
Stot=dcmplx(0d0,0d0)
SphaseL=dcmplx(0d0,0d0)
SphaseR=dcmplx(0d0,0d0)
forall (i=1:8) SphaseL(i,i)=dcmplx(1.0_8,0.0_8)
SphaseR=SphaseL
Stot=SphaseL


call init_phase_mat(Smag,Ssc,phi)


if (varphase.eq.'scL') then 
      SphaseL = Ssc
      call concat(SphaseL,Ssc,SphaseL,4)
      call invert(SphaseL,8)
elseif ( (varphase.eq.'scR').or.(varphase.eq.'sc2d') ) then 
      SphaseR = Ssc
      call concat(SphaseR,Ssc,SphaseR,4)
elseif (varphase.eq.'scLR') then  !!!!!!!!! our system configuration
      SphaseL = Ssc
      SphaseR = Ssc
elseif (varphase.eq.'magL') then
      SphaseL = Smag
      call concat(SphaseL,Smag,SphaseL,4)
      call invert(SphaseL,8)
elseif ( (varphase.eq.'magR').or.(varphase.eq.'mag2d') ) then
      SphaseR = Smag
      call concat(SphaseR,Smag,SphaseR,4)
elseif (varphase.eq.'magLR') then
      SphaseL = Smag
      SphaseR = Smag
else
      write(*,*) 'wrong parameter varphase'
      stop
end if


!concatenation of segments
     call concat(Stot,SphaseL,Stot,4)
     call concat(Stot,S2,Stot,4)
     call concat(Stot,SphaseR,Stot,4)
     call concat(Stot,S3,Stot,4)
     call concat(Stot,Sb,Stot,4)
     call concat(Stot,S1,Stot,4)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	phase matrices for:
!	* magnetic phase (not used)
!	* SC phase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_phase_mat(Smag,Ssc,phi)
  implicit none
  complex*16, intent(out)	:: Smag(8,8), Ssc(8,8)
  real*8, intent(in)		:: phi
  
  
Smag=dcmplx(0d0,0d0)
Ssc=dcmplx(0d0,0d0)
  
! scattering matrix for half the magnetic phase jump: 
! operator =exp(i phase/2 sigma_z tau_z)
! electron phase jump
  Smag(1,1) = exp( dcmplx(0.0_8,  -1.0_8 * pi * phi) )
  Smag(2,2) = exp( dcmplx(0.0_8,   1.0_8 * pi * phi) )
! hole phase jump
  Smag(3,3) = exp( dcmplx(0.0_8,   1.0_8 * pi * phi) )
  Smag(4,4) = exp( dcmplx(0.0_8,  -1.0_8 * pi * phi) )

  ! opposite b.c.
  Smag(5,5) = exp( dcmplx(0.0_8,   1.0_8 * pi * phi) )
  Smag(6,6) = exp( dcmplx(0.0_8,  -1.0_8 * pi * phi) )
  Smag(7,7) = exp( dcmplx(0.0_8,  -1.0_8 * pi * phi) )
  Smag(8,8) = exp( dcmplx(0.0_8,   1.0_8 * pi * phi) )
  
! scattering matrix for half the superconducting phase jump
! operator =exp(i phase/2 tau_z)
! electron phase jump
  Ssc(1,1) = exp( dcmplx(0.0_8,  -1.0_8 * pi * phi) )
  Ssc(2,2) = exp( dcmplx(0.0_8,  -1.0_8 * pi * phi) )
! hole phase jump
  Ssc(3,3) = exp( dcmplx(0.0_8,   1.0_8 * pi * phi) )
  Ssc(4,4) = exp( dcmplx(0.0_8,   1.0_8 * pi * phi) )

  ! opposite b.c.
  Ssc(5,5) = exp( dcmplx(0.0_8,   1.0_8 * pi * phi) )
  Ssc(6,6) = exp( dcmplx(0.0_8,   1.0_8 * pi * phi) )
  Ssc(7,7) = exp( dcmplx(0.0_8,  -1.0_8 * pi * phi) )
  Ssc(8,8) = exp( dcmplx(0.0_8,  -1.0_8 * pi * phi) )  

end subroutine




end module
