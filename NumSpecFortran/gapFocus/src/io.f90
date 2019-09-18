module io

implicit none

real*8          :: dL, dL0, dL1, dL2, epsmin,epsmax,dkfrac
real*8		:: Lmax,bmag,btheta,alpha,Deltamag,Deltaphase,mu0
real*8		:: Lmax1,bmag1,btheta1,alpha1,Deltamag1,Deltaphase1,mu01
real*8		:: Deltaphase2
real*8          :: potshift, pot, Lj, SCmass, mass
integer         :: Nseg,Nseg1,Ndata
integer		:: nmu,flatness
character(len=200) ::outfname,varphase

character(len=200) :: Aphilist,ABlist,AmuSClist,Amulist
integer :: Anstepphi,AnstepB,AnstepmuSC,Anstepmu
character(len=200) :: dataname


!!!! splitting subroutine split_read_real modified from https://groups.google.com/forum/#!topic/comp.lang.fortran/UFGS4c7UGKg by Thomas Koenig 12/02/2018 !!!!
  integer, private :: ii
  private :: is_sep


namelist /parameters/ dL, dL0, dL1, dL2, epsmin,epsmax,Lmax,bmag,btheta,alpha,Deltamag,Deltaphase,mu0,&
		Lmax1,bmag1,btheta1,alpha1,Deltamag1,Deltaphase1,mu01,Deltaphase2,potshift, pot,&
		Lj, SCmass, mass,Nseg,Nseg1,Ndata,nmu, dkfrac,flatness,outfname,varphase

namelist /newparameters/ dataname,Anstepphi,AnstepB,AnstepmuSC,Anstepmu,&
		Aphilist,ABlist,AmuSClist,Amulist


contains


subroutine get_infile_parameter()
    implicit none
    character(len=100) :: infile_name

    call get_command_argument(1,infile_name) ! on the command line ./JJ_bands a b c d e , infile_name takes the value of the 1st argument (0 is the invoked program): a

    if (len_trim(infile_name) .eq. 0) then
        write (*,*) "#Please insert name of input file"
        stop
    end if

    open(42,file=trim(infile_name))
    read(42,parameters)
    close(42)
end subroutine



subroutine get_infile_newparameters()
    implicit none
    character(len=100) :: infile_name

    call get_command_argument(2,infile_name)
    if (len_trim(infile_name) .eq. 0) then
        write (*,*) "#Please insert name of NEW input file"
        stop
    end if

    open(17,file=trim(infile_name))
    read(17,nml=newparameters)
    close(17)
end subroutine



subroutine definelooprange(rangemin,rangemax,nstep,step,nstepnew,stepnew,generatedarray,flag)
	!! Generates the array GENERATEDARRAY of a parameter 'flag'
	!! beginning with RANGEMIN (included)
	!! finishing with RANGEMAX (included)
	!! with NSTEP steps in total.
	!! If NSTEP==0, then the array begins with RANGEMIN, adding STEP until it reaches at most RANGEMAX.
	!! The updated values of NSTEP and STEP, describing the final returned array, are available in NSTEPNEW and STEPNEW.

	!! Constraints: RANGEMIN =< RANGEMAX & STEP >0

	implicit none
	real*8,intent(in)	:: rangemax, rangemin, step
	integer,intent(in)	:: nstep
	character(len=*),intent(in)	:: flag
	real*8,intent(out)	:: stepnew
	integer,intent(out)	:: nstepnew
	real*8,allocatable,dimension(:),intent(out)	:: generatedarray
	integer			:: i
	
	nstepnew=nstep
	stepnew=step

	if (rangemax.eq.rangemin) then
			nstepnew=1
			stepnew=0.
	else if (rangemax.lt.rangemin) then
		write(*,*) "#Positive parameters required, and MAX >= MIN, in parameter ", flag,'.'
		stop
	else
		if (nstep.gt.1) then
			stepnew=(rangemax-rangemin)/(nstep-1)
		else if (nstep.eq.1) then
			stepnew=0.
		else if ((step.ne.0d0).and.(step.le.(rangemax-rangemin))) then
			nstepnew=idnint((rangemax-rangemin)/step)+1 !round it to remove machine error, instead of int
		else
			write(*,*) "#Need a number of steps or a smaller step in parameter ", flag,'.'
			stop
		end if
	end if
	
	!allocate(generatedarray(nstepnew))
	generatedarray=(/ (rangemin+stepnew*i,i=0,nstepnew-1) /)
end subroutine definelooprange


subroutine makelist(arrayin,arrayout,flag)
	!! Generate a list of numerical values ARRAYOUT for the variable FLAG (character)
	!! reading the instructions in ARRAYIN (character):
	!! 	* a simple float isolated by commas indicates a single value
	!!		ex: arrayin='[... ,]2.5[, ...]'
	!!	* 3 floats(!) isolated by colons are read as (min:step:max) with min included
	!!		ex: arrayin='[... ,]1.:0.3:3.[, ...]'
	!!	* 2 floats around 1 integer(!) isolated by colons are read as (min:Nstep:max)
	!!		with min and max included. Ex: arrayin='[... ,]1.:3:3.[, ...]'
	!! Exemple: arrayin='0.:5:10.,12,13.5,0.:5.:10.' gives the following output,
	!!		arrayout=(/ 0. 2.5 5. 7.5 10. 12. 13.5 0. 5. 10. /)

	!! Attention: The central element of a generating sequence (:) is read as float or integer
	!!		depending on the presence of a dot (.) or not.
	!! Constraints: For the generating sequence (:), see subroutine definelooprange.


	implicit none
	character(len=*),intent(in) :: arrayin,flag
	real*8,allocatable,dimension(:),intent(out) :: arrayout
	character(len=20),allocatable,dimension(:) :: subarray,subsubarray
	character(len=20) :: seq
	integer :: n,i,nn,ntot,ind

	real*8 :: rangemin, rangemax, step, stepnew
	integer :: nstep,nstepnew,testfloat
	real*8, allocatable, dimension(:) :: prelist

	!read input - first level
	call stringsplit(arrayin,',',n,subarray)

	!compute the total length of the list
	ntot=0
	do i=1,n !loop over each sequence
		seq=subarray(i)
		call split_read_real(seq,',:',nn,subsubarray) !read sequence - second level
		if (nn.eq.1) then !isolated number
			ntot=ntot+1
		else if (nn.eq.3) then !list generated by the subroutine definelooprange
			read(subsubarray(1),*) rangemin
			read(subsubarray(3),*) rangemax
			testfloat=count(transfer(trim(subsubarray(2)),'a',len_trim(subsubarray(2)))== '.')
			if (testfloat.eq.0) then !analyse if the central number is float or integer
				read(subsubarray(2),*) nstep
				step=0.
			else if (testfloat.eq.1) then
				read(subsubarray(2),*) step
				nstep=0
			else
				print*,'Arrays of variable miswritten.'
				stop
			endif
			call definelooprange(rangemin, rangemax, nstep, step, nstepnew, stepnew, prelist, flag)
			ntot=ntot+nstepnew
		else
			print*,'Arrays of variable miswritten.'
			stop
		endif
		deallocate(subsubarray)
	enddo
	
	! fill in the list (copy-pasted from above with attribution in arrayout)
	allocate(arrayout(ntot))
	ind=1
	do i=1,n
		seq=subarray(i)
		call split_read_real(seq,',:',nn,subsubarray)
		if (nn.eq.1) then
			read(subsubarray(1),*) arrayout(ind)
			ind=ind+1
		else if (nn.eq.3) then
			read(subsubarray(1),*) rangemin
			read(subsubarray(3),*) rangemax
			testfloat=count(transfer(trim(subsubarray(2)),'a',len_trim(subsubarray(2)))== '.')
			if (testfloat.eq.0) then
				read(subsubarray(2),*) nstep
				step=0.
			else if (testfloat.eq.1) then
				read(subsubarray(2),*) step
				nstep=0
			else
				print*,'Arrays of variable miswritten.'
				stop
			endif
			call definelooprange(rangemin, rangemax, nstep, step, nstepnew, stepnew, prelist, flag)
			arrayout(ind:ind+nstepnew-1)=prelist
			ind=ind+nstepnew
		else
			print*,'Arrays of variable miswritten.'
			stop
		endif
		deallocate(subsubarray)
	enddo		
end subroutine makelist


subroutine stringsplit(str,separator,length,spstr) !string separator with blanks and commas
	character(len=*),intent(in) :: str,separator
	integer,intent(out) :: length
	character(len=20),allocatable,dimension(:),intent(out) :: spstr
	
	length=count(transfer(trim(str),'a',len_trim(str))== separator)+1
	allocate(spstr(length))
	read(str,*) spstr(1:length)
end subroutine stringsplit


!! copied/modified split routine
subroutine split_read_real(str,sep,length, aa)
    character(len=*), intent(in) :: str,sep
    integer,intent(out) :: length
    character(len=20), dimension(:), allocatable :: aa
    integer :: n,i
    integer :: from, to, mylen
    n = 1
    mylen = len_trim(str)
    do i=1, mylen
       if (is_sep(str(i:i),sep)) n = n + 1
    end do
    length=n
    allocate (aa(n))
    n = 1
    from = 0
    to = 1
    do while (to <= mylen)
       if (is_sep(str(to:to),sep)) then
          read (unit=str(from+1:to-1),fmt=*) aa(n)
          n = n + 1
          from = to
       end if
       to = to + 1
    end do
    read (unit=str(from+1:mylen),fmt=*) aa(n)
end subroutine split_read_real


logical function is_sep(c,sep)
    character(len=1), intent(in) :: c
    character(len=*),intent(in) :: sep
    logical, dimension(0:255) :: fld
    
    fld = merge(.true.,.false., [(index(sep,achar(ii))>0,ii=0,255)])

    is_sep = fld(ichar(c(1:1)))
end function is_sep



end module
