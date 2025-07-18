Module constants  
  implicit none 
  integer, parameter :: sample_size=100,simu=10000,ewma=50,level=4 
  real, parameter :: threshhold_a1=-2.5,threshhold_a2=0.3,threshhold_a3=2.7
  real, parameter :: lamda=0.1,co_lamda=0.9
End module constants

program main

    use IMSL
	use constants
	implicit none

	real Limit,arl,std
	integer erre
	real dClkStart, dClkFinish
	integer i,j
	real shiftm,shifts

    open(10,file="F:/Fortran_locationscale/TEXT/GpLogic_831-L1N100")
    call cpu_time(dClkStart)
  
	call limitB(Limit,arl,std)
	write(10,*),"Limit is",Limit
	write(10,*),"arl is",arl
	write(10,*),"std is",std
    
	shifts=0.0

	shiftm=0.02

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Location shift is ",shiftm
	write(10,*),"arl is",arl
	write(10,*),"std is",std

	shiftm=0.05

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Location shift is ",shiftm
	write(10,*),"arl is",arl
	write(10,*),"std is",std

	shiftm=0.10

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Location shift is ",shiftm
	write(10,*),"arl is",arl
	write(10,*),"std is",std

	shiftm=0.20

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Location shift is ",shiftm
	write(10,*),"arl is",arl
	write(10,*),"std is",std
	call cpu_time(dClkFinish)

	shiftm=-0.02

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Location shift is ",shiftm
	write(10,*),"arl is",arl
	write(10,*),"std is",std

	shiftm=-0.05

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Location shift is ",shiftm
	write(10,*),"arl is",arl
	write(10,*),"std is",std

	shiftm=-0.10

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Location shift is ",shiftm
	write(10,*),"arl is",arl
	write(10,*),"std is",std

	shiftm=-0.20

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Location shift is ",shiftm
	write(10,*),"arl is",arl
	write(10,*),"std is",std
	call cpu_time(dClkFinish)
	shiftm=0.0
	shifts=0.02

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Scale shift is ",shifts
	write(10,*),"arl is",arl
	write(10,*),"std is",std

	shifts=0.05

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Scale shift is ",shifts
	write(10,*),"arl is",arl
	write(10,*),"std is",std
	
	shifts=0.10

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Scale shift is ",shifts
	write(10,*),"arl is",arl
	write(10,*),"std is",std


	shifts=0.20

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Scale shift is ",shifts
	write(10,*),"arl is",arl
	write(10,*),"std is",std

	shifts=-0.02

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Scale shift is ",shifts
	write(10,*),"arl is",arl
	write(10,*),"std is",std

	shifts=-0.05

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Scale shift is ",shifts
	write(10,*),"arl is",arl
	write(10,*),"std is",std
	
	shifts=-0.10

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Scale shift is ",shifts
	write(10,*),"arl is",arl
	write(10,*),"std is",std
	write(10,*) "Time is", dClkFinish-dClkStart, "seconds"	

	shifts=-0.20

	call OCARL(Limit,shiftm,shifts,arl,std)
    write(10,*),"Scale shift is ",shifts
	write(10,*),"arl is",arl
	write(10,*),"std is",std	
	write(10,*) "Time is", dClkFinish-dClkStart, "seconds"	
    stop

end program main

subroutine limitB(Limit,arl,std)

	use IMSL
	use constants
	implicit none

    real simuarl(simu)
    real Limit,LimitL,LimitR
    integer i,j,k,erre,iseed
    real Test
	real arl,st,std
	real sample_ordinal(level),sample_ewma(level),sample(sample_size)
	real Prob(level)
	real Sigma(level,level)

    call LogICProb(Prob)
    call Sigmafun(Prob,Sigma)




	limitl=0.4
	limitR=1.4
    
    do while (abs(arl-370.0) >= 4)
	Limit = (LimitL + LimitR) / 2.0
    simuarl(:)=-1.0
    erre=0
	do i=1,simu
      sample_ewma=sample_size*Prob
	  Test=0.0
	  do j = 1, ewma, 1	
	  call RNGET (iseed)  
      call LogicRm(sample)
      call ordinal_sample(sample,sample_ordinal)
      sample_ewma=lamda*sample_ordinal+co_lamda*sample_ewma
	  !print*,sample_ewma
      call GpTest(Prob,sample_ewma,Sigma,Test)
	  !print*,Test
	  if (Test>Limit) then
	      exit
	  endif
	  end do
	  if (Test>Limit) then
	      erre=erre+1
	      cycle
	  endif
      simuarl(i)=0.0
	  do while (Test<Limit)
      call RNGET (iseed)  
      call LogicRm(sample)
      call ordinal_sample(sample,sample_ordinal)
      sample_ewma=lamda*sample_ordinal+co_lamda*sample_ewma
      call GpTest(Prob,sample_ewma,Sigma,Test)
	  !print*,Test
	  simuarl(i)=simuarl(i)+1.0
     end do
	print*,i
	print*,simuarl(i)
	print*,"****************"
    end do
    arl=sum(simuarl)/(simu-erre)
    write(10,*),"***********************"
	write(10,*),"Limit is",limit
	write(10,*),"IC ARL is", arl
	if (arl < 370.0) then
	    LimitL = Limit
     else
	    LimitR = Limit
     end if
	 enddo

	st=0.0                                                                                                                                                                                                                                                                                                                                                                            
	do j=1,simu
		if (simuarl(j)==-1.0) then
			cycle
		endif
	 	st=st+(simuarl(j)-arl)**2.0
	enddo 
	std=sqrt(st/(simu-1.0-erre))
	std=std/(sqrt(dble(simu)-1-dble(erre)))	
	return 
end subroutine limitB
subroutine OCARL(Limit,shiftm,shifts,arl,std)

	use IMSL
	use constants
	implicit none

    real simuarl(simu)
    real Limit,LimitL,LimitR
    integer i,j,k,erre,iseed
    real Test
	real arl,st,std
	real sample_ordinal(level),sample_ewma(level),sample(sample_size)
	real Prob(level)
	real Sigma(level,level)
	real shiftm,shifts

    call LogICProb(Prob)
    call Sigmafun(Prob,Sigma)

    simuarl(:)=-1.0
    erre=0
	do i=1,simu
      sample_ewma=sample_size*Prob
	  Test=0.0
	  do j = 1, ewma, 1	
	  call RNGET (iseed)  
      call LogicRm(sample) 
      call ordinal_sample(sample,sample_ordinal)
      sample_ewma=lamda*sample_ordinal+co_lamda*sample_ewma
      call GpTest(Prob,sample_ewma,Sigma,Test)
	  if (Test>Limit) then
	      exit
	  endif
	  end do
	  if (Test>Limit) then
	      erre=erre+1
	      cycle
	  endif
      simuarl(i)=0.0
	  do while (Test<Limit)
      call RNGET (iseed)  
      call LogicRm(sample)
	  sample=(sample+shiftm)/(1+shifts)
      call ordinal_sample(sample,sample_ordinal)
      sample_ewma=lamda*sample_ordinal+co_lamda*sample_ewma
      call GpTest(Prob,sample_ewma,Sigma,Test)
	  !print*,Test
	  simuarl(i)=simuarl(i)+1.0
     end do
    print*,i
	print*,simuarl(i)
	print*,"****************"
    end do
    arl=sum(simuarl)/(simu-erre)
	
	st=0.0                                                                                                                                                                                                                                                                                                                                                                            
	do j=1,simu
		if (simuarl(j)==-1.0) then
			cycle
		endif
	 	st=st+(simuarl(j)-arl)**2.0
	enddo 
	std=sqrt(st/(simu-1.0-erre))
	std=std/(sqrt(dble(simu)-1-dble(erre)))	
	return 
end subroutine OCARL

subroutine LogicRm(sample)
	use IMSL
	use constants
	implicit none
	real sample(sample_size)
	integer i

	call RNUN(sample_size,sample)
	do i = 1, sample_size, 1
		sample(i)=-log(1.0/sample(i)-1)
	end do
    return
end subroutine LogicRm
!Classify X^* as X

subroutine ordinal_sample(sample,sample_ordinal)
	use IMSL
	use constants
	implicit none
	real sample(sample_size)
	real sample_ordinal(level)
	integer i
    
    sample_ordinal(:)=0
    do i = 1, sample_size, 1
        if ( sample(i)<=threshhold_a1 ) then
	 	   sample_ordinal(1)=sample_ordinal(1)+1
	 	elseif( sample(i)>threshhold_a1.and.sample(i)<=threshhold_a2 ) then
           sample_ordinal(2)=sample_ordinal(2)+1
        elseif( sample(i)>threshhold_a2.and.sample(i)<=threshhold_a3 ) then
           sample_ordinal(3)=sample_ordinal(3)+1
        else
           sample_ordinal(4)=sample_ordinal(4)+1
	    end if 
    end do
	
    return
end subroutine ordinal_sample

function LogCDF(a)
	implicit none
	real a
	real LogCDF
	LogCDF=1.0/(1+exp(-a))
	return
end function LogCDF

subroutine LogICProb(Prob)
	use IMSL
	use constants
	implicit none
	real Prob(level)
	real interval(level)
	real LogCDF

	interval(1)=LogCDF(threshhold_a1)
	interval(2)=LogCDF(threshhold_a2)
	interval(3)=LogCDF(threshhold_a3)

	Prob(1)= interval(1)
	Prob(2)=interval(2)-interval(1)
	Prob(3)=interval(3)-interval(2)
	Prob(4)=1-interval(3)

	return
	
end subroutine LogICProb

subroutine Sigmafun(Prob,Sigma)
    use IMSL
	use constants
	implicit none

	real Sigma(level,level),Sigmap(level,level),Sigmaq(level,level)
	real Prob(level)
	integer i,j

	do i = 1, level, 1
	Sigmap(i,i)=Prob(i)
	end do
    
    do i = 1, level, 1
       do j = 1, level, 1
       	Sigmaq(i,j)=Prob(i)*Prob(j)
       end do
    end do
    
    Sigma=Sigmap-Sigmaq
	Sigma= .i. Sigma
	return	
end subroutine Sigmafun


subroutine GpTest(Prob,sample_ewma,Sigma,Test)
  	use IMSL
	use constants
	implicit none
    integer i,j
    real Prob(level),sample_ewma(level)
    real Probm(1,level),sample_ewmam(1,level)
    real var(1,level),tanvesvar(level,1)
	real Sigma(level,level)
	real Testm(1,1),Test 

    Probm(1,:)=Prob
    sample_ewmam(1,:)=sample_ewma
    var=sample_ewmam-sample_size*Probm
    tanvesvar= .t. var
    Testm=var .x. Sigma .x. tanvesvar
    Test=Testm(1,1)/(1.0*sample_size)
	return

		
end subroutine GpTest





