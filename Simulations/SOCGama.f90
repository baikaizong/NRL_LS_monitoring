Module constants  
  implicit none 
  integer, parameter :: sample_size=100,simu=10000,ewma=50,level=4 
  real, parameter :: threshhold_a1=1.2,threshhold_a2=4.5,threshhold_a3=9.6,DF=3
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

    open(10,file="F:/Fortran_locationscale/TEXT/SOCGama_906-L1N100")
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
	real score(level)

    call GamaICProb(Prob)
    call Locscore(Prob,score)

	limitl=2.0
	limitR=4.0
    
    do while (abs(arl-370.0) >= 4)
	Limit = (LimitL + LimitR) / 2.0
    simuarl(:)=-1.0
    erre=0
	do i=1,simu
      sample_ewma=sample_size*Prob
	  Test=0.0
	  do j = 1, ewma, 1	
	  call RNGET (iseed)  
      call RNGAM (sample_size, DF, sample)
	  !sample=sample+0.2
	  !print*,sample
      call ordinal_sample(sample,sample_ordinal)
      sample_ewma=lamda*sample_ordinal+co_lamda*sample_ewma
	  !print*,sample_ewma
      call SOCTest(Prob,sample_ewma,score,Test)
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
      call RNGAM (sample_size, DF, sample)
	  !sample=sample+0.2
      call ordinal_sample(sample,sample_ordinal)
      sample_ewma=lamda*sample_ordinal+co_lamda*sample_ewma
      call SOCTest(Prob,sample_ewma,score,Test)
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
	real Score(level)
	real shiftm,shifts

    call GamaICProb(Prob)
    call Locscore(Prob,score)

    simuarl(:)=-1.0
    erre=0
	do i=1,simu
      sample_ewma=sample_size*Prob
	  Test=0.0
	  do j = 1, ewma, 1	
	  call RNGET (iseed)  
      call RNGAM (sample_size, DF, sample) 
	  !print*,sample
      call ordinal_sample(sample,sample_ordinal)
      sample_ewma=lamda*sample_ordinal+co_lamda*sample_ewma
      call SOCTest(Prob,sample_ewma,score,Test)
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
      call RNGAM (sample_size, DF, sample)
	  sample=(sample+shiftm)/(1+shifts)
      call ordinal_sample(sample,sample_ordinal)
      sample_ewma=lamda*sample_ordinal+co_lamda*sample_ewma
      call SOCTest(Prob,sample_ewma,score,Test)
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

subroutine GamaICProb(Prob)
	use IMSL
	use constants
	implicit none
	real Prob(level)
	real interval(level)

	interval(1)=GAMDF(threshhold_a1,DF)
	interval(2)=GAMDF(threshhold_a2,DF)
	interval(3)=GAMDF(threshhold_a3,DF)

	Prob(1)= interval(1)
	Prob(2)=interval(2)-interval(1)
	Prob(3)=interval(3)-interval(2)
	Prob(4)=1-interval(3)

	return
	
end subroutine GamaICProb

subroutine Locscore(Prob,Score)
	use IMSL
	use constants
	implicit none 
	real Prob(level),score(level),Sumscore(level)
    Sumscore(1)=Prob(1)
    Sumscore(2)=Sumscore(1)+Prob(2)
    Sumscore(3)=Sumscore(2)+Prob(3)
    Sumscore(4)=1.0

    score(1)=Sumscore(1)-1.0
    score(2)=Sumscore(1)+Sumscore(2)-1.0
    score(3)=Sumscore(2)+Sumscore(3)-1.0
    score(4)=Sumscore(3)
    return
	
end subroutine Locscore

subroutine SOCTest(Prob,sample_ewma,score,Test)
  	use IMSL
	use constants
	implicit none
    integer i
    real Prob(level),sample_ewma(level)
	real score(level)
	real Test 
	Test=0.0
	do i = 1, level, 1
		Test=Test+sample_ewma(i)*score(i)
	end do
    Test=ABS(Test)

	return

		
end subroutine SOCTest





