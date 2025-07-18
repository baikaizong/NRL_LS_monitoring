Module constants  
  implicit none 
  integer, parameter :: sample_size=100,shapesquare=6,simu=10000,ewma=50,factor=3,cells=36 
  integer, parameter :: levelone=3,leveltwo=3,levelthree=4
  real, parameter :: lamda=0.1,co_lamda=0.9,shape=3.0
  real, parameter :: TOL=0.00001
  real level(3)
  data level /3,3,4/

End module constants
program main

  use IMSL
  use constants
  implicit none

  real Limit,arl,std
  integer erre
  real dClkStart, dClkFinish
  integer i,j
  real sample(cells),sampleIC(cells)
  real shift(3),Cov(3,3)
  real mProb(cells)
  real Sigmascore(cells,cells)
  real Prob(factor,levelthree)

  do i = 1, factor, 1
    do j = 1, factor, 1
       Cov(i,j)=0.5**(ABS(i-j))
    end do    
  end do 

  open(11,file="E:/Fortran_locationscale/TXT/Probgama.txt")
  do i = 1, cells, 1
  read(11,*),mProb(i)   
  end do
  close(11)
  Limit=0.9617188 
  call ICProb(Prob)
  call scoreSigma(Prob,mProb,Sigmascore)
  open(10,file="E:/Fortran_locationscale/TXT/MOCgamaOCcase_103-1221.txt")
  call cpu_time(dClkStart)

  write(10,*),"limit is",limit
  shift(:)=1.0
  Cov(1,2)=0.52
  Cov(2,1)=0.52
  call OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"Cov(1,2) is ",Cov(1,2)
  write(10,*),"arl is",arl
  write(10,*),"std is",std

  Cov(1,2)=0.55
  Cov(2,1)=0.55
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"Cov(1,2) is ",Cov(1,2)
  write(10,*),"arl is",arl
  write(10,*),"std is",std

  Cov(1,2)=0.60
  Cov(2,1)=0.60
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"Cov(1,2) is ",Cov(1,2)
  write(10,*),"arl is",arl
  write(10,*),"std is",std

  Cov(1,2)=0.70
  Cov(2,1)=0.70
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"Cov(1,2) is ",Cov(1,2)
  write(10,*),"arl is",arl
  write(10,*),"std is",std

 

  Cov(1,2)=0.48
  Cov(2,1)=0.48
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"Cov(1,2) is ",Cov(1,2)
  write(10,*),"arl is",arl
  write(10,*),"std is",std

   Cov(1,2)=0.45
   Cov(2,1)=0.45
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"Cov(1,2) is ",Cov(1,2)
  write(10,*),"arl is",arl
  write(10,*),"std is",std

  Cov(1,2)=0.40
  Cov(2,1)=0.40
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"Cov(1,2) is ",Cov(1,2)
  write(10,*),"arl is",arl
  write(10,*),"std is",std

  Cov(1,2)=0.30
  Cov(2,1)=0.30
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"Cov(1,2) is ",Cov(1,2)
  write(10,*),"arl is",arl
  write(10,*),"std is",std


  call cpu_time(dClkFinish)
  write(10,*) "Time is", dClkFinish-dClkStart, "seconds"  
  stop

end program main

subroutine msampleIC(sample)
  use IMSL
  use constants
  implicit none

  real mnsample(sample_size,factor),mnsamplenormal(shapesquare,factor),sample(cells)
  real Cov(factor,factor),RSIG(factor,factor),threshhold(3,4) 
  integer irank,iseed
  integer cell(sample_size,factor)
  real Wis(factor,factor)
  integer i,j,k,u,v,x,w

   threshhold(1,1)=3.0
   threshhold(1,2)=6.0
   threshhold(1,3)=1000000.0
   threshhold(2,1)=4.0
   threshhold(2,2)=7.0
   threshhold(2,3)=1000000.0
   threshhold(3,1)=3.0
   threshhold(3,2)=5.0
   threshhold(3,3)=7.0
   threshhold(3,4)=1000000.0

  do i = 1, factor, 1
    do j = 1, factor, 1
       Cov(i,j)=0.5**(ABS(i-j))
    end do    
  end do 
  
  call CHFAC (factor,Cov,factor,TOL,irank,RSIG,factor)  
  do i = 1, sample_size, 1
    call RNGET (iseed)   
    call RNMVN(shapesquare,factor,RSIG,factor,mnsamplenormal,shapesquare)
    Wis=(.t. mnsamplenormal) .x. mnsamplenormal 
    do j = 1, factor, 1
      mnsample(i,j)=wis(j,j)/2.0
    end do
  end do

  cell(:,:)=1
    do i = 1, sample_size, 1
    	do j = 1, factor, 1
    		do k = 1, level(j), 1
    			if ( mnsample(i,j)>threshhold(j,k) ) then
    				cell(i,j)=cell(i,j)+1
				else
				    exit
    			end if
    		end do    		
    	end do       
     end do 
    
     sample(:)=0.0
     k=0
     do u = 1, level(1), 1
     	do v = 1, level(2), 1
     		do x = 1, level(3), 1
			        k=k+1
     				do i = 1, sample_size, 1
     					if ( cell(i,1)==u .and. cell(i,2)==v .and. cell(i,3)==x ) then
     						sample(k)=sample(k)+1.0
     					end if	
     				end do
     		end do
     	end do
     end do

	 return
	 
end subroutine msampleIC
subroutine msample(shift,Cov,sample)
  use IMSL
  use constants
  implicit none

  real mnsample(sample_size,factor),mnsamplenormal(shapesquare,factor),sample(cells)
  real Cov(factor,factor),RSIG(factor,factor),threshhold(3,4) 
  integer irank,iseed
  integer cell(sample_size,factor)
  integer i,j,k,u,v,x,w
  real shift(factor)
  real Wis(factor,factor)

   threshhold(1,1)=3.0
   threshhold(1,2)=6.0
   threshhold(1,3)=1000000.0
   threshhold(2,1)=4.0
   threshhold(2,2)=7.0
   threshhold(2,3)=1000000.0
   threshhold(3,1)=3.0
   threshhold(3,2)=5.0
   threshhold(3,3)=7.0
   threshhold(3,4)=1000000.0

  
  call CHFAC (factor,Cov,factor,TOL,irank,RSIG,factor)  
  do i = 1, sample_size, 1
    call RNGET (iseed)   
    call RNMVN(shapesquare,factor,RSIG,factor,mnsamplenormal,shapesquare)
    Wis=(.t. mnsamplenormal) .x. mnsamplenormal 
    do j = 1, factor, 1
      mnsample(i,j)=(wis(j,j)/2.0)
    end do
	mnsample(i,:)=mnsample(i,:)
  end do

  cell(:,:)=1
    do i = 1, sample_size, 1
      do j = 1, factor, 1
        do k = 1, level(j), 1
          if ( mnsample(i,j)>threshhold(j,k) ) then
            cell(i,j)=cell(i,j)+1
        else
            exit
          end if
        end do        
      end do       
     end do 
    
     sample(:)=0.0
     k=0
     do u = 1, level(1), 1
      do v = 1, level(2), 1
        do x = 1, level(3), 1
              k=k+1
            do i = 1, sample_size, 1
              if ( cell(i,1)==u .and. cell(i,2)==v .and. cell(i,3)==x ) then
                sample(k)=sample(k)+1.0
              end if  
            end do
        end do
      end do
     end do

   return
   
end subroutine msample
subroutine ICProb(SProb)
	use IMSL
	use constants
	implicit none
	real SProb(factor,levelthree),threshhold(3,4) 
	real interval(factor,6)
	integer i,j

   threshhold(1,1)=3.0
   threshhold(1,2)=6.0
   threshhold(1,3)=1000000.0
   threshhold(2,1)=4.0
   threshhold(2,2)=7.0
   threshhold(2,3)=1000000.0
   threshhold(3,1)=3.0
   threshhold(3,2)=5.0
   threshhold(3,3)=7.0
   threshhold(3,4)=1000000.0
     
    do i = 1,factor, 1
    	j=1
    	interval(i,j)=0.0
    	do j = 2, level(i), 1
    		interval(i,j)=GAMDF(threshhold(i,j-1),shape)
    	end do
    	interval(i,j)=1.0  	
    end do  
    do i = 1, factor, 1
       do j = 1, level(i), 1
       	  Sprob(i,j)=interval(i,j+1)-interval(i,j)
       end do
    end do

	return
	
end subroutine ICProb


subroutine Locuniscore(Prob,Score)
	use IMSL
	use constants
	implicit none 
	real Prob(factor,levelthree),score(factor,levelthree),Sumscore(factor,levelthree)
	integer i,j
    Sumscore(:,:)=0.0
    do i = 1, factor, 1
    	j=1
    	Sumscore(i,j)=Prob(i,j)
    	do j = 2, level(i), 1
    		Sumscore(i,j)=Sumscore(i,j-1)+Prob(i,j)
    	end do 	
    end do
    do i = 1, factor, 1
    	j=1
    	score(i,j)=Sumscore(i,j)-1.0
    	do j = 2, level(i), 1
    		score(i,j)=Sumscore(i,j-1)+Sumscore(i,j)-1.0
    	end do 	
    end do
    return
	
end subroutine Locuniscore

function Scorescale(a)
	real a
	real Scorescale
	
	Scorescale=(a**2-a)*log(a/(1-a))
	return
end function Scorescale

subroutine Scauniscore(Prob,Score)
	use IMSL
	use constants
	implicit none 
	real Prob(factor,levelthree),score(factor,levelthree),Sumscore(factor,levelthree)
	real Scorescale
	integer i,j

    do i = 1, factor, 1
    	j=1
    	Sumscore(i,j)=Prob(i,j)
    	do j = 2, level(i), 1
    		Sumscore(i,j)=Sumscore(i,j-1)+Prob(i,j)
    	end do 	
    end do


	do i = 1, factor, 1
		j=1
		score(i,j)=Scorescale(Sumscore(i,j))/Prob(i,j)
		do j = 2, level(i)-1, 1
			score(i,j)=(Scorescale(Sumscore(i,j))-Scorescale(Sumscore(i,j-1)))/Prob(i,j)
		end do
			score(i,j)=(-Scorescale(Sumscore(i,j-1)))/Prob(i,j)		
	end do
    return	
end subroutine Scauniscore


subroutine Sigmafun(mProb,Sigma)
    use IMSL
	use constants
	implicit none

	real Sigma(cells,cells),Sigmap(cells,cells),Sigmaq(cells,cells)
	real mProb(cells)
	integer i,j

	do i = 1, cells, 1
	Sigmap(i,i)=mProb(i)
	end do
    
    do i = 1, cells, 1
       do j = 1, cells, 1
       	Sigmaq(i,j)=mProb(i)*mProb(j)
       end do
    end do
    
    Sigma=Sigmap-Sigmaq
	return	
end subroutine Sigmafun

subroutine scoreSigma(Prob,mProb,Sigmascore)
  use IMSL
  use constants
  implicit none 
  real mProb(cells)
  real Prob(factor,levelthree),Score(factor,levelthree)
  real ScoreLS(6,36)
  real Sigma(cells,cells),Sigmascore(cells,cells),Intermatrix(6,6)
  integer i,j,k,t
  real unitone(levelone),unittwo(leveltwo),unitthree(levelthree)

  call Locuniscore(Prob,Score(:,:))

  call Sigmafun(mProb,Sigma)
    unitone(:)=1.0
    unittwo(:)=1.0
    unitthree(:)=1.0
    j=0
    i=i+1
    call Kron(Score(1,1:level(1)),unittwo,unitthree,ScoreLS(i,:))
    i=i+1
    call Kron(unitone,Score(2,1:level(2)),unitthree,ScoreLS(i,:))
    i=i+1
    call Kron(unitone,unittwo,Score(3,1:level(3)),ScoreLS(i,:)) 
    i=i+1
    call Kron(Score(1,1:level(1)),Score(2,1:level(2)),unitthree,ScoreLS(i,:))
    i=i+1
    call Kron(Score(1,1:level(1)),unittwo,Score(3,1:level(3)),ScoreLS(i,:))
    i=i+1
    call Kron(unitone,Score(2,1:level(2)),Score(3,1:level(3)),ScoreLS(i,:))         

 

  Intermatrix=scoreLS .x. Sigma .x. ( .t. scoreLS)

  Intermatrix= .i. Intermatrix

  Sigmascore=( .t. scoreLS) .x. Intermatrix .x. scoreLS
    return
end subroutine scoreSigma

subroutine Kron(unitone,unittwo,unitthree,scoreunit)
	use IMSL
	use constants
	implicit none 

	integer i,j,k,t
	real unitone(levelone),unittwo(leveltwo),unitthree(levelthree)
	real scoreunit(cells)
	t=0
	do i = 1, level(1), 1
		do j = 1, level(2), 1
			do k = 1, level(3), 1
				t=t+1
				scoreunit(t)=unitone(i)*unittwo(j)*unitthree(k)
			end do			
		end do		
	end do
	return
	
end subroutine Kron

subroutine MOCtest(mProb,sample_ewma,Sigmascore,Test)
  	use IMSL
	use constants
	implicit none
    integer i,j
    real mProb(cells),sample_ewma(cells)
    real Probm(1,cells),sample_ewmam(1,cells)
    real var(1,cells),tanvesvar(cells,1)
	  real Sigmascore(cells,cells)
	  real Testm(1,1),Test 

    Probm(1,:)=mProb
    sample_ewmam(1,:)=sample_ewma
    var=sample_ewmam-sample_size*Probm
    tanvesvar= .t. var
    Testm=var .x. Sigmascore .x. tanvesvar
    Test=Testm(1,1)/(1.0*sample_size)
	return

	return

		
end subroutine MOCtest


subroutine OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)

use IMSL
use constants
implicit none

    real simuarl(simu)
    real Limit
    integer i,j,k,erre
    real Test
    real arl,st,std
    real sample(cells),sample_ewma(cells)
    real mProb(cells)
    real Sigmascore(cells,cells)
    real Prob(factor,levelthree),shift(3),Cov(3,3)


    simuarl(:)=-1.0
    erre=0
    do i=1,simu
      sample_ewma=sample_size*mProb
    Test=0.0
    do j = 1, ewma, 1 
      call msampleIC(sample)
      sample_ewma=lamda*sample+co_lamda*sample_ewma
      call MOCtest(mProb,sample_ewma,Sigmascore,Test) 
    if (Test>Limit) then
        exit
    endif
    end do
    if (Test>Limit) then
        erre=erre+1
        cycle
    endif
    simuarl(i)=0.0
  print*,i
  print*,simuarl(i)
    do while (Test<Limit)
      call msample(shift,Cov,sample)
      sample_ewma=lamda*sample+co_lamda*sample_ewma
      call MOCtest(mProb,sample_ewma,Sigmascore,Test) 
      simuarl(i)=simuarl(i)+1.0
     end do
    print*,i
    print*,simuarl(i)
    print*,"****************"
    end do
    arl=0.0   
    do j=1,simu
    if (simuarl(j)==-1.0) then
      cycle
    endif
    arl=arl+simuarl(j)
    enddo 
    arl=arl/(simu-erre)

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

