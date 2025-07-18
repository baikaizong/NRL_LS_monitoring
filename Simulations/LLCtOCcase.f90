Module constants  
  implicit none 
  integer, parameter :: sample_size=100,simu=10000,ewma=50,factor=3,cells=36 
  integer, parameter :: levelone=3,leveltwo=3,levelthree=4
  real, parameter :: lamda=0.1,co_lamda=0.9,DF=4.0
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
  real shift(3),Cov(3,3)
  real mProb(cells)
  real Sigmascore(cells,cells)
  real Prob(factor,levelthree)


    do i = 1, factor, 1
    do j = 1, factor, 1
       Cov(i,j)=0.5**(ABS(i-j))
    end do    
  end do 
  open(11,file="F:/Fortran_locationscale/TXT/Probt.txt")
  do i = 1, cells, 1
  read(11,*),mProb(i)   
  end do
  close(11)
  Limit=2.331250 
  call tDFmatrix(Prob)
  call scoreSigma(Prob,mProb,Sigmascore)

  open(10,file="F:/Fortran_locationscale/TXT/LLCtOCcase_924-33.txt")
  call cpu_time(dClkStart)
  shift(:)=1.0
  write(10,*),"limit is",limit
   shift(3)=1.02
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"shift is ",shift
  write(10,*),"arl is",arl
  write(10,*),"std is",std

   shift(3)=1.05
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"shift is ",shift
  write(10,*),"arl is",arl
  write(10,*),"std is",std

  shift(3)=1.10
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"shift is ",shift
  write(10,*),"arl is",arl
  write(10,*),"std is",std

  shift(3)=1.20
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"shift is ",shift
  write(10,*),"arl is",arl
  write(10,*),"std is",std

 

  shift(3)=0.98
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"shift is ",shift
  write(10,*),"arl is",arl
  write(10,*),"std is",std

   shift(3)=0.95
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"shift is ",shift
  write(10,*),"arl is",arl
  write(10,*),"std is",std

  shift(3)=0.90
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"shift is ",shift
  write(10,*),"arl is",arl
  write(10,*),"std is",std

  shift(3)=0.80
  call  OCARL(Limit,shift,Cov,mProb,Prob,Sigmascore,arl,std)
  write(10,*),"shift is ",shift
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

  real mnsample(sample_size,factor),mnsamplenormal(sample_size,factor),mnsamplechi(sample_size),sample(cells)
  real Cov(factor,factor),RSIG(factor,factor),threshhold(3,4) 
  integer irank,iseed
  integer cell(sample_size,factor)
  integer i,j,k,u,v,x,w

   threshhold(1,1)=-0.2
   threshhold(1,2)=1.6
   threshhold(1,3)=1000000.0
   threshhold(2,1)=-0.5
   threshhold(2,2)=0.8
   threshhold(2,3)=1000000.0
   threshhold(3,1)=-0.7
   threshhold(3,2)=1.0
   threshhold(3,3)=1.8
   threshhold(3,4)=1000000.0

    do i = 1, factor, 1
    do j = 1, factor, 1
       Cov(i,j)=0.5**(ABS(i-j))
    end do    
  end do 
  
  call CHFAC (factor,Cov,factor,TOL,irank,RSIG,factor)
  call RNGET (iseed)   
  call RNMVN(sample_size,factor,RSIG,factor,mnsamplenormal,sample_size)
  call RNGET (iseed)
  call RNCHI (sample_size, DF, mnsamplechi)
  do i = 1, sample_size, 1
    mnsample(i,:)=mnsamplenormal(i,:)/sqrt(mnsamplechi(i)/DF)
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

  real mnsample(sample_size,factor),mnsamplenormal(sample_size,factor),mnsamplechi(sample_size),sample(cells)
  real Cov(factor,factor),RSIG(factor,factor),threshhold(3,4) 
  integer irank,iseed
  integer cell(sample_size,factor)
  integer i,j,k,u,v,x,w
  real shift(factor)

   threshhold(1,1)=-0.2
   threshhold(1,2)=1.6
   threshhold(1,3)=1000000.0
   threshhold(2,1)=-0.5
   threshhold(2,2)=0.8
   threshhold(2,3)=1000000.0
   threshhold(3,1)=-0.7
   threshhold(3,2)=1.0
   threshhold(3,3)=1.8
   threshhold(3,4)=1000000.0

  
  call CHFAC (factor,Cov,factor,TOL,irank,RSIG,factor)
  call RNGET (iseed)   
  call RNMVN(sample_size,factor,RSIG,factor,mnsamplenormal,sample_size)
  call RNGET (iseed)
  call RNCHI (sample_size, DF, mnsamplechi)
  do i = 1, sample_size, 1
    mnsample(i,:)=(mnsamplenormal(i,:)/sqrt(mnsamplechi(i)/DF))*shift
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

subroutine tDFmatrix(SProb)
	use IMSL
	use constants
	implicit none
	real SProb(factor,levelthree),threshhold(3,4) 
	real interval(factor,6)
	integer i,j

   threshhold(1,1)=-0.2
   threshhold(1,2)=1.6
   threshhold(1,3)=1000000.0
   threshhold(2,1)=-0.5
   threshhold(2,2)=0.8
   threshhold(2,3)=1000000.0
   threshhold(3,1)=-0.7
   threshhold(3,2)=1.0
   threshhold(3,3)=1.8
   threshhold(3,4)=1000000.0
     
    do i = 1,factor, 1
    	j=1
    	interval(i,j)=0.0
    	do j = 2, level(i), 1
    		interval(i,j)=TDF(threshhold(i,j-1),DF)
    	end do
    	interval(i,j)=1.0  	
    end do  
    do i = 1, factor, 1
       do j = 1, level(i), 1
       	  Sprob(i,j)=interval(i,j+1)-interval(i,j)
       end do
    end do

	return
	
end subroutine tDFmatrix


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
	real Prob(factor,levelthree),Score(2,factor,levelthree)
	real ScoreLS(23,36)
	real Sigma(cells,cells),Sigmascore(cells,cells),Intermatrix(23,23)
	integer i,j,k,t
	real unitthree(3),unitfour(4),Jdesignthree(2,3),Jdesignfour(3,4)


    unitthree(:)=1.0
    unitfour(:)=1.0
    data Jdesignthree /1,0,0,1,-1,-1/
    data Jdesignfour /0,0,1,0,1,0,1,0,0,-1,-1,-1/
	call Sigmafun(mProb,Sigma)
    i=0
    i=i+1
    call Kron(Jdesignthree(1,:),unitthree,unitfour,ScoreLS(i,:))
    i=i+1
    call Kron(Jdesignthree(2,:),unitthree,unitfour,ScoreLS(i,:))
    i=i+1
    call Kron(unitthree,Jdesignthree(1,:),unitfour,ScoreLS(i,:))
    i=i+1
    call Kron(unitthree,Jdesignthree(2,:),unitfour,ScoreLS(i,:))
    i=i+1
    call Kron(unitthree,unitthree,Jdesignfour(1,:),ScoreLS(i,:))
    i=i+1
    call Kron(unitthree,unitthree,Jdesignfour(2,:),ScoreLS(i,:)) 
    i=i+1
    call Kron(unitthree,unitthree,Jdesignfour(3,:),ScoreLS(i,:))
    i=i+1
    call Kron(Jdesignthree(1,:),Jdesignthree(1,:),unitfour,ScoreLS(i,:))
    i=i+1
    call Kron(Jdesignthree(1,:),Jdesignthree(2,:),unitfour,ScoreLS(i,:))
    i=i+1
    call Kron(Jdesignthree(2,:),Jdesignthree(1,:),unitfour,ScoreLS(i,:))           	
    i=i+1
    call Kron(Jdesignthree(2,:),Jdesignthree(2,:),unitfour,ScoreLS(i,:))    
    i=i+1
    call Kron(Jdesignthree(1,:),unitthree,Jdesignfour(1,:),ScoreLS(i,:))     
    i=i+1
    call Kron(Jdesignthree(1,:),unitthree,Jdesignfour(2,:),ScoreLS(i,:))
    i=i+1
    call Kron(Jdesignthree(1,:),unitthree,Jdesignfour(3,:),ScoreLS(i,:))
    i=i+1
    call Kron(Jdesignthree(2,:),unitthree,Jdesignfour(1,:),ScoreLS(i,:))     
    i=i+1
    call Kron(Jdesignthree(2,:),unitthree,Jdesignfour(2,:),ScoreLS(i,:))
    i=i+1
    call Kron(Jdesignthree(2,:),unitthree,Jdesignfour(3,:),ScoreLS(i,:))
    i=i+1
    call Kron(unitthree,Jdesignthree(1,:),Jdesignfour(1,:),ScoreLS(i,:))     
    i=i+1
    call Kron(unitthree,Jdesignthree(1,:),Jdesignfour(2,:),ScoreLS(i,:))
    i=i+1
    call Kron(unitthree,Jdesignthree(1,:),Jdesignfour(3,:),ScoreLS(i,:))
    i=i+1
    call Kron(unitthree,Jdesignthree(2,:),Jdesignfour(1,:),ScoreLS(i,:))     
    i=i+1
    call Kron(unitthree,Jdesignthree(2,:),Jdesignfour(2,:),ScoreLS(i,:))
    i=i+1
    call Kron(unitthree,Jdesignthree(2,:),Jdesignfour(3,:),ScoreLS(i,:)) 


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

subroutine LLCTest(mProb,sample_ewma,Sigmascore,Test)
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

		
end subroutine LLCTest

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
      call LLCTest(mProb,sample_ewma,Sigmascore,Test) 
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
      call LLCTest(mProb,sample_ewma,Sigmascore,Test) 
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