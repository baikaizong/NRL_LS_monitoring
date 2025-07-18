Module constants  
  implicit none 
  integer, parameter :: sample_size=10,samplelevelseven=1059,samplelevelsix=2192,simu=10000,ewma=40,factor=3,cells=27
  integer, parameter :: levels=3
  real, parameter :: lamda=0.1,co_lamda=0.9
  real level(3)
  data level /3,3,3/

End module constants
program main

  use IMSL
  use constants
  implicit none

  call Onlinetest()

  stop

end program main

subroutine msample(threshhold,cell,sampleseven)
  use IMSL
  use constants
  implicit none

  real mnsamplelevelseven(samplelevelseven,factor),sampleseven(cells)
  real threshhold(levels,factor) 
  integer cell(samplelevelseven,factor)
  integer i,j,k,u,v,x,w

    open(unit=11, file="D:/Location-Scale-Ordinal/Fortran_locationscale/Case_wine/wine_IC.txt", status="old", action="read")
    do i = 1, samplelevelseven, 1
      read(11,*),mnsamplelevelseven(i,:)
    end do
    close(11)
    
    cell(:,:)=1
    do i = 1, samplelevelseven, 1
      do j = 1, factor, 1
        do k = 1, level(j), 1
          if ( mnsamplelevelseven(i,j)>threshhold(k,j) ) then
            cell(i,j)=cell(i,j)+1
        else
            exit
          end if
        end do        
      end do       
     end do 
    
     sampleseven(:)=0.0
     k=0
     do u = 1, level(1), 1
      do v = 1, level(2), 1
        do x = 1, level(3), 1
              k=k+1
            do i = 1, samplelevelseven, 1
              if ( cell(i,1)==u .and. cell(i,2)==v .and. cell(i,3)==x ) then
                sampleseven(k)=sampleseven(k)+1.0
              end if  
            end do
        end do
      end do
     end do
     return
  
end subroutine msample

subroutine sampleonline(cell,sample)
  use IMSL
  use constants
  implicit none
  integer cell(samplelevelseven,factor)
  real sample(cells),realindex(sample_size)
  integer indextem(sample_size),iseed
  integer i,j,k,u,v,x,w
    
    CALL RNUN (sample_size, realindex)
    realindex=realindex*samplelevelseven
    indextem=INT(realindex)+1

     sample(:)=0.0
     k=0
     do u = 1, level(1), 1
      do v = 1, level(2), 1
        do x = 1, level(3), 1
              k=k+1
            do i = 1, sample_size, 1
              if ( cell(indextem(i),1)==u .and. cell(indextem(i),2)==v .and. cell(indextem(i),3)==x ) then
                sample(k)=sample(k)+1.0
              end if  
            end do
        end do
      end do
     end do
     return 
end subroutine sampleonline

subroutine msampleOC(threshhold,cell)
	use IMSL
	use constants
	implicit none

	real mnsamplelevelsix(samplelevelsix,factor)
	real threshhold(levels,factor) 
	integer cell(samplelevelsix,factor)
	integer i,j,k,u,v,x,w

    open(unit=11, file="D:/Location-Scale-Ordinal/Fortran_locationscale/Case_wine/wine_OC.txt", status="old", action="read")
    do i = 1, samplelevelsix, 1
    	read(11,*),mnsamplelevelsix(i,:)
    end do
    close(11)
    
    cell(:,:)=1
    do i = 1, samplelevelsix, 1
    	do j = 1, factor, 1
    		do k = 1, level(j), 1
    			if ( mnsamplelevelsix(i,j)>threshhold(k,j) ) then
    				cell(i,j)=cell(i,j)+1
				else
				    exit
    			end if
    		end do    		
    	end do       
     end do 
  
     return
	
end subroutine msampleOC

subroutine sampleonlineOC(tem,cell,sample)
	use IMSL
	use constants
	implicit none
	integer cell(samplelevelsix,factor)
	real sample(cells),realindex(sample_size)
	integer indextem(sample_size),iseed
	integer i,j,k,u,v,x,w,tem
    
  sample(:)=0.0
     k=0
     do u = 1, level(1), 1
     	do v = 1, level(2), 1
     		do x = 1, level(3), 1
			        k=k+1
     				do i = tem, tem+sample_size-1, 1
     					if ( cell(i,1)==u .and. cell(i,2)==v .and. cell(i,3)==x ) then
     						sample(k)=sample(k)+1.0
     					end if	
     				end do
     		end do
     	end do
     end do
     return	
end subroutine sampleonlineOC

subroutine Locuniscore(Prob,Score)
	use IMSL
	use constants
	implicit none 
	real Prob(factor,levels),score(factor,levels),Sumscore(factor,levels)
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
	real Prob(factor,levels),score(factor,levels),Sumscore(factor,levels)
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
	real Prob(factor,levels),Score(2,factor,levels)
	real ScoreLS(9,cells)
	real Sigma(cells,cells),Sigmascore(cells,cells),Intermatrix(9,9)
	integer i,j,k,t
	real unitone(levels),unittwo(levels),unitthree(levels)

	call Locuniscore(Prob,Score(1,:,:))
	call Scauniscore(Prob,score(2,:,:))
	call Sigmafun(mProb,Sigma)
    unitone(:)=1.0
    unittwo(:)=1.0
    unitthree(:)=1.0
    j=0
    i=i+1
    call Kron(Score(1,1,1:level(1)),unittwo,unitthree,ScoreLS(i,:))
    i=i+1
    call Kron(unitone,Score(1,2,1:level(2)),unitthree,ScoreLS(i,:))
    i=i+1
    call Kron(unitone,unittwo,Score(1,3,1:level(3)),ScoreLS(i,:))
    i=i+1
    call Kron(Score(2,1,1:level(1)),unittwo,unitthree,ScoreLS(i,:))
    i=i+1
    call Kron(unitone,Score(2,2,1:level(2)),unitthree,ScoreLS(i,:))
    i=i+1
    call Kron(unitone,unittwo,Score(2,3,1:level(3)),ScoreLS(i,:))  
    i=i+1
    call Kron(Score(1,1,1:level(1)),Score(1,2,1:level(2)),unitthree,ScoreLS(i,:))
    i=i+1
    call Kron(Score(1,1,1:level(1)),unittwo,Score(1,3,1:level(3)),ScoreLS(i,:))
    i=i+1
    call Kron(unitone,Score(1,2,1:level(2)),Score(1,3,1:level(3)),ScoreLS(i,:))       	

 
    
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
	real unitone(levels),unittwo(levels),unitthree(levels)
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

subroutine MLSOTest(mProb,sample_ewma,Sigmascore,Test)
  	use IMSL
	use constants
	implicit none
    integer i,j
    real mProb(cells),sample_ewma(cells)
    real Probm(1,cells),sample_ewmam(1,cells)
    real var(1,cells),tanvesvar(cells,1)
	  real Sigmascore(cells,cells)
	  real*8 Testm(1,1),Test 

    Probm(1,:)=mProb
    sample_ewmam(1,:)=sample_ewma
    var=sample_ewmam-sample_size*Probm
    tanvesvar= .t. var
    Testm=var .x. Sigmascore .x. tanvesvar
    Test=Testm(1,1)/(1.0*sample_size)
	return
		
end subroutine MLSOTest

subroutine Onlinetest()

  use IMSL
  use constants
  implicit none

    real simuarl(simu)
    integer i,j,k,erre,iseed,tem
    real*8 Test
	 real threshhold(levels,factor)
    real sample(cells),sampleIC(cells),sample_ewma(cells)
    real mProb(cells)
    real Sigmascore(cells,cells)
    real Prob(factor,levels)
	 integer cell(samplelevelseven,factor),cellOC(samplelevelsix,factor)

   Prob(:,1)=0.25
   Prob(:,2)=0.5
   Prob(:,3)=0.25

   open(unit=12, file="D:/Location-Scale-Ordinal/Fortran_locationscale/Case_wine/thresholds20200228.txt", status="old", action="read")
	
   do i = 1, levels, 1
   	  read(12,*),threshhold(i,:)
   end do
   close(12)
   call msample(threshhold,cell,sampleIC)
   mProb=1.0*sampleIC/samplelevelseven
   call scoreSigma(Prob,mProb,Sigmascore)
   call msampleOC(threshhold,cellOC)
    sample_ewma=sample_size*mProb
    Test=0.0    
	open(unit=13, file="D:/Location-Scale-Ordinal/Fortran_locationscale/Case_wine/Onlinetest20200228.txt")
    do j = 1, ewma, 1 
      call sampleonline(cell,sample)
      sample_ewma=lamda*sample+co_lamda*sample_ewma
      call MLSOTest(mProb,sample_ewma,Sigmascore,Test)
	  print*,Test
	  write(13,*),Test
    end do 
    tem=1.0

    
    do i=1,10
      call sampleonlineOC(tem,cellOC,sample)
	  !print*,sum(sample)
      sample_ewma=lamda*sample+co_lamda*sample_ewma
      call MLSOTest(mProb,sample_ewma,Sigmascore,Test)
	  print*,Test
      write(13,*),Test
      tem=tem+40
	  !print*,Test 
    end do
  return 
end subroutine Onlinetest

