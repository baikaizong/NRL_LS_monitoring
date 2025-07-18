Module constants  
  implicit none 
  integer, parameter :: sample_size=100,simu=10000,ewma=50,factor=3,cells=36 ,simu_IC=100000
  integer, parameter :: levelone=3,leveltwo=3,levelthree=4
  !real, parameter :: lamda=0.1,co_lamda=0.9
  real, parameter :: TOL=0.00001
  real,parameter::CovIC(factor,factor)=reshape((/1.0,0.5,0.25,0.5,1.0,0.5,0.25,0.5,1.0/),(/factor,factor/))
  real level(3)
  data level /3,3,4/

End module constants




subroutine msample(shiftlocation,shiftscale,Cov,sample)
  use IMSL
  use constants
  implicit none

  real mnsample(sample_size,factor),sample(cells)
  real Cov(factor,factor),RSIG(factor,factor),threshhold(3,4) 
  integer irank,iseed
  integer cell(sample_size,factor)
  integer i,j,k,u,v,x,w
  real shiftlocation(factor),shiftscale(factor)

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
  call RNMVN(sample_size,factor,RSIG,factor,mnsample,sample_size)
  do j = 1, factor, 1
    mnsample(:,j)=(mnsample(:,j)+shiftlocation(j))*(1.0+shiftscale(j))
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

subroutine CDFmatrix(SProb)
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
   Sprob=0.0   
    do i = 1,factor, 1
      j=1
      interval(i,j)=0.0
      do j = 2, level(i), 1
        interval(i,j)=ANORDF(threshhold(i,j-1))
      end do
      interval(i,j)=1.0   
    end do  
    do i = 1, factor, 1
       do j = 1, level(i), 1
          Sprob(i,j)=interval(i,j+1)-interval(i,j)
       end do
    end do
  return
  
end subroutine CDFmatrix


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
  real Prob(factor,levelthree),Score(2,factor,levelthree)
  real ScoreLS(9,36),ScoreLStem(1,36)
  real Sigma(cells,cells),Sigmascore(9,cells,cells),Intermatrix(1,1)
  integer i,j,k,t
  real unitone(levelone),unittwo(leveltwo),unitthree(levelthree)

  call Locuniscore(Prob,Score(1,:,:))
  call Scauniscore(Prob,score(2,:,:))

  call Sigmafun(mProb,Sigma)
    unitone(:)=1.0
    unittwo(:)=1.0
    unitthree(:)=1.0
    i=0
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
  do i = 1, 9, 1
    ScoreLStem(1,:)=scoreLS(i,:)
    Intermatrix=ScoreLStem .x. Sigma .x. ( .t. ScoreLStem)
    Intermatrix= .i. Intermatrix
    Sigmascore(i,:,:)=( .t. ScoreLStem) .x. Intermatrix .x. ScoreLStem
  end do
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

subroutine MLSOTest(mProb,sample_ewma,Sigmascore,Testm)
    use IMSL
  use constants
  implicit none
    integer i,j
    real mProb(cells),sample_ewma(cells)
    real Probm(1,cells),sample_ewmam(1,cells)
    real var(1,cells),tanvesvar(cells,1)
    real Sigmascore(9,cells,cells)
    real Testm(1,9)

    Probm(1,:)=mProb
    sample_ewmam(1,:)=sample_ewma
    var=sample_ewmam-sample_size*Probm
    tanvesvar= .t. var
    do i = 1, 9, 1
      Testm(:,i:i)=var .x. Sigmascore(i,:,:) .x. tanvesvar
	end do
  return  
end subroutine MLSOTest

subroutine ExpSig(mProb,lamda,Exp_LRT,Sig_LRT)

  use IMSL
  use constants
  implicit none

    real simuarl(simu),zerovector(factor),Exp_LRT(9),Sig_LRT(9)
    real Limit,LimitL,LimitR
    integer i,j,k,iseed
    real Testm(1,9)
    real sample(cells),sample_ewma(cells)
    real mProb(cells)
    real Sigmascore(9,cells,cells)
    real Prob(factor,levelthree)
    real Cov(factor,factor)
    real lamda,co_lamda
    Exp_LRT=0.0
    Sig_LRT=0.0
    co_lamda=1.0-lamda
    zerovector=0.0
   call CDFmatrix(Prob)
   call scoreSigma(Prob,mProb,Sigmascore)

    sample_ewma=sample_size*mProb
    do j = 1, ewma, 1 
      call msample(zerovector,zerovector,CovIC,sample)
      sample_ewma=lamda*sample+co_lamda*sample_ewma
    end do
    do j = 1, simu_IC, 1
        print*,j
        print*,"----------------" 
      call msample(zerovector,zerovector,CovIC,sample)
      sample_ewma=lamda*sample+co_lamda*sample_ewma
      call MLSOTest(mProb,sample_ewma,Sigmascore,Testm) 
      do i = 1, 9, 1
        Exp_LRT(i)=Exp_LRT(i)+Testm(1,i)
      end do     
    end do
    Exp_LRT=Exp_LRT/simu_IC
    do j = 1, ewma, 1 
      call msample(zerovector,zerovector,CovIC,sample)
      sample_ewma=lamda*sample+co_lamda*sample_ewma
    end do
    do j = 1, simu_IC, 1
        print*,j
        print*,"----------------" 
      call msample(zerovector,zerovector,CovIC,sample)
      sample_ewma=lamda*sample+co_lamda*sample_ewma
      call MLSOTest(mProb,sample_ewma,Sigmascore,Testm) 
      do i = 1, 9, 1
        Sig_LRT(i)=Sig_LRT(i)+((Exp_LRT(i)-Testm(1,i))**2.0)
      end do     
    end do
    do i = 1, 9, 1
      Sig_LRT(i)=sqrt(Sig_LRT(i)/simu_IC)
    end do    
  return 
end subroutine ExpSig

subroutine Limitsearch(limitl,limitR,Limit,arl,std,mProb,lamda,Exp_LRT,Sig_LRT)

  use IMSL
  use constants
  implicit none

    real simuarl(simu),zerovector(factor),Exp_LRT(9),Sig_LRT(9)
    real Limit,LimitL,LimitR
    integer i,j,k,iseed
    real Test,Testm(1,9)
    real arl,st,std
    real sample(cells),sample_ewma(cells)
    real mProb(cells)
    real Sigmascore(9,cells,cells)
    real Prob(factor,levelthree)
    real Cov(factor,factor)
    real lamda,co_lamda
    arl=0.0
    std=0.0
    co_lamda=1.0-lamda
    zerovector=0.0
   call CDFmatrix(Prob)
   call scoreSigma(Prob,mProb,Sigmascore)

  do while (abs(arl-370.0) >= 4)
  Limit = (LimitL + LimitR) / 2.0
    simuarl(:)=-1.0
    i=1
  do while (i<=simu)
      sample_ewma=sample_size*mProb
    Test=0.0
    do j = 1, ewma, 1 
      call msample(zerovector,zerovector,CovIC,sample)
      sample_ewma=lamda*sample+co_lamda*sample_ewma
      call MLSOTest(mProb,sample_ewma,Sigmascore,Testm)
      Test=0
      do k = 1, 9, 1
        Testm(1,k)=(Testm(1,k)-Exp_LRT(k))/Sig_LRT(k)
         if ( Test<Testm(1,k) ) then
           Test=Testm(1,k) 
         end if
       end do 
    if (Test>Limit) then
        exit
    endif
    end do
    if (Test>Limit) then
        cycle
    endif
    simuarl(i)=0.0
    do while (Test<Limit)
      call msample(zerovector,zerovector,CovIC,sample)
      sample_ewma=lamda*sample+co_lamda*sample_ewma
      call MLSOTest(mProb,sample_ewma,Sigmascore,Testm)  
      do k = 1, 9, 1
        Testm(1,k)=(Testm(1,k)-Exp_LRT(k))/Sig_LRT(k)
         if ( Test<Testm(1,k) ) then
           Test=Testm(1,k) 
         end if
       end do  
    simuarl(i)=simuarl(i)+1.0
     end do
    print*,i
    print*,simuarl(i)
    print*,"****************" 
    i=i+1
    end do
    arl=(sum(simuarl))/(simu)
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
    st=st+(simuarl(j)-arl)**2.0
  enddo 
  std=sqrt(st/(simu-1.0))
  std=std/(sqrt(dble(simu)-1)) 
  return 
end subroutine Limitsearch

subroutine OCARL(Limit,shiftlocation,shiftscale,Cov,mProb,Prob,Sigmascore,arl,std,lamda,Exp_LRT,Sig_LRT)

use IMSL
use constants
implicit none

    real simuarl(simu),Exp_LRT(9),Sig_LRT(9)
    real Limit
    integer i,j,k
    real Testm(1,9),Test
    real arl,st,std
    real sample(cells),sample_ewma(cells)
    real mProb(cells)
    real Sigmascore(9,cells,cells)
    real Prob(factor,levelthree),shift(3),Cov(3,3)
    real zerovector(factor),shiftlocation(factor),shiftscale(factor)
    real lamda,co_lamda
    co_lamda=1.0-lamda
    zerovector=0.0
    simuarl(:)=-1.0
    arl=0.0
    std=0.0
    i=1
    do while(i<=simu)
    sample_ewma=sample_size*mProb
    Test=0.0
    do j = 1, ewma, 1 
      call msample(zerovector,zerovector,CovIC,sample)
      sample_ewma=lamda*sample+co_lamda*sample_ewma
      call MLSOTest(mProb,sample_ewma,Sigmascore,Testm)
      do k = 1, 9, 1
        Testm(1,k)=(Testm(1,k)-Exp_LRT(k))/Sig_LRT(k)
         if ( Test<Testm(1,k) ) then
           Test=Testm(1,k) 
         end if
       end do  
    if (Test>Limit) then
        exit
    endif
    end do
    if (Test>Limit) then
        cycle
    endif
    
    simuarl(i)=0.0
    do while (Test<Limit)
      call msample(shiftlocation,shiftscale,Cov,sample)
      sample_ewma=lamda*sample+co_lamda*sample_ewma
      call MLSOTest(mProb,sample_ewma,Sigmascore,Testm) 
      do k = 1, 9, 1
        Testm(1,k)=(Testm(1,k)-Exp_LRT(k))/Sig_LRT(k)
         if ( Test<Testm(1,k) ) then
           Test=Testm(1,k) 
         end if
       end do 
      simuarl(i)=simuarl(i)+1.0
     end do
    print*,i
    print*,simuarl(i)
    print*,"****************"
    i=i+1
    end do
    arl=0.0   
    do j=1,simu
    if (simuarl(j)==-1.0) then
      cycle
    endif
    arl=arl+simuarl(j)
    enddo 
    arl=arl/(simu)

   st=0.0                                                                                                                                                                                                                                                                                                                                                                            
   do j=1,simu
    st=st+(simuarl(j)-arl)**2.0
  enddo 
  std=sqrt(st/(simu-1.0))
  std=std/(sqrt(dble(simu)-1)) 
  return 
end subroutine OCARL

program main

  use IMSL
  use constants
  implicit none

  real Limit,arl,std,limitL,limitR,Exp_LRT(9),Sig_LRT(9)
  real dClkStart, dClkFinish
  integer i,j
  real sample(cells)
  real shiftlocation(factor),shiftscale(factor),Cov(factor,factor)
  real mProb(cells)
  real Sigmascore(9,cells,cells)
  real Prob(factor,levelthree)
  real shift_value(10)
  real lamda
  call cpu_time(dClkStart)
  Cov=CovIC  
  open(10,file="./Di_TXT/Di_MLSOnormal-230210.txt")
  open(12,file="./Di_TXT/Di_MLSOnormal-230210-A.txt")
  open(11,file="./MOCP/Probnormal_MLSO.txt")
  do i = 1, cells, 1
  read(11,*),mProb(i)   
  end do
  close(11)   
  lamda=0.10
  write(10,*),"lamda is",lamda
  write(10,*),"*****************************"
  !Limit=1.22875 
  call CDFmatrix(Prob)
  call scoreSigma(Prob,mProb,Sigmascore)
  call ExpSig(mProb,lamda,Exp_LRT,Sig_LRT)
  LimitL=5.0
  limitR=7.0
  !limit=0.621875
  call Limitsearch(limitl,limitR,Limit,arl,std,mProb,lamda,Exp_LRT,Sig_LRT)
      write(10,*),"lambda is",lamda
      write(10,*),"Limit is",limit
      write(10,*),"arl is",arl
      write(10,*),"std is",std 
  shift_value=(/-0.20,-0.10,-0.04,0.04,0.10,0.20,0.0,0.0,0.0,0.0/) 
  shiftscale=0.0
  shiftlocation=0.0
  do i = 1, 6, 1
      shiftlocation(1)=shift_value(i)
      call OCARL(Limit,shiftlocation,shiftscale,Cov,mProb,Prob,Sigmascore,arl,std,lamda,Exp_LRT,Sig_LRT)
      write(10,*),"shiftlocation(1) is",shiftlocation(1)
      write(10,*),"shiftscale(1) is",shiftscale(1)
      write(10,*),"arl is",arl
      write(10,*),"std is",std 
      if (arl>=100 ) then
        write(12,"(1X,I3,A1,1X,A1,F4.2,A1)"),int(arl),"&","(",std,")"
        else if (arl>=10 .and. arl<100) then
        write(12,"(1X,F4.1,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")" 
        else 
        write(12,"(1X,F4.2,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")"  
      end if
  end do
  shiftscale=0.0
  shiftlocation=0.0
  do i = 1, 6, 1
      shiftscale(1)=shift_value(i)
      call OCARL(Limit,shiftlocation,shiftscale,Cov,mProb,Prob,Sigmascore,arl,std,lamda,Exp_LRT,Sig_LRT)
      write(10,*),"shiftlocation(1) is",shiftlocation(1)
      write(10,*),"shiftscale(1) is",shiftscale(1)
      write(10,*),"arl is",arl
      write(10,*),"std is",std 
      if (arl>=100 ) then
        write(12,"(1X,I3,A1,1X,A1,F4.2,A1)"),int(arl),"&","(",std,")"
        else if (arl>=10 .and. arl<100) then
        write(12,"(1X,F4.1,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")" 
        else 
        write(12,"(1X,F4.2,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")"  
      end if
  end do

  shift_value=(/-0.10,-0.05,0.05,0.10,0.0,0.0,0.0,0.0,0.0,0.0/)
  shiftscale=0.0
  shiftlocation=0.0
  shiftscale(1)=0.02
  do i = 1, 4, 1
      shiftlocation(1)=shift_value(i)
      call OCARL(Limit,shiftlocation,shiftscale,Cov,mProb,Prob,Sigmascore,arl,std,lamda,Exp_LRT,Sig_LRT)
      write(10,*),"shiftlocation(1) is",shiftlocation(1)
      write(10,*),"shiftscale(1) is",shiftscale(1)
      write(10,*),"arl is",arl
      write(10,*),"std is",std 
      if (arl>=100 ) then
        write(12,"(1X,I3,A1,1X,A1,F4.2,A1)"),int(arl),"&","(",std,")"
        else if (arl>=10 .and. arl<100) then
        write(12,"(1X,F4.1,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")" 
        else 
        write(12,"(1X,F4.2,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")"  
      end if
  end do
  shiftscale=0.0
  shiftlocation=0.0
  shiftlocation(1)=0.02
  do i = 1, 4, 1
      shiftscale(1)=shift_value(i)
      call OCARL(Limit,shiftlocation,shiftscale,Cov,mProb,Prob,Sigmascore,arl,std,lamda,Exp_LRT,Sig_LRT)
      write(10,*),"shiftlocation(1) is",shiftlocation(1)
      write(10,*),"shiftscale(1) is",shiftscale(1)
      write(10,*),"arl is",arl
      write(10,*),"std is",std 
      if (arl>=100 ) then
        write(12,"(1X,I3,A1,1X,A1,F4.2,A1)"),int(arl),"&","(",std,")"
        else if (arl>=10 .and. arl<100) then
        write(12,"(1X,F4.1,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")" 
        else 
        write(12,"(1X,F4.2,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")"  
      end if
  end do

  shift_value=(/-0.06,-0.03,0.03,0.06,0.0,0.0,0.0,0.0,0.0,0.0/) 
  shiftscale=0.0
  shiftlocation=0.0
  shiftscale(1)=0.05
  do i = 1, 4, 1
      shiftlocation(1)=shift_value(i)
      call OCARL(Limit,shiftlocation,shiftscale,Cov,mProb,Prob,Sigmascore,arl,std,lamda,Exp_LRT,Sig_LRT)
      write(10,*),"shiftlocation(1) is",shiftlocation(1)
      write(10,*),"shiftscale(1) is",shiftscale(1)
      write(10,*),"arl is",arl
      write(10,*),"std is",std 
      if (arl>=100 ) then
        write(12,"(1X,I3,A1,1X,A1,F4.2,A1)"),int(arl),"&","(",std,")"
        else if (arl>=10 .and. arl<100) then
        write(12,"(1X,F4.1,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")" 
        else 
        write(12,"(1X,F4.2,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")"  
      end if
  end do
  shiftscale=0.0
  shiftlocation=0.0
  shiftlocation(1)=0.05
  do i = 1, 4, 1
      shiftscale(1)=shift_value(i)
      call OCARL(Limit,shiftlocation,shiftscale,Cov,mProb,Prob,Sigmascore,arl,std,lamda,Exp_LRT,Sig_LRT)
      write(10,*),"shiftlocation(1) is",shiftlocation(1)
      write(10,*),"shiftscale(1) is",shiftscale(1)
      write(10,*),"arl is",arl
      write(10,*),"std is",std 
      if (arl>=100 ) then
        write(12,"(1X,I3,A1,1X,A1,F4.2,A1)"),int(arl),"&","(",std,")"
        else if (arl>=10 .and. arl<100) then
        write(12,"(1X,F4.1,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")" 
        else 
        write(12,"(1X,F4.2,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")"  
      end if
  end do


  write(10,*),"++++++++++++++++++++++++++++++++++++++++++++++++++"

  shift_value=(/-0.20,-0.10,-0.04,0.04,0.10,0.20,0.0,0.0,0.0,0.0/) 
  shiftscale=0.0
  shiftlocation=0.0
  do i = 1, 6, 1
      shiftscale(2)=shift_value(i)
      call OCARL(Limit,shiftlocation,shiftscale,Cov,mProb,Prob,Sigmascore,arl,std,lamda,Exp_LRT,Sig_LRT)
      write(10,*),"shiftscale(2) is",shiftscale(2)
      write(10,*),"shiftCov(1,2) is",0
      write(10,*),"arl is",arl
      write(10,*),"std is",std 
      if (arl>=100 ) then
        write(12,"(1X,I3,A1,1X,A1,F4.2,A1)"),int(arl),"&","(",std,")"
        else if (arl>=10 .and. arl<100) then
        write(12,"(1X,F4.1,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")" 
        else 
        write(12,"(1X,F4.2,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")"  
      end if
  end do
  shiftscale=0.0
  shiftlocation=0.0
  do i = 1, 6, 1
      Cov(1,2)=CovIC(1,2)+shift_value(i)
      call OCARL(Limit,shiftlocation,shiftscale,Cov,mProb,Prob,Sigmascore,arl,std,lamda,Exp_LRT,Sig_LRT)
      write(10,*),"shiftscale(2) is",shiftscale(2)
      write(10,*),"shiftCov(1,2) is",shift_value(i)
      write(10,*),"arl is",arl
      write(10,*),"std is",std 
      if (arl>=100 ) then
        write(12,"(1X,I3,A1,1X,A1,F4.2,A1)"),int(arl),"&","(",std,")"
        else if (arl>=10 .and. arl<100) then
        write(12,"(1X,F4.1,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")" 
        else 
        write(12,"(1X,F4.2,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")"  
      end if
  end do

  shift_value=(/-0.10,-0.05,0.05,0.10,0.0,0.0,0.0,0.0,0.0,0.0/)
  shiftscale=0.0
  shiftlocation=0.0
  Cov(1,2)=CovIC(1,2)+0.02
  do i = 1, 4, 1
      shiftscale(2)=shift_value(i)
      call OCARL(Limit,shiftlocation,shiftscale,Cov,mProb,Prob,Sigmascore,arl,std,lamda,Exp_LRT,Sig_LRT)
      write(10,*),"shiftscale(2) is",shiftlocation(2)
      write(10,*),"shiftCov(1,2) is",0.02
      write(10,*),"arl is",arl
      write(10,*),"std is",std 
      if (arl>=100 ) then
        write(12,"(1X,I3,A1,1X,A1,F4.2,A1)"),int(arl),"&","(",std,")"
        else if (arl>=10 .and. arl<100) then
        write(12,"(1X,F4.1,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")" 
        else 
        write(12,"(1X,F4.2,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")"  
      end if
  end do
  shiftscale=0.0
  shiftlocation=0.0
  shiftscale(2)=0.02
  do i = 1, 4, 1
      Cov(1,2)=CovIC(1,2)+shift_value(i)
      call OCARL(Limit,shiftlocation,shiftscale,Cov,mProb,Prob,Sigmascore,arl,std,lamda,Exp_LRT,Sig_LRT)
      write(10,*),"shiftscale(2) is",0.02
      write(10,*),"shiftCov(1,2) is",shift_value(i)
      write(10,*),"arl is",arl
      write(10,*),"std is",std 
      if (arl>=100 ) then
        write(12,"(1X,I3,A1,1X,A1,F4.2,A1)"),int(arl),"&","(",std,")"
        else if (arl>=10 .and. arl<100) then
        write(12,"(1X,F4.1,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")" 
        else 
        write(12,"(1X,F4.2,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")"  
      end if
  end do

  shift_value=(/-0.06,-0.03,0.03,0.06,0.0,0.0,0.0,0.0,0.0,0.0/)
  shiftscale=0.0
  shiftlocation=0.0
  Cov(1,2)=CovIC(1,2)+0.05
  do i = 1, 4, 1
      shiftscale(2)=shift_value(i)
      call OCARL(Limit,shiftlocation,shiftscale,Cov,mProb,Prob,Sigmascore,arl,std,lamda,Exp_LRT,Sig_LRT)
      write(10,*),"shiftscale(2) is",shiftlocation(2)
      write(10,*),"shiftCov(1,2) is",0.05
      write(10,*),"arl is",arl
      write(10,*),"std is",std 
      if (arl>=100 ) then
        write(12,"(1X,I3,A1,1X,A1,F4.2,A1)"),int(arl),"&","(",std,")"
        else if (arl>=10 .and. arl<100) then
        write(12,"(1X,F4.1,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")" 
        else 
        write(12,"(1X,F4.2,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")"  
      end if
  end do
  shiftscale=0.0
  shiftlocation=0.0
  shiftscale(2)=0.05
  do i = 1, 4, 1
      Cov(1,2)=CovIC(1,2)+shift_value(i)
      call OCARL(Limit,shiftlocation,shiftscale,Cov,mProb,Prob,Sigmascore,arl,std,lamda,Exp_LRT,Sig_LRT)
      write(10,*),"shiftscale(2) is",0.05
      write(10,*),"shiftCov(1,2) is",shift_value(i)
      write(10,*),"arl is",arl
      write(10,*),"std is",std 
      if (arl>=100 ) then
        write(12,"(1X,I3,A1,1X,A1,F4.2,A1)"),int(arl),"&","(",std,")"
        else if (arl>=10 .and. arl<100) then
        write(12,"(1X,F4.1,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")" 
        else 
        write(12,"(1X,F4.2,A1,1X,A1,F4.2,A1)"),arl,"&","(",std,")"  
      end if
  end do  
  call cpu_time(dClkFinish)
  write(10,*) "Time is", dClkFinish-dClkStart, "seconds"  
  stop

end program main
