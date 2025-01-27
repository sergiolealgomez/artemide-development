!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			harpy
!
!	An interface for artemide
!
!				A.Vladimirov (13.01.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module harpy
use TMDs
use TMDs_inKT
use TMDX_DY
use TMDX_DY_mT
use TMDX_SIDIS
use aTMDe_control
use uTMDPDF
use TMDR_model
use TMDF

!!! this flag is requared to guaranty that artemide is not started twice (it lead to the crush)
logical::started=.false.

 contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GENERAL
  subroutine Initialize(file)
    character(len=*)::file
    if(started) then
      write(*,*) 'artemide already runs'
    else
      call artemide_Initialize(file)
      started=.true.
    end if
  end subroutine Initialize
  
  
  subroutine ShowStatistics()
   call artemide_ShowStatistics()
  end  subroutine ShowStatistics
  
  !! call for parameters from the model
  subroutine SetReplica_TMDR(num)
  integer:: num
  call artemide_SetReplica_TMDR(num)
  end subroutine SetReplica_TMDR
  
  !!
  subroutine SetReplica_uTMDPDF(num)
  integer:: num
  call artemide_SetReplica_uTMDPDF(num)
  end subroutine SetReplica_uTMDPDF
  
  !!
  subroutine SetReplica_uTMDFF(num)
  integer:: num
  call artemide_SetReplica_uTMDFF(num)
  end subroutine SetReplica_uTMDFF
  
!!
  subroutine SetReplica_lpTMDPDF(num)
  integer:: num
  call artemide_SetReplica_lpTMDPDF(num)
  end subroutine SetReplica_lpTMDPDF
  
  !!
  subroutine SetReplica_SiversTMDPDF(num)
  integer:: num
  call artemide_SetReplica_SiversTMDPDF(num)
  end subroutine SetReplica_SiversTMDPDF
  
  !!!It sets the W mass and decay width
  subroutine SetEWparameters_Main(EWparm)
  real*8,intent(in)::EWparm(:)
  call artemide_SetEWparameters(M_W=EWparm(1),G_W=EWparm(2))
  end subroutine SetEWparameters_Main
  
  !!!It sets the DELTA in the W boson propagator
  subroutine SetDELTA_Main(DELTANew)
  real*8,intent(in)::DELTANew
  call artemide_SetDELTA(DELTANew)
  end subroutine SetDELTA_Main
  
  !!!It sets the mass and the width of new boson
  subroutine SetNewBoson_Main(NewBoson)
  real*8,intent(in)::NewBoson(:)
  call artemide_SetNewBoson(MNew=NewBoson(1),GammaNew=NewBoson(2),RatioNew=NewBoson(3))
  end subroutine SetNewBoson_Main
  
  
  !!!Sets the non-pertrubative parameters lambda
  subroutine SetLambda_Main(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    call artemide_SetNPparameters(lambdaIN)
  end subroutine SetLambda_Main  
  
    !!!Sets the non-pertrubative parameters lambda
  subroutine SetLambda_TMDR(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    call artemide_SetNPparameters_TMDR(lambdaIN)
  end subroutine SetLambda_TMDR
  
       !!!Sets the non-pertrubative parameters lambda
   subroutine SetLambda_uTMDPDF(lambdaIN)
     real*8,intent(in)::lambdaIN(:)
     call artemide_SetNPparameters_uTMDPDF(lambdaIN)
   end subroutine SetLambda_uTMDPDF
  
      !!!Sets the non-pertrubative parameters lambda flavour dependent
  subroutine SetLambda_uTMDPDF_f(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    real*8,dimension(:,:),allocatable::lambdaIN_f
    integer::i,j
    
    if (allocated(lambdaIN_f)) then
        deallocate(lambdaIN_f)
    end if
    
    allocate(lambdaIN_f(-5:5,1:size(lambdaIN,1)/11))
    
    do i=1,11
        do j=1,size(lambdaIN,1)/11
        lambdaIN_f(i-6,j)=lambdaIN(11*(j-1)+i)
        end do
    end do

    call artemide_SetNPparameters_uTMDPDF_f(lambdaIN_f)
  end subroutine SetLambda_uTMDPDF_f
  
  !!!Sets the non-pertrubative parameters lambda
  subroutine SetLambda_uTMDFF(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    call artemide_SetNPparameters_uTMDFF(lambdaIN)
  end subroutine SetLambda_uTMDFF
  
  !!!Sets the non-pertrubative parameters lambda
  subroutine SetLambda_lpTMDPDF(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    call artemide_SetNPparameters_lpTMDPDF(lambdaIN)
  end subroutine SetLambda_lpTMDPDF
  
  !!!Sets the non-pertrubative parameters lambda
  subroutine SetLambda_SiversTMDPDF(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    call artemide_SetNPparameters_SiversTMDPDF(lambdaIN)
  end subroutine SetLambda_SiversTMDPDF
  
  
  !!!! this routine set the variations of scales
  !!!! it is used for the estimation of errors
  subroutine SetScaleVariation(c1_in,c2_in,c3_in,c4_in)
    real*8::c1_in,c2_in,c3_in,c4_in
    
    call artemide_SetScaleVariations(c1_in,c2_in,c3_in,c4_in)
    
  end subroutine SetScaleVariation
  
  !! reset the number for PDF replica for uTMDPDF
  subroutine SetPDFreplica(rep)
    integer::rep
    call uTMDPDF_SetPDFreplica(rep)
  end subroutine SetPDFreplica
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!! DNP
  function getDNP(b,mu,f)
    real*8::getDNP
    real*8::b,mu
    integer::f
    
    getDNP=DNP(mu,b,f)
  end function getDNP
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!uTMDPDF
  !!!!!!!! upolarized TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDPDF_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDPDF_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDPDF_5_Evolved=uTMDPDF_5(x,bt,muf,zetaf,hadron)
    
  end function uTMDPDF_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDPDF_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDPDF_50_Evolved=uTMDPDF_50(x,bt,muf,zetaf,hadron)
    
  end function uTMDPDF_50_Evolved

  !!!!!!!! upolarized TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDPDF_5_Optimal(x,bt,hadron)
    real*8:: uTMDPDF_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDPDF_5_Optimal=uTMDPDF_5(x,bt,hadron)
    
  end function uTMDPDF_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_50_Optimal(x,bt,hadron)
    real*8:: uTMDPDF_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDPDF_50_Optimal=uTMDPDF_50(x,bt,hadron)
    
  end function uTMDPDF_50_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!uTMDFF
  !!!!!!!! upolarized TMDFF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDFF_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDFF_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDFF_5_Evolved=uTMDFF_5(x,bt,muf,zetaf,hadron)
    
  end function uTMDFF_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDFF_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDFF_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDFF_50_Evolved=uTMDFF_50(x,bt,muf,zetaf,hadron)
    
  end function uTMDFF_50_Evolved

    !!!!!!!! upolarized TMDFF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDFF_5_Optimal(x,bt,hadron)
    real*8:: uTMDFF_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDFF_5_Optimal=uTMDFF_5(x,bt,hadron)
    
  end function uTMDFF_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDFF_50_Optimal(x,bt,hadron)
    real*8:: uTMDFF_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDFF_50_Optimal=uTMDFF_50(x,bt,hadron)
    
  end function uTMDFF_50_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SiversTMDPDF
  !!!!!!!! Sivers TMDFF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function SiversTMDPDF_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: SiversTMDPDF_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  SiversTMDPDF_5_Evolved=SiversTMDPDF_5(x,bt,muf,zetaf,hadron)
    
  end function SiversTMDPDF_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function SiversTMDPDF_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: SiversTMDPDF_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  SiversTMDPDF_50_Evolved=SiversTMDPDF_50(x,bt,muf,zetaf,hadron)
    
  end function SiversTMDPDF_50_Evolved

    !!!!!!!! Sivers TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function SiversTMDPDF_5_Optimal(x,bt,hadron)
    real*8:: SiversTMDPDF_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  SiversTMDPDF_5_Optimal=SiversTMDPDF_5(x,bt,hadron)
    
  end function SiversTMDPDF_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function SiversTMDPDF_50_Optimal(x,bt,hadron)
    real*8:: SiversTMDPDF_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  SiversTMDPDF_50_Optimal=SiversTMDPDF_50(x,bt,hadron)
    
  end function SiversTMDPDF_50_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!lpTMDPDF
  !!!!!!!! linearly polarized TMDFF
 
  ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function lpTMDPDF_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: lpTMDPDF_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  lpTMDPDF_50_Evolved=lpTMDPDF_50(x,bt,muf,zetaf,hadron)
    
  end function lpTMDPDF_50_Evolved
  
  ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function lpTMDPDF_50_Optimal(x,bt,hadron)
    real*8:: lpTMDPDF_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  lpTMDPDF_50_Optimal=lpTMDPDF_50(x,bt,hadron)
    
  end function lpTMDPDF_50_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TMDs IN KT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!! upolarized TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDPDF_kT_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDPDF_kT_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDPDF_kT_5_Evolved=uTMDPDF_kT_5(x,bt,muf,zetaf,hadron)
    
  end function uTMDPDF_kT_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_kT_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDPDF_kT_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDPDF_kT_50_Evolved=uTMDPDF_kT_5(x,bt,muf,zetaf,hadron)
    
  end function uTMDPDF_kT_50_Evolved

  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDPDF_kT_5_Optimal(x,bt,hadron)
    real*8:: uTMDPDF_kT_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDPDF_kT_5_Optimal=uTMDPDF_kT_5(x,bt,hadron)
    
  end function uTMDPDF_kT_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_kT_50_Optimal(x,bt,hadron)
    real*8:: uTMDPDF_kT_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDPDF_kT_50_Optimal=uTMDPDF_kT_5(x,bt,hadron)
    
  end function uTMDPDF_kT_50_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!! upolarized TMDFF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDFF_kT_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDFF_kT_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDFF_kT_5_Evolved=uTMDFF_kT_5(x,bt,muf,zetaf,hadron)
    
  end function uTMDFF_kT_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDFF_kT_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDFF_kT_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDFF_kT_50_Evolved=uTMDFF_kT_5(x,bt,muf,zetaf,hadron)
    
  end function uTMDFF_kT_50_Evolved

  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDFF_kT_5_Optimal(x,bt,hadron)
    real*8:: uTMDFF_kT_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDFF_kT_5_Optimal=uTMDFF_kT_5(x,bt,hadron)
    
  end function uTMDFF_kT_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDFF_kT_50_Optimal(x,bt,hadron)
    real*8:: uTMDFF_kT_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDFF_kT_50_Optimal=uTMDFF_kT_5(x,bt,hadron)
    
  end function uTMDFF_kT_50_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!! linearly polarized gluon TMDPDF

  ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function lpTMDPDF_kT_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: lpTMDPDF_kT_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  lpTMDPDF_kT_50_Evolved=lpTMDPDF_kT_50(x,bt,muf,zetaf,hadron)
    
  end function lpTMDPDF_kT_50_Evolved
  
  ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function lpTMDPDF_kT_50_Optimal(x,bt,hadron)
    real*8:: lpTMDPDF_kT_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  lpTMDPDF_kT_50_Optimal=lpTMDPDF_kT_50(x,bt,hadron)
    
  end function lpTMDPDF_kT_50_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! Sivers TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function SiversTMDPDF_kT_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: SiversTMDPDF_kT_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  SiversTMDPDF_kT_5_Evolved=SiversTMDPDF_kT_5(x,bt,muf,zetaf,hadron)
    
  end function SiversTMDPDF_kT_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function SiversTMDPDF_kT_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: SiversTMDPDF_kT_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  SiversTMDPDF_kT_50_Evolved=SiversTMDPDF_kT_5(x,bt,muf,zetaf,hadron)
    
  end function SiversTMDPDF_kT_50_Evolved

  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function SiversTMDPDF_kT_5_Optimal(x,bt,hadron)
    real*8:: SiversTMDPDF_kT_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  SiversTMDPDF_kT_5_Optimal=SiversTMDPDF_kT_5(x,bt,hadron)
    
  end function SiversTMDPDF_kT_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function SiversTMDPDF_kT_50_Optimal(x,bt,hadron)
    real*8:: SiversTMDPDF_kT_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  SiversTMDPDF_kT_50_Optimal=SiversTMDPDF_kT_5(x,bt,hadron)
    
  end function SiversTMDPDF_kT_50_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DY CROSS-SECTION
  
    function DY_xSec_SingleN(process,s,qT,Q,y,includeCuts,CutParameters,Num)
    integer,intent(in),dimension(1:3)::process		!the number of process
    real*8,intent(in)::s				!Mandelshtam s
    real*8,intent(in),dimension(1:2)::qT		!(qtMin,qtMax)
    real*8,intent(in),dimension(1:2)::Q			!(Qmin,Qmax)
    real*8,intent(in),dimension(1:2)::y			!(ymin,ymax)
    logical,intent(in)::includeCuts			!include cuts
    real*8,intent(in),dimension(1:4)::CutParameters	!(p1,p2,eta1,eta2)
    integer,intent(in)::Num				!number of sections
    real*8::DY_xSec_Single
    real*8::X
    
    call xSec_DY(X,process,s,qT,Q,y,includeCuts,CutParameters,Num)
    DY_xSec_SingleN=X
  
  end function DY_xSec_SingleN
  
  function DY_xSec_Single(process,s,qT,Q,y,includeCuts,CutParameters)
    integer,intent(in),dimension(1:3)::process		!the number of process
    real*8,intent(in)::s				!Mandelshtam s
    real*8,intent(in),dimension(1:2)::qT		!(qtMin,qtMax)
    real*8,intent(in),dimension(1:2)::Q			!(Qmin,Qmax)
    real*8,intent(in),dimension(1:2)::y			!(ymin,ymax)
    logical,intent(in)::includeCuts			!include cuts
    real*8,intent(in),dimension(1:4)::CutParameters	!(p1,p2,eta1,eta2)
    real*8::DY_xSec_Single
    real*8::X
    
    call xSec_DY(X,process,s,qT,Q,y,includeCuts,CutParameters)
    DY_xSec_Single=X
  
  end function DY_xSec_Single
  
  function DY_xSec_List(process,s,qT,Q,y,includeCuts,CutParameters,ListLength)
    integer,intent(in)::ListLength
    integer,intent(in),dimension(:,:)::process			!the number of process
    real*8,intent(in),dimension(:)::s				!Mandelshtam s
    real*8,intent(in),dimension(:,:)::qT			!(qtMin,qtMax)
    real*8,intent(in),dimension(:,:)::Q				!(Qmin,Qmax)
    real*8,intent(in),dimension(:,:)::y				!(ymin,ymax)
    logical,intent(in),dimension(:)::includeCuts		!include cuts
    real*8,intent(in),dimension(:,:)::CutParameters	!(p1,p2,eta1,eta2)
    real*8,dimension(1:ListLength)::DY_xSec_List
    real*8,dimension(1:ListLength)::DY_xSec_List_test
    integer::i
    
    call xSec_DY_List(DY_xSec_List_test,process,s,qT,Q,y,includeCuts,CutParameters)
    if (TMDF_IsconvergenceLost()) then
       do i=1,ListLength
       DY_xSec_List(i)=1d9
       end do
    else
       do i=1,ListLength
       DY_xSec_List(i)=DY_xSec_List_test(i)
       end do
    end if
  
  end function DY_xSec_List
  
  function DY_xSec_mT_List(process,s,qT,Q,y,mT,includeCuts,CutParameters,ListLength)
    integer,intent(in)::ListLength
    integer,intent(in),dimension(:,:)::process			!the number of process
    real*8,intent(in),dimension(:)::s				!Mandelshtam s
    real*8,intent(in),dimension(:,:)::qT			!(qtMin,qtMax)
    real*8,intent(in),dimension(:,:)::Q				!(Qmin,Qmax)
    real*8,intent(in),dimension(:,:)::y				!(ymin,ymax)
    real*8,intent(in),dimension(:,:)::mT				!(ymin,ymax)
    logical,intent(in),dimension(:)::includeCuts		!include cuts
    real*8,intent(in),dimension(:,:)::CutParameters	!(p1,p2,eta1,eta2)
    real*8,dimension(1:ListLength)::DY_xSec_mT_List
    integer,dimension(1:ListLength)::n
    integer::i,NQ
    real*8::deltaQ
    
    deltaQ=Q(1,2)-Q(1,1)
    
    do i=1,ListLength
        if ((qT(i,2)-qT(i,1)-int(qT(i,2)-qT(i,1)))>0d0) then
           n(i)=4*int(qT(i,2)-qT(i,1)+1d0)
        else
           n(i)=4*int(qT(i,2)-qT(i,1))
        end if
    end do
    
    if (deltaQ>60d0.and.deltaQ<=120d0) then
        NQ=8
    else if (deltaQ>120d0) then
        NQ=16
    else
        NQ=4
    end if
    
    call xSec_DY_mT_List(DY_xSec_mT_List,process,s,qT,Q,y,mT,includeCuts,CutParameters,n,NQ)
  
  end function DY_xSec_mT_List
  
  function DY_xSec_BINLESS_List(process,s,qT,Q,y,includeCuts,CutParameters,ListLength)
    integer,intent(in)::ListLength
    integer,intent(in),dimension(:,:)::process			!the number of process
    real*8,intent(in),dimension(:)::s				!Mandelshtam s
    real*8,intent(in),dimension(:)::qT			!(qtMin,qtMax)
    real*8,intent(in),dimension(:)::Q				!(Qmin,Qmax)
    real*8,intent(in),dimension(:)::y				!(ymin,ymax)
    logical,intent(in),dimension(:)::includeCuts		!include cuts
    real*8,intent(in),dimension(:,:)::CutParameters	!(p1,p2,eta1,eta2)
    real*8,dimension(1:ListLength)::DY_xSec_BINLESS_List
    
    call xSec_DY_List_BINLESS(DY_xSec_BINLESS_List,process,s,qT,Q,y,includeCuts,CutParameters)
  
  end function DY_xSec_BINLESS_List
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SIDIS CROSS-SECTION
  
  function SIDIS_xSec_Single(process,s,pT,z,x,Q,doCut,Cuts)
    integer,intent(in),dimension(1:3)::process			!the number of process
    real*8,intent(in)::s					!Mandelshtam s
    real*8,intent(in),dimension(1:2)::pT			!(qtMin,qtMax)
    real*8,intent(in),dimension(1:2)::z				!(zmin,zmax)
    real*8,intent(in),dimension(1:2)::x				!(xmin,xmax)
    real*8,intent(in),dimension(1:2)::Q				!(Qmin,Qmax)    
    logical,intent(in)::doCut					!triger cuts
    real*8,intent(in),dimension(1:4)::Cuts			!(ymin,yMax,W2min,W2max)
    real*8::SIDIS_xSec_Single
    
    call xSec_SIDIS(SIDIS_xSec_Single,process,s,pT,z,x,Q,doCut,Cuts)
  
  end function SIDIS_xSec_Single
  
  function SIDIS_xSec_Single_withMasses(process,s,pT,z,x,Q,doCut,Cuts,masses)
    integer,intent(in),dimension(1:3)::process			!the number of process
    real*8,intent(in)::s					!Mandelshtam s
    real*8,intent(in),dimension(1:2)::pT			!(qtMin,qtMax)
    real*8,intent(in),dimension(1:2)::z				!(zmin,zmax)
    real*8,intent(in),dimension(1:2)::x				!(xmin,xmax)
    real*8,intent(in),dimension(1:2)::Q				!(Qmin,Qmax)    
    logical,intent(in)::doCut					!triger cuts
    real*8,intent(in),dimension(1:4)::Cuts			!(ymin,yMax,W2min,W2max)
    real*8,intent(in),dimension(1:2)::masses			!(mTARGET,mPRODUCT)
    real*8::SIDIS_xSec_Single_withMasses
    
    call xSec_SIDIS(SIDIS_xSec_Single_withMasses,process,s,pT,z,x,Q,doCut,Cuts,masses)
  
  end function SIDIS_xSec_Single_withMasses
  
  function SIDIS_xSec_List(process,s,pT,z,x,Q,doCut,Cuts,masses,ListLength)
    integer,intent(in)::ListLength
    integer,intent(in),dimension(:,:)::process			!the number of process
    real*8,intent(in),dimension(:)::s				!Mandelshtam s
    real*8,intent(in),dimension(:,:)::pT			!(qtMin,qtMax)
    real*8,intent(in),dimension(:,:)::z				!(zmin,zmax)
    real*8,intent(in),dimension(:,:)::x				!(xmin,xmax)
    real*8,intent(in),dimension(:,:)::Q				!(Qmin,Qmax)        
    logical,intent(in),dimension(:)::doCut			!triger cuts
    real*8,intent(in),dimension(:,:)::Cuts			!(ymin,yMax,W2min,W2max)
    real*8,intent(in),dimension(:,:)::masses			!(mTARGET,mPRODUCT)
    real*8,dimension(1:ListLength)::SIDIS_xSec_List
    
    call xSec_SIDIS_List_forharpy(SIDIS_xSec_List,process,s,pT,z,x,Q,doCut,Cuts,masses)
  
  end function SIDIS_xSec_List
  
  function SIDIS_xSec_BINLESS_List(process,s,pT,z,x,Q,masses,ListLength)
    integer,intent(in)::ListLength
    integer,intent(in),dimension(:,:)::process			!the number of process
    real*8,intent(in),dimension(:)::s				!Mandelshtam s
    real*8,intent(in),dimension(:)::pT			        !(qt)
    real*8,intent(in),dimension(:)::z				!(z)
    real*8,intent(in),dimension(:)::x				!(x)
    real*8,intent(in),dimension(:)::Q				!(Q)        
    real*8,intent(in),dimension(:,:)::masses			!(mTARGET,mPRODUCT)
    real*8,dimension(1:ListLength)::SIDIS_xSec_BINLESS_List
    
    call xSec_SIDIS_BINLESS_List_forharpy(SIDIS_xSec_BINLESS_List,process,s,pT,z,x,Q,masses)
  
  end function SIDIS_xSec_BINLESS_List
  
  
end module harpy
