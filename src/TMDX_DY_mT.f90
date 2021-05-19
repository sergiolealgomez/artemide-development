!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.01
!
!	Evaluation of the TMD cross-section for DY-like cross-sections
!	
!	if you use this module please, quote 1706.01473
!
!	ver 1.0: release (AV, 10.05.2017)
!	ver 1.1: multiple updates (AV, 5.10.2017)
!	ver 1.2: module is renamed, and multiple renaming of functions (AV, 15.10.2017)
!	ver 1.31: part of functions migrated to TMDF, rest updated (AV, 1.06.2018)
!	ver 1.4: encapsulation of cuts, and process,+ multiple updates (AV, 18.01.2019)
!	ver 2.01: Added Higgs xSec, piresum option, and coefficient function moved to separate file (AV, 17.06.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDX_DY_mT
use IO_functions
use TMDF
use LeptonCutsDY_mT
use QCDinput
use EWinput

implicit none
  private
  
   !Current version of module
 character (len=7),parameter :: moduleName="TMDX-DY-mT"
 character (len=5),parameter :: version="v2.01"
 !Last appropriate verion of constants-file
  integer,parameter::inputver=6
  
  real*8 :: tolerance=0.0005d0
  
  integer::outputlevel
  integer::messageTrigger
  
  !!!!
  !!!! in the module the kinematic is stored in the varibles "kinematic" real*8,dimension(1:6)
  !!!! which is (qT,s,Q,Q^2,x0,y,exp[y])
  !!!! where x0=sqrt[(Q^2+q_T^2)/s]   (if exactX1X2) or x0=Q^2/s (otherwise)
  !!!!
  
  !The variables for all parameters of the model!
  real*8:: s_global,Q_global,y_global,mT_global
  !! Set of process definition, for Prefactor 1, Prefactor 2, structure function, etc
  !! = (/p1,p2,p3/)
  integer,dimension(1:3)::process_global
  !!cut parameters
  !!!!! this variable = (pT1,pT2,etaMIN,etaMAX)
  !real*8,dimension(1:4)::CutParameters_global
  real*8,allocatable::CutParameters_global(:)
  
  !!other global parameters see SetXParameters  
  integer:: orderH_global
  logical:: includeCuts_global
  logical:: includeCuts_global_neutrino=.false.
  logical::usePIresum
  integer:: exactX1X2    !!!=1 if exact x's=true, =0 otherwise
  
  !!! number of sections for PT-integral by default
  integer::NumPTdefault=4
  
  real*8::c2_global!,muHard_global
  

  
  integer::GlobalCounter
  integer::CallCounter
  integer::messageCounter
  
  real*8::hc2
  
  logical::started=.false.
  
  public::alphaQCD
  
  public::TMDX_DY_mT_XSetup,TMDX_DY_mT_Initialize,TMDX_DY_mT_SetCuts,&
  TMDX_DY_mT_SetScaleVariation,TMDX_DY_mT_setProcess,TMDX_DY_mT_ShowStatistic,&
  TMDX_DY_mT_ResetCounters,TMDX_DY_mT_IsInitialized
  
  public::CalcXsec_DY_mT,CalcXsec_DY_mT_Yint,CalcXsec_DY_mT_Qint_Yint,&
  CalcXsec_DY_mT_PTint_Qint_Yint,CalcXsec_DY_mT_Qint,xSec_DY_mT,xSec_DY_mT_List,&
  CalcXsec_DY_PTint_Qint_Yint_mTint
  
 interface TMDX_DY_mT_SetCuts
    module procedure SetCuts_sym,SetCuts_asym,SetCuts_asym_leppT,SetCuts_asym_neppT,&
    SetCuts_asym_leppT_neupT
 end interface
  
 interface CalcXsec_DY_mT
    module procedure xSecSingle,xSecList
  end interface
  
 interface TMDX_DY_mT_setProcess
    module procedure TMDX_setProcess1,TMDX_setProcess3,TMDX_setProcess30
 end interface
 
 interface CalcXsec_DY_mT_Yint
    module procedure xSecSingle_Yint ,xSecList_Yint,xSecSingle_Ycomplete ,xSecList_Ycomplete
 end interface
 
  interface CalcXsec_DY_mT_Qint
    module procedure xSecSingle_Qint ,xSecList_Qint
 end interface
 
 interface CalcXsec_DY_mT_Qint_Yint
    module procedure xSecSingle_Qint_Yint, xSecList_Qint_Yint,&
    xSecSingle_Qint_Ycomplete,xSecList_Qint_Ycomplete
 end interface
 

 interface CalcXsec_DY_mT_PTint_Qint_Yint
    module procedure xSecSingle_PTint_Qint_Yint, xSecList_PTint_Qint_Yint, &
	    xSecSingle_PTint_Qint_Ycomplete, xSecList_PTint_Qint_Ycomplete,&
	    xSecSingle_PTintN_Qint_Yint, xSecList_PTintN_Qint_Yint,&
	    xSecSingle_PTintN_Qint_Ycomplete, xSecList_PTintN_Qint_Ycomplete,&
	    xSecListList_PTint_Qint_Yint,xSecListList_PTint_Qint_Ycomplete,&
	    xSecListList_PTintN_Qint_Yint,xSecListList_PTintN_Qint_Ycomplete,&
	    xSecListPY_PTint_Qint_Yint,xSecListPY_PTintN_Qint_Yint
    
 end interface
 
 interface CalcXsec_DY_PTint_Qint_Yint_mTint
    module procedure xSecSingle_PTint_Qint_Yint_mTint, xSecList_PTint_Qint_Yint_mTint,&
	    xSecSingle_PTintN_Qint_Yint_mTint, xSecList_PTintN_Qint_Yint_mTint,&
	    xSecListList_PTint_Qint_Yint_mTint, xSecListList_PTintN_Qint_Yint_mTint,&
	    xSecListPY_PTint_Qint_Yint_mTint, xSecListPY_PTintN_Qint_Yint_mTint,&
	    xSecSingle_PTintN_Qint_Ycomplete_mTint,&
	    xSecSingle_PTint_Qint_Ycomplete_mTint,&
	    xSecSingle_PTint_Qcomplete_Ycomplete_mTint,&
	    xSecSingle_PTintN_Qcomplete_Ycomplete_mTint,&
	    xSecSingle_PTint_Qcomplete_Ycomplete_mTint_omp,&
	    xSecSingle_PTint_Ycomplete_mTint_NRW,&
	    xSecSingle_PTintN_Ycomplete_mTint_NRW
    
 end interface
 
 interface xSec_DY_mT
    module procedure MainInterface_AsAAAloo,MainInterface_isAAAloo
 end interface
 
    contains
  
    function TMDX_DY_mT_IsInitialized()
        logical::TMDX_DY_mT_IsInitialized
        TMDX_DY_mT_IsInitialized=started
    end function TMDX_DY_mT_IsInitialized

   !! Initialization of the package
  subroutine TMDX_DY_mT_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path,line
    logical::initRequared,dummyLogical
    character(len=8)::orderMain
    integer::i,FILEver
    !$ integer:: omp_get_thread_num
    
    if(started) return
    
    if(present(prefix)) then
      path=trim(adjustl(prefix))//trim(adjustr(file))
    else
      path=trim(adjustr(file))
    end if
  
    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    !!! Search for output level
    call MoveTO(51,'*0   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) FILEver
    if(FILEver<inputver) then
      write(*,*) 'artemide.'//trim(moduleName)//': const-file version is too old.'
      write(*,*) '		     Update the const-file with artemide.setup'
      write(*,*) '  '
      stop
    end if
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>2) write(*,*) '--------------------------------------------- '
    if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger
    
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) hc2
    
    call MoveTO(51,'*13   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequared
    if(.not.initRequared) then
      if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not requared. '
      started=.false.
      return
    end if
    
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) orderMain
    SELECT CASE(orderMain)
      CASE ("LO")
	orderH_global=0
      CASE ("LO+")
	orderH_global=0
      CASE ("NLO")
	orderH_global=1
      CASE ("NLO+")
	orderH_global=1
      CASE ("NNLO")
	orderH_global=2
      CASE ("NNLO+")
	orderH_global=2
      CASE ("NNNLO")
	orderH_global=3
      CASE DEFAULT
	if(outputLevel>0) write(*,*) 'WARNING arTeMiDe.TMDX_DY_mT:try to set unknown order. Switch to NLO.'
	orderH_global=1
     END SELECT
     if(outputLevel>2) write(*,*) '	artemide.TMDX_DY_mT: the used order is ',trim(orderMain)
     
     !!exact values of x1x2
     call MoveTO(51,'*p2   ')
     read(51,*) dummyLogical
     if(dummyLogical) then 
      exactX1X2=1
     else
      exactX1X2=0
     end if
     if(outputLevel>2 .and. dummyLogical) write(*,*) '	artemide.TMDX_DY_mT: qT/Q correction for x1 and x2 variables are included.'
     call MoveTO(51,'*p3   ')
     read(51,*) usePIresum
     if(outputLevel>2 .and. usePIresum) write(*,*) '	artemide.TMDX_DY_mT: pi-resummation in coef.function included.'
     
     call MoveTO(51,'*B   ')
     call MoveTO(51,'*p1  ')
     read(51,*) tolerance
     call MoveTO(51,'*p2  ')
     read(51,*) NumPTdefault
     
!$    if(outputLevel>1) write(*,*) '	artemide.TMDX_DY_mT: parallel evaluation of cross-sections is to be used'
!$    call MoveTO(51,'*C   ')
!$    call MoveTO(51,'*p1  ')
!$    read(51,*) i
!$    call OMP_set_num_threads(i)
!$    if(outputLevel>1) write(*,*) '	artemide.TMDX_DY_mT: number of threads for parallel evaluation is set to ', i	

!$     if(outputLevel>2) write(*,*) '------TEST OF PARALLEL PROCESSING ----------'
!$OMP PARALLEL
!$     if(outputLevel>2) write(*,*) '   artemide.TMDX_DY_mT:thread num ',  omp_get_thread_num(), ' ready.'
!$OMP END PARALLEL
    CLOSE (51, STATUS='KEEP')
    
     if(.not.EWinput_IsInitialized()) then
	if(outputLevel>1) write(*,*) '.. initializing EWinput (from ',moduleName,')'
	if(present(prefix)) then
	  call EWinput_Initialize(file,prefix)
	else
	  call EWinput_Initialize(file)
	end if
      end if
      
      if(.not.TMDF_IsInitialized()) then
	if(outputLevel>1) write(*,*) '.. initializing TMDF (from ',moduleName,')'
	if(present(prefix)) then
	  call TMDF_Initialize(file,prefix)
	else
	  call TMDF_Initialize(file)
	end if
      end if
    
     includeCuts_global=.false.
     c2_global=1d0
     
     GlobalCounter=0
     CallCounter=0
     messageCounter=0
     
     started=.true.
    write(*,*)  color('----- arTeMiDe.TMD_DY_mT '//trim(version)//': .... initialized',c_green)
  end subroutine TMDX_DY_mT_Initialize

  subroutine TMDX_DY_mT_ResetCounters()
  if(outputlevel>2) call TMDX_DY_mT_ShowStatistic
  GlobalCounter=0
  CallCounter=0
  messageCounter=0
  end subroutine TMDX_DY_mT_ResetCounters
  
  subroutine TMDX_DY_mT_ShowStatistic()
  
      write(*,'(A,ES12.3)') 'TMDX DY-mT statistics      total calls of point xSec  :  ',Real(GlobalCounter)
      write(*,'(A,ES12.3)') '                              total calls of xSecF :  ',Real(CallCounter)
      write(*,'(A,F12.3)')  '                                         avarage M :  ',Real(GlobalCounter)/Real(CallCounter)
  end subroutine TMDX_DY_mT_ShowStatistic
  
  
  
  !!!!Call this after TMD initializetion but before NP, and X parameters
  subroutine TMDX_DY_mT_SetScaleVariation(c2_in)
    real*8::c1_in,c2_in,c3_in,c4_in
    
    if(outputLevel>1) write(*,*) 'TMDX_DY_mT: c2 scale reset:',c2_in
    
    if(c2_in<0.1d0 .or. c2_in>10.d0) then
    if(outputLevel>0) write(*,*) 'TMDX_DY_mT WARNING: variation in c2 is enourmous. c2 is set to 2'
     c2_global=2d0
    else
    c2_global=c2_in
    end if
    
  end subroutine TMDX_DY_mT_SetScaleVariation

  !!sets the cuts
  !! argument includeCuts_global_in=logical, if .true. will add to calculation the evaluation of cut leptonic tensor
  !! call BEFORE SetXParameters
  subroutine SetCuts_sym(include_arg,pT_arg,eta_min,eta_max)
  logical:: include_arg
  real*8:: pT_arg,eta_max,eta_min
  
  includeCuts_global=include_arg
  
  
  if (allocated(CutParameters_global)) then
  
  deallocate(CutParameters_global)
  
  end if
  
  allocate(CutParameters_global(1:4))
  
  CutParameters_global=(/pT_arg,pT_arg,eta_min,eta_max/)
  
  end subroutine SetCuts_sym
  
  !!sets the cuts (asymetric)
  !! argument includeCuts_global_in=logical, if .true. will add to calculation the evaluation of cut leptonic tensor
  !! call BEFORE SetXParameters
  subroutine SetCuts_asym(include_arg,pT1_arg,pT2_arg,eta_min,eta_max)
  logical:: include_arg
  real*8:: pT1_arg,pT2_arg,eta_max,eta_min
  
  includeCuts_global=include_arg
  
  
  if (allocated(CutParameters_global)) then
  
  deallocate(CutParameters_global)
  
  end if
  
  allocate(CutParameters_global(1:4))
  
  CutParameters_global=(/pT1_arg,pT2_arg,eta_min,eta_max/)
  
  end subroutine SetCuts_asym
  
  subroutine SetCuts_asym_leppT(include_arg,pT1_arg1,pT1_arg2,pT2_arg,eta_min,eta_max)
  logical:: include_arg
  real*8:: pT1_arg1,pT1_arg2,pT2_arg,eta_max,eta_min
  
  includeCuts_global=include_arg
  
  
  if (allocated(CutParameters_global)) then
  
  deallocate(CutParameters_global)
  
  end if
  
  allocate(CutParameters_global(1:5))
  
  CutParameters_global=(/pT1_arg1,pT1_arg2,pT2_arg,eta_min,eta_max/)
  
  end subroutine SetCuts_asym_leppT
  
  subroutine SetCuts_asym_neppT(include_arg,bin_neutrino,pT1_arg,pT2_arg1,pT2_arg2,eta_min,eta_max)
  logical:: include_arg
  logical:: bin_neutrino
  real*8:: pT1_arg,pT2_arg1,pT2_arg2,eta_max,eta_min
  
  includeCuts_global=include_arg
  includeCuts_global_neutrino=bin_neutrino
  
  
  if (allocated(CutParameters_global)) then
  
  deallocate(CutParameters_global)
  
  end if
  
  allocate(CutParameters_global(1:5))
  
  CutParameters_global=(/pT1_arg,pT2_arg1,pT2_arg2,eta_min,eta_max/)
  
  end subroutine SetCuts_asym_neppT
  
  
  subroutine SetCuts_asym_leppT_neupT(include_arg,pT1_arg1,pT1_arg2,pT2_arg1,pT2_arg2,eta_min,eta_max)
  logical:: include_arg
  real*8:: pT1_arg1,pT1_arg2,pT2_arg1,pT2_arg2,eta_max,eta_min
  
  includeCuts_global=include_arg
  
  
  if (allocated(CutParameters_global)) then
  
  deallocate(CutParameters_global)
  
  end if
  
  allocate(CutParameters_global(1:6))
    
  
  CutParameters_global=(/pT1_arg1,pT1_arg2,pT2_arg1,pT2_arg2,eta_min,eta_max/)
  

  
  end subroutine SetCuts_asym_leppT_neupT

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!! PROCESS DEFINITION
  
  function processArrayFromInteger(p)
    integer,intent(in)::p
    integer,dimension(1:3)::processArrayFromInteger
    SELECT CASE(p)
      case(1)
	processArrayFromInteger=(/1,1,5/) !!! p + p -> Z + gamma^*   (e.g. ATLAS, CMS, LHCb)
      case(2)
	processArrayFromInteger=(/1,1,6/) !!! p + pbar -> Z + gamma^*   (e.g. CDF,D0)
      case(4)
	processArrayFromInteger=(/1,2,1001/)
      case(5)
	processArrayFromInteger=(/1,1,5/)
      case(7)
	processArrayFromInteger=(/1,1,6/)
      case default
	write(*,*) 'ERROR: arTeMiDe_DY_mT: unknown process is called. p=',p
	write(*,*) 'Evaluation stop'
	stop
      end SELECT
  end function processArrayFromInteger
  
    !!!set variables for process definition
  subroutine TMDX_setProcess30(p0)
  integer,dimension(1:3)::p0
  
  process_global(1)=p0(1)
  process_global(2)=p0(2)
  process_global(3)=p0(3)
  end subroutine TMDX_setProcess30
  
  !!!set variables for process definition
  subroutine TMDX_setProcess3(p1,p2,p3)
  integer::p1,p2,p3
  
  process_global(1)=p1
  process_global(2)=p2
  process_global(3)=p3
  end subroutine TMDX_setProcess3
  
  !!!set variables for process definition
  subroutine TMDX_setProcess1(p)
  integer::p
  
  call TMDX_setProcess30(processArrayFromInteger(p))
  end subroutine TMDX_setProcess1
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS FOR OPERATION WITH KINEMATICS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !sets the main parameters of cross-section (x,zeta,etc)
  !the variables process defines the type of process
  subroutine TMDX_DY_mT_XSetup(s,Q,y,mT)
    real*8::s,Q,y,mT
    
    if(.not.started) then
    write(*,*) 'ERROR: arTeMiDe.TMDX_DY_mT is not initialized. Evaluation terminated'
    stop
    end if
    
    s_global=s
    Q_global=Q
    y_global=y
    mT_global=mT
    
!     if(includeCuts_global) then
!       call SetCutParameters(pT1_global,pT2_global,eta_min_global,eta_max_global)
!     end if
    
  end subroutine TMDX_DY_mT_XSetup
  
  !!! function makes kinematic array from the given set of qT,s,Q,y
  !!! array has 7 often appearing entries
  function kinematicArray(qT,s,Q,y,mT)
  real*8,dimension(1:8)::kinematicArray
  real*8::qT,s,Q,y,mT
  
  kinematicArray=(/qT,s,Q,Q**2,sqrt((Q**2+exactX1X2*qT**2)/s),y,exp(y),mT/)
  
  end function kinematicArray
  
  !!!intrinsic change the value of Q within kinematic array var
  subroutine SetQ(Q,var)
    real*8,dimension(1:8)::var
    real*8::Q
   
    var(3)=Q
    var(4)=Q**2
    var(5)=sqrt((var(4)+exactX1X2*var(1)**2)/var(2))
   
  end subroutine SetQ
  
  !!!intrinsic change the value of y within kinematic array var
  subroutine SetY(y,var)
    real*8,dimension(1:8)::var
    real*8::y
    
    var(6)=y
    var(7)=exp(y)
    
  end subroutine SetY
  
  !!!intrinsic change the value of qT within kinematic array var
  subroutine SetQT(qT_in,var)
    real*8,dimension(1:8)::var
    real*8::qT_in
    
    var(1)=qT_in
    var(5)=sqrt((var(4)+exactX1X2*var(1)**2)/var(2))
    
  end subroutine SetQT
  
  subroutine SetmT(mT,var)
    real*8,dimension(1:8)::var
    real*8::mT
    
    var(8)=mT
    
  end subroutine SetmT


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS FOR PREFACTORS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!! for easier read coeff-functions are split into separate file
  INCLUDE 'DYcoeff-func.f90'
  
  !!!!! Prefactor 2 is (universal part) x (cuts) x H
  function PreFactor2(kin,process, includeCuts_in,CutParam)
    real*8,dimension(1:8),intent(in)::kin
    logical,intent(in)::includeCuts_in
    real*8::PreFactor2,cutPrefactor,uniPart
    !real*8,dimension(1:4),intent(in)::CutParam
    real*8,intent(in)::CutParam(:)
    integer,dimension(1:3),intent(in)::process
    integer::length
    
    
    length=size(CutParam)
  
  !!!!! cut part
    if(includeCuts_in) then
       !!! here include cuts of lepton tensor
       if (process(1)==1) then
       !massless electron and muon
       cutPrefactor=CutFactor4(qT=kin(1),Q_in=kin(3),y_in=kin(6),mT=kin(8),CutParameters=CutParam,&
       neutrino=includeCuts_global_neutrino)
       else if (process(1)==3) then
       !massive electron
       cutPrefactor=CutFactor4(qT=kin(1),Q_in=kin(3),y_in=kin(6),mT=kin(8),CutParameters=CutParam,&
       neutrino=includeCuts_global_neutrino,lep_type=1)
       else if (process(1)==4) then
       !massive muon
       cutPrefactor=CutFactor4(qT=kin(1),Q_in=kin(3),y_in=kin(6),mT=kin(8),CutParameters=CutParam,&
       neutrino=includeCuts_global_neutrino,lep_type=2)
       else if (process(1)==5) then
       !massive electron+massive muon
       cutPrefactor=CutFactor4(qT=kin(1),Q_in=kin(3),y_in=kin(6),mT=kin(8),CutParameters=CutParam,&
       neutrino=includeCuts_global_neutrino,lep_type=2)+&
       CutFactor4(qT=kin(1),Q_in=kin(3),y_in=kin(6),mT=kin(8),CutParameters=CutParam,&
       neutrino=includeCuts_global_neutrino,lep_type=1)
       else
       !massless case
       cutPrefactor=CutFactor4(qT=kin(1),Q_in=kin(3),y_in=kin(6),mT=kin(8),CutParameters=CutParam,&
       neutrino=includeCuts_global_neutrino)
       end if
    
    else
	!!! this is uncut lepton tensor
       cutPrefactor=(1+0.5d0*(kin(1)/kin(3))**2)
    
    end if  
  
    
   !!!! universal part

  SELECT CASE(process(2))
    case(-10221191)
	uniPart=1d0
    CASE(1)
	!8 aEm^2/3 /Nc/Q^4/s !!! 
	uniPart=8d0*(alphaEM(kin(3))**2)/(3d0*kin(2)*kin(4)**2)*&
	    HardCoefficientDY(kin(3))*&
	    hc2*1d9!from GeV to pb
	!!! IsySymmetric=.true.  !!! state is IsySymmetric-function
    CASE(2)
	!4 pi aEm^2/3 /Nc/Q^2/s
	uniPart=1.3962634015954636d0*(alphaEM(kin(3))**2)/(kin(2)*kin(4))*&
	    HardCoefficientDY(kin(3))*&
	    hc2*1d9!from GeV to pb
	!!!IsySymmetric=.false.	!!! state is IsySymmetric-function
    CASE (3) !Zboson in the narrow-width approximation
	!4 pi^2 aem/Ns/s Br(z->ee+mumu)
	uniPart=13.15947253478581d0*alphaEM(kin(3))/kin(2)*&
	    HardCoefficientDY(kin(3))*&
	    hc2*1d9*&!from GeV to pb
	    0.03645d0!Br from PDG, ee+mumu
    CASE (4) !Wboson in the narrow-width approximation
	!(12 pi aEm/Nc) 1/(Q2**2*s) Br(z->enue+munumu) !!! 
	uniPart=12.566370614359172d0*(alphaEM(kin(3)))/(kin(2)*kin(4)**2)*&
	    HardCoefficientDY(kin(3))*&
	    hc2*1d9*&!from GeV to pb
	    0.1086d0!Br from PDG, enue+munumu
    !!! IsySymmetric=.true.  !!! state is IsySymmetric-function
    CASE (5) !exclusive HIGGSboson production
	! (2\pi) *pi Mh^2 as(mu)/36/s/vev^2 * H*cT^2
	! (1.033)^2 is correction for mT mass in Ct at LO.
	uniPart=(1d0/18d0)*MH2*(As(c2_global*kin(3))/VEVH)**2/kin(2)*&
	   HardCoefficientHIGGS(kin(3))*(EffCouplingHFF(kin(3))**2)*1.0677023627519822d0*&
	    hc2*1d9!from GeV to pb
	    
	cutPrefactor=1d0 !!! cut-prefactor is different in this case! 
    CASE DEFAULT 
      write(*,*) 'ERROR: arTeMiDe.TMDX_DY_mT: unknown process p2=',process(2),' .Evaluation stop.'
      stop
  END SELECT
  
  !!! this is case of xF integration the weight is 2sqrt[(Q^2+q_T^2)/s] Cosh[y]
  if(process(1)==2) then
    uniPart=2d0*kin(5)*cosh(kin(6))*uniPart
  end if

  PreFactor2=uniPart*cutPrefactor
  
  end function PreFactor2
  
  !!!! Set Prefactor1
  !!!! it multiplies the cross-section as a whole, does not participate in the integration
  function PreFactor1(p1)
  real*8::Prefactor1
  integer::p1
  SELECT CASE(p1)
    CASE(1)
	PreFactor1=1d0
    CASE(2)
	PreFactor1=1d0
	CASE(3)
	PreFactor1=1d0
	CASE(4)
	PreFactor1=1d0
	CASE(5) 
	PreFactor1=1d0
    CASE DEFAULT 
      write(*,*) 'ERROR: arTeMiDe.TMDX_DY_mT: unknown process p1=',p1,' .Evaluation stop.'
      stop
  END SELECT
  end function PreFactor1
  
  !!! check is the process y-symmetric
  function IsySymmetric(p2)
  logical::IsySymmetric
  integer::p2
  if(p2==1 .or. p2==3 ) then 
    IsySymmetric=.true.
  else
    IsySymmetric=.false.
  end if
  end function IsySymmetric
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS CALCULATING CROSS-SECTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !---------------------------------INTEGRATED------------------------------------------------------------------
  
  !!! this is help function which evaluate xSec at single qt (without lists) with only prefactor 2
  !!!! this is extended (and default) version of xSec, which include all parameters
  function xSec(var,process,incCut,CutParam)
    real*8:: xSec,FF
    real*8::x1,x2
    real*8,dimension(1:8),intent(in)::var
    logical,intent(in)::incCut
    !real*8,dimension(1:4),intent(in)::CutParam
    real*8,intent(in)::CutParam(:)
    integer,dimension(1:3),intent(in)::process
    real*8::X
    GlobalCounter=GlobalCounter+1
   
    x1=var(5)*var(7)
    x2=var(5)/var(7)
    
    X=PreFactor2(var,process,incCut,CutParam)
    
    if (X>0d0) then
    FF=TMDF_F(var(4),var(1),x1,x2,var(3)*c2_global,var(4),var(4),process(3))
    xSec=X*FF
    else
    xSec=0d0
    end if
    
  end function xSec
  
!---------------------------------INTEGRATED over mT---------------------------------------------------------------
  
function xSec_mTint_1(var,process,incCut,CutParam,mT_min,mT_max)
    real*8:: xSec_mTint_1
    real*8::x1,x2,FF
    real*8,dimension(1:8),intent(in)::var
    logical,intent(in)::incCut
    !real*8,dimension(1:4),intent(in)::CutParam
    real*8,intent(in)::CutParam(:)
    integer,dimension(1:3),intent(in)::process
    real*8::mT_min,mT_max,X_mT
    real*8::mT_b_max,mT_b_min
    integer::N_mT,j,i
    GlobalCounter=GlobalCounter+1
    
    
    X_mT=0d0
    if ((mT_max-mT_min)/60d0-real(int((mT_max-mT_min)/60d0))>0d0) then
        N_mT=int((mT_max-mT_min)/60d0)+1
        do j=1,N_mT
            mT_b_min=mT_min+real(j-1)*(mT_max-mT_min)/real(N_mT)
            mT_b_max=mT_min+real(j)*(mT_max-mT_min)/real(N_mT)
            X_mT=X_mT+xSec_mTint_1_S(var,process,incCut,CutParam,mT_b_min,mT_b_max)
        end do
    else
        N_mT=int((mT_max-mT_min)/60d0)
        do j=1,N_mT
            mT_b_min=mT_min+real(j-1)*(mT_max-mT_min)/real(N_mT)
            mT_b_max=mT_min+real(j)*(mT_max-mT_min)/real(N_mT)
            X_mT=X_mT+xSec_mTint_1_S(var,process,incCut,CutParam,mT_b_min,mT_b_max)
        end do
    end if
    
    x1=var(5)*var(7)
    x2=var(5)/var(7)
    
    if (X_mT>0d0) then
    FF=TMDF_F(var(4),var(1),x1,x2,var(3)*c2_global,var(4),var(4),process(3))
    xSec_mTint_1=X_mT*FF
    else
    xSec_mTint_1=0d0
    end if
    
  
end function xSec_mTint_1

function xSec_mTint_1_S(var,process,incCut,CutParam,mT_min,mT_max)
    real*8,dimension(1:8),intent(in)::var
    integer,dimension(1:3),intent(in)::process
    logical,intent(in)::incCut
    real*8,intent(in)::CutParam(:)
    real*8::mT_min,mT_max,delta_mT
    real*8::X1,X2,X3,X4,X5
    real*8::mT2,mT3,mT4
    real*8::value_Max,xSec_mTint_1_S
    
    delta_mT=(mT_max-mT_min)/4d0
    mT2=mT_min+delta_mT
    mT3=mT_min+2d0*delta_mT
    mT4=mT_max-delta_mT
    
    
    call SetmT(mT_min,var)
    X1=2d0*mT_min*PreFactor2(var,process,incCut,CutParam)

    call SetmT(mT2,var)
    X2=2d0*mT2*PreFactor2(var,process,incCut,CutParam)
    
    call SetmT(mT3,var)
    X3=2d0*mT3*PreFactor2(var,process,incCut,CutParam)

    call SetmT(mT4,var)
    X4=2d0*mT4*PreFactor2(var,process,incCut,CutParam)

    call SetmT(mT_max,var)
    X5=2d0*mT_max*PreFactor2(var,process,incCut,CutParam)
    
    value_Max=delta_mT*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/3d0
    
    xSec_mTint_1_S=xSec_mTint_1_S_Rec(var,process,incCut,cutParam,mT_min,mT3,X1,X2,X3,value_Max)+&
    xSec_mTint_1_S_Rec(var,process,incCut,cutParam,mT3,mT_max,X3,X4,X5,value_Max)
    
end function xSec_mTint_1_S

recursive function xSec_mTint_1_S_Rec(var,process,incCut,CutParam,mT_min,mT_max,X1,X3,X5,&
                    value_Max) result(Inter_mT_Rec)
    real*8,dimension(1:8),intent(in)::var
    integer,dimension(1:3),intent(in)::process
    real*8,intent(in)::CutParam(:)
    logical,intent(in)::incCut
    real*8,intent(in)::value_Max
    real*8::mT_min,mT_max,delta_mT
    real*8::X1,X2,X3,X4,X5,Inter_mT_Rec
    real*8::mT2,mT3,mT4
    real*8::valueAB,valueACB
    
    delta_mT=(mT_max-mT_min)/4d0
    mT2=mT_min+delta_mT
    mT3=mT_min+2d0*delta_mT
    mT4=mT_max-delta_mT
    
    call SetmT(mT2,var)
    X2=2d0*mT2*PreFactor2(var,process,incCut,CutParam)

    call SetmT(mT4,var)
    X4=2d0*mT4*PreFactor2(var,process,incCut,CutParam)

    valueAB=(mT_max-mT_min)*(X1+4d0*X3+X5)/6d0
    valueACB=delta_mT*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/3d0
    
    if (value_Max>0d0) then
        if (ABS((valueACB-valueAB)/value_Max)>tolerance) then
        Inter_mT_Rec=xSec_mTint_1_S_Rec(var,process,incCut,cutParam,mT_min,mT3,X1,X2,X3,value_Max)+&
        xSec_mTint_1_S_Rec(var,process,incCut,cutParam,mT3,mT_max,X3,X4,X5,value_Max)
        else
        Inter_mT_Rec=valueACB
        end if
    else
        if (ABS((valueACB-valueAB))>0.000001d0) then
        Inter_mT_Rec=xSec_mTint_1_S_Rec(var,process,incCut,cutParam,mT_min,mT3,X1,X2,X3,value_Max)+&
        xSec_mTint_1_S_Rec(var,process,incCut,cutParam,mT3,mT_max,X3,X4,X5,value_Max)
        else
        Inter_mT_Rec=valueACB
        end if
    end if
    
end function xSec_mTint_1_S_Rec
!------------------------------------------------------------------------------------------------
  
  !---------------------------------INTEGRATED over Y---------------------------------------------------------------
  
  !!! function determines the best value of PT-sections from PT-bin size, and Q
  !!! it is determined by formula Q/PT< val/ (2 k) => def+2K
  function NumPT_auto(dPT,Q)
    real,parameter::val=40.
    real::dPT,Q,rat
    integer::i,NumPT_auto
    rat=Q/dPT
    
    if(rat>40.) then
        NumPT_auto=NumPTdefault
        return
    else
        do i=1,5
            if(rat>(40./2./i)) then
                NumPT_auto=NumPTdefault+2*i
                return
            end if
        end do
    end if
    if(outputlevel>1) write(*,*) 'arTeMiDe_DY_mT:WARNING! Fail to automatically determine number of Pt-section for a bin. &
                                                Possibly Pt-bin is too large', dPT
    NumPT_auto=NumPTdefault+12
    
  end function NumPT_auto
      
  function yFromXF(xF,var)
  real*8,dimension(1:8)::var
  real*8:: yFromXF,xF
!     yFromXF=asinh(Sqrt(s_global/(Q2_global+qt**2))*xF/2d0)
    yFromXF=asinh(xF/2d0/var(5))
  end function yFromXF
  
  !!!
  function Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
    real*8,dimension(1:8) :: var
    logical,intent(in)::incCut
    !real*8,dimension(1:4),intent(in)::CutParam
    real*8,intent(in)::CutParam(:)
    real*8 :: Xsec_Yint
    real*8 :: ymin, ymax,ymin_in,ymax_in
    real*8 :: ymin_Check,ymax_Check
    integer,dimension(1:3),intent(in)::process
  
    
    !!! evaluate correspnding y's
    !!! in the case process=2 the integral is over xF
    if(process(1)==2) then
      ymin=yFromXF(ymin_in,var)
      ymax=yFromXF(ymax_in,var)
    else
      ymin=ymin_in
      ymax=ymax_in
   end if
    
    ymin_Check=log(var(5))+0.000000001d0
    ymax_Check=-log(var(5))-0.000000001d0
    
    if(IsySymmetric(process(2)) .and. (ABS(ymax+ymin)<tolerance)) then!!! symetric integral
    if(ymax > ymax_check) then
        ymax=ymax_Check
    end if!!!!! else case: automatically taken into account
    
    Xsec_Yint=2d0*integralOverYpoint_S(var,process,incCut,CutParam,0d0,ymax)!!! 2 since symmetric
      
    else !!!non-symmetric integral!!!!!!!!
    
      if(ymax<ymin_check .or. ymin>ymax_check) then !!! the case then y is outside physicsl region
	  Xsec_Yint=0d0
      else
	  if(ymax > ymax_check) then
	  ymax=yMax_check
	  end if!!!!! else case: automatically taken into account
	  if(ymin < ymin_check) then
	  ymin=ymin_check
	  end if!!!!! else case: automatically taken into account
	
      Xsec_Yint=integralOverYpoint_S(var,process,incCut,CutParam,ymin,ymax)
      
      end if
  end if
  
  
  end function Xsec_Yint
  !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the approximate value of integral to weight the tolerance.
  !!!! evaluation is done by adaptive simpson
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  function integralOverYpoint_S(var,process,incCut,CutParam,yMin_in,yMax_in)
   real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   !real*8,dimension(1:4),intent(in)::CutParam
   real*8,intent(in)::CutParam(:)
   integer,dimension(1:3),intent(in)::process
   real*8 ::integralOverYpoint_S
   real*8 :: X1,X2,X3,X4,X5
   real*8 :: y2,y3,y4,deltay
   real*8 :: yMin_in,yMax_in
   real*8::valueMax
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
      
   call SetY(yMin_in,var)
   X1= xSec(var,process,incCut,CutParam)   
   call SetY(y2,var)
   X2= xSec(var,process,incCut,CutParam)   
   call SetY(y3,var)
   X3= xSec(var,process,incCut,CutParam)   
   call SetY(y4,var)
   X4= xSec(var,process,incCut,CutParam)   
   call SetY(yMax_in,var)
   X5= xSec(var,process,incCut,CutParam)
   
   !!approximate integral value
   valueMax=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   integralOverYpoint_S=IntegralOverYpoint_S_Rec(var,process,incCut,CutParam,yMin_in,y3,X1,X2,X3,valueMax)+&
	  IntegralOverYpoint_S_Rec(var,process,incCut,CutParam,y3,yMax_in,X3,X4,X5,valueMax)
  end function integralOverYpoint_S
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function integralOverYpoint_S_Rec(var,process,incCut,CutParam,yMin_in,yMax_in,X1,X3,X5,valueMax) result(interX)
   real*8,dimension(1:8) ::var
   logical,intent(in)::incCut
   !real*8,dimension(1:4),intent(in)::CutParam
   real*8,intent(in)::CutParam(:)
   integer,dimension(1:3),intent(in)::process
   real*8 :: interX,X1,X2,X3,X4,X5
   real*8 :: value,valueAB,valueACB
   real*8 :: yMin_in,yMax_in,y2,y3,y4,deltay
   real*8,intent(in)::valueMax
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
   
   
   call SetY(y2,var)
   X2= xSec(var,process,incCut,CutParam)
   
   call SetY(y4,var)
   X4= xSec(var,process,incCut,CutParam)
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
  
   
   !write(*,*) yMin_in,yMax_in,(valueACB-valueAB)/valueMax,ABS((valueACB-valueAB)/valueMax)>tolerance
   !if(GlobalCounter>200) stop
!    if(abs(valueMax)<1d-10) write(*,*) valueMax
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=integralOverYpoint_S_Rec(var,process,incCut,CutParam,yMin_in,y3,X1,X2,X3,valueMax)&
	  +integralOverYpoint_S_Rec(var,process,incCut,CutParam,y3,yMax_in,X3,X4,X5,valueMax)
   else
    interX=valueACB
   end if
   
  end function integralOverYpoint_S_Rec
  
  
  
  !!!
  function Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,mT_min,mT_max)
    real*8,dimension(1:8) :: var
    logical,intent(in)::incCut
    !real*8,dimension(1:4),intent(in)::CutParam
    real*8,intent(in)::CutParam(:)
    real*8 :: Xsec_Yint_mTint_1
    real*8 :: ymin, ymax,ymin_in,ymax_in,mT_min,mT_max
    real*8 :: ymin_Check,ymax_Check
    integer,dimension(1:3),intent(in)::process
  
    
    !!! evaluate correspnding y's
    !!! in the case process=2 the integral is over xF
    if(process(1)==2) then
      ymin=yFromXF(ymin_in,var)
      ymax=yFromXF(ymax_in,var)
    else
      ymin=ymin_in
      ymax=ymax_in
   end if
    
    ymin_Check=log(var(5))+0.000000001d0
    ymax_Check=-log(var(5))-0.000000001d0
    
    if(IsySymmetric(process(2)) .and. (ABS(ymax+ymin)<tolerance)) then!!! symetric integral
    if(ymax > ymax_check) then
        ymax=ymax_Check
    end if!!!!! else case: automatically taken into account

    Xsec_Yint_mTint_1=2d0*integralOverYpointmTint_S_1(var,process,incCut,CutParam,0d0,ymax,mT_min,mT_max)!!! 2 since symmetric

    else !!!non-symmetric integral!!!!!!!!
    
      if(ymax<ymin_check .or. ymin>ymax_check) then !!! the case then y is outside physicsl region
	  Xsec_Yint_mTint_1=0d0
      else
	  if(ymax > ymax_check) then
	  ymax=yMax_check
	  end if!!!!! else case: automatically taken into account
	  if(ymin < ymin_check) then
	  ymin=ymin_check
	  end if!!!!! else case: automatically taken into account
	
      Xsec_Yint_mTint_1=integralOverYpointmTint_S_1(var,process,incCut,CutParam,ymin,ymax,mT_min,mT_max)
      
      end if
  end if

  end function Xsec_Yint_mTint_1
  
  !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the approximate value of integral to weight the tolerance.
  !!!! evaluation is done by adaptive simpson
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  function integralOverYpointmTint_S_1(var,process,incCut,CutParam,yMin_in,yMax_in,mT_min,mT_max)
   real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   !real*8,dimension(1:4),intent(in)::CutParam
   real*8,intent(in)::CutParam(:)
   integer,dimension(1:3),intent(in)::process
   real*8 ::integralOverYpointmTint_S_1
   real*8 :: X1,X2,X3,X4,X5
   real*8 :: y2,y3,y4,deltay
   real*8 :: yMin_in,yMax_in,mT_min,mT_max
   real*8::valueMax
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
      
   call SetY(yMin_in,var)
   X1= xSec_mTint_1(var,process,incCut,CutParam,mT_min,mT_max)   
   call SetY(y2,var)
   X2= xSec_mTint_1(var,process,incCut,CutParam,mT_min,mT_max)   
   call SetY(y3,var)
   X3= xSec_mTint_1(var,process,incCut,CutParam,mT_min,mT_max)   
   call SetY(y4,var)
   X4= xSec_mTint_1(var,process,incCut,CutParam,mT_min,mT_max)   
   call SetY(yMax_in,var)
   X5= xSec_mTint_1(var,process,incCut,CutParam,mT_min,mT_max)
   
   !!approximate integral value
   valueMax=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   integralOverYpointmTint_S_1=IntegralOverYpointmTint_S_Rec_1(var,process,incCut,CutParam,&
   yMin_in,y3,mT_min,mT_max,X1,X2,X3,valueMax)+&
   IntegralOverYpointmTint_S_Rec_1(var,process,incCut,CutParam,y3,yMax_in,&
   mT_min,mT_max,X3,X4,X5,valueMax)
  end function integralOverYpointmTint_S_1
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function IntegralOverYpointmTint_S_Rec_1(var,process,incCut,CutParam,&
                    yMin_in,yMax_in,mT_min,mT_max,X1,X3,X5,valueMax) result(interX)
   real*8,dimension(1:8) ::var
   logical,intent(in)::incCut
   !real*8,dimension(1:4),intent(in)::CutParam
   real*8,intent(in)::CutParam(:)
   integer,dimension(1:3),intent(in)::process
   real*8 :: interX,X1,X2,X3,X4,X5
   real*8 :: value,valueAB,valueACB
   real*8 :: yMin_in,yMax_in,y2,y3,y4,deltay,mT_min,mT_max
   real*8,intent(in)::valueMax
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
   
   
   call SetY(y2,var)
   X2= xSec_mTint_1(var,process,incCut,CutParam,mT_min,mT_max)
   
   call SetY(y4,var)
   X4= xSec_mTint_1(var,process,incCut,CutParam,mT_min,mT_max)
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
  
   
   !write(*,*) yMin_in,yMax_in,(valueACB-valueAB)/valueMax,ABS((valueACB-valueAB)/valueMax)>tolerance
   !if(GlobalCounter>200) stop
!    if(abs(valueMax)<1d-10) write(*,*) valueMax
   
   if (abs(valueMax).gt.0d0) then
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=IntegralOverYpointmTint_S_Rec_1(var,process,incCut,CutParam,yMin_in,y3,&
    mT_min,mT_max,X1,X2,X3,valueMax)+&
    IntegralOverYpointmTint_S_Rec_1(var,process,incCut,CutParam,y3,yMax_in,&
    mT_min,mT_max,X3,X4,X5,valueMax)
   else
    interX=valueACB
   end if
   
   else
   
    If(ABS(valueACB-valueAB)>0.000001d0) then
    interX=IntegralOverYpointmTint_S_Rec_1(var,process,incCut,CutParam,yMin_in,y3,&
    mT_min,mT_max,X1,X2,X3,valueMax)+&
    IntegralOverYpointmTint_S_Rec_1(var,process,incCut,CutParam,y3,yMax_in,&
    mT_min,mT_max,X3,X4,X5,valueMax)
   else
    interX=valueACB
   end if
   
   end if
   
  end function IntegralOverYpointmTint_S_Rec_1
  
  !---------------------------------INTEGRATED over Q---------------------------------------------------------------
  function Xsec_Qint(var,process,incCut,CutParam,Q_min,Q_max)
    real*8,dimension(1:8)::var
    logical,intent(in)::incCut
    real*8,dimension(1:4),intent(in)::CutParam
    integer,dimension(1:3),intent(in)::process
    real*8:: Xsec_Qint
    real*8:: Q_min,Q_max
    
    Xsec_Qint=integralOverQpoint_S(var,process,incCut,CutParam,Q_min,Q_max)
  end function Xsec_Qint
  
  !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  !!!! taking into account minimum calls of y-integral we have  =81 points
  function integralOverQpoint_S(var,process,incCut,CutParam,QMin_in,QMax_in)
  real*8,dimension(1:8)::var
  logical,intent(in)::incCut
   real*8,dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
    real*8 ::integralOverQpoint_S
   real*8 :: X1,X2,X3,X4,X5
   real*8 :: QMin_in,QMax_in
   real*8::valueMax,Q2,Q3,Q4,deltaQ
   
   deltaQ=QMax_in-QMin_in
   Q2=QMin_in+deltaQ/4d0
   Q3=QMin_in+deltaQ/2d0
   Q4=QMax_in-deltaQ/4d0
   
   call SetQ(QMin_in,var)
   X1=2*QMin_in*xSec(var,process,incCut,CutParam)
   
   call SetQ(Q2,var)
   X2=2*Q2*xSec(var,process,incCut,CutParam)
   
   call SetQ(Q3,var)
   X3=2*Q3*xSec(var,process,incCut,CutParam)
   
   call SetQ(Q4,var)
   X4=2*Q4*xSec(var,process,incCut,CutParam)
   
   call SetQ(QMax_in,var)
   X5=2*Qmax_in*xSec(var,process,incCut,CutParam)
   
      !!approximate integral value
   valueMax=deltaQ*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   integralOverQpoint_S=IntegralOverQpoint_S_Rec(var,process,incCut,CutParam,QMin_in,Q3,X1,X2,X3,valueMax)+&
	IntegralOverQpoint_S_Rec(var,process,incCut,CutParam,Q3,QMax_in,X3,X4,X5,valueMax)
  end function integralOverQpoint_S
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function integralOverQpoint_S_Rec(var,process,incCut,CutParam,QMin_in,QMax_in,X1,X3,X5,valueMax) result(interX)
   real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   real*8,dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real*8 :: interX,X1,X2,X3,X4,X5
   real*8 :: valueAB,valueACB
   real*8 :: QMin_in,QMax_in,Q2,Q3,Q4,deltaQ
   real*8,intent(in)::valueMax
   
   deltaQ=QMax_in-QMin_in
   Q2=QMin_in+deltaQ/4d0
   Q3=QMin_in+deltaQ/2d0
   Q4=QMax_in-deltaQ/4d0
   
   call SetQ(Q2,var)
   X2=2*Q2*xSec(var,process,incCut,CutParam)
      
   call SetQ(Q4,var)
   X4=2*Q4*xSec(var,process,incCut,CutParam)
   
   valueAB=deltaQ*(X1+4d0*X3+X5)/6d0
   valueACB=deltaQ*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=integralOverQpoint_S_Rec(var,process,incCut,CutParam,QMin_in,Q3,X1,X2,X3,valueMax)&
	  +integralOverQpoint_S_Rec(var,process,incCut,CutParam,Q3,Qmax_in,X3,X4,X5,valueMax)
   else
    interX=valueACB
   end if
  end function integralOverQpoint_S_Rec
  

  
  !---------------------------------INTEGRATED over Y over Q  16---------------------------------------------------------------
  !!!! No need for check over Y they take a place within y-integration
    !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  !!!! taking into account minimum calls of y-integral we have  =81 points
  function Xsec_Qint_Yint_16(var,process,incCut,CutParam,QMin_in,QMax_in,ymin_in,ymax_in)
  real*8,dimension(1:8)::var
  logical,intent(in)::incCut
   real*8,dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real*8:: Xsec_Qint_Yint_16
   real*8 :: ymin_in,ymax_in,QMin_in,QMax_in
   real*8::valueMax,deltaQ,medQ
   real*8,dimension(1:8)::Qinter8
   real*8::Qinter,Qinter_inc,Qinter_25,Qinter_med,Qinter_75,Qinter_fin
   real*8::Qinter_125,Qinter_375,Qinter_625,Qinter_875,Qvalue
   integer::i
   
   
   
    deltaQ=(QMax_in-QMin_in)/16d0
    medQ=QMin_in+8d0*deltaQ

    !!approximate integral value
      
    call SetQ(QMin_in,var)
    Qinter_inc=2d0*QMin_in*Xsec_Yint(var,process,incCut,CutParam,&
    ymin_in,ymax_in)
    Qinter=Qinter_inc
    

    
    do i=1,8
        Qvalue=QMin_in+(2d0*real(i)-1d0)*deltaQ
        call SetQ(Qvalue,var)
        Qinter8(i)=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        ymin_in,ymax_in)
        Qinter=Qinter+4d0*Qinter8(i)
    end do
    
  
        Qvalue=QMin_in+2d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_125=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        ymin_in,ymax_in)
        Qinter=Qinter+2d0*Qinter_125
        
        Qvalue=QMin_in+4d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_25=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        ymin_in,ymax_in)
        Qinter=Qinter+2d0*Qinter_25
        
        Qvalue=QMin_in+6d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_375=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        ymin_in,ymax_in)
        Qinter=Qinter+2d0*Qinter_375
        
        Qvalue=QMin_in+8d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_med=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        ymin_in,ymax_in)
        Qinter=Qinter+2d0*Qinter_med
        
        Qvalue=QMin_in+10d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_625=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        ymin_in,ymax_in)
        Qinter=Qinter+2d0*Qinter_625
        
        Qvalue=QMin_in+12d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_75=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        ymin_in,ymax_in)
        Qinter=Qinter+2d0*Qinter_75
        
        Qvalue=QMin_in+14d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_875=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        ymin_in,ymax_in)
        Qinter=Qinter+2d0*Qinter_875
        


    
    call SetQ(QMax_in,var)
    Qinter_fin=2d0*QMax_in*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
   
    valueMax=deltaQ*(Qinter+Qinter_fin)/3d0
   
   Xsec_Qint_Yint_16=IntegralOverQYpoint_S_Rec_16(var,process,incCut,CutParam,QMin_in,medQ,&
   ymin_in,ymax_in,Qinter_inc,Qinter8(1),Qinter_125,Qinter8(2),Qinter_25,Qinter8(3),Qinter_375,Qinter8(4),&
   Qinter_med,valueMax)+IntegralOverQYpoint_S_Rec_16(var,process,incCut,CutParam,&
   medQ,QMax_in,ymin_in,ymax_in,Qinter_med,Qinter8(5),Qinter_625,Qinter8(6),Qinter_75,Qinter8(7),&
   Qinter_875,Qinter8(8),Qinter_fin,valueMax)
   
   
  end function Xsec_Qint_Yint_16
  

  recursive function IntegralOverQYpoint_S_Rec_16(var,process,incCut,CutParam,&
			      QMin_in,QMax_in,yMin_in,yMax_in,Qinter_inc,Qinter_125,Qinter_25,Qinter_375,Qinter_med,Qinter_625,&
			      Qinter_75,Qinter_875,Qinter_fin,valueMax) result(interX)
   real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   real*8,dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real*8 :: interX
   real*8 :: valueAB,valueACB
   real*8 :: yMin_in,yMax_in,QMin_in,QMax_in,deltaQ
   real*8,intent(in)::valueMax
   real*8,dimension(1:8)::Qinter8
   real*8::Qinter,Qinter_inc,Qinter_25,Qinter_med,Qinter_75,Qinter_fin,medQ
   real*8::Qinter_125,Qinter_375,Qinter_625,Qinter_875,Qvalue
   integer::i
   
   deltaQ=(QMax_in-QMin_in)/16d0
   medQ=QMin_in+8d0*deltaQ
   
   
    Qinter=Qinter_inc+Qinter_fin
    
    do i=1,8
        Qvalue=QMin_in+(2d0*real(i)-1d0)*deltaQ
        call SetQ(Qvalue,var)
        Qinter8(i)=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        yMin_in,yMax_in)
        Qinter=Qinter+4d0*Qinter8(i)
    end do
    
    
   
   valueAB=((QMax_in-QMin_in)/6d0)*(Qinter_inc+4d0*Qinter_med+Qinter_fin)
   valueACB=deltaQ*(Qinter+2d0*Qinter_125+2d0*Qinter_25+2d0*Qinter_375+2d0*Qinter_med+2d0*Qinter_625+&
   2d0*Qinter_75+2d0*Qinter_875)/3d0
   
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=IntegralOverQYpoint_S_Rec_16(var,process,incCut,CutParam,QMin_in,medQ,yMin_in,yMax_in,&
    Qinter_inc,Qinter8(1),Qinter_125,Qinter8(2),Qinter_25,Qinter8(3),Qinter_375,Qinter8(4),Qinter_med,valueMax)&
	  +IntegralOverQYpoint_S_Rec_16(var,process,incCut,CutParam,medQ,QMax_in,yMin_in,yMax_in,&
	  Qinter_med,Qinter8(5),Qinter_625,Qinter8(6),Qinter_75,Qinter8(7),Qinter_875,Qinter8(8),Qinter_fin,valueMax)
   else
    interX=valueACB
   end if
  end function IntegralOverQYpoint_S_Rec_16
  
  
  
  !---------------------------------INTEGRATED over Y over Q  8---------------------------------------------------------------
  !!!! No need for check over Y they take a place within y-integration
    !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  !!!! taking into account minimum calls of y-integral we have  =81 points
  function Xsec_Qint_Yint_8(var,process,incCut,CutParam,QMin_in,QMax_in,ymin_in,ymax_in)
  real*8,dimension(1:8)::var
  logical,intent(in)::incCut
   real*8,dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real*8:: Xsec_Qint_Yint_8
   real*8 :: ymin_in,ymax_in,QMin_in,QMax_in
   real*8::valueMax,deltaQ,medQ
   real*8,dimension(1:4)::Qinter8
   real*8::Qinter,Qinter_inc,Qinter_25,Qinter_med,Qinter_75,Qinter_fin
   real*8::Qvalue
   integer::i
   
   
   
    deltaQ=(QMax_in-QMin_in)/8d0
    medQ=QMin_in+4d0*deltaQ

    !!approximate integral value
      
    call SetQ(QMin_in,var)
    Qinter_inc=2d0*QMin_in*Xsec_Yint(var,process,incCut,CutParam,&
    ymin_in,ymax_in)
    Qinter=Qinter_inc
    

    
    do i=1,4
        Qvalue=QMin_in+(2d0*real(i)-1d0)*deltaQ
        call SetQ(Qvalue,var)
        Qinter8(i)=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        ymin_in,ymax_in)
        Qinter=Qinter+4d0*Qinter8(i)
    end do
    
        
        Qvalue=QMin_in+2d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_25=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        ymin_in,ymax_in)
        Qinter=Qinter+2d0*Qinter_25
        
        Qvalue=QMin_in+4d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_med=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        ymin_in,ymax_in)
        Qinter=Qinter+2d0*Qinter_med
        
        
        Qvalue=QMin_in+6d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_75=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        ymin_in,ymax_in)
        Qinter=Qinter+2d0*Qinter_75


    
    call SetQ(QMax_in,var)
    Qinter_fin=2d0*QMax_in*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
   
    valueMax=deltaQ*(Qinter+Qinter_fin)/3d0
   
   Xsec_Qint_Yint_8=IntegralOverQYpoint_S_Rec_8(var,process,incCut,CutParam,QMin_in,medQ,&
   ymin_in,ymax_in,Qinter_inc,Qinter8(1),Qinter_25,Qinter8(2),&
   Qinter_med,valueMax)+IntegralOverQYpoint_S_Rec_8(var,process,incCut,CutParam,&
   medQ,QMax_in,ymin_in,ymax_in,Qinter_med,Qinter8(3),Qinter_75,Qinter8(4),Qinter_fin,valueMax)
   
   
  end function Xsec_Qint_Yint_8
  

  recursive function IntegralOverQYpoint_S_Rec_8(var,process,incCut,CutParam,&
			      QMin_in,QMax_in,yMin_in,yMax_in,Qinter_inc,Qinter_25,Qinter_med,&
			      Qinter_75,Qinter_fin,valueMax) result(interX)
   real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   real*8,dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real*8 :: interX
   real*8 :: valueAB,valueACB
   real*8 :: yMin_in,yMax_in,QMin_in,QMax_in,deltaQ
   real*8,intent(in)::valueMax
   real*8,dimension(1:4)::Qinter8
   real*8::Qinter,Qinter_inc,Qinter_25,Qinter_med,Qinter_75,Qinter_fin,medQ
   real*8::Qvalue
   integer::i
   
   deltaQ=(QMax_in-QMin_in)/8d0
   medQ=QMin_in+4d0*deltaQ
   
   
    Qinter=Qinter_inc+Qinter_fin
    
    do i=1,4
        Qvalue=QMin_in+(2d0*real(i)-1d0)*deltaQ
        call SetQ(Qvalue,var)
        Qinter8(i)=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        yMin_in,yMax_in)
        Qinter=Qinter+4d0*Qinter8(i)
    end do
    
    
   
   valueAB=((QMax_in-QMin_in)/6d0)*(Qinter_inc+4d0*Qinter_med+Qinter_fin)
   valueACB=deltaQ*(Qinter+2d0*Qinter_25+2d0*Qinter_med+2d0*Qinter_75)/3d0
   
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=IntegralOverQYpoint_S_Rec_8(var,process,incCut,CutParam,QMin_in,medQ,yMin_in,yMax_in,&
    Qinter_inc,Qinter8(1),Qinter_25,Qinter8(3),Qinter_med,valueMax)&
	  +IntegralOverQYpoint_S_Rec_8(var,process,incCut,CutParam,medQ,QMax_in,yMin_in,yMax_in,&
	  Qinter_med,Qinter8(3),Qinter_75,Qinter8(4),Qinter_fin,valueMax)
   else
    interX=valueACB
   end if
  end function IntegralOverQYpoint_S_Rec_8  
  
  
  
  
  !---------------------------------INTEGRATED over Y over Q 4---------------------------------------------------------------
  !!!! No need for check over Y they take a place within y-integration
    !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  !!!! taking into account minimum calls of y-integral we have  =81 points
  function Xsec_Qint_Yint(var,process,incCut,CutParam,QMin_in,QMax_in,ymin_in,ymax_in)
  real*8,dimension(1:8)::var
  logical,intent(in)::incCut
   real*8,dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real*8:: Xsec_Qint_Yint
   real*8 :: X1,X2,X3,X4,X5,X6,X7,X8,X9
   real*8 :: yMin_in,yMax_in,QMin_in,QMax_in
   real*8::valueMax,Q2,Q3,Q4,deltaQ,Q5,Q6,Q7,Q8
   real*8::value1,value2
    
    deltaQ=QMax_in-QMin_in
   Q2=QMin_in+deltaQ/8d0
   Q3=QMin_in+2d0*deltaQ/8d0
   Q4=QMin_in+3d0*deltaQ/8d0
   Q5=QMin_in+4d0*deltaQ/8d0
   Q6=QMin_in+5d0*deltaQ/8d0
   Q7=QMin_in+6d0*deltaQ/8d0
   Q8=QMin_in+7d0*deltaQ/8d0
   
   call SetQ(QMin_in,var)
   X1=2*QMin_in*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
   
   call SetQ(Q2,var)
   X2=2*Q2*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
   
   call SetQ(Q3,var)
   X3=2*Q3*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
   
   call SetQ(Q4,var)
   X4=2*Q4*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
   
   call SetQ(Q5,var)
   X5=2*Q5*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
   
   call SetQ(Q6,var)
   X6=2*Q6*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
   
   call SetQ(Q7,var)
   X7=2*Q7*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
   
   call SetQ(Q8,var)
   X8=2*Q8*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
   
   call SetQ(QMax_in,var)
   X9=2*QMax_in*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
   
      !!approximate integral value
    value1=2d0*(deltaQ/4d0)*(7d0*X1+32d0*X3+12d0*X5+32d0*X7+7d0*X9)/45d0
    value2=2d0*(deltaQ/8d0)*(7d0*X1+32d0*X2+12d0*X3+32d0*X4+14d0*X5+32d0*X6+&
   12d0*X7+32d0*X8+7d0*X9)/45d0
   valueMax=(2d0**6d0*value2-value1)/(2d0**6d0-1d0)
   
   Xsec_Qint_Yint=IntegralOverQYpoint_S_Rec(var,process,incCut,CutParam,QMin_in,Q5,yMin_in,yMax_in,&
   X1,X2,X3,X4,X5,valueMax)+IntegralOverQYpoint_S_Rec(var,process,incCut,CutParam,Q5,QMax_in,yMin_in,yMax_in,&
   X5,X6,X7,X8,X9,valueMax)
  end function Xsec_Qint_Yint
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function IntegralOverQYpoint_S_Rec(var,process,incCut,CutParam,&
			      QMin_in,QMax_in,yMin_in,yMax_in,X1,X3,X5,X7,X9,valueMax) result(interX)
   real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   real*8,dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real*8 :: interX,X1,X2,X3,X4,X5,X6,X7,X8,X9
   real*8 :: valueAB,valueACB,valueR
   real*8 :: yMin_in,yMax_in,QMin_in,QMax_in,Q2,Q3,Q4,deltaQ,Q5,Q6,Q7,Q8
   real*8,intent(in)::valueMax
   
   deltaQ=QMax_in-QMin_in
   Q2=QMin_in+deltaQ/8d0
   Q3=QMin_in+2d0*deltaQ/8d0
   Q4=QMin_in+3d0*deltaQ/8d0
   Q5=QMin_in+4d0*deltaQ/8d0
   Q6=QMin_in+5d0*deltaQ/8d0
   Q7=QMin_in+6d0*deltaQ/8d0
   Q8=QMax_in-deltaQ/8d0
   
   call SetQ(Q2,var)
   X2=2*Q2*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
      
   call SetQ(Q4,var)
   X4=2*Q4*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
   
   call SetQ(Q6,var)
   X6=2*Q6*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
   
   call SetQ(Q8,var)
   X8=2*Q8*Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
   
   valueAB=2d0*(deltaQ/4d0)*(7d0*X1+32d0*X3+12d0*X5+32d0*X7+7d0*X9)/45d0
   valueACB=2d0*(deltaQ/8d0)*(7d0*X1+32d0*X2+12d0*X3+32d0*X4+14d0*X5+&
   32d0*X6+12d0*X7+32d0*X8+7d0*X9)/45d0
   valueR=(2d0**6d0*valueACB-valueAB)/(2d0**6d0-1d0)
   
   If(ABS((valueACB-valueR)/valueMax)>tolerance) then
    interX=integralOverQYpoint_S_Rec(var,process,incCut,CutParam,QMin_in,Q5,yMin_in,yMax_in,&
    X1,X2,X3,X4,X5,valueMax)+integralOverQYpoint_S_Rec(var,process,incCut,CutParam,Q5,Qmax_in,&
    yMin_in,yMax_in,X5,X6,X7,X8,X9,valueMax)
   else
    interX=valueR
   end if
  end function IntegralOverQYpoint_S_Rec
  
  !---------------------------------INTEGRATED over Y over Q 4 diff in mT---------------------------------------------------------------
  !!!! No need for check over Y they take a place within y-integration
    !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  !!!! taking into account minimum calls of y-integral we have  =81 points
  function Xsec_Qint_Yint_mTint_1(var,process,incCut,CutParam,QMin_in,QMax_in,yMin_in,yMax_in,mT_min,mT_max)
    real*8,dimension(1:8)::var
    logical,intent(in)::incCut
    real*8,dimension(1:4),intent(in)::CutParam
    integer,dimension(1:3),intent(in)::process
    real*8:: Xsec_Qint_Yint_mTint_1
    real*8 :: X1,X2,X3,X4,X5
    real*8 :: yMin_in,yMax_in,QMin_in,QMax_in,mT_min,mT_max
    real*8::valueMax,Q2,Q3,Q4,deltaQ
    
    deltaQ=(QMax_in-QMin_in)/4d0
    Q2=QMin_in+deltaQ
    Q3=QMin_in+2d0*deltaQ
    Q4=QMax_in-deltaQ
   
    call SetQ(QMin_in,var)
    if (QMin_in.le.mT_min) then
    X1=0d0
    else if ((mT_min<QMin_in).and.(QMin_in.le.mT_max)) then
    X1=2*QMin_in*Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,&
    mT_min,QMin_in-0.000001d0)
    else if (mT_max<QMin_in) then
    X1=2*QMin_in*Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,mT_min,mT_max)
    end if

   
    call SetQ(Q2,var)
    if (Q2.le.mT_min) then
    X2=0d0
    else if ((mT_min<Q2).and.(Q2.le.mT_max)) then
    X2=2*Q2*Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,&
    mT_min,Q2-0.000001d0)
    else if (mT_max<Q2) then
    X2=2*Q2*Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,mT_min,mT_max)
    end if
    
   
    call SetQ(Q3,var)
    if (Q3.le.mT_min) then
    X3=0d0
    else if ((mT_min<Q3).and.(Q3.le.mT_max)) then
    X3=2*Q3*Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,&
    mT_min,Q3-0.000001d0)
    else if (mT_max<Q3) then
    X3=2*Q3*Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,mT_min,mT_max)
    end if
   
    call SetQ(Q4,var)
    if (Q4.le.mT_min) then
    X4=0d0
    else if ((mT_min<Q4).and.(Q4.le.mT_max)) then
    X4=2*Q4*Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,&
    mT_min,Q4-0.00001d0)
    else if (mT_max<Q4) then
    X4=2*Q4*Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,mT_min,mT_max)
    end if
   
    call SetQ(QMax_in,var)
    if (QMax_in.le.mT_min) then
    X5=0d0
    else if ((mT_min<QMax_in).and.(QMax_in.le.mT_max)) then
    X5=2*QMax_in*Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,&
    mT_min,QMax_in-0.000001d0)
    else if (mT_max<QMax_in) then
    X5=2*QMax_in*Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,mT_min,mT_max)
    end if
   

   
    !!approximate integral value
    valueMax=deltaQ*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/3d0
   
    Xsec_Qint_Yint_mTint_1=IntegralOverQYpointmTint_S_Rec_1(var,process,incCut,CutParam,QMin_in,Q3,yMin_in,yMax_in,&
    mT_min,mT_max,X1,X2,X3,valueMax)+IntegralOverQYpointmTint_S_Rec_1(var,process,incCut,CutParam,Q3,QMax_in,&
    yMin_in,yMax_in,mT_min,mT_max,X3,X4,X5,valueMax)
  end function Xsec_Qint_Yint_mTint_1
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function IntegralOverQYpointmTint_S_Rec_1(var,process,incCut,CutParam,&
			      QMin_in,QMax_in,yMin_in,yMax_in,mT_min,mT_max,X1,X3,X5,valueMax) result(interX)
   real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   real*8,dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real*8 :: interX,X1,X2,X3,X4,X5
   real*8 :: valueAB,valueACB
   real*8 :: yMin_in,yMax_in,QMin_in,QMax_in,Q2,Q3,Q4,deltaQ,mT_min,mT_max
   real*8,intent(in)::valueMax
   
    deltaQ=(QMax_in-QMin_in)/4d0
    Q2=QMin_in+deltaQ
    Q3=QMin_in+2d0*deltaQ
    Q4=QMax_in-deltaQ
   
    call SetQ(Q2,var)
    if (Q2.le.mT_min) then
    X2=0d0
    else if ((mT_min<Q2).and.(Q2.le.mT_max)) then
    X2=2*Q2*Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,&
    mT_min,Q2-0.000001d0)
    else if (mT_max<Q2) then
    X2=2*Q2*Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,mT_min,mT_max)
    end if
      
    call SetQ(Q4,var)
    if (Q4.le.mT_min) then
    X4=0d0
    else if ((mT_min<Q4).and.(Q4.le.mT_max)) then
    X4=2*Q4*Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,&
    mT_min,Q4-0.000001d0)
    else if (mT_max<Q4) then
    X4=2*Q4*Xsec_Yint_mTint_1(var,process,incCut,CutParam,ymin_in,ymax_in,mT_min,mT_max)
    end if
   
    valueAB=(QMax_in-QMin_in)*(X1+4d0*X3+X5)/6d0
    valueACB=deltaQ*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/3d0
 
    if (abs(valueMax).gt.0d0)then
    
    if(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=IntegralOverQYpointmTint_S_Rec_1(var,process,incCut,CutParam,QMin_in,Q3,yMin_in,yMax_in,&
    mT_min,mT_max,X1,X2,X3,valueMax)+IntegralOverQYpointmTint_S_Rec_1(var,process,incCut,CutParam,&
    Q3,QMax_in,yMin_in,yMax_in,mT_min,mT_max,X3,X4,X5,valueMax)
    else
    interX=valueACB
    end if
    
    else
    
    if(ABS(valueACB-valueAB)>0.000001d0) then
    interX=IntegralOverQYpointmTint_S_Rec_1(var,process,incCut,CutParam,QMin_in,Q3,yMin_in,yMax_in,&
    mT_min,mT_max,X1,X2,X3,valueMax)+IntegralOverQYpointmTint_S_Rec_1(var,process,incCut,CutParam,&
    Q3,QMax_in,yMin_in,yMax_in,mT_min,mT_max,X3,X4,X5,valueMax)
    else
    interX=valueACB
    end if
    
    end if
    
  end function IntegralOverQYpointmTint_S_Rec_1
  
  
  
  !---------------------------------INTEGRATED over Y complete over Q complete S16---------------------------------------------------------------
  !!!! No need for check over Y they take a place within y-integration
    !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  !!!! taking into account minimum calls of y-integral we have  =81 points
  function Xsec_Qcomplete_Ycomplete_S16(var,process,incCut,CutParam,Qmin_in,Qmax_in)
  real*8,dimension(1:8)::var
  logical,intent(in)::incCut
   !real*8,dimension(1:4),intent(in)::CutParam
   real*8,intent(in)::CutParam(:)
   integer,dimension(1:3),intent(in)::process
   real*8:: Xsec_Qcomplete_Ycomplete_S16
   real*8 :: QMin_in,QMax_in
   real*8::valueMax,deltaQ,medQ,Qvalue
   real*8,dimension(1:8)::Qinter8
   real*8::Qinter,Qinter_inc,Qinter_25,Qinter_med,Qinter_75,Qinter_fin
   real*8::Qinter_125,Qinter_375,Qinter_625,Qinter_875
   integer::i
    
    deltaQ=(QMax_in-QMin_in)/16d0
    medQ=Qmin_in+8d0*deltaQ

      !!approximate integral value
      
    call SetQ(QMin_in,var)
    Qinter_inc=2d0*QMin_in*Xsec_Yint(var,process,incCut,CutParam,&
    log(var(5)),-log(var(5)))
    Qinter=Qinter_inc
    

    
    do i=1,8
        Qvalue=QMin_in+(2d0*real(i)-1d0)*deltaQ
        call SetQ(Qvalue,var)
        Qinter8(i)=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        log(var(5)),-log(var(5)))
        Qinter=Qinter+4d0*Qinter8(i)
    end do
    
  
        Qvalue=QMin_in+2d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_125=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        log(var(5)),-log(var(5)))
        Qinter=Qinter+2d0*Qinter_125
        
        Qvalue=QMin_in+4d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_25=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        log(var(5)),-log(var(5)))
        Qinter=Qinter+2d0*Qinter_25
        
        Qvalue=QMin_in+6d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_375=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        log(var(5)),-log(var(5)))
        Qinter=Qinter+2d0*Qinter_375
        
        Qvalue=QMin_in+8d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_med=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        log(var(5)),-log(var(5)))
        Qinter=Qinter+2d0*Qinter_med
        
        Qvalue=QMin_in+10d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_625=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        log(var(5)),-log(var(5)))
        Qinter=Qinter+2d0*Qinter_625
        
        Qvalue=QMin_in+12d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_75=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        log(var(5)),-log(var(5)))
        Qinter=Qinter+2d0*Qinter_75
        
        Qvalue=QMin_in+14d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_875=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        log(var(5)),-log(var(5)))
        Qinter=Qinter+2d0*Qinter_875
        


    
    call SetQ(QMax_in,var)
    Qinter_fin=2d0*QMax_in*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
   
    valueMax=deltaQ*(Qinter+Qinter_fin)/3d0
   
   Xsec_Qcomplete_Ycomplete_S16=IntegralOverQYcomplete_S16_Rec(var,process,incCut,CutParam,QMin_in,medQ,&
   Qinter_inc,Qinter8(1),Qinter_125,Qinter8(2),Qinter_25,Qinter8(3),Qinter_375,Qinter8(4),Qinter_med,valueMax)+&
   IntegralOverQYcomplete_S16_Rec(var,process,incCut,CutParam,&
   medQ,QMax_in,Qinter_med,Qinter8(5),Qinter_625,Qinter8(6),Qinter_75,Qinter8(7),Qinter_875,Qinter8(8),&
   Qinter_fin,valueMax)
  end function Xsec_Qcomplete_Ycomplete_S16
  
  
  recursive function IntegralOverQYcomplete_S16_Rec(var,process,incCut,CutParam,&
			      QMin_in,QMax_in,Qinter_inc,Qinter_125,Qinter_25,Qinter_375,Qinter_med,Qinter_625,&
			      Qinter_75,Qinter_875,Qinter_fin,valueMax) result(interX)
   real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   !real*8,dimension(1:4),intent(in)::CutParam
   real*8,intent(in)::CutParam(:)
   integer,dimension(1:3),intent(in)::process
   real*8 :: interX
   real*8 :: valueAB,valueACB,val
   real*8 :: QMin_in,QMax_in,deltaQ,medQ,Qvalue
   real*8::Qinter,Qinter_inc,Qinter_25,Qinter_med,Qinter_75,Qinter_fin
   real*8,dimension(1:8)::Qinter8
   real*8,intent(in)::valueMax
   real*8::Qinter_125,Qinter_375,Qinter_625,Qinter_875
   integer::i
   
   deltaQ=(QMax_in-QMin_in)/16d0
   medQ=QMin_in+8d0*deltaQ
   
   
    Qinter=Qinter_inc+Qinter_fin
    
    do i=1,8
        Qvalue=QMin_in+(2d0*real(i)-1d0)*deltaQ
        call SetQ(Qvalue,var)
        Qinter8(i)=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        log(var(5)),-log(var(5)))
        Qinter=Qinter+4d0*Qinter8(i)
    end do
    
    
   
   valueAB=((QMax_in-QMin_in)/6d0)*(Qinter_inc+4d0*Qinter_med+Qinter_fin)
   valueACB=deltaQ*(Qinter+2d0*Qinter_125+2d0*Qinter_25+2d0*Qinter_375+2d0*Qinter_med+2d0*Qinter_625+&
   2d0*Qinter_75+2d0*Qinter_875)/3d0
   
  
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=IntegralOverQYcomplete_S16_Rec(var,process,incCut,CutParam,QMin_in,medQ,&
   Qinter_inc,Qinter8(1),Qinter_125,Qinter8(2),Qinter_25,Qinter8(3),Qinter_375,Qinter8(4),Qinter_med,valueMax)+&
   IntegralOverQYcomplete_S16_Rec(var,process,incCut,CutParam,medQ,QMax_in,&
   Qinter_med,Qinter8(5),Qinter_625,Qinter8(6),Qinter_75,Qinter8(7),Qinter_875,Qinter8(8),Qinter_fin,valueMax)
   else
    interX=valueACB
   end if
  end function IntegralOverQYcomplete_S16_Rec
  
  
  
  
  !---------------------------------INTEGRATED over Y complete over Q complete S8---------------------------------------------------------------
  !!!! No need for check over Y they take a place within y-integration
    !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  !!!! taking into account minimum calls of y-integral we have  =81 points
  function Xsec_Qcomplete_Ycomplete_S8(var,process,incCut,CutParam,Qmin_in,Qmax_in)
  real*8,dimension(1:8)::var
  logical,intent(in)::incCut
   !real*8,dimension(1:4),intent(in)::CutParam
   real*8,intent(in)::CutParam(:)
   integer,dimension(1:3),intent(in)::process
   real*8:: Xsec_Qcomplete_Ycomplete_S8
   real*8 :: QMin_in,QMax_in
   real*8::valueMax,deltaQ,medQ,Qvalue
   real*8,dimension(1:4)::Qinter8
   real*8::Qinter,Qinter_inc,Qinter_25,Qinter_med,Qinter_75,Qinter_fin
   integer::i
    
    deltaQ=(QMax_in-QMin_in)/8d0
    medQ=Qmin_in+4d0*deltaQ
   
      !!approximate integral value
      
    call SetQ(QMin_in,var)
    Qinter_inc=2d0*QMin_in*Xsec_Yint(var,process,incCut,CutParam,&
    log(var(5)),-log(var(5)))
    Qinter=Qinter_inc
    

    
    do i=1,4
        Qvalue=QMin_in+(2d0*real(i)-1d0)*deltaQ
        call SetQ(Qvalue,var)
        Qinter8(i)=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        log(var(5)),-log(var(5)))
        Qinter=Qinter+4d0*Qinter8(i)
    end do
    
        
        Qvalue=QMin_in+2d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_25=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        log(var(5)),-log(var(5)))
        Qinter=Qinter+2d0*Qinter_25
        

        
        Qvalue=QMin_in+4d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_med=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        log(var(5)),-log(var(5)))
        Qinter=Qinter+2d0*Qinter_med
        

        
        Qvalue=QMin_in+6d0*deltaQ
        call SetQ(Qvalue,var)
        Qinter_75=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        log(var(5)),-log(var(5)))
        Qinter=Qinter+2d0*Qinter_75
        

        


    
    call SetQ(QMax_in,var)
    Qinter_fin=2d0*QMax_in*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
   
    valueMax=deltaQ*(Qinter+Qinter_fin)/3d0
   
   Xsec_Qcomplete_Ycomplete_S8=IntegralOverQYcomplete_S8_Rec(var,process,incCut,CutParam,QMin_in,medQ,&
   Qinter_inc,Qinter8(1),Qinter_25,Qinter8(2),Qinter_med,valueMax)+&
   IntegralOverQYcomplete_S8_Rec(var,process,incCut,CutParam,&
   medQ,QMax_in,Qinter_med,Qinter8(3),Qinter_75,Qinter8(4),Qinter_fin,valueMax)
  end function Xsec_Qcomplete_Ycomplete_S8
  
  recursive function IntegralOverQYcomplete_S8_Rec(var,process,incCut,CutParam,&
			      QMin_in,QMax_in,Qinter_inc,Qinter_25,Qinter_med,&
			      Qinter_75,Qinter_fin,valueMax) result(interX)
   real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   !real*8,dimension(1:4),intent(in)::CutParam
   real*8,intent(in)::CutParam(:)
   integer,dimension(1:3),intent(in)::process
   real*8 :: interX
   real*8 :: valueAB,valueACB,val
   real*8 :: QMin_in,QMax_in,deltaQ,medQ,Qvalue
   real*8::Qinter,Qinter_inc,Qinter_25,Qinter_med,Qinter_75,Qinter_fin
   real*8,dimension(1:4)::Qinter8
   real*8,intent(in)::valueMax
   integer::i
   
   deltaQ=(QMax_in-QMin_in)/8d0
   medQ=QMin_in+4d0*deltaQ
   
    Qinter=Qinter_inc+Qinter_fin
    
    do i=1,4
        Qvalue=QMin_in+(2d0*real(i)-1d0)*deltaQ
        call SetQ(Qvalue,var)
        Qinter8(i)=2d0*Qvalue*Xsec_Yint(var,process,incCut,CutParam,&
        log(var(5)),-log(var(5)))
        Qinter=Qinter+4d0*Qinter8(i)
    end do
    
    
   
   valueAB=((QMax_in-QMin_in)/6d0)*(Qinter_inc+4d0*Qinter_med+Qinter_fin)
   valueACB=deltaQ*(Qinter+2d0*Qinter_25+2d0*Qinter_med+2d0*Qinter_75)/3d0
   
  
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=IntegralOverQYcomplete_S8_Rec(var,process,incCut,CutParam,QMin_in,medQ,&
   Qinter_inc,Qinter8(1),Qinter_25,Qinter8(2),Qinter_med,valueMax)+&
   IntegralOverQYcomplete_S8_Rec(var,process,incCut,CutParam,medQ,QMax_in,&
   Qinter_med,Qinter8(3),Qinter_75,Qinter8(4),Qinter_fin,valueMax)
   else
    interX=valueACB
   end if
  end function IntegralOverQYcomplete_S8_Rec
  
  
  
  !---------------------------------INTEGRATED over Y complete over Q complete S4---------------------------------------------------------------
  !!!! No need for check over Y they take a place within y-integration
    !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  !!!! taking into account minimum calls of y-integral we have  =81 points
  function Xsec_Qcomplete_Ycomplete_S4(var,process,incCut,CutParam,QMin_in,QMax_in)
  real*8,dimension(1:8)::var
  logical,intent(in)::incCut
   !real*8,dimension(1:4),intent(in)::CutParam
   real*8,intent(in)::CutParam(:)
   integer,dimension(1:3),intent(in)::process
   real*8:: Xsec_Qcomplete_Ycomplete_S4
   real*8 :: X1,X2,X3,X4,X5
   real*8 :: QMin_in,QMax_in
   real*8::valueMax,Q2,Q3,Q4,deltaQ
   integer::i
    
    deltaQ=(QMax_in-QMin_in)

   Q2=QMin_in+deltaQ/4d0
   Q3=QMin_in+deltaQ/2d0
   Q4=QMax_in-deltaQ/4d0
   
   call SetQ(QMin_in,var)
   X1=2*QMin_in*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
   
   call SetQ(Q2,var)
   X2=2*Q2*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
   
   call SetQ(Q3,var)
   X3=2*Q3*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
   
   call SetQ(Q4,var)
   X4=2*Q4*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
   
   call SetQ(QMax_in,var)
   X5=2*QMax_in*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
   
      !!approximate integral value
    
    valueMax=deltaQ*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   Xsec_Qcomplete_Ycomplete_S4=IntegralOverQYcomplete_S4_Rec(var,process,incCut,CutParam,QMin_in,Q3,&
   X1,X2,X3,valueMax)+IntegralOverQYcomplete_S4_Rec(var,process,incCut,CutParam,&
   Q3,QMax_in,X3,X4,X5,valueMax)
  end function Xsec_Qcomplete_Ycomplete_S4
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function IntegralOverQYcomplete_S4_Rec(var,process,incCut,CutParam,&
			      QMin_in,QMax_in,X1,X3,X5,valueMax) result(interX)
   real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   !real*8,dimension(1:4),intent(in)::CutParam
   real*8,intent(in)::CutParam(:)
   integer,dimension(1:3),intent(in)::process
   real*8 :: X1,X2,X3,X4,X5
   real*8 :: interX
   real*8 :: valueAB,valueACB,val
   real*8 :: QMin_in,QMax_in,Q2,Q3,Q4,deltaQ
   real*8,intent(in)::valueMax
   integer::i
   
   deltaQ=(QMax_in-QMin_in)
   
   Q2=QMin_in+deltaQ/4d0
   Q3=QMin_in+deltaQ/2d0
   Q4=QMax_in-deltaQ/4d0
   
   call SetQ(Q2,var)
   X2=2*Q2*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
      
   call SetQ(Q4,var)
   X4=2*Q4*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
    
   valueAB=deltaQ*(X1+4d0*X3+X5)/6d0
   valueACB=deltaQ*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
  
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=IntegralOverQYcomplete_S4_Rec(var,process,incCut,CutParam,QMin_in,Q3,&
   X1,X2,X3,valueMax)+&
   IntegralOverQYcomplete_S4_Rec(var,process,incCut,CutParam,Q3,QMax_in,&
   X3,X4,X5,valueMax)
   else
    interX=valueACB
   end if
  end function IntegralOverQYcomplete_S4_Rec
  
  
  
  
  !---------------------------------INTEGRATED over Y over Q over mT---------------------------------------------------------------
  !!!! No need for check over Y they take a place within y-integration
    !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  !!!! taking into account minimum calls of y-integral we have  =81 points
  function Xsec_Qint_Yint_mTint(var,process,incCut,CutParam,Qmin_in,Qmax_in,ymin_in,ymax_in,&
    mT_min,mT_max,NQ)
    real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   real*8,dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real*8:: Xsec_Qint_Yint_mTint
   real*8 :: yMin_in,yMax_in,QMin_in,QMax_in,mT_min,mT_max
   real*8::valueMax_mT,deltamT,med_mT
   integer,intent(in),optional::NQ
   real*8::mTinter,mTinter_inc,mTinter_fin,mTvalue,med_mY
   !real*8::mTinter_125,mTinter_25,mTinter_375,mTinter_med,mTinter_625,mTinter_75,mTinter_875
   !real*8,dimension(1:8)::mTinter8
   integer::i,j,Qbin_l
   real*8::Qbin_min,Qbin_max,Qbin_v
   real*8::X1,X2,X3,X4,X5,X6,X7,X8,X9
   real*8::mT2,mT3,mT4,mT5,mT6,mT7,mT8
   real*8::value1,value2
   

    
!      deltamT=(mT_max-mT_min)/16
!      med_mT=mT_min+8d0*deltamT
    deltamT=mT_max-mT_min
    mT2=mT_min+deltamT/8d0
    mT3=mT_min+2d0*deltamT/8d0
    mT4=mT_min+3d0*deltamT/8d0
    mT5=mT_min+4d0*deltamT/8d0
    mT6=mT_min+5d0*deltamT/8d0
    mT7=mT_min+6d0*deltamT/8d0
    mT8=mT_max-deltamT/8d0
    
    call SetmT(mT_min,var)
    if (Qmin_in>mT_min) then
        if ((Qmax_in-Qmin_in)/60d0-real(int((Qmax_in-Qmin_in)/60d0))>0d0) then
        X1=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT_min*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X1=X1+Qbin_v
        end do
        else
        X1=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT_min*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X1=X1+Qbin_v
        end do
        end if
    else
        if ((Qmax_in-mT_min-0.000001d0)/60d0-real(int((Qmax_in-mT_min-0.000001d0)/60d0))>0d0) then
        X1=0d0
        Qbin_l=int((Qmax_in-mT_min-0.000001d0)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=mT_min+0.000001d0+real(j-1)*(Qmax_in-mT_min-0.000001d0)/real(Qbin_l)
        Qbin_max=mT_min+0.000001d0+real(j)*(Qmax_in-mT_min-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT_min*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X1=X1+Qbin_v
        end do
        else
        X1=0d0
        Qbin_l=int((Qmax_in-mT_min-0.000001d0)/60d0)
        do j=1,Qbin_l
        Qbin_min=mT_min+0.000001d0+real(j-1)*(Qmax_in-mT_min-0.000001d0)/real(Qbin_l)
        Qbin_max=mT_min+0.000001d0+real(j)*(Qmax_in-mT_min-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT_min*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X1=X1+Qbin_v
        end do
        end if
    end if
    
    call SetmT(mT2,var)
    if (Qmin_in>mT2) then
        if ((Qmax_in-Qmin_in)/60d0-real(int((Qmax_in-Qmin_in)/60d0))>0d0) then
        X2=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT2*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X2=X2+Qbin_v
        end do
        else
        X2=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT2*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X2=X2+Qbin_v
        end do
        end if
    else
        if ((Qmax_in-mT2-0.000001d0)/60d0-real(int((Qmax_in-mT2-0.000001d0)/60d0))>0d0) then
        X2=0d0
        Qbin_l=int((Qmax_in-mT2-0.000001d0)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=mT2+0.000001d0+real(j-1)*(Qmax_in-mT2-0.000001d0)/real(Qbin_l)
        Qbin_max=mT2+0.000001d0+real(j)*(Qmax_in-mT2-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT2*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X2=X2+Qbin_v
        end do
        else
        X2=0d0
        Qbin_l=int((Qmax_in-mT2-0.000001d0)/60d0)
        do j=1,Qbin_l
        Qbin_min=mT2+0.000001d0+real(j-1)*(Qmax_in-mT2-0.000001d0)/real(Qbin_l)
        Qbin_max=mT2+0.000001d0+real(j)*(Qmax_in-mT2-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT2*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X2=X2+Qbin_v
        end do
        end if
    end if
    
    call SetmT(mT3,var)
    if (Qmin_in>mT3) then
        if ((Qmax_in-Qmin_in)/60d0-real(int((Qmax_in-Qmin_in)/60d0))>0d0) then
        X3=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT3*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X3=X3+Qbin_v
        end do
        else
        X3=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT3*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X3=X3+Qbin_v
        end do
        end if
    else
        if ((Qmax_in-mT3-0.000001d0)/60d0-real(int((Qmax_in-mT3-0.000001d0)/60d0))>0d0) then
        X3=0d0
        Qbin_l=int((Qmax_in-mT3-0.000001d0)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=mT3+0.000001d0+real(j-1)*(Qmax_in-mT3-0.000001d0)/real(Qbin_l)
        Qbin_max=mT3+0.000001d0+real(j)*(Qmax_in-mT3-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT3*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X3=X3+Qbin_v
        end do
        else
        X3=0d0
        Qbin_l=int((Qmax_in-mT3-0.000001d0)/60d0)
        do j=1,Qbin_l
        Qbin_min=mT3+0.000001d0+real(j-1)*(Qmax_in-mT3-0.000001d0)/real(Qbin_l)
        Qbin_max=mT3+0.000001d0+real(j)*(Qmax_in-mT3-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT3*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X3=X3+Qbin_v
        end do
        end if
    end if
    
    call SetmT(mT4,var)
    if (Qmin_in>mT4) then
        if ((Qmax_in-Qmin_in)/60d0-real(int((Qmax_in-Qmin_in)/60d0))>0d0) then
        X4=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT4*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X4=X4+Qbin_v
        end do
        else
        X4=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT4*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X4=X4+Qbin_v
        end do
        end if
    else
        if ((Qmax_in-mT4-0.000001d0)/60d0-real(int((Qmax_in-mT4-0.000001d0)/60d0))>0d0) then
        X4=0d0
        Qbin_l=int((Qmax_in-mT4-0.000001d0)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=mT4+0.000001d0+real(j-1)*(Qmax_in-mT4-0.000001d0)/real(Qbin_l)
        Qbin_max=mT4+0.000001d0+real(j)*(Qmax_in-mT4-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT4*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X4=X4+Qbin_v
        end do
        else
        X4=0d0
        Qbin_l=int((Qmax_in-mT4-0.000001d0)/60d0)
        do j=1,Qbin_l
        Qbin_min=mT4+0.000001d0+real(j-1)*(Qmax_in-mT4-0.000001d0)/real(Qbin_l)
        Qbin_max=mT4+0.000001d0+real(j)*(Qmax_in-mT4-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT4*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X4=X4+Qbin_v
        end do
        end if
    end if
    
    call SetmT(mT5,var)
    if (Qmin_in>mT5) then
        if ((Qmax_in-Qmin_in)/60d0-real(int((Qmax_in-Qmin_in)/60d0))>0d0) then
        X5=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT5*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X5=X5+Qbin_v
        end do
        else
        X5=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT5*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X5=X5+Qbin_v
        end do
        end if
    else
        if ((Qmax_in-mT5-0.000001d0)/60d0-real(int((Qmax_in-mT5-0.000001d0)/60d0))>0d0) then
        X5=0d0
        Qbin_l=int((Qmax_in-mT5-0.000001d0)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=mT5+0.000001d0+real(j-1)*(Qmax_in-mT5-0.000001d0)/real(Qbin_l)
        Qbin_max=mT5+0.000001d0+real(j)*(Qmax_in-mT5-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT5*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X5=X5+Qbin_v
        end do
        else
        X5=0d0
        Qbin_l=int((Qmax_in-mT5-0.000001d0)/60d0)
        do j=1,Qbin_l
        Qbin_min=mT5+0.000001d0+real(j-1)*(Qmax_in-mT5-0.000001d0)/real(Qbin_l)
        Qbin_max=mT5+0.000001d0+real(j)*(Qmax_in-mT5-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT5*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X5=X5+Qbin_v
        end do
        end if
    end if
    
    call SetmT(mT6,var)
    if (Qmin_in>mT6) then
        if ((Qmax_in-Qmin_in)/60d0-real(int((Qmax_in-Qmin_in)/60d0))>0d0) then
        X6=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT6*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X6=X6+Qbin_v
        end do
        else
        X6=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT6*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X6=X6+Qbin_v
        end do
        end if
    else
        if ((Qmax_in-mT6-0.000001d0)/60d0-real(int((Qmax_in-mT6-0.000001d0)/60d0))>0d0) then
        X6=0d0
        Qbin_l=int((Qmax_in-mT6-0.000001d0)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=mT6+0.000001d0+real(j-1)*(Qmax_in-mT6-0.000001d0)/real(Qbin_l)
        Qbin_max=mT6+0.000001d0+real(j)*(Qmax_in-mT6-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT6*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X6=X6+Qbin_v
        end do
        else
        X6=0d0
        Qbin_l=int((Qmax_in-mT6-0.000001d0)/60d0)
        do j=1,Qbin_l
        Qbin_min=mT6+0.000001d0+real(j-1)*(Qmax_in-mT6-0.000001d0)/real(Qbin_l)
        Qbin_max=mT6+0.000001d0+real(j)*(Qmax_in-mT6-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT6*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X6=X6+Qbin_v
        end do
        end if
    end if
    
    call SetmT(mT7,var)
    if (Qmin_in>mT7) then
        if ((Qmax_in-Qmin_in)/60d0-real(int((Qmax_in-Qmin_in)/60d0))>0d0) then
        X7=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT7*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X7=X7+Qbin_v
        end do
        else
        X7=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT7*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X7=X7+Qbin_v
        end do
        end if
    else
        if ((Qmax_in-mT7-0.000001d0)/60d0-real(int((Qmax_in-mT7-0.000001d0)/60d0))>0d0) then
        X7=0d0
        Qbin_l=int((Qmax_in-mT7-0.000001d0)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=mT7+0.000001d0+real(j-1)*(Qmax_in-mT7-0.000001d0)/real(Qbin_l)
        Qbin_max=mT7+0.000001d0+real(j)*(Qmax_in-mT7-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT7*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X7=X7+Qbin_v
        end do
        else
        X7=0d0
        Qbin_l=int((Qmax_in-mT7-0.000001d0)/60d0)
        do j=1,Qbin_l
        Qbin_min=mT7+0.000001d0+real(j-1)*(Qmax_in-mT7-0.000001d0)/real(Qbin_l)
        Qbin_max=mT7+0.000001d0+real(j)*(Qmax_in-mT7-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT7*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X7=X7+Qbin_v
        end do
        end if
    end if
    
    call SetmT(mT8,var)
    if (Qmin_in>mT8) then
        if ((Qmax_in-Qmin_in)/60d0-real(int((Qmax_in-Qmin_in)/60d0))>0d0) then
        X8=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT8*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X8=X8+Qbin_v
        end do
        else
        X8=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT8*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X8=X8+Qbin_v
        end do
        end if
    else
        if ((Qmax_in-mT8-0.000001d0)/60d0-real(int((Qmax_in-mT8-0.000001d0)/60d0))>0d0) then
        X8=0d0
        Qbin_l=int((Qmax_in-mT8-0.000001d0)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=mT8+0.000001d0+real(j-1)*(Qmax_in-mT8-0.000001d0)/real(Qbin_l)
        Qbin_max=mT8+0.000001d0+real(j)*(Qmax_in-mT8-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT8*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X8=X8+Qbin_v
        end do
        else
        X8=0d0
        Qbin_l=int((Qmax_in-mT8-0.000001d0)/60d0)
        do j=1,Qbin_l
        Qbin_min=mT8+0.000001d0+real(j-1)*(Qmax_in-mT8-0.000001d0)/real(Qbin_l)
        Qbin_max=mT8+0.000001d0+real(j)*(Qmax_in-mT8-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT8*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X8=X8+Qbin_v
        end do
        end if
    end if
    
    call SetmT(mT_max,var)
    if (Qmin_in>mT_max) then
        if ((Qmax_in-Qmin_in)/60d0-real(int((Qmax_in-Qmin_in)/60d0))>0d0) then
        X9=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT_max*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X9=X9+Qbin_v
        end do
        else
        X9=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT_max*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X9=X9+Qbin_v
        end do
        end if
    else
        if ((Qmax_in-mT_max-0.000001d0)/60d0-real(int((Qmax_in-mT_max-0.000001d0)/60d0))>0d0) then
        X9=0d0
        Qbin_l=int((Qmax_in-mT_max-0.000001d0)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=mT_max+0.000001d0+real(j-1)*(Qmax_in-mT_max-0.000001d0)/real(Qbin_l)
        Qbin_max=mT_max+0.000001d0+real(j)*(Qmax_in-mT_max-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT_max*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X9=X9+Qbin_v
        end do
        else
        X9=0d0
        Qbin_l=int((Qmax_in-mT_max-0.000001d0)/60d0)
        do j=1,Qbin_l
        Qbin_min=mT_max+0.000001d0+real(j-1)*(Qmax_in-mT_max-0.000001d0)/real(Qbin_l)
        Qbin_max=mT_max+0.000001d0+real(j)*(Qmax_in-mT_max-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT_max*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X9=X9+Qbin_v
        end do
        end if
    end if

    !We use Richardson's extrapolation
    value1=2d0*(deltamT/4d0)*(1d0*X1+32d0*X3+12d0*X5+32d0*X7+7d0*X9)/45d0
    value2=2d0*(deltamT/8d0)*(7d0*X1+32d0*X2+12d0*X3+32d0*X4+14d0*X5+32d0*X6+&
    12d0*X7+32d0*X8+7d0*X9)/45d0
    valueMax_mT=(2d0**6d0*value2-value1)/(2d0*6d0-1d0)
    
    Xsec_Qint_Yint_mTint=IntegralOverQYmTpoint_S_Rec(var,process,incCut,CutParam,QMin_in,QMax_in,yMin_in,yMax_in,&
    mT_min,mT5,X1,X2,X3,X4,X5,valueMax_mT,NQ)+&
    IntegralOverQYmTpoint_S_Rec(var,process,incCut,CutParam,QMin_in,QMax_in,yMin_in,yMax_in,&
    mT5,mT_max,X5,X6,X7,X8,X9,valueMax_mT,NQ)
   
  end function Xsec_Qint_Yint_mTint
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function IntegralOverQYmTpoint_S_Rec(var,process,incCut,CutParam,&
			      QMin_in,QMax_in,yMin_in,yMax_in,mT_min,mT_max,&
			      X1,X3,X5,X7,X9,valueMax_mT,NQ) result(interX_mT)
   real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   real*8,dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real*8 :: interX_mT
   real*8 :: valueAB,valueACB,valueR
   real*8 :: yMin_in,yMax_in,QMin_in,QMax_in,mT_min,mT_max,deltamT
   real*8,intent(in)::valueMax_mT
   integer,intent(in),optional::NQ
   real*8,dimension(1:8)::mTinter8
   real*8::mTinter_125,mTinter_25,mTinter_375,mTinter_med,mTinter_625,mTinter_75,mTinter_875
   real*8::med_mT,mTvalue,mTinter,mTinter_inc,mTinter_fin
   integer::i,j,Qbin_l
   real*8::Qbin_min,Qbin_max,Qbin_v
   real*8::X1,X2,X3,X4,X5,X6,X7,X8,X9
   real*8::mT2,mT3,mT4,mT5,mT6,mT7,mT8
   
    deltamT=mT_max-mT_min
    mT2=mT_min+deltamT/8d0
    mT3=mT_min+2d0*deltamT/8d0
    mT4=mT_min+3d0*deltamT/8d0
    mT5=mT_min+4d0*deltamT/8d0
    mT6=mT_min+5d0*deltamT/8d0
    mT7=mT_min+6d0*deltamT/8d0
    mT8=mT_max-deltamT/8d0
    
    call SetmT(mT2,var)
    if (Qmin_in>mT2) then
        if ((Qmax_in-Qmin_in)/60d0-real(int((Qmax_in-Qmin_in)/60d0))>0d0) then
        X2=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT2*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X2=X2+Qbin_v
        end do
        else
        X2=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT2*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X2=X2+Qbin_v
        end do
        end if
    else
        if ((Qmax_in-mT2-0.000001d0)/60d0-real(int((Qmax_in-mT2-0.000001d0)/60d0))>0d0) then
        X2=0d0
        Qbin_l=int((Qmax_in-mT2-0.000001d0)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=mT2+0.000001d0+real(j-1)*(Qmax_in-mT2-0.000001d0)/real(Qbin_l)
        Qbin_max=mT2+0.000001d0+real(j)*(Qmax_in-mT2-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT2*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X2=X2+Qbin_v
        end do
        else
        X2=0d0
        Qbin_l=int((Qmax_in-mT2-0.000001d0)/60d0)
        do j=1,Qbin_l
        Qbin_min=mT2+0.000001d0+real(j-1)*(Qmax_in-mT2-0.000001d0)/real(Qbin_l)
        Qbin_max=mT2+0.000001d0+real(j)*(Qmax_in-mT2-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT2*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X2=X2+Qbin_v
        end do
        end if
    end if
   
   
    call SetmT(mT4,var)
    if (Qmin_in>mT4) then
        if ((Qmax_in-Qmin_in)/60d0-real(int((Qmax_in-Qmin_in)/60d0))>0d0) then
        X4=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT4*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X4=X4+Qbin_v
        end do
        else
        X4=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT4*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X4=X4+Qbin_v
        end do
        end if
    else
        if ((Qmax_in-mT4-0.000001d0)/60d0-real(int((Qmax_in-mT4-0.000001d0)/60d0))>0d0) then
        X4=0d0
        Qbin_l=int((Qmax_in-mT4-0.000001d0)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=mT4+0.000001d0+real(j-1)*(Qmax_in-mT4-0.000001d0)/real(Qbin_l)
        Qbin_max=mT4+0.000001d0+real(j)*(Qmax_in-mT4-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT4*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X4=X4+Qbin_v
        end do
        else
        X4=0d0
        Qbin_l=int((Qmax_in-mT4-0.000001d0)/60d0)
        do j=1,Qbin_l
        Qbin_min=mT4+0.000001d0+real(j-1)*(Qmax_in-mT4-0.000001d0)/real(Qbin_l)
        Qbin_max=mT4+0.000001d0+real(j)*(Qmax_in-mT4-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT4*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X4=X4+Qbin_v
        end do
        end if
    end if
    
    
    call SetmT(mT6,var)
    if (Qmin_in>mT6) then
        if ((Qmax_in-Qmin_in)/60d0-real(int((Qmax_in-Qmin_in)/60d0))>0d0) then
        X6=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT6*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X6=X6+Qbin_v
        end do
        else
        X6=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT6*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X6=X6+Qbin_v
        end do
        end if
    else
        if ((Qmax_in-mT6-0.000001d0)/60d0-real(int((Qmax_in-mT6-0.000001d0)/60d0))>0d0) then
        X6=0d0
        Qbin_l=int((Qmax_in-mT6-0.000001d0)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=mT6+0.000001d0+real(j-1)*(Qmax_in-mT6-0.000001d0)/real(Qbin_l)
        Qbin_max=mT6+0.000001d0+real(j)*(Qmax_in-mT6-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT6*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X6=X6+Qbin_v
        end do
        else
        X6=0d0
        Qbin_l=int((Qmax_in-mT6-0.000001d0)/60d0)
        do j=1,Qbin_l
        Qbin_min=mT6+0.000001d0+real(j-1)*(Qmax_in-mT6-0.000001d0)/real(Qbin_l)
        Qbin_max=mT6+0.000001d0+real(j)*(Qmax_in-mT6-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT6*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X6=X6+Qbin_v
        end do
        end if
    end if
    
    call SetmT(mT8,var)
    if (Qmin_in>mT8) then
        if ((Qmax_in-Qmin_in)/60d0-real(int((Qmax_in-Qmin_in)/60d0))>0d0) then
        X8=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT8*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X8=X8+Qbin_v
        end do
        else
        X8=0d0
        Qbin_l=int((Qmax_in-Qmin_in)/60d0)
        do j=1,Qbin_l
        Qbin_min=Qmin_in+real(j-1)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_max=Qmin_in+real(j)*(Qmax_in-Qmin_in)/real(Qbin_l)
        Qbin_v=2*mT8*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X8=X8+Qbin_v
        end do
        end if
    else
        if ((Qmax_in-mT8-0.000001d0)/60d0-real(int((Qmax_in-mT8-0.000001d0)/60d0))>0d0) then
        X8=0d0
        Qbin_l=int((Qmax_in-mT8-0.000001d0)/60d0)+1
        do j=1,Qbin_l
        Qbin_min=mT8+0.000001d0+real(j-1)*(Qmax_in-mT8-0.000001d0)/real(Qbin_l)
        Qbin_max=mT8+0.000001d0+real(j)*(Qmax_in-mT8-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT8*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X8=X8+Qbin_v
        end do
        else
        X4=0d0
        Qbin_l=int((Qmax_in-mT8-0.000001d0)/60d0)
        do j=1,Qbin_l
        Qbin_min=mT8+0.000001d0+real(j-1)*(Qmax_in-mT8-0.000001d0)/real(Qbin_l)
        Qbin_max=mT8+0.000001d0+real(j)*(Qmax_in-mT8-0.000001d0)/real(Qbin_l)
        Qbin_v=2*mT8*Xsec_Qint_Yint(var,process,incCut,CutParam,&
        Qbin_min,Qbin_max,yMin_in,yMax_in) 
        X8=X8+Qbin_v
        end do
        end if
    end if
   
   
    valueAB=2d0*(deltamT/4d0)*(7d0*X1+32d0*X3+12d0*X5+32d0*X7+7d0*X9)/45d0
    valueACB=2d0*(deltamT/8d0)*(7d0*X1+32d0*X2+12d0*X3+32d0*X4+14d0*X5+&
    32d0*X6+12d0*X7+32d0*X8+7d0*X9)/45d0
    valueR=(2d0**6d0*valueACB-valueAB)/(2d0**6d0-1d0)
   
   If(ABS((valueACB-valueR)/valueMax_mT)>tolerance) then
    interX_mT=integralOverQYmTpoint_S_Rec(var,process,incCut,CutParam,QMin_in,QMax_in,yMin_in,yMax_in,mT_min,mT5,&
    X1,X2,X3,X4,X5,valueMax_mT,NQ)+integralOverQYmTpoint_S_Rec(var,process,incCut,CutParam,Qmin_in,Qmax_in,yMin_in,yMax_in,&
    mT5,mT_max,X5,X6,X7,X8,X9,valueMax_mT,NQ)
   else
    interX_mT=valueR
   end if
   
  end function IntegralOverQYmTpoint_S_Rec
  
  
  
  !---------------------------------INTEGRATED over Y complete over Q complete over mT---------------------------------------------------------------
  !!!! No need for check over Y they take a place within y-integration
    !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  !!!! taking into account minimum calls of y-integral we have  =81 points
  function Xsec_Qcomplete_Ycomplete_mTint(var,process,incCut,CutParam,QMax_in,&
    mT_min,mT_max,n)
    real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   !real*8,dimension(1:4),intent(in)::CutParam
   real*8,intent(in)::CutParam(:)
   integer,dimension(1:3),intent(in)::process
   real*8:: Xsec_Qcomplete_Ycomplete_mTint
   real*8 ::X1,X2,X3,X4,X5,mT2,mT3,mT4
   real*8 :: QMax_in,mT_min,mT_max
   real*8::valueMax_mT,deltamT,medmT
   integer::i
   integer,intent(in)::n
    
   deltamT=(mT_max-mT_min)
   mT2=mT_min+deltamT/4d0
   mT3=mT_min+deltamT/2d0
   mT4=mT_max-deltamT/4d0
   
   call SetmT(mT_min,var)
   if(n==16) then
   X1=2*mT_min*Xsec_Qcomplete_Ycomplete_S16(var,process,incCut,CutParam,mT_min+0.000001d0,Qmax_in)
   else if (n==8) then
   X1=2*mT_min*Xsec_Qcomplete_Ycomplete_S8(var,process,incCut,CutParam,mT_min+0.000001d0,Qmax_in)
   else
   X1=2*mT_min*Xsec_Qcomplete_Ycomplete_S4(var,process,incCut,CutParam,mT_min+0.000001d0,Qmax_in)
   end if
   
   call SetmT(mT2,var)
   if (n==16) then
   X2=2*mT2*Xsec_Qcomplete_Ycomplete_S16(var,process,incCut,CutParam,mT2+0.000001d0,Qmax_in)
   else if (n==8) then
   X2=2*mT2*Xsec_Qcomplete_Ycomplete_S8(var,process,incCut,CutParam,mT2+0.000001d0,Qmax_in)
   else
   X2=2*mT2*Xsec_Qcomplete_Ycomplete_S4(var,process,incCut,CutParam,mT2+0.000001d0,Qmax_in)
   end if
   
   call SetmT(mT3,var)
   if (n==16) then
   X3=2*mT3*Xsec_Qcomplete_Ycomplete_S16(var,process,incCut,CutParam,mT3+0.000001d0,Qmax_in)
   else if (n==8) then
   X3=2*mT3*Xsec_Qcomplete_Ycomplete_S8(var,process,incCut,CutParam,mT3+0.000001d0,Qmax_in)
   else 
   X3=2*mT3*Xsec_Qcomplete_Ycomplete_S4(var,process,incCut,CutParam,mT3+0.000001d0,Qmax_in)
   end if
   
   call SetmT(mT4,var)
   if (n==16) then
   X4=2*mT4*Xsec_Qcomplete_Ycomplete_S16(var,process,incCut,CutParam,mT4+0.000001d0,Qmax_in)
   else if (n==8) then
   X4=2*mT4*Xsec_Qcomplete_Ycomplete_S8(var,process,incCut,CutParam,mT4+0.000001d0,Qmax_in)
   else
   X4=2*mT4*Xsec_Qcomplete_Ycomplete_S4(var,process,incCut,CutParam,mT4+0.000001d0,Qmax_in)
   end if
   
   call SetmT(mT_max,var)
   if (n==16) then
   X5=2*mT_max*Xsec_Qcomplete_Ycomplete_S16(var,process,incCut,CutParam,mT_max+0.000001d0,Qmax_in)
   else if (n==8) then
   X5=2*mT_max*Xsec_Qcomplete_Ycomplete_S8(var,process,incCut,CutParam,mT_max+0.000001d0,Qmax_in)
   else
   X5=2*mT_max*Xsec_Qcomplete_Ycomplete_S4(var,process,incCut,CutParam,mT_max+0.000001d0,Qmax_in)
   end if
   

   
    
      !!approximate integral value
   valueMax_mT=deltamT*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   Xsec_Qcomplete_Ycomplete_mTint=IntegralOverQYcompletemTpoint_S_Rec(var,process,incCut,CutParam,QMax_in,mT_min,mT3,&
   X1,X2,X3,valueMax_mT,n)+IntegralOverQYcompletemTpoint_S_Rec(var,process,incCut,&
   CutParam,QMax_in,mT3,mT_max,X3,X4,X5,valueMax_mT,n)
  end function Xsec_Qcomplete_Ycomplete_mTint
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function IntegralOverQYcompletemTpoint_S_Rec(var,process,incCut,CutParam,&
			      QMax_in,mT_min,mT_max,X1,X3,X5,&
			      valueMax_mT,n) result(interX_mT)
   real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   !real*8,dimension(1:4),intent(in)::CutParam
   real*8,intent(in)::CutParam(:)
   integer,dimension(1:3),intent(in)::process
   real*8 :: interX_mT,X1,X2,X3,X4,X5
   real*8 :: valueAB,valueACB
   real*8 :: QMax_in,mT_min,mT_max,mT2,mT3,mT4,deltamT
   integer:: i
   real*8,intent(in)::valueMax_mT
   integer,intent(in)::n
   
    deltamT=(mT_max-mT_min)
    !medmT=mT_min+deltamT*4d0
   
   mT2=mT_min+deltamT/4d0
   mT3=mT_min+deltamT/2d0
   mT4=mT_max-deltamT/4d0
   
   call SetmT(mT2,var)
   if (n==16) then
   X2=2*mT2*Xsec_Qcomplete_Ycomplete_S16(var,process,incCut,CutParam,mT2+0.000001d0,Qmax_in)
   else if (n==8) then
   X2=2*mT2*Xsec_Qcomplete_Ycomplete_S8(var,process,incCut,CutParam,mT2+0.000001d0,Qmax_in)
   else 
   X2=2*mT2*Xsec_Qcomplete_Ycomplete_S4(var,process,incCut,CutParam,mT2+0.000001d0,Qmax_in)
   end if
      
   call SetmT(mT4,var)
   if (n==16) then
   X4=2*mT4*Xsec_Qcomplete_Ycomplete_S16(var,process,incCut,CutParam,mT4+0.000001d0,Qmax_in)
   else if (n==8) then
   X4=2*mT4*Xsec_Qcomplete_Ycomplete_S8(var,process,incCut,CutParam,mT4+0.000001d0,Qmax_in)
   else
   X4=2*mT4*Xsec_Qcomplete_Ycomplete_S4(var,process,incCut,CutParam,mT4+0.000001d0,Qmax_in)
   end if
   
   
   
   
   valueAB=(deltamT/6d0)*(X1+4d0*X3+X5)
   valueACB=deltamT*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   If(ABS((valueACB-valueAB)/valueMax_mT)>tolerance) then
    interX_mT=IntegralOverQYcompletemTpoint_S_Rec(var,process,incCut,CutParam,Qmax_in,mT_min,mT3,&
    X1,X2,X3,valueMax_mT,n)+&
    IntegralOverQYcompletemTpoint_S_Rec(var,process,incCut,CutParam,Qmax_in,&
    mT3,mT_max,X3,X4,X5,valueMax_mT,n)
   else
    interX_mT=valueACB
   end if
  end function IntegralOverQYcompletemTpoint_S_Rec
  
  
  
  !---------------------------------INTEGRATED over Y complete over mT in NRW---------------------------------------------------------------
  !!!! No need for check over Y they take a place within y-integration
    !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  !!!! taking into account minimum calls of y-integral we have  =81 points
  function Xsec_Ycomplete_mTint_NRW(var,process,incCut,CutParam,&
    mT_min,mT_max)
    real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   real*8,intent(in)::CutParam(:)
   !real*8,dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real*8:: Xsec_Ycomplete_mTint_NRW
   real*8 :: X1,X2,X3,X4,X5
   real*8 :: QMax_in,mT_min,mT_max
   real*8::valueMax_mT,mT2,mT3,mT4,deltamT
    
    deltamT=mT_max-mT_min
   mT2=mT_min+deltamT/4d0
   mT3=mT_min+deltamT/2d0
   mT4=mT_max-deltamT/4d0
   
   call SetmT(mT_min,var)
   X1=2*mT_min*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
   
   call SetmT(mT2,var)
   X2=2*mT2*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
   
   call SetmT(mT3,var)
   X3=2*mT3*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
   
   call SetmT(mT4,var)
   X4=2*mT4*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
   
   call SetmT(mT_max,var)
   X5=2*mT_max*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
   
      !!approximate integral value
   valueMax_mT=deltamT*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   Xsec_Ycomplete_mTint_NRW=IntegralOverYcompletemTpoint_NRW_S_Rec(var,process,incCut,CutParam,mT_min,mT3,&
   X1,X2,X3,valueMax_mT)+IntegralOverYcompletemTpoint_NRW_S_Rec(var,process,incCut,CutParam,&
   mT3,mT_max,X3,X4,X5,valueMax_mT)
  end function Xsec_Ycomplete_mTint_NRW
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function IntegralOverYcompletemTpoint_NRW_S_Rec(var,process,incCut,CutParam,&
			      mT_min,mT_max,X1,X3,X5,valueMax_mT) result(interX_mT)
   real*8,dimension(1:8)::var
   logical,intent(in)::incCut
   real*8,intent(in)::CutParam(:)
   !real*8,dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real*8 :: interX_mT,X1,X2,X3,X4,X5
   real*8 :: valueAB,valueACB
   real*8 :: mT_min,mT_max,mT2,mT3,mT4,deltamT
   real*8,intent(in)::valueMax_mT
   
   deltamT=mT_max-mT_min
   mT2=mT_min+deltamT/4d0
   mT3=mT_min+deltamT/2d0
   mT4=mT_max-deltamT/4d0
   
   call SetmT(mT2,var)
   X2=2*mT2*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
      
   call SetmT(mT4,var)
   X4=2*mT4*Xsec_Yint(var,process,incCut,CutParam,log(var(5)),-log(var(5)))
   
   valueAB=deltamT*(X1+4d0*X3+X5)/6d0
   valueACB=deltamT*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   If(ABS((valueACB-valueAB)/valueMax_mT)>tolerance) then
    interX_mT=IntegralOverYcompletemTpoint_NRW_S_Rec(var,process,incCut,CutParam,mT_min,mT3,&
    X1,X2,X3,valueMax_mT)+IntegralOverYcompletemTpoint_NRW_S_Rec(var,process,incCut,CutParam,&
    mT3,mT_max,X3,X4,X5,valueMax_mT)
   else
    interX_mT=valueACB
   end if
  end function IntegralOverYcompletemTpoint_NRW_S_Rec
  
  
  
  !---------------------------------INTEGRATED over Y over Q over pT-------------------------------------------------------------
  !!!integration over PT is made by Num-sections
  !!!N even  
  function Xsec_PTint_Qint_Yint(process,incCut,CutParam,s_in,qt_min,qt_max,Q_min,Q_max,ymin_in,ymax_in,mT_in,Num)
    real*8,dimension(1:8)::var
    logical,intent(in)::incCut
    real*8,dimension(1:4),intent(in)::CutParam
    integer,dimension(1:3),intent(in)::process
    real*8:: Xsec_PTint_Qint_Yint,X0,Xfin
    real*8:: ymin_in,ymax_in,Q_min,Q_max,qt_min,qt_max,s_in,mT_in
    integer :: Num
    
    if(qt_min<1d-3) then
      var=kinematicArray(1d-3,s_in,(Q_min+Q_max)/2d0,(ymin_in+ymax_in)/2d0,mT_in)
    else
      var=kinematicArray(qt_min,s_in,(Q_min+Q_max)/2d0,(ymin_in+ymax_in)/2d0,mT_in)
    end if
    
    X0=2d0*qt_min*Xsec_Qint_Yint(var,process,incCut,CutParam,Q_min,Q_max,ymin_in,ymax_in)
    
    call Xsec_PTint_Qint_Yint_0(process,incCut,CutParam,s_in,qt_min,qt_max,Q_min,Q_max,ymin_in,ymax_in,mT_in,Num,Xfin,X0)
    Xsec_PTint_Qint_Yint=Xfin
    
  end function Xsec_PTint_Qint_Yint
  
  
  !!!integration over PT is made by Num-sections
  !!!N even
  !!! X0 is value of the function at qt_min input
  !!! X0 is value of the function at qt_max output
  !!! !!! Xfin is value of the cross-section
  subroutine Xsec_PTint_Qint_Yint_0(process,incCut,CutParam,s_in,qt_min,qt_max,Q_min,Q_max,ymin_in,ymax_in,mT_in,Num,Xfin,X0)
    real*8,dimension(1:8)::var
    logical,intent(in)::incCut
    real*8,dimension(1:4),intent(in)::CutParam
    integer,dimension(1:3),intent(in)::process
    real*8:: Xfin,X0
    real*8:: ymin_in,ymax_in,Q_min,Q_max,qt_min,qt_max,s_in,mT_in
    integer :: i,Num
    
    real*8::deltaQT,qT_cur,inter
    
    if(mod(num,2)>0) then 
      write(*,*) 'ERROR: arTeMiDe_DY: number of Simpson sections is odd. Evaluation stop.'
      stop
    end if
    
    deltaQT=(qt_max-qt_min)/Num
    inter=X0!!!first term is calculated eqarlier
    
    var=kinematicArray(qt_min,s_in,(Q_min+Q_max)/2d0,(ymin_in+ymax_in)/2d0,mT_in)
    
    !!!! even terms
    do i=1,Num-1,2
    qT_cur=qt_min+i*deltaQT
    call SetQT(qT_cur,var)
    inter=inter+8d0*qt_cur*Xsec_Qint_Yint(var,process,incCut,CutParam,Q_min,Q_max,ymin_in,ymax_in)
    end do
    
    if(Num>2) then
    !!!! odd terms
    do i=2,Num-2,2
    qT_cur=qt_min+i*deltaQT
    call SetQT(qT_cur,var)
    inter=inter+4d0*qt_cur*Xsec_Qint_Yint(var,process,incCut,CutParam,Q_min,Q_max,ymin_in,ymax_in)
    end do
    end if
    
    call SetQT(qT_max,var)
    X0=2d0*qt_max*Xsec_Qint_Yint(var,process,incCut,CutParam,Q_min,Q_max,ymin_in,ymax_in)!!!! last term
    inter=inter+X0  
    
    Xfin=deltaQT/3d0*inter
    
  end subroutine Xsec_PTint_Qint_Yint_0
  
  
  
  
  !---------------------------------INTEGRATED over Y over Q over pT over mT-------------------------------------------------------------
  !!!integration over PT is made by Num-sections
  !!!N even  
  function Xsec_PTint_Qint_Yint_mTint(process,incCut,CutParam,s_in,qt_min,qt_max,Q_min,Q_max,ymin_in,ymax_in,&
    mT_min,mT_max,Num,NQ)
    real*8,dimension(1:8)::var
    logical,intent(in)::incCut
    real*8,dimension(1:4),intent(in)::CutParam
    integer,dimension(1:3),intent(in)::process
    real*8:: Xsec_PTint_Qint_Yint_mTint,X0_mT,Xfin_mT
    real*8:: ymin_in,ymax_in,Q_min,Q_max,qt_min,qt_max,s_in,mT_min,mT_max
    integer,intent(in) :: Num
    integer,intent(in),optional::NQ
    real*8::Q_v,Q_b_max,Q_b_min
    integer::Q_l,j
    
    if(qt_min<1d-3) then
      var=kinematicArray(1d-3,s_in,(Q_min+Q_max)/2d0,(ymin_in+ymax_in)/2d0,(mT_min+mT_max)/2d0)
    else
      var=kinematicArray(qt_min,s_in,(Q_min+Q_max)/2d0,(ymin_in+ymax_in)/2d0,(mT_min+mT_max)/2d0)
    end if
    
    
    if ((Q_max-Q_min)/60d0-real(int((Q_max-Q_min)/60d0))>0d0) then
        X0_mT=0d0
        Q_l=int((Q_max-Q_min)/60d0)+1
        do j=1,Q_l
        Q_b_min=Q_min+real(j-1)*(Q_max-Q_min)/real(Q_l)
        Q_b_max=Q_min+real(j)*(Q_max-Q_min)/real(Q_l)
        Q_v=2d0*qt_min*Xsec_Qint_Yint_mTint_1(var,process,incCut,CutParam,&
        Q_b_min,Q_b_max,ymin_in,ymax_in,mT_min,mT_max)
        X0_mT=X0_mT+Q_v
        end do
    else
        X0_mT=0d0
        Q_l=int((Q_max-Q_min)/60d0)
        do j=1,Q_l
        Q_b_min=Q_min+real(j-1)*(Q_max-Q_min)/real(Q_l)
        Q_b_max=Q_min+real(j)*(Q_max-Q_min)/real(Q_l)
        Q_v=2d0*qt_min*Xsec_Qint_Yint_mTint_1(var,process,incCut,CutParam,&
        Q_b_min,Q_b_max,ymin_in,ymax_in,mT_min,mT_max)
        X0_mT=X0_mT+Q_v
        end do
    end if
    
    call Xsec_PTint_Qint_Yint_mTint_0(process,incCut,CutParam,s_in,qt_min,qt_max,Q_min,Q_max,ymin_in,ymax_in,&
    mT_min,mT_max,Num,Xfin_mT,X0_mT,NQ)
    Xsec_PTint_Qint_Yint_mTint=Xfin_mT
    
  end function Xsec_PTint_Qint_Yint_mTint
  
  
  !!!integration over PT is made by Num-sections
  !!!N even
  !!! X0 is value of the function at qt_min input
  !!! X0 is value of the function at qt_max output
  !!! !!! Xfin is value of the cross-section
  subroutine Xsec_PTint_Qint_Yint_mTint_0(process,incCut,CutParam,s_in,qt_min,qt_max,Q_min,Q_max,ymin_in,ymax_in,&
    mT_min,mT_max,Num,Xfin_mT,X0_mT,NQ)
    real*8,dimension(1:8)::var
    logical,intent(in)::incCut
    real*8,dimension(1:4),intent(in)::CutParam
    integer,dimension(1:3),intent(in)::process
    real*8:: Xfin_mT,X0_mT
    real*8:: ymin_in,ymax_in,Q_min,Q_max,qt_min,qt_max,s_in,mT_min,mT_max
    integer :: i,Num
    integer,intent(in),optional::NQ
    real*8::Q_v,Q_b_max,Q_b_min
    integer::Q_l,j
    
    real*8::deltaQT,qT_cur,inter_mT
    
    if(mod(num,2)>0) then 
      write(*,*) 'ERROR: arTeMiDe_DY: number of Simpson sections is odd. Evaluation stop.'
      stop
    end if
    
    deltaQT=(qt_max-qt_min)/Num
    inter_mT=X0_mT!!!first term is calculated eqarlier
    
    var=kinematicArray(qt_min,s_in,(Q_min+Q_max)/2d0,(ymin_in+ymax_in)/2d0,(mT_min+mT_max)/2d0)
 
    !!!! even terms
    do i=1,Num-1,2
    qT_cur=qt_min+i*deltaQT
    call SetQT(qT_cur,var)
    if ((Q_max-Q_min)/60d0-real(int((Q_max-Q_min)/60d0))>0d0) then
        Q_l=int((Q_max-Q_min)/60d0)+1
        do j=1,Q_l
        Q_b_min=Q_min+real(j-1)*(Q_max-Q_min)/real(Q_l)
        Q_b_max=Q_min+real(j)*(Q_max-Q_min)/real(Q_l)
        inter_mT=inter_mT+8d0*qT_cur*Xsec_Qint_Yint_mTint_1(var,process,incCut,CutParam,&
        Q_b_min,Q_b_max,ymin_in,ymax_in,mT_min,mT_max)
        end do
    else
        Q_l=int((Q_max-Q_min)/60d0)
        do j=1,Q_l
        Q_b_min=Q_min+real(j-1)*(Q_max-Q_min)/real(Q_l)
        Q_b_max=Q_min+real(j)*(Q_max-Q_min)/real(Q_l)
        inter_mT=inter_mT+8d0*qT_cur*Xsec_Qint_Yint_mTint_1(var,process,incCut,CutParam,&
        Q_b_min,Q_b_max,ymin_in,ymax_in,mT_min,mT_max)
        end do
    end if
    end do
    
    if(Num>2) then
    !!!! odd terms
    do i=2,Num-2,2
    qT_cur=qt_min+i*deltaQT
    call SetQT(qT_cur,var)

    if ((Q_max-Q_min)/60d0-real(int((Q_max-Q_min)/60d0))>0d0) then
        Q_l=int((Q_max-Q_min)/60d0)+1
        do j=1,Q_l
        Q_b_min=Q_min+real(j-1)*(Q_max-Q_min)/real(Q_l)
        Q_b_max=Q_min+real(j)*(Q_max-Q_min)/real(Q_l)
        inter_mT=inter_mT+4d0*qt_cur*Xsec_Qint_Yint_mTint_1(var,process,incCut,CutParam,&
        Q_b_min,Q_b_max,ymin_in,ymax_in,mT_min,mT_max)
        end do
    else
        Q_l=int((Q_max-Q_min)/60d0)
        do j=1,Q_l
        Q_b_min=Q_min+real(j-1)*(Q_max-Q_min)/real(Q_l)
        Q_b_max=Q_min+real(j)*(Q_max-Q_min)/real(Q_l)
        inter_mT=inter_mT+4d0*qt_cur*Xsec_Qint_Yint_mTint_1(var,process,incCut,CutParam,&
        Q_b_min,Q_b_max,ymin_in,ymax_in,mT_min,mT_max)
        end do
    end if 
    
    end do
    end if
    
    call SetQT(qT_max,var)
    if ((Q_max-Q_min)/60d0-real(int((Q_max-Q_min)/60d0))>0d0) then
        X0_mT=0d0
        Q_l=int((mT_max-mT_min)/60d0)+1
        do j=1,Q_l
        Q_b_min=Q_min+real(j-1)*(Q_max-Q_min)/real(Q_l)
        Q_b_max=Q_min+real(j)*(Q_max-Q_min)/real(Q_l)
        Q_v=2d0*qt_max*Xsec_Qint_Yint_mTint_1(var,process,incCut,CutParam,&
        Q_b_min,Q_b_max,ymin_in,ymax_in,mT_min,mT_max)!!!! last term
        X0_mT=X0_mT+Q_v
        end do
    else
        X0_mT=0d0
        Q_l=int((Q_max-Q_min)/60d0)
        do j=1,Q_l
        Q_b_min=Q_min+real(j-1)*(Q_max-Q_min)/real(Q_l)
        Q_b_max=Q_min+real(j)*(Q_max-Q_min)/real(Q_l)
        Q_v=2d0*qt_max*Xsec_Qint_Yint_mTint_1(var,process,incCut,CutParam,&
        Q_b_min,Q_b_max,ymin_in,ymax_in,mT_min,mT_max)!!!! last term
        X0_mT=X0_mT+Q_v
        end do
    end if

    inter_mT=inter_mT+X0_mT 
    
    Xfin_mT=(deltaQT/3d0)*inter_mT
    
  end subroutine Xsec_PTint_Qint_Yint_mTint_0
  
  
  
  !---------------------------------INTEGRATED over Y complete over Q complete over pT over mT-------------------------------------------------------------
  !!!integration over PT is made by Num-sections
  !!!N even  
  function Xsec_PTint_Qcomplete_Ycomplete_mTint(process,incCut,CutParam,s_in,qT_min,qT_max,Q_max,y_in,&
    mT_min,mT_max,Num,n)
    real*8,dimension(1:8)::var
    logical,intent(in)::incCut
    !real*8,dimension(1:4),intent(in)::CutParam
    real*8,intent(in)::CutParam(:)
    integer,dimension(1:3),intent(in)::process
    real*8:: Xsec_PTint_Qcomplete_Ycomplete_mTint,X0_mT,Xfin_mT
    real*8:: y_in,Q_max,qT_min,qT_max,s_in,mT_min,mT_max
    integer :: Num
    integer,intent(in)::n
    
    if(qt_min<1d-3) then
      var=kinematicArray(1d-3,s_in,(mT_min+Q_max)/2d0,y_in,(mT_min+mT_max)/2d0)
    else
      var=kinematicArray(qt_min,s_in,(mT_min+Q_max)/2d0,y_in,(mT_min+mT_max)/2d0)
    end if
    
    X0_mT=2d0*qt_min*Xsec_Qcomplete_Ycomplete_mTint(var,process,incCut,CutParam,Q_max,mT_min,mT_max,n)
    
    call Xsec_PTint_Qcomplete_Ycomplete_mTint_0(process,incCut,CutParam,s_in,qT_min,qT_max,Q_max,y_in,&
    mT_min,mT_max,Num,n,Xfin_mT,X0_mT)
    Xsec_PTint_Qcomplete_Ycomplete_mTint=Xfin_mT
    
  end function Xsec_PTint_Qcomplete_Ycomplete_mTint
  
  
  !!!integration over PT is made by Num-sections
  !!!N even
  !!! X0 is value of the function at qt_min input
  !!! X0 is value of the function at qt_max output
  !!! !!! Xfin is value of the cross-section
  subroutine Xsec_PTint_Qcomplete_Ycomplete_mTint_0(process,incCut,CutParam,s_in,qT_min,qT_max,Q_max,y_in,&
    mT_min,mT_max,Num,n,Xfin_mT,X0_mT)
    real*8,dimension(1:8)::var
    logical,intent(in)::incCut
    !real*8,dimension(1:4),intent(in)::CutParam
    real*8,intent(in)::CutParam(:)
    integer,dimension(1:3),intent(in)::process
    real*8:: Xfin_mT,X0_mT
    real*8:: y_in,Q_max,qT_min,qT_max,s_in,mT_min,mT_max
    integer :: i,Num
    integer,intent(in)::n
    
    real*8::deltaQT,qT_cur,inter_mT
    
    if(mod(num,2)>0) then 
      write(*,*) 'ERROR: arTeMiDe_DY: number of Simpson sections is odd. Evaluation stop.'
      stop
    end if
    
    deltaQT=(qT_max-qT_min)/Num
    inter_mT=X0_mT!!!first term is calculated eqarlier
    
    var=kinematicArray(qt_min,s_in,(mT_min+Q_max)/2d0,y_in,(mT_min+mT_max)/2d0)
    
    !!!! even terms
    do i=1,Num-1,2
    qT_cur=qT_min+real(i)*deltaQT
    call SetQT(qT_cur,var)
    inter_mT=inter_mT+8d0*qT_cur*Xsec_Qcomplete_Ycomplete_mTint(var,process,incCut,CutParam,Q_max,&
    mT_min,mT_max,n)
    end do
    
    if(Num>2) then
    !!!! odd terms
    do i=2,Num-2,2
    qT_cur=qT_min+real(i)*deltaQT
    call SetQT(qT_cur,var)
    inter_mT=inter_mT+4d0*qT_cur*Xsec_Qcomplete_Ycomplete_mTint(var,process,incCut,CutParam,Q_max,&
    mT_min,mT_max,n)
    end do
    end if
    
    call SetQT(qT_max,var)
    X0_mT=2d0*qT_max*Xsec_Qcomplete_Ycomplete_mTint(var,process,incCut,CutParam,Q_max,&
    mT_min,mT_max,n)!!!! last term
    inter_mT=inter_mT+X0_mT 
    
    Xfin_mT=(deltaQT/3d0)*inter_mT
    
  end subroutine Xsec_PTint_Qcomplete_Ycomplete_mTint_0
  
  
  
  !---------------------------------INTEGRATED over Y complete Q=Mw complete over pT over mT-------------------------------------------------------------
  !!!integration over PT is made by Num-sections
  !!!N even  
  function Xsec_PTint_Ycomplete_mTint_NRW(process,incCut,CutParam,s_in,qt_min,qt_max,Q,y_in,&
    mT_min,Num)
    real*8,dimension(1:8)::var
    logical,intent(in)::incCut
    real*8,intent(in)::CutParam(:)
    !real*8,dimension(1:4),intent(in)::CutParam
    integer,dimension(1:3),intent(in)::process
    real*8:: Xsec_PTint_Ycomplete_mTint_NRW,X0_mT,Xfin_mT
    real*8:: y_in,Q,qt_min,qt_max,s_in,mT_min
    integer :: Num
    
    if(qt_min<1d-3) then
      var=kinematicArray(1d-3,s_in,Q,y_in,(mT_min+Q)/2d0)
    else
      var=kinematicArray(qt_min,s_in,Q,y_in,(mT_min+Q)/2d0)
    end if
    
    X0_mT=2d0*qt_min*Xsec_Ycomplete_mTint_NRW(var,process,incCut,CutParam,mT_min,Q)
    
    call Xsec_PTint_Ycomplete_mTint_NRW_0(process,incCut,CutParam,s_in,qt_min,qt_max,Q,y_in,&
    mT_min,Num,Xfin_mT,X0_mT)
    Xsec_PTint_Ycomplete_mTint_NRW=Xfin_mT
    
  end function Xsec_PTint_Ycomplete_mTint_NRW
  
  
  !!!integration over PT is made by Num-sections
  !!!N even
  !!! X0 is value of the function at qt_min input
  !!! X0 is value of the function at qt_max output
  !!! !!! Xfin is value of the cross-section
  subroutine Xsec_PTint_Ycomplete_mTint_NRW_0(process,incCut,CutParam,s_in,qt_min,qt_max,Q,y_in,&
    mT_min,Num,Xfin_mT,X0_mT)
    real*8,dimension(1:8)::var
    logical,intent(in)::incCut
    real*8,intent(in)::CutParam(:)
    !real*8,dimension(1:4),intent(in)::CutParam
    integer,dimension(1:3),intent(in)::process
    real*8:: Xfin_mT,X0_mT
    real*8:: y_in,Q,qt_min,qt_max,s_in,mT_min
    integer :: i,Num
    
    real*8::deltaQT,qT_cur,inter_mT
    
    if(mod(num,2)>0) then 
      write(*,*) 'ERROR: arTeMiDe_DY: number of Simpson sections is odd. Evaluation stop.'
      stop
    end if
    
    deltaQT=(qt_max-qt_min)/Num
    inter_mT=X0_mT!!!first term is calculated eqarlier
    
    var=kinematicArray(qt_min,s_in,Q,y_in,(mT_min+Q)/2d0)
    
    !!!! even terms
    do i=1,Num-1,2
    qT_cur=qt_min+i*deltaQT
    call SetQT(qT_cur,var)
    inter_mT=inter_mT+8d0*qT_cur*Xsec_Ycomplete_mTint_NRW(var,process,incCut,CutParam,&
    mT_min,Q)
    end do
    
    if(Num>2) then
    !!!! odd terms
    do i=2,Num-2,2
    qT_cur=qt_min+i*deltaQT
    call SetQT(qT_cur,var)
    inter_mT=inter_mT+4d0*qT_cur*Xsec_Ycomplete_mTint_NRW(var,process,incCut,CutParam,&
    mT_min,Q)
    end do
    end if
    
    call SetQT(qT_max,var)
    X0_mT=2d0*qt_max*Xsec_Ycomplete_mTint_NRW(var,process,incCut,CutParam,&
    mT_min,Q)!!!! last term
    inter_mT=inter_mT+X0_mT 
    
    Xfin_mT=deltaQT/3d0*inter_mT
    
  end subroutine Xsec_PTint_Ycomplete_mTint_NRW_0
  
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!INTERFACES TO CALCULATING CROSS-SECTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !---------------------------------INTEGRATED------------------------------------------------------------------
  
  !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList(X_list,qt_List)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    integer :: i,length
    real*8,dimension(1:8)::var
    length=size(qt_list)
    CallCounter=CallCounter+length
     do i=1,length
       var=kinematicArray(qt_List(i),s_global,Q_global,y_global,mT_global)
       X_List(i)=PreFactor1(process_global(1))*xSec(var,process_global,includeCuts_global,CutParameters_global)
     end do
  end subroutine xSecList
  
  !!!!Evaluate differential xSec at single point
  subroutine xSecSingle(X,qT_in)
   real*8,dimension(1:8)::var
   real*8:: X,qT_in
   CallCounter=CallCounter+1
   var=kinematicArray(qt_in,s_global,Q_global,y_global,mT_global)
   X=PreFactor1(process_global(1))*xSec(var,process_global,includeCuts_global,CutParameters_global)
  end subroutine xSecSingle
  
  
!---------------------------------INTEGRATED over Y---------------------------------------------------------------
    !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_Yint(X_list,qt_List,yMin_in,yMax_in)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::yMin_in,yMax_in
    real*8,dimension(1:8)::var
    integer :: i,length
    length=size(qt_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(var)
     do i=1,length
	var=kinematicArray(qt_List(i),s_global,Q_global,y_global,mT_global)
       X_List(i)=PreFactor1(process_global(1))*Xsec_Yint(var,process_global,includeCuts_global,CutParameters_global,yMin_in,yMax_in)
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecList_Yint
  
  !!
  subroutine xSecSingle_Yint(X,qt,yMin_in,yMax_in)
    real*8::X,qT
    real*8::yMin_in,yMax_in
    real*8,dimension(1:8)::var
    
   CallCounter=CallCounter+1
   var=kinematicArray(qt,s_global,Q_global,y_global,mT_global)
   X=PreFactor1(process_global(1))*Xsec_Yint(var,process_global,includeCuts_global,CutParameters_global,yMin_in,yMax_in)
  end subroutine xSecSingle_Yint
  
      !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_Ycomplete(X_list,qt_List)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    integer :: i,length
    real*8,dimension(1:8)::var
    length=size(qt_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(var)
     do i=1,length
       var=kinematicArray(qt_List(i),s_global,Q_global,y_global,mT_global)
       X_List(i)=PreFactor1(process_global(1))*Xsec_Yint(var,process_global,includeCuts_global,CutParameters_global,&
		    log(var(5)),-log(var(5)))
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecList_Ycomplete
  
  subroutine xSecSingle_Ycomplete(X,qt)
    real*8::X,qT
    real*8,dimension(1:8)::var
    
   CallCounter=CallCounter+1
   var=kinematicArray(qt,s_global,Q_global,y_global,mT_global)
   X=PreFactor1(process_global(1))*Xsec_Yint(var,process_global,includeCuts_global,CutParameters_global,log(var(5)),-log(var(5)))
  end subroutine xSecSingle_Ycomplete
  
 
  !---------------------------------INTEGRATED over Q---------------------------------------------------------------
  !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_Qint(X_list,qt_List,Q_min,Q_max)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max
    integer :: i,length
    real*8,dimension(1:8)::var
    length=size(qt_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(var)
     do i=1,length
       var=kinematicArray(qt_List(i),s_global,Q_global,y_global,mT_global)
       X_List(i)=PreFactor1(process_global(1))*Xsec_Qint(var,process_global,includeCuts_global,CutParameters_global,Q_min,Q_max)
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecList_Qint
  
  !!
  subroutine xSecSingle_Qint(X,qt,Q_min,Q_max)
    real*8::X,qT
    real*8::Q_min,Q_max
    real*8,dimension(1:8)::var
    
   CallCounter=CallCounter+1
   var=kinematicArray(qt,s_global,Q_global,y_global,mT_global)
   X=PreFactor1(process_global(1))*Xsec_Qint(var,process_global,includeCuts_global,CutParameters_global,Q_min,Q_max)
  end subroutine xSecSingle_Qint
  
  
!---------------------------------INTEGRATED over Y over Q---------------------------------------------------------------
  !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_Qint_Yint(X_list,qt_List,Q_min,Q_max,yMin_in,yMax_in)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::yMin_in,yMax_in,Q_min,Q_max
    integer :: i,length
    real*8,dimension(1:8)::var
    length=size(qt_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(var)
     do i=1,length
       var=kinematicArray(qt_List(i),s_global,Q_global,y_global,mT_global)
       X_List(i)=PreFactor1(process_global(1))*Xsec_Qint_Yint(var,process_global,includeCuts_global,CutParameters_global,&
		Q_min,Q_max,yMin_in,yMax_in)
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecList_Qint_Yint
  
  !!
  subroutine xSecSingle_Qint_Yint(X,qt,Q_min,Q_max,yMin_in,yMax_in)
    real*8::X,qT
    real*8::yMin_in,yMax_in,Q_min,Q_max
    real*8,dimension(1:8)::var
   CallCounter=CallCounter+1
   var=kinematicArray(qt,s_global,Q_global,y_global,mT_global)
   X=PreFactor1(process_global(1))*Xsec_Qint_Yint(var,process_global,includeCuts_global,CutParameters_global,&
	Q_min,Q_max,yMin_in,yMax_in)
  end subroutine xSecSingle_Qint_Yint
  
      !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_Qint_Ycomplete(X_list,qt_List,Q_min,Q_max)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max
    integer :: i,length
    real*8,dimension(1:8)::var
    length=size(qt_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(var)
     do i=1,length
       var=kinematicArray(qt_List(i),s_global,Q_global,y_global,mT_global)
       X_List(i)=PreFactor1(process_global(1))*Xsec_Qint_Yint(var,process_global,includeCuts_global,CutParameters_global,&
			      Q_min,Q_max,log(var(5)),-log(var(5)))
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecList_Qint_Ycomplete
  
  subroutine xSecSingle_Qint_Ycomplete(X,qt,Q_min,Q_max)
      real*8::X,qT
      real*8::Q_min,Q_max
      real*8,dimension(1:8)::var
    
   CallCounter=CallCounter+1
   var=kinematicArray(qt,s_global,Q_global,y_global,mT_global)
   X=PreFactor1(process_global(1))*Xsec_Qint_Yint(var,process_global,includeCuts_global,CutParameters_global,&
	      Q_min,Q_max,log(var(5)),-log(var(5)))
  end subroutine xSecSingle_Qint_Ycomplete
  
  !---------------------------------INTEGRATED over Y over Q  over PT----------------------------------------------------------
  !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_PTintN_Qint_Yint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_in,yMax_in,num)
    real*8, intent(in) :: qtMIN_list(:),qtMAX_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::yMin_in,yMax_in,Q_min,Q_max
    integer :: i,length,num
    length=size(qtMIN_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED)
     do i=1,length
       X_List(i)=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint(process_global,includeCuts_global,CutParameters_global,&
				s_global,qtMIN_List(i),qtMAX_list(i),Q_min,Q_max,yMin_in,yMax_in,mT_global,num)
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecList_PTintN_Qint_Yint
  
  subroutine xSecList_PTint_Qint_Yint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_in,yMax_in)
    real*8, intent(in) :: qtMIN_list(:),qtMAX_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::yMin_in,yMax_in,Q_min,Q_max
    
    call xSecList_PTintN_Qint_Yint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_in,yMax_in,NumPTdefault)

  end subroutine xSecList_PTint_Qint_Yint
  
  subroutine xSecListList_PTintN_Qint_Yint(X_list,qt_List,Q_min,Q_max,yMin_in,yMax_in,num)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max,X0,Xfin,yMin_in,yMax_in
    integer :: i,length,length2,num
    real*8,dimension(1:8)::var
    length2=size(qt_list)
    length=size(X_list)
    if( (length2-length) .ne. 1) then    
      write(*,*) 'ERROR: arTeMiDe_DY: pt integration : sizes of lists (for cross-seciont and pt-bins) are inconsistent'
      write(*,*) 'Evaluation stop'
      stop
    end if 
    
    !-----------------parallel version-------------
    !$ call xSecList_PTintN_Qint_Yint(X_list,qt_List(1:length2-1),qt_list(2:length2),Q_min,Q_max,yMin_in,yMax_in,num)
    !$ return
    
    !-----------------single-thread version-------------
    CallCounter=CallCounter+length
    
    var=kinematicArray(qt_list(1),s_global,Q_global,y_global,mT_global)
    X0=2d0*qt_list(1)*Xsec_Qint_Yint(var,process_global,includeCuts_global,CutParameters_global,Q_min,Q_max,yMin_in,yMax_in)
    do i=1,length
       call Xsec_PTint_Qint_Yint_0(process_global,includeCuts_global,CutParameters_global,&
		    s_global,qt_list(i),qt_list(i+1),Q_min,Q_max,yMin_in,yMax_in,mT_global,Num,Xfin,X0)
       X_List(i)=PreFactor1(process_global(1))*Xfin
     end do
  end subroutine xSecListList_PTintN_Qint_Yint
  
  subroutine xSecListList_PTint_Qint_Yint(X_list,qt_List,Q_min,Q_max,yMin_in,yMax_in)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max,yMin_in,yMax_in
    
    call xSecListList_PTintN_Qint_Yint(X_list,qt_List,Q_min,Q_max,yMin_in,yMax_in,NumPTdefault)
    
  end subroutine xSecListList_PTint_Qint_Yint
  
  !!
  subroutine xSecSingle_PTintN_Qint_Yint(X,qt_Min,qt_Max,Q_min,Q_max,yMin_in,yMax_in,num)
    real*8::X,qt_Min,qt_Max
    real*8::yMin_in,yMax_in,Q_min,Q_max
    integer::num
    
   CallCounter=CallCounter+1
   X=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint(process_global,includeCuts_global,CutParameters_global,&
				  s_global,qt_Min,qt_Max,Q_min,Q_max,yMin_in,yMax_in,mT_global,num)
  end subroutine xSecSingle_PTintN_Qint_Yint
  
  subroutine xSecSingle_PTint_Qint_Yint(X,qt_Min,qt_Max,Q_min,Q_max,yMin_in,yMax_in)
    real*8::X,qt_Min,qt_Max
    real*8::yMin_in,yMax_in,Q_min,Q_max
   X=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint(process_global,includeCuts_global,CutParameters_global,&
			    s_global, qt_Min,qt_Max,Q_min,Q_max,yMin_in,yMax_in,mT_global,NumPTdefault)
  end subroutine xSecSingle_PTint_Qint_Yint
  
  subroutine xSecListPY_PTint_Qint_Yint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_List,yMax_List,num)
    real*8, intent(in) :: qtMIN_list(:),qtMAX_list(:),yMin_List(:),yMax_List(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max
    integer :: i,length,num
    length=size(qtMIN_list)
    if(size(qtMAX_list)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: pt integration (2) : sizes of lists (pt min-max bins) are inconsistent'
      write(*,*) 'Evaluation stop'
      stop
    end if
    
    if(size(yMin_List)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: pt integration (3) : sizes of lists (pt-bins vs yMin) are inconsistent'
      write(*,*) 'Evaluation stop'
      stop
    end if
    
    if(size(yMax_List)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: pt integration (4) : sizes of lists (pt-bins vs yMax) are inconsistent'
      write(*,*) 'Evaluation stop'
      stop
    end if
    
    
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED)
     do i=1,length
       X_List(i)=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint(process_global,includeCuts_global,CutParameters_global,&
			  s_global,qtMIN_List(i),qtMAX_list(i),Q_min,Q_max,yMin_List(i),yMax_List(i),mT_global,num)
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecListPY_PTint_Qint_Yint
  
  subroutine xSecListPY_PTintN_Qint_Yint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_List,yMax_List)
    real*8, intent(in) :: qtMIN_list(:),qtMAX_list(:),yMin_List(:),yMax_List(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max
    call xSecListPY_PTint_Qint_Yint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_List,yMax_List,NumPTdefault)
  end subroutine xSecListPY_PTintN_Qint_Yint
  
  !---------------------------------INTEGRATED over Y (complete) over Q  over PT----------------------------------------------------------!
  
  !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  !! I set y in (-1000,1000) since the check is made in the integration routine
  subroutine xSecList_PTintN_Qint_Ycomplete(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,num)
    real*8, intent(in) :: qtMIN_list(:),qtMAX_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max
    integer :: i,length,num
    length=size(qtMIN_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED)
     do i=1,length
       X_List(i)=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint(process_global,includeCuts_global,CutParameters_global,&
				s_global,qtMIN_list(i),qtMAX_list(i),Q_min,Q_max,-1000d0,1000d0,mT_global,num)
     end do
    !$OMP END PARALLEL DO
  end subroutine xSecList_PTintN_Qint_Ycomplete
  
  subroutine xSecList_PTint_Qint_Ycomplete(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max)
    real*8, intent(in) :: qtMIN_list(:),qtMAX_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max
    
    call xSecList_PTintN_Qint_Ycomplete(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,NumPTdefault)
  end subroutine xSecList_PTint_Qint_Ycomplete
  
  subroutine xSecSingle_PTintN_Qint_Ycomplete(X,qt_min,qt_max,Q_min,Q_max,num)
      real*8::X,qT_min,qT_max
      real*8::Q_min,Q_max
      integer::num
    
   CallCounter=CallCounter+1
   X=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint(process_global,includeCuts_global,CutParameters_global,&
			    s_global,qt_min,qt_max,Q_min,Q_max,-1000d0,1000d0,mT_global,num)
  end subroutine xSecSingle_PTintN_Qint_Ycomplete
  
  
  subroutine xSecSingle_PTint_Qint_Ycomplete(X,qt_min,qt_max,Q_min,Q_max)
      real*8::X,qT_min,qT_max
      real*8::Q_min,Q_max
    
   CallCounter=CallCounter+1
   X=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint(process_global,includeCuts_global,CutParameters_global,&
			    s_global,qt_min,qt_max,Q_min,Q_max,-1000d0,1000d0,mT_global,NumPTdefault)
  end subroutine xSecSingle_PTint_Qint_Ycomplete
  
  
  subroutine xSecListList_PTintN_Qint_Ycomplete(X_list,qt_List,Q_min,Q_max,num)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max,X0,Xfin
    integer :: i,length,length2,num
    real*8,dimension(1:8)::var
    
    length2=size(qt_list)
    length=size(X_list)
    if( (length2-length) .ne. 1) then    
      write(*,*) 'ERROR: arTeMiDe_DY: pt integration : sizes of lists (for cross-seciont and pt-bins) are inconsistent'
      write(*,*) 'Evaluation stop'
      stop
    end if  
    
    !-----------------parallel version-------------
    !$ call xSecList_PTintN_Qint_Ycomplete(X_list,qt_List(1:length2-1),qt_list(2:length2),Q_min,Q_max,num)
    !$ return
    
    !-----------------single-thread version-------------
    CallCounter=CallCounter+length
    
    var=kinematicArray(qt_list(1),s_global,Q_global,y_global,mT_global)
    X0=2d0*qt_list(1)*Xsec_Qint_Yint(var,process_global,includeCuts_global,CutParameters_global,Q_min,Q_max,-1000d0,1000d0)
    do i=1,length
       call Xsec_PTint_Qint_Yint_0(process_global,includeCuts_global,CutParameters_global,&
		      s_global, qt_list(i),qt_list(i+1),Q_min,Q_max,-1000d0,1000d0,mT_global,Num,Xfin,X0)
       X_List(i)=PreFactor1(process_global(1))*Xfin
     end do
  end subroutine xSecListList_PTintN_Qint_Ycomplete
  
  
  
  subroutine xSecListList_PTint_Qint_Ycomplete(X_list,qt_List,Q_min,Q_max)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max
    
    call xSecListList_PTintN_Qint_Ycomplete(X_list,qt_List,Q_min,Q_max,NumPTdefault)
    
  end subroutine xSecListList_PTint_Qint_Ycomplete
  
  
  
  !---------------------------------INTEGRATED over Y over Q  over PT over mT----------------------------------------------------------
  !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_PTintN_Qint_Yint_mTint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_in,yMax_in,&
  mT_min,mT_max,num)
    real*8, intent(in) :: qtMIN_list(:),qtMAX_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::yMin_in,yMax_in,Q_min,Q_max, mT_min,mT_max
    integer :: i,length,num
    length=size(qtMIN_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED)
     do i=1,length
       X_List(i)=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint_mTint(process_global,includeCuts_global,CutParameters_global,&
				s_global,qtMIN_List(i),qtMAX_list(i),Q_min,Q_max,yMin_in,yMax_in,mT_min,mT_max,num)
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecList_PTintN_Qint_Yint_mTint
  
  subroutine xSecList_PTint_Qint_Yint_mTint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_in,yMax_in)
    real*8, intent(in) :: qtMIN_list(:),qtMAX_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::yMin_in,yMax_in,Q_min,Q_max,mT_min,mT_max
    
    call xSecList_PTintN_Qint_Yint_mTint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_in,yMax_in,&
    mT_min,mT_max,NumPTdefault)

  end subroutine xSecList_PTint_Qint_Yint_mTint
  
  subroutine xSecListList_PTintN_Qint_Yint_mTint(X_list,qt_List,Q_min,Q_max,yMin_in,yMax_in,&
  mT_min,mT_max,num)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max,X0,Xfin,yMin_in,yMax_in,mT_min,mT_max
    integer :: i,length,length2,num
    real*8,dimension(1:8)::var
    length2=size(qt_list)
    length=size(X_list)
    if( (length2-length) .ne. 1) then    
      write(*,*) 'ERROR: arTeMiDe_DY: pt integration : sizes of lists (for cross-seciont and pt-bins) are inconsistent'
      write(*,*) 'Evaluation stop'
      stop
    end if 
    
    !-----------------parallel version-------------
    !$ call xSecList_PTintN_Qint_Yint(X_list,qt_List(1:length2-1),qt_list(2:length2),Q_min,Q_max,yMin_in,yMax_in,num)
    !$ return
    
    !-----------------single-thread version-------------
    CallCounter=CallCounter+length
    
    var=kinematicArray(qt_list(1),s_global,Q_global,y_global,mT_global)
    X0=2d0*qt_list(1)*Xsec_Qint_Yint_mTint(var,process_global,includeCuts_global,CutParameters_global,Q_min,Q_max,yMin_in,&
    yMax_in,mT_min,mT_max)
    do i=1,length
       call Xsec_PTint_Qint_Yint_mTint_0(process_global,includeCuts_global,CutParameters_global,&
		    s_global,qt_list(i),qt_list(i+1),Q_min,Q_max,yMin_in,yMax_in,mT_min,mT_max,Num,Xfin,X0)
       X_List(i)=PreFactor1(process_global(1))*Xfin
     end do
  end subroutine xSecListList_PTintN_Qint_Yint_mTint
  
  subroutine xSecListList_PTint_Qint_Yint_mTint(X_list,qt_List,Q_min,Q_max,yMin_in,yMax_in,mT_min,mT_max)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max,yMin_in,yMax_in,mT_min,mT_max
    
    call xSecListList_PTintN_Qint_Yint_mTint(X_list,qt_List,Q_min,Q_max,yMin_in,yMax_in,&
    mT_min,mT_max,NumPTdefault)
    
  end subroutine xSecListList_PTint_Qint_Yint_mTint
  
  !!
  subroutine xSecSingle_PTintN_Qint_Yint_mTint(X,qt_Min,qt_Max,Q_min,Q_max,yMin_in,yMax_in,&
  mT_min,mT_max,num)
    real*8::X,qt_Min,qt_Max
    real*8::yMin_in,yMax_in,Q_min,Q_max,mT_min,mT_max
    integer::num
    
   CallCounter=CallCounter+1
   X=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint_mTint(process_global,includeCuts_global,CutParameters_global,&
				  s_global,qt_Min,qt_Max,Q_min,Q_max,yMin_in,yMax_in,mT_min,mT_max,num)
  end subroutine xSecSingle_PTintN_Qint_Yint_mTint
  
  subroutine xSecSingle_PTint_Qint_Yint_mTint(X,qt_Min,qt_Max,Q_min,Q_max,yMin_in,yMax_in,mT_min,mT_max)
    real*8::X,qt_Min,qt_Max
    real*8::yMin_in,yMax_in,Q_min,Q_max,mT_min,mT_max
   X=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint_mTint(process_global,includeCuts_global,CutParameters_global,&
			    s_global, qt_Min,qt_Max,Q_min,Q_max,yMin_in,yMax_in,mT_min,mT_max,NumPTdefault)
  end subroutine xSecSingle_PTint_Qint_Yint_mTint
  
  subroutine xSecListPY_PTint_Qint_Yint_mTint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_List,yMax_List,&
  mT_min,mT_max,num)
    real*8, intent(in) :: qtMIN_list(:),qtMAX_list(:),yMin_List(:),yMax_List(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max,mT_min,mT_max
    integer :: i,length,num
    length=size(qtMIN_list)
    if(size(qtMAX_list)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: pt integration (2) : sizes of lists (pt min-max bins) are inconsistent'
      write(*,*) 'Evaluation stop'
      stop
    end if
    
    if(size(yMin_List)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: pt integration (3) : sizes of lists (pt-bins vs yMin) are inconsistent'
      write(*,*) 'Evaluation stop'
      stop
    end if
    
    if(size(yMax_List)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: pt integration (4) : sizes of lists (pt-bins vs yMax) are inconsistent'
      write(*,*) 'Evaluation stop'
      stop
    end if
    
    
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED)
     do i=1,length
       X_List(i)=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint_mTint(process_global,includeCuts_global,CutParameters_global,&
			  s_global,qtMIN_List(i),qtMAX_list(i),Q_min,Q_max,yMin_List(i),yMax_List(i),mT_min,mT_max,num)
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecListPY_PTint_Qint_Yint_mTint
  
  subroutine xSecListPY_PTintN_Qint_Yint_mTint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_List,yMax_List,mT_min,mT_max)
    real*8, intent(in) :: qtMIN_list(:),qtMAX_list(:),yMin_List(:),yMax_List(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max,mT_min,mT_max
    call xSecListPY_PTint_Qint_Yint_mTint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_List,yMax_List,&
    mT_min,mT_max,NumPTdefault)
  end subroutine xSecListPY_PTintN_Qint_Yint_mTint
  
  subroutine xSecSingle_PTintN_Qint_Ycomplete_mTint(X,qt_Min,qt_Max,Q_min,Q_max,&
  mT_min,mT_max,num)
    real*8::X,qt_Min,qt_Max
    real*8::yMin_in,yMax_in,Q_min,Q_max,mT_min,mT_max
    integer::num
    
   CallCounter=CallCounter+1
   X=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint_mTint(process_global,includeCuts_global,CutParameters_global,&
				  s_global,qt_Min,qt_Max,Q_min,Q_max,-1000d0,1000d0,mT_min,mT_max,num)
				  
  end subroutine xSecSingle_PTintN_Qint_Ycomplete_mTint
  
  subroutine xSecSingle_PTint_Qint_Ycomplete_mTint(X,qt_Min,qt_Max,Q_min,Q_max,mT_min,mT_max)
    real*8::X,qt_Min,qt_Max
    real*8::yMin_in,yMax_in,Q_min,Q_max,mT_min,mT_max
    
   X=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint_mTint(process_global,includeCuts_global,CutParameters_global,&
			    s_global, qt_Min,qt_Max,Q_min,Q_max,-1000d0,1000d0,mT_min,mT_max,NumPTdefault)
  end subroutine xSecSingle_PTint_Qint_Ycomplete_mTint
  
  
    subroutine xSecSingle_PTintN_Qcomplete_Ycomplete_mTint(X,qt_Min,qt_Max,Q_max,mT_min,mT_max,n,num)
        real*8,intent(in)::qt_Min,qt_Max,Q_max,mT_min,mT_max
        real*8,intent(out)::X
        integer,intent(in)::num,n
        
        X=PreFactor1(process_global(1))*Xsec_PTint_Qcomplete_Ycomplete_mTint(process_global,includeCuts_global,&
        CutParameters_global,s_global,qt_Min,qt_Max,Q_max,y_global,mT_min,mT_max,num,n)
			    
    end subroutine xSecSingle_PTintN_Qcomplete_Ycomplete_mTint
    
    
    subroutine xSecSingle_PTint_Qcomplete_Ycomplete_mTint(X,qt_Min,qt_Max,Q_max,mT_min,mT_max,n)
        real*8,intent(in)::qt_Min,qt_Max,Q_max,mT_min,mT_max
        real*8,intent(out)::X
        integer,intent(in)::n
      
        X=PreFactor1(process_global(1))*Xsec_PTint_Qcomplete_Ycomplete_mTint(process_global,includeCuts_global,&
        CutParameters_global,s_global,qt_Min,qt_Max,Q_max,y_global,mT_min,mT_max,NumPTdefault,n)
			    
    end subroutine xSecSingle_PTint_Qcomplete_Ycomplete_mTint
    
    subroutine xSecSingle_PTint_Qcomplete_Ycomplete_mTint_omp(X,qT,Q_max,mT_min,mT_max,n)
        real*8,intent(in)::Q_max,mT_min,mT_max
        !real*8,intent(out)::X
        integer,intent(in)::n
        integer::i,length
        real*8,intent(in),dimension(:)::qT
        real*8,intent(out),dimension(:)::X
        
        length=size(qT,1)
        
        if (size(X,1) /= length-1 ) then        
            write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_List: sizes of qT and X lists are not equal.'
            write(*,*) 'Evaluation stop'
            stop
        end if
        
        !$OMP PARALLEL DO DEFAULT(SHARED)
        do i=1,length-1
        X(i)=PreFactor1(process_global(1))*Xsec_PTint_Qcomplete_Ycomplete_mTint(process_global,includeCuts_global,&
        CutParameters_global,s_global,qT(i),qT(i+1),Q_max,y_global,mT_min,mT_max,NumPTdefault,n)        
        end do
        !$OMP END PARALLEL DO 
			    
    end subroutine xSecSingle_PTint_Qcomplete_Ycomplete_mTint_omp
    
    
    subroutine xSecSingle_PTint_Ycomplete_mTint_NRW(X,qt_Min,qt_Max,mT_min)
        real*8,intent(in)::qt_Min,qt_Max,mT_min
        real*8,intent(out)::X
        
        X=PreFactor1(process_global(1))*Xsec_PTint_Ycomplete_mTint_NRW(process_global,includeCuts_global,&
        CutParameters_global,s_global,qt_Min,qt_Max,Q_global,y_global,mT_min,NumPTdefault)
			    
    end subroutine xSecSingle_PTint_Ycomplete_mTint_NRW
    
    
    subroutine xSecSingle_PTintN_Ycomplete_mTint_NRW(X,qt_Min,qt_Max,mT_min,num)
        real*8,intent(in)::qt_Min,qt_Max,mT_min
        real*8,intent(out)::X
        integer,intent(in)::num
        
        X=PreFactor1(process_global(1))*Xsec_PTint_Ycomplete_mTint_NRW(process_global,includeCuts_global,&
        CutParameters_global,s_global,qt_Min,qt_Max,Q_global,y_global,mT_min,num)
			    
    end subroutine xSecSingle_PTintN_Ycomplete_mTint_NRW
    
    
    
  
  
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!THE MAIN INTERFACE TO CROSS-SECTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! interface for integer,s,array,array,array,logical,optional, optional
  subroutine MainInterface_isAAAloo(X,process,s,qT,Q,y,mT,includeCuts,CutParameters,Num)
    integer,intent(in)::process					!the number of process
    real*8,intent(in)::s					!Mandelshtam s
    real*8,intent(in),dimension(1:2)::qT			!(qtMin,qtMax)
    real*8,intent(in),dimension(1:2)::Q				!(Qmin,Qmax)
    real*8,intent(in),dimension(1:2)::y				!(ymin,ymax)
    real*8,intent(in),dimension(1:2)::mT            !(mTmin,mTmax)
    logical,intent(in)::includeCuts				!include cuts
    real*8,intent(in),dimension(1:4),optional::CutParameters	!(p1,p2,eta1,eta2)
    integer,intent(in),optional::Num				!number of sections
    
    real*8::X
  
  integer::nn
  real*8,dimension(1:4)::CutParam
  integer,dimension(1:3)::ppp
  
  !!! determine number of sections
  if(present(Num)) then
    nn=Num
  else
    nn=NumPTdefault
  end if
  
  !!! determine cut parameters
  if(includeCuts) then
    if(present(CutParameters)) then
      CutParam=CutParameters
    else
      write(*,*) 'ERROR: arTeMiDe_DY: called includeCuts=true, while CutParameters are undefined'
      write(*,*) 'Evaluation stop'
      stop
    end if
  else
    CutParam=(/0d0,0d0,0d0,0d0/)
  end if
  
  ppp=processArrayFromInteger(process)
  
  !!!! evaluation
  CallCounter=CallCounter+1
  X=PreFactor1(ppp(1))*Xsec_PTint_Qint_Yint_mTint(ppp,includeCuts,CutParameters,&
				  s,qT(1),qT(2),Q(1),Q(2),y(1),y(2),mT(1),mT(2),nn)
  
  end subroutine MainInterface_isAAAloo
  
  !!!! interface for array,s,array,array,array,logical,optional, optional
  subroutine MainInterface_AsAAAloo(X,process,s,qT,Q,y,mT,includeCuts,CutParameters,Num)
!   function xSec_DY(process,s,qT,Q,y,includeCuts,CutParameters,Num)
    integer,intent(in),dimension(1:3)::process			!the number of process
    real*8,intent(in)::s					        !Mandelshtam s
    real*8,intent(in),dimension(1:2)::qT			!(qtMin,qtMax)
    real*8,intent(in),dimension(1:2)::Q				!(Qmin,Qmax)
    real*8,intent(in),dimension(1:2)::y				!(ymin,ymax)
    real*8,intent(in),dimension(1:2)::mT            !(mTmin,mTmax)
    logical,intent(in)::includeCuts				!include cuts
    real*8,intent(in),dimension(1:4),optional::CutParameters	!(p1,p2,eta1,eta2)
    integer,intent(in),optional::Num				!number of sections
    
    real*8::X
  
  integer::nn
  real*8,dimension(1:4)::CutParam
  
  
  !! determine number of sections
  if(present(Num)) then
    nn=Num
  else
    nn=NumPT_auto(real(qT(2)-qT(1)),real((Q(2)+Q(1))/2.))
  end if
    
  !!! determine cut parameters
  if(includeCuts) then
    if(present(CutParameters)) then
      CutParam=CutParameters
    else
      write(*,*) 'ERROR: arTeMiDe_DY_mT: called includeCuts=true, while CutParameters are undefined'
      write(*,*) 'Evaluation stop'
      stop
    end if
  else
    CutParam=(/0d0,0d0,0d0,0d0/)
  end if
  
  !!!! evaluation
  CallCounter=CallCounter+1
  X=PreFactor1(process(1))*Xsec_PTint_Qint_Yint_mTint(process,includeCuts,CutParameters,&
				  s,qT(1),qT(2),Q(1),Q(2),y(1),y(2),mT(1),mT(2),nn)
  
  end subroutine MainInterface_AsAAAloo
  
  subroutine xSec_DY_mT_List(X,process,s,qT,Q,y,mT,includeCuts,CutParameters,Num,NQ)
    integer,intent(in),dimension(:,:)::process			!the number of process
    real*8,intent(in),dimension(:)::s				!Mandelshtam s
    real*8,intent(in),dimension(:,:)::qT			!(qtMin,qtMax)
    real*8,intent(in),dimension(:,:)::Q				!(Qmin,Qmax)
    real*8,intent(in),dimension(:,:)::y				!(ymin,ymax)
    real*8,intent(in),dimension(:,:)::mT            !(mTmin,mTmax)
    logical,intent(in),dimension(:)::includeCuts		!include cuts
    real*8,intent(in),dimension(:,:)::CutParameters	        !(p1,p2,eta1,eta2)
    integer,intent(in),dimension(:),optional::Num		!number of sections
    integer,intent(in)::NQ
    real*8,dimension(:),intent(out)::X
    real*8,allocatable::XqT(:),NqT(:,:)
    integer :: i,length,lenNqT,lenXqT,j,k
    integer,allocatable::nn(:),NumqT(:)
    real*8,allocatable::yN(:,:),QN(:,:),mTN(:,:),sN(:),CutParametersN(:,:)
    integer,allocatable::processN(:,:)
    logical,allocatable::includeCutsN(:)
    real*8::delqT,Sec,XSec_qT
    
    length=size(s)
    allocate(NumqT(1:size(Num,1)))
    do i=1,size(Num,1)
        NumqT(i)=Num(i)/4
    end do
    lenNqT=sum(NumqT(:))+1
    lenXqT=sum(NumqT(:))
    
     allocate(NqT(1:lenXqT,1:2))
     allocate(XqT(1:lenXqT))
     allocate(sN(1:lenXqT))
     allocate(includeCutsN(1:lenXqT))
     allocate(CutParametersN(1:lenXqT,1:4))
     allocate(processN(1:lenXqT,1:3))
     allocate(yN(1:lenXqT,1:2))
     allocate(QN(1:lenXqT,1:2))
     allocate(mTN(1:lenXqT,1:2))
    
    j=0
    do i=1,length
        do k=1,NumqT(i)
            j=j+1
            sN(j)=s(i)
            includeCutsN(j)=includeCuts(i)
            yN(j,1:2)=y(i,1:2)
            mTN(j,1:2)=mT(i,1:2)
            QN(j,1:2)=Q(i,1:2)
            CutParametersN(j,1:4)=CutParameters(i,1:4)
            processN(j,1:3)=process(i,1:3)
        end do
    end do
    

     
    
    
    !!! cheking sizes
    if(size(X)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY_mT: xSec_DY_mT_List: sizes of xSec and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(process,1)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_mT_List: sizes of process and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(qT,1)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY_mT: xSec_DY_mT_List: sizes of qT and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(y,1)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY_mT: xSec_DY_mT_List: sizes of y and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if (size(mT,1)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY_mT: xSec_DY_mT_List: sizes of mT and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(Q,1)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY_mT: xSec_DY_mT_List: sizes of Q and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(includeCuts)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY_mT: xSec_DY_mT_List: sizes of includeCuts and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(CutParameters,1)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY_mT: xSec_DY_mT_List: sizes of CutParameters and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(process,2)/=3) then
      write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_mT_List: process list must be (:,1:3).'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(qT,2)/=2) then
      write(*,*) 'ERROR: arTeMiDe_DY_mT: xSec_DY_mT_List: qt list must be (:,1:2).'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(y,2)/=2) then
      write(*,*) 'ERROR: arTeMiDe_DY_mT: xSec_DY_mT_List: y list must be (:,1:2).'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(Q,2)/=2) then
      write(*,*) 'ERROR: arTeMiDe_DY_mT: xSec_DY_mT_List: Q list must be (:,1:2).'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(mT,2)/=2) then
      write(*,*) 'ERROR: arTeMiDe_DY_mT: xSec_DY_mT_List: mT list must be (:,1:2).'
      write(*,*) 'Evaluation stop'
      stop
    end if
    
    
    CallCounter=CallCounter+length
    
    
     allocate(nn(1:lenXqT))
     if(present(Num)) then
         if(size(Num,1)/=length) then
 	  write(*,*) 'ERROR: arTeMiDe_DY_mT: xSec_DY_mT_List: sizes of Num and s lists are not equal.'
 	  write(*,*) 'Evaluation stop'
         stop
         end if
         do i=1,lenXqT
            nn(i)=4
         end do
     else
         do i=1,length
             nn(i)=NumPT_auto(real(qT(i,2)-qT(i,1)),real((Q(i,2)+Q(i,1))/2.))
         end do
     end if
    


    j=0
    do i=1,length
        delqT=(qT(i,2)-qT(i,1))/real(NumqT(i))
        do k=1,NumqT(i)
        j=j+1
        NqT(j,1)=qT(i,1)+real(k-1)*delqT
        NqT(j,2)=qT(i,1)+real(k)*delqT
        end do
     end do
    

    !$OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,lenXqT
        XqT(i)=PreFactor1(processN(i,1))*Xsec_PTint_Qint_Yint_mTint(processN(i,1:3),includeCutsN(i),&
        CutParametersN(i,1:4),sN(i),NqT(i,1),NqT(i,2),QN(i,1),QN(i,2),yN(i,1),yN(i,2),&
        mTN(i,1),mTN(i,2),nn(i),4)
    end do
    !$OMP END PARALLEL DO
    
     j=0
     do i=1,length
         Sec=0d0
         do k=1,NumqT(i)
             j=j+1
             Sec=XqT(j)+Sec
         end do
         X(i)=Sec
     end do
    
    
    deallocate(NqT)
    deallocate(XqT)
    deallocate(CutParametersN)
    deallocate(yN)
    deallocate(QN)
    deallocate(mTN)
    deallocate(sN)
    deallocate(includeCutsN)
    deallocate(processN)
    deallocate(nn)
    
    
  end subroutine xSec_DY_mT_List
  
    subroutine alphaQCD(alphaS,mu)
        real*8,intent(in)::mu
        real*8,intent(out)::alphaS
        
        alphaS=AsPythia(mu)
    end subroutine alphaQCD
  
  
  
end module TMDX_DY_mT
