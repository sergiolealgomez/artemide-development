!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.4
!
!	Evaluation of the leptonic cuts for the DrellYan
!
!	if you use this module please, quote 1706.01473
!
!	ver 1.0: release (AV, 10.05.2017)
!	ver 1.32 update nessacary for parallelisation (AV,17.09.2018)
!	ver 1.32 CutFactor4 added, asymetric cuts are introduced (AV,03.12.2018)
!	ver.1.4  Deleted old routines. Incapsulated variables (AV. 18.01.2019(
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module LeptonCutsDY_mT
use IO_functions
implicit none
private
! public

!!!Parameters of cut and cross-section
!!!!! this is array =(/  pT1lim,pT2lim,etaMin,etaMax,Exp(2etaMin),exp(2etaMax) /)
real*8::cutParam_global(1:6)
real*8::leppT_upper,neupT_upper,leppT_lower
logical::include_leppT_upper=.false.
logical::include_neupT_upper=.false.
logical::Sec_diff_lT=.false.
logical::Sec_diff_nT=.false.


!!! number of divisions in Simpsons
integer,parameter::num=128!64!24!!
!! Tolerance (absolute)
real*8,parameter::tolerance=0.000001d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! the kinematic variables are passed via an aray
!!! (q_T,Q,y, Q^2 , qt^2/Q^2+qt^2 , Sqrt[Q^2+qT^2] )

public:: SetCutParameters
public:: CutFactor4,CutFactorA
public:: Set_neutrino_pTupper,Set_lepton_pTupper


  interface SetCutParameters
    module procedure SetCutParameters_sym,SetCutParameters_asym
  end interface

 contains


!SEt parameters of cut
subroutine SetCutParameters_sym(pT_in,etaMin_in,etaMax_in)
  real*8::pT_in,etaMax_in,etaMin_in

  cutParam_global=(/ pT_in**2,&
	      pT_in**2,&
	      etaMin_in,&
	      etaMax_in,&
	      EXP(2*etaMin_in),&
	      EXP(2*etaMax_in) /)

end subroutine SetCutParameters_sym

!SEt parameters of cut
! with asymetric cuts for pT
subroutine SetCutParameters_asym(pT1_in,pT2_in,etaMin_in,etaMax_in)
  real*8::pT1_in,pT2_in,etaMax_in,etaMin_in

  !! for definetines we order pt1>pt2
  if(pT1_in>=pT2_in) then
  cutParam_global=(/ pT1_in**2,&
	      pT2_in**2,&
	      etaMin_in,&
	      etaMax_in,&
	      EXP(2d0*etaMin_in),&
	      EXP(2d0*etaMax_in) /)
  else
  cutParam_global=(/ pT2_in**2,&
	      pT1_in**2,&
	      etaMin_in-tolerance,&
	      etaMax_in+tolerance,&
	      EXP(2d0*etaMin_in),&
	      EXP(2d0*etaMax_in) /)
  end if

end subroutine SetCutParameters_asym


subroutine Set_lepton_pTupper(diff_lT,lepton_lower_pT,lepton_upper_pT)
    real*8::lepton_upper_pT,lepton_lower_pT
    logical:: diff_lT
    Sec_diff_lT=diff_lT
    leppT_upper=lepton_upper_pT
    leppT_lower=lepton_lower_pT
    include_leppT_upper=.true.
    
    !write(*,*) "Hallo 4"


end subroutine Set_lepton_pTupper

subroutine Set_neutrino_pTupper(diff_nT,neutrino_pT)
    real*8::neutrino_pT
    logical:: diff_nT
    Sec_diff_nT=diff_nT
    neupT_upper=neutrino_pT
    include_neupT_upper=.true.


end subroutine Set_neutrino_pTupper

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	This function survives after many modification of module. For different version of integral evaluation see /history
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Integration is done by more accurate method. Secting the phi-plane, and define the boundaries of eta
!!!! The integal over eta is done exactly
!!!! Unfortunately works only for large enough Q.
!!!! Gives very precise result 10^-5 accuracy
!!!!
!!!! CutParameters = (/ kT1, kT2, etaMin, etaMax /)
function CutFactor4(qT,Q_in,y_in,mT,CutParameters,neutrino,lep_type)
  real*8:: qT,CutFactor4,Q_in,y_in,mT
  real*8:: dphi,leppT2_upper,neupT2_upper
  integer :: i
  real*8,dimension(1:6)::varC
  real*8,dimension(1:7)::var
  !real*8,dimension(1:4),optional,intent(in)::CutParameters
  real*8,optional,intent(in)::CutParameters(:)
  logical,intent(in)::neutrino
  integer,optional,intent(in)::lep_type
  
  if(present(CutParameters)) then
  
    if (size(CutParameters)==4) then

    varC(1)=CutParameters(1)**2
    varC(2)=CutParameters(2)**2
    varC(3)=CutParameters(3)
    varC(4)=CutParameters(4)
    varC(5)=exp(2d0*CutParameters(3))
    varC(6)=exp(2d0*CutParameters(4))
    
    else if ((size(CutParameters)==5).and.(.not.neutrino)) then
    
    varC(1)=CutParameters(1)**2d0
    leppT2_upper=CutParameters(2)**2d0
    varC(2)=CutParameters(3)**2d0
    varC(3)=CutParameters(4)
    varC(4)=CutParameters(5)
    varC(5)=exp(2d0*CutParameters(4))
    varC(6)=exp(2d0*CutParameters(5))
    
    else if ((size(CutParameters)==5).and.(neutrino)) then
    
    varC(1)=CutParameters(1)**2d0
    varC(2)=CutParameters(2)**2d0
    neupT2_upper=CutParameters(3)**2d0
    varC(3)=CutParameters(4)
    varC(4)=CutParameters(5)
    varC(5)=exp(2d0*CutParameters(4))
    varC(6)=exp(2d0*CutParameters(5))
    
    else if (size(CutParameters)==6) then
    
    varC(1)=CutParameters(1)**2d0
    leppT2_upper=CutParameters(2)**2d0
    varC(2)=CutParameters(3)**2d0
    neupT2_upper=CutParameters(4)**2d0
    varC(3)=CutParameters(5)
    varC(4)=CutParameters(6)
    varC(5)=exp(2d0*CutParameters(5))
    varC(6)=exp(2d0*CutParameters(6))
    
    end if
    
  else
    varC=cutParam_global
  end if

  var=(/qT,Q_in,y_in,Q_in**2,qT**2/(Q_in**2+qT**2),SQRT(Q_in**2+qT**2),mT/)

  if(varC(3)<varC(4)) then
  if(y_in<varC(3) .or. y_in>varC(4)) then
    CutFactor4=0d0
  else
    CutFactor4=0d0
    dphi=2d0*3.14159265358979d0/num
    
    if (size(CutParameters)==4) then
    
    do i=0,num
        if (present(lep_type)) then
      CutFactor4=CutFactor4+Simp(i,num)*IntegralOverEtaFixedPhiEXACT(var,varC,i*dphi,&
      lep_type=lep_type)
        else
      CutFactor4=CutFactor4+Simp(i,num)*IntegralOverEtaFixedPhiEXACT(var,varC,real(i)*dphi)
        end if
    end do
    
    else if ((size(CutParameters)==5).and.(.not.neutrino)) then
    
    do i=0,num
      CutFactor4=CutFactor4+Simp(i,num)*IntegralOverEtaFixedPhiEXACT(var,varC,i*dphi,lpT2=leppT2_upper)
    end do
    
    else if ((size(CutParameters)==5).and.(neutrino)) then
    
    do i=0,num
      CutFactor4=CutFactor4+Simp(i,num)*IntegralOverEtaFixedPhiEXACT(var,varC,i*dphi,npT2=neupT2_upper)
    end do
    
    else if (size(CutParameters)==6) then
    
    do i=0,num
      CutFactor4=CutFactor4+Simp(i,num)*IntegralOverEtaFixedPhiEXACT(var,varC,i*dphi,lpT2=leppT2_upper,&
      npT2=neupT2_upper)
    end do
    
    end if
    
      CutFactor4=CutFactor4*dphi/3d0
  end if
  else
    CutFactor4=0d0
  end if

end function CutFactor4

!!!! integral for anti-symmetric part of lepton tensor
function CutFactorA(qT,Q_in,y_in,mT,CutParameters)
  real*8:: qT,CutFactorA,Q_in,y_in,mT
  real*8:: dphi
  integer :: i
  real*8,dimension(1:6)::varC
  real*8,dimension(1:7)::var
  real*8,dimension(1:4),optional,intent(in)::CutParameters

  if(present(CutParameters)) then
    varC(1)=CutParameters(1)**2
    varC(2)=CutParameters(2)**2
    varC(3)=CutParameters(3)
    varC(4)=CutParameters(4)
    varC(5)=exp(2d0*CutParameters(3))
    varC(6)=exp(2d0*CutParameters(4))
  else
    varC=cutParam_global
  end if

  var=(/qT,Q_in,y_in,Q_in**2,qT**2/(Q_in**2+qT**2),SQRT(Q_in**2+qT**2),mT/)

  if(varC(3)<varC(4)) then
  if(y_in<varC(3) .or. y_in>varC(4)) then
    CutFactorA=0d0
  else
    CutFactorA=0d0
    dphi=2d0*3.14159265358979d0/num
    do i=0,num
      CutFactorA=CutFactorA+Simp(i,num)*IntegralOverEtaFixedPhiEXACT_A(var,varC,i*dphi)
    end do
      CutFactorA=CutFactorA*dphi/3d0
  end if
  else
    CutFactorA=0d0
  end if

end function CutFactorA



!it is the same function as IntegralOverEtaFixedPhi, but the integration is done exactly.
! the boundaries are defined as before
!	       1     2    3       4          5                       6                 7
! var =(/ qT,  Q_in, y_in,   Q_in**2,  qT**2/(Q_in**2+qT**2), SQRT(Q_in**2+qT**2), mT/)
! varC=(/ pT1, pT2,  etaMin, etaMax,   EXP(2*etaMin),         EXP(2*etaMax) /)
function IntegralOverEtaFixedPhiEXACT(var,varC,phi,lpT2,npT2,lep_type)
  real*8,dimension(1:6)::varC
  real*8,dimension(1:7)::var
  real*8:: IntegralOverEtaFixedPhiEXACT,phi
  real*8::cosp,chhy,l1square,l1,chhynum,chhyden,term1,term2,term3
  real*8,optional,intent(in)::lpT2,npT2
  real*8::m
  integer,optional,intent(in)::lep_type
  
    !if lep_type=1 m is the mass of a electron or anti-electron
    !if lep_type=2 m is the mass of a muon or anti-muon
    
    if (present(lep_type)) then
    if (lep_type==1) then
        m=0.510998928d0*1d-3 !GeV
    else if (lep_type==2) then
        m=105.6583715d0*1d-3 !GeV
    end if
    end if
  
    if (present(lep_type)) then !mass if
    
    if ((present(lpT2)).and.(present(npT2)).and.&
    (Integrand2THETA(var,varC,phi,lpT2=lpT2,npT2=npT2,lep_type=lep_type)==0)) then
    !!!! Here is the point of assumption!! because if Q~qT the integration region is difficult and the boundaries may be inproper
    IntegralOverEtaFixedPhiEXACT=0d0
    else if ((present(lpT2)).and.(.not.present(npT2)).and.&
    (Integrand2THETA(var,varC,phi,lpT2=lpT2,lep_type=lep_type)==0)) then
    !!!! Here is the point of assumption!! because if Q~qT the integration region is difficult and the boundaries may be inproper
    IntegralOverEtaFixedPhiEXACT=0d0
    else if ((.not.present(lpT2)).and.(present(npT2)).and.&
    (Integrand2THETA(var,varC,phi,npT2=npT2,lep_type=lep_type)==0)) then
    !!!! Here is the point of assumption!! because if Q~qT the integration region is difficult and the boundaries may be inproper
    IntegralOverEtaFixedPhiEXACT=0d0
    else if ((.not.present(lpT2)).and.(.not.present(npT2)).and.&
    (Integrand2THETA(var,varC,phi,lep_type=lep_type)==0)) then
    !!!! Here is the point of assumption!! because if Q~qT the integration region is difficult and the boundaries may be inproper
    IntegralOverEtaFixedPhiEXACT=0d0
    else !cuts if
    
    cosp=COS(phi)
    
    l1=(var(7)**2d0*SQRT(var(7)**2d0+var(1)**2d0)+var(7)**2d0*var(1)*cosp)&
    /(2d0*var(7)**2d0+(2d0*var(1)**2)*(1d0-cosp**2))
    l1square=l1**2d0

    chhynum=(2d0*l1**2-var(7)**2d0-m**2d0+var(2)**2d0+2d0*l1*SQRT(l1**2+var(1)**2d0-2d0*cosp*l1*var(1)))
    chhyden=(2d0*l1*SQRT(1+(m/l1)**2d0)*SQRT(var(2)**2d0+var(1)**2d0))
    chhy=chhynum/chhyden

    term1=1/(4d0*SQRT(1+(m/l1)**2d0)*SQRT(var(2)**2d0+var(1)**2d0)*SQRT(-1d0+chhy**2d0+0.000001d0))
    term2=-l1**2d0-m**2d0+chhy*l1*SQRT(1+(m/l1)**2d0)*SQRT(var(2)**2d0+var(1)**2d0)
    term3=ABS(-(4d0*l1**2-6d0*l1*var(1)*cosp+2d0*var(1)**2d0)/(SQRT(l1**2d0+var(1)**2d0-2d0*l1*var(1)*cosp))&
    -4d0*l1+2d0*var(1)*cosp)
    

    IntegralOverEtaFixedPhiEXACT=term1*term2/term3
    
    end if !cuts if
    
    else !mass if

    if ((present(lpT2)).and.(present(npT2)).and.&
    (Integrand2THETA(var,varC,phi,lpT2=lpT2,npT2=npT2)==0)) then
    !!!! Here is the point of assumption!! because if Q~qT the integration region is difficult and the boundaries may be inproper
    IntegralOverEtaFixedPhiEXACT=0d0
    else if ((present(lpT2)).and.(.not.present(npT2)).and.&
    (Integrand2THETA(var,varC,phi,lpT2=lpT2)==0)) then
    !!!! Here is the point of assumption!! because if Q~qT the integration region is difficult and the boundaries may be inproper
    IntegralOverEtaFixedPhiEXACT=0d0
    else if ((.not.present(lpT2)).and.(present(npT2)).and.&
    (Integrand2THETA(var,varC,phi,npT2=npT2)==0)) then
    !!!! Here is the point of assumption!! because if Q~qT the integration region is difficult and the boundaries may be inproper
    IntegralOverEtaFixedPhiEXACT=0d0
    else if ((.not.present(lpT2)).and.(.not.present(npT2)).and.&
    (Integrand2THETA(var,varC,phi)==0)) then
    !!!! Here is the point of assumption!! because if Q~qT the integration region is difficult and the boundaries may be inproper
    IntegralOverEtaFixedPhiEXACT=0d0
    else !cuts if

    cosp=COS(phi)
    
    l1=(var(7)**2d0*SQRT(var(7)**2d0+var(1)**2d0)+var(7)**2d0*var(1)*cosp)&
    /(2d0*var(7)**2d0+(2d0*var(1)**2)*(1d0-cosp**2))
    l1square=l1**2d0

    chhynum=(2d0*l1**2-var(7)**2d0+var(2)**2d0+2d0*l1*SQRT(l1**2+var(1)**2d0-2d0*cosp*l1*var(1)))
    chhyden=(2d0*l1*SQRT(var(2)**2d0+var(1)**2d0))
    chhy=chhynum/chhyden

    term1=1/(4d0*SQRT(var(2)**2d0+var(1)**2d0)*SQRT(-1d0+chhy**2d0+0.000001d0))
    term2=-l1**2d0+chhy*l1*SQRT(var(2)**2d0+var(1)**2d0)
    term3=ABS(-(4d0*l1**2-6d0*l1*var(1)*cosp+2d0*var(1)**2d0)/(SQRT(l1**2d0+var(1)**2d0-2d0*l1*var(1)*cosp))&
    -4d0*l1+2d0*var(1)*cosp)

    IntegralOverEtaFixedPhiEXACT=term1*term2/term3

    end if !end cuts if
    
    end if !end mass if
  
  
end function IntegralOverEtaFixedPhiEXACT


!it is the same function as IntegralOverEtaFixedPhi, but the integration is done exactly.
! the boundaries are defined as before
!	   1     2    3       4          5                       6
! var =(/ qT,  Q_in, y_in,   Q_in**2,  qT**2/(Q_in**2+qT**2), SQRT(Q_in**2+qT**2), mT/)
! varC=(/ pT1, pT2,  etaMin, etaMax,   EXP(2*etaMin),         EXP(2*etaMax) /)
!!! ANTI SYMMETRIC ONE!
function IntegralOverEtaFixedPhiEXACT_A(var,varC,phi)
  real*8,dimension(1:6)::varC
  real*8,dimension(1:7)::var
  real*8:: IntegralOverEtaFixedPhiEXACT_A,phi


  if(Integrand2THETA(var,varC,phi)==0) then
  !!!! Here is the point of assumption!! because if Q~qT the integration region is difficult and the boundaries may be inproper
  IntegralOverEtaFixedPhiEXACT_A=0d0
  else


  IntegralOverEtaFixedPhiEXACT_A=0

  end if
end function IntegralOverEtaFixedPhiEXACT_A



!!!!the theta (0 or 1) function of the integrand in the coordinates rapidity-angle
!! Search for the boundary of integration between t1 and t2 at fixed phi
!! by devision on 2
!	   1     2    3       4          5                       6
! var =(/ qT,  Q_in, y_in,   Q_in**2,  qT**2/(Q_in**2+qT**2), SQRT(Q_in**2+qT**2), mT/)
! varC=(/ pT1, pT2,  etaMin, etaMax,   EXP(2*etaMin),         EXP(2*etaMax) /)
function Integrand2THETA(var,varC,p1,lpT2,npT2,lep_type)
  real*8,dimension(1:6)::varC
  real*8,dimension(1:7)::var
  real*8:: h1max,h1min,p1
  integer::Integrand2THETA
  real*8::cosp1,chhy,l1square,l2square,exp2h2,l1,chhynum,chhyden
  real*8,optional,intent(in)::lpT2,npT2
  integer::theta
  integer,optional,intent(in)::lep_type
  real*8::m
  logical::zero

    cosp1=COS(p1)

    if (present(lep_type)) then !mass if
        if(lep_type==1) then
            m=0.510998928d0*1d-3 !GeV
        else
            m=105.6583715d0*1d-3 !GeV
        end if
        
        l1=(var(7)**2d0*SQRT(var(7)**2d0+var(1)**2d0)+var(7)**2d0*var(1)*cosp1)&
        /(2d0*var(7)**2d0+(2d0*var(1)**2)*(1d0-cosp1**2))
        l1square=l1**2d0


        chhynum=(2d0*l1**2-var(7)**2d0-m**2d0+var(2)**2d0+2d0*l1*SQRT(l1**2+var(1)**2d0-2d0*cosp1*l1*var(1)))
        chhyden=(2d0*l1*SQRT(1+(m/l1)**2d0)*SQRT(var(2)**2d0+var(1)**2d0))
        chhy=chhynum/chhyden
        
        if (chhy>=1d0) then
        h1max=var(3)+ACOSH(chhy)
        h1min=var(3)-ACOSH(chhy)
        zero=.false.
        else
        h1max=1000d0
        h1min=-1000d0
        zero=.true.
        end if

        l2square=var(1)**2+l1square-2d0*var(1)*l1*cosp1
        !exp2h2=(EXP(2d0*var(3)+h1max)*var(6)-EXP(var(3)+2d0*h1max)*l1)/(EXP(h1max)*var(6)-EXP(var(3))*l1)
        
    
    else !mass if
    
        l1=(var(7)**2d0*SQRT(var(7)**2d0+var(1)**2d0)+var(7)**2d0*var(1)*cosp1)/&
        (2d0*var(7)**2d0+(2d0*var(1)**2d0)*(1d0-cosp1**2d0))
        l1square=l1**2d0

        chhynum=(2d0*l1**2d0-var(7)**2d0+var(2)**2d0+2d0*l1*SQRT(l1**2d0+var(1)**2d0-2d0*cosp1*l1*var(1)))
        chhyden=(2d0*l1*SQRT(var(2)**2d0+var(1)**2d0))
        chhy=chhynum/chhyden

        if (chhy>=1d0) then
        h1max=var(3)+ACOSH(chhy)
        h1min=var(3)-ACOSH(chhy)
        zero=.false.
        else
        h1max=1000d0
        h1min=-1000d0
        zero=.true.
        end if

        l2square=var(1)**2+l1square-2d0*var(1)*l1*cosp1
        !exp2h2=(EXP(2d0*var(3)+h1max)*var(6)-EXP(var(3)+2d0*h1max)*l1)/(EXP(h1max)*var(6)-EXP(var(3))*l1)
        
        
    end if !end mass if
    
    if ((present(lpT2)).and.(present(npT2))) then
    
    if ((l1square>varC(1)).and.(l2square>varC(2)).and.((varC(4)>h1max).and.(h1min>varC(3))).and.&
    (.not.zero).and.(lpT2>l1square).and.(npT2>l2square)) then
    Integrand2THETA=1
    else
    Integrand2THETA=0
    end if
    
    else if ((present(lpT2)).and.(.not.present(npT2))) then
    
    if ((l1square>varC(1)).and.(l2square>varC(2)).and.((varC(4)>h1max).and.(h1min>varC(3))).and.&
    (.not.zero).and.(lpT2>l1square)) then
    Integrand2THETA=1
    else
    Integrand2THETA=0
    end if
    
    else if ((.not.present(lpT2)).and.(present(npT2))) then
    
    if ((l1square>varC(1)).and.(l2square>varC(2)).and.((varC(4)>h1max).and.(h1min>varC(3))).and.&
    (.not.zero).and.(npT2>l2square)) then
    Integrand2THETA=1
    else
    Integrand2THETA=0
    end if
    
    else if ((.not.present(lpT2)).and.(.not.present(npT2))) then
    
    if ((l1square>varC(1)).and.(l2square>varC(2)).and.((varC(4)>h1max).and.(h1min>varC(3))).and.&
    (.not.zero)) then
    Integrand2THETA=1
    else
    Integrand2THETA=0
    end if
    
    end if
   
end function Integrand2THETA



function Simp(i,n)
  integer::i,n
  real*8::Simp
  if((i==0).or.(i==n)) then
  Simp=1d0
  else
    if(MOD(i,2)==1) then
    Simp=4d0
    else
    Simp=2d0
    end if
  end if
end function Simp



end module LeptonCutsDY_mT
