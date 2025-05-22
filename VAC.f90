! *** 5 May -> CRS and VAC for the given Atmosphere ***
SUBROUTINE SK_prof_forATM
CHARACTER HEAD*30,GAS*3,ATMO*30
PARAMETER (JM=65,MCO2=28,MH2O=6,MSO2=4,NKTR=32)  ! JM - number of STANDARD CRS sections
DIMENSION KT_CO2(NKTR), KT_H2O(NKTR),KT_SO2(NKTR),VAC(200,NKTR)
DIMENSION GAS(3),ZA(200),WA(200),TA(200),ROCO2(200),ROH2O(200),ROSO2(200) &
                        ,ACO2(200,MCO2),AH2O(200,MH2O),ASO2(200,MSO2)
DIMENSION W(JM),TV(JM),SCO2(JM,MCO2),SH2O(JM,MH2O),SSO2(JM,MSO2),SCO2T(JM,MCO2)
COMMON/IGA/ICO2,IH2O,ISO2
DATA KT_SO2/23*0,1,2,3,4,5*0/,KT_H2O/1,2,3,4,22*0,5,3*0,6,0/,&
 KT_CO2/4*0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28/

DATA IBG/0/
SAVE IBG
OPEN(1,FILE='ATMOSPHERE.INP')
READ(1,2)ATMO
CLOSE(1)
2 FORMAT(A30)
! *** the given Atmosphere *** !
CALL NEWATM(ZA,WA,TA,JAM,NG,GAS,ROCO2,ROH2O,ROSO2,ATMO)  ! JAM - number LEVELs it the given atmospheric profile
 TV=0.
! ------------------------ !
IF(IBG==0)THEN
!*** Standard CRS reading (TH(...) see in CO2_T_COR ***)!
IBG=1
OPEN(9,FILE='./F_LW/W_STAND')
DO J=1,JM
READ(9,*)W(J)
END DO
CLOSE(9)

!*** TV(...) ***!
DO J=1,JM
WWWWW=W(J)
DO JJ=2,JAM
W1=WA(JJ-1) ; W2=WA(JJ)
IF(W2==0.)EXIT
IF(W1>=WWWWW.AND.W2<=WWWWW)EXIT
END DO
C2=(W1-WWWWW)/(W1-W2) ; C1=1.-C2
TV(J)=C1*TA(JJ-1)+C2*TA(JJ)
END DO 

! *** CO2 *** !
  IF(ICO2>0)THEN
OPEN(10,FILE='./F_LW/S.CO2')
DO J=1,JM
READ(10,*)SCO2(J,:)
END DO
CLOSE(10)
  END IF
! *** H2O *** !
  IF(IH2O>0)THEN
OPEN(10,FILE='./F_LW/S.H2O')
DO J=1,JM
READ(10,*)SH2O(J,:)
END DO
CLOSE(10)
  END IF
! *** SO2 *** !
  IF(ISO2>0)THEN
OPEN(10,FILE='./F_LW/S.SO2')
DO J=1,JM
READ(10,*)SSO2(J,:)
END DO
CLOSE(10)
  END IF
END IF !*** Standard CRS reading OK! ***!
! ------------------------ !

! ==== Atmosphera -> CRS ============ !

! ***** CO2 T-correction!  ***** !
CALL CO2_COR_T(W,SCO2,SCO2T,TV)

DO J=1,JAM
WWWWW=WA(J) 
DO JJ=2,JM
W1=W(JJ-1) ; W2=W(JJ)
IF(W2==0.)EXIT
IF(W2<=WWWWW)EXIT
END DO
C1=1. ; C2=0.
IF(W1>=WWWWW)THEN
C2=(W1-WWWWW)/(W1-W2) ; C1=1.-C2
ENDIF
IF(JJ<=JM)THEN
  IF(ISO2>0)ASO2(J,:)=C1*SSO2(JJ-1,:)+C2*SSO2(JJ,:) ! *** SO2 *** !
  IF(IH2O>0)AH2O(J,:)=C1*SH2O(JJ-1,:)+C2*SH2O(JJ,:) ! *** H2O *** !
 !###  IF(ICO2>0)ACO2(J,:)=C1*SCO2(JJ-1,:)+C2*SCO2(JJ,:) ! *** CO2 no T *** !
 IF(ICO2>0)ACO2(J,:)=C1*SCO2T(JJ-1,:)+C2*SCO2T(JJ,:) ! *** CO2-T *** !
ELSE
ASO2(J,:)=0. ; AH2O(J,:)=0. ; ACO2(J,:)=0.
END IF
END DO

! ****** VAC (Vol. Absorpt. Coeff.) profiles in EACH Gas ****** !
DO J=1,JAM
ACO2(J,:)=EXP(ACO2(J,:))*ROCO2(J)
AH2O(J,:)=EXP(AH2O(J,:))*ROH2O(J)
ASO2(J,:)=EXP(ASO2(J,:))*ROSO2(J)
END DO
! ************************** !
! ****** VAC profiles in KD-Term ****** !VAC(200,NKTR) KT_SO2/2
VAC=0.
DO KK=1,NKTR
K1=KT_CO2(KK) ; K2=KT_H2O(KK) ; K3=KT_SO2(KK)
IF(K1>0)VAC(:,KK)=VAC(:,KK)+ACO2(:,K1)
IF(K2>0)VAC(:,KK)=VAC(:,KK)+AH2O(:,K2)
IF(K3>0)VAC(:,KK)=VAC(:,KK)+ASO2(:,K3)
END DO
 
  OPEN(40,FILE='V_'//ATMO)
DO J=1,JAM
WRITE(40,32)ZA(J),WA(J),VAC(J,:)
32 FORMAT(F7.3,F9.5,32E12.4)
END  DO
CLOSE(40)

!**** TOTAL VALIDATION ***!
!### OPEN(40,FILE='TV')
!### DO J=1,JM
!### WRITE(40,*)W(J),TV(J)
!### ENDDO
!### CLOSE(40)

!###OPEN(40,FILE='KD-SO2.ATM')
!###DO J=1,JAM
!###WRITE(40,39)WA(J),ASO2(J,:)
!###39 FORMAT(F9.5,6F11.5)
!### ENDDO


!###OPEN(40,FILE='KD-H2O.ATM')
!###DO J=1,JAM
!###WRITE(40,39)WA(J),AH2O(J,:)
!###39 FORMAT(F9.5,6F11.5)
!### ENDDO

!###OPEN(40,FILE='KD-CO2_T.ATM')
!###DO J=1,JAM
!###WRITE(40,39)WA(J),ACO2(J,:)
!###39 FORMAT(F9.5,28F11.5)
!###END DO

WRITE(*,*)' ***  VAC Profiles are ready  *** '

END


SUBROUTINE CO2_COR_T(W,SCO2,SCO2T,TV)  !////////////////////////////////////////
! **** 26 Apr.2025 -> T-correction in the initial 'S.CO2' .  **** !
! ***  3-9 colums in 'S.CO2'  *** !   
PARAMETER (J_S=65, I_M=7,J_T=26)
CHARACTER ATMO*30
DIMENSION S0(J_S,I_M),S_(J_S,I_M),W(J_S),COR_T(J_T,I_M),WT(J_T),TH(J_S),TV(J_S) &
  ,SCO2(J_S,J_T),SCO2T(J_S,J_T) 
DIMENSION WATM(200),TATM(200) ! <*** Carrent Atmospheric Profile *** 
DATA TH/725.00,709.00,693.00,677.00,661.00,645.00,630.00,616.00,602.00,587.00,570.50,553.50,537.00,520.50,&
503.50,486.50,469.50,453.50,438.00,423.00,409.50,397.00,385.00,372.50,357.50,340.50,322.00,301.50,&
282.50,268.50,257.50,248.00,240.00,236.50,235.00,230.50,225.50,220.00,214.00,207.00,198.50,189.00,&
180.00,173.50,169.00,165.50,163.00,162.00,162.50,166.50,167.50,162.50,157.50,152.50,147.50,142.50,&
134.50,123.00,121.00,130.00,138.50,144.00,146.50,143.50,137.00/
DATA IBG/0/
SAVE IBG,TH
IF(IBG==0)THEN
IBG=1
! ***   Reading ***
!  Initial Cross_Section !
DO J=1,J_S 
S0(J,1:7)=SCO2(J,3:9)
END DO
! --- !
OPEN(10,FILE='./F_LW/S(T)') ; OPEN(11,FILE='./F_LW/W(T)_STAND')
DO J=1,J_T 
READ(10,*)COR_T(J,:)
READ(11,*)WT(J)
END DO
CLOSE(10) ; CLOSE(11)
END IF
! ================= new Atmosphere ========================= !
! TV is defined above !
! ********** T- correction *********** !
S_=S0
DO J=36,61
RAT=(TV(J)-TH(J))/TH(J)
  DO I=1,I_M
  S_(J,I)=S0(J,I)+COR_T(J-35,I)*RAT 
  END DO ! ***DO I=1,I_M
END DO ! ***DO J=36,61

! **** OTOBRAJENIE NA SCO2T *** !
SCO2T=SCO2
DO J=36,61
SCO2T(J,3:9)=S_(J,1:7)
END DO
! *** Finish T-correction *** !  
END

SUBROUTINE NEWATM(ZA,WA,TA,JM,NG,GAS,ROCO2,ROH2O,ROSO2,ATMO)
CHARACTER HEAD*30,GAS*3,ATMO*30
DIMENSION GAS(3),ZA(200),WA(200),TA(200),ROCO2(200),ROH2O(200),ROSO2(200),RORO(200,3) 
COMMON/IGA/ICO2,IH2O,ISO2
OPEN(140,FILE='./Atmospheres/'//ATMO)
READ(140,2)HEAD ; WRITE(*,2)HEAD
2 FORMAT(A30)
READ(140,*)NG,JM
GAS='0'
ICO2=0 ; IH2O=0 ; ISO2=0
DO JG=1,NG 
READ(140,4)GAS(JG)
IF(GAS(JG)=='CO2')ICO2=JG
IF(GAS(JG)=='H2O')IH2O=JG
IF(GAS(JG)=='SO2')ISO2=JG
4 FORMAT(A3)
END DO

WA=0. ; TA=0. ; ROCO2=0. ; ROH2O=0. ; ROSO2=0.
RORO=0.
DO J=1,JM
READ(140,*)ZA(J),PA,TA(J),RORO(J,1:NG) 
 WA(J)=ALOG(PA*1013.16)
END DO
CLOSE(140)
IF(ICO2>0)ROCO2(:)=RORO(:,ICO2)
IF(IH2O>0)ROH2O(:)=RORO(:,IH2O)
IF(ISO2>0)ROSO2(:)=RORO(:,ISO2)
END
 