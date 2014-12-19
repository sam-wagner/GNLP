!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS MODULE IS USED TO DECLARE ALL THE GLOBAL VARIABLES.  THE 
!  PLANET RADII ARE IN KM AND THE GRAVITATIONAL PARAMETERS ARE IN 
!  KM^2/S^3.  AU IS IN KM.  THE J2000_ELEMENTS AND CENT_RATE ARRAYS
!  ARE TAKEN FROM ORBITAL MECHANICS FOR ENGINEERING STUDENTS.
!
MODULE CONSTANTS
IMPLICIT NONE

public
double precision, parameter :: mu=132712440018.d0
double precision, parameter :: pi=4.d0*atan(1.d0)
double precision, parameter :: rad=pi/180.d0
double precision, parameter :: AU=149599650.d0
double precision, parameter :: arc_sec=1.d0/3600.d0
double precision, parameter :: re=6378.14d0
double precision, parameter :: mue=398600.5d0
double precision :: R_planet(9), mu_planet(9) 
DOUBLE PRECISION :: J2000_ELEMENTS(9,6), CENT_RATES(9,6)
parameter(R_planet=(/2440.d0, 6052.d0, 6378.137d0, 3396.d0, & 
					71490.d0, 60270.d0, 25560.d0, 24730.d0, &
                    1195.d0/))
parameter(mu_planet= (/22032.d0, 324859.d0, 398600.4418d0, &
                       42828.d0, 126686534.d0, 37931187.d0, &
                       5793939.d0, 6836529.d0, 871.d0/))


PARAMETER(J2000_ELEMENTS=RESHAPE((/  &
  0.38709893000000001D0,       0.72333199000000004D0,   &
   1.0000001100000000D0,        1.5236623100000000D0,   &
   5.2033630100000003D0,        9.5370703199999998D0,   &
   19.191263930000002D0,        30.068963480000001D0,   &
   39.481686770000003D0,       0.20563069000000000D0,   &
  6.7732299999999999D-3,       1.6710220000000001D-2,   &
  9.3412330000000002D-2,       4.8392659999999997D-2,   &
  5.4150600000000000D-2,       4.7167710000000002D-2,   &
  8.5858700000000007D-3,       0.24880766000000001D0,   &
   7.0048700000000004D0,        3.3947099999999999D0,   &
  5.0000000000000002D-5,        1.8506100000000001D0,   &
   1.3052999999999999D0,        2.4844599999999999D0,   &
  0.76985999999999999D0,        1.7691699999999999D0,   &
   17.141749999999998D0,        48.331670000000003D0,   &
   76.680689999999998D0,       -11.260640000000000D0,   &
   49.578539999999997D0,        100.55615000000000D0,   &
   113.71504000000000D0,        74.229879999999994D0,   &
   131.72169000000000D0,        110.30347000000000D0,   &
   77.456450000000004D0,        131.53298000000001D0,   &
   102.94719000000001D0,        336.04084000000000D0,   &
   14.753850000000000D0,        92.431939999999997D0,   &
   170.96423999999999D0,        44.971350000000001D0,   &
   224.06675999999999D0,        252.25084000000001D0,   &
   181.97972999999999D0,        100.46435000000000D0,   &
   355.45332000000002D0,        34.404380000000003D0,   &
   49.944319999999998D0,        313.23218000000003D0,  &
   304.88002999999998D0,        238.92881000000000D0 /), &
    SHAPE(J2000_ELEMENTS)))

PARAMETER(CENT_RATES=RESHAPE((/&
   6.6000000000000003D-7,    9.1999999999999998D-7,   &
  -4.9999999999999998D-8,   -7.2210000000000002D-5,   &
   6.0736999999999998D-4,   -3.0152999999999998D-3,   &
   1.5202499999999999D-3,   -1.2519600000000001D-3,   &
  -7.6912000076845288D-4,    2.5270000000000000D-5,   &
  -4.9379999999999998D-5,   -3.8040000000000002D-5,   &
   1.1902000000000000D-4,   -1.2879999999999999D-4,   &
  -3.6761999999999998D-4,   -1.9149999999999999D-4,   &
   2.5140000000000000D-5,    6.4649997511878610D-5,   &
   -23.510000000000002D0,    -2.8599999999999999D0,   &
   -46.939999999999998D0,    -25.469999999999999D0,   &
   -4.1500000000000004D0,     6.1100000000000003D0,   &
   -2.0899999999999999D0,    -3.6400000000000001D0,   &
    11.069999694824219D0,    -446.30000000000001D0,   &
   -996.88999999999999D0,    -18228.250000000000D0,   &
   -1020.1900000000001D0,     1217.1700000000001D0,   &
   -1591.0500000000000D0,    -1681.4000000000001D0,   &
   -151.25000000000000D0,    -37.330001831054688D0,   &
    573.57000000000005D0,    -108.80000000000000D0,   &
    1198.2800000000000D0,     1560.7800000000000D0,   &
    839.92999999999995D0,    -1948.8900000000001D0,   &
    1312.5599999999999D0,    -844.42999999999995D0,   &
   -132.25000000000000D0,     538101628.28999996D0,   &
    210664136.06000000D0,     129597740.63000000D0,   &
    68905103.780000001D0,     10925078.350000000D0,   &
    4401052.9500000002D0,     1542547.7900000000D0,   &
    786449.20999999996D0,     522747.90625000000D0 /), &
    SHAPE(CENT_RATES)))
END MODULE CONSTANTS


MODULE MATH_MODULE
USE CONSTANTS
IMPLICIT NONE

CONTAINS
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
FUNCTION NORM(VECTOR)
DOUBLE PRECISION :: NORM  
DOUBLE PRECISION, INTENT (IN) :: VECTOR(3)
DOUBLE PRECISION :: V(3)

V=VECTOR
NORM=SQRT(V(1)*V(1)+V(2)*V(2)+V(3)*V(3))
  
END FUNCTION NORM
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
FUNCTION CROSS_PRODUCT(VEC1,VEC2)
IMPLICIT NONE

DOUBLE PRECISION :: CROSS_PRODUCT(3)
DOUBLE PRECISION, INTENT(IN) :: VEC1(3),VEC2(3)

CROSS_PRODUCT=(/ VEC1(2)*VEC2(3)-VEC2(2)*VEC1(3) ,&
                 VEC1(3)*VEC2(1)-VEC2(3)*VEC1(1) ,&
                 VEC1(1)*VEC2(2)-VEC2(1)*VEC1(2) /)

END FUNCTION CROSS_PRODUCT
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
FUNCTION ZERO_TO_2PI(THETA)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: THETA
DOUBLE PRECISION :: ZERO_TO_2PI, DUM

DUM=THETA

IF (THETA>=2.d0*PI) THEN
    DUM=THETA-FLOOR(THETA/(2.d0*PI))*2.d0*PI
ELSE IF (theta<0.d0) THEN
    DUM=THETA-(FLOOR(THETA/(2.d0*PI)))*2.d0*PI
END IF
  
ZERO_TO_2PI=DUM;

END FUNCTION ZERO_TO_2PI
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
FUNCTION ASINH1(X)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: X
DOUBLE PRECISION :: ASINH1
ASINH1=LOG(X+SQRT(X*X+1))
END FUNCTION ASINH1
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!

FUNCTION ACOTH1(X)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: X
DOUBLE PRECISION :: ACOTH1
ACOTH1=0.5D0*(LOG(1.D0+1.D0/X)-LOG(1.D0-1.D0/X))
END FUNCTION ACOTH1    
    
END MODULE MATH_MODULE




MODULE ORBITAL_FUNCTIONS
USE CONSTANTS
USE MATH_MODULE
IMPLICIT NONE



!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!    FUNCTIONS THAT ARE INCLUDED IN THIS MODEL ARE DESCRIBED IN DETAIL BELOW.  TO USE THIS M-
!        ODULE THE CONSTANTS AND MATH_MODULE MUST BE DECLARED AND COMPILED IN THE APPROPRIATE 
!        MANNER.
!
!********************************************************************************************!
!    JDATE(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
!        THIS FUNCTION TAKES THE DATE, IN INTEGER VALUES AND OUTPUTS THE JULIAN DATES
!
!        YEAR- INTEGER
!        MONTH- INTEGER
!        DAY- INTEGER
!        HOUR- INTEGER
!        MINUTE- INTEGER
!        SECOND- INTEGER
!
!********************************************************************************************!
!    GDATE(JD,YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)
!        THIS SUBROUTINE CONVERTS THE JULIAN DATE TO GREGORIAN DATE.
!
!        JD-     INPUT JULIAN DATE, DOUBLE PRECISION
!        YEAR-   OUTPUT YEAR (1900-2100), INTEGER
!        MONTH-  OUTPUT MONTH (1-12), INTEGER
!        DAY-    OUTPUT DAY (1-31), INTEGER
!        HOUR-   OUTPUT HOUR (1-23), INTEGER
!        MINUTE- OUTPUT MINUTE(1-59), INTEGER
!        SECOND- OUTPUT SECOND (1-59), INTEGER
!
!********************************************************************************************!
!    DAYS_TO_DATE(DAY_DBLE, YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
!        THIS SUBROUTINE CONVERTS THE DAYS OF THE YEAR TO MONTH, DAY, HOUR, MINUTE, SECOND
!
!        DAY_DBLE- NUMBER OF DAYS INTO THE CURRENT YEAR, DOUBLE PRECISION
!        YEAR-   INTEGER
!        MONTH-  INTEGER
!        DAY-    INTEGER
!        HOUR-   INTEGER
!        MINUTE- INTEGER
!        SECOND- INTEGER
!  
!********************************************************************************************!
!    DATESTRING(JD,DATESTRING)
!        THIS SUBROUTINE TAKES THE JULIAN DATE AND CONVERTS IT TO A STRING WITH ALL THE TIME
!            INFORMATION, MOSTLY FOR USE WHEN OUTPUTTING/WRITING DATA FILES.
!
!        JD-         JULIAN DATE, DOUBLE PRECISION
!        DATESTRING- OUTPUT STRING CONTAINING DATE INFORMATION, CHARACTER ARRAY
!
!********************************************************************************************!
!    OE2SV(a,e,i,RA,w,NU)
!        THIS FUNCTION CONVERTS THE ORBITAL ELEMENTS TO THE STATE VECTOR, RADIUS THEN VELOCI-
!            TY.  OUTPUT UNITS ARE THE SAME AS THE INPUT UNITS
!
!        a-  ORBIT SEMI-MAJOR AXIS, USUALLY IN KM OR AU
!        e-  ORBIT ECCENTRICITY
!        i-  ORBIT INCLINATION, RADIAN
!        RA- ASCENDING NODE, RADIAN
!        w-  ARGUMENT OF PERIAPSIS, RADIAN
!        NU- TRUE ANOMALY ANGLE, RADIAN 
!
!********************************************************************************************!
!    SV2OE(R,V)
!        THIS FUNCTION CONVERTS THE STATE VECTOR TO CLASSICAL ORBITAL ELEMENTS. UNITS USUALLY 
!            IN KM OR AU
!
!        R- RADIUS VECTOR OF LENGTH 3, DOUBLE PRECISION
!        V- VELOCITY VECTOR OF LENGTH 3, DOUBLE PRECISION
!
!********************************************************************************************!
!    PLAN_ELEM(JD,P,ASTEROID,N1,N2, OE)
!        THIS SUBROUTINE OUTPUTS THE ORBITAL ELEMENTS OF THE PLANET OR USER DEFINED EPHEMERIS
!            FILE, ASTEROID. THE TIME FRAME FOR THE PLANETS IS 1900-2100.  THE TIME FRAME FOR 
!            CUSTOM  BODIES IS DEPENDENT  ON THE INPUT FILE.   FOR CUSTOM BODIES THE ORBIT IS 
!            PROPOGATED,  USING  A  SOLUTION  TO  KEPLER'S  PROBLEM,  FROM  THE NEAREST POINT 
!            PREVIOUS TO THE INPUT JULIAN DATE.
!
!        JD-       JULIAN DATE CORRESPONDING TO OUTPUT ORBITAL ELEMENTS
!        P-        PLANET NUMBER OF THE PLANET WHO'S ORBITAL ELEMENTS ARE DESIRED. PLANETS 
!                      1-9, USE 10 TO USE EPHEMERIS FILE
!        ASTEROID- EPHEMERIS DATA FOR CUSTOM BODY.  IF THE PLANET NUMBER IS 9 OR LESS YOU CAN 
!                      INPUT A DUMMY ARRAY.  SIZE IS N1xN2 DOUBLE PRECISION
!        N1-       WIDTH OF ASTEROID ARRAY, SHOULD BE 13, INTEGER
!        N2-       LENGTH OF EPHEMERIS DATA, INTEGER
!        OE-       OUTPUT ARRAY OF SIZE 6 CONTAINING THE OUTPUT ORBITAL ORBITAL ELEMENTS, 
!                      DOUBLE PRECISION
!
!********************************************************************************************!
!    KEPLER(e,M,EC)
!        THIS FUNCTION COMPUTES SOLUTIONS TO KEPLER'S EQUATIONS
!            M=EC-e*SIN(EC)
!
!        e-      ORBIT ECCENTRICITY, DOUBLE PRECISION
!        M-      FINAL MEAN ANOMALY ANGLE AFTER DESIRED TIME PASSAGE IN RADIANS, DOUBLE 
!                    PRECISION
!        EC-     INITIAL GUESS FOR ECCENTRIC ANOMALY ANGLE IN RADIANS, DOUBLE PRECISION
!        KEPLER- OUTPUT ECCENTRIC ANOMALY ANGLE IN RADIANS, DOUBLE PRECISION
!
!********************************************************************************************!
!    LAMBERT_BATTIN(R1,R2,T,OT)
!        THIS FUNCTION USING BATTIN'S SOLUTION  TO  SOLVER  LAMBERT'S PROBLEM.  IT IS ASSUMED
!            THAT THE ORBIT BEING SOLVED  IS  SUN-CENTERED,  IF NOT THE VALUE OF MU SHOULD BE
!            CHANGED  IN THE  CONSTANTS MODULE.  THE  UNITS  REQUIRED  FOR  THIS FUNCTION ARE
!            STANDARD SI UNITS (KILOMETERS AND SECONDS).
!
!        R1- INITIAL ORBIT RADIUS VECTOR OF LENGTH 3 (KM), DOUBLE PRECISION
!        R2- FINAL ORBIT RADIUS VECTOR OF LENGTH 3 (KM), DOUBLE PRECISION
!        T-  TRANSFER TIME IN SECONDS, DOUBLE PRECISION
!        OT- TRANSFER ORBIT TYPE. 1-PROGRADE, 2-RETROGRADE
!
!********************************************************************************************!
!    XI_BATTIN(X)
!        THIS  FUNCTION COMPUTES  THE CONTINUED FRACTION EXPANSION FOR THE XI FUNCTION.  THIS 
!            FUNCTION WILL REALLY ONLY EVER BE USED BY THE BATTIN LAMBERT SOLUTION.
!
!        X- ARGUMENT, DOUBLE PRECISION
!
!********************************************************************************************!
!    K_BATTIN(U)
!        THIS  FUNCTION COMPUTES  THE CONTINUED  FRACTION EXPANSION  FOR THE  K  FUNCTION  IN 
!            BATTIN'S LAMBERT SOLUTION.  THIS FUNCTION WILL ONLY  EVER BE USED BY  THE BATTIN 
!            LAMBERT'S SOLUTION FUNCTION.
!
!        U- ARGUMENT, DOUBLE PRECISION
!
!********************************************************************************************!
!    FG_BATTIN(a,S,C,NU,T,R1,R2)
!        THIS FUNCTION COMPUTER THE LAGRANGE MULTIPLIER COEFFICIENTS USING INFORMATION ALREADY
!            CALCULATED IN THE LAMBERT SOLUTION FUNCTION
!
!        a- SEMI-MAJOR AXIS COMPUTED IN LAMBERT_BATTIN, DOUBLE PRECISION
!        S-
!        C-  CHORD LENGTH BETWEEN THE R1 AND R2
!        NU- ORBIT TRUE ANOMALY ANGLE IN RADIANS, DOUBLE PRECISION
!        T-  ORBIT TRANSFER TIME IN SECONDS, DOUBLE PRECISION
!        R1- INITIAL RADIUS VECTOR IN KILOMETERS, DOUBLE PRECISION
!        R2- FINAL RADIUS VECTOR IN KILOMETERS, DOUBLE PRECISION
!
!********************************************************************************************!
!    TWO_IMPULSE(J0, OT, N1, N2, PNUM, QNUM, P1, P2, dT, J1, J2, V0_INF, VF_INF)
!        THIS FUNCTION COMPUTES THE DEPARTURE AND ARRIVAL V_INF ARRAYS FOR THE STANDARD TWO
!            IMPULSE MISSION.  THE TWO DELA-V'S CAN THEM BE CALCULATED USING THE REND_POST_P-
!            ROC SUBROUTINE.
!        J0-     INITIAL STARTING JULIAN DATE, DOUBLE PRECISION
!        OT-     ORBIT TYPE, INTEGER
!        N1-     EPHEMERIS FILE LENGTH, INTEGER
!        N2-     NUMBER OF ENTRIES PER LINE IN THE EPHEMERIS FILE, INTEGER
!        PNUM-   LENGTH OF OUTPUT ARRAYS, COMPUTED BY: INT(FLOOR((JF-J0)/dT))+1, INTEGER
!                (IN THIS CASE JF IS THE FINAL JULIAN DATE FOR THE LAST LAUNCH DATE TO BE
!                 SEARCHED)
!        QNUM-   WIDTH OF OUTPUT ARRAYS, COMPUTED BY: INT(FLOOR(DBLE(M_LENGTH)/dT)), INTEGER
!        P1-     INITIAL BODY, 1-9 PLANETS, 10 EPHEMERIS FILE, INTEGER
!        P2-     FINAL BODY, INTEGER
!        dT-     SEARCH TIME STEP, DOUBLE PRECISION
!        J1-     ARRAY OF LAUNCH DATES IN JULIAN DAYS, LENGTH OF PNUM, DOUBLE PRECISION
!        J2-     ARRAY OF ARRIVAL DATES, SIZE IS PNUMxQNUM, DOUBLE PRECISION
!        V0_INF- DEPARTURE V_INF ARRAY, SIZE IS PNUMxQNUM, DOUBLE PRECISION
!        VF_INF- ARRIVAL V_INF ARRAY, SIZE IS PNUMxQNUM, DOUBLE PRECISION
!  
!********************************************************************************************!
!    RETURN_DV
!        THIS SUBROUTINE DOES ALL THE POST PROCESSING FOR SIMPLE 4 IMPULSE RETURN MISSIONS.
!
!********************************************************************************************!
!    REND_POST_PROC(P1,P2,PER_DEP, ECC_DEP, PER_ARR, ECC_ARR, V0_INF, VF_INF, J1,J2,&
!					N1,N2,DV1,DV2,DV_TOT)
!        THIS FUNCTION COMPUTES AND WRITES ALL THE OUTPUT INFORMATION THAT MAY BE NEEDED FOR
!            ROBOTIC RENDEZVOUS AND/OR DIRECT INTERCEPT MISSIONS.
!
!        P1-      DEPARTURE PLANET 1-9(USUALLY 3), INTEGER
!        P2-      ARRIVAL PLANET/BODY 1-10.  IF PLANET 2 IS 10 IT IS ASSUMED TO BE A SMALL
!                     ASTEROID, SO DV2 IS JUST THE ARRIVAL V_INF. INTEGER
!        PER_DEP- DEPARTURE PERIGEE ALTITUDE, DOUBLE PRECISION
!        ECC_DEP- DEPARTURE ECCENTRICITY, DOUBLE PRECISION
!        PER_ARR- ARRIVAL PERIGEE ALTITUDE, DOUBLE PRECISION
!        ECC_ARR- ARRIVAL ECCENTRICITY, DOUBLE PRECISION
!        V0_INF-  DEPARTURE V_INF, SIZE IS N1xN2, DOUBLE PRECISION
!        VF_INF-  ARRIVAL V_INF, SIZE IS N1xN2, DOUBLE PRECISION
!        J1-      DEPARTURE JULIAN DATES, LENGTH IS N1, DOUBLE PRECISION
!        J2-      ARRIVAL JULIAN DATES, SIZE IS N1xN2, DOUBLE PRECISION
!        N1-      NUMBER OF DEPARTURE DATES, INTEGER
!        N2-      NUMBER OF ARRIVAL COMBINATIONS, INTEGER
!        DV1-     DEPARTURE DELTA-V, SIZE N1xN2, DOUBLE PRECISION
!        DV2-     ARRIVAL DELTA-V, SIZE N1xN2, DOUBLE PRECISION
!        DV_TOT-  TOTAL DELTA-V, SIZE N1xN2, DOUBLE PRECISION
!********************************************************************************************!
!    ORBIT_PROPOGATION(R,V,T)
!        THIS FUNCTION PROPOGATES THE ORBTI OUT BY THE THE TIME T, WHICH IS IN SECONDS. R AND
!        V SHOULD BE IN SI UNITS (KILMOETERS AND KILOMETERS/SECOND).
!
!        
!
!********************************************************************************************!
!    S_STUMPFF
!
!********************************************************************************************!
!    C_STUMPFF
!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
    CONTAINS
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS FUNCTION TAKES THE DATE, IN INTEGER VALUES AND OUTPUTS THE 
!  JULIAN DATES.  THIS ALGORITHM IS VALID FROM MARCH 1, 1900 TO 
!  FEB. 28, 2100.
!
!    INPUTS:
!      YR = YEAR, INTEGER
!      M0 = MONTH, INTEGER
!      D  = DAY, INTEGER
!      M  = MONTH, INTEGER
!      S = SECOND, INTEGER
!
!    OUTPUT:
!      JDATE = JULIAN DATE, DOUBLE PRECISION
!    
FUNCTION JDATE(YR,MO,D,H,M,S)
IMPLICIT NONE

DOUBLE PRECISION :: JDATE
INTEGER, INTENT(IN) :: yr,mo,d,h,m,s

JDATE=367.d0*DBLE(YR)-FLOOR(7.d0*(DBLE(YR)+FLOOR((DBLE(MO)+&
      9.d0)/12.d0))/4.d0) + FLOOR(275.d0*DBLE(MO)/9.d0)+DBLE(D)+ &
      1.7210135d6+(((DBLE(S)/60.d0+DBLE(M)))/60.d0+dble(H))/24.d0

END FUNCTION JDATE
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
SUBROUTINE EPHEMERIS_READER(EPHEMERIS,N3,N2,NAME,CORE_NUMBER)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N3,N2, CORE_NUMBER
CHARACTER(LEN=350) :: NAME
DOUBLE PRECISION, INTENT(INOUT) :: EPHEMERIS(N3,N2)
CHARACTER(LEN=350) DUM1, ASTEROID_NAME
DOUBLE PRECISION :: JDATE, DUM
INTEGER :: N1,I, EXP, COUNT, FILE_ID

FILE_ID=CORE_NUMBER+100
N1=N3+75
ASTEROID_NAME=TRIM(ADJUSTL(NAME))//'.txt'

OPEN(FILE_ID,FILE=TRIM(ASTEROID_NAME))


COUNT=0
DO I=1,N1,1
    READ(FILE_ID,'(a)') DUM1
    IF(I.GT.35 .AND. COUNT.LT.N3)THEN
        COUNT=COUNT+1
        ! VARIABLE 1
        READ(DUM1(1:17),'(F17.9)') EPHEMERIS(COUNT,1)
        !WRITE(*,*) EPHEMERIS(COUNT,1), NAME

        !VARIABLE 2
        READ(DUM1(52:69),"(F18.16)") DUM
        READ(DUM1(71:73),"(I3)") EXP
        EPHEMERIS(COUNT,2)=DUM*10.D0**EXP

        !VARIABLE 3
        READ(DUM1(76:93),"(F18.16)") DUM
        READ(DUM1(95:97),"(I3)") EXP
        EPHEMERIS(COUNT,3)=DUM*10.D0**EXP

        !VARIABLE 4
        READ(DUM1(100:117),"(F18.16)") DUM
        READ(DUM1(119:121),"(I3)") EXP
        EPHEMERIS(COUNT,4)=DUM*10.D0**EXP
        
        !VARIABLE 5
        READ(DUM1(124:141),"(F18.16)") DUM
        READ(DUM1(143:145),"(I3)") EXP
        EPHEMERIS(COUNT,5)=DUM*10.D0**EXP        

        !VARIABLE 6
        READ(DUM1(148:165),"(F18.16)") DUM
        READ(DUM1(167:169),"(I3)") EXP
        EPHEMERIS(COUNT,6)=DUM*10.D0**EXP 
        
        !VARIABLE 7
        READ(DUM1(172:189),"(F18.16)") DUM
        READ(DUM1(191:193),"(I3)") EXP
        EPHEMERIS(COUNT,7)=DUM*10.D0**EXP 

        !VARIABLE 8
        READ(DUM1(196:213),"(F18.16)") DUM
        READ(DUM1(215:217),"(I3)") EXP
        EPHEMERIS(COUNT,8)=DUM*10.D0**EXP 

        !VARIABLE 9
        READ(DUM1(220:237),"(F18.16)") DUM
        READ(DUM1(239:241),"(I3)") EXP
        EPHEMERIS(COUNT,9)=DUM*10.D0**EXP 

        !VARIABLE 10
        READ(DUM1(244:261),"(F18.16)") DUM
        READ(DUM1(263:265),"(I3)") EXP
        EPHEMERIS(COUNT,10)=DUM*10.D0**EXP 
        
        !VARIABLE 11
        READ(DUM1(268:285),"(F18.16)") DUM
        READ(DUM1(287:289),"(I3)") EXP
        EPHEMERIS(COUNT,11)=(DUM*10.D0**EXP)
        
        !VARIABLE 12
        READ(DUM1(292:309),"(F18.16)") DUM
        READ(DUM1(311:313),"(I3)") EXP
        EPHEMERIS(COUNT,12)=DUM*10.D0**EXP         

        !VARIABLE 13
        READ(DUM1(316:333),"(F18.16)") DUM
        READ(DUM1(335:337),"(I3)") EXP
        EPHEMERIS(COUNT,13)=DUM*10.D0**EXP   
    END IF
END DO
CLOSE(FILE_ID)
END SUBROUTINE EPHEMERIS_READER


!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
SUBROUTINE GDATE(JD,YEAR,MONTH,DAY,HOUR,MINUTE,SECOND)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: JD
INTEGER, INTENT(INOUT) :: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
DOUBLE PRECISION :: DAY_DBLE, TU, TEMP
INTEGER :: LEAPYEARS

!**********************FIND THE YEAR AND DAYS OF THE YEAR HERE*********************!
TEMP=JD-2415019.5d0
TU=TEMP/365.25d0
YEAR=1900+INT(TU)
LEAPYEARS=INT((YEAR-1901)*0.25d0)
DAY_DBLE=TEMP-((YEAR-1900)*365.d0+DBLE(LEAPYEARS))

!***************CHECK TO SEE IF IT'S THE BEGINNING OF THE YEAR*********************!
IF (DAY_DBLE .LT. 1.0d0) THEN
    YEAR=YEAR-1
    LEAPYEARS=INT((YEAR-1901)*0.25d0)
    DAY_DBLE=TEMP-((YEAR-1900)*365.d0+DBLE(LEAPYEARS))
END IF

!CALL THE FUNCTIONT TO DETERMINE THE MONTH, DAY, HOUR, MINUTE, AND SECOND
CALL DAYS_TO_DATE(DAY_DBLE, YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
END SUBROUTINE GDATE
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!




!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
SUBROUTINE DAYS_TO_DATE(DAY_DBLE, YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: DAY_DBLE
INTEGER, INTENT(INOUT) :: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
DOUBLE PRECISION :: TEMP
INTEGER :: TEMP_INT, I, DAYOFYEAR, LMONTH(12)
LMONTH=31
LMONTH(2)=28
LMONTH(4)=30
LMONTH(6)=30
LMONTH(9)=30
LMONTH(11)=30
IF(MOD(YEAR,4).EQ.0) LMONTH(2)=29

!**********************DETERMINE MONTH AND DAY**********************!
DAYOFYEAR=INT(DAY_DBLE)
I=1
TEMP_INT=0

DO WHILE((DAYOFYEAR .GT. TEMP_INT+LMONTH(I)).AND.(I.LT.12))
    TEMP_INT=TEMP_INT+LMONTH(I)
    I=I+1
END DO

MONTH=I
DAY=DAYOFYEAR-TEMP_INT

!*****************DETERMINE HOURS MINUTES AND SECONDS***************!
TEMP=(DAY_DBLE-DBLE(DAYOFYEAR))*24.d0
HOUR=INT(TEMP)
TEMP=(TEMP-DBLE(HOUR))*60.d0
MINUTE=INT(TEMP)
TEMP=(TEMP-DBLE(MINUTE))*60.d0
SECOND=INT(TEMP)

END SUBROUTINE DAYS_TO_DATE
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!




!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
SUBROUTINE DATESTRING(JD, DATE_STRING)
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: JD
CHARACTER(LEN=*), INTENT(INOUT) :: DATE_STRING
INTEGER :: YR, MO, DY, HR, MN,SEC, YEAR_REDUCED
CHARACTER(LEN=3) :: MONTH_STR
CHARACTER(LEN=2) :: DAY_STR, YEAR_STR, HR_STR, MN_STR, SEC_STR
CHARACTER(LEN=1) :: DAY_STR_DUM, HR_STR_DUM, MN_STR_DUM, SEC_STR_DUM

CALL GDATE(JD,YR,MO,DY,HR,MN,SEC)

!YEAR STRING ASSIGNED HERE
YEAR_REDUCED=MOD(YR,100)
WRITE(YEAR_STR,'(I2)') YEAR_REDUCED

!MONTH STRING ASSIGNED HERE
IF (MO==1) THEN 
    MONTH_STR='JAN'
ELSEIF (MO==2) THEN
    MONTH_STR='FEB'
ELSEIF (MO==3) THEN
    MONTH_STR='MAR'
ELSEIF (MO==4) THEN
    MONTH_STR='APR'
ELSEIF (MO==5) THEN
    MONTH_STR='MAY'
ELSEIF (MO==6) THEN
    MONTH_STR='JUN'
ELSEIF (MO==7) THEN
    MONTH_STR='JUL'
ELSEIF (MO==8) THEN
    MONTH_STR='AUG'
ELSEIF (MO==9) THEN
    MONTH_STR='SEP'
ELSEIF (MO==10) THEN
    MONTH_STR='OCT'
ELSEIF (MO==11) THEN
    MONTH_STR='NOV'
ELSE
    MONTH_STR='DEC'
END IF

!DAY STRING ASSIGNED HERE
IF (DY<10.D0) THEN
    WRITE(DAY_STR_DUM,'(I1)') DY
    DAY_STR='0'//DAY_STR_DUM
ELSE
    WRITE(DAY_STR,'(I2)') DY
END IF

!HOUR STRING ASSIGNED HERE
IF (HR<10.D0) THEN
    WRITE(HR_STR_DUM,'(I1)') HR
    HR_STR='0'//HR_STR_DUM
ELSE
    WRITE(HR_STR,'(I2)') HR
END IF

!MINUTE STRING ASSIGNED HERE    
IF (MN<10.D0) THEN
    WRITE(MN_STR_DUM,'(I1)') MN
    MN_STR='0'//MN_STR_DUM
ELSE
    WRITE(MN_STR,'(I2)') MN
END IF

!SECOND STRING ASSIGNED HERE
IF (HR<10.D0) THEN
    WRITE(SEC_STR_DUM,'(I1)') SEC
    SEC_STR='0'//SEC_STR_DUM
ELSE
    WRITE(SEC_STR,'(I2)') SEC
END IF

DATE_STRING=DAY_STR//'-'//MONTH_STR//'-'//YEAR_STR//' A.D. '&
                //HR_STR//':'//MN_STR//':'//SEC_STR

END SUBROUTINE DATESTRING
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!




!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!THIS ALGORITHM CONVERTS THE ORBITAL ELEMENTS TO THE STATE VECTOR
! INPUTS ARE ALL DOUBLE PRECISION :: 
!     a = SEMI-MAJOR AXIS (KM)
!     e = ORBIT ECCENTRICITY
!     i = ORBIT ECCENTRICITY  (RAD)
!     RA= LONGITUDE OF THE ASCENDING NODE  (RAD)
!     w = ARGUMENT OF PERIAPSIS  (RAD)
!     nu= TRUE ANOMALY  (RAD)
! OUTPUT IS A DOUBLE PRECISION ARRAY WITH A LENGTH OF 6 :: 
!     SV(1:3)=R
!     SV(4:6)=V
!
! THE SUNS GRAVITATIONAL PARAMETER(mu) IS DEFINED AS A CONSTANT.  
! FOR OTHER BODIES THIS CAN BE CHANGED.
!
FUNCTION OE2SV(a,e,i,RA,w,NU)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: a,e,i,RA,w,NU
DOUBLE PRECISION :: OE2SV(6), P,H,R, R_PQW(3),V_PQW(3), &
    R_VEC(3), V_VEC(3), RMAT(3,3)

P=a*(1.d0-e*e)
H=sqrt(P*MU)
R=P/(1.d0+e*COS(NU))

R_PQW=(/ R*COS(NU), R*SIN(NU), 0.d0 /)
V_PQW=(/ -SIN(NU), e+COS(NU), 0.d0 /)
V_PQW=MU/H*V_PQW

RMAT(1,1)=COS(RA)*COS(w)-SIN(RA)*SIN(w)*COS(i)
RMAT(1,2)=-COS(RA)*SIN(w)-SIN(RA)*COS(w)*COS(i)
RMAT(1,3)=SIN(RA)*SIN(i)
RMAT(2,1)=SIN(RA)*COS(w)+COS(RA)*SIN(w)*COS(i)
RMAT(2,2)=-SIN(RA)*SIN(w)+COS(RA)*COS(w)*COS(i)
RMAT(2,3)=-COS(RA)*SIN(i)
RMAT(3,1)=SIN(w)*SIN(i)
RMAT(3,2)=COS(w)*SIN(i)
RMAT(3,3)=COS(i)

R_VEC=MATMUL(RMAT,R_PQW)
V_VEC=MATMUL(RMAT,V_PQW)

OE2SV(1:3)=R_VEC
OE2SV(4:6)=V_VEC

END FUNCTION OE2SV
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!




!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!  THIS FUNCTION CONVERTS THE CLASSICAL ORBITAL ELEMENTS TO THE EQUINOCTIAL ORBITAL ELEMENTS
!
!  LAST UPDATED: 2/18/2013
!
!
FUNCTION OE2EQ(OE)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: OE(6)
DOUBLE PRECISION :: OE2EQ(6), a,e,i,RA,w,NU
DOUBLE PRECISION :: h, k, p, q, LAMBDA, M, ECC, EPS, HA, B
EPS=1.D-8

a=OE(1)
e=OE(2)
i=OE(3)
RA=OE(4)
w=OE(5)
NU=OE(6)

h=e*SIN(w+RA)
k=e*SIN(w+RA)
p=TAN(i/2.D0)*SIN(RA)
q=TAN(I/2.D0)*COS(RA)

IF (e.LT.(1.D0-EPS)) THEN
    ! ELLIPTICAL ORBITS
    ECC=2.D0*ATAN2(SQRT(1.D0-e)*TAN(NU/2.D0),SQRT(1.D0+e))
    M=ECC-e*SIN(ECC)
ELSE IF (e.GT.(1.D0+EPS)) THEN
    !HYPERBOLIC ORBITS
    HA=2.D0*ATANH(SQRT((e-1.D0)/(e+1.D0))*TAN(NU/2.D0))
    M=e*SINH(HA)-HA
ELSE
    !PARABOLIC ORBITS
    B=TAN(NU/2.D0)
    M=B+B**3/3.D0
END IF
LAMBDA=M+w+RA

OE2EQ(1)=a
OE2EQ(2)=h
OE2EQ(3)=k
OE2EQ(4)=p
OE2EQ(5)=q
OE2EQ(6)=LAMBDA

END FUNCTION OE2EQ
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!




FUNCTION EQ2SV(EQ, MU_P)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: EQ(6), MU_P
DOUBLE PRECISION :: EQ2SV(6)
DOUBLE PRECISION :: e, w, RA, MA, ECC, THETA, L, F, R_MAG, V_MAG, P_SEMI, h_ANG_MOM
DOUBLE PRECISION :: R(3), V(3), EQU2SV(3,3), RF(3), VF(3)
DOUBLE PRECISION :: a, h, k, p, q, LAMBDA

a=EQ(1)
h=EQ(2)
k=EQ(3)
p=EQ(4)
q=EQ(5)
LAMBDA=EQ(6)


e= SQRT(H**2+k**2)
w=ATAN2(h,k)-ATAN2(p,q)
RA=ATAN2(p,q)
MA=ZERO_TO_2PI(LAMBDA- ATAN2(h,k))
ECC=kepler(e, MA);
THETA=2.D0*ATAN2(SQRT(1.D0+e)*TAN(ECC/2.D0),SQRT(1.D0-e))
L=THETA+ATAN2(h,k)
F=ECC+ATAN2(h,k)
R_MAG=a*(1.D0-e*COS(F-(w+RA)))
V_MAG=SQRT(2.D0*MU_P/R_MAG-MU_P/a)
p_SEMI=a*(1.D0-e**2)
h_ANG_MOM=SQRT(MU_P*p_SEMI)
R(1)=R_MAG*COS(L)
R(2)=R_MAG*SIN(L)
R(3)=R_MAG*0.D0
V(1)=h_ANG_MOM/P_SEMI*(-h-SIN(L))
V(1)=h_ANG_MOM/P_SEMI*(k+COS(L))
V(1)=h_ANG_MOM/P_SEMI*0.D0
EQU2SV(1,1)=1.D0/(1+p**2+q**2)*(1.D0-p**2+q**2)
EQU2SV(1,2)=1.D0/(1+p**2+q**2)*(2.D0*p*q)
EQU2SV(1,3)=1.D0/(1+p**2+q**2)*(2.D0*p)
EQU2SV(2,1)=1.D0/(1+p**2+q**2)*(2.D0*p*q)
EQU2SV(2,2)=1.D0/(1+p**2+q**2)*(1.D0+p**2-q**2)
EQU2SV(2,3)=1.D0/(1+p**2+q**2)*(-2.D0*q)
EQU2SV(3,1)=1.D0/(1+p**2+q**2)*(-2.D0*p)
EQU2SV(3,2)=1.D0/(1+p**2+q**2)*(2.D0*q)
EQU2SV(3,3)=1.D0/(1+p**2+q**2)*(1.D0-p**2-q**2)
RF=MATMUL(R,EQU2SV)
VF=MATMUL(V,EQU2SV)
EQ2SV(1:3)=RF
EQ2SV(4:6)=VF

END FUNCTION EQ2SV

!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!  THIS FUNCTION CONVERTS THE CLASSICAL THE EQUINOCTIAL ORBITAL ELEMENTS TO CLASSICAL ORBITAL
!   ELEMENTS.
!
!  LAST UPDATED: 2/18/2013
!
!
FUNCTION EQ2OE(EQ)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: EQ(6)
DOUBLE PRECISION :: EQ2OE(6), a,e,i,RA,w,NU, M
DOUBLE PRECISION :: h, k, p, q, LAMBDA, ECC, EPS
EPS=1.D-8

a=EQ(1)
h=EQ(2)
k=EQ(3)
p=EQ(4)
q=EQ(5)
LAMBDA=EQ(6)

e=SQRT(h**2+k**2)
i=2.D0*ATAN(SQRT(p**2+q**2))
RA=ATAN2(p,q)
w=ATAN2(h,k)-ATAN2(p,q)
M=LAMBDA-ATAN2(h,k)
ECC=KEPLER(M,e)

IF (e.LT.(1.D0-EPS)) THEN
    ! ELLIPTICAL ORBITS
    NU=2.D0*ATAN2(SQRT(1.D0+e)*TAN(ECC/2.D0),SQRT(1.D0-e))
ELSE IF (e.GT.(1.D0+EPS)) THEN
    !HYPERBOLIC ORBITS
    NU=2.D0*ATAN2(SQRT(e+1.D0)*TANH(M/2.D0), SQRT(e-1.D0))
ELSE
    !PARABOLIC ORBITS
    NU=2.D0*ATAN(ECC)
END IF


EQ2OE(1)=a
EQ2OE(2)=e
EQ2OE(3)=i
EQ2OE(4)=w
EQ2OE(5)=RA
EQ2OE(6)=NU

END FUNCTION EQ2OE
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!



!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
!THIS ALGORITHM CONVERTS THE ORBITAL ELEMENTS TO THE STATE VECTOR
! INPUTS ARE ALL DOUBLE PRECISION :: 
!     R(3) = RADIUS VECTOR (KM)
!     V(3) = VELOCITY VECTOR (KM/S)
! OUTPUT IS A DOUBLE PRECISION ARRAY WITH A LENGTH OF 6 :: 
!     SV2OE(1) = a = SEMI-MAJOR AXIS (KM)
!     SV2OE(2) = e = ORBIT ECCENTRICITY
!     SV2OE(3) = i = ORBIT ECCENTRICITY (RAD)
!     SV2OE(4) = RA= LONGITUDE OF THE ASCENDING NODE (RAD)
!     SV2OE(5) = w = ARGUMENT OF PERIAPSIS (RAD)
!     SV2OE(6) = nu= TRUE ANOMALY (RAD)
!
! THE SUNS GRAVITATIONAL PARAMETER IS DEFINED AS A CONSTANT.  FOR 
! OTHER BODIES THIS CAN BE CHANGED.
FUNCTION SV2OE(R,V)
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: R(3),V(3)
DOUBLE PRECISION :: SV2OE(6),RMAG,VMAG,H, HVEC(3),EVEC(3),NVEC(3), &
    a,e,i,N, RA, w, NU, ENERGY

RMAG=NORM(R); VMAG=NORM(V)

HVEC=CROSS_PRODUCT(R,V); H=NORM(HVEC)

NVEC=CROSS_PRODUCT((/0.d0,0.d0,1.d0/),HVEC)
N=NORM(NVEC)
EVEC=((VMAG*VMAG-MU/RMAG)*R-DOT_PRODUCT(R,V)*V)/MU
e=NORM(EVEC)

ENERGY=VMAG*VMAG/2.d0-MU/RMAG

a=-MU/(2*ENERGY)

IF (e==1.d0) a=2.d8

i=ACOS(HVEC(3)/h)
RA=ACOS(NVEC(1)/N)

IF (NVEC(2)<0.d0) RA=2*PI-RA

w=ACOS(DOT_PRODUCT(NVEC,EVEC)/(N*e))

IF(EVEC(3)<0.d0) w=2*PI-w

NU=ACOS(DOT_PRODUCT(EVEC,R)/(e*RMAG))

IF (DOT_PRODUCT(R,V)<0.d0) NU=2*PI-NU
  
SV2OE(1)=a
SV2OE(2)=e
SV2OE(3)=i
SV2OE(4)=RA
SV2OE(5)=w
SV2OE(6)=NU

END FUNCTION SV2OE
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!



!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS SUBROUTINE RETURNS THE ORBITAL ELEMENTS OF THE SPECIFIED 
!  PLANET.  IF A CUSTOM BODY IS USED THE EPHEMERIS DATA IS INPUT AS 
!  THE ASTEROID ARRAY.  THIS SUBROUTINE IMPLEMENTS THE PLANETARY 
!  EPHEMERIS METHOD FROM:
!  TITLE: ORBITAL MECHANICS FOR ENGINEERING STUDENTS
!  AUTHOR: HOWARD D. CURTIS
!  YEAR: 2005
!
!    INPUTS: 
!      JD = JULIAN DATE, DOUBLE PRECISION
!      P  = PLANERTY NUMBER, 1-8 PLANETS, 9 PLUTO, 10 CUSTOM BODY IN 
!           ASTEROID, INTEGER
!      ASTEROID (N1,N2) = EPHEMERIS DATA IN THE ORDER IT COMES FROM 
!           NASA'S HORIZON SYSTEM
!      N1               = LENGTH OF THE EPHEMERIS DATA, INTEGER
!      N2               = WIDTHS OF THE EPHEMRIS DATA (SHOULD BE 13), 
!                         INTEGER
!    OUTPUTS:
!      PLAN_ELEM(1) = a,  SEMI-MAJOR AXIS (KM), DOUBLE PRECISION
!      PLAN_ELEM(2) = e,  ORBIT ECCENTRICITY, DOUBLE PRECISION
!      PLAN_ELEM(3) = i,  ORBIT ECCENTRICITY (RAD), DOUBLE PRECISION
!      PLAN_ELEM(4) = RA, LONGITUDE OF THE ASCENDING NODE (RAD), 
!                     DOUBLE PRECISION
!      PLAN_ELEM(5) = w,  ARGUMENT OF PERIAPSIS (RAD), 
!                     DOUBLE PRECISION
!      PLAN_ELEM(6) = nu, TRUE ANOMALY (RAD), DOUBLE PRECISION
!
SUBROUTINE PLAN_ELEM(JD,P,ASTEROID,N1,N2, OE)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N1, N2, P
DOUBLE PRECISION, INTENT(IN) :: JD, ASTEROID(N1,N2)
DOUBLE PRECISION, INTENT(INOUT) :: OE(6)
 
DOUBLE PRECISION :: T, a, e, i, RA, w_bar, L, w, M, E_0, EA, NU
DOUBLE PRECISION :: dJ, dT, M_0, dM
INTEGER :: II

T=(JD-2451545.d0)/36525.d0
!write(*,*) jd
IF (P<10) THEN
    a=J2000_ELEMENTS(P,1)+CENT_RATES(P,1)*T
    a=a*AU
    e=J2000_ELEMENTS(P,2)+CENT_RATES(P,2)*T
    i=(J2000_ELEMENTS(P,3)+CENT_RATES(P,3)*T*ARC_SEC)*RAD
    RA=ZERO_TO_2PI((J2000_ELEMENTS(P,4)+CENT_RATES(P,4)*T*ARC_SEC)*RAD)
    w_bar=zero_to_2PI((J2000_ELEMENTS(P,5)+CENT_RATES(P,5)*&
          T*arc_sec)*rad)
    L=zero_to_2PI((J2000_ELEMENTS(P,6)+CENT_RATES(P,6)*T*ARC_SEC)*RAD)                          
    w=W_BAR-RA
    w=ZERO_TO_2PI(w)
    M=L-W_BAR
    EA=KEPLER(e,M)
    EA=ZERO_TO_2PI(EA)
    NU=((2.d0*ATAN(TAN(EA/2.d0)*SQRT((1.d0+e)/(1.d0-e)))))
    NU=ZERO_TO_2PI(NU)
ELSE
    !write(*,*) n1, n2, asteroid(2,1), asteroid(1,1)
    dJ=ASTEROID(2,1)-ASTEROID(1,1)
    II=FLOOR((JD-ASTEROID(1,1))/dJ)+1

    IF (NINT(JD)<NINT(ASTEROID(II,1)))THEN
        II=II-1
    ELSEIF (NINT(JD)>NINT(ASTEROID(II,1)))THEN
        II=II+1
    END IF
    IF (II.GT.N1) THEN
        WRITE(*,*), N1, II
    END IF
    e=ASTEROID(II,2)
    i=ASTEROID(II,4)*RAD
    RA=ASTEROID(II,5)*RAD
    w=ASTEROID(II,6)*RAD
    NU=ASTEROID(II,10)*RAD
    a=ASTEROID(II,11)*AU

    dT=JD-ASTEROID(II,1)
    dT=dT*24.d0*60.d0*60.d0

    IF (dT>1.d-8) THEN
        E_0=2.d0*ATAN2(TAN(NU/2.d0)*SQRT((1.d0 - e)),&
            SQRT((1.d0 + e)))
        E_0=ZERO_TO_2PI(E_0)
        M_0=E_0-e*SIN(E_0)
        dM=SQRT(MU/a**3)*dT
        M=M_0+dM
        EA=KEPLER(e,M)
        NU=2.d0*ATAN2(TAN(EA/2.d0)*SQRT((1.d0+e)),SQRT((1.d0-e)))
        NU=ZERO_TO_2PI(NU); 
    END IF
END IF

OE(1)=a; OE(2)=e; OE(3)=i; OE(4)=RA; OE(5)=w; OE(6)=NU

RETURN
END SUBROUTINE PLAN_ELEM
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!


!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS FUNCTION IS THE FINAL FUNCTION TO ITERATIVELY SOLVE KEPLERS 
!  PROBLEM.  IF AN ECCENTRICIC ANOMALY ANGLE ISN'T CONVERGED UPON, 
!  A VALUE OF 1.D12 IS RETURNED.  THIS ALGORITHM IS AN IMPLEMENTATION
!  OF THE KEPLER'S SOLUTION ALGORITHMS DEVELOPED IN:
!  TITLE: FUNDAMENTALS OF ASTRODYNAMICS AND APPLICATIONS, 3RD EDITION
!  AUTHOR: DAVID A. VALLADO
!  YEAR: 2007
!
!    INPUTS:
!      e = ORBIT ECCENTRICIITY, DOUBLE PRECISION 
!      M_DUM = MEAN ANOMALY ANGLE (RAD), DOUBLE PRECISION
!
!    OUTPUTS:
!      KEPLER = CONVERGED ECCENTRICIC ANOMALY ANGLE (RAD), 
!               DOUBLE PRECISION
!
FUNCTION KEPLER(e,M_DUM)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: e,M_DUM
DOUBLE PRECISION :: KEPLER, TOL, ERROR, ENEW, EA, H, HNEW, M
INTEGER :: ITER, ITERMAX

TOL=1.d-8
ERROR=1.d0
ITERMAX=100
ITER=0
M=M_DUM

IF (e.LT.1.d0) THEN
    EA=M+e*SIN(M)+e**2/2.d0*SIN(2.d0*M)
    DO WHILE (ERROR.GE.TOL .AND. ITER.LE.ITERMAX)
        ENEW=EA+(M-EA+e*SIN(EA))/(1.d0-e*COS(EA))
        ERROR=ABS(EA-ENEW)
        EA=ENEW
        ITER=ITER+1
    END DO
    
    IF (ITER.GE.ITERMAX) THEN
        KEPLER=1.D12
    ELSE
        KEPLER=EA
    END IF    
        
ELSE IF (e>=1.d0) THEN
    IF (e<1.6d0)THEN
        IF(M.LT.0.D0 .AND. M.GT.-PI)THEN
            H=M-e
        ELSEIF(M.GT.PI)THEN
            H=M-e
        ELSE
            H=M+e
        END IF    
    ELSE
        IF(e.LT.3.6D0 .AND. ABS(M).GT.PI)THEN
            H=M-M/ABS(M)*e
        ELSE
            H=M/(e-1.D0)
        END IF
    END IF
    DO WHILE(ERROR.GE.TOL .AND. ITER.LE.ITERMAX)
        HNEW=H+(M-e*SINH(H)+H)/(e*COSH(H)-1.D0)
        ERROR=ABS(H-HNEW)
        H=HNEW
        ITER=ITER+1
    END DO

    IF (ITER.GE.ITERMAX) THEN
        KEPLER=1.D12
    ELSE
        KEPLER=H
    END IF 
END IF        

END FUNCTION KEPLER

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!

! THIS ALGORITHM SOLVES KEPLER'S PROBLEM WITH UNIVERSAL VARIABLES AND OUTPUTS THE FINAL\
!   RADIUS AND VELOCITY VECTORS.  THIS APPROACH IS BASE ON "ORBITAL MECHANICS" BY CHOBOTOV 
!   AND "ORBITAL MECHANICS" BY PRUSSING AND CONWAY.  THIS ALGORITHM IS EXTREMELY ROBUST DUE 
!   TO THE LAGEURRE ROOT FINDING METHOD AND THE INITIAL GUESS ESTIMATE GIVEN IN PRUSSINGS
!   BOOK.
!
! as of right now this algorithm isn't working correctly
SUBROUTINE KEPLER_UNIVERSAL(M,e,R0,V0,RF,VF)
IMPLICIT NONE
INTEGER :: N_LAG, ITER, ITERMAX
DOUBLE PRECISION, INTENT(INOUT) :: RF(3), VF(3)
DOUBLE PRECISION, INTENT(IN) :: M,e,R0(3), V0(3)
DOUBLE PRECISION :: X, F, FP, FPP, X_UP, a, ALPHA, N, DT, DELTA, G, F_DOT, G_DOT,R
DOUBLE PRECISION :: R0_MAG, V0_MAG, RP, P, LP, Z, C, S, SIGMA0, EPSILON, DZ, X_NEW

R0_MAG=NORM(R0)
V0_MAG=NORM(V0)
a=R0_MAG*MU/(2.D0*MU-R0_MAG*V0_MAG**2)
ALPHA=1/a
n=SQRT(MU/a**3)
DT=M/n
RP=a*(1.D0-e)
SIGMA0=DOT_PRODUCT(R0,V0)/SQRT(MU)

!THIS CORRECTS DT FOR MULTIPLE REVOLUTION CASES
IF (ALPHA.GT.0.D0)THEN
P=2.D0*PI*SQRT(a**3/MU)
LP=A*(1.D0-e**2)
DT=DT-DT/ABS(DT)*INT(ABS(DT)/P)*LP
END IF
!DETERMINE INITIAL GUESS FROM PRUSSING HERE BY USING A SECANT APPROXIMATION OF F USING THE 
!UPPER BOUND OF X, WITH IS AT THE PERIGEE RADIUS
X_UP=SQRT(MU)*DT/RP
Z=X_UP**2/a
C=C_STUMPFF(Z)
S=S_STUMPFF(Z)
F=(1.D0-R0_MAG*ALPHA)*S*X_UP**3
X=MU*DT**2/(RP*(F+sqrt(MU)*DT))

EPSILON=1.D-8
DZ=1.D8
ITER=0
N_LAG=5
ITERMAX=50
! THE LAGUERRE NUMERICAL METHOD IS USED TO SOLVE FOR THE UNIVERSAL VARIABLE X. 
! IN THIS CASE N IS SET TO 5

DO WHILE(ITER.LT.ITERMAX .AND. DZ.GT.EPSILON)
    ITER=ITER+1
    Z=X**2*ALPHA
    S=S_STUMPFF(Z)
    C=C_STUMPFF(Z)
    F=(1.D0-R0_MAG*ALPHA)*S*X**3+SIGMA0*C*X**2+R0_MAG*X-SQRT(MU)*DT
    FP=C*X**2+SIGMA0*(1.D0-S*Z)*X+R0_MAG*(1.D0-C*Z)
    FPP=(1.D0-R0_MAG/a)*(1.D0-S*Z)*X+SIGMA0*(1.D0-C*Z)
    DELTA=SQRT(ABS(DBLE(N_LAG-1)**2*FP**2-DBLE(N_LAG*(N_LAG-1))*F*FPP))
    X_NEW=X-N_LAG*F/(FP+FP/ABS(FP)*DELTA)
    DZ=ABS(X_NEW**2*ALPHA-Z)
    X=X_NEW
END DO

!THE F AND G FUNCTIONS ARE USED TO FIND THE ENDING RADIUES AND VELOCITY VECTORS

Z=X**2*ALPHA
R=a+a*(DOT_PRODUCT(R0,V0)/SQRT(MU*a)*SIN(X/SQRT(a))-(1.D0-R0_MAG/a)*COS(X/SQRT(a)))
C=C_STUMPFF(Z)
S=S_STUMPFF(Z)
F=1.D0-X**2/R0_MAG*C
G=DT-X**3/SQRT(MU)*S
F_DOT=SQRT(MU)/(R*R0_MAG)*(S*Z-1.D0)*X
G_DOT=1.D0-X**2/R*C
RF=F*R0+G*V0
VF=F_DOT*R0+G_DOT*V0

END SUBROUTINE KEPLER_UNIVERSAL

!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!





!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!THIS LAMBERT SOLUTION ALGORITHM

!THIS ALGORITHM CURRENTLY DOESN'T WORK, PLEASE DON'T TRY AND USE IT.
SUBROUTINE LAMBERT_AHNLEE(R1,R2,DT,OT, V)
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: R1(3), R2(3), DT
INTEGER, INTENT(IN) :: OT
DOUBLE PRECISION, INTENT(INOUT) :: V(6)
DOUBLE PRECISION :: R1MAG, R2MAG,C12(3), TOL, NU, RATIO, GAMMA_MIN, GAMMA_MAX, GAMMA_0
DOUBLE PRECISION :: GAMMA1, PHI, F1, F2, VR1, VT1, h, v1, V2, RHO, p, e, ECC1, ECC2, a, ALPHA
DOUBLE PRECISION :: TF, n, EPS, H1, H2, NH, W1, W2, NP, DTF, TFR, GAMMA2, DTFDG1, GVR2, GVR1 
DOUBLE PRECISION :: GA1, GA2, GAN, ZI1, BGAMMA, q, EPSIL, BTHETA, ETA, BOMEGA1, DELTA_GAMMA1
DOUBLE PRECISION :: BOMEGA2, BLAMBDA, ALPHA_MIN, ALPHA_MAX, S, C, FG(3)
INTEGER :: ITERMAX, ITER


TOL=1.d-8
R1MAG=NORM(R1)
R2MAG=NORM(R2)
C12=CROSS_PRODUCT(R1,R2)
NU=ACOS(DOT_PRODUCT(R1,R2)/(R1MAG*R2MAG))

!DETERMINE THE TRUE ANOMALY ANGLE USING THE ORBIT TYPE
!1 IS PROGRADE, 2 IS RETROGRADE
IF (OT==1) THEN
    IF (C12(3)<=0.d0) NU=2.d0*PI-NU;
END IF

IF (OT==2) THEN
    IF (C12(3)>=0.d0) NU=2.d0*PI-NU;
END IF

PHI=NU
 
RATIO=R1MAG/R2MAG

IF (NU.LT.PI) THEN
    GAMMA_MIN=ATAN2(COS(NU)-R1MAG/R2MAG,SIN(NU))
ELSE
    GAMMA_MIN= -PI/2.D0
END IF

GAMMA_MAX=ATAN2(SIN(PHI)+SQRT(2.D0*R1MAG/R2MAG*(1.D0-COS(PHI))),1.D0-COS(PHI))

GAMMA_0=(GAMMA_MIN+GAMMA_MAX)/2.D0

GAMMA1=GAMMA_0

EPS=1.D-6

ALPHA_MIN=0.95D0
ALPHA_MAX=0.85D0


!WRITE(*,*) ALPHA_MIN, ALPHA_MAX

TF=0.D0
ITER=0
ITERMAX=50
DTF=1.D8
TFR=DT


DO WHILE(DTF.GT.TOL .AND. ITER.LE.ITERMAX)
    ITER=ITER+1
    RHO=TAN(GAMMA1)
    VT1=SQRT(mu*(1.D0-COS(PHI))/(R1MAG*(RATIO-COS(PHI)+RHO*SIN(PHI))))
    VR1=RHO*VT1
    V1=SQRT(VT1*VT1+VR1*VR1)
    
    h=R1MAG*VT1
    
    F1=ATAN2(h*VR1/MU,h*VT1/MU-1.D0)
    F2=F1+PHI

    a=1.D0/(2.D0/R1MAG-V1**2/MU)
    p=h**2/MU
    e=SQRT(1.D0-p/a)
    ALPHA=1.D0/a
    ZI1=R1MAG*V1**3/(2.D0*MU*(1.D0-COS(PHI)))*&
        (-SIN(PHI)*COS(2.D0*GAMMA1)+(RATIO-COS(PHI))*SIN(2.D0*GAMMA1))
    
    V2=SQRT(MU*(2.D0/R1MAG-ALPHA))
    
    ! IN THIS PART THE TOF IS CALCULATED FORM THE CURRENT GAMMA GUESS
    IF(e<1.D0-EPS) THEN ! ELLIPTICAL ORBIT
        ECC1=2.D0*ATAN2(TAN(F1/2.D0)*SQRT(1.D0-e),SQRT(1.D0+e))
        ECC2=2.D0*ATAN2(TAN(F2/2.D0)*SQRT(1.D0-e),SQRT(1.D0+e))
        n=SQRT(MU/a**3)
        TF=(ECC2-ECC1-e*(SIN(ECC2)-SIN(ECC1)))/n
    ELSE IF(e>1.D0+EPS) THEN ! HYPERBOLIC ORBIT
        H1=2.D0*ATAN2(TAN(F1/2.D0)*SQRT(e-1.D0),SQRT(1.D0+e))
        H2=2.D0*ATAN2(TAN(F2/2.D0)*SQRT(e-1.D0),SQRT(1.D0+e))
        NH=SQRT(MU/(-(a**3)))
        TF=(e*(SINH(H2)-SINH(H1))-H2+H1)/NH
    ELSE !PARABOLIC ORBITA
        W1=TAN(F1/2.D0)
        W2=TAN(F2/2.D0)
        NP=2.D0*SQRT(MU/P**3)
        TF=((W2**3/3.D0+W2)-(W1**3/3.D0+W1))/NP
    END IF

    DTF=ABS(TF-TFR)
    
    GAMMA2=ATAN((1.D0/RATIO-1.D0)/TAN(PHI/2.0)-TAN(GAMMA1)/RATIO)
    BGAMMA=-(R2MAG*COS(GAMMA2)**2)/(R1MAG*COS(GAMMA1)**2)
    
    IF (DTF.GT.TOL) THEN
        IF(e<1.D0-EPS)THEN !ELLIPTICAL ORBIT
            GAN=3.D0*N*TF/(2.D0*ALPHA)
            GVR1=(R1MAG*SQRT(ALPHA)*(COS(ECC1)-e))/(SQRT(MU)*e)
            GVR2=(R2MAG*SQRT(ALPHA)*(COS(ECC2)-e))/(SQRT(MU)*e)
            GA1=SIN(ECC1)*((COS(ECC1)-e)/(2.D0*ALPHA)+R1MAG/e)
            GA2=SIN(ECC2)*((COS(ECC2)-e)/(2.D0*ALPHA)+R2MAG/e)
            DTFDG1=1.D0/n*(GVR2*(SIN(GAMMA2)*V1/V2*ZI1+V2*COS(GAMMA2)*BGAMMA) &
                   -GVR1*(SIN(GAMMA1)*ZI1+V1*COS(GAMMA1))  &
                   +(GA2-GA1-GAN)*(-2.D0*V1/MU)*ZI1  )
            
        ELSE IF (e>1.D0+EPS) THEN !HYPERBOLIC ORBIT
            GAN=3*NH*TF/(2.D0*ALPHA)
            GVR1=R1MAG*SQRT(-ALPHA)*(e-COSH(ECC1))/(SQRT(MU)*e)
            GVR2=R2MAG*SQRT(-ALPHA)*(e-COSH(ECC2))/(SQRT(MU)*e)
            GA1=SINH(ECC1)*((e-COSH(ECC1))/(2.D0*ALPHA)-R1MAG/e)
            GA2=SINH(ECC2)*((e-COSH(ECC2))/(2.D0*ALPHA)-R2MAG/e)
            DTfDG1=1.D0/NH*(GVR2*(SIN(GAMMA2)*V1/V2*ZI1+V2*COS(GAMMA2)*BGAMMA)&
                   -GVR1*(SIN(GAMMA1)*ZI1+V1*COS(GAMMA1))&
                   +(GA2-GA1-GAN)*(-2.D0*V1/MU)*ZI1)

        ELSE !PARABOLIC ORBIT
            q=p/(1.D0+e)
            EPSIL=p*V1/MU*ZI1
            ETA=R1MAG*(ZI1*COS(GAMMA1)-V1*SIN(GAMMA1))
            BOMEGA1=(1.D0+W1**2)*((SIN(2.D0*GAMMA1-F1)-SIN(F1))/V1*ZI1+COS(2.D0*GAMMA1-F1))
            BOMEGA2=(1.D0+W2**2)*((SIN(2.D0*GAMMA2-F2)-SIN(F2))/V2*V1/V2*ZI1 &
                     +COS(2.D0*GAMMA2-F2)*BGAMMA)
            BLAMBDA=-EPSIL/2.D0
            BTHETA=(-p/4.D0*EPSIL+h/MU*ETA)
            DTFDG1=1.D0/h*(2.D0*h/q*BTHETA-ETA)*TF +&
                   2.D0*q**2/h*((1.D0+W2**2)*BOMEGA2-(1.D0+W1**2)*BOMEGA1 -&
                   2.D0*BLAMBDA*((W2**3/3.D0+W2**5/5.D0)-(W1**3/3.D0+W1**5/5.D0)))                                      

        END IF
        
        DELTA_GAMMA1=(TFR-TF)/DTFDG1
        
        IF ((GAMMA1+DELTA_GAMMA1).LE.GAMMA_MIN)THEN
            GAMMA1=ALPHA_MIN*GAMMA_MIN+(1.D0-ALPHA_MIN)*GAMMA1
        ELSE IF ((GAMMA1+DELTA_GAMMA1).GE.GAMMA_MAX) THEN
            GAMMA1=ALPHA_MAX*GAMMA_MAX+(1.D0-ALPHA_MAX)*GAMMA1
        ELSE
            GAMMA1=GAMMA1+DELTA_GAMMA1
        END IF

    END IF
END DO

!NOW THE INITIAL AND FINAL VELOCITIES MUST BE FOUND.  I JUST USED THE F AND G FUNCTIONS FROM 
!THE BATTIN SOLVER, AS IT ONLY REQUIRES THE SEMI-MAJOR AXIS, WHICH CAN BE CALCULATED FROM THE
!FLIGHT PATH ANGLE
IF (ITER.LE.ITERMAX) THEN
    RHO=TAN(GAMMA1)
    VT1=SQRT(mu*(1.D0-COS(PHI))/(R1MAG*(RATIO-COS(PHI)+RHO*SIN(PHI))))
    VR1=RHO*VT1
    V1=SQRT(VT1*VT1+VR1*VR1)
    h=R1MAG*VT1

    a=1.D0/(2.D0/R1MAG-V1**2/MU)
    
    !WRITE(*,*) "NEW SOLVER SEMI-MAJOR AXIS", a
    C=SQRT(R1MAG*R1MAG+R2MAG*R2MAG-2*R1MAG*R2MAG*COS(NU))
    S=(R1MAG+R2MAG+C)/2.d0
    FG=FG_BATTIN(a,S,C,NU,DT,R1MAG,R2MAG)
    V(1:3)=(R2-FG(1)*R1)/FG(2)
    V(4:6)=(FG(3)*R2-R1)/FG(2)
ELSE
    V=1.D8
END IF


END SUBROUTINE LAMBERT_AHNLEE

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS SUBROUTINE CONTAINS THE BATTIN LAMBERT SOLUTION METHOD, AS 
!  DESCRIBED IN HIS BOOK.  IF A SOLUTION ISN'T CONVERGED UPON THEN 
!  EACH MEMBER OF THE OUTPUT VELOCITY VECTORS ARE SET TO A VALUE OF 
!  1.D12.  THE SUN GRAVITATIONAL CONSTANT IS ASSUMED, BUT THAT CAN 
!  BE CHANGED BY MODIFYING THE CONSTANT MU PARAMETER.
!
!    REFERENCE:
!      TITLE: AN INTRODUCTION TO THE MATHEMATICS AND METHODS OF
!             ASTRODYNAMICS, REVISED EDITION
!      AUTHOR: RICHARD H. BATTIN
!      YEAR: 1999
!
!    INPUTS:
!      R1(3) = INITIAL RADIUS VECTOR (KM), DOUBLE PRECISION
!      R2(3) = FINAL RADIUS VECTOR (KM), DOUBLE PRECISION
!      DT    = REQUIRED TIME-OF-FLIGHT (SEC), DOUBLE PRECISION 
!      OT    = ORBIT TYPE, 1 FOR PROGRADE, 2 FOR RETROGRADE, INTEGER
!
!    OUTPUTS:
!      V(6) = INITIAL AND FINAL VELOCITIES OF THE DETERMINED SOLUTION

SUBROUTINE LAMBERT_BATTIN(R1,R2,DT,OT, V)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: R1(3), R2(3), DT
DOUBLE PRECISION, INTENT(INOUT) :: V(6)
INTEGER, INTENT(IN) :: OT
DOUBLE PRECISION :: R1MAG, R2MAG,C12(3), TOL, NU, C, S, EPS, LTOP
DOUBLE PRECISION :: TANSQ2W, M, a, X, V1(3), V2(3), Y,Y1
DOUBLE PRECISION :: XNEW, dX, ROP, L, ETA
DOUBLE PRECISION :: DENOM, H1, H2, B, U,K, XI, FG(3), LAM,T_P,T
INTEGER :: ITERMAX, ITER

! CONVERGENCE TOLERANCE
TOL=1.d-8

R1MAG=NORM(R1)
R2MAG=NORM(R2)

! DETERMINE TRUE ANOMALY ANGLE HERE
C12=CROSS_PRODUCT(R1,R2)
NU=ACOS(DOT_PRODUCT(R1,R2)/(R1MAG*R2MAG))

!DETERMINE THE TRUE ANOMALY ANGLE USING THE ORBIT TYPE
!1 IS PROGRADE, 2 IS RETROGRADE
IF (OT==1) THEN
    IF (C12(3)<=0.d0) NU=2.d0*PI-NU
END IF

IF (OT==2) THEN
    IF (C12(3)>=0.d0) NU=2.d0*PI-NU
END IF


C=SQRT(R1MAG*R1MAG+R2MAG*R2MAG-2*R1MAG*R2MAG*COS(NU))
S=(R1MAG+R2MAG+C)/2.d0
EPS=(R2MAG-R1MAG)/R1MAG

LAM=SQRT(R1MAG*R2MAG)*COS(NU*.5D0)/S
T=SQRT(8.D0*MU/S**3)*DT
T_P=4.D0/3.D0*(1.D0-LAM**3)
M=T**2/(1.D0+LAM)**6

TANSQ2W=(EPS*EPS*0.25d0)/(SQRT(R2MAG/R1MAG)+R2MAG/R1MAG*&
        (2.d0+SQRT(R2MAG/R1MAG)))
ROP=SQRT(R2MAG*R1MAG)*(COS(NU*0.25d0)*COS(NU*0.25d0)+tansq2w)

IF (NU<PI) THEN
    LTOP=(SIN(NU*25.d-2)*SIN(NU*25.d-2)+TANSQ2W)
    L=LTOP/(LTOP+COS(NU*5.d-1))
ELSE
    LTOP=COS(NU*25.d-2)*COS(NU*25.d-2)+TANSQ2W
    L=(LTOP-COS(NU*5.d-1))/LTOP
END IF

!INITIAL GUESS IS SET HERE
IF(T.LE.T_P)THEN
    X=0.D0
    XNEW=0.D0
ELSE
    X=L
    XNEW=L
END IF

DX=1.d0
ITER=0
ITERMAX=20
! THIS TO LOOP DOES THE SUCCESSIVE SUBSTITUTION 
DO WHILE(DX>=TOL .and. ITER<=ITERMAX)
    XI=XI_BATTIN(X)
    DENOM=(1.d0+2.d0*X+L)*(4.d0*X+XI*(3.d0+X))
    H1=(L+X)**2*(1.d0+3.d0*X+XI)/DENOM
    H2=(M*(X-L+XI))/DENOM
    B=27.d0*H2*25.d-2/(1.d0+H1)**3
    U=B/(2.d0*(SQRT(1.d0+B)+1.d0))
    K=K_BATTIN(U)
    Y=(1.d0+H1)/3.d0*(2.d0+SQRT(1.d0+B)/(1.d0+2.d0*U*K*K))
    XNEW=SQRT(((1.d0-L)/2.d0)**2+M/(Y*Y))-(1.d0+L)/2.d0
    Y1=SQRT(M/(L+X)*(1.d0+X));
    DX=ABS(X-XNEW);
    X=XNEW
    ITER=ITER+1;  
END DO


IF (ITER >= ITERMAX .or. isnan(x) ) THEN
    ! IF A SOLUTION ISN'T FOUND THE FINAL VELOCITIES ARE
    !  OUTPUT AS A LARGE NUMBER
    V1(1)=1.D12
    V1(2)=1.D12
    V1(3)=1.D12
    V2(1)=1.D12
    V2(2)=1.D12
    V2(3)=1.D12
ELSE
    a=MU*DT*DT/(16.d0*ROP*ROP*X*Y*Y)
    FG=FG_BATTIN(a,S,C,NU,DT,R1MAG,R2MAG)
    V1=(R2-FG(1)*R1)/FG(2)
    V2=(FG(3)*R2-R1)/FG(2)
END IF

V(1:3)=V1
V(4:6)=V2

END SUBROUTINE LAMBERT_BATTIN
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!



!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS FUNCTION COMPUTES THE FIRST CONTINUED FRACTION FOR THE BATTIN
!  ALGORITHM.  
!
!    INPUT:
!      X = CURRENT X ESTIMATE, DOUBLE PRECISION
!
!    OUTPUT:
!      XI_BATTIN = VALUE FOR THE CONTINUED FRACTION, DOUBLE PRECISION
FUNCTION XI_BATTIN(X)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: X
DOUBLE PRECISION :: C(20), ETA, XI_BATTIN

C=(/ &
0.25396825396825395D0, 0.25252525252525254D0, 0.25174825174825177D0, &
0.25128205128205128D0, 0.25098039215686274D0, 0.25077399380804954D0, &
0.25062656641604009D0, 0.25051759834368531D0, 0.25043478260869567D0, &
0.25037037037037035D0, 0.25031928480204341D0, 0.25027808676307006D0, &
0.25024437927663734D0, 0.25021645021645023D0, 0.25019305019305021D0, &
0.25017325017325015D0, 0.25015634771732331D0, 0.25014180374361883D0, &
0.25012919896640828D0, 0.25011820330969264D0  /)

ETA=X/(SQRT(1.d0+X)+1.d0)**2;

XI_BATTIN=8.d0*(SQRT(1.d0+X)+1.d0)/(3.d0+1.d0/ &
(5.d0+ETA+9.d0/7.d0*ETA/(1.d0+C(1)*ETA/&
(1.d0+C(2)*ETA/(1.d0+C(3)*ETA/(1.d0+C(4)*ETA/(1.d0+C(5)*ETA/&
(1.d0+C(6)*ETA/(1.d0+C(7)*ETA/(1.d0+C(8)*ETA/(1.d0+C(9)*ETA/&
(1.d0+C(10)*ETA/(1.d0+C(11)*ETA/(1.d0+C(12)*ETA/(1.d0+C(13)*ETA/&
(1.d0+C(14)*ETA/(1.d0+C(15)*ETA/(1.d0+C(16)*ETA/(1.d0+C(17)*ETA/&
(1.d0+C(18)*ETA/(1.d0+C(19)*ETA/(1.d0+C(20)*ETA))))))))))))))))))))))
END FUNCTION XI_BATTIN
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!




!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS FUNCTION COMPUTED THE K CONTINUED FRACTION FOR THE BATTIN
!  LAMBERT SOLUTION ALGORITHM.  
!
!    INPUT:
!      U = SAME AS THE U VARIABLE IN THE BATTIN DESCRIPTTION,
!          DOUBLE PRECISION
!
!    OUTPUT:
!      K_BATTIN = CONTINUED FRACTION VALUE, DOUBLE PRECISION
FUNCTION K_BATTIN(U)
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: U
DOUBLE PRECISION :: K_BATTIN, D(21)

D=(/ &
0.33333333333333331D0, 0.14814814814814814D0, 0.29629629629629628D0, &   
0.22222222222222221D0, 0.27160493827160492D0, 0.23344556677890010D0, &
0.26418026418026419D0, 0.23817663817663817D0, 0.26056644880174290D0, &
0.24079807361541108D0, 0.25842383737120578D0, 0.24246606855302508D0, &
0.25700483091787441D0, 0.24362139917695474D0, 0.25599545906059318D0, &
0.24446916326782844D0, 0.25524057782122300D0, 0.24511784511784512D0, &
0.25465465465465464D0, 0.24563024563024563D0, 0.25418664443054689D0/)

K_BATTIN= D(1)/(1.d0+D(2)*U/(1.d0+D(3)*U/(1.d0+D(4)*U/&
  (1.d0+D(5)*U/(1.d0+D(6)*U/(1.d0+D(7)*U/(1.d0+D(8)*U/&
  (1.d0+D(9)*U/(1.d0+D(10)*U/(1.d0+D(11)*U/(1.d0+D(12)*U/&
  (1.d0+D(13)*U/(1.d0+D(14)*U/(1.d0+D(15)*U/(1.d0+D(16)*U/&
  (1.d0+D(17)*U/(1.d0+D(18)*U/(1.d0+D(19)*U/(1.d0+D(20)*U/&
  (1.d0+D(21)*U))))))))))))))))))))
END FUNCTION K_BATTIN
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!




!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS FUNCTION COMPUTES THE LAGRANGE COEFFICIENT VALUES TO COMPUTE
!  THE FINAL 2 VELOCITY VECTORS.  THE SUN'S GRAVITATIONAL PARAMETER
!  IS SPECIFIED BY THE SET CONSTANT VALUE, MU.  THIS CAN BE CHANGED 
!  FOR OTHER BODIES
!
!    INPUT:
!      a =
!      S = SEMIPARAMETER, DOUBLE PRECISION
!      C = CHORD, DOUBLE PRECISION
!      NU = TRUE ANOMALY ANGLE (RAD), DOUBLE PRECISION
!      T = SCALED TIME-OF-FLIGHT PARAMETER, DOUBLE PRECISION
!      R1 = INITIAL RADIUS MAGNITUDE (KM), DOUBLE PRECISION
!      R2 = FINAL RADIUS MAGNITUDE (KM), DOUBLE PRECISION
!
!    OUTPUT:
!      FG_BATTIN(1) = F, DOUBLE PRECISION
!      FG_BATTIN(2) = G, DOUBLE PRECISION
!      FG_BATTIN(3) = G_DOT, DOUBLE PRECISION
!
FUNCTION FG_BATTIN(A,S,C,NU,T,R1,R2)
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: A, S, C, NU, T, R1, R2
DOUBLE PRECISION :: FG_BATTIN(3), SMALL_NUMBER, BE, A_MIN, T_MIN, &
    DUM, AE, DE, F, G, GDOT, AH, BH, DH

SMALL_NUMBER= 1.D-3

IF (A>SMALL_NUMBER) THEN
    BE=2.D0*ASIN((SQRT((S-C)/(2*A))))
    IF (NU>PI) BE=-BE

    A_MIN=S*5.D-1
    T_MIN=SQRT(A_MIN**3/MU)*(PI-BE+SIN(BE))
    DUM=(SQRT(S/(2.D0*A)))
    AE=2.D0*ASIN(DUM)
    IF (T>T_MIN) AE=2.D0*PI-AE

    DE=AE-BE
    F=1.D0-A/R1*(1-COS(DE))
    G=T-SQRT(A*A*A/MU)*(DE-SIN(DE))
    GDOT=1.D0-A/R2*(1.D0-COS(DE))

ELSE IF (A<-SMALL_NUMBER) THEN
    AH=2.D0*ASINH1(SQRT(S/(-2.D0*A)))
    BH=2.D0*ASINH1(SQRT((S-C)/(-2.D0*A)))
    DH=AH-BH
    F=1.D0-A/R1*(1.D0-COSH(DH))
    G=T-SQRT(-A**3/MU)*(SINH(DH)-DH)
    GDOT=1.D0-A/R2*(1.D0-COSH(DH))

ELSE
    F=0.D0
    G=0.D0
    GDOT=0.D0
END IF

FG_BATTIN(1)=F
FG_BATTIN(2)=G
FG_BATTIN(3)=GDOT

END FUNCTION FG_BATTIN
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!

!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
SUBROUTINE LAMBERT_UNIVERSAL(R1,R2,DT,OT,V)
IMPLICIT NONE
INTEGER, INTENT(IN) :: OT
DOUBLE PRECISION, INTENT(IN) :: R1(3), R2(3), DT
DOUBLE PRECISION, INTENT(INOUT) :: V(6)

INTEGER :: ITER, ITERMAX
DOUBLE PRECISION :: TOL, R1MAG, R2MAG, C12(3), NU, A, Z0, Z, ZNEW
DOUBLE PRECISION :: C, S, Y, X, TP, F, DZ, DF, D2F, DC, DS, D2C, &
    D2S, N, DY, D2Y, DX, D2X, G, GDOT, SIGN_DF

TOL=1.D-8
R1MAG=NORM(R1)
R2MAG=NORM(R2)
C12=CROSS_PRODUCT(R1,R2)
NU=ACOS(DOT_PRODUCT(R1,R2)/(R1MAG*R2MAG))
!Determine the true anomaly angle using the orbit type
!1 is prograde, 2 is retrograde
IF (OT==1) THEN
    IF (C12(3) <= 0.d0) THEN
        NU=2.D0*PI-NU;
    END IF
END IF
 
IF (OT==2) THEN
    IF (c12(3)>=0.d0) THEN
        NU=2.D0*PI-NU;
    END IF
END IF

A=SQRT(R1MAG*R2MAG)*SIN(NU)/SQRT(1.D0-COS(NU))

!COMPUTE THE PARABOLIC TOF TO DETERMINE THE INITIAL GUESS FOR Z
! AT THE PARABOLIC TOF Z=0
Z=0.0D0
C=C_STUMPFF(Z)
S=S_STUMPFF(Z)
Y=R1MAG+R2MAG-A*(1.D0-Z*S)/SQRT(C)
X=SQRT(Y/C)
TP=(X**3*S+A*SQRT(Y))/SQRT(MU)

IF ((DT-TP)<-.001D0*DT) THEN
    Z0=-PI/32.D0 ! HYPERBOLIC INITIAL GUESS
ELSE IF ((DT-TP)>0.001D0*DT) THEN
    Z0=2.d0*PI**2 ! ELLIPTICAL INITIAL GUESS
ELSE 
    Z0=0.D0 ! INITIAL GUESS FOR NEAR PARABOLIC ORBITS   
END IF

ITER=0
ITERMAX=30
DZ=1.D0

Z=Z0
N=2.D0
DO WHILE(DZ>=tol .and. iter<=itermax)
    ITER=ITER+1
    S=S_STUMPFF(Z)
    C=C_STUMPFF(Z)
    Y=R1MAG+R2MAG-A*(1.D0-Z*S)/SQRT(C)
    IF (Y<0.D0) THEN
        ITER=100
    END IF
    X=SQRT(Y/C)
    DS=1.D0/(2.D0*Z)*(C-3*S)
    DC=(1.D0-Z*S - 2.D0*C)/(2.D0*Z)
    D2S=(DC-3.D0*DS)/(2.D0*Z)-(C-3.D0*S)/(2.D0*Z**2)
    D2C=(2.D0*C+Z*S-1.D0)/(2.D0*Z**2)-(S+Z*DS+2.D0*DC)/(2.D0*z)
   
    ! OMES ALTERNATIVE FOR Y
    DY=A/4*SQRT(C)
    D2Y=A*DC/(8.D0*SQRT(C))
    
    DX=(DY/C-(Y*DC)/C**2)/(2.D0*SQRT(Y/C))
    D2X= (D2Y/C-(2.D0*DY*DC)/C**2+(2.D0*Y*DC**2)/C**3-(Y*D2C)/C**2)/&
         (2.D0*SQRT(Y/C))-(DY/C-(Y*DC)/C**2)**2/(4.D0*SQRT(Y/C)**3)
    
    F=S*X**3-SQRT(MU)*DT+A*SQRT(Y)
    DF=X**3*DS+3.D0*S*X**2*DX+A*DY/(2.D0*SQRT(Y))
    D2F=X**3*D2S-(A*DY**2)/(4.D0*SQRT(Y)**3)+(A*D2Y)/(2.D0*SQRT(Y))+&
        6.D0*X**2*DX*DS +6.D0*S*X*DX**2+3.D0*S*X**2*D2X
    
    SIGN_DF=ABS(DF)/DF
    !LAGEURRE METHOD
    !ZNEW=Z-N*F/(DF+SIGN_DF*SQRT(ABS((N-1.D0)**2*DF**2-&
    !     N*(N-1.D0)*F*D2F)))
    
    !HALLEY'S METHOD
    ZNEW=Z-2.D0*F*DF/(2.D0*DF**2-F*D2F)
    DZ=ABS(Z-ZNEW)
    Z=ZNEW
END DO

IF (ITER >= ITERMAX .or. isnan(z)) THEN
    IF (ITER == 100) THEN
        CALL LAMBERT_GOODING(R1,R2,DT,OT, V)
    ELSE
        V(1)=1.D12; V(2)=1.D12; V(3)=1.D12
        V(4)=1.D12; V(5)=1.D12; V(6)=1.D12
    END IF
ELSE
    S=S_STUMPFF(Z)
    C=C_STUMPFF(Z)
    Y=R1MAG+R2MAG-A*(1.D0-Z*S)/SQRT(C)
    F=1.D0-Y/R1MAG
    G=A*SQRT(Y/MU)  
    GDOT=1.D0-Y/R2MAG
    V(1:3)=(R2-F*R1)/G
    V(4:6)=(GDOT*R2-R1)/G
END IF

END SUBROUTINE LAMBERT_UNIVERSAL
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!



!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
! THIS FUNCTION CALCULATES THE C STUMPFF FUNCTION
!
!    INPUTS:
!      Z = CURRENT Z VARIABLE ESTIMATE, DOUBLE PRECISION
!
!    OUTPUTS:
!      C_STUMPFF = C STUMPFF FUNCTION VALUE, DOUBLE PRECISION
FUNCTION C_STUMPFF(Z)
IMPLICIT NONE

DOUBLE PRECISION :: C_STUMPFF,DUM, SMALL
DOUBLE PRECISION, INTENT(IN) :: Z

SMALL=1.d-6

    IF (Z>SMALL) THEN
        DUM=(1.d0-COS(SQRT(Z)))/Z
    ELSEIF (Z<-SMALL) THEN
        DUM=(COSH(SQRT(-Z))-1.d0)/(-Z)
    ELSE
        DUM=0.5d0    
    END IF

    C_STUMPFF=DUM

END FUNCTION C_STUMPFF

FUNCTION dC_STUMPFF(Z)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: Z
DOUBLE PRECISION :: dC_STUMPFF, DUM

DUM=1.D0/(2.D0*Z)*(1.D0-Z*S_STUMPFF(Z)-2.D0*C_STUMPFF(Z))
DC_STUMPFF=DUM
END FUNCTION dC_STUMPFF



!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
! THIS FUNCTION CALCULATES THE S STUMPFF FUNCTION
!
!    INPUTS:
!      Z = CURRENT Z VARIABLE ESTIMATE, DOUBLE PRECISION
!
!    OUTPUTS:
!      S_STUMPFF = S STUMPFF FUNCTION VALUE, DOUBLE PRECISION
FUNCTION S_STUMPFF(Z)
IMPLICIT NONE

DOUBLE PRECISION :: S_STUMPFF, DUM, SMALL
DOUBLE PRECISION, INTENT(IN) :: Z
SMALL=1.d-6

IF(Z>SMALL) THEN
    DUM=(SQRT(Z)-SIN(SQRT(Z)))/SQRT(Z)**3
ELSEIF (Z<-SMALL) THEN
    DUM=(SINH(SQRT(-Z))-SQRT(-Z))/SQRT(-Z)**3
ELSE
    DUM=1.d0/6.d0    
END IF

S_STUMPFF=DUM

END FUNCTION S_STUMPFF

! THIS FUNCTION CALCULATES THE DERIVATIVE OF THE S STUMPFF FUNCTION
!
!    INPUT:
!      Z
FUNCTION dS_STUMPFF(Z)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: Z
DOUBLE PRECISION :: dS_STUMPFF, DUM

DUM=1.D0/(2.D0*Z)*(C_STUMPFF(Z)-3.D0*S_STUMPFF(Z))
DS_STUMPFF=DUM
END FUNCTION dS_STUMPFF
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!




!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!  THIS SECTION CONTAINS THE SUBROUTINES FOR THE GOODING LAMBERT 
!    LAMBERT SOLUTION
!
!    INPUTS:
!      R1 = INITIAL RADIUS VECTOR, DOUBLE PRECISION
!      R2 = FINAL RADIUS VECTOR, DOUBLE PRECISION
!      DT = TRANSFER TIME, MUST AGREE WITH MU, DOUBLE PRECISION
!      OT = ORBIT TYPE, 1=PROGRADE, 2=RETROGRADE, INTEGER
!
!    OUTPUTS:
!      V(1:3) = INITIAL VELOCITY VECTOR, DOUBLE PRECISION
!      V(4:6) = FINAL VELOCITY VECTORS, DOUBLE PRECISION
!
SUBROUTINE LAMBERT_GOODING(R1,R2,DT,OT, V)
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: R1(3), R2(3), DT
DOUBLE PRECISION, INTENT(INOUT) :: V(6)
INTEGER, INTENT(IN) :: OT

DOUBLE PRECISION :: TOL, R1MAG, R2MAG, C12(3), NU, C12MAG, &
    TAN1U(3), TAN2U(3), C, S, T, Q, T0, R1UNIT(3), R2UNIT(3), &
    C12UNIT(3), C1, C2, C3, X0, X01, X02, X03, W_BIG, X_NEW, &
    W_LITTLE, LAM, GAMMA, V1(3), V2(3), X, SIG, RHO, Z, &
    VR1(3), VR2(3), VT1(3), VT2(3), DUM, T2, TD1, TD2, TDIFF  
INTEGER :: SIGNX, ERROR1, ITER, ITERMAX, NREV
NREV=0
ITERMAX=20
TOL=1.D-8
ERROR1=0

R1MAG=NORM(R1)
R2MAG=NORM(R2)
R1UNIT=R1/R1MAG
R2UNIT=R2/R2MAG

C12=CROSS_PRODUCT(R1, R2)
C12MAG=NORM(C12)
C12UNIT=C12/C12MAG

NU=ACOS(DOT_PRODUCT(R1,R2)/(R1MAG*R2MAG))

!DETERMINE THE TRUE ANOMALY ANGLE USING THE ORBIT TYPE
!1 IS PROGRADE, 2 IS RETROGRADE
IF (OT==1) THEN
    IF (C12(3) <= 0.D0) THEN
        NU=2.D0*PI-NU;
    END IF
END IF
IF (OT==2) THEN
    IF(C12(3)>=0.D0) THEN
        NU=2.D0*PI-NU;
    END IF
END IF

!TANGENTIAL UNIT VECTORS
TAN1U=CROSS_PRODUCT(C12UNIT, R1UNIT)
TAN2U=CROSS_PRODUCT(C12UNIT, R2UNIT)
IF (NU >=PI) THEN
   TAN1U=-TAN1U
   TAN2U=-TAN2U
END IF

!============CONSTANTS DETERMINED FORM PROBLEM  GEOMETRY===================
C=SQRT(R1MAG*R1MAG+R2MAG*R2MAG - 2.D0*R1MAG*R2MAG*COS(NU))
S=(R1MAG+R2MAG+C)/2.D0
T=SQRT(8.D0*MU/S**3)*DT
Q=SQRT(R1MAG*R2MAG)/S * COS(NU*.5D0)
IF (NU<=PI) THEN
    Q=ABS(Q)
ELSE
    Q=-ABS(Q)
END IF

!INITIAL CONDITIONS
DUM=0.D0
CALL TIME_DERIVATIVES(DUM,Q,NREV, T0,DUM,DUM)

IF ((T0-T)<0.D0) THEN
    SIGNX=-1
ELSE 
    SIGNX=1
END IF  
!=================CONSTANTS FOR INITIAL CONDITIONS SPLICING================
C1=0.5D0 
C2=0.03D0 

!INITIAL VALUES FOR SINGLE REVOLUTION CASE
X0=0.D0

IF (SIGNX==1) THEN
    X0=T0*(T0-T)/(4.D0*T)
ELSE
    X01=-(T-T0)/(T-T0+4.D0);
    X02=-SQRT((T-T0)/(T+0.5D0*T0))
    W_BIG=X01+C1*SQRT(2.D0-NU/PI)
    IF (W_BIG>=0.D0) THEN
        X03=X01
    ELSE
        W_LITTLE=SQRT(SQRT(SQRT(SQRT(-W_BIG))))
        X03=X01+W_LITTLE*(X02-X01)
    END IF
        LAM=1.D0+C1*X03*(4.D0/(4.D0-T-T0))-C2*X03**2*SQRT(1.D0+X01)
        X0=LAM*X03
END IF

IF (X0<=-1) THEN 
    ERROR1=1
END IF

!============SOLVE EACH CASE GIVEN THE INITIAL CONDITIONS ABOVE============
!FOR 0-REVOLUTION CASE THE VL'S ARE RETURNED AS NAN'S AND THE VH'S RETURN
!THE CALCULATED VELOCITIES
GAMMA=SQRT(MU*S/2.D0)
IF (ERROR1==1) THEN
    V1=1.D12
    V2=1.D12
ELSE    

!X=HALLEY(X0,T,TOL,Q,NREV)
ITER=0
X=X0

DO WHILE (ABS(TDIFF)>=TOL .AND. ITER <=ITERMAX)
    ITER=ITER+1
    CALL TIME_DERIVATIVES( DUM, Q, NREV, T2, TD1, TD2)
    TDIFF=T2-T
    X_NEW=X-2.D0*TDIFF*TD1/(2.D0*TD1*TD1-TDIFF*TD2)
    X=X_NEW
END DO  
        
IF (ITER.GE.ITERMAX .OR. ISNAN(X)) THEN
    V1=1.D12
    V2=1.D12
ELSE  
    IF (C==0) THEN
        SIG=1.D0
        RHO=0.D0
        Z=ABS(X)
    ELSE
        SIG= 2.D0*SQRT(R1MAG*R2MAG/C**2)*SIN(.5D0*NU)
        RHO= (R1MAG-R2MAG)/C
        Z=   SQRT(1.D0+Q**2*(X**2-1.D0))
    END IF
    !RADIAL COMPONENTS OF VELOCITY
    DUM= GAMMA*((Q*Z-X)-RHO*(Q*Z+X))/R1MAG
    VR1=DUM*R1UNIT
    DUM=-GAMMA*((Q*Z-X)+RHO*(Q*Z+X))/R2MAG
    VR2=DUM*R2UNIT
    !TANGENTIAL COMPONENTS OF VELOCITY
    DUM=GAMMA*SIG*(Z+Q*X)
    VT1=DUM/R1MAG*TAN1U
    VT2=DUM/R2MAG*TAN2U
    V1=VR1+VT1
    V2=VR2+VT2
    END IF
END IF
V(1:3)=V1
V(4:6)=V2
END SUBROUTINE LAMBERT_GOODING
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
! THIS FUNCTION CALCULATES THE TIME-OF-FLIGHT EQUATIONS AND THE
!  DERIVATIVES
!
!    INPUTS:
!      X = CURRENT X ESTIMATE, DOUBLE PRECISION
!      Q = LAMBERT PARAMETER, DOUBLE PRECISION
!      NREV = NUMBER OF COMPLETE REVOLUTIONS, INTEGER
!      
!    OUTPUTS:
!      T = CALCULATED TIME-OF-FLIGHT, DOUBLE PRECISION
!      TD1 = TIME-OF-FLIGHT FIRST DERIVATE, DOUBLE PRECISION
!      TD2 = TIME-OF-FLIGHT SECOND DERIVATIVE, DOUBLE PRECISION
!
SUBROUTINE TIME_DERIVATIVES(X, Q, NREV, T, TD1, TD2)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: Q
DOUBLE PRECISION, INTENT (IN OUT) :: T, TD1, TD2, X
INTEGER, INTENT(IN) :: NREV
DOUBLE PRECISION :: K, E, DELT, SIG1, DSIG1, D2SIG1, D3SIG1
DOUBLE PRECISION :: SIG2, DSIG2, D2SIG2, D3SIG2, Y, Z, F ,G, D,XX
DOUBLE PRECISION :: U, BETA, TOL

TOL=1.D-8
DELT=1.D-1
XX=X
IF (XX<-1.D0) THEN
    XX=ABS(X) - 2.D0
ELSEIF (X==-1) THEN
    XX=X+DELT    
END IF
    
K=Q*Q
E=XX*XX-1.D0
IF (XX==1) THEN
    ! EXACT PARABOLIC SOLUTION, PROBABLY NEVER ACTUALLY USED
    T  = 4.D0/3.D0*(1-Q**3)   
    TD1= 4.D0/5.D0*(Q**5 - 1.D0)
    TD2= TD1+120.D0/70.D0*(1.D0-Q**7)
ELSE IF (ABS(XX-1.D0)<TOL) THEN
    ! NEAR PARABOLIC SOLUTION, USE THE TRANS/SERIES REPRESENTATION
    CALL SIGMA( -E, SIG1, DSIG1, D2SIG1)!, D3SIG1)
    CALL SIGMA( -E*Q*Q, SIG2, DSIG2, D2SIG2)!, D3SIG2)
    T=SIG1- Q**3*SIG2
    TD1=2.D0*X*(Q**5.D0*DSIG2-DSIG1)
    TD2=TD1/X+4.D0*X*X*(D2SIG1-Q**7.D0*D2SIG2)
    
ELSE
    ! ALL CASES NOT EXACTLY PARABLOLIC OR CLOSE TO PARABOLIC
    Y=SQRT(ABS(E))
    Z=SQRT(1.D0-Q**2+Q**2*X**2)
    F=Y*(Z-Q*X)
    U=-E
    BETA=Q*Z-X
    G=(X**2-Q**2*(-E))/(X*Z-Q*(-E))
    IF (E<0) THEN
        D=ATAN2(F,G)+PI*DBLE(NREV)
    ELSE
        D=ATANH(F/G)
        !D=LOG(F+G)
    END IF
    
    T=2.D0*(X-Q*Z-D/Y)/E
    TD1=(3.D0*X*T+4.D0*Q**3*X/Z-4.D0)/U
    TD2=(3.D0*T+5.0*X*TD1+4.D0*(Q/Z)**3*(1.D0-Q**2))/U 
END IF
END SUBROUTINE TIME_DERIVATIVES


! THIS FUNCTION COMPUTES THE SIGMA AND SIGMA DERIVATIVES FOR THE 
!  GOODING LAMBERT SOLUTION
!
!    INPUT:
!      U = INITIAL U VALUE FOR THE SIGMA FUNCTION, DOUBLE PRECISION
!
!    OUTPUT:
!      SIGM = SIGMA(U), DOUBLE PRECISION
!      DSIGMA = D(SIGMA)/DU, DOUBLE PRECISION
!      D2SIGMA = D2(SIGMA)/DU2, DOUBLE PRECISION
!
SUBROUTINE sigma(u, sigm, dsigma, d2sigma)
IMPLICIT NONE
DOUBLE PRECISION , INTENT(IN) :: U
DOUBLE PRECISION, INTENT(INOUT) :: SIGM, DSIGMA, D2SIGMA
         
SIGM=-2.D0*(SQRT(U)*SQRT(1.D0-U)-ASIN(SQRT(U)))/SQRT(U)**3
DSIGMA=-(U/SQRT(1.D0-U)+(3.D0*asin(SQRT(U)))/SQRT(U)-&
       3.D0/SQRT((1.D0-U)))/U**2
D2SIGMA=(15.D0*ASIN(SQRT(U))*SQRT(1.D0-U)**3-15.D0*SQRT(U)+&
        20.D0*SQRT(U)**3-3.D0*SQRT(U)**5 )/&
        (2.D0*SQRT(U)**7*SQRT((1.D0-U))**3)
END SUBROUTINE SIGMA


!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
! THIS SUBROUTINE CONTAINS THE ALGORITHM FOR THE SUN LAMBERT 
!  SOLUTION.
!
!    INPUTS:
!      R1(3) = INITIAL RADIUS VECTOR (KM), DOUBLE PRECISION
!      R2(3) = FINAL RADIUS VECTOR (KM), DOUBLE PRECISION
!      DT    = TIME-OF-FLIGHT (SEC), DOUBLE PRECISION
!      OT    = ORBIT TYPE, 1=PROGRADE, 2=RETROGRADE, INTEGER
!
!    OUTPUTS:
!      V(1:3) = INITIAL VELOCITY VECTOR (KM/S), DOUBLE PRECISION
!      V(4:6) = FINAL VELOCITY VECTOR (KM/S), DOUBLE PRECISION
SUBROUTINE lambert_sun(R1,R2,DT,OT, V)
IMPLICIT NONE
INTEGER, INTENT(IN) :: OT
DOUBLE PRECISION, INTENT(IN) :: R1(3), R2(3), DT
DOUBLE PRECISION, INTENT(INOUT) :: V(6)
DOUBLE PRECISION :: TOL, R1MAG, R2MAG, C12(3), NU, SIGMA, C, M, N, &
    TAU, TAU_P, TAU_ME, X, XNEW, DX, Y, F, FP, FPP, N_F, V1(3), &
    V2(3), VR, VC, EC(3), ER1(3), ER2(3), sign_fp
INTEGER :: ITERMAX, ITER 
ITERMAX=20
TOL=1.D-8

R1MAG=NORM(R1)
R2MAG=NORM(R2)
C12=cross_product(R1, R2)
NU=ACOS(DOT_PRODUCT(R1,R2)/(r1mag*r2mag))
 
!Determine the true anomaly angle using the orbit type
!1 is prograde, 2 is retrograde
IF (OT==1) THEN
    IF(C12(3) <= 0.D0) THEN
        NU=2.D0*PI-NU;
   END IF
END IF
 
IF (OT==2) THEN
    IF(C12(3)>=0.D0) THEN
        NU=2.D0*PI-NU;
    END IF
END IF

C=SQRT(R1MAG*R1MAG+R2MAG*R2MAG-2*R1MAG*R2MAG*COS(NU))
M=R1MAG+R2MAG+C
N=R1MAG+R2MAG-C
SIGMA=ABS(SQRT(4.D0*R1MAG*R2MAG/M**2*COS(0.5D0*NU)**2))

IF(NU<PI-TOL) THEN
    SIGMA=SIGMA
ELSEIF (NU>PI+TOL) THEN
    SIGMA=-SIGMA
ELSE
    SIGMA=0
END IF

TAU=4.D0*DT*SQRT(MU/M**3)
TAU_P=2.D0/3.D0*(1.D0-SIGMA**3)

 
IF (TAU>TAU_P+TOL) THEN !ELLIPTICAL ORBIT
    !INITIAL VALUES FOR ELIPTICAL AND PARABOLIC ORBITS
    TAU_ME=ACOS(SIGMA)+SIGMA*SQRT(1.D0-SIGMA**2)
    IF(TAU<TAU_ME-TOL) THEN
        X=0.5D0
    ELSEIF (TAU>TAU_ME+TOL)THEN
        X=-0.5D0
    ELSE
        X=0.D0
    END IF
ELSEIF(TAU<TAU_P-TOL)THEN
    !HYPERBOLIC ORBIT
    X=1.5D0
ELSE
    !PARABOLIC ORBIT
    X=1.D0
END IF

DX=1.D0
ITER=0
N_F=4.D0

DO WHILE(DX>=TOL .AND. ITER.LE.ITERMAX)
    ITER=ITER+1
    CALL lambert_sun_F(X,SIGMA,F,FP,FPP,TAU,NU)
    SIGN_FP=-FP/ABS(FP)
    !LAGUERRE
    !XNEW=X-N_F*F/(FP+SIGN_FP*FP/ABS(FP)*SQRT(ABS((N_F-1.D0)**2*&
    !     FP**2-N_F*(N_F-1)*F*FPP)))
    
    !HALLEY
    XNEW=X-2.D0*F*FP/(2.D0*FP**2-F*FPP)
    DX=ABS(X-XNEW)
    X=XNEW
END DO

IF (ITER >= ITERMAX .or. isnan(x)) THEN
    V1(1)=1.D12; V1(2)=1.D12; V1(3)=1.D12
    V2(1)=1.D12; V2(2)=1.D12; V2(3)=1.D12
ELSE
    IF(NU<PI-TOL) THEN
        Y=SQRT(1.D0-SIGMA*SIGMA*(1.D0-X*X))
    ELSEIF (NU>PI+TOL)THEN
        Y=-SQRT(1.D0-SIGMA*SIGMA*(1.D0-X*X))
    ELSE
        Y=1
    END IF

    VC=SQRT(MU)*(Y/SQRT(N)+X/SQRT(M))
    VR=SQRT(MU)*(Y/SQRT(N)-X/SQRT(M))
    EC=(R2-R1)/C
    ER1=R1/R1MAG
    ER2=R2/R2MAG
    V1=VC*EC+VR*ER1
    V2=VC*EC-VR*ER2
END IF

V(1:3)=V1
V(4:6)=V2

END SUBROUTINE lambert_sun
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!


!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
! THIS SUBROUTINE CALCULATES THE TIME-OF-FLIGHT EQUATION AND THE 
!  DERIVATIVES.
!
!    INPUTS:
!      X     = CURRENT X VALUE, DOUBLE PRECISION
!      SIGMA = LAMBERT PARAMETER, DOUBLE PRECISION
!      TAU   = SPECIFIED TIME PARAMETER, DOUBLE PRECISION
!      NU    = TRUE ANOMALY (RAD), DOUBLE PRECISION
!
!    OUTPUTS:
!      F   = TIME-OF-FLIGHT EQUATION, DOUBLE PRECISION
!      FP  = FIRST DER.E OF THE TIME EQUATION, DOUBLE PRECISION
!      FPP = SECOND DER. OF THE TIME EQUATION, DOUBLE PRECISION
!
SUBROUTINE lambert_sun_F(X, SIGMA, F, FP, FPP, TAU, NU)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: SIGMA, TAU, NU
DOUBLE PRECISION, INTENT(INOUT)  :: X, F, FP, FPP
DOUBLE PRECISION  :: TOL, DUM1,DUM2,Y, COTX, COTY
TOL=1.D-12

IF(NU<PI-TOL) THEN
    Y=SQRT(1.D0-SIGMA*SIGMA*(1.D0-X*X))
ELSEIF (NU>PI+TOL)THEN
    Y=-SQRT(1.D0-SIGMA*SIGMA*(1.D0-X*X))
ELSE
    Y=1.d0
END IF

IF (X.LE.(1.D0-TOL))THEN
    COTX=DBLE(ACOS(X))
    COTY=ATAN(SQRT(1.D0-Y*Y)/Y)
    F=1.D0/SQRT((1.D0-X*X)**3)*(COTX-COTY-X*SQRT(1.D0-X*X)+&
      Y*SQRT(1.D0-Y*Y))
ELSE IF (X.GE.(1.D0+TOL)) THEN
    DUM1=X/SQRT(X*X-1.D0)
    DUM2=Y/SQRT(Y*Y-1.D0)
    F=1.D0/SQRT((X*X-1.D0)**3)*(-ACOTH1(DUM1)+ACOTH1(DUM2)+&
      X*SQRT(X*X-1.D0)-Y*SQRT(Y*Y-1.D0))    
ELSE
    F=2.D0/3.D0*(1.D0-SIGMA**3)
END IF
    
FP=1.D0/(1.D0-X*X)*(3.D0*X*F-2.D0*(1.D0-SIGMA**3*X/ABS(Y)))
FPP=1.D0/(X*(1.D0-X*X))*((1.D0+4.D0*X*X)*FP+2.D0*&
    (1.D0-SIGMA**5*X**3/ABS(Y)**3))
F=F-TAU
END SUBROUTINE lambert_sun_F
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!




!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!*******************************************************************!
!SUBROUTINE TWO_IMPULSE(J0, OT, N1, N2, PNUM, QNUM, P1, P2, dT, ASTEROID, J1, J2, V0_INF,&
!                       NREV,SOL_TYPE)
SUBROUTINE TWO_IMPULSE(J0, OT, N1, N2, PNUM, QNUM, P1, P2, dT, ASTEROID, V0_INF,&
                       NREV, iter)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N1, N2, PNUM, QNUM, P1, P2, OT, NREV
integer, intent(inout) :: iter(pnum, qnum)
DOUBLE PRECISION, INTENT(IN) :: J0, dT, ASTEROID(N1,N2)
!DOUBLE PRECISION, INTENT(INOUT) :: J1(PNUM), J2(PNUM,QNUM) 
DOUBLE PRECISION, INTENT(INOUT) ::V0_INF(PNUM,QNUM)!, VF_INF(PNUM,QNUM)
DOUBLE PRECISION :: OE1(6), OE2(6), SV1(6), SV2(6), V(6), T, j1, j2
INTEGER :: P, Q, CON_COUNT, NON_CON_COUNT, iteration
CON_COUNT=0
NON_CON_COUNT=0
!$OMP PARALLEL private(oe1,sv1,oe2,sv2,t,V,j1,j2, iteration) 
!$OMP DO
DO  P=1,PNUM,1
    J1=J0+DBLE(P-1)*dT
    CALL PLAN_ELEM(J1,P1,ASTEROID,N1,N2,OE1)
    SV1=OE2SV(OE1(1),OE1(2),OE1(3),OE1(4),OE1(5),OE1(6))

    DO Q=1,QNUM,1
        J2=J1+DBLE(Q)*dT
        CALL PLAN_ELEM(J2,P2,ASTEROID,N1,N2,OE2)
        SV2=OE2SV(OE2(1),OE2(2),OE2(3),OE2(4),OE2(5),OE2(6))
        T=(J2-J1)*86400.d0
        
        iteration=0
        CALL LAMBERT_BATTIN(SV1(1:3),SV2(1:3),T,OT, V)
        
        !CALL lambert_sun(SV1(1:3), SV2(1:3),T,OT,V)
        !CALL LAMBERT_AHNLEE(SV1(1:3), SV2(1:3),T,OT,V)
        !CALL LAMBERT_GOODING(SV1(1:3),SV2(1:3),T,OT, V)
        !CALL lambert_UNIVERSAL(SV1(1:3), SV2(1:3),T,OT,V)
        !CALL LAMBERT_GOODING2(SV1(1:3),SV2(1:3),T,OT,NREV, V)
        iter(p,q)=iteration
        !IF (ABS(V(1)-1.D8).GT.1.D0) THEN
        !    CON_COUNT=CON_COUNT+1
        !ELSE
        !    NON_CON_COUNT=NON_CON_COUNT+1        
        !END IF
        IF (abs(V(1)-1.d12).LE.1.D-8) THEN
            V0_INF(P,Q)=1.d12
        ELSE
            V0_INF(p,q)=NORM(V(1:3)-SV1(4:6))
        END IF
        !VF_INF(p,q)=NORM(V(4:6)-SV2(4:6))
     END DO
    !WRITE(*,*) P,"/",PNUM
END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
!WRITE(*,*) "CON COUNT", CON_COUNT
!WRITE(*,*) "NON CON COUNT", NON_CON_COUNT
!WRITE(*,*) "PERCENT CONVERGED",  DBLE(CON_COUNT)/DBLE(CON_COUNT+NON_CON_COUNT)*100.D0
!WRITE(*,*) "PERCENT NON-CONVERGED", DBLE(NON_CON_COUNT)/DBLE(CON_COUNT+NON_CON_COUNT)*100.D0

END SUBROUTINE TWO_IMPULSE
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!




!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
FUNCTION RETURN_DV(LENGTH,NL,NS,EARTH_OA,MAX_EV,V1_INF, V2_INF, V3_INF, V4_INF, J0, dT, QNUM,&
                   C3_LIM)

INTEGER, INTENT(IN) :: LENGTH, NL, NS, QNUM
DOUBLE PRECISION :: RETURN_DV(NL,13), MINIMUM, C3
DOUBLE PRECISION, INTENT(IN) :: EARTH_OA, MAX_EV,J0,dT, C3_LIM, V1_INF(NL,QNUM)
DOUBLE PRECISION, INTENT(IN) :: V3_INF(NL,QNUM),V4_INF(NL,QNUM),V2_INF(NL,QNUM)
DOUBLE PRECISION :: VC,DUM,DV1(NL,QNUM),DV2(NL,QNUM),DV3(NL,QNUM), DV_TOT(NL,QNUM), J1(NL)
DOUBLE PRECISION :: DV4(NL,QNUM), J2(NL), J3(NL), J4(NL), DV1_MIN(NL), DV2_MIN(NL)
DOUBLE PRECISION :: DV3_MIN(NL), DV4_MIN(NL), C3_1(NL), C3_2(NL), C3_3(NL), C3_4(NL)
DOUBLE PRECISION :: DV_TOT_MIN(NL), VP4
INTEGER :: P, Q, I(NL), INDEX

VC=SQRT(MUE/(RE+EARTH_OA))
DUM=2.d0*MUE/(RE+EARTH_OA)
DV2=V2_INF
DV3=V3_INF

DO P=1,NL,1
    J1(P)=J0+DBLE(P-1)*DT
    J4(p)=J1(p) + DBLE(LENGTH)
    
    DO Q=1,QNUM,1

    dV1(P,Q)=SQRT(V1_INF(P,Q)*V1_INF(P,Q)+DUM)-VC
    VP4=sqrt(V4_INF(P,Q)*V4_INF(P,Q)+2.d0*MUE/(RE+100.d0))
    dV4(P,Q)=VP4-MAX_EV
        
    IF (dV4(P,Q)<0.d0) dV4(P,Q)=0.d0
        
    END DO
END DO

dV_TOT=dV1+dV2+dV3+dV4

DO P=1,NL,1
    MINIMUM=1.d8
    DO Q=1,QNUM,1
        C3=V1_INF(P,Q)*V1_INF(P,Q)
        IF (C3>C3_LIM) dV_TOT(P,Q)=5.d8
        
        IF (dV_TOT(P,Q) < MINIMUM) THEN
            MINIMUM=dV_TOT(P,Q)
            INDEX=Q
        END IF
    END DO
    I(P)=INDEX
END DO  

DO P=1,NL,1
    dV1_MIN(P)=dV1(P,I(P))
    dV2_MIN(P)=dV2(p,I(P))
    dV3_MIN(P)=dV3(p,I(P))
    dV4_MIN(P)=dV4(p,I(P))
    dV_TOT_MIN(P)=dV_TOT(p,I(P))
    J2(P)=J1(P)+DBLE(I(P))*dT
    J3(P)=J2(P)+DBLE(NS)
    C3_1(P)=V1_INF(P,I(P))*V1_INF(P,I(P))
    C3_2(P)=V2_INF(P,I(P))*V2_INF(P,I(P))
    C3_3(P)=V3_INF(P,I(P))*V3_INF(p,I(P))
    C3_4(P)=V4_INF(P,I(P))*V4_INF(p,I(P))
END DO

RETURN_DV(1:NL,1)=J1
RETURN_DV(1:NL,2)=J2
RETURN_DV(1:NL,3)=J3
RETURN_DV(1:NL,4)=J4
RETURN_DV(1:NL,5)=dV1_MIN
RETURN_DV(1:NL,6)=dV2_MIN
RETURN_DV(1:NL,7)=dV3_MIN
RETURN_DV(1:NL,8)=dV4_MIN
RETURN_DV(1:NL,9)=dV_TOT_MIN
RETURN_DV(1:NL,10)=C3_1
RETURN_DV(1:NL,11)=C3_2
RETURN_DV(1:NL,12)=C3_3
RETURN_DV(1:NL,13)=C3_4
        
END FUNCTION RETURN_DV
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!




!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
SUBROUTINE REND_POST_PROC(P1, P2, PER_DEP, ECC_DEP, PER_ARR, ECC_ARR, C3_MAX, V0_INF, &
                          VF_INF, J1, J2, J_LIM, N1, N2, DV1, DV2, DV_TOT, MIN_DATA)
IMPLICIT NONE
INTEGER, INTENT(IN) :: P1,P2,N1,N2
DOUBLE PRECISION, INTENT(IN) :: PER_DEP, ECC_DEP, PER_ARR, ECC_ARR, V0_INF(N1,N2)
DOUBLE PRECISION, INTENT(IN) :: VF_INF(N1,N2), J1(N1), J2(N1,N2), J_LIM, C3_MAX
DOUBLE PRECISION, INTENT(INOUT) :: DV1(N1,N2), DV2(N1,N2), DV_TOT(N1,N2), MIN_DATA(N1,7)
DOUBLE PRECISION :: C3_DEP(N1,N2), C3_ARR(N1,N2), VP_DEP, VP_ARR, a_ARR, a_DEP
DOUBLE PRECISION :: DV_MIN_TEMP, RP_DEP, RP_ARR, ARG, XX
INTEGER :: MIN_INDEX, P, Q

C3_DEP=V0_INF**2
C3_ARR=VF_INF**2
! CALCULATE DEPARTURE ORBIT INFORMATION
RP_DEP=R_PLANET(P1)+PER_DEP
a_DEP=RP_DEP/(1.d0-ECC_DEP)
VP_DEP=SQRT(2.d0*MU_PLANET(P1)/RP_DEP-MU_PLANET(P1)/a_DEP)

! CALCULATE ARRIVAL ORBIT INFORMATION
DV1=ABS(SQRT((C3_DEP)+(2.D0*MU_PLANET(P1)/RP_DEP))-VP_DEP)

IF (P2==10) THEN
    DV2=VF_INF
ELSE
    RP_ARR=R_PLANET(P2)+PER_ARR
    a_ARR=RP_ARR/(1.d0-ECC_ARR)
    VP_ARR=SQRT(2.D0*MU_PLANET(P2)/RP_ARR-MU_PLANET(P2)/a_ARR)
    DV2=ABS(SQRT((C3_ARR)+(2.D0*MU_PLANET(P2)/RP_ARR))-VP_ARR)
END IF
DV_TOT=DV1+DV2

ARG=-1.d0
XX=SQRT(ARG)
DO P=1,N1,1
    DO Q=1,N2,1
        IF (J2(P,Q) .GE. J_LIM) THEN
            DV1(P,Q)=XX
            DV2(P,Q)=XX
            DV_TOT(P,Q)=XX
        END IF
        IF(C3_DEP(P,Q) .GE. C3_MAX) THEN
            DV1(P,Q)=XX
            DV2(P,Q)=XX
            DV_TOT(P,Q)=XX
        END IF
    END DO
END DO


DO P=1,N1,1
    DV_MIN_TEMP=1.D8
    DO Q=1,N2,1
        IF(DV_TOT(P,Q)<DV_MIN_TEMP)THEN
            DV_MIN_TEMP=DV_TOT(P,Q)
            MIN_INDEX=Q
        END IF
    END DO
    MIN_DATA(P,1)=J1(P)
    MIN_DATA(P,2)=J2(P,MIN_INDEX)
    MIN_DATA(P,3)=DV1(P,MIN_INDEX)
    MIN_DATA(P,4)=DV2(P,MIN_INDEX)
    MIN_DATA(P,5)=DV_TOT(P,MIN_INDEX)
    MIN_DATA(P,6)=V0_INF(P,MIN_INDEX)
    MIN_DATA(P,7)=VF_INF(P,MIN_INDEX)
END DO  

END SUBROUTINE REND_POST_PROC                    
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!




!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
FUNCTION ORBIT_PROPOGATION(R,V,T)
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: R(3),V(3),T
DOUBLE PRECISION :: OE(6), E0, M0, dM, M, E, NU0, NUF, ECC, ORBIT_PROPOGATION(6)
DOUBLE PRECISION  :: N0, dN, N, H0,H
OE=SV2OE(R,V)
IF (T==0) THEN
    ORBIT_PROPOGATION(1:3)=R
    ORBIT_PROPOGATION(4:6)=V
ELSE
    IF (OE(1)>0.D0)THEN    
        NU0=OE(6)
        ECC=OE(2)
        E0=2.d0*ATAN2(TAN(NU0/2.d0)*SQRT(1.d0-ECC),SQRT(1.d0+ECC))
        M0=E0-ECC*SIN(E0)
        dM=sqrt(MU/OE(1)**3)*T
        M=M0+dM
        E=KEPLER(ECC,M)
        NUF=2*ATAN2(TAN(E/2.d0)*SQRT(1.d0+ECC),SQRT(1.d0-ECC))
        ORBIT_PROPOGATION=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),NUF)
    ELSE
        NU0=OE(6)
        ECC=OE(2)
        H0=2.D0*ATANH(TAN(NU0/2.D0)*SQRT((ECC-1.D0)/(ECC+1.D0)))
        N0=ECC*SINH(H0)-H0
        dN=SQRT(MU/(-OE(1)**3))*T
        N=N0+dN
        H=KEPLER(ECC,N)
        NUF=2.D0*ATAN2(TANH(H/2.D0)*SQRT(ECC+1.D0),SQRT(ECC-1.D0))
        ORBIT_PROPOGATION=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),NUF)
    END IF
END IF
END FUNCTION ORBIT_PROPOGATION
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!





!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
SUBROUTINE GA_MGA(V_INF_IN, V_INF_OUT, P_NUM, R_P, DV_GA, G)
IMPLICIT NONE
INTEGER, INTENT(IN)  ::  P_NUM
DOUBLE PRECISION, INTENT(IN)  ::  V_INF_IN(3), V_INF_OUT(3)
DOUBLE PRECISION, INTENT(INOUT)  ::  R_P, DV_GA, G
DOUBLE PRECISION  ::  V_INF_IN_MAG, V_INF_OUT_MAG, a1, a2, d
DOUBLE PRECISION  :: F, dF, e2_NEW, e2, a_R, de2, TOL
INTEGER :: ITER, ITERMAX
!write(*,*) "ga_mga called"
ITERMAX=20
TOL=1.D-8

V_INF_IN_MAG=NORM(V_INF_IN)
V_INF_OUT_MAG=NORM(V_INF_OUT)
a1=-MU_PLANET(P_NUM)/V_INF_IN_MAG**2
a2=-MU_PLANET(P_NUM)/V_INF_OUT_MAG**2
a_R=a2/a1
d=ACOS(DOT_PRODUCT(V_INF_IN, V_INF_OUT)/(V_INF_IN_MAG*V_INF_OUT_MAG))
!NEWTON SOLVER TO DETERMINE e2 AND IN TURN R_P
e2=1.5D0
de2=1.D0
ITER=0

DO WHILE(ITER.LE.ITERMAX .AND. de2.GE.TOL)
    ITER=ITER+1
    F=(a_R*(e2-1.D0)+1.D0)*SIN(d-ASIN(1.D0/e2))-1.D0
    dF=(a_R*e2-a_R+1.D0)*COS(d-ASIN(1.D0/e2))/(e2**2*SQRT(1.D0-1.D0/e2**2))+&
        a_R*SIN(d-ASIN(1.D0/e2))
    e2_NEW=e2-F/dF
    de2=ABS(e2_NEW-e2)
    e2=e2_NEW
END DO
!write(*,*) "e2", ITER, e2
!IF (ISNAN(E2)) THEN
!    WRITE(*,*) V_INF_IN_MAG
!    WRITE(*,*) V_INF_OUT_MAG
!END IF
IF(ITER.LE.ITERMAX) THEN
    R_P=a2*(1.D0-e2)
    dV_GA=ABS(SQRT(V_INF_IN_MAG**2+2.D0*MU_PLANET(P_NUM)/R_P)-SQRT(V_INF_OUT_MAG**2+2.D0*MU_PLANET(P_NUM)/R_P))   
    !Penalty function
    IF (R_P<(R_PLANET(P_NUM)*1.047))THEN
        G=-2.D0*LOG10(R_P/(1.047*R_PLANET(P_NUM)))
    ELSE
        G=0
    END IF
ELSE
    R_P=0
    dV_GA=1.D8
    G=0.D0
END IF

IF (ISNAN(DV_GA) .OR. ISNAN(E2) .OR. ISNAN(R_P))THEN
    DV_GA=1.D8
    G=1.D8
END IF

END SUBROUTINE GA_MGA
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!




!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
SUBROUTINE RESONANT_ORBIT(T1,T2,P1,P2, V1_INF_IN,V2_INF_OUT, V1_INF_OUT, V2_INF_IN)
IMPLICIT NONE
INTEGER, INTENT(IN) :: P1, P2
DOUBLE PRECISION, INTENT(IN) :: V1_INF_IN(3), V2_INF_OUT(3), T1, T2
DOUBLE PRECISION, INTENT(INOUT) :: V1_INF_OUT(3), V2_INF_IN(3)
DOUBLE PRECISION :: a, OE(6), V_SC_SUN, V_INF_MAG, VE1(3), RE1(3), VE2(3), RE2(3), SV(6)
DOUBLE PRECISION :: ASTEROID(1,1), THETA, PHI, dPHI, G
DOUBLE PRECISION :: V_HAT(3), N_HAT(3), C_HAT(3), T_ECL_VNC(3,3) 
DOUBLE PRECISION :: T_VNC_ECL(3,3), VE1_MAG, RE1_MAG,P, R_P_TOT, R_P_MAX,R_P1_F, R_P2_F
DOUBLE PRECISION :: V1_INF_OUT_F(3), V2_INF_IN_F(3), R_P1, R_P2, ARG, XX, DV_DUM
INTEGER :: ITER, PNUM, N1, N2, TEST

!WRITE(*,*) V1_INF_IN
!WRITE(*,*) V2_INF_OUT
ASTEROID=0.D0
N1=1
N2=1

P=(T2-T1)*86400.D0
a=((P/(2.D0*PI))**2*MU)**(1.D0/3.D0)

CALL PLAN_ELEM(T1, P1,ASTEROID,N1,N2,OE)
SV=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))

RE1=SV(1:3)
VE1=SV(4:6)

CALL PLAN_ELEM(T2,P2,ASTEROID,N1,N2,OE)
SV=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))
RE2=SV(1:3)
VE2=SV(4:6)

RE1_MAG=NORM(RE1)
V_SC_SUN=SQRT(MU*(2.D0/RE1_MAG-1.D0/a))
V_INF_MAG=NORM(V1_INF_IN)
VE1_MAG=NORM(VE1)
THETA=ACOS((V_INF_MAG**2+VE1_MAG**2-V_SC_SUN**2)/(2.D0*V_INF_MAG*VE1_MAG))

V_HAT=VE1/VE1_MAG
N_HAT=CROSS_PRODUCT(RE1,VE1)/NORM(CROSS_PRODUCT(RE1,VE1))
C_HAT=CROSS_PRODUCT(V_HAT, N_HAT)
T_ECL_VNC(1:3,1)=V_HAT
T_ECL_VNC(1:3,2)=N_HAT
T_ECL_VNC(1:3,3)=C_HAT
T_VNC_ECL(1,1:3)=V_HAT
T_VNC_ECL(2,1:3)=N_HAT
T_VNC_ECL(3,1:3)=C_HAT

dPHI=0.01
PHI=0.D0
PNUM=FLOOR(2.D0*PI/dPHI)
ITER=0
TEST=0

!WRITE(*,*) THETA
R_P_TOT=0.D0
R_P_MAX=0.D0
DO WHILE(ITER.LE.PNUM)
    V1_INF_OUT(1)=V_INF_MAG*COS(PI-THETA)
    V1_INF_OUT(2)=V_INF_MAG*SIN(PI-THETA)*COS(PHI)
    V1_INF_OUT(3)=V_INF_MAG*(-SIN(PI-THETA)*SIN(PHI))
    V1_INF_OUT=MATMUL(T_ECL_VNC,V1_INF_OUT)
    V2_INF_IN=V1_INF_OUT+VE1-VE2
    !CHECK THE FIRST GRAVITY ASSIST PERIGEE RADIUS FIRST
    CALL GA_MGA(V1_INF_IN, V1_INF_OUT, P1, R_P1, DV_DUM, g)

!    WRITE(*,*) V1_INF_OUT

    CALL GA_MGA(V2_INF_IN, V2_INF_OUT, P2, R_P2, DV_DUM, g)
    IF(R_P1 .GE. R_PLANET(P1)) THEN
        IF (R_P2 .GE. R_PLANET(P2)) THEN
            TEST=1
            R_P_TOT=R_P1+R_P2-2*R_PLANET(P1)
            IF (R_P_TOT.GT.R_P_MAX) THEN
                R_P_MAX=R_P_TOT
                R_P1_F=R_P1
                R_P2_F=R_P2
                V1_INF_OUT_F=V1_INF_OUT
                V2_INF_IN_F=V2_INF_IN            
            END IF
        END IF
    END IF

!    write(*,*) R_P1, R_P2, ITER

    PHI=PHI+dPHI
    ITER=ITER+1
END DO
R_P1=R_P1_F
R_P2=R_P2_F
IF(TEST.EQ.1) THEN
    V1_INF_OUT=V1_INF_OUT_F
    V2_INF_IN=V2_INF_IN_F
ELSE
    
    
    V1_INF_OUT=1.D12
    V2_INF_IN=1.D12
END IF

!WRITE(*,*) "a        =", a/AU
!WRITE(*,*) "V_SC_SUN =", V_SC_SUN
!WRITE(*,*) "V_INF_MAG=", V_INF_MAG
!WRITE(*,*) "VE1_MAG  =", VE1_MAG
!WRITE(*,*) "THETA    =", THETA
!WRITE(*,*) "PHI      =", PHI
!WRITE(*,*) "R_P1     =", R_P1/R_PLANET(P1)
!WRITE(*,*) "R_P2     =", R_P2/R_PLANET(P2)


END SUBROUTINE RESONANT_ORBIT
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!




!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
SUBROUTINE DEPARTURE_MGADSM(P1, P2, J1, J2, C3, ALPHA, BETA, ETA, DSM, V_INF, ASTEROID,N1,N2)
IMPLICIT NONE
INTEGER, INTENT(IN) :: P1,P2,N1,N2
DOUBLE PRECISION, INTENT(IN) :: C3, ALPHA, BETA, ETA, J1, J2, ASTEROID(N1,N2)
DOUBLE PRECISION, INTENT(INOUT) :: DSM, V_INF(3)

DOUBLE PRECISION :: OE(6), SV(6), V_SC(3), DT1, DT2, T
DOUBLE PRECISION :: R1(3), R2(3), V1(3), V2(3)
DOUBLE PRECISION :: SV_MID(6), V(6)
INTEGER :: OT

OT=1
!ASTEROID=0.D0
!N1=1
!N2=1

T=J2-J1
DT1=T*ETA*86400.D0
DT2=T*(1.D0-ETA)*86400.D0
!write(*,*) p1, j1
!write(*,*) p2, j2

CALL PLAN_ELEM(J1, P1,ASTEROID,N1,N2,OE)
SV=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))
R1=SV(1:3)
V1=SV(4:6)

CALL PLAN_ELEM(J2, P2,ASTEROID,N1,N2,OE)
SV=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))
R2=SV(1:3)
V2=SV(4:6)

V_SC(1)=V1(1)+SQRT(C3)*COS(ALPHA)*COS(BETA)
V_SC(2)=V1(2)+SQRT(C3)*SIN(ALPHA)*COS(BETA)
V_SC(3)=V1(3)+SQRT(C3)*SIN(BETA)

SV_MID=ORBIT_PROPOGATION(R1,V_SC,DT1)
CALL LAMBERT_BATTIN(SV_MID(1:3),R2,DT2,OT, V)
!CALL LAMBERT_SUN(SV_MID(1:3),R2,DT2,OT, V)
!CALL LAMBERT_GOODING(SV_MID(1:3),R2,DT2,OT, V)
DSM=NORM(V(1:3)-SV_MID(4:6))
V_INF=V(4:6)-V2
END SUBROUTINE DEPARTURE_MGADSM
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!




!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
SUBROUTINE GA_MGADSM(P1,P2,J1,J2,RP,GAMMA,ETA,V_INF1, V_INF2, DSM, G, ASTEROID, N1, N2)
IMPLICIT NONE
INTEGER, INTENT(IN) :: P1, P2, N1, N2
DOUBLE PRECISION, INTENT(IN) :: J1, J2, RP, GAMMA, ETA, V_INF1(3), ASTEROID(N1,N2)
DOUBLE PRECISION, INTENT(INOUT) :: V_INF2(3), DSM, G
DOUBLE PRECISION :: DT1, DT2, T, V_INF_MAG, DELTA, OE(6)
DOUBLE PRECISION :: SV(6), R1(3), R2(3), V1(3), V2(3), v(6), e
DOUBLE PRECISION :: I(3), J(3), K(3), V_INF_OUT(3), V_SC(3), SV_MID(6)
INTEGER :: OT

OT=1
T=J2-J1
DT1=T*ETA*86400.D0
DT2=T*(1.D0-ETA)*86400.D0
CALL PLAN_ELEM(J1, P1,ASTEROID,N1,N2,OE)
SV=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))
R1=SV(1:3)
V1=SV(4:6)

CALL PLAN_ELEM(J2, P2,ASTEROID,N1,N2,OE)
SV=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))
R2=SV(1:3)
V2=SV(4:6)

V_INF_MAG=NORM(V_INF1)

e=1.D0+RP*V_INF_MAG**2/MU_PLANET(P1)
DELTA=2.D0*ASIN(1.D0/e)

!write(*,*) "e", e
I=V_INF1/V_INF_MAG
J=CROSS_PRODUCT(I,V1)/NORM(CROSS_PRODUCT(I,V1))
K=CROSS_PRODUCT(I,J)
V_INF_OUT=COS(DELTA)*I+COS(GAMMA)*SIN(DELTA)*J+SIN(GAMMA)*SIN(DELTA)*K
V_INF_OUT=V_INF_OUT*V_INF_MAG
V_SC=V_INF_OUT+V1

SV_MID=ORBIT_PROPOGATION(R1,V_SC,DT1)

!write(*,*) "r_mid", norm(sv_mid(1:3))
!write(*,*) "v_mid", norm(sv_mid(4:6))

!CALL LAMBERT_BATTIN(SV_MID(1:3),R2,DT2,OT, V)
CALL LAMBERT_GOODING(SV_MID(1:3),R2,DT2,OT, V)
!write(*,*) "v(1)", v(1)
if (v(1).gt. 9.9d11) then
    !write(*,*) "called gooding"
    CALL LAMBERT_sun(SV_MID(1:3),R2,DT2,OT, V)
    if (v(1).gt.9.9d11) then
        !write(*,*) "called battin"
        CALL lambert_battin(sv_mid(1:3), r2, dt2, ot, V)
        if (v(1).gt.9.9d11) then
            !write(*,*) "called universal"
            call LAMBERT_UNIVERSAL(SV_MID(1:3),R2,DT2,OT, V)
        end if
    end if
    
end if
DSM=NORM(V(1:3)-SV_MID(4:6))
!write(*,*) "dsm", dsm
V_INF2=V(4:6)-V2
!write(*,*) "v2_inf", norm(v_inf2)
IF (RP<(R_PLANET(P1)*1.05D0))THEN
    G=-2.D0*LOG10(RP/(1.05D0*R_PLANET(P1)))
ELSE
    G=0
END IF

END SUBROUTINE GA_MGADSM
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!




!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
SUBROUTINE ARRIVAL_MGADSM(P1, e, RP, V_INF, dV)
IMPLICIT NONE
INTEGER, INTENT(IN) :: P1
DOUBLE PRECISION, INTENT(IN) :: e, RP, V_INF(3)
DOUBLE PRECISION, INTENT(INOUT) :: dV
DOUBLE PRECISION :: a, VP1, VP2, V_INF_MAG
V_INF_MAG=NORM(V_INF)
a=RP/(1.D0-e)
VP1=SQRT(V_INF_MAG**2+2.D0*MU_PLANET(P1)/RP)
VP2=SQRT(2.D0*MU_PLANET(P1)/RP-MU_PLANET(P1)/a)
dV=ABS(VP1-VP2)
END SUBROUTINE ARRIVAL_MGADSM
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!









!test gooding


!SUBROUTINE LAMBERT_GOODING2(R1,R2,DT,OT,NREV, V)
!IMPLICIT NONE
!
!DOUBLE PRECISION, INTENT(IN) :: R1(3), R2(3), dt
!DOUBLE PRECISION, INTENT(INOUT) :: V(6)
!INTEGER, INTENT(IN) :: OT, NREV
!
!DOUBLE PRECISION :: tol, r1mag, r2mag, C12(3), nu, c12mag
!double precision :: R1UNIT(3), R2UNIT(3), C12UNIT(3), TAN1U(3), TAN2U(3), GM
!DOUBLE PRECISION :: vr11,vt11, Vr12, vt12, vr21, vt21, vr22, vt22
!INTEGER :: N
!TOL=1.D-8
!
!R1MAG=NORM(R1)
!R2MAG=NORM(R2)
!
!R1UNIT=R1/r1mag
!R2UNIT=R2/r2mag
!
!C12=CROSS_PRODUCT(R1, R2)
!C12MAG=NORM(C12)
!C12UNIT=C12/C12MAG
!
!NU=ACOS(DOT_PRODUCT(R1,R2)/(R1MAG*R2MAG))
!
!!Determine the true anomaly angle using the orbit type
!!1 is prograde, 2 is retrograde
!IF (OT==1) THEN
!    IF (C12(3) <= 0.D0) THEN
!        NU=2.D0*PI-NU;
!    END IF
!END IF
!
!IF (OT==2) THEN
!    IF(C12(3)>=0.D0) THEN
!        NU=2.D0*PI-NU;
!    END IF
!END IF
!
!!write(*,*) "calling vlamb"
!!N=1
!CALL vlamb(gm,r1MAG,r2MAG ,NU, DT,N ,vr11,vt11, Vr12, vt12, vr21, vt21, vr22, vt22)
!
!
!
!
!!Tangential unit vectors
!TAN1U=CROSS_PRODUCT(C12UNIT, R1UNIT)
!TAN2U=CROSS_PRODUCT(C12UNIT, R2UNIT)
!IF (NU >=PI) THEN
!   TAN1U=-TAN1U
!   TAN2U=-TAN2U
!END IF
!
!
!!write(*,*) "exited vlamb", n, v(1:3)
!
!V(1:3)=VR11*R1UNIT+VT11*TAN1U
!V(4:6)=VR12*R2UNIT+VT12*TAN2U
!
!
!END SUBROUTINE LAMBERT_GOODING2
!
!
!SUBROUTINE tlamb(m, q, qsqfm1, x, n, t, dt, d2t, d3t)
!
!! Code converted using TO_F90 by Alan Miller
!! Date: 2014-03-05  Time: 00:30:37
!
!INTEGER, INTENT(IN)                      :: m
!DOUBLE PRECISION, INTENT(IN)             :: q
!DOUBLE PRECISION, INTENT(IN)             :: qsqfm1
!DOUBLE PRECISION, INTENT(IN)             :: x
!INTEGER, INTENT(IN)                      :: n
!DOUBLE PRECISION, INTENT(INOUT)            :: t
!DOUBLE PRECISION, INTENT(INOUT)            :: dt
!DOUBLE PRECISION, INTENT(INOUT)            :: d2t
!DOUBLE PRECISION, INTENT(INOUT)            :: d3t
!!IMPLICIT DOUBLE PRECISION (a-h, o-z)
!LOGICAL :: lm1, l1, l2, l3
!!DATA pi, sw /3.14159265358979D0, 0.4D0/
!!DATA sw /3.14159265358979D0, 0.4D0/
!DOUBLE PRECISION :: SW
!DOUBLE PRECISION :: A, AA, B, BB, F, FG1, FG1SQ, G, QSQ, QX, QZ, QZ2
!DOUBLE PRECISION :: TERM, TOLD, TQ, TQSUM, TQTERM, TTERM, TTMOLD, TWOI1, U, U01, U1I
!DOUBLE PRECISION :: U2I, U3I, XSQ, Y, Z
!INTEGER :: I, P
!SW=0.4D0
!
!lm1 = n == -1
!l1=n >= 1
!l2=n >= 2
!l3=n >= 3
!qsq=q*q
!xsq=x*x
!u=(1.d0-x)*(1.d0+x)
!IF(.NOT.lm1) THEN
!!         NEEDED IF SERIES AND OTHERWISE USEFUL WHEN Z=0
!  dt = 0.d0
!  d2t = 0.d0
!  d3t = 0.d0
!END IF
!IF(lm1.OR. m > 0 .OR. x < 0.d0 .OR. DABS(u) > sw) THEN
!!         DIRECT COMPUTATION (NOT SERIES)
!  y = DSQRT(DABS(u))
!  z = DSQRT(qsqfm1+qsq*xsq)
!  qx = q*x
!  IF (qx <= 0.d0) THEN
!    a = z-qx
!    b=q*z-x
!  END IF
!  IF (qx < 0.d0 .AND. lm1) THEN
!    aa=qsqfm1/a
!    bb=qsqfm1*(qsq*u-xsq)/b
!  END IF
!  IF(qx == 0.d0 .AND.lm1 .OR.qx > 0.d0) THEN
!    aa=z + qx
!    bb=q*z+x
!  END IF
!  IF (qx > 0.d0) THEN
!    a=qsqfm1/a
!    b=qsqfm1*(qsq*u-xsq)/bb
!  END IF
!  IF(.NOT.lm1) THEN
!    IF(qx*u >= 0.d0) THEN
!      g=x*z+q*u
!    ELSE
!      g= (xsq-qsq*u)/(x*z-q*u)
!    END IF
!    f= a*y
!    IF (x <= 1.d0) THEN
!      t=m*pi*DATAN2(f,g)
!    ELSE
!      IF (f > sW) THEN
!        t=DLOG(f+g)
!      ELSE
!        fg1 = f/(g+1.d0)
!        term = 2.d0*fg1
!        fg1sq = fg1*fg1
!        t = term
!        twoi1 = 1.d0
!        1                  twoi1=twoi1 + 2.d0
!        term = term*fg1sq
!        told = t
!        t = t + term/twoi1
!        IF (t /= told) GO TO 1
!!                         CONTINUE LOOPING FOR INVERSE TANH
!      END IF
!    END IF
!    t = 2.d0*(t/y + b)/u
!    IF (l1 .AND. z /= 0.d0) THEN
!      qz=q/z
!      qz2=qz*qz
!      qz=qz*qz2
!      dt= (3.d0*x*t -4.d0*(a + qx*qsqfm1)/z)/u
!      IF (l2) d2t=(3.d0*t+5.d0*x*dt+4.d0*qz*qz2*x*qsqfm1)/u
!      IF (l3) d3t=(8.d0*dt+7.d0*x*d2t-12.d0*qz*qz2*x*qsqfm1)/u
!    END IF
!  ELSE
!    dt=b
!    d2t=bb
!    d3t=aa
!  END IF
!ELSE
!!         COMPUTE BY SERIES
!  u01=1.d0
!  IF (l1) u1i = 1.d0
!  IF (l2) u2i = 1.d0
!  IF (l3) u3i = 1.d0
!  term = 4.d0
!  
!  tq = q*qsqfm1
!  i = 0
!  IF (q < 5.d-1) tqsum = 1.d0 - q*qsq
!  IF (q >= 5.d-1) tqsum = (1.d0/(1.d0+q)+q)*qsqfm1
!  ttmold = term/3.d0
!  t = ttmold*tqsum
!!             START OF THE LOOP
!  2       i = i + 1
!  p = i
!  u01 = u01*u
!  IF (l1 .AND. i > 1) u1i = u1i*u
!  IF (l2 .AND. i > 2) u2i = u2i*u
!  IF (l3 .AND. i > 3) u3i = u3i*u
!  term = term*(p-0.5D0)/p
!  tq=tq*qsq
!  tqsum = tqsum + tq
!  told = t
!  tterm = TERM/(2.D0+P+3.D0)
!  TQTERM=tterm*Tqsum
!  t=t-u01*((1.5D0*p+0.25D0)*tqterm/(p*p-0.25D0)-ttmold*tq)
!  ttmold = tterm
!  tqterm = tqterm*p
!  IF (l1) dt=dt+tqterm*u1i
!  IF (l2) d2t=d2t+tqterm*u2i*(p-1.d0)
!  IF (l3) d3t=d3t+tqterm*u3i*(p-1.d0)*(p-2.d0)
!  IF (i < n .OR. t /= told) GO TO 2
!!             END OF LOOP
!  IF (l3) d3t=8.d0*x*(1.5D0*d2t-xsq*d3t)
!  IF (l2) d2t=2.d0*(2.d0*xsq*d2t-dt)
!  IF (l1) dt=-2.d0*x*dt
!  t=t/xsq
!END IF
!RETURN
!END SUBROUTINE tlamb
!
!
!
!SUBROUTINE xlamb(m,q,qsqfm1,tin,n,x,xpl)
!
!INTEGER, INTENT(IN)                      :: m
!DOUBLE PRECISION, INTENT(IN OUT)         :: q
!DOUBLE PRECISION, INTENT(IN OUT)         :: qsqfm1
!DOUBLE PRECISION, INTENT(IN)             :: tin
!INTEGER, INTENT(IN OUT)                     :: n
!DOUBLE PRECISION, INTENT(OUT)            :: x
!DOUBLE PRECISION, INTENT(OUT)            :: xpl
!!IMPLICIT DOUBLE PRECISION (a-h,o-z)
!!DOUBLE PRECISION, PARAMETER :: pi=3.141592653589793D0
!!DOUBLE PRECISION, PARAMETER :: tol=3.d-7
!!DOUBLE PRECISION, PARAMETER :: c0=1.7D0
!!DOUBLE PRECISION, PARAMETER :: c1=0.5D0
!!DOUBLE PRECISION, PARAMETER :: c2=0.03D0
!!DOUBLE PRECISION, PARAMETER :: c3=0.15D0
!!DOUBLE PRECISION, PARAMETER :: c41=1.d0
!!DOUBLE PRECISION, PARAMETER :: c42=0.24D0
!
!DOUBLE PRECISION :: TOL, C0 , C1, C2, C3, C41, C42
!DOUBLE PRECISION :: THR2, T0, DT, D2T, D3T, TDIFF, W, XM, XMOLD, XTEST
!DOUBLE PRECISION :: TDIFFM, TMIN, D2T2, T, X0, TDIFF0
!INTEGER :: i
!
!!d8rt(x)=DSQRT(DSQRT(DSQRT(x)))
!
!
!tol=3.d-7
!c0=1.7D0
!c1=0.5D0
!c2=0.03D0
!c3=0.15D0
!c41=1.d0
!c42=0.24D0
!
!
!
!thr2 = DATAN2(qsqfm1, 2.d0*q)/pi
!IF (m == 0) THEN
!!         SINGLE-REV STARTER FROM T( AT X=0) &BILINEAR
!  n = 1
!  CALL tlamb(m,q,qsqfm1, 0.d0, 0, t0, dt, d2t, d3t)
!  tdiff=tin-t0
!  IF (tdiff <= 0.d0) THEN
!    x=t0*tdiff/(-4.d0*tin)
!!           -4 IS THE VALUE OF DT, FOR X=0
!  ELSE
!    x=-tdiff/(tdiff+4.d0)
!    w=x+c0*DSQRT(2.d0*(1.d0-thr2))
!    IF (w < 0.d0) x=x-DSQRT(DSQRT(DSQRT(DSQRT(-w))))* (x+DSQRT(tdiff/(tdiff+1.5D0)))
!    w=4.d0/(4.d0+tdiff)
!    x=x*(1.d0+x*(c1*w-c2*x*DSQRT(w)))
!  END IF
!ELSE
!!    WITH MILTI-REVS, FIRST GET T(MIN) AS BASIS FOR STARTER
!  xm=1.d0/(1.5D0*(m+5.d-1)*pi)
!  IF (thr2 < 5.d-1) xm=DSQRT(DSQRT(DSQRT(2.D0*THR2)))*xm
!  IF(thr2 > 5.d-1) xm=(2.d0-DSQRT(DSQRT(DSQRT(2.d0-2.d0*thr2))))*xm
!!         STARTER FOR TMIN
!  DO  i=1,12
!    CALL tlamb(m,q,qsqfm1,xm,3,tmin,dt,d2t,d3t)
!    IF (d2t == 0.d0) GO TO 2
!    xmold=xm
!    xm=xm-dt*d2t/(d2t*d2t-dt*d3t/2.d0)
!    xtest=DABS(xmold/xm-1.d0)
!    IF (xtest <= tol) GO TO 2
!  END DO
!  n=-1
!  RETURN
!!         BREAK OFF & EXIT IF TMIN NOT LOCATED - SHOULD NEVER HAPPER
!!         NOW PROCEED FROM T(MIN) TO FULL STARTER
!  2    CONTINUE
!  tdiffm=tin-tmin
!  IF (tdiffm < 0.d0) THEN
!    n=0
!    RETURN
!!         EXIT IF NO SOLUTION WITH THIS M
!  ELSE IF (tdiffm == 0.d0) THEN
!    x=xm
!    n=1
!    RETURN
!!           EXIT IF UNIQUE SOLUTION ALREADY FROM T(MIN)
!  ELSE
!    n=3
!    IF (d2t == 0.d0) d2t=6.d0*m*pi
!    x=DSQRT(tdiffm/(d2t/2.d0+tdiffm/(1.d0-xm)**2))
!    w=xm+x
!    w=w*4.d0/(4.d0+tdiffm)+(1.d0-w)**2
!    x=x*(1.d0-(1.d0+m+c41*(thr2-0.5D0))/(1.d0+c3*m)*  &
!        x*(c1*w+c2*x*DSQRT(w)))+xm
!    d2t2=d2t/2.d0
!    IF (x >= 1.d0) THEN
!      n=1
!      GO TO 3
!    END IF
!!           NO INFINITE SOLUTIONS WITH X>XM
!  END IF
!!         WE HAVE THE STARTERS, SO PROCEED BY HALLEY
!END IF
!5  CONTINUE
!DO  i=1,3
!  CALL tlamb(m,q,qsqfm1,x,2,t,dt,d2t,d3t)
!  t=tin-t
!  IF (dt /= 0.d0) x=x+t*dt/(dt*dt+t*d2t/2.d0)
!  END DO
!  IF (n /= 3) RETURN
!!       EXIT IF ONLY ONE SOLUTION, NORMALLY WHEN M=0
!  n=2
!  xpl=x
!!       SECOND MULTI-REV STARTER
!  3  CALL tlamb(m,q,qsqfm1,0.d0,0,t0,dt,d2t,d3t)
!  tdiff0=t0-tmin
!  tdiff=tin-t0
!  IF (tdiff <= 0) THEN
!    x=xm-DSQRT(tdiffm/(d2t2-tdiffm*(d2t2/tdiff0-1.d0/xm**2)))
!  ELSE
!    x=-tdiff/(tdiff+4.d0)
!    w=x0+c0*DSQRT(2.d0*(1.d0-thr2))
!    IF (w < 0.d0) x=x-DSQRT(DSQRT(DSQRT(DSQRT(-W))))* (x+DSQRT(tdiff/(tdiff+1.5D0*t0)))
!    w=4.d0/(4.d0+tdiff)
!    x=x*(1.d0+(1.d0+m+c42*(thr2-0.5D0))/(1.d0+c3*m)*x *(c1*w-c2*x*DSQRT(w)))
!    IF (x <= -1.d0) THEN
!      n=n-1
!!           NO INFITE SOLUTIONS WITH X<XM
!      IF (n == 1) x=xpl
!    END IF
!  END IF
!  GO TO 5
!END SUBROUTINE xlamb
!
!SUBROUTINE vlamb(gm,r1,r2,th,tdelt,n,vr11,vt11,  &
!    vr12, vt12, vr21, vt21, vr22, vt22)
!
!DOUBLE PRECISION, INTENT(IN OUT)         :: gm
!DOUBLE PRECISION, INTENT(IN)             :: r1
!DOUBLE PRECISION, INTENT(IN)             :: r2
!DOUBLE PRECISION, INTENT(IN)             :: th
!DOUBLE PRECISION, INTENT(IN)             :: tdelt
!INTEGER, INTENT(INOUT)                      :: n
!DOUBLE PRECISION, INTENT(IN OUT)            :: vr11
!DOUBLE PRECISION, INTENT(IN OUT)            :: vt11
!DOUBLE PRECISION, INTENT(INOUT)            :: vr12
!DOUBLE PRECISION, INTENT(INOUT)            :: vt12
!DOUBLE PRECISION, INTENT(INOUT)            :: vr21
!DOUBLE PRECISION, INTENT(INOUT)            :: vt21
!DOUBLE PRECISION, INTENT(INOUT)            :: vr22
!DOUBLE PRECISION, INTENT(INOUT)            :: vt22
!!IMPLICIT DOUBLE PRECISION (a-h,o-z)
!!DOUBLE PRECISION, PARAMETER :: pi=2.141592653589793D0
!!DOUBLE PRECISION, PARAMETER :: twopi=2.d0*pi
!DOUBLE PRECISION :: TWOPI, THR2, DR,R1R2,R1R2TH, CSQ, C, S, GMS, QSQFM1, Q
!DOUBLE PRECISION :: RHO, SIG, T, X, X1, X2, UNUSED, QZMINX, QZPLX, ZPLQX, VT2, VR1, VT1, VR2
!
!INTEGER :: I, M
!
!!write(*,*) "entering vlamb"
!
!TWOPI=2.D0*PI
!m=FLOOR(th/twopi)
!thr2=th/2.d0-m*pi
!dr=r1-r2
!r1r2=r1*r2
!r1r2th=4.d0*r1r2*dsIn(thr2)**2
!csq=dr**2+r1r2th
!c=DSQRT(csq)
!s=(r1+r2+c)/2.d0
!gms=DSQRT(gm*s/2.d0)
!qsqfm1=c/s
!q=DSQRT(r1r2)*DCOS(thr2)/s
!IF (c /= 0.d0) THEN
!  rho=dr/c
!  sig=r1r2th/csq
!ELSE
!  rho=0.d0
!  sig=1.d0
!END IF
!t=4.d0*gms*tdelt/s**2
!CALL xlamb(m,q,qsqfm1, T, n, x1, x2)
!!         PROCEED FOR SINGLE SOLUTION, OR A PAIR
!DO  i=1,n
!  IF (i == 1) THEN
!    x=x1
!  ELSE
!    x=x2
!  END IF
!  CALL tlamb(m,q,qsqfm1,x,-1,unused,qzminx,qzplx,zplqx)
!  vt2=gms*zplqx*DSQRT(sig)
!  vr1=gms*(qzminx-qzplx*rho)/r1
!  vt1=vt2/r1
!  vr2=-gms*(qzminx-qzplx*rho)/r2
!  vt2=vt2/r2
!  IF (i == 1) THEN
!    vr11=vr1
!    vt11=vt1
!    vr12=vr2
!    vt12=vt2
!  ELSE
!    vr21=vr1
!    vt21=vt1
!    vr22=vr2
!    vt22=vt2
!  END IF
!END DO
!RETURN
!END SUBROUTINE vlamb


END MODULE ORBITAL_FUNCTIONS
