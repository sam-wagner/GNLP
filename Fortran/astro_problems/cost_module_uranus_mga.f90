MODULE COST_MODULE
USE ORBITAL_FUNCTIONS
USE MATH_MODULE
USE CONSTANTS
IMPLICIT NONE
 ! THIS MODULE MINIMIZES AN EARTH-URANUS TRANSFER
CONTAINS


! THIS OBJECTIVE FUNCTION IS THE SAME AS THE CASSINI AND GALILEO MGA COST FUNCTIONS.
! SO FAR I HAVEN'T USED THIS TO TEST ANYMORE THAN TEST FOR POSSIBLE URANUS TRAJECTORIES.
! FROM WHAT I'VE FOUND IT IS SOMEWHAT DIFFICULT, MOSTLY DUE TO LARGE ARRIVAL INSERTION BURNS
! OR HIGH LAUNCH C3 VALUES WHICH ARE OFTEN REQUIRED.  IT MAY BE LIKELY THE ORBIT INSERTION 
! TERMS NEED TO BE ADJUSTED AS WELL.

SUBROUTINE COST(N_DOUBLE, N_INT, N1, N2, X, CHROM_INT, FITNESS, ARRAY, G_CON, NCON)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N_DOUBLE, N_INT, N1, N2, CHROM_INT(N_INT), NCON
DOUBLE PRECISION, INTENT(IN) :: X(N_DOUBLE), ARRAY(N1,N2)
DOUBLE PRECISION, INTENT(INOUT) :: FITNESS, G_CON(NCON)

DOUBLE PRECISION :: OE(6), R1(3), V1(3), R2(3), V2(3), V(6), V1_INF(3), V1_INF_MAG,&
	C3_LIM, DV1, DV2, DVF, DT, SV(6), VF_INF(3), &
	ARR_PEN, VP1, VP2, RPL, eF, rpf, GA_CHECK, asteroid(N1,N2), DV_GA(CHROM_INT(2)), &
	G(CHROM_INT(2)), R_P(CHROM_INT(2)), V_INF_IN(CHROM_INT(2),3), &
	V_INF_OUT(CHROM_INT(2),3), J(CHROM_INT(2)+2), LENGTH_PENALTY, TOF_TOT

INTEGER :: PL, PF, NGA, I, K, OT, NREV, n11, n22, P(CHROM_INT(2)+1)

!SET THE ORBIT TYPE TO PROGRADE FOR THE LAMBERT SOLVER
OT=1
!NREV=0 !NO LONGER REVELANT AS THE CURRENT LAMBERT SOLVERS ARE ALL 0-REVOLUTION
NGA=CHROM_INT(2)
PL=CHROM_INT(1)

! PROPOSED EARTH DEPARTURE LIMITS.  LARGER C3'S ARE ALLOWED BUT RESULT IN A REQUIRE DELTA-V 
! FOR THE SPACECRAFT.  THE DEPARTURE ORBIT PERIGEE RADIUS IS ASSUMED TO BE 300 KM. 
C3_LIM=20.d0
RPL=6378.137D0+300.D0
! PROPOSED URANUS INSERTION ORBIT
RPF=46200.D0
eF=.98D0

!ASSIGN THE VALUES FOR THE PLANETARY ARRAY ORDER FROM THE INPUT VARIABLES
P(NGA+1)=CHROM_INT(3)
P(1:NGA)=CHROM_INT(4:NGA+3)



!SET ALL THE RELEVANT DATES FOR THE MIDDION FROM THE REAL VALUES CHROMOSOME
J(1)=X(1)
DO I=2,NGA+1,1
    J(I)=X(I)+J(I-1)
END DO
J(NGA+2)=X(N_DOUBLE)+J(NGA+1)


! DEPARTURE PHASE
CALL PLAN_ELEM(J(1), PL,ARRAY,N1,N2,OE)
SV=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))
R1=SV(1:3) 
V1=SV(4:6)
CALL PLAN_ELEM(J(2), P(1),ARRAY,N1,N2,OE)
SV=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))
R2=SV(1:3)
V2=SV(4:6)


DT=(J(2)-J(1))*86400.D0
!CALL lambert_GOODING(R1,R2,DT,OT,V)
CALL lambert_sun(R1,R2,DT,OT,V)
! JUST IN CASE THE SUN LAMBERT SOLVER FAILS.
IF(ISNAN(NORM(V)) .or. v(1).gt.9.9d11)THEN
    call LAMBERT_battin(R1,R2,DT,OT, V)
END IF
V1_INF=V1-V(1:3)
V1_INF_MAG=NORM(V1_INF)


R1=R2
V1=V2
V_INF_IN(1,1:3)=V(4:6)-V1

!GRAVITY ASSISTS
DO I=1,NGA,1
    DT=(J(I+2)-J(I+1))*86400.D0
    CALL PLAN_ELEM(J(I+2), P(I+1),ARRAY,N1,N2,OE)
    SV=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))  
    R2=SV(1:3)
    V2=SV(4:6)
    CALL lambert_sun(R1,R2,DT,OT,V)
    IF(ISNAN(NORM(V)) .or. V(1).gt.9.9d11)THEN
        CALL LAMBERT_BATTIN(R1,R2,DT,OT, V)
    END IF
    V_INF_OUT(I,1:3)=V(1:3)-V1
    IF (I.EQ.NGA) THEN
        VF_INF=V(4:6)-V2
    ELSE
        V_INF_IN(I+1,1:3)=V(4:6)-V2
    END IF
    V1=V2
    R1=R2
END DO


! TEST FOR POSSIBLE RESONANT ORBITS AROUND THE EARTH OR VENUS (SO FAR THIS HAS REALLY ONLY BEEN 
! USED FOR THE GALILEO MISSION BUT YOU COULD EASILY TEST FOR RESONANT ORBITS AROUND OTHER PLANETS 
! BY ADJUSTING THE PERIOD THAT IS CHECK.  RIGHT NOW IT SIMPLY CHECKS IF THE TRAJECTORY IS CLOSE TO
! A 2 YEAR EARTH RESONANT ORBIT.

!DO I=1,NGA,1
!    !RESONANT_ORBIT(T1,T2,P1,P2, V1_INF_IN,V2_INF_OUT, V1_INF_OUT, V2_INF_IN)
!    IF (I.GE.1 .and. i.lt.nga) THEN
!        IF(P(I).EQ.3 .AND. P(I+1).EQ.3) THEN
!            IF(ABS(J(I+2)-J(I+1)-2.D0*365.242189D0).LE.5.D0) THEN
!                CALL RESONANT_ORBIT(J(I+1),J(I+2),P(I),P(I+1), V_INF_IN(I,1:3),V_INF_OUT(I+1,1:3), &
!                                   V_INF_OUT(I,1:3), V_INF_IN(I+1,1:3))
!            END IF
!        ELSE IF (P(I).EQ.2 .AND. P(I+1).EQ.2) THEN
!            IF(ABS(J(I+2)-J(I+1)-2.D0*224.7D0).LE.5.D0) THEN
!                CALL RESONANT_ORBIT(J(I+1),J(I+2),P(I),P(I+1), V_INF_IN(I,1:3),V_INF_OUT(I+1,1:3), &
!                                   V_INF_OUT(I,1:3), V_INF_IN(I+1,1:3))
!            END IF
!        END IF
!    END IF 
!END DO

DO I=1,NGA,1
    CALL GA_MGA(V_INF_IN(I,1:3), V_INF_OUT(I,1:3), P(I), R_P(I), DV_GA(I), G(I))
    !WRITE(*,*) "DV_GA", I, DV_GA(I)
END DO


!ARRIVAL PHASE
CALL ARRIVAL_MGADSM(P(NGA+1), eF, RPF, VF_INF, dVF)


DV1=SQRT(C3_LIM+2.D0*MU_PLANET(PL)/RPL)
DV2=SQRT(V1_INF_MAG**2+2.D0*MU_PLANET(PL)/RPL)
IF (V1_INF_MAG.GT.SQRT(C3_LIM)) THEN
    ARR_PEN=DV2-DV1
    
!ELSE IF (V1_INF_MAG.LT.SQRT(2.D0)) THEN
!    ARR_PEN=SQRT(2.D0)-V1_INF_MAG
ELSE
    ARR_PEN=0.D0
END IF

!Penalize flights longer than 10 years
TOF_TOT=J(NGA+2)-J(1)
LENGTH_PENALTY=0.D0
!IF (TOF_TOT.GT.(365.25D0*15.D0)) THEN
!    LENGTH_PENALTY=0.1D0*(TOF_TOT-365.25D0*15.D0)
!END IF

FITNESS=ARR_PEN+SUM(DV_GA)+SUM(G)+DVF+LENGTH_PENALTY

!write(*,*) "c3lim", c3_lim
!write(*,*) "launch c3", v1_inf_mag**2
!write(*,*) "v1_inf", v1_inf_mag
!WRITE(*,*) "DV_DEP", ARR_PEN
!DO I=1,NGA,1
!    write(*,*) "dv_ga", I, dv_ga(I)
!    write(*,*) "r_P", I, r_p(I)
!    write(*,*) "alt", r_p(I)-r_planet(p(I))
!    WRITE(*,*) "G", I, G(I)
!END DO 
!write(*,*) "dvF", dvF
!WRITE(*,*) "FINAL FITNESS VALUE", FITNESS



!P1=CHROM_INT(1)
!P2=CHROM_INT(2)
!P3=CHROM_INT(3)
!P4=CHROM_INT(4)
!P5=CHROM_INT(5)
!
!T1=X(1)
!T2=X(2)+T1
!T3=X(3)+T2
!T4=X(4)+T3
!T5=X(5)+T4
!
!DT12=(T2-T1)*86400.D0
!DT23=(T3-T2)*86400.D0
!DT34=(T4-T3)*86400.D0
!DT45=(T5-T4)*86400.D0
!
!!OT=1
!!ASTEROID=ARRAY
!
!CALL PLAN_ELEM(T1, P1,ASTEROID,N1,N2,OE)
!SV1=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))
!CALL PLAN_ELEM(T2, P2,ASTEROID,N1,N2,OE)
!SV2=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))
!CALL PLAN_ELEM(T3, P3,ASTEROID,N1,N2,OE)
!SV3=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))
!CALL PLAN_ELEM(T4, P4,ASTEROID,N1,N2,OE)
!SV4=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))
!CALL PLAN_ELEM(T5, P5,ASTEROID,N1,N2,OE)
!SV5=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))
!
!R1=SV1(1:3)
!V1=SV1(4:6)
!R2=SV2(1:3)
!V2=SV2(4:6)
!R3=SV3(1:3)
!V3=SV3(4:6)
!R4=SV4(1:3)
!V4=SV4(4:6)
!R5=SV5(1:3)
!V5=SV5(4:6)
!
!CALL lambert_GOODING(R1,R2,DT12,OT, 0, V)
!IF(ISNAN(NORM(V)))THEN
!    V=LAMBERT_BATTIN(R1,R2,DT12,OT)
!
!END IF
!V1_INF=V(1:3)-V1
!V21_INF=V(4:6)-V2
!
!CALL lambert_GOODING(R2,R3,DT23,OT, 0, V)
!IF(ISNAN(NORM(V)))THEN
!    V=LAMBERT_BATTIN(R2,R3,DT23,OT)
!    !CALL LAMBERT_SUN(R2,R3,DT23,OT,V)
!!    CALL lambert_universal(R2,R3,DT23,OT, V)
!END IF
!V22_INF=V(1:3)-V2
!V31_INF=V(4:6)-V3
!
!CALL lambert_GOODING(R3,R4,DT34,OT,0,V)
!IF(ISNAN(NORM(V)))THEN
!    V=LAMBERT_BATTIN(R3,R4,DT34,OT)
!!    CALL lambert_universal(R2,R3,DT23,OT, V)
!END IF
!V32_INF=V(1:3)-V3
!V41_INF=V(4:6)-V4
!
!CALL lambert_GOODING(R4,R5,DT45,OT, 0, v)
!IF(ISNAN(NORM(V)))THEN
!    V=LAMBERT_BATTIN(R4,R5,DT45,OT)
!!    CALL lambert_universal(R2,R3,DT23,OT, V)
!END IF
!V42_INF=V(1:3)-V4
!V5_INF=V(4:6)-V5
!
!
!
!
!!CHECK HERE FOR RESONANT ORBITS BETWEEN EARTH-EARTH TRANSFER
!if (p3.eq.3 .and. p4.eq.3) then
!    IF (ABS(T4-T3-2.D0*365.242189D0).LE.5.D0)THEN
!        !WRITE(*,*) "T3     ", T3
!        !WRITE(*,*) "T4     ", T4
!        !WRITE(*,*) "P3     ", P3
!        !WRITE(*,*) "P4     ", P4
!        !WRITE(*,*) "V31_INF", V31_INF
!        !WRITE(*,*) "v42_INF", V42_INF
!        CALL RESONANT_ORBIT(T3,T4,P3,P4, V31_INF,V42_INF, V32_INF, V41_INF)
!    END IF
!end if
!
!if (p3.eq.2 .and. p4.eq.2) then
!    IF (ABS(T4-T3-2.D0*224.7D0).LE.5.D0)THEN
!        CALL RESONANT_ORBIT(T3,T4,P3,P4, V31_INF,V42_INF, V32_INF, V41_INF)
!    END IF
!end if
!
!!WRITE(*,*) "V_INF_IN ", V21_INF
!!WRITE(*,*) "V_INF_OUT", V22_INF
!!WRITE(*,*) "V_INF_IN ", V31_INF
!!WRITE(*,*) "V_INF_OUT", V32_INF
!!WRITE(*,*) "V_INF_IN ", V41_INF
!!WRITE(*,*) "V_INF_OUT", V42_INF
!
!
!V1_INF_MAG=NORM(V1_INF)
!V21_INF_MAG=NORM(V21_INF)
!V22_INF_MAG=NORM(V22_INF)
!V31_INF_MAG=NORM(V31_INF)
!V32_INF_MAG=NORM(V32_INF)
!V41_INF_MAG=NORM(V41_INF)
!V42_INF_MAG=NORM(V42_INF)
!V5_INF_MAG=NORM(V5_INF)
!
!!WRITE(*,*) "V_INF_IN ", V21_INF
!!WRITE(*,*) "V_INF_OUT", V22_INF
!!WRITE(*,*) "V_INF_IN ", V31_INF
!!WRITE(*,*) "V_INF_OUT", V32_INF
!!WRITE(*,*) "V_INF_IN ", V41_INF
!!WRITE(*,*) "V_INF_OUT", V42_INF
!
!!FIRST GRAVITY ASSIST
!CALL GA_MGA(V21_INF, V22_INF, P2, R_P1, DV_GA1, G1)
!
!
!!SECOND GRAVITY ASSIST
!CALL GA_MGA(V31_INF, V32_INF, P3, R_P2, DV_GA2, G2)
!
!!THIRD GRAVITY ASSIST
!CALL GA_MGA(V41_INF, V42_INF, P4, R_P3, DV_GA3, G3)
!
!
!
! 
!    
!C3_LIM=17.D0
!
!!IF (V1_INF_MAG.LT.SQRT(2.D0)) THEN
!!    DV1=SQRT(2.D0)-V1_INF_MAG
!!ELSE
!!    DV1=0.D0
!!END IF
!
!IF (V1_INF_MAG.GT.SQRT(C3_LIM)) THEN
!    ARR_PEN=V1_INF_MAG-SQRT(C3_LIM)
!ELSE IF (V1_INF_MAG.LT.SQRT(2.D0)) THEN
!    ARR_PEN=SQRT(2.D0)-V1_INF_MAG
!ELSE
!    ARR_PEN=0.D0
!END IF
!
!CALL ARRIVAL_MGADSM(P5, eF, RPF, V5_INF, dV2)
!
!!FITNESS=dV1+dV2+dV_GA1+dV_GA2+dV_GA3+G1+G2+G3
!FITNESS=ARR_PEN+dV2+dV_GA1+dV_GA2+dV_GA3+G1+G2+G3
!!write(*,*) "c3lim", c3_lim
!!write(*,*) "launch c3", v1_inf_mag**2
!!write(*,*) "dv1", dv1
!!write(*,*) "dv2", dv2
!!write(*,*) "dv_ga1", dv_ga1
!!write(*,*) "dv_ga2", dv_ga2
!!write(*,*) "dv_ga3", dv_ga3
!!write(*,*) "r_P1", r_p1
!!write(*,*) "r_P2", r_p2
!!write(*,*) "r_P3", r_p3
!!
!!write(*,*) "alt", r_p1-r_planet(p2)
!!write(*,*) "alt", r_p2-r_planet(p3)
!!write(*,*) "alt", r_p3-r_planet(p4)
!!
!!write(*,*) "g1", g1
!!write(*,*) "g2", g2
!!write(*,*) "g2", g3
!
!!FITNESS=V1_INF_MAG+V5_INF_MAG+dV_GA1+dV_GA2+dV_GA3+G1+G2+G3
!WRITE(*,*) FITNESS
! CHECK TO SEE THAT THIS FLYBY ORDER HASN'T ALREADY BEEN OPTIMIZED
do i=1,n1,1
    !WRITE(*,*) I
    NGA=NINT(ARRAY(I,2))
    IF (NGA.GT.0) THEN
        GA_CHECK=0.D0
        DO K=1,NGA,1
            GA_CHECK=GA_CHECK+ABS(ARRAY(I,3+K)-DBLE(P(K)))
            !WRITE(*,*) ARRAY(I,3+K), P(K)
            !WRITE(*,*) GA_CHECK
        END DO
    
        IF (GA_CHECK.LT.1.D-8 .AND. FITNESS.GT.ARRAY(I,N2)) THEN
            !WRITE(*,*) GA_CHECK, FITNESS
            FITNESS=1.D30
        
        END IF
    END IF
    !if(NINT(ARRAY(I,2)).EQ.CHROM_INT(2) .AND. NINT(ARRAY(I,3)).EQ.CHROM_INT(3) .AND. &
    !   NINT(ARRAY(I,4)).EQ.CHROM_INT(4) .and. fitness.gt.array(i,6)) then
    !    fitness=1.d30
    !end if
end do
nga=chrom_int(2)
!MAKE SURE SOLUTIONS AREN'T PICKED THAT GO TO AN OUTER PLANET AND BACK TO AN INNER PLANET
DO I = 1,NGA-1,1
    !WRITE(*,*) I, P(I)
    !WRITE(*,*) P(I)
    IF (P(I).EQ.6) THEN
        FITNESS=1.D30
        !WRITE(*,*) "VIOLATED"
    END IF
END DO

IF (NGA.GT.2) THEN
    DO I = 1,NGA-2,1
        IF (P(I).GT.4 .AND. P(I+1).LT.P(I)) THEN
            FITNESS=1.D30
            !WRITE(*,*) "VIOLATED"
        END IF
    END DO
END IF

!DEALLOCATE (DV_GA, G, R_P,V_INF_IN, V_INF_OUT, P, J)
if (fitness.gt.9.9d7 .or. isnan(fitness) ) fitness=1.d30
END SUBROUTINE COST


END MODULE COST_MODULE
