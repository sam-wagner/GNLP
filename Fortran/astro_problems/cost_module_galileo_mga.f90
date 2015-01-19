MODULE COST_MODULE
USE ORBITAL_FUNCTIONS
USE MATH_MODULE
USE CONSTANTS
IMPLICIT NONE
 ! THIS MODULE MINIMIZES AN EARTH-SATURN TYPE TRANSFER
CONTAINS



SUBROUTINE COST(N_DOUBLE, N_INT, N1, N2, X, CHROM_INT, FITNESS, ARRAY, G_CON, NCON)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N_DOUBLE, N_INT, N1, N2, CHROM_INT(N_INT), NCON
DOUBLE PRECISION, INTENT(IN) :: X(N_DOUBLE), ARRAY(N1,N2)
DOUBLE PRECISION, INTENT(INOUT) :: FITNESS, G_CON(NCON)
DOUBLE PRECISION :: OE(6), R1(3), V1(3), R2(3), V2(3), R3(3), V3(3)
DOUBLE PRECISION :: R4(3), V4(3), R5(3), V5(3)
DOUBLE PRECISION :: SV1(6), SV2(6), SV3(6), SV4(6), SV5(6) 
DOUBLE PRECISION :: V1_INF(3),V21_INF(3)
DOUBLE PRECISION :: V22_INF(3), V31_INF(3), V32_INF(3), V41_INF(3), V42_INF(3), V5_INF(3)
DOUBLE PRECISION :: V1_INF_MAG, V3_INF_MAG, V21_INF_MAG, V22_INF_MAG, V31_INF_MAG
DOUBLE PRECISION :: V32_INF_MAG, V41_INF_MAG, V42_INF_MAG, V5_INF_MAG
DOUBLE PRECISION :: V(6), R_P1, R_P2, R_P3, T1, T2, T3,T4,T5, DT12, DT23, DT34,DT45
DOUBLE PRECISION :: dV_GA1, dV_GA2, dV_GA3, G1, G2, G3
DOUBLE PRECISION :: dV1, dV2, Vc, VPF, af, VP1_C3_LIM, C3_LIM, dV1_C3_LIM 
INTEGER :: P1,P2,P3,P4,P5, COUNT

DOUBLE PRECISION, ALLOCATABLE :: DV_GA(:), G(:), R_P(:),V_INF_IN(:,:), V_INF_OUT(:,:), J(:)
DOUBLE PRECISION :: DVF, DT, SV(6), VF_INF(3)!, V_INF_IN(3), V_INF_OUT(3)
DOUBLE PRECISION :: ARR_PEN, VP1, VP2, RPL, eF, rpf, GA_CHECK, asteroid(1,1)
INTEGER :: PL, PF, NGA, I, K, OT, NREV, n11, n22
INTEGER, ALLOCATABLE :: P(:)

OT=1
NREV=0
ASTEROID=0.D0
NGA=CHROM_INT(2)
!WRITE(*,*) "NGA", NGA
PL=CHROM_INT(1)
!FINAL ORBIT INSERTION PARAMETERS
eF=0.998
RPF=286000.D0


!WRITE(*,*) "MADE IT TO COST MODULE"
ALLOCATE (DV_GA(NGA), G(NGA), R_P(NGA),V_INF_IN(NGA,3), V_INF_OUT(NGA,3), P(NGA+1), J(NGA+2))

P(NGA+1)=CHROM_INT(3)
P(1:NGA)=CHROM_INT(4:NGA+3)


!WRITE(*,*) "NGA", NGA
!WRITE(*,*) "P", P


J(1)=X(1)
DO I=2,NGA+1,1
    J(I)=X(I)+J(I-1)
    !J(3)=X(3)+J(2)
    !J(4)=X(4)+J(3)
END DO
J(NGA+2)=X(N_DOUBLE)+J(NGA+1)

!WRITE(*,*) "NGA", NGA
!DO I=1,NGA+2,1
!    WRITE(*,*) "J", I, J(I)
!END DO

! DEPARTURE PHASE

CALL PLAN_ELEM(J(1), PL,ASTEROID,N11,N22,OE)
SV=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))
R1=SV(1:3) 
V1=SV(4:6)
CALL PLAN_ELEM(J(2), P(1),ASTEROID,N11,N22,OE)
SV=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))
R2=SV(1:3)
V2=SV(4:6)

DT=(J(2)-J(1))*86400.D0
!CALL lambert_GOODING(R1,R2,DT,OT,V)

CALL lambert_sun(R1,R2,DT,OT,V)
IF(ISNAN(NORM(V)) .or. v(1).gt.9.9d11)THEN
    CALL LAMBERT_battin(R1,R2,DT,OT, V)
END IF
V1_INF=V1-V(1:3)
V1_INF_MAG=NORM(V1_INF)


R1=R2
V1=V2
V_INF_IN(1,1:3)=V(4:6)-V1

!GRAVITY ASSISTS
DO I=1,NGA,1
    DT=(J(I+2)-J(I+1))*86400.D0
    CALL PLAN_ELEM(J(I+2), P(I+1),ASTEROID,N11,N22,OE)
    SV=OE2SV(OE(1),OE(2),OE(3),OE(4),OE(5),OE(6))  
    R2=SV(1:3)
    V2=SV(4:6)
    !CALL lambert_GOODING(R1,R2,DT,OT,V)
    
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


! TEST FOR POSSIBLE RESONANT ORBITS AROUND THE EARTH OR VENUS
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
!WRITE(*,*) "CALCULATING GA'S"
DO I=1,NGA,1
    CALL GA_MGA(V_INF_IN(I,1:3), V_INF_OUT(I,1:3), P(I), R_P(I), DV_GA(I), G(I))
    !WRITE(*,*) "DV_GA", I, DV_GA(I)
END DO


!ARRIVAL
!WRITE(*,*) "VF_INF", NORM(VF_INF)
CALL ARRIVAL_MGADSM(P(NGA+1), eF, RPF, VF_INF, dVF)


C3_LIM=17.D0
RPL=6378.137D0+300.D0
DV1=SQRT(C3_LIM+2.D0*MU_PLANET(PL)/RPL)
DV2=SQRT(V1_INF_MAG**2+2.D0*MU_PLANET(PL)/RPL)
IF (V1_INF_MAG.GT.SQRT(C3_LIM)) THEN
    ARR_PEN=DV2-DV1
    
!ELSE IF (V1_INF_MAG.LT.SQRT(2.D0)) THEN
!    ARR_PEN=SQRT(2.D0)-V1_INF_MAG
ELSE
    ARR_PEN=0.D0
END IF


FITNESS=ARR_PEN+SUM(DV_GA)+SUM(G)+DVF

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



! CHECK TO SEE THAT THIS FLYBY ORDER HASN'T ALREADY BEEN OPTIMIZED MORE THAN 4 TIMES
COUNT=0
do i=1,n1,1
    NGA=NINT(ARRAY(I,2))
    IF (NGA.GT.0) THEN
        GA_CHECK=0.D0
        DO K=1,NGA,1
            GA_CHECK=GA_CHECK+ABS(ARRAY(I,3+K)-DBLE(P(K)))
            !WRITE(*,*) ARRAY(I,3+K), P(K)
            !WRITE(*,*) GA_CHECK
        END DO
        IF (GA_CHECK.LT.1.D-8 .AND. COUNT.LT.3) THEN
            COUNT=COUNT+1
            GA_CHECK=1.D8
        END IF
    
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

! CHECK TO MAKE SURE THE FITNESS VALUES ISN'T A NAN SO THAT THE OBJECT FUNCTION 
! WON'T CAUSE PROBLEMS WITH THE SOLVERS.

if (fitness.gt.9.9d7 .or. isnan(fitness) ) fitness=1.d30
END SUBROUTINE COST


END MODULE COST_MODULE
