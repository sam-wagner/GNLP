MODULE COST_MODULE
IMPLICIT NONE


CONTAINS

! THE COST FUNCTION SHOULD BE INSERTED HERE

SUBROUTINE COST(N_DOUBLE, N_INT, N1, N2, CHROM_DOUBLE, CHROM_INT, FITNESS, ARRAY, G_CON, NCON)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N_INT, N_DOUBLE, N1, N2, CHROM_INT(N_INT), NCON
DOUBLE PRECISION, INTENT(IN) :: CHROM_DOUBLE(N_DOUBLE)
DOUBLE PRECISION, INTENT(INOUT) :: FITNESS, G_CON(NCON), ARRAY(N1,N2)

DOUBLE PRECISION :: X(N_DOUBLE), DUM, PI
INTEGER :: I

pi=4.d0*atan(1.d0)
ARRAY(1,2)=ARRAY(1,2)+1.D0

X=CHROM_DOUBLE
DUM=0.D0
DO I=1,N_DOUBLE,1
    DUM=DUM+(X(I))**2-10.D0*COS(2.D0*PI*(X(I)))
END DO    
    

FITNESS=10.D0*dble(n_double)+DUM

IF (ncon.gt.0) then
    g_con=-1.d0
end if

END SUBROUTINE COST

END MODULE COST_MODULE
