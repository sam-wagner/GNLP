MODULE COST_MODULE
IMPLICIT NONE


CONTAINS

! THE COST FUNCTION SHOULD BE INSERTED HERE

SUBROUTINE COST(N_DOUBLE, N_INT, N1, N2, CHROM_DOUBLE, CHROM_INT, FITNESS, ARRAY, G_CON, NCON)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N_INT, N_DOUBLE, N1, N2, CHROM_INT(N_INT), NCON
DOUBLE PRECISION, INTENT(IN) :: CHROM_DOUBLE(N_DOUBLE)
DOUBLE PRECISION, INTENT(INOUT) :: FITNESS, G_CON(NCON), ARRAY(N1,N2)

DOUBLE PRECISION :: X(N_DOUBLE), X1, X2, SUM1, SUM2
INTEGER :: I


ARRAY(1,2)=ARRAY(1,2)+1.D0

X=CHROM_DOUBLE

X1=X(1)
X2=X(2)


sum1 = 0.D0
sum2 = 0.D0

DO i = 1,5,1
	SUM1=SUM1 + DBLE(I) * COS(DBLE(I+1)*x1+DBLE(I));
	SUM2=SUM2 + DBLE(I) * COS(DBLE(I+1)*x2+DBLE(I));
END DO

FITNESS=SUM1*SUM2+200.D0

END SUBROUTINE COST

END MODULE COST_MODULE
