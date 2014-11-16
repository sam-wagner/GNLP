MODULE UNCMIN_MOD
USE COST_MODULE
IMPLICIT NONE
CONTAINS

!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!************************                                             ***********************!
!************************         UNCMIN NLP SOLVER ALGORITM          ***********************!
!************************                                             ***********************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!

!   THIS VERSION OF UNCMIN WAS MODIFIED TO ALLOWS THE INTEGER CHROMOSOME TO BE PASSED TO THE
!   COST FUNCTION.  SOME UPDATING WAS DONE IN THE HOPES OF MAKING DEVELOPING A CUDA OR OPENCL
!   VERSION OF UNCMIN EASIER WHEN THE TIME COMES.  
!   LAST UPDATED 9/2013 BY:
!       SAMUEL WAGNER
!           PHD CANDIDATE AT IOWA STATE UNIVERSITY OF SCIENCE AND TECHNOLOGY
!           DEPARTMENT OF AEROSPACE ENGINEERING
!           THEWAGS@IASTATE.EDU OR THEWAGS.05@OUTLOOK.COM
!           www.linkedin.com/in/samuelawagner/
!
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
! BAKSLV solves A'*x=b where A is a lower triangular matrix.
!
!  Discussion:
!
!    BAKSLV solves the linear equations A'*X=B, where A is a
!    lower triangular matrix and A' is the transpose of A.
!
!    This routine is required by the UNCMIN minimization program.
!
!    If B is no longer required by calling routine, then vectors B and
!    X may share the same storage, and the output value of X will
!    overwrite the input value of B.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the number of rows and columns in A.
!
!    Input, real ( kind = 8 ) A(NR,N), the N by N matrix, containing the lower
!    triangular matrix.  A is not altered by this routine.
!
!    Output, real ( kind = 8 ) X(N), the solution vector.
!
!    Input, real ( kind = 8 ) B(N), the right hand side vector.
!
SUBROUTINE BAKSLV (NR, N, A, X, B)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, NR
DOUBLE PRECISION, INTENT(IN) :: A(NR,N), B(N)
DOUBLE PRECISION, INTENT(INOUT) :: X(N)
INTEGER :: I, IP1

!  Solve L' * x = b.

I = N
X(I)=B(I)/A(I,I)

IF (n==1) THEN
    RETURN
END IF

DO
    IP1=I
    I=I-1
    X(I)=(B(I)-DOT_PRODUCT(X(IP1:N), A(IP1:N,I)))/A(I,I)
    IF ( I == 1 ) THEN
        EXIT
    END IF
END DO

END SUBROUTINE BAKSLV



!*****************************************************************************80
!
!! CHLHSN finds the L*L' decomposition of the perturbed model hessian matrix.
!
!  Discussion:
!
!    The perturbed model Hessian matrix has the form
!
!      A + MU * I
!
!    (where 0 <= MU and I is the identity matrix) which is safely
!    positive definite.
!
!    If A is safely positive definite upon entry, then MU=0.
!
!    1. If A has any negative diagonal elements, then choose 0 < MU
!    such that the diagonal of A:=A+MU*I is all positive
!    with the ratio of its smallest to largest element on the
!    order of sqrt ( EPSM ).
!
!    2. A undergoes a perturbed Cholesky decomposition which
!    results in an LL+ decomposition of A+D, where D is a
!    non-negative diagonal matrix which is implicitly added to
!    A during the decomposition if A is not positive definite.
!    A is retained and not changed during this process by
!    copying L into the upper triangular part of A and the
!    diagonal into UDIAG.  Then the Cholesky decomposition routine
!    is called.  On return, ADDMAX contains the maximum element of D.
!
!    3. If ADDMAX=0, A was positive definite going into step 2
!    and return is made to calling program.  Otherwise,
!    the minimum number SDD which must be added to the
!    diagonal of A to make it safely strictly diagonally dominant
!    is calculated.  Since A + ADDMAX * I and A + SDD * I are safely
!    positive definite, choose MU = min ( ADDMAX, SDD ) and decompose
!    A + MU * I to obtain L.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input/output, real A(NR,N), contains an N by N matrix.
!    On input, A is the model hessian.  Only the lower triangular part and
!    diagonal are stored.  On output, A contains the factor L of the
!    LL+ decomposition of the perturbed model hessian in the lower triangular
!    part and diagonal, and contains the hessian in the upper triangular part
!    and UDIAG.
!
!    Input, real ( kind = 8 ) EPSM, the machine epsilon.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Output, real ( kind = 8 ) UDIAG(N), the diagonal of the hessian.
!
!  Local variables:
!
!    tol              tolerance
!    diagmn           minimum element on diagonal of a
!    diagmx           maximum element on diagonal of a
!    offmax           maximum off-diagonal element of a
!    offrow           sum of off-diagonal elements in a row of a
!    evmin            minimum eigenvalue of a
!    evmax            maximum eigenvalue of a
!
SUBROUTINE CHLHSN ( nr, n, a, epsm, sx, udiag )
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, NR
DOUBLE PRECISION, INTENT(IN) :: EPSM, SX(N)
DOUBLE PRECISION, INTENT(INOUT) :: A(NR,N), UDIAG(N)
DOUBLE PRECISION :: ADDMAX, AMU, DIAGMX, DIAGMN, EVMAX, EVMIN, OFFMAX, OFFROW
DOUBLE PRECISION :: POSMAX, SDD, TOL
INTEGER :: I,J

!  Scale the hessian.
DO j = 1, n
    DO i = j, n
        A(i,j) = A(i,j) / ( SX(i) * SX(j) )
    END DO
END DO

!  Step1
tol = sqrt ( epsm )

diagmx = a(1,1)
diagmn = a(1,1)

DO i = 2, n
    IF ( a(i,i) < diagmn ) THEN
        diagmn = a(i,i)
    END IF
    IF(diagmx<A(i,i))THEN
        diagmx=A(i,i)
    END IF
END DO

posmax = max ( diagmx, 0.0D+00 )

IF ( diagmn <= posmax * tol ) THEN

    amu = tol * ( posmax - diagmn ) - diagmn
!
!  Find the largest off-diagonal element of A.
!
    IF ( amu == 0.0D0 ) THEN

      offmax = 0.0D0

        DO i = 2, n
            DO j = 1, i-1
                IF(offmax<abs(A(i,j))) THEN
                    offmax=abs ( a(i,j) )
                END IF
            END DO
        END DO

        amu = offmax

        IF ( amu == 0.0D+00 ) THEN
            amu = 1.0D+00
        ELSE
            amu = amu * ( 1.0D+00 + tol )
        END IF

    END IF

!  A = A + MU*I
    DO i = 1, n
        A(i,i) = A(i,i) + AMU
    END DO

    diagmx = diagmx + amu

END IF
!
!  Step2
!
!  Copy lower triangular part of A to upper triangular part
!  and diagonal of A to udiag
!
  do j = 1, n
    udiag(j) = a(j,j)
    do i = j + 1, n
      a(j,i) = a(i,j)
    end do
  end do

  call choldc ( nr, n, a, diagmx, tol, addmax )
!
!  Step3
!
!  If ADDMAX=0, A was positive definite going into step 2,
!  the ll+ decomposition has been done, and we return.
!
!  Otherwise, 0 < ADDMAX.  perturb A so that it is safely
!  diagonally dominant and find ll+ decomposition
!
  if ( 0.0D+00 < addmax ) then
!
!  Restore original A (lower triangular part and diagonal)
!
    do j = 1, n
      a(j,j) = udiag(j)
      do i = j+1, n
        a(i,j) = a(j,i)
      end do
    end do
!
!  Find SDD such that A+sdd*i is safely positive definite
!  note:  evmin<0 since A is not positive definite;
!
    evmin = 0.0D+00
    evmax = a(1,1)

    do i = 1, n

      offrow = sum ( abs ( a(i,1:i-1) ) ) + sum ( abs ( a(i+1:n,i) ) )
      evmin = min ( evmin, a(i,i)-offrow )
      evmax = max ( evmax, a(i,i)+offrow )

    end do

    sdd = tol * ( evmax - evmin ) - evmin
!
!  Perturb A and decompose again.
!
    amu = min ( sdd, addmax )

    do i = 1, n
      a(i,i) = a(i,i) + amu
      udiag(i) = a(i,i)
    end do
!
!  A is now guaranteed safely positive definite
!
    call choldc ( nr, n, a, 0.0D+00, tol, addmax )

  end if
!
!  Unscale the hessian and Cholesky decomposition matrix.
!
  do j = 1, n

    a(j:n,j) = sx(j:n) * a(j:n,j)

    do i = 1, j - 1
      a(i,j) = sx(i) * sx(j) * a(i,j)
    end do

    udiag(j) = udiag(j) * sx(j) * sx(j)

  end do

  return
end subroutine chlhsn






!*****************************************************************************80
!
!! CHOLDC finds the perturbed L*L' decomposition of A+D.
!
!  Discussion:
!
!    D is a non-negative diagonal matrix added to A if
!    necessary to allow the Cholesky decomposition to continue.
!
!    The normal Cholesky decomposition is performed.  However, if at any
!    point the algorithm would attempt to set
!      L(I,I) = sqrt ( TEMP )
!    with
!      TEMP < TOL * DIAGMX,
!    then L(I,I) is set to sqrt ( TOL * DIAGMX )
!    instead.  This is equivalent to adding TOL * DIAGMX-TEMP to A(I,I)
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input/output, real ( kind = 8 ) A(NR,N), the N by N matrix.
!    On input, the matrix for which to find the perturbed
!    Cholesky decomposition.
!    On output, the lower triangular part contains the L factor,
!    and the diagonal of A.
!
!    Input, real ( kind = 8 ) DIAGMX, the maximum diagonal element of A.
!
!    Input, real ( kind = 8 ) TOL, a tolerance.
!
!    Output, real ( kind = 8 ) ADDMAX, the maximum amount implicitly added to
!    the diagonal of A in forming the Cholesky decomposition of A+D.
!
!  Local variables:
!
!    aminl    smallest element allowed on diagonal of L.
!
!    amnlsq   =aminl**2
!
!    offmax   maximum off-diagonal element in column of a
!
SUBROUTINE CHOLDC ( nr, n, a, diagmx, tol, addmax )
IMPLICIT NONE
INTEGER, INTENT(IN) :: NR, N
DOUBLE PRECISION :: DIAGMX, TOL
DOUBLE PRECISION, INTENT(INOUT) :: A(NR,N), ADDMAX
DOUBLE PRECISION :: AMINL, AMNLSQ, OFFMAX, SUM2, TEMP
INTEGER :: I, J, K

ADDMAX=0.0D0
AMINL=sqrt(DIAGMX*TOL)
AMNLSQ=AMINL**2

!  Form column J of L.

DO J = 1, n
!  Find diagonal elements of L.
    SUM2=SUM(A(j,1:j-1)**2)
    TEMP=A(j,j)-SUM2
    IF(amnlsq<=temp)THEN
        A(j,j)=SQRT(TEMP)
! Find maximum off-diagonal element in column.
    ELSE
        OFFMAX=0.0D0

        DO i = j+1, n
            IF (OFFMAX<abs(A(i,j))) THEN
                OFFMAX=abs(A(i,j))
            END IF
        END DO

        IF (OFFMAX<=AMNLSQ) THEN
            OFFMAX=AMNLSQ
        END IF
!  Add to diagonal element to allow Cholesky decomposition to continue
        A(j,j)=SQRT(OFFMAX)
        ADDMAX=MAX(addmax, offmax-temp)

    END IF
!  Find (I,J) element of lower triangular matrix.
    DO I = J+1, N
        SUM2=0.0D0
        DO K=1,J-1
            SUM2=SUM2+A(I,K)*A(J,K)
        END DO
        A(I,J)=(A(I,J)-SUM2)/A(J,J)
    END DO
END DO

END SUBROUTINE CHOLDC


!*****************************************************************************80
!
!  DFAULT sets default values for the optimization algorithm.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), an initial guess for the solution,
!    used to compute a maximum stepsize.
!
!    Output, real ( kind = 8 ) TYPSIZ(N), a typical size for each component
!    of X.
!
!    Output, real ( kind = 8 ) FSCALE, an estimate of the scale of the
!    minimization function.
!
!    Output, integer ( kind = 4 ) METHOD, specifies the algorithm to use to
!    solve the minimization problem.
!
!    Output, integer ( kind = 4 ) IEXP, set to 0 if minimization function
!    not expensive to evaluate.
!
!    Output, integer ( kind = 4 ) MSG, a message to inhibit certain automatic
!    checks and output.
!
!    Output, integer ( kind = 4 ) NDIGIT, the number of good digits in
!    minimization function.
!
!    Output, integer ( kind = 4 ) ITNLIM, the maximum number of allowable
!    iterations.
!
!    Output, integer ( kind = 4 ) IAGFLG, set to 0, meaning the analytic
!    gradient is not supplied.
!
!    Output, integer ( kind = 4 ) IAHFLG, set to 0, meaning the analytic hessian is
!    not supplied.
!
!    Output, integer ( kind = 4 ) IPR, the device to which to send output.
!
!    Output, real ( kind = 8 ) DLT, the trust region radius.
!
!    Output, real ( kind = 8 ) GRADTL, a tolerance at which the gradient is
!    considered close enough to zero to terminate algorithm.
!
!    Output, real ( kind = 8 ) STEPMX, the maximum stepsize, set to 0.0 to trip
!    the default maximum in OPTCHK.
!
!    Output, real ( kind = 8 ) STEPTL, a tolerance at which successive
!    iterates are considered close enough to terminate the algorithm.
!
SUBROUTINE DFAULT(N, X, typsiz, fscale, method, iexp, msg, ndigit, itnlim, &
                  iagflg, iahflg, ipr, dlt, gradtl, stepmx, steptl )
IMPLICIT NONE
INTEGER, INTENT(IN) :: N
INTEGER, INTENT(INOUT) :: METHOD, IEXP, MSG, NDIGIT, ITNLIM, IAGFLG, IAHFLG, IPR
DOUBLE PRECISION, INTENT(IN) :: X(N)
DOUBLE PRECISION, INTENT(INOUT) :: TYPSIZ(N), FSCALE, DLT, GRADTL, STEPMX, STEPTL
DOUBLE PRECISION :: EPSM

!  Typical size of X and minimization function.
TYPSIZ(1:N)=1.0D0
FSCALE=1.0D0

!  Tolerances.
DLT=-1.0D0
EPSM=EPSILON(EPSM)
GRADTL=EPSM**(1.D0/3.D0)
STEPMX=0.0D0
STEPTL=SQRT(EPSM)

!  Flags.
METHOD=1
IEXP=1
MSG=9
NDIGIT=-1
ITNLIM=150
IAGFLG=0
IAHFLG=0
IPR=6

END SUBROUTINE DFAULT






!*****************************************************************************80
!
!! DNRM2 returns the euclidean norm of a vector.
!
!  Discussion:
!
!     DNRM2 ( X ) = sqrt ( X' * X )
!
!  Author:
!
!    Sven Hammarling
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input, real ( kind = 8 ) X(*), the vector whose norm is to be computed.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries of X.
!
!    Output, real ( kind = 8 ) DNRM2, the Euclidean norm of X.
!
FUNCTION DNRM2(N, X, INCX)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, INCX
DOUBLE  PRECISION, INTENT(IN) :: X(N)
DOUBLE PRECISION :: ABSXI, DNRM2, NORM, SCALE, SSQ
INTEGER :: IX

IF (N<1 .or. INCX<1) THEN
    NORM=0.D0
ELSE IF (n==1) THEN
    NORM=ABS(X(1))
ELSE
    SCALE=0.D0
    SSQ=1.D0

    DO IX=1,1+(N-1)*INCX, INCX
        IF(X(IX)/=0.D0) THEN
            ABSXI=ABS(X(IX))
            IF (SCALE<absxi) THEN
                ssq = 1.0D+00 + ssq * ( scale / absxi )**2
                scale = absxi
            ELSE
                ssq = ssq + ( absxi / scale )**2
            END IF
        END IF
    END DO
    NORM=SCALE*SQRT(SSQ)
END IF

DNRM2=NORM

END FUNCTION DNRM2


!*****************************************************************************80
!
!  DOGDRV finds the next Newton iterate by the double dogleg method.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the old iterate, "X[K-1]".
!
!    Input, real ( kind = 8 ) F, the function value at the old iterate, "F(X)".
!
!    Input, real ( kind = 8 ) G(N), the gradient at the old iterate.
!
!    Input, real ( kind = 8 ) A(N,N), the Cholesky decomposition of the
!    Hessian matrix in lower triangular part and diagonal.
!
!    Input, real ( kind = 8 ) P(N), the Newton step.
!
!    Output, real ( kind = 8 ) XPLS(N), the new iterate "X[K]".
!
!    Output, real ( kind = 8 ) FPLS, the function value at the new iterate,
!    F(XPLS).
!
!    Input, external COST, the name of the subroutine to evaluate the function,
!    of the form
!
!      subroutine COST ( n, x, f )
!      integer ( kind = 4 ) n
!      real x(n)
!      real f
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, real ( kind = 8 ) STEPMX, the maximum allowable step size.
!
!    Input, real ( kind = 8 ) STEPTL, the relative step size at which
!    successive iterates are considered close enough to terminate algorithm.
!
!    Input/output, real ( kind = 8 ) DLT, the trust region radius.
!    [retain value between successive calls].
!
!    Output, integer ( kind = 4 ) IRETCD, the return code.
!    0, satisfactory XPLS found
!    1, failed to find satisfactory XPLS sufficiently distinct from X.
!
!    Output, logical MXTAKE, TRUE if a step of maximum length was used.
!
!    Workspace, real ( kind = 8 ) SC(N), holds the current step.
!
!    Workspace, real ( kind = 8 ) WRK1(N).
!
!    Workspace, real ( kind = 8 ) WRK2(N).
!
!    Workspace, real ( kind = 8 ) WRK3(N).
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
SUBROUTINE DOGDRV( NR, N, X, F, G, A, P, XPLS, FPLS, SX, STEPMX, STEPTL, &
                   DLT, IRETCD, MXTAKE, SC, WRK1, WRK2, WRK3, IPR,&
                   N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
IMPLICIT NONE

INTEGER, INTENT(IN) :: NR, N, IPR, NCON
INTEGER, INTENT(INOUT) :: IRETCD
DOUBLE PRECISION, INTENT(IN) :: X(N), SX(N), F, G(N), A(NR,N), P(N)
DOUBLE PRECISION, INTENT(IN) :: STEPMX, STEPTL
DOUBLE PRECISION, INTENT(INOUT) :: XPLS(N), FPLS, DLT, G_CON(NCON)
DOUBLE PRECISION, INTENT(INOUT) :: SC(N), WRK1(N), WRK2(N), WRK3(N)
LOGICAL, INTENT(INOUT) :: MXTAKE
!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER, INTENT(IN) :: N_INT, N1, N2, CHROM_INT(N_INT), INFO
DOUBLE PRECISION, INTENT(IN) :: INPUT_ARRAY(N1,N2)
!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION :: CLN, ETA, FPLSP, RNWTLN
LOGICAL :: FSTDOG, NWTAKE

IRETCD=4
FSTDOG=.TRUE.

RNWTLN=sqrt(sum(sx(1:n)**2 * p(1:n)**2))

DO
    !
    !  Find new step by double dogleg algorithm.
    !
    CALL DOGSTP(NR,N,G,A,P,SX,RNWTLN,DLT,NWTAKE,FSTDOG,WRK1,WRK2,CLN,ETA,SC,IPR,STEPMX)

    !
    !  Check new point and update trust region.
    !

    CALL TREGUP(NR, N, X, F, G, A, SC, SX, NWTAKE, STEPMX, STEPTL, DLT, IRETCD, &
                WRK3, FPLSP, XPLS, FPLS, MXTAKE, IPR, 2, WRK1, N_INT, N1, N2, &
                CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)

    IF(iretcd<=1)THEN
        EXIT
    END IF

END DO

RETURN

END SUBROUTINE DOGDRV


!*****************************************************************************80
!
!! DOGSTP finds a new step by the double dogleg algorithm.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) G(N), the gradient at the current iterate.
!
!    Input, real ( kind = 8 ) A(NR,N), the Cholesky decomposition of the
!    hessian in the lower triangle and diagonal.
!
!    Input, real ( kind = 8 ) P(N), the Newton step.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, real ( kind = 8 ) RNWTLN, the Newton step length.
!
!    Input/output, real ( kind = 8 ) DLT, the trust region radius.
!
!    Input/output, logical NWTAKE, TRUE if a Newton step was taken.
!
!    Input/output, logical FSTDOG, TRUE if on first leg of dogleg.
!
!    Input/output, real ( kind = 8 ) SSD(N), workspace [cauchy step to
!    the minimum of the quadratic model in the scaled steepest descent
!    direction] [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) V(N), workspace  [retain value
!    between successive calls]
!
!    Workspace, real ( kind = 8 ) CLN, the cauchy length.
!    [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) ETA, [retain value between successive calls]
!
!    Output, real ( kind = 8 ) SC(N), the current step.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
!    Input, real ( kind = 8 ) STEPMX, the maximum allowable step size.
!
!  Local variables:
!
!    CLN, the length of cauchy step
!
SUBROUTINE dogstp ( nr, n, g, a, p, sx, rnwtln, dlt, nwtake, fstdog, ssd, v, &
  cln, eta, sc, ipr, stepmx )
IMPLICIT NONE
INTEGER, INTENT(IN) :: NR, N, IPR
DOUBLE PRECISION, INTENT(IN):: G(N), P(N), SX(N), A(NR,N), RNWTLN, STEPMX
DOUBLE PRECISION, INTENT(INOUT) :: SC(N), SSD(N), V(N), DLT, CLN, ETA
LOGICAL, INTENT(INOUT) :: FSTDOG, NWTAKE
DOUBLE PRECISION :: ALAM, ALPHA, BETA, DOT1, DOT2, TMP
INTEGER :: I, J

!  Can we take a Newton step?

IF( rnwtln <= dlt )THEN

    nwtake = .true.
    sc(1:n) = p(1:n)
    dlt = rnwtln

ELSE
!
!  The Newton step is too long.
!  The Cauchy step is on double dogleg curve.
!
    nwtake = .false.

    IF(fstdog) THEN
!
!  Calculate double dogleg curve, SSD.
!
        fstdog = .false.
        alpha = sum ( ( g(1:n) / sx(1:n) )**2 )
        beta = 0.0D+00
        DO i = 1, n
            tmp = 0.0D0
            DO j = i, n
                tmp = tmp + ( a(j,i) * g(j) ) / ( sx(j) * sx(j) )
            END DO
            beta = beta + tmp * tmp
        END DO

        ssd(1:n) = - ( alpha / beta ) * g(1:n) / sx(1:n)

        cln = alpha * sqrt ( alpha ) / beta

        eta = 0.2D+00 + ( 0.8D+00 * alpha * alpha ) / &
             ( - beta * dot_product ( g, p ) )

        v(1:n) = eta * sx(1:n) * p(1:n) - ssd(1:n)

        IF ( dlt == - 1.0D+00 ) THEN
            dlt = min ( cln, stepmx )
        END IF

    END IF
!
!  Take a partial step in the Newton direction.
!
    IF ( eta * rnwtln <= dlt ) THEN

        sc(1:n) = ( dlt / rnwtln ) * p(1:n)
!
!  Take a step in steepest descent direction.
!
    ELSE IF ( dlt <= cln ) THEN

        sc(1:n) = ( dlt / cln ) * ssd(1:n) / sx(1:n)
!
!  Convex combination of SSD and eta*p which has scaled length DLT.
!
    ELSE

        dot1 = dot_product ( v, ssd )
        dot2 = dot_product ( v, v )
        alam = ( -dot1 + sqrt ( ( dot1 * dot1 ) &
                 - dot2 * ( cln * cln - dlt * dlt ) ) ) / dot2
        sc(1:n) = ( ssd(1:n) + alam * v(1:n) ) / sx(1:n)
    END IF
END IF

end subroutine dogstp






!*****************************************************************************80
!
! FORSLV solves A*x=b where A is lower triangular matrix.
!
SUBROUTINE forslv ( nr, n, a, x, b )
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, NR
DOUBLE PRECISION, INTENT(IN) :: A(NR,N), B(N)
DOUBLE PRECISION, INTENT(INOUT) :: X(N)
INTEGER :: I

X(1) = B(1) / A(1,1)

DO i = 2, n
    X(i) = (B(i) - dot_product( A(i,1:i-1), X(1:i-1) ) ) / A(i,i)
END DO

END SUBROUTINE forslv






!*****************************************************************************80
!
!! FSTOCD approximates the gradient of a function using central differences.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the point at which the gradient is to
!    be approximated.
!
!    Input, real ( kind = 8 ) SX(N), the scaling factors for X.
!
!    Input, real ( kind = 8 ) RNOISE, the relative noise in COST [F(X)].
!
!    Output, real ( kind = 8 ) G(N), a central difference approximation
!    to the gradient.
!
SUBROUTINE fstocd ( n, x, sx, rnoise, g, N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON )
IMPLICIT NONE

INTEGER, INTENT(IN) :: NCON
DOUBLE PRECISION, INTENT(INOUT) :: G_CON(NCON)

INTEGER, INTENT(IN) :: N
DOUBLE PRECISION, INTENT(IN) :: SX(N), RNOISE
DOUBLE PRECISION, INTENT(IN OUT) :: G(N), X(N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ADDITIONAL VARIABLE CHANGES
INTEGER, INTENT(IN) :: N_INT, N1, N2, CHROM_INT(N_INT), INFO
DOUBLE PRECISION, INTENT(IN) :: INPUT_ARRAY(N1,N2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION :: FMINUS, FPLUS, STEPI, THIRD, XTEMPI
INTEGER :: I


!  Find I-th stepsize, evaluate two neighbors in direction of I-th
!  unit vector, and evaluate I-th component of gradient.

  third = 1.0D0 / 3.0D0

DO i = 1, n
    stepi = rnoise**third * max ( abs ( x(i) ), 1.0D+00 / sx(i) )
    xtempi = x(i)
    x(i) = xtempi + stepi
    !call COST ( n, x, fplus )
    CALL COST(N, N_INT, N1, N2, X, CHROM_INT, FPLUS, INPUT_ARRAY, G_CON, NCON)
    x(i) = xtempi - stepi
    !call COST ( n, x, fminus )
    CALL COST(N, N_INT, N1, N2, X, CHROM_INT, FMINUS, INPUT_ARRAY, G_CON, NCON)
    x(i) = xtempi
    g(i) = ( fplus - fminus ) / ( 2.0D+00 * stepi )
END DO

END SUBROUTINE fstocd






!*****************************************************************************80
!
!! FSTOFD approximates a derivative by a first order approximation.
!
!  Discussion:
!
!    The routine finds the first order forward finite difference
!    approximation A to the first derivative of the function defined
!    by the subprogram "fname" evaluated at the new iterate "xpls".
!
!    For optimization, use this routine to estimate:
!
!    * the first derivative (gradient) of the optimization function "COST"
!      if no analytic user routine has been supplied;
!
!    * the second derivative (hessian) of the optimization function
!      if no analytic user routine has been supplied for the hessian but
!      one has been supplied for the gradient ("COST") and if the
!      optimization function is inexpensive to evaluate.
!
!    m=1 (optimization) algorithm estimates the gradient of the function
!    (COST).   COST(x) # f: r(n)-->r(1)
!
!    m=n (systems) algorithm estimates the jacobian of the function
!    COST(x) # f: r(n)-->r(n).
!
!    m=n (optimization) algorithm estimates the hessian of the optimization
!    function, where the hessian is the first derivative of "COST"
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix A.
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A, and the dimension
!    of the problem.
!
!    Input, real ( kind = 8 ) XPLS(N), the point at which the derivative is
!    to be estimated.
!
!    Input, external COST, the name of the subroutine to evaluate
!    the function, of the form
!
!      subroutine COST ( n, x, f )
!      integer ( kind = 4 ) n
!      real x(n)
!      real f
!
!    Input, real ( kind = 8 ) FPLS(M).
!    If M is 1, (optimization), then this is the function value at
!    the new iterate.
!    If M = N for optimization, then this is the value of the first
!    derivative of the function.
!    If M = N for nonlinear systems, then this is the value
!    of the associated minimization function.
!
!    Output, real ( kind = 8 ) A(NR,N), the N by N finite difference
!    approximation.  Only the lower triangular matrix and diagonal are
!    computed.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, real ( kind = 8 ) RNOISE, the relative noise or inaccuracy in the
!    function value COST.
!
!    Workspace, real ( kind = 8 ) FHAT(M).
!
!    Input, integer ( kind = 4 ) ICASE, problem specifier:
!    1, optimization (gradient)
!    2, systems
!    3, optimization (hessian)
!
!  Local variables:
!
!    real STEPSZ, the stepsize in the J-th variable direction
!
SUBROUTINE FSTOFD(NR,M,N,XPLS,FPLS,A,SX,RNOISE,FHAT,ICASE, N_INT, N1, N2, &
                  CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
IMPLICIT NONE

INTEGER, INTENT(IN) :: NCON
DOUBLE PRECISION, INTENT(INOUT) ::G_CON(NCON)

INTEGER, INTENT(IN) :: M, N, NR, ICASE
DOUBLE PRECISION, INTENT(IN) :: FPLS(M), SX(N), RNOISE
DOUBLE PRECISION, INTENT(INOUT) :: A(NR,N), XPLS(N), FHAT(M)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NEW INPUT VARIABLES
INTEGER, INTENT(IN) :: N_INT, N1, N2, INFO, CHROM_INT(N_INT)
DOUBLE PRECISION :: INPUT_ARRAY(N1,N2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER :: I, J
DOUBLE PRECISION :: STEPSZ, XTMPJ
DOUBLE PRECISION :: F_DUM
!  Find the J-th column of A.
!  Each column is the derivative of f(COST) with respect to xpls(j).
!
DO j = 1, n
    stepsz = sqrt ( rnoise ) * max ( abs ( xpls(j) ), 1.0D0 / sx(j) )
    xtmpj = xpls(j)
    xpls(j) = xtmpj + stepsz
    !CALL COST( n, xpls, fhat(1) )


    !!WRITE(*,*) "made it to COST inside fstofd"
    !!WRITE(*,*) N, N_INT,N1,N2,INFO
    !!WRITE(*,*) XPLS, CHROM_INT
    !!WRITE(*,*) FHAT(1)
    !!WRITE(*,*) INPUT_ARRAY

    CALL COST(N, N_INT, N1, N2, XPLS, CHROM_INT, FHAT(1), INPUT_ARRAY, G_CON, NCON)
    !!WRITE(*,*) "MADE IT OUT OF COST IN FSTOFD"
    !!WRITE(*,*) "F_DUM", F_DUM
    !FHAT(1)=F_DUM

    !!WRITE(*,*) "m", m
    !!WRITE(*,*) "fhat", fhat(1)

    fhat=fhat(1)

    xpls(j) = xtmpj
    a(1:m,j) = ( fhat(1:m) - fpls(1:m) ) / stepsz
    !!WRITE(*,*) a
END DO

IF( icase /= 3 ) then
    RETURN
END IF
!
!  If computing the hessian, A must be symmetric.
!
DO j = 1, n-1
    DO i = j+1, m
        A(i,j) = ( A(i,j) + A(j,i) ) / 2.0D0
    END DO
END DO

END SUBROUTINE fstofd






!*****************************************************************************80
!
!! GRDCHK checks an analytic gradient against an estimated gradient.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), a point at which the gradient is
!       to be checked.
!
!    Input, external COST, the name of the subroutine that evaluates
!    the optimization function, of the form
!
!      subroutine COST ( n, x, f )
!      integer ( kind = 4 ) n
!      real f
!      real x(n)
!
!    Input, real ( kind = 8 ) F, the function value at X.
!
!    Input, real ( kind = 8 ) G(N), the gradient value at X.
!
!    Input, real ( kind = 8 ) TYPSIZ(N), a typical size for each component of X.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling values:
!    SX(1:N)=1.0/TYPSIZ(1:N)
!
!    Input, real ( kind = 8 ) FSCALE, an estimate of the scale of the
!    objective function COST.
!
!    Input, real ( kind = 8 ) RNF, the relative noise in the optimization
!    function COST.
!
!    Input, real ( kind = 8 ) ANALTL, a tolerance for comparison of estimated
!    and analytical gradients
!
!    Workspace, real ( kind = 8 ) WRK1(N).
!
!    Output, integer ( kind = 4 ) MSG, message or error code.
!      0: no error detected.
!    -21: probable coding error of gradient
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
SUBROUTINE GRDCHK(N,X,F,G,TYPSIZ,SX,FSCALE,RNF,ANALTL,WRK1,MSG,IPR, &
                  N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
IMPLICIT NONE

INTEGER, INTENT(IN) ::N, NCON
INTEGER, INTENT(INOUT) :: IPR, MSG
DOUBLE PRECISION, INTENT(IN OUT) :: WRK1(N), X(N), G_CON(NCON)
DOUBLE PRECISION, INTENT(IN) :: F, G(N), TYPSIZ(N), SX(N), FSCALE
DOUBLE PRECISION, INTENT(IN) :: RNF, ANALTL
!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER, INTENT(IN) :: N_INT, N1, N2, CHROM_INT(N_INT), INFO
DOUBLE PRECISION, INTENT(IN) :: INPUT_ARRAY(N1,N2)
!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION :: GS, VALUE(1), WRK(1)
INTEGER :: I, KER

MSG=0

!
!  Compute the first order finite difference gradient;
!  compare it to the analytic gradient.
!
  value(1) = f
  !call fstofd ( 1, 1, n, x, value, wrk1, sx, rnf, wrk, 1 )
  CALL FSTOFD(1, 1, N, X, VALUE, WRK1, SX, RNF, WRK, 1, &
              N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
  ker = 0

  do i = 1, n

    gs = max ( abs ( f ), fscale ) / max ( abs ( x(i) ), typsiz(i) )

    if ( max ( abs ( g(i) ), gs ) * analtl < abs ( g(i) - wrk1(i) ) ) then
      ker = 1
    end if

  end do

  if ( ker /= 0 ) then
    !write ( ipr, * ) ' '
    !write ( ipr, * ) 'GRDCHK - probable error in analytic gradient.'
    !write ( ipr, * ) ' '
    !write ( ipr, * ) ' grdchk     comp            analytic            est'
    !write ( ipr, 902 ) ( i, g(i), wrk1(i), i = 1, n )
    msg = -21
  end if

  return

  902 format(' grdchk    ',i5,3x,e20.13,3x,e20.13)
END SUBROUTINE GRDCHK






!*****************************************************************************80
!
!! HESCHK checks an analytic hessian against a computed estimate.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), a point at which the Hessian is to
!    be checked.
!
!    Input, external COST, the name of the subroutine that evaluates
!    the optimization function, of the form:
!
!      subroutine COST ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!
!    Input, external D1FCN, the name of the subroutine to evaluate the
!    gradient of the function, of the form:
!
!      subroutine d1fcn ( n, x, g )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) g(n)
!      real ( kind = 8 ) x(n)
!
!    Input, external D2FCN, the name of the subroutine to evaluate the
!    Hessian of the function, of the form:
!
!      subroutine d2fcn ( nr, n, x, h )
!      integer ( kind = 4 ) nr
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) h(nr,n)
!      real ( kind = 8 ) x(n)
!
!    Input, real ( kind = 8 ) F, the function value at X.
!
!    Output, real ( kind = 8 ) G(N), the value of the gradient at X.
!
!    Output, real ( kind = 8 ) A(NR,N), the analytic Hessian matrix will
!    be stored in the lower triangle and diagonal.
!
!    Input, real ( kind = 8 ) TYPSIZ(N), a typical size for each component of X.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix.
!
!    Input, real ( kind = 8 ) RNF, the relative noise in the optimization
!    function COST.
!
!    Input, real ( kind = 8 ) ANALTL, a tolerance for comparison of the
!    estimated and analytic gradients.
!
!    Input, integer ( kind = 4 ) IAGFLG, is 1 if the analytic gradient is supplied.
!
!    Workspace, real ( kind = 8 ) UDIAG(N).
!
!    Workspace, real ( kind = 8 ) WRK1(N).
!
!    Workspace, real ( kind = 8 ) WRK2(N).
!
!    Input/output, integer ( kind = 4 ) MSG, message or error code
!    on input : if =1xx do not compare analytic + estimated hessian.
!    on output: =-22, probable coding error of hessian.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
SUBROUTINE HESCHK(nr, n, x, f, g, a, typsiz, sx, rnf, &
                  analtl, iagflg, udiag, wrk1, wrk2, msg, ipr,&
                  N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N, NR, IAGFLG, IPR, NCON
INTEGER, INTENT(INOUT) :: MSG
DOUBLE PRECISION, INTENT(IN) :: ANALTL, F, TYPSIZ(N), SX(N)
DOUBLE PRECISION, INTENT(IN) :: RNF
DOUBLE PRECISION, INTENT(INOUT) :: A(NR,N), G(N), UDIAG(N), WRK1(N), WRK2(N), X(N), G_CON(NCON)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NEW INPUT VARIABLES
INTEGER, INTENT(IN) :: N_INT, N1, N2, CHROM_INT(N_INT), INFO
DOUBLE PRECISION, INTENT(IN) :: INPUT_ARRAY(N1, N2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION :: HS
INTEGER :: I, J, KER
!  integer ( kind = 4 ) n
!  integer ( kind = 4 ) nr

!  real ( kind = 8 ) a(nr,n)
!  real ( kind = 8 ) analtl
  !external d1fcn
  !external d2fcn
!  real ( kind = 8 ) f
  !external COST
!  real ( kind = 8 ) g(n)
!  real ( kind = 8 ) hs
!  integer ( kind = 4 ) i
!  integer ( kind = 4 ) iagflg
!  integer ( kind = 4 ) ipr
!  integer ( kind = 4 ) j
!  integer ( kind = 4 ) ker
!  integer ( kind = 4 ) msg
!  real ( kind = 8 ) rnf
!  real ( kind = 8 ) sx(n)
!  real ( kind = 8 ) typsiz(n)
!  real ( kind = 8 ) udiag(n)
!  real ( kind = 8 ) wrk1(n)
!  real ( kind = 8 ) wrk2(n)
!  real ( kind = 8 ) x(n)
!
!  Compute the finite difference approximation A to the hessian.
!
  if ( iagflg == 1 ) then
    !call fstofd ( nr, n, n, x, g, a, sx, rnf, wrk1, 3 )
    CALL FSTOFD(NR, N, N, X, G, A, SX, RNF, WRK1, 3, &
                N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
  else
    !call sndofd ( nr, n, x, f, a, sx, rnf, wrk1, wrk2 )
    CALL SNDOFD(NR, N, X, F, A, SX, RNF, WRK1, WRK2, N_INT, &
                N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
  end if

  ker = 0
!
!  Copy lower triangular part of A to upper triangular part
!  and diagonal of A to UDIAG.
!
  do j = 1, n
    udiag(j) = a(j,j)
    do i = j+1, n
      a(j,i) = a(i,j)
    end do
  end do
!
!  Compute analytic hessian and compare to finite difference approximation.
!
  call d2fcn ( nr, n, x, a )

  do j = 1, n

    hs = max ( abs ( g(j) ), 1.0D+00 ) / max ( abs ( x(j) ), typsiz(j) )

    if ( max ( abs ( udiag(j) ), hs ) * analtl &
       < abs ( a(j,j) - udiag(j) ) ) then
      ker = 1
    end if

    do i = j+1, n
      if ( max ( abs ( a(i,j) ), hs ) * analtl &
        < abs ( a(i,j) - a(j,i) ) ) then
        ker = 1
      end if
    end do

  end do

  if ( ker /= 0 ) then

    !write ( ipr, '(a)' ) ' '
    !write ( ipr, '(a)' ) 'HESCHK:'
    !write ( ipr, '(a)' ) '  Probable error in coding of analytic hessian.'
    !write ( ipr, '(a)' ) '            row  col              analytic' &
    !  // '              (estimate)'
    !write ( ipr, '(a)' ) ' '

    do i = 1, n
      do j = 1, i-1
        !write(ipr,902) i, j, a(i,j), a(j,i)
      end do
      !write(ipr,902) i, i, a(i,i), udiag(i)
    end do

    msg = -22

  end if

  return
  902 format('heschk    ',2i5,2x,e20.13,2x,'(',e20.13,')')
END SUBROUTINE HESCHK






!*****************************************************************************80
!
!! HOOKDR finds the next Newton iterate by the More-Hebdon method.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the old iterate, "X[K-1]".
!
!    Input, real ( kind = 8 ) F, the function value at the old iterate.
!
!    Input, real ( kind = 8 ) G(N), the gradient or an approximation, at
!    the old iterate.
!
!    Input, real ( kind = 8 ) A(NR,N), the Cholesky decomposition of hessian
!    in lower triangular part and diagonal.  Hessian in upper triangular
!    part and UDIAG.
!
!    Input, real ( kind = 8 ) UDIAG(N), the diagonal of the hessian matrix.
!
!    Input, real ( kind = 8 ) P(N), the Newton step.
!
!    Output, real ( kind = 8 ) XPLS(N), the new iterate X[K].
!
!    Output, real ( kind = 8 ) FPLS, the function value at the new iterate.
!
!    Input, external COST, the name of the subroutine to evaluate the function,
!    of the form
!
!      subroutine COST ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) x(n)
!      real ( kind = 8 ) f
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, real ( kind = 8 ) STEPMX, the maximum allowable step size.
!
!    Input, real ( kind = 8 ) STEPTL, the relative step size at which successive
!    iterates are considered close enough to terminate algorithm.
!
!    Input/output, real ( kind = 8 ) DLT, the trust region radius.
!
!    Output, integer ( kind = 4 ) IRETCD, return code
!    0, satisfactory xpls found
!    1, failed to find satisfactory xpls sufficiently distinct from x.
!
!    Output, logical MXTAKE, is TRUE if a step of maximum length was used.
!
!    Workspace, real ( kind = 8 ) AMU, [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) DLTP, [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) PHI, [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) PHIP0, [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) SC(N).
!
!    Workspace, real ( kind = 8 ) XPLSP(N).
!
!    Workspace, real ( kind = 8 ) WRK0(N).
!
!    Input, real ( kind = 8 ) EPSM, the machine epsilon
!
!    Input, integer ( kind = 4 ) ITNCNT, the iteration count.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
SUBROUTINE HOOKDR( nr, n, x, f, g, a, udiag, p, xpls, fpls, sx, stepmx, &
  steptl, dlt, iretcd, mxtake, amu, dltp, phi, phip0, sc, xplsp, wrk0, epsm, &
  itncnt, ipr, N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N, NR, ITNCNT, IPR, NCON
INTEGER, INTENT(INOUT) :: IRETCD
LOGICAL, INTENT(INOUT) :: MXTAKE
DOUBLE PRECISION, INTENT(IN) :: X(N), F, G(N), UDIAG(N), P(N)
DOUBLE PRECISION, INTENT(IN) :: SX(N), STEPMX, STEPTL
DOUBLE PRECISION, INTENT(IN) :: EPSM
DOUBLE PRECISION, INTENT(INOUT) ::  DLT, AMU, DLTP, PHI, PHIP0, SC(N), FPLS
DOUBLE PRECISION, INTENT(INOUT) :: XPLSP(N), WRK0(N), A(NR,N), XPLS(N), G_CON(NCON)
!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER, INTENT(IN) :: N_INT, N1, N2, CHROM_INT(N_INT), INFO
DOUBLE PRECISION, INTENT(IN) :: INPUT_ARRAY(N1,N2)

!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION :: ALPHA, BETA, FPLSP, RNWTLN, TMP
LOGICAL :: FSTIME, NWTAKE
INTEGER :: I, J

iretcd = 4
fstime = .true.

rnwtln = sqrt ( sum ( sx(1:n)**2 * p(1:n)**2 ) )
!
!  If first iteration and trust region not provided by user,
!  compute initial trust region.
!
  if ( itncnt <= 1 ) then

    amu = 0.0D+00

    if ( dlt == -1.0D+00 ) then

      alpha = sum ( ( g(1:n) / sx(1:n) )**2 )

      beta = 0.0D+00
      do i = 1, n
        tmp = 0.0D+00
        do j = i, n
          tmp = tmp + ( a(j,i) * g(j) ) / ( sx(j) * sx(j) )
        end do
        beta = beta + tmp * tmp
      end do

      dlt = alpha * sqrt ( alpha ) / beta
      dlt = min ( dlt, stepmx )

    end if

  end if
!
!  Find the new step by More-Hebdon algorithm.
!
  do

    call hookst ( nr, n, g, a, udiag, p, sx, rnwtln, dlt, amu, dltp, phi, &
      phip0, fstime, sc, nwtake, wrk0, epsm, ipr )

    dltp = dlt
!
!  Check the new point and update trust region.
!

    CALL TREGUP(NR, N, X, F, G, A, SC, SX, NWTAKE, STEPMX, STEPTL, DLT, IRETCD, &
                XPLSP, FPLSP, XPLS, FPLS, MXTAKE, IPR, 3, UDIAG, N_INT, N1, N2, &
                CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)

    if ( iretcd <= 1 ) then
      exit
    end if

  end do

  return
END SUBROUTINE HOOKDR






!*****************************************************************************80
!
!! HOOKST finds the new step by the More-Hebdon algorithm.
!
!  Modified:
!
!    15 May 2005
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) G(N), the gradient at the current iterate.
!
!    Input, real ( kind = 8 ) A(NR,N), an N by N array.  It contains the
!    Cholesky decomposition of the hessian in the lower triangular
!    part and diagonal; the hessian or approximation in the upper
!    triangular part.
!
!    Input, real ( kind = 8 ) UDIAG(N), the diagonal of the hessian. whose lower
!    triangular part is stored in A.
!
!    Input, real ( kind = 8 ) P(N), the Newton step.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, real ( kind = 8 ) RNWTLN, the Newton step length.
!
!    Input/output, real ( kind = 8 ) DLT, the trust region radius.
!
!    Workspace, real ( kind = 8 ) AMU, [retain value between successive calls]
!
!    Input, real ( kind = 8 ) DLTP, the trust region radius at last exit
!    from this routine.
!
!    Workspace, real ( kind = 8 ) PHI, [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) PHIP0, [retain value between successive calls]
!
!    Input/output, logical FSTIME, TRUE if first entry to this routine
!    during k-th iteration.
!
!    Output, real ( kind = 8 ) SC(N), the current step.
!
!    Output, logical NWTAKE, is TRUE if a Newton step taken.
!
!    Workspace, real ( kind = 8 ) WRK0(N).
!
!    Input, real ( kind = 8 ) EPSM, the machine epsilon.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
SUBROUTINE HOOKST( nr, n, g, a, udiag, p, sx, rnwtln, dlt, amu, &
  dltp, phi, phip0, fstime, sc, nwtake, wrk0, epsm, ipr )
IMPLICIT NONE

INTEGER, INTENT(IN) :: N, NR, IPR
!INTEGER, INTENT(INOUT) ::
DOUBLE PRECISION, INTENT(IN) :: G(N), UDIAG(N), P(N), SX(N)
DOUBLE PRECISION, INTENT(IN) :: RNWTLN, DLTP, EPSM
DOUBLE PRECISION, INTENT(INOUT) :: AMU, PHI, PHIP0, SC(N), WRK0(N), DLT
DOUBLE PRECISION, INTENT(INOUT) :: A(NR,N)
LOGICAL, INTENT(INOUT) :: FSTIME, NWTAKE
DOUBLE PRECISION :: ADDMAX, AMULO, AMUUP, PHIP, STEPLN
INTEGER :: I, J
DOUBLE PRECISION, PARAMETER :: hi = 1.50D+00, alo = 0.75D+00

!
!  Take a Newton step?
!
  if ( rnwtln <= hi * dlt ) then
    nwtake = .true.
    sc(1:n) = p(1:n)
    dlt = min ( dlt, rnwtln )
    amu = 0.0D+00
    return
  end if
!
!  Newton step not taken.
!
  nwtake = .false.

  if ( 0.0D+00 < amu ) then
    amu = amu - ( phi + dltp ) * ( ( dltp - dlt ) + phi ) / ( dlt * phip )
  end if

  phi = rnwtln - dlt

  if ( fstime ) then

    wrk0(1:n) = sx(1:n) * sx(1:n) * p(1:n)
!
!  Solve L * Y = (SX**2)*P
!
    call forslv ( nr, n, a, wrk0, wrk0 )

    phip0 = -dnrm2 ( n, wrk0, 1 )**2 / rnwtln
    fstime = .false.

  end if

  phip = phip0
  amulo = -phi / phip
  amuup = 0.0D+00
  do i = 1, n
    amuup = amuup + ( g(i) * g(i) ) / ( sx(i) * sx(i) )
  end do
  amuup = sqrt ( amuup ) / dlt
!
!  Test the value of amu; generate next amu if necessary.
!
  do

    if ( amu < amulo .or. amuup < amu ) then
      amu = max ( sqrt ( amulo * amuup ), amuup * 1.0D-03 )
    end if
!
!  Copy (h,udiag) to L
!  where h <-- h + amu*(sx**2) [do not actually change (h,udiag)]
!
    do j = 1, n
      a(j,j) = udiag(j) + amu * sx(j) * sx(j)
      a(j+1:n,j) = a(j,j+1:n)
    end do
!
!  Factor h=l(l+)
!
    call choldc ( nr, n, a, 0.0D+00, sqrt ( epsm ), addmax )
!
!  Solve h*p = l(l+) * sc = -g
!
    wrk0(1:n) = -g(1:n)

    call lltslv ( nr, n, a, sc, wrk0 )
!
!  Reset H.  Note since UDIAG has not been destroyed, we need do
!  nothing here.  H is in the upper part and in UDIAG, still intact
!
    stepln = sqrt ( dot_product ( sx(1:n)**2, sc(1:n)**2 ) )

    phi = stepln - dlt

    wrk0(1:n) = sx(1:n)**2 * sc(1:n)

    call forslv ( nr, n, a, wrk0, wrk0 )

    phip = -dnrm2 ( n, wrk0, 1 )**2 / stepln
!
!  If SC not acceptable hookstep, then select new AMU.
!
    if ( ( stepln < alo * dlt .or. hi * dlt < stepln ) .and. &
      ( 0.0D+00 < amuup - amulo ) ) then

      amulo = max ( amulo, amu - ( phi / phip ) )

      if ( phi < 0.0D+00 ) then
        amuup = min ( amuup, amu )
      end if

      amu = amu - ( stepln * phi ) / ( dlt * phip )
!
!  SC is acceptable hookstep.
!
    else

      exit

    end if

  end do

  return
END SUBROUTINE HOOKST






!*****************************************************************************80
!
!! HSNINT provides initial hessian when using secant updates.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Output, real ( kind = 8 ) A(NR,N), the initial N by N Hessian.  Only the
!    lower triangle of the matrix is assigned values.
!
!    Input, real ( kind = 8 ) SX(N), the scaling factors for X.
!
!    Input, integer ( kind = 4 ) METHOD, specifies the algorithm to use to solve
!    the minimization problem.
!    1 or 2: factored secant method used
!    3:  unfactored secant method used
!
SUBROUTINE HSNINT(nr, n, a, sx, method )
IMPLICIT NONE
INTEGER, INTENT(IN) :: NR, N, METHOD
DOUBLE PRECISION, INTENT(IN) ::SX(N)
DOUBLE PRECISION, INTENT(INOUT) :: A(NR,N)
INTEGER :: J
!  integer ( kind = 4 ) n
!  integer ( kind = 4 ) nr

!  real ( kind = 8 ) a(nr,n)
!  integer ( kind = 4 ) j
!  integer ( kind = 4 ) method
!  real ( kind = 8 ) sx(n)

DO j = 1, n

    IF ( method == 3 ) THEN
        A(j,j)=SX(j)**2
    ELSE
        A(j,j)=SX(j)
    END IF

    a(j+1:n,j) = 0.0D+00

END DO

RETURN
END SUBROUTINE HSNINT






!*****************************************************************************80
!
!! LLTSLV solves A*x=b where A = L * L'.
!
!  Discussion:
!
!    L is a lower triangular matrix.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) A(NR,N), contains the lower triangular matrix L.
!
!    Output, real X(N), the solution vector.
!
!    Input, real B(N), the right hand side vector.  If B is not required by
!    the calling program, then B and X may share the same storage.
!
SUBROUTINE LLTSLV( nr, n, a, x, b )
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, NR
DOUBLE PRECISION, INTENT(IN) :: A(NR,N), B(N)
DOUBLE PRECISION, INTENT(INOUT) :: X(N)

!
!  Forward solve, result in X.
!
  call forslv ( nr, n, a, x, b )
!
!  Back solve, result in X.
!
  call bakslv ( nr, n, a, x, x )

  RETURN
  END SUBROUTINE LLTSLV




!*****************************************************************************80
!
!! LNSRCH finds a next Newton iterate by line search.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the old iterate, sometimes called X[K-1].
!
!    Input, real ( kind = 8 ) F, the function value at the old iterate.
!
!    Input, real ( kind = 8 ) G(N), the gradient at the old iterate, or an
!    approximation to that value.
!
!    Input, real ( kind = 8 ) P(N), the (non-zero) Newton step.
!
!    Output, real ( kind = 8 ) XPLS(N), the new iterate.
!
!    Output, real ( kind = 8 ) FPLS, the function value at the new iterate.
!
!    Input, external COST, the name of the subroutine to evaluate the function,
!    of the form
!
!      subroutine COST ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) x(n)
!      real ( kind = 8 ) f
!
!    Output, integer ( kind = 4 ) IRETCD, the return code.
!
!    Output, logical MXTAKE, is TRUE if a step of maximum size was used.
!
!    Input, real ( kind = 8 ) STEPMX, the maximum allowable step size.
!
!    Input, real ( kind = 8 ) STEPTL, the relative step size at which successive
!    iterates are considered close enough to terminate algorithm.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
!  Local variables:
!
!    sln, the Newton length.
!
!    rln, the relative length of Newton step
!
SUBROUTINE lnsrch ( n, x, f, g, p, xpls, fpls, mxtake, iretcd, stepmx, &
  steptl, sx, ipr, N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON )
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, IPR, NCON
INTEGER, INTENT(INOUT) :: IRETCD
DOUBLE PRECISION, INTENT(IN) :: X(N), F, G(N), STEPMX, STEPTL, SX(N)
DOUBLE PRECISION, INTENT(INOUT) :: XPLS(N), FPLS, P(N), G_CON(NCON)
LOGICAL, INTENT(INOUT) :: MXTAKE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NEW INPUT VARIABLES
INTEGER, INTENT(IN) :: N_INT, N1, N2, CHROM_INT(N_INT), INFO
DOUBLE PRECISION, INTENT(IN) :: INPUT_ARRAY(N1,N2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION :: A, ALMBDA, B, DISC, PFPLS, PLMBDA, RLN, RMNLMB, SCL
DOUBLE PRECISION :: SLN, SLP, T1, T2, T3, TLMBDA
INTEGER :: I

  mxtake = .false.
  iretcd = 2

  sln = sqrt ( sum ( ( sx(1:n) * p(1:n) )**2 ) )
!
!  Newton step longer than maximum allowed.
!
  if ( stepmx < sln ) then
    scl = stepmx / sln
    p(1:n) = p(1:n) * stepmx / sln
    sln = stepmx
  end if

  slp = dot_product ( g, p )

  rln = 0.0D+00
  do i = 1, n
    rln = max ( rln, abs ( p(i) ) / max ( abs ( x(i) ), 1.0D+00 / sx(i) ) )
  end do

  rmnlmb = steptl / rln
  almbda = 1.0D+00
!
!  Check if new iterate satisfactory.  Generate new lambda if necessary.
!
  do

    if ( iretcd < 2 ) then
      exit
    end if

    xpls(1:n) = x(1:n) + almbda * p(1:n)

    !call COST ( n, xpls, fpls )
    CALL COST(N, N_INT, N1, N2, XPLS, CHROM_INT, FPLS, INPUT_ARRAY, G_CON, NCON)

    if ( f + slp * 0.0001D+00 * almbda < fpls ) then
      go to 130
    end if
!
!  Solution found.
!
    iretcd = 0

    if ( almbda == 1.0D+00 .and. 0.99D+00 * stepmx < sln ) then
      mxtake = .true.
    end if

    cycle
!
!  Solution not (yet) found.
!
130  continue
!
!  No satisfactory XPLS found sufficiently distinct from X.
!
    if ( almbda < rmnlmb ) then
      iretcd = 1
      cycle
    end if
!
!  Calculate new lambda.
!
!  First backtrack: quadratic fit.
!
    if ( almbda == 1.0D+00 ) then
      tlmbda = -slp / ( 2.0D+00 * ( fpls - f - slp ) )
      go to 170
    end if
!
!  All subsequent backtracks: cubic fit.
!
150 continue

    t1 = fpls - f - almbda * slp
    t2 = pfpls - f - plmbda * slp
    t3 = 1.0D+00 / ( almbda - plmbda )
    a = t3 * ( t1 / ( almbda * almbda ) - t2 / ( plmbda * plmbda ) )
    b = t3 * ( t2 *  almbda / ( plmbda * plmbda ) &
      - t1 * plmbda / ( almbda * almbda ) )
    disc = b * b - 3.0D+00 * a * slp

    if ( disc <= b * b ) then
      go to 160
    end if
!
!  Only one positive critical point, must be minimum.
!
    tlmbda = ( - b + sign ( 1.0D+00, a ) * sqrt ( disc ) ) / ( 3.0D+00 * a )
    go to 165
!
!  Both critical points positive, first is minimum.
!
160 continue

    tlmbda = ( -b - sign ( 1.0D+00, a ) * sqrt ( disc ) ) / ( 3.0D+00 * a )

165 continue

    if ( 0.5D+00 * almbda < tlmbda ) then
      tlmbda = 0.5D+00 * almbda
    end if

170 continue

    plmbda = almbda
    pfpls = fpls

    if ( almbda * 0.1D+00 <= tlmbda ) then
      almbda = tlmbda
    else
      almbda = almbda * 0.1D+00
    end if

  end do

RETURN
END SUBROUTINE LNSRCH






!*****************************************************************************80
!
!! MVMLTL computes y = L * x where L is a lower triangular matrix stored in A.
!
!  Discussion:
!
!    Note that X and Y cannot share storage.
!
!  Modified:
!
!    29 May 2001
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) A(NR,N), the N by N lower triangular matrix.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
!
!    Output, real ( kind = 8 ) Y(N), the result.
!

SUBROUTINE MVMLTL( nr, n, a, x, y )
IMPLICIT NONE
INTEGER, INTENT(IN) :: NR, N
DOUBLE PRECISION, INTENT(IN) :: A(NR,N), X(N)
DOUBLE PRECISION, INTENT(INOUT):: Y(N)
INTEGER :: I

DO i = 1, n
    Y(I)= dot_product(A(I,1:I), X(1:I))
END DO

RETURN
END SUBROUTINE MVMLTL



!*****************************************************************************80
!
!! MVMLTS computes y = A * x where A is a symmetric matrix.
!
!  Discussion:
!
!    A is a symmetric N by N matrix stored in its lower triangular part
!    and X and Y are N vectors.
!
!    X and Y cannot share storage.
!
!  Modified:
!
!    25 August 2001
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) A(NR,N), the symmetric N by N matrix.  The entries
!    of A are stored in the lower half of the array.
!
!    Input, real ( kind = 8 ) X(N), the vector to be multiplied by A.
!
!    Output, real ( kind = 8 ) Y(N), the result.
!
SUBROUTINE MVMLTS( nr, n, a, x, y)
IMPLICIT NONE
INTEGER, INTENT(IN) :: NR, N
DOUBLE PRECISION, INTENT(IN) ::A(NR,N), X(N)
DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
INTEGER :: I
!  integer ( kind = 4 ) n
!  integer ( kind = 4 ) nr

!  real ( kind = 8 ) a(nr,n)
!  integer ( kind = 4 ) i
!  real ( kind = 8 ) x(n)
!  real ( kind = 8 ) y(n)

DO i = 1, n

    Y(i) = dot_product(A(i,1:i), X(1:i) ) &
         + dot_product(A(i+1:n,i), X(i+1:n) )
END DO

RETURN
END SUBROUTINE MVMLTS






!*****************************************************************************80
!
!! MVMLTU computes y = L' * x where L is a lower triangular matrix.
!
!  Discussion:
!
!    Note that X and Y cannot share storage.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) A(NR,N), the N by N lower triangular matrix,
!
!    Input, real ( kind = 8 ) X(N), the matrix to be multiplied.
!
!    Output, real ( kind = 8 ) Y(N), the result vector.
!
SUBROUTINE MVMLTU( nr, n, a, x, y )
IMPLICIT NONE
INTEGER, INTENT(IN) :: NR, N
DOUBLE PRECISION, INTENT(IN) :: A(NR,N), X(N)
DOUBLE PRECISION, INTENT(INOUT) :: Y(N)
INTEGER :: I

DO i = 1, n
    y(i) = dot_product ( x(i:n), a(i:n,i) )
END DO

RETURN
END SUBROUTINE MVMLTU


!*****************************************************************************80
!
!! OPTCHK checks the input to the optimization routine.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), an approximate solution of the problem.
!
!    Input/output, real ( kind = 8 ) TYPSIZ(N), a typical size for each
!    component of X.  If TYPSIZ(I) is zero, it is reset to 1.
!
!    Input, real ( kind = 8 ) SX(N), the  diagonal scaling matrix for X.
!
!    Input/output, real ( kind = 8 ) FSCALE, an estimate of the scale of
!    the objective function COST.
!
!    Input, real ( kind = 8 ) GRADTL, the tolerance at which the gradient
!    is considered close enough to zero to terminate the algorithm.
!
!    Input/output, integer ( kind = 4 ) ITNLIM, the maximum number of allowable
!    iterations.
!
!    Input/output, integer ( kind = 4 ) NDIGIT, the number of good digits in
!    optimization function COST.
!
!    Input, real ( kind = 8 ) EPSM, the machine epsilon.
!
!    Input/output, real ( kind = 8 ) DLT, the trust region radius.
!
!    Input/output, integer ( kind = 4 ) METHOD, the algorithm indicator.
!
!    Input/output, integer ( kind = 4 ) IEXP, the expense flag.
!
!    Input/output, integer ( kind = 4 ) IAGFLG, = 1 if analytic gradient supplied.
!
!    Input/output, integer ( kind = 4 ) IAHFLG, = 1 if analytic hessian supplied.
!
!    Input/output, real ( kind = 8 ) STEPMX, the maximum step size.
!
!    Input/output, integer ( kind = 4 ) MSG, the message and error code.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
SUBROUTINE optchk ( n, x, typsiz, sx, fscale, gradtl, itnlim, ndigit, epsm, &
  dlt, method, iexp, iagflg, iahflg, stepmx, msg, ipr )
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, IPR
INTEGER, INTENT(INOUT) :: ITNLIM, NDIGIT, METHOD, IEXP, IAGFLG, IAHFLG, MSG
DOUBLE PRECISION, INTENT(IN) :: X(N), GRADTL, EPSM
DOUBLE PRECISION, INTENT(INOUT) :: TYPSIZ(N), FSCALE, DLT, STEPMX, SX(N)
DOUBLE PRECISION :: STPSIZ
INTEGER :: I

!
!  Check that parameters only take on acceptable values.
!  if not, set them to default values.
!
  if ( method < 1 .or. 3 < method ) then
    method = 1
  end if

  if ( iagflg /= 1 ) then
    iagflg = 0
  end if

  if ( iahflg /= 1 ) then
    iahflg = 0
  end if

  if ( iexp /= 0 ) then
    iexp = 1
  end if

  if ( mod ( msg/2, 2 ) == 1 .and. iagflg == 0 ) then
    !write ( ipr, 906 ) msg, iagflg
    msg = -6
    return
  end if

  if ( mod ( msg/4, 2 ) == 1 .and. iahflg == 0 ) then
    !write ( ipr, 907 ) msg, iahflg
    msg = -7
    return
  end if
!
!  Check N.
!
  if ( n <= 0 ) then
    !write ( ipr, * ) ' '
    !write ( ipr, * ) 'OPTCHK - Fatal error!'
    !write ( ipr, * ) '  Illegal nonpositive value of N = ', n
    msg = -1
    return
  end if

  if ( n == 1 .and. mod ( msg, 2 ) == 0 ) then
    !write ( ipr, 902 )
    msg = -2
    return
  end if
!
!  Compute the scale matrix.
!
  do i = 1, n
    if ( typsiz(i) == 0.0D+00 ) then
      typsiz(i) = 1.0D+00
    end if
  end do

  typsiz(1:n) = abs ( typsiz(1:n) )
  sx(1:n) = 1.0D+00 / typsiz(1:n)
!
!  Check maximum step size.
!
  if ( stepmx <= 0.0D+00 ) then

    stpsiz = sqrt ( sum ( x(1:n)**2 * sx(1:n)**2 ) )

    stepmx = max ( 1.0D+03 * stpsiz, 1.0D+03 )

  end if
!
!  Check the function scale.
!
  if ( fscale == 0.0D+00 ) then
    fscale = 1.0D+00
  end if

  if ( fscale < 0.0D+00 ) then
    fscale = -fscale
  end if
!
!  Check gradient tolerance
!
  if ( gradtl < 0.0D+00 ) then
    !write ( ipr, 903 ) gradtl
    msg = -3
    return
  end if
!
!  Check iteration limit
!
  if ( itnlim <= 0 ) then
    !write ( ipr, 904 ) itnlim
    msg = -4
    return
  end if
!
!  Check number of digits of accuracy in function COST.
!
  if ( ndigit == 0 ) then
    !write ( ipr, 905 ) ndigit
    msg = -5
    return
  end if

  if ( ndigit < 0 ) then
    ndigit = -log10 ( epsm )
  end if
!
!  Check trust region radius.
!
  if ( dlt <= 0.0D+00 ) then
    dlt = -1.0D+00
  end if

  if ( stepmx < dlt ) then
    dlt = stepmx
  end if

  902 format(' optchk    +++ warning +++  this package is inefficient', &
    'for problems of size n=1.'/ &
    ' optchk    check installation libraries for more appropriate routines.'/ &
    ' optchk    if none, set msg and resubmit.')
  903 format(' optchk    illegal tolerance.  gradtl=',e20.13)
  904 format(' optchk    illegal iteration limit.  itnlim=',i5)
  905 format(' optchk    minimization function has no good digits.', &
     'ndigit=',i5)
  906 format(' optchk    user requests that analytic gradient be', &
     ' accepted as properly coded (msg=',i5, '),'/ &
     ' optchk    but analytic gradient not supplied (iagflg=',i5, ').')
  907 format(' optchk    user requests that analytic hessian be', &
     ' accepted as properly coded (msg=',i5, '),'/ &
     ' optchk    but analytic hessian not supplied(iahflg=',i5, ').')

RETURN
END SUBROUTINE optchk


!*****************************************************************************80
!
!!  is a driver for the nonlinear optimization package.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, an rough solution of
!    the problem.  On output, the computed solution.
!
!    Input, external COST, the name of the subroutine that evaluates
!    the optimization function, of the form:
!      subroutine COST ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!
!    Input, external D1FCN, the name of the subroutine to evaluate gradient
!    of COST, of the form:
!      subroutine d1fcn ( n, x, g )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) g(n)
!      real ( kind = 8 ) x(n)
!
!    Input, external D2FCN, the name of the subroutine to evaluate the
!    Hessian of the function, of the form:
!      subroutine d2fcn ( nr, n, x, h )
!      integer ( kind = 4 ) nr
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) h(nr,n)
!      real ( kind = 8 ) x(n)
!
!    Input, real ( kind = 8 ) TYPSIZ(N), a typical size for each component of X.
!
!    Input, real ( kind = 8 ) FSCALE, an estimate of the scale of the objective
!    function.
!
!    Input, integer ( kind = 4 ) METHOD, the algorithm to use to solve
!    minimization problem:
!    1, line search
!    2, double dogleg
!    3, More-Hebdon
!
!    Input, integer ( kind = 4 ) IEXP, function expense flag.
!    Set IEXP to 1 if optimization function COST is expensive to
!    evaluate,  and 0 otherwise.  If set then hessian will
!    be evaluated by secant update instead of
!    analytically or by finite differences.
!
!    Input/output, integer ( kind = 4 ) MSG.
!    On input, set it positive to inhibit certain automatic checks
!    On output. < 0 indicates an error occurred.
!
!    Input, integer ( kind = 4 ) NDIGIT, the number of good digits in optimization
!    function COST.
!
!    Input, integer ( kind = 4 ) ITNLIM, the maximum number of allowable iterations.
!
!    Input, integer ( kind = 4 ) IAGFLG, is 1 if analytic gradient supplied.
!
!    Input, integer ( kind = 4 ) IAHFLG, is 1 if analytic hessian supplied.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
!    Input, real ( kind = 8 ) DLT, the trust region radius.
!
!    Input, real ( kind = 8 ) GRADTL, the tolerance at which gradient
!    considered close enough to zero to terminate algorithm.
!
!    Input, real ( kind = 8 ) STEPMX, the maximum allowable step size.
!
!    Input, real ( kind = 8 ) STEPTL, the relative step size at which
!    successive iterates are considered close enough to terminate algorithm
!
!    Input/output, real ( kind = 8 ) XPLS(N); on exit, XPLS is the local
!    minimizer.
!
!    Input/output, real ( kind = 8 ) FPLS; on exit, the function value at XPLS.
!
!    Input/output, real ( kind = 8 ) GPLS(N); on exit, the gradient at XPLS.
!
!    Output, integer ( kind = 4 ) ITRMCD, the termination code.
!
!    Workspace, real ( kind = 8 ) A(NR,N), workspace for hessian (or estimate)
!    and its Cholesky decomposition.
!
!    Workspace, real ( kind = 8 ) UDIAG(N), workspace for diagonal of hessian.
!
!    Workspace, real ( kind = 8 ) G(N), workspace for gradient at current
!    iterate.
!
!    Workspace, real ( kind = 8 ) P(N), workspace for the step.
!
!    Workspace, real ( kind = 8 ) SX(N), workspace for diagonal scaling matrix.
!
!    Workspace, real ( kind = 8 ) WRK0(N), WRK1(N), WRK2(N), WRK3(N).
!
!  Local variables:
!
!    analtl, tolerance for gradient and hessian checking.
!
!    epsm, machine epsilon.
!
!    f, function value: COST(x).
!
!    itncnt, current iteration, k
!
!    rnf, relative noise in optimization function COST.
!
!    noise=10.**(-ndigit)
!
SUBROUTINE  optdrv( nr, n, x, typsiz, fscale, method, &
    iexp, msg, ndigit, itnlim, iagflg, iahflg, ipr, dlt, gradtl, stepmx, &
    steptl, xpls, fpls, gpls, itrmcd, a, udiag, g, p, sx, wrk0, wrk1, wrk2, &
    wrk3, N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON )
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, NR, NCON
INTEGER, INTENT(INOUT) :: IAGFLG, IAHFLG
INTEGER, INTENT(INOUT) :: MSG, ITRMCD, ITNLIM, NDIGIT, METHOD, IPR, IEXP
DOUBLE PRECISION, INTENT(IN) :: GRADTL
DOUBLE PRECISION, INTENT(IN) :: STEPTL
DOUBLE PRECISION, INTENT(INOUT) :: FSCALE, DLT
DOUBLE PRECISION, INTENT(INOUT) :: WRK0(N), WRK1(N), WRK2(N), WRK3(N), G_CON(NCON)
DOUBLE PRECISION, INTENT(INOUT) :: XPLS(N), X(N), FPLS, GPLS(N), A(NR,N)
DOUBLE PRECISION, INTENT(INOUT) :: UDIAG(N), G(N), P(N), SX(N), TYPSIZ(N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ADDITIONAL VARIABLE CHANGES
INTEGER, INTENT(IN) :: N_INT, N1, N2, CHROM_INT(N_INT), INFO
DOUBLE PRECISION, INTENT(IN) :: INPUT_ARRAY(N1,N2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER :: ICSCMX, IRETCD, ITNCNT
DOUBLE PRECISION :: AMU, AMUSAV, ANALTL, DLPSAV, DLTP, DLTSAV, EPSM, F
DOUBLE PRECISION ::PHI, PHIP0, PHISAV, PHPSAV, RNF, VALUE(1), WRK(1), STEPMX
LOGICAL :: MXTAKE, NOUPDT

!
!  Initialization.
!
  p(1:n) = 0.0D+00
  itncnt = 0
  iretcd = -1
  epsm = epsilon ( epsm )

!WRITE(*,*) "CALL OPTCHK"
  call optchk ( n, x, typsiz, sx, fscale, gradtl, itnlim, ndigit, epsm, &
    dlt, method, iexp, iagflg, iahflg, stepmx, msg, ipr )

  if ( msg < 0 ) then
    return
  end if

  rnf = max ( 10.0D+00**(-ndigit), epsm )

  analtl = max ( 1.0D-02, sqrt ( rnf ) )

  if ( mod ( msg / 8, 2 ) == 0 ) then
    !write ( ipr, 901 )
    !write ( ipr, 900 ) typsiz(1:n)
    !write ( ipr, 902 )
    !write ( ipr, 900 ) sx(1:n)
    !write ( ipr, 903 ) fscale
    !write ( ipr, 904 ) ndigit,iagflg,iahflg,iexp,method,itnlim,epsm
    !write ( ipr, 905 ) stepmx,steptl,gradtl,dlt,rnf,analtl
  end if
!
!  Evaluate COST(x)
!
  !call COST ( n, x, f )
  !WRITE(*,*) "made it to COST"
  CALL COST(N, N_INT, N1, N2, X, CHROM_INT, F, INPUT_ARRAY, G_CON, NCON)
  !WRITE(*,*) "made it out of COST"
!
!  Evaluate analytic or finite difference gradient and check analytic
!  gradient, if requested.
!
  if ( iagflg /= 1 ) then

    value(1) = f

    CALL FSTOFD(1, 1, N, X, VALUE, G, SX, RNF, WRK, 1, &
                N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)

  else

    call d1fcn ( n, x, g )

    if ( mod ( msg/2, 2 ) /= 1 ) then

      call grdchk ( n, x, f, g, typsiz, sx, fscale, rnf, analtl, wrk1, &
        msg, ipr, N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON )

      if ( msg < 0 ) then
        return
      end if

    end if

  end if
  !WRITE(*,*) "made if to optstp"
  call optstp ( n, x, f, g, wrk1, itncnt, icscmx, itrmcd, gradtl, steptl, &
    sx, fscale, itnlim, iretcd, mxtake, ipr, msg )
  !WRITE(*,*) "MADE IT OUT OF OPTSTP"
  if ( itrmcd /= 0 ) then
    go to 700
  end if

  if ( iexp /= 1 ) then
    go to 80
  end if
!
!  If optimization function expensive to evaluate (iexp=1), then
!  hessian will be obtained by secant updates.  Get initial hessian.
!
  call hsnint ( nr, n, a, sx, method )
  go to 90

80 continue
!
!  Evaluate analytic or finite difference hessian and check analytic
!  hessian if requested (only if user-supplied analytic hessian
!  routine d2fcn fills only lower triangular part and diagonal of a).
!
  if ( iahflg == 1 ) then
    go to 82
  end if

  if ( iagflg == 1 ) then

    CALL FSTOFD(NR, N, N, X, G, A, SX, RNF, WRK1, 3, &
                N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)

  else

    CALL SNDOFD(NR, N, X, F, A, SX, RNF, WRK1, WRK2, N_INT, N1, N2, &
                CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
  end if

  go to 88

82 continue

  if ( mod ( msg / 4, 2 ) == 0 ) then
    go to 85
  end if

  call d2fcn ( nr, n, x, a )
  go to 88

85 continue

  !call heschk ( nr, n, x, f, g, a, typsiz, &
  !  sx, rnf, analtl, iagflg, udiag, wrk1, wrk2, msg, ipr )
  !WRITE(*,*) "CALL HESCHK"
  CALL HESCHK(NR, N, X, F, G, A, TYPSIZ, SX, RNF, ANALTL, IAGFLG, &
              UDIAG, WRK1, WRK2, MSG, IPR, N_INT, N1, N2, &
              CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
!
!  HESCHK evaluates d2fcn and checks it against the finite
!  difference hessian which it calculates by calling fstofd
!  (if iagflg == 1) or sndofd (otherwise).
!
  if ( msg < 0 ) then
    return
  end if

88 continue

90 continue

  if ( mod ( msg / 8, 2 ) == 0 ) then
    call result ( nr, n, x, f, g, a, p, itncnt, 1, ipr )
  end if
!
!  iteration
!
  100 continue

  itncnt = itncnt + 1
!
!  Find perturbed local model hessian and its ll+ decomposition
!  (skip this step if line search or dogstep techniques being used with
!  secant updates.  Cholesky decomposition l already obtained from
!  secfac.)
!
  if ( iexp == 1 .and. method /= 3 ) then
    go to 105
  end if

  103   continue
    !WRITE(*,*) "CALL CHLHSN"
  call chlhsn ( nr, n, a, epsm, sx, udiag )
  105 continue
!
!  Solve for Newton step:  ap = -g
!
  wrk1(1:n) = - g(1:n)
    !WRITE(*,*) "CALL LLTSLV"
  call lltslv ( nr, n, a, p, wrk1 )
!
!  Decide whether to accept Newton step  xpls = x + p
!  or to choose xpls by a global strategy.
!
  if ( iagflg == 0 .and. method /= 1 ) then

    dltsav = dlt

    if ( method /= 2 ) then
      amusav = amu
      dlpsav = dltp
      phisav = phi
      phpsav = phip0
    end if

  end if

  if ( method == 1 ) then

    !call lnsrch ( n, x, f, g, p, xpls, fpls, mxtake, iretcd, &
    !  stepmx, steptl, sx, ipr )
    !WRITE(*,*) "CALL LNSRCH"
    CALL LNSRCH(N, X, F, G, P, XPLS, FPLS, MXTAKE, IRETCD, STEPMX, STEPTL, &
                SX, IPR, N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)

  else if ( method == 2 ) then

    call dogdrv ( nr, n, x, f, g, a, p, xpls, fpls, sx, stepmx, &
      steptl, dlt, iretcd, mxtake, wrk0, wrk1, wrk2, wrk3, ipr, &
      N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON )

  else if ( method == 3 ) then

    call hookdr ( nr, n, x, f, g, a, udiag, p, xpls, fpls, sx, stepmx, &
      steptl, dlt, iretcd, mxtake, amu, dltp, phi, phip0, wrk0, &
      wrk1, wrk2, epsm, itncnt, ipr, N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)

  end if
!
!  If could not find satisfactory step and forward difference
!  gradient was used, retry using central difference gradient.
!
  if ( iretcd /= 1 .or. iagflg /= 0 ) then
    go to 112
  end if
!
!  Set iagflg for central differences.
!
     iagflg = -1
     !write(ipr,906) itncnt
     !call fstocd ( n, x, sx, rnf, g )
     CALL FSTOCD(N,X,SX,RNF, G,  N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
     if ( method == 1 ) then
       go to 105
     end if

     dlt = dltsav

     if ( method == 2 ) then
       go to 105
     end if

     amu = amusav
     dltp = dlpsav
     phi = phisav
     phip0 = phpsav
     go to 103
!
!  Calculate step for output
!
  112 continue

  p(1:n) = xpls(1:n) - x(1:n)
!
!  Calculate the gradient at XPLS.
!
  if ( iagflg == -1 ) then
    !call fstocd ( n, xpls, sx, rnf, gpls )
    CALL FSTOCD(N, XPLS, SX, RNF, GPLS, N_INT, N1, N2, &
                CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
  else if ( iagflg == 0 ) then
    value(1) = fpls

    CALL FSTOFD(1, 1, N, XPLS, VALUE, GPLS, SX, RNF, WRK, 1, &
                N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
  else
    call d1fcn ( n, xpls, gpls )
  end if
!
!  Check whether stopping criteria satisfied.
!
  call optstp ( n, xpls, fpls, gpls, x, itncnt, icscmx, itrmcd, gradtl, &
    steptl, sx, fscale, itnlim, iretcd, mxtake, ipr, msg )

  if ( itrmcd /= 0 ) then
    go to 690
  end if
!
!  Evaluate hessian at xpls
!
  if ( iexp == 0 ) then
    go to 130
  end if

  if ( method == 3 ) then
     call secunf ( nr, n, x, g, a, udiag, xpls, gpls, epsm, itncnt, rnf, &
       iagflg, noupdt, wrk1, wrk2, wrk3 )
  else
    call secfac ( nr, n, x, g, a, xpls, gpls, epsm, itncnt, rnf, iagflg, &
      noupdt, wrk0, wrk1, wrk2, wrk3 )
  end if

  go to 150

  130 continue

  if ( iahflg == 1 ) then
    go to 140
  end if

  if ( iagflg == 1 ) then
    !call fstofd ( nr, n, n, xpls, gpls, a, sx, rnf, wrk1, 3 )
    CALL FSTOFD(NR, N, N, XPLS, GPLS, A, SX, RNF, WRK1, 3, &
                N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
  else
    !call sndofd ( nr, n, xpls, fpls, a, sx, rnf, wrk1, wrk2 )
    CALL SNDOFD(NR, N, XPLS, FPLS, A, SX, RNF, WRK1, WRK2, N_INT, N1, N2, &
                CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
  end if

  go to 150

  140 continue

  call d2fcn ( nr, n, xpls, a )

  150 continue

  if ( mod ( msg / 16, 2 ) == 1 ) then
    call result ( nr, n, xpls, fpls, gpls, a, p, itncnt, 1, ipr )
  end if
!
!  x <-- xpls
!  g <-- gpls
!  f <-- fpls
!
  f = fpls
  x(1:n) = xpls(1:n)
  g(1:n) = gpls(1:n)

  go to 100
!
!  Termination.
!
!  Reset XPLS, FPLS, GPLS, if previous iterate solution
!
  690 if ( itrmcd /= 3 ) then
    go to 710
  end if

  700 continue

  fpls = f
  xpls(1:n) = x(1:n)
  gpls(1:n) = g(1:n)
!
!  Print results
!
  710 continue

  if ( mod ( msg / 8, 2 ) == 0) then
    call result ( nr, n, xpls, fpls, gpls, a, p, itncnt, 0, ipr )
  end if

  msg = 0

  900 format('        ',5(e20.13,3x))
  901 format('0    typical x')
  902 format('     diagonal scaling matrix for x')
  903 format('     typical f =',e20.13)
  904 format('0    number of good digits in COST=',i5/ &
             '     gradient flag  =',i5,'   (=1 if analytic', &
             ' gradient supplied)'/ &
             '     hessian flag   =',i5,'   (=1 if analytic', &
             ' hessian supplied)'/ &
             '     expense flag   =',i5, '   (=1 if', &
             ' minimization function expensive to evaluate)'/ &
             '     method to use  =',i5,'   (=1,2,3 for line', &
             ' search, double dogleg, more-hebdon respectively)'/ &
             '     iteration limit=',i5/ &
             '     machine epsilon=',e20.13)

  905 format('0    maximum step size =',e20.13/ &
             '     step tolerance    =',e20.13/ &
             '     gradient tolerance=',e20.13/ &
             '     trust reg radius  =',e20.13/ &
             '     rel noise in COST  =',e20.13/ &
             '     anal-fd tolerance =',e20.13)

  906 format('     shift from forward to central differences', &
     ' in iteration ', i5)

  return
END SUBROUTINE



!*****************************************************************************80
!
!! OPTIF0 provides a simple interface to the minimization package.
!
!  Modified:
!
!    29 May 2001
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, an rough solution of
!    the problem.  On output, the computed solution.
!
!    Input, external COST, the name of the subroutine that evaluates
!    the optimization function, of the form
!
!      subroutine COST ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x(n)
!
!    Output, real ( kind = 8 ) XPLS(N), estimated local minimizer of
!    the function.
!
!    Output, real ( kind = 8 ) FPLS, the function value at XPLS.
!
!    Output, real ( kind = 8 ) GPLS(N), the gradient at XPLS.
!
!    Output, integer ( kind = 4 ) ITRMCD, the termination code.
!    1, relative gradient close to zero.
!       The current iterate is probably solution.
!    2, successive iterates within tolerance.
!       The current iterate is probably solution.
!    3, the last global step failed to locate a point lower than X.
!       Either x is an approximate local minimum of the function,
!       the function is too non-linear for this algorithm,
!       or STEPTL is too large.
!    4, iteration limit exceeded.  The algorithm failed.
!    5, maximum step size exceeded 5 consecutive times.
!       Either the function is unbounded below, becomes asymptotic to a
!       finite value from above in some direction, or STEPMX is too small.
!
SUBROUTINE OPTIF0 (n, x, xpls, fpls, gpls, itrmcd, N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, NCON
INTEGER, INTENT(INOUT) :: ITRMCD
DOUBLE PRECISION, INTENT(INOUT) :: X(N), XPLS(N), FPLS, GPLS(N), G_CON(NCON)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER,INTENT(IN) :: N_INT, N1, N2, CHROM_INT(N_INT), INFO
DOUBLE PRECISION, INTENT(IN) :: INPUT_ARRAY(N1,N2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION :: A(N,N), DLT, FSCALE, GRADTL, STEPMX, STEPTL, WRK(N,9)
INTEGER :: IAGFLG, IAHFLG, IEXP, IPR, ITNLIM, METHOD, MSG, NDIGIT, NR

! equivalence wrk(n,1) = udiag(n)
!             wrk(n,2) = g(n)
!             wrk(n,3) = p(n)
!             wrk(n,4) = typsiz(n)
!             wrk(n,5) = sx(n)
!             wrk(n,6) = wrk0(n)
!             wrk(n,7) = wrk1(n)
!             wrk(n,8) = wrk2(n)
!             wrk(n,9) = wrk3(n)
!
  nr = n

  call dfault ( n, x, wrk(1,4), fscale, method, iexp, msg, ndigit, &
    itnlim, iagflg, iahflg, ipr, dlt, gradtl, stepmx, steptl )

  call  OPTDRV( nr, n, x, wrk(1,4), fscale, &
    method, iexp, msg, ndigit, itnlim, iagflg, iahflg, ipr, &
    dlt, gradtl, stepmx, steptl, xpls, fpls, gpls, itrmcd, &
    a, wrk(1,1), wrk(1,2), wrk(1,3), wrk(1,5), wrk(1,6), &
    wrk(1,7), wrk(1,8), wrk(1,9) , N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON)

RETURN
END SUBROUTINE OPTIF0


!*****************************************************************************80
!
!! OPTSTP: unconstrained minimization stopping criteria
!
!  Discussion:
!
!    OPSTP determines whether the optimization algorithm should terminate,
!    due to any of the following:
!    1) the problem has been solved to the user's tolerance;
!    2) convergence within user tolerance;
!    3) iteration limit reached;
!    4) divergence or too restrictive maximum step (stepmx) suspected;
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) XPLS(N), the new iterate, X[K].
!
!    Input, real ( kind = 8 ) FPLS, the function value at the new iterate.
!
!    Input, real ( kind = 8 ) GPLS(N), the gradient at the new iterate, or an
!    approximation of that value.
!
!    Input, real ( kind = 8 ) X(N), the old iterate X[K-1].
!
!    Input, integer ( kind = 4 ) ITNCNT, the current iteration K.
!
!    Input/output, integer ( kind = 4 ) ICSCMX, the number of consecutive steps
!    greater than or equal to STEPMX.
!    [retain value between successive calls].
!
!    Output, integer ( kind = 4 ) ITRMD, the termination code.
!
!    Input, real ( kind = 8 ) GRADTL, the tolerance at which relative gradient
!    considered close enough to zero to terminate algorithm.
!
!    Input, real ( kind = 8 ) STEPTL, the relative step size at which
!    successive iterates are considered close enough to terminate algorithm.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, real ( kind = 8 ) FSCALE, an estimate of the scale of objective
!    function.
!
!    Input, integer ( kind = 4 ) ITNLIM, the maximum number of allowable iterations.
!
!    Output, integer ( kind = 4 ) IRETCD, the return code.
!
!    Input, logical MXTAKE, TRUE if a step of maximum length was used.
!
!    Output, integer ( kind = 4 ) IPR, the device to which to send output.
!
!    Input, integer ( kind = 4 ) MSG, if includes a term 8, suppress output.
!
SUBROUTINE optstp (n, xpls, fpls, gpls, x, itncnt, icscmx, &
  itrmcd, gradtl, steptl, sx, fscale, itnlim, iretcd, mxtake, ipr, msg )
IMPLICIT NONE

INTEGER, INTENT(IN) :: N, ITNCNT, ITNLIM, MSG
INTEGER, INTENT(INOUT) :: ICSCMX, ITRMCD, IRETCD, IPR
DOUBLE PRECISION, INTENT(IN) :: XPLS(N), FPLS, GPLS(N), X(N), GRADTL, STEPTL
DOUBLE PRECISION, INTENT(IN) :: SX(N), FSCALE
LOGICAL, INTENT(IN) :: MXTAKE
DOUBLE PRECISION :: D, RELGRD, RELSTP, RGX, RSX
INTEGER :: I, JTRMCD

itrmcd = 0

!
!  Last global step failed to locate a point lower than X.
!
  if ( iretcd == 1 ) then
    jtrmcd = 3
    go to 600
  end if
!
!  Find direction in which relative gradient maximum.
!  Check whether within tolerance
!
  d = max ( abs ( fpls ), fscale )

  rgx = 0.0D+00
  do i = 1, n
    relgrd = abs ( gpls(i) ) * max ( abs ( xpls(i) ), 1.0D+00 / sx(i) ) / d
    rgx = max ( rgx, relgrd )
  end do

  jtrmcd = 1
  if ( rgx <= gradtl ) then
    go to 600
  end if

  if ( itncnt == 0 ) then
    return
  end if
!
!  Find direction in which relative stepsize is maximum.
!  Check whether within tolerance.
!
  rsx = 0.0D+00
  do i = 1, n
    relstp = abs ( xpls(i) - x(i) ) / max ( abs ( xpls(i) ), 1.0D+00 / sx(i) )
    rsx = max ( rsx, relstp )
  end do

  jtrmcd = 2
  if ( rsx <= steptl ) then
    go to 600
  end if
!
!  Check iteration limit.
!
  jtrmcd = 4
  if ( itnlim <= itncnt ) then
    go to 600
  end if
!
!  Check number of consecutive steps \ stepmx
!
  if ( .not. mxtake ) then
    icscmx = 0
    return
  else
    if ( mod ( msg / 8, 2 ) == 0 ) then
      !write ( ipr, 900 )
    end if
    icscmx = icscmx + 1
    if ( icscmx < 5 ) then
      return
    end if
    jtrmcd = 5
  end if
!
!  Print termination code
!
  600 continue

  itrmcd = jtrmcd

  if ( itrmcd == 1 ) then
    !write ( ipr, 901 )
  else if ( itrmcd == 2 ) then
    !write(ipr,902)
  else if ( itrmcd == 3 ) then
    !write(ipr,903)
  else if ( itrmcd == 4 ) then
    !write(ipr,904)
  else if ( itrmcd == 5 ) then
    !write(ipr,905)
  end if

  !900 format('0optstp    step of maximum length (stepmx) taken')
  !901 format('0optstp    relative gradient close to zero.'/ &
  !           ' optstp    current iterate is probably solution.')
  !902 format('0optstp    successive iterates within tolerance.'/ &
  !           ' optstp    current iterate is probably solution')
  !903 format('0optstp    last global step failed to locate a point', &
  !           ' lower than x.'/ &
  !           ' optstp    either x is an approximate local minimum', &
  !           ' of the function',/ &
  !           ' optstp    the function is too non-linear for this algorithm,'/ &
  !           ' optstp    or steptl is too large.')
  !904 format('optstp    iteration limit exceeded.'/'optstp    algorithm failed.')
  !905 format('0optstp    maximum step size exceeded 5 consecutive times.'/ &
  !           ' optstp    either the function is unbounded below',/ &
  !           ' optstp    becomes asymptotic to a finite value', &
  !           ' from above in some direction',/ &
  !           ' optstp    or stepmx is too small')

RETURN
END SUBROUTINE OPTSTP


!*****************************************************************************80
!
!! QRAUX1 interchanges two rows of an upper Hessenberg matrix.
!
!  Discussion:
!
!    QRAUX1 interchanges rows I and I+1 of the upper Hessenberg matrix
!    R, columns I to N.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the leading dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the matrix.
!
!    Input/output, real ( kind = 8 ) R(NR,N), the N by N upper Hessenberg
!    matrix.
!
!    Input, integer ( kind = 4 ) I, the index of the first row to interchange.
!
SUBROUTINE QRAUX1( nr, n, r, i )
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, NR, I
DOUBLE PRECISION, INTENT(INOUT) :: R(NR,N)
INTEGER :: J
!  integer ( kind = 4 ) n
!  integer ( kind = 4 ) nr

!  integer ( kind = 4 ) i
!  integer ( kind = 4 ) j
!  real ( kind = 8 ) r(nr,n)

DO j = i, n
    CALL R8_SWAP (R(I,J), R(I+1,J))
END DO

RETURN
END SUBROUTINE QRAUX1


!*****************************************************************************80
!
!! QRAUX2 pre-multiplies an upper Hessenberg matrix by a Jacobi rotation.
!
!  Discussion:
!
!    QRAUX2 pre-multiplies an upper Hessenberg matrix by a Jacobi rotation
!    J(I,I+1,A,B)
!
!  Modified:
!
!    15 December 2001
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) R(NR,N), the N by N upper Hessenberg
!    matrix.
!
!    Input, integer ( kind = 4 ) I, the index of the row.
!
!    Input, real ( kind = 8 ) A, B, scalars that define the rotation.
!
SUBROUTINE QRAUX2(nr, n, r, i, a, b)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, NR, I
DOUBLE PRECISION, INTENT(IN) :: A, B
DOUBLE PRECISION, INTENT(INOUT) :: R(NR,N)
DOUBLE PRECISION :: C, DEN, S, Y, Z
INTEGER :: J

DEN=sqrt(A*A+B*B)
C=A/DEN
S=B/DEN

DO J = I, N
    Y=R(I,J)
    Z=R(I+1,J)
    R(I,J)=C*Y-S*Z
    R(I+1,J)=S*Y+C*Z
END DO

RETURN
END SUBROUTINE QRAUX2





!*****************************************************************************80
!
!! QRUPDT updates a QR factorization.
!
!  Discussion:
!
!    The routine finds an orthogonal N by N matrix Q* and an upper triangular
!    N by N matrix R* such that (Q*)(R*) = R + U*V'
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input/output, real ( kind = 8 ) A(NR,N), on input, contains the original QR
!    factorization.  On output, contains the revised factorization.
!
!    Input, real ( kind = 8 ) U(N), V(N), vectors that describe the rank
!    one update applied to the original matrix A.
!
SUBROUTINE QRUPDT(NR, N, A, U, V)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, NR
DOUBLE PRECISION, INTENT(IN) :: V(N)
DOUBLE PRECISION, INTENT(INOUT) :: A(NR,N), U(N)
DOUBLE PRECISION :: T1, T2
INTEGER :: I, K

!
!  Determine the last non-zero in U.
!
  k = n

  do while ( u(k) == 0.0D+00 .and. 1 < k )
    k = k - 1
  end do
!
!  (k-1) Jacobi rotations transform
!    r + u(v+) --> (r*) + ( u(1) * e1 ) (v+)
!  which is upper Hessenberg
!
  if ( 1 < k ) then

    do i = k-1, 1, -1

      if ( u(i) == 0.0D+00 ) then
        call qraux1 ( nr, n, a, i )
        u(i) = u(i+1)
      else
        call qraux2 ( nr, n, a, i, u(i), -u(i+1) )
        u(i) = sqrt ( u(i) * u(i) + u(i+1) * u(i+1) )
      end if

    end do

  end if
!
!  R <-- R + ( u(1) * e1 ) (v+)
!
  a(1,1:n) = a(1,1:n) + u(1) * v(1:n)
!
!  (k-1) Jacobi rotations transform upper Hessenberg R
!  to upper triangular (R*)
!
    do i = 1, k-1

      if ( a(i,i) == 0.0D+00 ) then
        call qraux1 ( nr, n, a, i )
      else
        t1 = a(i,i)
        t2 = -a(i+1,i)
        call qraux2 ( nr, n, a, i, t1, t2 )
      end if

    end do

RETURN
END SUBROUTINE QRUPDT




!*****************************************************************************80
!
!! R8_SWAP swaps two R8's.
!
!  Modified:
!
!    22 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
!    Y have been interchanged.
!
SUBROUTINE r8_swap ( x, y )
IMPLICIT NONE
DOUBLE PRECISION, INTENT(INOUT) :: X,Y
DOUBLE PRECISION :: Z

Z=X
X=Y
Y=Z

RETURN
END SUBROUTINE R8_SWAP



!*****************************************************************************80
!
!! RESULT prints information about the optimization process.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the current iterate.
!
!    Input, real ( kind = 8 ) F, the function value at X.
!
!    Input, real ( kind = 8 ) G(N), the gradient at X.
!
!    Input, real ( kind = 8 ) A(NR,N), the N by N Hessian matrix at X.
!
!    Input, real ( kind = 8 ) P(N), the step taken.
!
!    Input, integer ( kind = 4 ) ITNCNT, the iteration number.
!
!    Input, integer ( kind = 4 ) IFLG, the flag controlling the amount of printout.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
SUBROUTINE RESULT( nr, n, x, f, g, a, p, itncnt, iflg, ipr )
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, NR, ITNCNT, IFLG, IPR
DOUBLE PRECISION, INTENT(IN) :: X(N), F, G(N), A(NR,N), P(N)
INTEGER :: I

  !write ( ipr, 903 ) itncnt

  if ( iflg /= 0 ) then
    !write ( ipr, * ) ' result       step'
    !write ( ipr,905) p(1:n)
  end if

  !write ( ipr, * ) ' result       x(k)'
  !write ( ipr, 905) x(1:n)
  !write ( ipr, * ) ' result     function at x(k)'
  !write ( ipr, 905) f
  !write ( ipr, * ) ' result       gradient at x(k)'
  !write ( ipr, 905) g(1:n)

  if ( iflg /= 0 ) then

    !write ( ipr, * ) ' result       Hessian at x(k)'
    do i = 1, n
      !write ( ipr, 900) i
      !write ( ipr, 902) a(i,1:i)
    end do

  end if

RETURN

  900 format(' result     row',i5)
  902 format(' result       ',5(2x,e20.13))
  903 format(/'0result    iterate k=',i5)
  905 format(' result               ',5(2x,e20.13) )
END SUBROUTINE RESULT



!*****************************************************************************80
!
!! SECFAC updates the hessian by the BFGS factored method.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the old iterate, X[K-1].
!
!    Input, real ( kind = 8 ) G(N), the gradient or an approximation,
!    at the old iterate.
!
!    Input/output, real ( kind = 8 ) A(NR,N).
!    On input, the Cholesky decomposition of hessian in lower part and diagonal.
!    On output, the updated Cholesky decomposition of hessian
!    in lower triangular part and diagonal
!
!    Input, real ( kind = 8 ) XPLS(N), the new iterate, X[K].
!
!    Input, real ( kind = 8 ) GPLSN(N), gradient, or an approximation,
!    at the new iterate.
!
!    Input, real ( kind = 8 ) EPSM, the machine epsilon.
!
!    Input, integer ( kind = 4 ) ITNCNT, the iteration count.
!
!    Input, real ( kind = 8 ) RNF, the relative noise in the optimization
!    function COST.
!
!    Input, integer ( kind = 4 ) IAGFLG, 1 if analytic gradient supplied.
!
!    Input/output, logical NOUPDT, is TRUE if there has been no update
!    yet.  The user should retain the output value between successive
!    calls.
!
!    Workspace, real ( kind = 8 ) S(N).
!
!    Workspace, real ( kind = 8 ) Y(N).
!
!    Workspace, real ( kind = 8 ) U(N).
!
!    Workspace, real ( kind = 8 ) W(N).
!
SUBROUTINE SECFAC(nr, n, x, g, a, xpls, gpls, epsm, itncnt, rnf, &
                  iagflg, noupdt, s, y, u, w )
IMPLICIT NONE
INTEGER, INTENT(IN) :: NR, N, ITNCNT, IAGFLG
DOUBLE PRECISION, INTENT(IN) :: X(N), G(N), XPLS(N), GPLS(N), EPSM, RNF
DOUBLE PRECISION, INTENT(INOUT) :: A(NR,N), S(N), Y(N), U(N), W(N)
LOGICAL, INTENT(INOUT) :: NOUPDT
DOUBLE PRECISION :: ALP, DEN1, DEN2, RELTOL, SNORM2, YNRM2
INTEGER :: I, J
LOGICAL :: SKPUPD

  if ( itncnt == 1 ) then
    noupdt = .true.
  end if

  s(1:n) = xpls(1:n) - x(1:n)
  y(1:n) = gpls(1:n) - g(1:n)

  den1 = dot_product ( s, y )

  snorm2 = dnrm2 ( n, s, 1)

  ynrm2 = dnrm2 ( n, y, 1)

  if ( den1 < sqrt ( epsm ) * snorm2 * ynrm2 ) then
    return
  end if

  call mvmltu ( nr, n, a, s, u )

  den2 = dot_product ( u, u )
!
!  L <-- sqrt ( den1 / den2 ) * L
!
  alp = sqrt ( den1 / den2 )

  if ( noupdt ) then

    u(1:n) = alp * u(1:n)

    do j = 1, n
      do i = j, n
        a(i,j) = alp * a(i,j)
      end do
    end do

    noupdt = .false.
    den2 = den1
    alp = 1.0D+00

  end if

  skpupd = .true.
!
!  W = l(l+)s = hs
!
  call mvmltl ( nr, n, a, u, w )
  i = 1

  if ( iagflg == 0 ) then
    reltol = sqrt ( rnf )
  else
    reltol = rnf
  end if

60  continue

  if ( i <= n .and. skpupd ) then

    if ( abs ( y(i) - w(i) ) < reltol * &
      max ( abs ( g(i) ), abs ( gpls(i) ) ) ) then
      i = i + 1
    else
      skpupd = .false.
    end if
    go to 60
  end if

  if ( skpupd ) then
    return
  end if
!
!  W = y-alp*l(l+)s
!
  w(1:n) = y(1:n) - alp * w(1:n)
!
!  ALP = 1 / sqrt ( den1 * den2 )
!
  alp = alp / den1
!
!  U = (l+) / sqrt ( den1 * den2 ) = (l+)s/ sqrt ( ( y+ ) s * (s+) l (l+) s )
!
  u(1:n) = alp * u(1:n)
!
!  Copy L into upper triangular part.  Zero L.
!
  do i = 2, n
    do j = 1, i-1
      a(j,i) = a(i,j)
      a(i,j) = 0.0D+00
    end do
  end do
!
!  Find Q, (l+) such that  q(l+) = (l+) + u(w+)
!
  call qrupdt ( nr, n, a, u, w )
!
!  Upper triangular part and diagonal of a now contain updated
!  Cholesky decomposition of hessian.  Copy back to lower triangular part.
!
  do i = 2, n
    do j = 1, i-1
      a(i,j) = a(j,i)
    end do
  end do

RETURN
END SUBROUTINE SECFAC




!*****************************************************************************80
!
!! SECUNF updates a Hessian matrix by the BFGS unfactored method.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the old iterate, X[K-1].
!
!    Input, real ( kind = 8 ) G(N), the gradient, or an approximate value,
!    at the  old iterate.
!
!    Input/output, real ( kind = 8 ) A(NR,N).
!    on entry: approximate hessian at old iterate
!    in upper triangular part (and udiag)
!    on exit:  updated approx hessian at new iterate
!    in lower triangular part and diagonal
!    [lower triangular part of symmetric matrix]
!
!    Input, real ( kind = 8 ) UDIAG(N), the diagonal entries of the hessian.
!
!    Input, real ( kind = 8 ) XPLS(N), the new iterate, X[K].
!
!    Input, real ( kind = 8 ) GPLS(N), the gradient or an approximate value, at
!    the new iterate
!
!    Input, real ( kind = 8 ) EPSM, the machine epsilon.
!
!    Input, integer ( kind = 4 ) ITNCNT, the iteration count.
!
!    Input, real ( kind = 8 ) RNF, the relative noise in the optimization
!    function.
!
!    Input, integer ( kind = 4 ) IAGFLG, =1 if analytic gradient supplied, =0 otherwise
!
!    Input/output, logical NOUPDT, TRUE if no update yet.
!    [retain value between successive calls]
!
!    Workspace, real ( kind = 8 ) S(N).
!
!    Workspace, real ( kind = 8 ) Y(N).
!
!    Workspace, real ( kind = 8 ) T(N).
!
SUBROUTINE SECUNF( nr, n, x, g, a, udiag, xpls, gpls, epsm, itncnt, &
                   rnf, iagflg, noupdt, s, y, t )
IMPLICIT NONE
INTEGER :: N, NR, ITNCNT, IAGFLG
DOUBLE PRECISION, INTENT(IN) :: X(N), G(N), UDIAG(N), XPLS(N), GPLS(N), EPSM
DOUBLE PRECISION, INTENT(IN) :: RNF
DOUBLE PRECISION, INTENT(INOUT) :: A(NR,N), S(N), Y(N), T(N)
LOGICAL, INTENT(INOUT) :: NOUPDT
DOUBLE PRECISION :: DEN1, DEN2, GAM, SNORM2, TOL, YNRM2
INTEGER :: I, J
LOGICAL :: SKPUPD
!
!  Copy hessian in upper triangular part and UDIAG to
!  lower triangular part and diagonal.
!
  do j = 1, n
    a(j,j) = udiag(j)
    do i = j+1, n
      a(i,j) = a(j,i)
    end do
  end do

  if ( itncnt == 1 ) then
    noupdt = .true.
  end if

  s(1:n) = xpls(1:n) - x(1:n)
  y(1:n) = gpls(1:n) - g(1:n)

  den1 = dot_product ( s, y )

  snorm2 = dnrm2 ( n, s, 1 )

  ynrm2 = dnrm2 ( n, y, 1 )

  if ( den1 < sqrt ( epsm ) * snorm2 * ynrm2 ) then
    return
  end if

  call mvmlts ( nr, n, a, s, t )

  den2 = dot_product ( s, t )

  if ( noupdt ) then
!
!  H <-- [(s+)y/(s+)hs]h
!
    gam = den1 / den2
    den2 = gam * den2
    do j = 1, n
      t(j) = gam * t(j)
      do i = j, n
        a(i,j) = gam * a(i,j)
      end do
    end do
    noupdt = .false.

  end if

  skpupd = .true.
!
!  Check update condition on row I.
!
  do i = 1, n

    tol = rnf * max ( abs ( g(i) ), abs ( gpls(i) ) )
    if ( iagflg == 0 ) then
      tol = tol / sqrt ( rnf )
    end if

    if ( tol <= abs ( y(i) - t(i) ) ) then
      skpupd = .false.
      exit
    end if

  end do

  if ( skpupd ) then
    return
  end if
!
!  BFGS update
!
  do j = 1, n
    do i = j, n
      a(i,j) = a(i,j) + y(i) * y(j) / den1 - t(i) * t(j) / den2
    end do
  end do

RETURN
END SUBROUTINE SECUNF



!*****************************************************************************80
!
!! SNDOFD approximates a Hessian with a second order finite difference.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) XPLS(N), the new iterate, X[K].
!
!    Input, external COST, the name of the subroutine to evaluate the function,
!    of the form
!
!      subroutine COST ( n, x, f )
!      integer ( kind = 4 ) n
!      real x(n)
!      real f
!
!    Input, real ( kind = 8 ) FPLS, the function value at the new iterate.
!
!    Output, real ( kind = 8 ) A(NR,N), the N by N  finite difference
!    approximation to the hessian matrix.  Only the lower triangular matrix and
!    diagonal are returned.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling for X.
!
!    Input, real ( kind = 8 ) RNOISE, the relative noise in the function.
!
!    Workspace, real ( kind = 8 ) STEPSZ(N), used for the stepsize.
!
!    Workspace, real ( kind = 8 ) ANBR(N), holds neighbors.
!
SUBROUTINE sndofd ( nr, n, xpls, fpls, a, sx, rnoise, stepsz, anbr, N_INT, N1, N2, &
                    CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON )
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, NR, NCON
DOUBLE PRECISION, INTENT(IN) :: FPLS, SX(N), RNOISE
DOUBLE PRECISION, INTENT(INOUT) :: A(NR,N),  STEPSZ(N), ANBR(N), XPLS(N), G_CON(NCON)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NEW INPUT VARIABLES
INTEGER, INTENT(IN) :: N_INT, N1, N2, CHROM_INT(N_INT), INFO
DOUBLE PRECISION, INTENT(IN) :: INPUT_ARRAY(N1,N2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION :: FHAT, OV3, XTMPI, XTMPJ
INTEGER :: I, J

!
!  Find I-th stepsize and evaluate neighbor in direction of I-th unit vector.
!
  ov3 = 1.0D+00 / 3.0D+00

  do i = 1, n
    stepsz(i) = rnoise**ov3 * max ( abs ( xpls(i) ), 1.0D+00 / sx(i) )
    xtmpi = xpls(i)
    xpls(i) = xtmpi + stepsz(i)
    !call COST ( n, xpls, anbr(i) )
    CALL COST(N, N_INT, N1, N2, XPLS, CHROM_INT, ANBR(I), INPUT_ARRAY, G_CON, NCON)
    xpls(i) = xtmpi
  end do
!
!  Calculate column I of A.
!
  do i = 1, n

    xtmpi = xpls(i)
    xpls(i) = xtmpi + 2.0D+00 * stepsz(i)
    !call COST ( n, xpls, fhat )
    CALL COST(N, N_INT, N1, N2, XPLS, CHROM_INT, FHAT, INPUT_ARRAY, G_CON, NCON)
    a(i,i) = ( ( fpls - anbr(i) ) &
      + ( fhat - anbr(i) ) ) / ( stepsz(i) * stepsz(i) )
!
!  Calculate sub-diagonal elements of column.
!
    if ( i /= n ) then

      xpls(i) = xtmpi + stepsz(i)

      do j = i + 1, n
        xtmpj = xpls(j)
        xpls(j) = xtmpj + stepsz(j)
        !call COST ( n, xpls, fhat )
        CALL COST(N, N_INT, N1, N2, XPLS, CHROM_INT, FHAT, INPUT_ARRAY, G_CON, NCON)
        a(j,i) = ( ( fpls - anbr(i) ) + ( fhat - anbr(j) ) ) &
          / ( stepsz(i) * stepsz(j) )
        xpls(j) = xtmpj
      end do

    end if

    xpls(i) = xtmpi

  end do

RETURN
END SUBROUTINE SNDOFD


!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
SUBROUTINE timestamp ( )
IMPLICIT NONE
INTEGER :: D, H, M, MM, N, S, VALUES(8), Y
CHARACTER(LEN=8) :: AMPM, DATE
CHARACTER(LEN=10) :: TIME
CHARACTER(LEN=5) :: ZONE
CHARACTER(LEN=9), PARAMETER, DIMENSION(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)

call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

RETURN
END SUBROUTINE TIMESTAMP






!*****************************************************************************80
!
!! TREGUP decides whether to accept the next optimization iterate.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NR, the row dimension of the matrix.
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X(N), the old iterate X[K-1].
!
!    Input, real ( kind = 8 ) F, the function value at the old iterate.
!
!    Input, real ( kind = 8 ) G(N), the gradient at the old iterate, or
!    an approximate value of it.
!
!    Input, real ( kind = 8 ) A(NR,N), the Cholesky decomposition of hessian in
!    lower triangular part and diagonal.  Hessian or approximation in
!    upper triangular part.
!
!    Input, external COST, the name of the subroutine to evaluate the function,
!    of the form
!
!      subroutine COST ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) x(n)
!      real ( kind = 8 ) f
!
!    Input, real ( kind = 8 ) SC(N), the current step.
!
!    Input, real ( kind = 8 ) SX(N), the diagonal scaling matrix for X.
!
!    Input, logical NWTAKE, is TRUE if a Newton step is taken.
!
!    Input, real ( kind = 8 ) STEPMX, the maximum allowable step size.
!
!    Input, real ( kind = 8 ) STEPTL, the relative step size at which successive
!    iterates are considered close enough to terminate algorithm.
!
!    Input/output, real ( kind = 8 ) DLT, the trust region radius.
!
!    Input/output, integer ( kind = 4 ) IRETCD, the status code.
!    0, xpls accepted as next iterate;  dlt trust region for next iteration.
!    1, xpls unsatisfactory but accepted as next iterate because
!      xpls-x < smallest allowable step length.
!    2, f(xpls) too large.  continue current iteration with new reduced dlt.
!    3, f(xpls) sufficiently small, but quadratic model predicts f(xpls)
!      sufficiently well to continue current iteration with new doubled dlt.
!
!    Workspace, real ( kind = 8 ) XPLSP(N), [value needs to be retained between
!    succesive calls of k-th global step].
!
!    Worskpace, real ( kind = 8 ) FPLSP, [retain value between successive
!    calls].
!
!    Output, real ( kind = 8 ) XPLS(N), the new iterate x[k].
!
!    Output, real ( kind = 8 ) FPLS, the function value at the new iterate.
!
!    Output, logical MXTAKE, is TRUE if a step of maximum length was taken.
!
!    Input, integer ( kind = 4 ) IPR, the device to which to send output.
!
!    Input, integer ( kind = 4 ) METHOD, the algorithm to use.
!    1, line search,
!    2, double dogleg,
!    3, More-Hebdon.
!
!    Input, real ( kind = 8 ) UDIAG(N), the diagonal of hessian in a(.,.)
!
SUBROUTINE TREGUP(nr, n, x, f, g, a, sc, sx, nwtake, stepmx, steptl, &
  dlt, iretcd, xplsp, fplsp, xpls, fpls, mxtake, ipr, method, udiag, &
  N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, INFO, G_CON, NCON )
IMPLICIT NONE

INTEGER, INTENT(IN) :: N, NR, IPR, METHOD, NCON
INTEGER, INTENT(INOUT) :: IRETCD
DOUBLE PRECISION, INTENT(IN) :: X(N), F, G(N), A(NR,N), SC(N), SX(N)
DOUBLE PRECISION, INTENT(IN) ::STEPMX, STEPTL
DOUBLE PRECISION, INTENT(INOUT) :: DLT, XPLSP(N), FPLSP, XPLS(N), FPLS, G_CON(NCON)
LOGICAL, INTENT(IN) :: NWTAKE
LOGICAL, INTENT(INOUT) :: MXTAKE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NEW INPUT VARIABLES
INTEGER, INTENT(IN) :: N_INT, N1, N2, CHROM_INT(N_INT), INFO
DOUBLE PRECISION, INTENT(IN) :: INPUT_ARRAY(N1, N2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION :: DLTF, DLTFP, DLTMP, RLN, SLP, TEMP, UDIAG(N)
INTEGER :: I, J

mxtake = .false.
xpls(1:n) = x(1:n) + sc(1:n)
!call COST ( n, xpls, fpls )

CALL COST(N, N_INT, N1, N2, XPLS, CHROM_INT, FPLS, INPUT_ARRAY, G_CON, NCON)
dltf = fpls - f
slp = dot_product ( g(1:n), sc(1:n) )
!
!  Next statement added for case of compilers which do not optimize
!  evaluation of next "if" statement (in which case fplsp could be
!  undefined).
!
  if ( iretcd == 4 ) then
    fplsp = 0.0D+00
  end if
!
!  Reset XPLS to XPLSP and terminate global step.
!
  if ( iretcd == 3 .and. ( fplsp <= fpls .or. 1.0D-04 * slp < dltf ) ) then
    iretcd = 0
    xpls(1:n) = xplsp(1:n)
    fpls = fplsp
    dlt = 0.5D+00 * dlt
    return
  end if

  if ( dltf <= 1.0D-04 * slp ) then
    go to 170
  end if

  rln = 0.0D+00

  do i = 1, n

    rln = max (                                        &
                rln,                                   &
                abs ( sc(i) ) / max (                  &
                                      abs ( xpls(i) ), &
                                      1.0D+00 / sx(i)  &
                                    )                  &
              )
  end do
!
!  Cannot find satisfactory xpls sufficiently distinct from x
!
  if ( rln < steptl ) then
    iretcd = 1
    return
  end if
!
!  Reduce trust region and continue global step
!
        iretcd = 2
        dltmp = -slp * dlt / ( 2.0D+00 * ( dltf - slp ) )

        if ( 0.1D+00 * dlt <= dltmp ) then
          go to 155
        end if

          dlt = 0.1D+00 * dlt
          go to 160

  155   continue
        dlt = dltmp

  160   continue
        return
!
!  FPLS sufficiently small
!
  170     continue

      dltfp = 0.0D+00

      if ( method == 2 ) then

        do i = 1, n
          temp = dot_product ( sc(i:n), a(i:n,i) )
          dltfp = dltfp + temp**2
        end do

      else

        do i = 1, n
          dltfp = dltfp + udiag(i) * sc(i) * sc(i)

          temp = 0.0D+00
          do j = i+1, n
            temp = temp + a(i,j) * sc(i) * sc(j)
          end do
          dltfp = dltfp + 2.0D+00 * temp
        end do

      end if

      dltfp = slp + dltfp / 2.0D+00

      if ( iretcd == 2 .or. &
         0.1D+00 * abs ( dltf ) < abs ( dltfp - dltf ) .or. &
         nwtake .or. &
         0.99D+00 * stepmx < dlt ) then
        go to 210
      end if
!
!  Double trust region and continue global step
!
        iretcd = 3
        xplsp(1:n) = xpls(1:n)
        fplsp = fpls
        dlt = min ( 2.0D+00 * dlt, stepmx )
        return
!
!  Accept XPLS as the next iterate.  Choose the new trust region.
!
  210       continue

  iretcd = 0

  if ( 0.99D+00 * stepmx < dlt ) then
    mxtake = .true.
  end if

  if ( dltf < 0.1D+00 * dltfp ) then
    if ( dltf <= 0.75D+00 * dltfp ) then
      dlt = min ( 2.0D+00 * dlt, stepmx )
    end if
  else
    dlt = 0.5D+00 * dlt
  end if

RETURN
END SUBROUTINE TREGUP







!*****************************************************************************80
!
!! UNCMIN minimizes a smooth nonlinear function of N variables.
!
!  Discussion:
!
!    A subroutine that computes the function value at any point
!    must be supplied.  Derivative values are not required.
!    This subroutine provides the simplest interface to the uncmin
!    minimization package.  The user has no control over options.
!
!    This routine uses a quasi-Newton algorithm with line search
!    to minimize the function represented by the subroutine COST.
!    At each iteration, the nonlinear function is approximated
!    by a quadratic function derived from a taylor series.
!    The quadratic function is minimized to obtain a search direction,
!    and an approximate minimum of the nonlinear function along
!    the search direction is found using a line search.  The
!    algorithm computes an approximation to the second derivative
!    matrix of the nonlinear function using quasi-Newton techniques.
!
!    The uncmin package is quite general, and provides many options
!    for the user.  However, this subroutine is designed to be
!    easy to use, with few choices allowed.  For example:
!
!    1.  only function values need be computed.  first derivative
!    values are obtained by finite differencing.  this can be
!    very costly when the number of variables is large.
!
!    2.  it is assumed that the function values can be obtained
!    accurately (to an accuracy comparable to the precision of
!    the computer arithmetic).
!
!    3.  at most 150 iterations are allowed.
!
!    4.  it is assumed that the function values are well-scaled,
!    that is, that the optimal function value is not pathologically
!    large or small.
!
!  Reference:
!
!    John Dennis, Robert Schnabel,
!    Numerical Methods for Unconstrained Optimization
!    and Nonlinear Equations,
!    SIAM, 1996,
!    ISBN13: 978-0-898713-64-0,
!    LC: QA402.5.D44.
!
!    Robert Schnabel, John Koontz, B E Weiss,
!    A modular system of algorithms for unconstrained minimization,
!    Report cu-cs-240-82,
!    Computer Science Department,
!    University of Colorado at Boulder, 1982.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the dimension of the problem.
!
!    Input, real ( kind = 8 ) X0(N), an initial estimate of the minimizer.
!
!    Input, external COST, the name of the routine to evaluate the minimization
!    function, of the form
!
!      subroutine COST ( n, x, f )
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) x(n)
!      real ( kind = 8 ) f
!
!    Output, real ( kind = 8 ) X(N), the local minimizer.
!
!    Output, real ( kind = 8 ) F, the function value at X.
!
!    Output, integer ( kind = 4 ) INFO, termination code.
!    0:  optimal solution found.
!    1:  terminated with gradient small, X is probably optimal.
!    2:  terminated with stepsize small, X is probably optimal.
!    3:  lower point cannot be found, X is probably optimal.
!    4:  iteration limit (150) exceeded.
!    5:  too many large steps, function may be unbounded.
!    -1:  insufficient workspace.
!
!    Workspace, real ( kind = 8 ) W(LW).
!
!    Input, integer ( kind = 4 ) LW, the size of the workspace vector W, which
!    must be at least N * ( N + 10 ).
!
!SUBROUTINE UNCMIN( n, x0, x, f, info, w, lw )
SUBROUTINE UNCMIN(N, N_INT, N1, N2, ITNLIM, INFO, X0, X, CHROM_INT, &
                 F, INPUT_ARRAY, G_CON, NCON)

IMPLICIT NONE
INTEGER, INTENT(IN) :: N, NCON
INTEGER, INTENT(INOUT) :: INFO, ITNLIM
DOUBLE PRECISION, INTENT(INOUT) :: X(N), F, X0(N), G_CON(NCON)
DOUBLE PRECISION :: DLT, EPSM, FSCALE, GRADTL, STEPMX, STEPTL
INTEGER :: IA, IAGFLG, IAHFLG, IEXP, IG, IPR, IT, IW1, IW2, IW3
INTEGER :: IW4, IW5, IW6, IW7, IW8, LWMIN, METHOD, MSG, NDIGIT, NR
DOUBLE PRECISION :: W(N*(N+10))
INTEGER :: LW

! NEW PARTS NEEDED FOR THE COST FUNCTION AND COMBINED ALGORITHM
INTEGER, INTENT(IN) :: N_INT, N1, N2, CHROM_INT(N_INT)
DOUBLE PRECISION, INTENT(IN) :: INPUT_ARRAY(N1,N2)
!
!  Subdivide workspace
!

! ADDED W AND LW HERE INSTEAD OF PASSING THE ARRAY AND LENGTH INTO UNCMIN
LW=N*(N+10)


  ig  = 1
  it  = ig  + n
  iw1 = it  + n
  iw2 = iw1 + n
  iw3 = iw2 + n
  iw4 = iw3 + n
  iw5 = iw4 + n
  iw6 = iw5 + n
  iw7 = iw6 + n
  iw8 = iw7 + n
  ia  = iw8 + n
  lwmin = ia + n*n-1

  if ( lw < lwmin ) then
    info = -1
    !write ( *, '(a)' ) ' '
    !write ( *, '(a)' ) 'UNCMIN - Fatal error!'
    !write ( *, '(a)' ) '  Insufficient workspace.'
    !write ( *, '(a)' ) '  LW < LWMIN.'
    !write ( *, '(a,i6)' ) '  LW = ', lw
    !write ( *, '(a,i6)' ) '  LWMIN = ', lwmin
    stop
  end if
!
!  Set up parameters for .
!
!  parameters that should not be changed when using condensed code
!
! nr     = parameter used to divide workspace
! method = 1 (line search) -- do not change
! msg    = 9 => no printing, n=1 allowed
! iahflg = 1 => analytic hessian  supplied (0 otherwise)
! ipr    = device for output (irrelevant in current version)
! dlt    = (irrelevant parameter for method = 1)
! epsm   = machine epsilon
!
NR=N
METHOD=1
MSG=9
IAHFLG=0
IPR=6
DLT=-1.0D+00
EPSM=epsilon(EPSM)
!
! parameters that may be changed:
!
! iexp   = 1 means function expensive to evaluate (iexp = 0 => cheap)
! iagflg = 1 means analytic gradient supplied (0 otherwise)
! ndigit = -1 means  assumes f is fully accurate
! itnlim = 150 = maximum number of iterations allowed
! gradtl = zero tolerance for gradient, for convergence tests
! stepmx = maximum allowable step size
! steptl = zero tolerance for step, for convergence tests
! fscale = typical order-of-magnitude size of function
! typsiz = typical order-of-magnitude size of x (stored in w(lt))
!
IEXP=1
IAGFLG=0
NDIGIT=-1
!ITNLIM=5
GRADTL=epsm**(1.0D+00 / 3.0D+00 )
  stepmx = 0.0D+00
  steptl = sqrt ( epsm )
  fscale = 1.0D+00
  w(it:it+n-1) = 1.0D+00
!
!  Minimize function
!

!!WRITE(*,*) N, N_INT, N1, N2, INFO
!!WRITE(*,*) X0
!!WRITE(*,*) X
!!WRITE(*,*) CHROM_INT, F
!WRITE(*,*) "call OPTDRV"
  call optdrv ( nr, n, x0, w(it), fscale, method, iexp, &
    msg, ndigit, itnlim, iagflg, iahflg,  ipr, dlt, gradtl, stepmx, steptl, &
    x, f, w(ig), info, w(ia), w(iw1), w(iw2), w(iw3), w(iw4), &
    w(iw5), w(iw6), w(iw7), w(iw8), N_INT, N1, N2, CHROM_INT, INPUT_ARRAY, &
    INFO, G_CON, NCON)

!WRITE(*,*) F

!!WRITE(*,*)f
!  if ( info == 1 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'UNCMIN - Note!'
!    write ( *, '(a)' ) '  INFO = 1.'
!    write ( *, '(a)' ) '  The iteration probably converged.'
!    write ( *, '(a)' ) '  The gradient is very small.'
!    return
!  end if
!
!  if ( info == 2 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'UNCMIN - Note!'
!    write ( *, '(a)' ) '  INFO = 2.'
!    write ( *, '(a)' ) '  The iteration probably converged.'
!    write ( *, '(a)' ) '  The stepsize is very small.'
!    return
!  end if
!
!  if ( info == 3 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'UNCMIN - Warning!'
!    write ( *, '(a)' ) '  INFO = 3.'
!    write ( *, '(a)' ) '  Cannot find a point with lower value.'
!    write ( *, '(a)' ) '  (But not completely happy with the current value.)'
!    return
!  end if
!
!  if ( info == 4 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'UNCMIN - Warning!'
!    write ( *, '(a)' ) '  INFO = 4.'
!    write ( *, '(a)' ) '  Too many iterations.'
!    return
!  end if
!
!  if ( info == 5 ) then
!    write ( *, '(a)' ) ' '
!    write ( *, '(a)' ) 'UNCMIN - Warning!'
!    write ( *, '(a)' ) '  INFO = 5.'
!    write ( *, '(a)' ) '  Too many large steps.'
!    write ( *, '(a)' ) '  The function may be unbounded.'
!    return
!  end if

RETURN
END SUBROUTINE UNCMIN    

!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!  THIS IS A DUMMY FUCNTION TO EVALUATE THE GRADIENT OF THE FUNCTION
!
!   G(I) = d F/d X(I).
!
SUBROUTINE D1FCN( N, X, G )
IMPLICIT NONE
INTEGER, INTENT(IN) :: N
DOUBLE PRECISION, INTENT(IN) :: X(N)
DOUBLE PRECISION, INTENT(INOUT) :: G(N)

write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'D1FCN - Fatal error!'
write ( *, '(a)' ) '  This is a dummy routine.'
write ( *, '(a)' ) '  The user is required to replace it with a'
write ( *, '(a)' ) '  routine that computes the gradient of F.'

RETURN
END SUBROUTINE D1FCN



!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!  THIS IS A DUMMY FUCNTION TO EVALUATE THE HESSIAN OF THE FUNCTION
!
!   H(I,J) = d d F/d X(I) d X(J)
!
!   IN THIS CASE:
!     H=A(NR,N) AND NR SHOULD BE EQUAL TO N
!
SUBROUTINE D2FCN( NR, N, X, A )
IMPLICIT NONE
INTEGER, INTENT(IN) :: N, NR
DOUBLE PRECISION, INTENT(IN) :: X(N)
DOUBLE PRECISION, INTENT(INOUT) :: A(NR,N)

write ( *, '(a)' ) ' '
write ( *, '(a)' ) 'D2FCN - Fatal error!'
write ( *, '(a)' ) '  This is a dummy routine.'
write ( *, '(a)' ) '  The user is required to replace it with a'
write ( *, '(a)' ) '  routine that computes the Hessian matrix of F.'

RETURN
END SUBROUTINE D2FCN                 

END MODULE UNCMIN_MOD

MODULE cobyla2
USE COST_MODULE
! The Fortran 77 version of this code was by Michael Powell
! (M.J.D.Powell @ damtp.cam.ac.uk)

! This Fortran 90 version by Alan Miller
! This version is in a subset of F90 and is compatible with version 3.0a
! of Lahey's free ELF90 compiler.
! Alan.Miller @ vic.cmis.csiro.au
! Latest revision - 28 June 2000

IMPLICIT NONE

!INTEGER, PARAMETER, PUBLIC :: dp = SELECTED_REAL_KIND(14, 60)

    CONTAINS
!
!********************************************************************!
!********************************************************************! 
!********************************************************************! 
!********************************************************************! 
! THIS SUBROUTINE CONTAINS THE DRIVER FOR THE COBYLA OPTIMIZATION 
!  PACKAGE.  THE ONLY MODIFICATIONS REQUIRED FOR THE COBYLA SOLVER IS 
!  MODIFYING IT TO CALL THE COST FUNCTION AS REQURIED FOR THE GA-NLP 
!  ALGORITHM
SUBROUTINE COBYLA_DRIVER(N, N_INT, N1, N2, x, CHROM_INT, FITNESS, &
    ARRAY, NCON, NLP_ITER_MAX )
IMPLICIT NONE

INTEGER, INTENT(IN) :: N, N_INT, N1, N2, CHROM_INT(N_INT), NCON, &
    NLP_ITER_MAX
DOUBLE PRECISION, INTENT(INOUT) :: FITNESS, X(N)

DOUBLE PRECISION, INTENT(IN) :: ARRAY(N1,N2)

DOUBLE PRECISION :: RHOBEG, RHOEND

INTEGER :: IPRINT, MAXFUN

MAXFUN=NLP_ITER_MAX
RHOBEG = 0.5D0
RHOEND = 1.D-6
IPRINT = 0

CALL COBYLA (N, NCON, X, RHOBEG, RHOEND, IPRINT, MAXFUN, N1, N2, &
    N_INT, CHROM_INT, ARRAY, FITNESS)

END SUBROUTINE COBYLA_DRIVER
!------------------------------------------------------------------------
SUBROUTINE cobyla (n, m, x, rhobeg, rhoend, iprint, maxfun, N1, N2, N_INT, CHROM_INT, ARRAY, f)

INTEGER, INTENT(IN)        :: n, m, N1, N2, N_INT, CHROM_INT(N_INT)
DOUBLE PRECISION, INTENT(IN) :: ARRAY(N1,N2)
DOUBLE PRECISION, INTENT(IN OUT)  :: x(n), rhobeg, rhoend, f 
!DOUBLE PRECISION, INTENT(IN OUT)  :: rhobeg
!DOUBLE PRECISION, INTENT(IN OUT)  :: rhoend
INTEGER, INTENT(IN)        :: iprint
INTEGER, INTENT(IN OUT)    :: maxfun

!  This subroutine minimizes an objective function F(X) subject to M
!  inequality constraints on X, where X is a vector of variables that has
!  N components.  The algorithm employs linear approximations to the
!  objective and constraint functions, the approximations being formed by
!  linear interpolation at N+1 points in the space of the variables.
!  We regard these interpolation points as vertices of a simplex.  The
!  parameter RHO controls the size of the simplex and it is reduced
!  automatically from RHOBEG to RHOEND.  For each RHO the subroutine tries
!  to achieve a good vector of variables for the current size, and then
!  RHO is reduced until the value RHOEND is reached.  Therefore RHOBEG and
!  RHOEND should be set to reasonable initial changes to and the required
!  accuracy in the variables respectively, but this accuracy should be
!  viewed as a subject for experimentation because it is not guaranteed.
!  The subroutine has an advantage over many of its competitors, however,
!  which is that it treats each constraint individually when calculating
!  a change to the variables, instead of lumping the constraints together
!  into a single penalty function.  The name of the subroutine is derived
!  from the phrase Constrained Optimization BY Linear Approximations.

!  The user must set the values of N, M, RHOBEG and RHOEND, and must
!  provide an initial vector of variables in X.  Further, the value of
!  IPRINT should be set to 0, 1, 2 or 3, which controls the amount of
!  printing during the calculation. Specifically, there is no output if
!  IPRINT=0 and there is output only at the end of the calculation if
!  IPRINT=1.  Otherwise each new value of RHO and SIGMA is printed.
!  Further, the vector of variables and some function information are
!  given either when RHO is reduced or when each new value of F(X) is
!  computed in the cases IPRINT=2 or IPRINT=3 respectively. Here SIGMA
!  is a penalty parameter, it being assumed that a change to X is an
!  improvement if it reduces the merit function
!             F(X)+SIGMA*MAX(0.0, - C1(X), - C2(X),..., - CM(X)),
!  where C1,C2,...,CM denote the constraint functions that should become
!  nonnegative eventually, at least to the precision of RHOEND. In the
!  printed output the displayed term that is multiplied by SIGMA is
!  called MAXCV, which stands for 'MAXimum Constraint Violation'.  The
!  argument MAXFUN is an integer variable that must be set by the user to a
!  limit on the number of calls of CALCFC, the purpose of this routine being
!  given below.  The value of MAXFUN will be altered to the number of calls
!  of CALCFC that are made.

!  In order to define the objective and constraint functions, we require
!  a subroutine that has the name and arguments
!             SUBROUTINE CALCFC (N,M,X,F,CON)
!             DIMENSION X(:),CON(:)  .
!  The values of N and M are fixed and have been defined already, while
!  X is now the current vector of variables. The subroutine should return
!  the objective and constraint functions at X in F and CON(1),CON(2),
!  ...,CON(M).  Note that we are trying to adjust X so that F(X) is as
!  small as possible subject to the constraint functions being nonnegative.

!  N.B. If the starting value for any x(i) is set to zero, that value will
!       not be changed.   This can be a useful feature in comparing
!       nested models.   If all the x(i)'s are set to zero, an error message
!       will result.

! Local variable

INTEGER :: mpp

mpp = m + 2
CALL cobylb (n, m, mpp, x, rhobeg, rhoend, iprint, maxfun, N1, N2, N_INT, CHROM_INT, ARRAY,f)

RETURN
END SUBROUTINE cobyla


!------------------------------------------------------------------------------

SUBROUTINE cobylb (n, m, mpp, x, rhobeg, rhoend, iprint, maxfun, N1, N2, N_INT, CHROM_INT, ARRAY,f)

! N.B. Arguments CON, SIM, SIMI, DATMAT, A, VSIG, VETA, SIGBAR, DX, W & IACT
!   have been removed.

INTEGER, INTENT(IN)        :: n, m, mpp, N1, N2, N_INT, CHROM_INT(N_INT)
DOUBLE PRECISION, INTENT(IN OUT)  :: x(n), f
DOUBLE PRECISION, INTENT(IN)      :: rhobeg, ARRAY(N1,N2)
DOUBLE PRECISION, INTENT(IN)      :: rhoend
INTEGER, INTENT(IN)        :: iprint
INTEGER, INTENT(OUT)       :: maxfun

!INTERFACE
!  SUBROUTINE calcfc (n, m, x, f, con)
!    IMPLICIT NONE
!    !INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)
!    INTEGER, INTENT(IN)    :: n, m
!    DOUBLE PRECISION, INTENT(IN)  :: x(n)
!    DOUBLE PRECISION, INTENT(OUT) :: f
!    DOUBLE PRECISION, INTENT(OUT) :: con(m)
!  END SUBROUTINE calcfc
!END INTERFACE

!  Set the initial values of some parameters. The last column of SIM holds
!  the optimal vertex of the current simplex, and the preceding N columns
!  hold the displacements from the optimal vertex to the other vertices.
!  Further, SIMI holds the inverse of the matrix that is contained in the
!  first N columns of SIM.

! Local variables

DOUBLE PRECISION :: con(mpp), sim(n,n+1), simi(n,n), datmat(mpp,n+1), a(n,m+1),      &
             vsig(n), veta(n), sigbar(n), dx(n), w(n)
DOUBLE PRECISION :: alpha, barmu, beta, cmin, cmax, cvmaxm, cvmaxp, delta, denom,    &
             dxsign, edgmax, error, gamma, pareta, parmu, parsig, phi,     &
             phimin, prerec, prerem, ratio, resmax, resnew, rho, temp, tempa, &
             total, trured, vmnew, vmold, weta, wsig
INTEGER   :: i, ibrnch, iflag, ifull, iptem, iptemp, j, jdrop, k, l, mp,  &
             nbest, nfvals, np

iptem = MIN(n,5)
iptemp = iptem + 1
np = n + 1
mp = m + 1
alpha = 0.25d0
beta = 2.1d0
gamma = 0.5D0
delta = 1.1D0
rho = rhobeg
parmu = 0.0D0
IF (iprint >= 2) WRITE(*, 10) rho
10 FORMAT (/'   The initial value of RHO is', G13.6,   &
           '  and PARMU is set to zero.')
nfvals = 0
temp = 1.0D0/rho
DO i=1,n
  sim(i,np) = x(i)
  DO j=1,n
    sim(i,j) = 0.0D0
    simi(i,j) = 0.0D0
  END DO
  sim(i,i) = rho
  simi(i,i) = temp
END DO
jdrop = np
ibrnch = 0

!  Make the next call of the user-supplied subroutine CALCFC. These
!  instructions are also used for calling CALCFC during the iterations of
!  the algorithm.

40 IF (nfvals >= maxfun .AND. nfvals > 0) THEN
  IF (iprint >= 1) WRITE(*, 50)
  50 FORMAT (/'   Return from subroutine COBYLA because the ',  &
             'MAXFUN limit has been reached.')
  GO TO 600
END IF
nfvals = nfvals + 1

CALL COST(n, N_INT, N1, N2, X, CHROM_INT, f, ARRAY, CON(1:m), m )
!CALL COST (n, m, x, f, con(1:m))
resmax = 0.0D0
IF (m > 0) THEN
  DO k=1,m
    resmax = MAX(resmax, - con(k))
  END DO
END IF
IF (nfvals == iprint-1 .OR. iprint == 3) THEN
  WRITE(*, 70) nfvals, f, resmax, x(1:iptem)
  70 FORMAT (/'   NFVALS = ', i5, '   F = ', G13.6, '    MAXCV = ',  &
             G13.6/ ('   X = ', 5G14.6))
  IF (iptem < n) WRITE(*, 80) x(iptemp:n)
  80 FORMAT (G19.6, G15.6)
END IF
con(mp) = f
con(mpp) = resmax
IF (ibrnch == 1) GO TO 440

!  Set the recently calculated function values in a column of DATMAT. This
!  array has a column for each vertex of the current simplex, the entries of
!  each column being the values of the constraint functions (if any)
!  followed by the objective function and the greatest constraint violation
!  at the vertex.

DO k=1,mpp
  datmat(k,jdrop) = con(k)
END DO
IF (nfvals > np) GO TO 130

!  Exchange the new vertex of the initial simplex with the optimal vertex if
!  necessary. Then, if the initial simplex is not complete, pick its next
!  vertex and calculate the function values there.

IF (jdrop <= n) THEN
  IF (datmat(mp,np) <= f) THEN
    x(jdrop) = sim(jdrop,np)
  ELSE
    sim(jdrop,np) = x(jdrop)
    DO k=1,mpp
      datmat(k,jdrop) = datmat(k,np)
      datmat(k,np) = con(k)
    END DO
    DO k=1,jdrop
      sim(jdrop,k) = -rho
      temp = -SUM( simi(k:jdrop, k) )
      simi(jdrop,k) = temp
    END DO
  END IF
END IF
IF (nfvals <= n) THEN
  jdrop = nfvals
  x(jdrop) = x(jdrop) + rho
  GO TO 40
END IF
130 ibrnch = 1

!  Identify the optimal vertex of the current simplex.

140 phimin = datmat(mp,np) + parmu*datmat(mpp,np)
nbest = np
DO j=1,n
  temp = datmat(mp,j) + parmu*datmat(mpp,j)
  IF (temp < phimin) THEN
    nbest = j
    phimin = temp
  ELSE IF (temp == phimin .AND. parmu == 0.0D0) THEN
    IF (datmat(mpp,j) < datmat(mpp,nbest)) nbest = j
  END IF
END DO

!  Switch the best vertex into pole position if it is not there already,
!  and also update SIM, SIMI and DATMAT.

IF (nbest <= n) THEN
  DO i=1,mpp
    temp = datmat(i,np)
    datmat(i,np) = datmat(i,nbest)
    datmat(i,nbest) = temp
  END DO
  DO i=1,n
    temp = sim(i,nbest)
    sim(i,nbest) = 0.0D0
    sim(i,np) = sim(i,np) + temp
    tempa = 0.0D0
    DO k=1,n
      sim(i,k) = sim(i,k) - temp
      tempa = tempa - simi(k,i)
    END DO
    simi(nbest,i) = tempa
  END DO
END IF

!  Make an error return if SIGI is a poor approximation to the inverse of
!  the leading N by N submatrix of SIG.

error = 0.0D0
DO i=1,n
  DO j=1,n
    temp = 0.0D0
    IF (i == j) temp = temp - 1.0D0
    temp = temp + DOT_PRODUCT( simi(i,1:n), sim(1:n,j) )
    error = MAX(error, ABS(temp))
  END DO
END DO
IF (error > 0.1D0) THEN
  IF (iprint >= 1) WRITE(*, 210)
  210 FORMAT (/'   Return from subroutine COBYLA because ',  &
              'rounding errors are becoming damaging.')
  GO TO 600
END IF

!  Calculate the coefficients of the linear approximations to the objective
!  and constraint functions, placing minus the objective function gradient
!  after the constraint gradients in the array A. The vector W is used for
!  working space.

DO k=1,mp
  con(k) = -datmat(k,np)
  DO j=1,n
    w(j) = datmat(k,j) + con(k)
  END DO
  DO i=1,n
    temp = DOT_PRODUCT( w(1:n), simi(1:n,i) )
    IF (k == mp) temp = -temp
    a(i,k) = temp
  END DO
END DO

!  Calculate the values of sigma and eta, and set IFLAG = 0 if the current
!  simplex is not acceptable.

iflag = 1
parsig = alpha*rho
pareta = beta*rho
DO j=1,n
  wsig = SUM( simi(j,1:n)**2 )
  weta = SUM( sim(1:n,j)**2 )
  vsig(j) = 1.0D0/SQRT(wsig)
  veta(j) = SQRT(weta)
  IF (vsig(j) < parsig .OR. veta(j) > pareta) iflag = 0
END DO

!  If a new vertex is needed to improve acceptability, then decide which
!  vertex to drop from the simplex.

IF (ibrnch == 1 .OR. iflag == 1) GO TO 370
jdrop = 0
temp = pareta
DO j=1,n
  IF (veta(j) > temp) THEN
    jdrop = j
    temp = veta(j)
  END IF
END DO
IF (jdrop == 0) THEN
  DO j=1,n
    IF (vsig(j) < temp) THEN
      jdrop = j
      temp = vsig(j)
    END IF
  END DO
END IF

!  Calculate the step to the new vertex and its sign.

temp = gamma*rho*vsig(jdrop)
dx(1:n) = temp*simi(jdrop,1:n)
cvmaxp = 0.0D0
cvmaxm = 0.0D0
DO k=1,mp
  total = DOT_PRODUCT( a(1:n,k), dx(1:n) )
  IF (k < mp) THEN
    temp = datmat(k,np)
    cvmaxp = MAX(cvmaxp, -total - temp)
    cvmaxm = MAX(cvmaxm, total - temp)
  END IF
END DO
dxsign = 1.0D0
IF (parmu*(cvmaxp - cvmaxm) > total + total) dxsign = -1.0D0

!  Update the elements of SIM and SIMI, and set the next X.

temp = 0.0D0
DO i=1,n
  dx(i) = dxsign*dx(i)
  sim(i,jdrop) = dx(i)
  temp = temp + simi(jdrop,i)*dx(i)
END DO
simi(jdrop,1:n) = simi(jdrop,1:n) / temp
DO j=1,n
  IF (j /= jdrop) THEN
    temp = DOT_PRODUCT( simi(j,1:n), dx(1:n) )
    simi(j,1:n) = simi(j,1:n) - temp*simi(jdrop,1:n)
  END IF
  x(j) = sim(j,np) + dx(j)
END DO
GO TO 40

!  Calculate DX = x(*)-x(0).
!  Branch if the length of DX is less than 0.5*RHO.

370 CALL trstlp (n, m, a, con, rho, dx, ifull)
IF (ifull == 0) THEN
  temp = SUM( dx(1:n)**2 )
  IF (temp < 0.25D0*rho*rho) THEN
    ibrnch = 1
    GO TO 550
  END IF
END IF

!  Predict the change to F and the new maximum constraint violation if the
!  variables are altered from x(0) to x(0) + DX.

resnew = 0.0D0
con(mp) = 0.0D0
DO k=1,mp
  total = con(k) - DOT_PRODUCT( a(1:n,k), dx(1:n) )
  IF (k < mp) resnew = MAX(resnew, total)
END DO

!  Increase PARMU if necessary and branch back if this change alters the
!  optimal vertex. Otherwise PREREM and PREREC will be set to the predicted
!  reductions in the merit function and the maximum constraint violation
!  respectively.

barmu = 0.0D0
prerec = datmat(mpp,np) - resnew
IF (prerec > 0.0D0) barmu = total/prerec
IF (parmu < 1.5D0*barmu) THEN
  parmu = 2.0D0*barmu
  IF (iprint >= 2) WRITE(*, 410) parmu
  410 FORMAT (/'   Increase in PARMU to', G13.6)
  phi = datmat(mp,np) + parmu*datmat(mpp,np)
  DO j=1,n
    temp = datmat(mp,j) + parmu*datmat(mpp,j)
    IF (temp < phi) GO TO 140
    IF (temp == phi .AND. parmu == 0.0) THEN
      IF (datmat(mpp,j) < datmat(mpp,np)) GO TO 140
    END IF
  END DO
END IF
prerem = parmu*prerec - total

!  Calculate the constraint and objective functions at x(*).
!  Then find the actual reduction in the merit function.

x(1:n) = sim(1:n,np) + dx(1:n)
ibrnch = 1
GO TO 40

440 vmold = datmat(mp,np) + parmu*datmat(mpp,np)
vmnew = f + parmu*resmax
trured = vmold - vmnew
IF (parmu == 0.0D0 .AND. f == datmat(mp,np)) THEN
  prerem = prerec
  trured = datmat(mpp,np) - resmax
END IF

!  Begin the operations that decide whether x(*) should replace one of the
!  vertices of the current simplex, the change being mandatory if TRURED is
!  positive. Firstly, JDROP is set to the index of the vertex that is to be
!  replaced.

ratio = 0.0D0
IF (trured <= 0.0) ratio = 1.0
jdrop = 0
DO j=1,n
  temp = DOT_PRODUCT( simi(j,1:n), dx(1:n) )
  temp = ABS(temp)
  IF (temp > ratio) THEN
    jdrop = j
    ratio = temp
  END IF
  sigbar(j) = temp*vsig(j)
END DO

!  Calculate the value of ell.

edgmax = delta*rho
l = 0
DO j=1,n
  IF (sigbar(j) >= parsig .OR. sigbar(j) >= vsig(j)) THEN
    temp = veta(j)
    IF (trured > 0.0D0) THEN
      temp = SUM( (dx(1:n) - sim(1:n,j))**2 )
      temp = SQRT(temp)
    END IF
    IF (temp > edgmax) THEN
      l = j
      edgmax = temp
    END IF
  END IF
END DO
IF (l > 0) jdrop = l
IF (jdrop == 0) GO TO 550

!  Revise the simplex by updating the elements of SIM, SIMI and DATMAT.

temp = 0.0D0
DO i=1,n
  sim(i,jdrop) = dx(i)
  temp = temp + simi(jdrop,i)*dx(i)
END DO
simi(jdrop,1:n) = simi(jdrop,1:n) / temp
DO j=1,n
  IF (j /= jdrop) THEN
    temp = DOT_PRODUCT( simi(j,1:n), dx(1:n) )
    simi(j,1:n) = simi(j,1:n) - temp*simi(jdrop,1:n)
  END IF
END DO
datmat(1:mpp,jdrop) = con(1:mpp)

!  Branch back for further iterations with the current RHO.

IF (trured > 0.0D0 .AND. trured >= 0.1D0*prerem) GO TO 140
550 IF (iflag == 0) THEN
  ibrnch = 0
  GO TO 140
END IF

!  Otherwise reduce RHO if it is not at its least value and reset PARMU.

IF (rho > rhoend) THEN
  rho = 0.5D0*rho
  IF (rho <= 1.5D0*rhoend) rho = rhoend
  IF (parmu > 0.0D0) THEN
    denom = 0.0D0
    DO k=1,mp
      cmin = datmat(k,np)
      cmax = cmin
      DO i=1,n
        cmin = MIN(cmin, datmat(k,i))
        cmax = MAX(cmax, datmat(k,i))
      END DO
      IF (k <= m .AND. cmin < 0.5D0*cmax) THEN
        temp = MAX(cmax,0.0D0) - cmin
        IF (denom <= 0.0D0) THEN
          denom = temp
        ELSE
          denom = MIN(denom,temp)
        END IF
      END IF
    END DO
    IF (denom == 0.0D0) THEN
      parmu = 0.0D0
    ELSE IF (cmax - cmin < parmu*denom) THEN
      parmu = (cmax - cmin)/denom
    END IF
  END IF
  IF (iprint >= 2) WRITE(*, 580) rho,parmu
  580 FORMAT (/'   Reduction in RHO to ', G13.6, '  and PARMU = ', G13.6)
  IF (iprint == 2) THEN
    WRITE(*, 70) nfvals, datmat(mp,np), datmat(mpp,np), sim(1:iptem,np)
    IF (iptem < n) WRITE(*, 80) x(iptemp:n)
  END IF
  GO TO 140
END IF

!  Return the best calculated values of the variables.

IF (iprint >= 1) WRITE(*, 590)
590 FORMAT (/'   Normal return from subroutine COBYLA')
IF (ifull == 1) GO TO 620

600 x(1:n) = sim(1:n,np)
f = datmat(mp,np)
resmax = datmat(mpp,np)
620 IF (iprint >= 1) THEN
  WRITE(*, 70) nfvals, f, resmax, x(1:iptem)
  IF (iptem < n) WRITE(*, 80) x(iptemp:n)
END IF
maxfun = nfvals

RETURN
END SUBROUTINE cobylb
!------------------------------------------------------------------------------

SUBROUTINE trstlp (n, m, a, b, rho, dx, ifull)

! N.B. Arguments Z, ZDOTA, VMULTC, SDIRN, DXNEW, VMULTD & IACT have been removed.

INTEGER, INTENT(IN)     :: n
INTEGER, INTENT(IN)     :: m
DOUBLE PRECISION, INTENT(IN)   :: a(:,:)
DOUBLE PRECISION, INTENT(IN)   :: b(:)
DOUBLE PRECISION, INTENT(IN)   :: rho
DOUBLE PRECISION, INTENT(OUT)  :: dx(:)
INTEGER, INTENT(OUT)    :: ifull

!  This subroutine calculates an N-component vector DX by applying the
!  following two stages. In the first stage, DX is set to the shortest
!  vector that minimizes the greatest violation of the constraints
!    A(1,K)*DX(1)+A(2,K)*DX(2)+...+A(N,K)*DX(N) .GE. B(K), K = 2,3,...,M,
!  subject to the Euclidean length of DX being at most RHO. If its length is
!  strictly less than RHO, then we use the resultant freedom in DX to
!  minimize the objective function
!           -A(1,M+1)*DX(1) - A(2,M+1)*DX(2) - ... - A(N,M+1)*DX(N)
!  subject to no increase in any greatest constraint violation. This
!  notation allows the gradient of the objective function to be regarded as
!  the gradient of a constraint. Therefore the two stages are distinguished
!  by MCON .EQ. M and MCON .GT. M respectively. It is possible that a
!  degeneracy may prevent DX from attaining the target length RHO. Then the
!  value IFULL = 0 would be set, but usually IFULL = 1 on return.

!  In general NACT is the number of constraints in the active set and
!  IACT(1),...,IACT(NACT) are their indices, while the remainder of IACT
!  contains a permutation of the remaining constraint indices.  Further, Z
!  is an orthogonal matrix whose first NACT columns can be regarded as the
!  result of Gram-Schmidt applied to the active constraint gradients.  For
!  J = 1,2,...,NACT, the number ZDOTA(J) is the scalar product of the J-th
!  column of Z with the gradient of the J-th active constraint.  DX is the
!  current vector of variables and here the residuals of the active
!  constraints should be zero. Further, the active constraints have
!  nonnegative Lagrange multipliers that are held at the beginning of
!  VMULTC. The remainder of this vector holds the residuals of the inactive
!  constraints at DX, the ordering of the components of VMULTC being in
!  agreement with the permutation of the indices of the constraints that is
!  in IACT. All these residuals are nonnegative, which is achieved by the
!  shift RESMAX that makes the least residual zero.

!  Initialize Z and some other variables. The value of RESMAX will be
!  appropriate to DX = 0, while ICON will be the index of a most violated
!  constraint if RESMAX is positive. Usually during the first stage the
!  vector SDIRN gives a search direction that reduces all the active
!  constraint violations by one simultaneously.

! Local variables

DOUBLE PRECISION :: z(n,n), zdota(m+1), vmultc(m+1), sdirn(n), dxnew(n), vmultd(m+1)
DOUBLE PRECISION :: acca, accb, alpha, beta, dd, optnew, optold, ratio, resmax,   &
             resold, sd, sp, spabs, ss, step, stpful, sumabs, temp, tempa, &
             tot, total, vsave, zdotv, zdotw, zdvabs, zdwabs
INTEGER   :: i, iact(m+1), icon, icount, isave, k, kk, kl, kp, kw, mcon,   &
             nact, nactx

ifull = 1
mcon = m
nact = 0
resmax = 0.0D0
DO i=1,n
  z(i,1:n) = 0.0D0
  z(i,i) = 1.0D0
  dx(i) = 0.0D0
END DO
IF (m >= 1) THEN
  DO k=1,m
    IF (b(k) > resmax) THEN
      resmax = b(k)
      icon = k
    END IF
  END DO
  DO k=1,m
    iact(k) = k
    vmultc(k) = resmax - b(k)
  END DO
END IF
IF (resmax == 0.0D0) GO TO 480
sdirn(1:n) = 0.0D0

!  End the current stage of the calculation if 3 consecutive iterations
!  have either failed to reduce the best calculated value of the objective
!  function or to increase the number of active constraints since the best
!  value was calculated. This strategy prevents cycling, but there is a
!  remote possibility that it will cause premature termination.

60 optold = 0.0D0
icount = 0
70 IF (mcon == m) THEN
  optnew = resmax
ELSE
  optnew = - DOT_PRODUCT( dx(1:n), a(1:n,mcon) )
END IF
IF (icount == 0 .OR. optnew < optold) THEN
  optold = optnew
  nactx = nact
  icount = 3
ELSE IF (nact > nactx) THEN
  nactx = nact
  icount = 3
ELSE
  icount = icount - 1
  IF (icount == 0) GO TO 490
END IF

!  If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to
!  the active set. Apply Givens rotations so that the last N-NACT-1 columns
!  of Z are orthogonal to the gradient of the new constraint, a scalar
!  product being set to zero if its nonzero value could be due to computer
!  rounding errors. The array DXNEW is used for working space.

IF (icon <= nact) GO TO 260
kk = iact(icon)
dxnew(1:n) = a(1:n,kk)
tot = 0.0D0
k = n
100 IF (k > nact) THEN
  sp = 0.0D0
  spabs = 0.0D0
  DO i=1,n
    temp = z(i,k)*dxnew(i)
    sp = sp + temp
    spabs = spabs + ABS(temp)
  END DO
  acca = spabs + 0.1D0*ABS(sp)
  accb = spabs + 0.2D0*ABS(sp)
  IF (spabs >= acca .OR. acca >= accb) sp = 0.0D0
  IF (tot == 0.0D0) THEN
    tot = sp
  ELSE
    kp = k + 1
    temp = SQRT(sp*sp + tot*tot)
    alpha = sp/temp
    beta = tot/temp
    tot = temp
    DO i=1,n
      temp = alpha*z(i,k) + beta*z(i,kp)
      z(i,kp) = alpha*z(i,kp) - beta*z(i,k)
      z(i,k) = temp
    END DO
  END IF
  k = k - 1
  GO TO 100
END IF

!  Add the new constraint if this can be done without a deletion from the
!  active set.

IF (tot /= 0.0D0) THEN
  nact = nact + 1
  zdota(nact) = tot
  vmultc(icon) = vmultc(nact)
  vmultc(nact) = 0.0D0
  GO TO 210
END IF

!  The next instruction is reached if a deletion has to be made from the
!  active set in order to make room for the new active constraint, because
!  the new constraint gradient is a linear combination of the gradients of
!  the old active constraints.  Set the elements of VMULTD to the multipliers
!  of the linear combination.  Further, set IOUT to the index of the
!  constraint to be deleted, but branch if no suitable index can be found.

ratio = -1.0D0
k = nact
130 zdotv = 0.0D0
zdvabs = 0.0D0
DO i=1,n
  temp = z(i,k)*dxnew(i)
  zdotv = zdotv + temp
  zdvabs = zdvabs + ABS(temp)
END DO
acca = zdvabs + 0.1D0*ABS(zdotv)
accb = zdvabs + 0.2D0*ABS(zdotv)
IF (zdvabs < acca .AND. acca < accb) THEN
  temp = zdotv/zdota(k)
  IF (temp > 0.0D0 .AND. iact(k) <= m) THEN
    tempa = vmultc(k)/temp
    IF (ratio < 0.0D0 .OR. tempa < ratio) THEN
      ratio = tempa
    END IF
  END IF
  IF (k >= 2) THEN
    kw = iact(k)
    dxnew(1:n) = dxnew(1:n) - temp*a(1:n,kw)
  END IF
  vmultd(k) = temp
ELSE
  vmultd(k) = 0.0D0
END IF
k = k - 1
IF (k > 0) GO TO 130
IF (ratio < 0.0D0) GO TO 490

!  Revise the Lagrange multipliers and reorder the active constraints so
!  that the one to be replaced is at the end of the list. Also calculate the
!  new value of ZDOTA(NACT) and branch if it is not acceptable.

DO k=1,nact
  vmultc(k) = MAX(0.0D0,vmultc(k) - ratio*vmultd(k))
END DO
IF (icon < nact) THEN
  isave = iact(icon)
  vsave = vmultc(icon)
  k = icon
  170 kp = k + 1
  kw = iact(kp)
  sp = DOT_PRODUCT( z(1:n,k), a(1:n,kw) )
  temp = SQRT(sp*sp + zdota(kp)**2)
  alpha = zdota(kp)/temp
  beta = sp/temp
  zdota(kp) = alpha*zdota(k)
  zdota(k) = temp
  DO i=1,n
    temp = alpha*z(i,kp) + beta*z(i,k)
    z(i,kp) = alpha*z(i,k) - beta*z(i,kp)
    z(i,k) = temp
  END DO
  iact(k) = kw
  vmultc(k) = vmultc(kp)
  k = kp
  IF (k < nact) GO TO 170
  iact(k) = isave
  vmultc(k) = vsave
END IF
temp = DOT_PRODUCT( z(1:n,nact), a(1:n,kk) )
IF (temp == 0.0D0) GO TO 490
zdota(nact) = temp
vmultc(icon) = 0.0D0
vmultc(nact) = ratio

!  Update IACT and ensure that the objective function continues to be
!  treated as the last active constraint when MCON>M.

210 iact(icon) = iact(nact)
iact(nact) = kk
IF (mcon > m .AND. kk /= mcon) THEN
  k = nact - 1
  sp = DOT_PRODUCT( z(1:n,k), a(1:n,kk) )
  temp = SQRT(sp*sp + zdota(nact)**2)
  alpha = zdota(nact)/temp
  beta = sp/temp
  zdota(nact) = alpha*zdota(k)
  zdota(k) = temp
  DO i=1,n
    temp = alpha*z(i,nact) + beta*z(i,k)
    z(i,nact) = alpha*z(i,k) - beta*z(i,nact)
    z(i,k) = temp
  END DO
  iact(nact) = iact(k)
  iact(k) = kk
  temp = vmultc(k)
  vmultc(k) = vmultc(nact)
  vmultc(nact) = temp
END IF

!  If stage one is in progress, then set SDIRN to the direction of the next
!  change to the current vector of variables.

IF (mcon > m) GO TO 320
kk = iact(nact)
temp = DOT_PRODUCT( sdirn(1:n), a(1:n,kk) )
temp = temp - 1.0D0
temp = temp/zdota(nact)
sdirn(1:n) = sdirn(1:n) - temp*z(1:n,nact)
GO TO 340

!  Delete the constraint that has the index IACT(ICON) from the active set.

260 IF (icon < nact) THEN
  isave = iact(icon)
  vsave = vmultc(icon)
  k = icon
  DO
    kp = k + 1
    kk = iact(kp)
    sp = DOT_PRODUCT( z(1:n,k), a(1:n,kk) )
    temp = SQRT(sp*sp + zdota(kp)**2)
    alpha = zdota(kp)/temp
    beta = sp/temp
    zdota(kp) = alpha*zdota(k)
    zdota(k) = temp
    DO i=1,n
      temp = alpha*z(i,kp) + beta*z(i,k)
      z(i,kp) = alpha*z(i,k) - beta*z(i,kp)
      z(i,k) = temp
    END DO
    iact(k) = kk
    vmultc(k) = vmultc(kp)
    k = kp
    IF (k >= nact) EXIT
  END DO
  iact(k) = isave
  vmultc(k) = vsave
END IF
nact = nact - 1

!  If stage one is in progress, then set SDIRN to the direction of the next
!  change to the current vector of variables.

IF (mcon > m) GO TO 320
temp = DOT_PRODUCT( sdirn(1:n), z(1:n,nact+1) )
sdirn(1:n) = sdirn(1:n) - temp*z(1:n,nact+1)
GO TO 340

!  Pick the next search direction of stage two.

320 temp = 1.0D0/zdota(nact)
sdirn(1:n) = temp*z(1:n,nact)

!  Calculate the step to the boundary of the trust region or take the step
!  that reduces RESMAX to zero. The two statements below that include the
!  factor 1.0E-6 prevent some harmless underflows that occurred in a test
!  calculation. Further, we skip the step if it could be zero within a
!  reasonable tolerance for computer rounding errors.

340 dd = rho*rho
sd = 0.0D0
ss = 0.0D0
DO i=1,n
  IF (ABS(dx(i)) >= 1.0E-6*rho) dd = dd - dx(i)**2
  sd = sd + dx(i)*sdirn(i)
  ss = ss + sdirn(i)**2
END DO
IF (dd <= 0.0D0) GO TO 490
temp = SQRT(ss*dd)
IF (ABS(sd) >= 1.0E-6*temp) temp = SQRT(ss*dd + sd*sd)
stpful = dd/(temp + sd)
step = stpful
IF (mcon == m) THEN
  acca = step + 0.1D0*resmax
  accb = step + 0.2D0*resmax
  IF (step >= acca .OR. acca >= accb) GO TO 480
  step = MIN(step,resmax)
END IF

!  Set DXNEW to the new variables if STEP is the steplength, and reduce
!  RESMAX to the corresponding maximum residual if stage one is being done.
!  Because DXNEW will be changed during the calculation of some Lagrange
!  multipliers, it will be restored to the following value later.

dxnew(1:n) = dx(1:n) + step*sdirn(1:n)
IF (mcon == m) THEN
  resold = resmax
  resmax = 0.0D0
  DO k=1,nact
    kk = iact(k)
    temp = b(kk) - DOT_PRODUCT( a(1:n,kk), dxnew(1:n) )
    resmax = MAX(resmax,temp)
  END DO
END IF

!  Set VMULTD to the VMULTC vector that would occur if DX became DXNEW. A
!  device is included to force VMULTD(K) = 0.0 if deviations from this value
!  can be attributed to computer rounding errors. First calculate the new
!  Lagrange multipliers.

k = nact
390 zdotw = 0.0D0
zdwabs = 0.0D0
DO i=1,n
  temp = z(i,k)*dxnew(i)
  zdotw = zdotw + temp
  zdwabs = zdwabs + ABS(temp)
END DO
acca = zdwabs + 0.1D0*ABS(zdotw)
accb = zdwabs + 0.2D0*ABS(zdotw)
IF (zdwabs >= acca .OR. acca >= accb) zdotw = 0.0D0
vmultd(k) = zdotw / zdota(k)
IF (k >= 2) THEN
  kk = iact(k)
  dxnew(1:n) = dxnew(1:n) - vmultd(k)*a(1:n,kk)
  k = k - 1
  GO TO 390
END IF
IF (mcon > m) vmultd(nact) = MAX(0.0D0,vmultd(nact))

!  Complete VMULTC by finding the new constraint residuals.

dxnew(1:n) = dx(1:n) + step*sdirn(1:n)
IF (mcon > nact) THEN
  kl = nact + 1
  DO k=kl,mcon
    kk = iact(k)
    total = resmax - b(kk)
    sumabs = resmax + ABS(b(kk))
    DO i=1,n
      temp = a(i,kk)*dxnew(i)
      total = total + temp
      sumabs = sumabs + ABS(temp)
    END DO
    acca = sumabs + 0.1*ABS(total)
    accb = sumabs + 0.2*ABS(total)
    IF (sumabs >= acca .OR. acca >= accb) total = 0.0
    vmultd(k) = total
  END DO
END IF

!  Calculate the fraction of the step from DX to DXNEW that will be taken.

ratio = 1.0D0
icon = 0
DO k=1,mcon
  IF (vmultd(k) < 0.0D0) THEN
    temp = vmultc(k)/(vmultc(k) - vmultd(k))
    IF (temp < ratio) THEN
      ratio = temp
      icon = k
    END IF
  END IF
END DO

!  Update DX, VMULTC and RESMAX.

temp = 1.0D0 - ratio
dx(1:n) = temp*dx(1:n) + ratio*dxnew(1:n)
DO k=1,mcon
  vmultc(k) = MAX(0.0D0,temp*vmultc(k) + ratio*vmultd(k))
END DO
IF (mcon == m) resmax = resold + ratio*(resmax - resold)

!  If the full step is not acceptable then begin another iteration.
!  Otherwise switch to stage two or end the calculation.

IF (icon > 0) GO TO 70
IF (step == stpful) GO TO 500
480 mcon = m + 1
icon = mcon
iact(mcon) = mcon
vmultc(mcon) = 0.0D0
GO TO 60

!  We employ any freedom that may be available to reduce the objective
!  function before returning a DX whose length is less than RHO.

490 IF (mcon == m) GO TO 480
ifull = 0

500 RETURN
END SUBROUTINE trstlp

END MODULE cobyla2
    
    
MODULE OPTIMIZATION_MODULE
USE COST_MODULE
USE Constrained_minimization
USE cobyla2
USE omp_lib
USE UNCMIN_MOD
IMPLICIT NONE



    CONTAINS

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS SUBROUTINE IS THE PROGRAM THAT CONTROLS THE FLOW OF THE GENETIC
!  ALGORITHM.  THIS IS THE SUBROUTINE THAT SHOULD BE CALL BE THE 
!  CALLED TO RUN THE GENETIC ALGORITHM.  
!    
!  INPUTS:
!    IPRINT        =   A VALUE GREATER THAN 0 PRINT INFORMATION FOR
!                      EACH GENERATION AS IT IS FINISH, INTEGER
!    N_POP         =   SIZE OF THE POPULATION, INTEGER
!    N_GEN         =   MAXIMUM NUMBER OF GENERATIONS TO RUN, INTEGER
!    N_INT         =   NUMBER OF INTEGER VARIABLES, INTEGER
!    N_DOUBLE      =   NUMBER OF REAL VALUED VARIABLES, INTEGER
!    N1            =   FIRST DIMENSION OF THE INPUT ARRAY, INTEGER
!    N2            =   SECOND DIMENSION OF THE INPUT ARRAY, INTEGER 
!    ITER_MAX_NLP  =   MAXIMUM NUMBER OF ITERATIONS ALLOWED FOR THE
!                      SPECIFIED NLP SOLVER, INTEGER
!    N_CON         =   NUMBER OF CONSTRAINTS, CAN ONLY BE USED WITH 
!                      COBYLA AND CONMIN, INTEGER
!    INTEGER_UPPER =   UPPER BOUNDS FOR INTEGER VARIABLES, 
!                      INTEGER(N_INT)
!    INTEGER_LOWER =   LOWER BOUNDS FOR INTEGER VARIABLES, 
!                      INTEGER(N_INT)
!    P_CROSS       =   PROBABILITY THAT A CROSSOVER WILL OCCUR, VALUES
!                      SHOULD TYPICALLY BE AROUND 0.9, 
!                      DOUBLE PRECISION
!    P_REP         =   PROBABILITY THAT REPRODUCTIONS WILL OCCUR, 
!                      VALUES SHOULD BE AROUND 0.1, DOUBLE PRECISION
!    P_MUT         =   PROBABILITY THAT A MUTATION WILL OCCUR, VALUES
!                      SHOULD BE TYPICALLY BE KEPT LOWER THAN 0.1,
!                      DOUBLE PRECISION
!    DOUBLE_UPPER  =   REAL VALUED VARIABLES UPPER BOUNDS, 
!                      DOUBLE PRECISION(N_DOUBLE)
!    DOUBLE_LOWER  =   REAL VALUED VARIABLES LOWER BOUNDS, 
!                      DOUBLE PRECISION(N_DOUBLE)
!    INPUT_ARRAY   =   INPUT ARRAY TO BE USED FOR ADDITIONAL INPUTS 
!                      THAT THE COST FUNCTION MAY NEED, 
!                      DOUBLE PRECISION (N1,N2)
!    CROSS_TYPE    =   TYPE OF CROSSOVER TO BE USED, OPTIONS ARE:
!                      UNIFORM, SINGLE_POINT, DOUBLE_POINT, 
!                      ARITHMETIC,  AND HEURISTION, 
!                      CHARACTER WITH A LENGTH OF 30
!    MUT_TYPE      =   TYPE OF MUTATION TO BE USED, OPTIONS ARE:
!                      UNIFORM, SLIDING, AND BOUNDARY, 
!                      CHARACTER WITH A LENGTH OF 30
!    SEL_TYP       =   SELECTION TYPE TO BE USE, OPTIONS ARE:
!                      ROULETTE AND TOURNAMENT, 
!                      CHARACTER WITH A LENGTH OF 30
!    OPT_TYPE      =   OPTIMIZATION TYPE TO BE USED, OPTIONS ARE:
!                      GEN, HYB_COBYLA, HYB_CONMIN, HYB_UNCMIN,
!                      CHARACTER WITH A LENGTH OF 30
!    SEED          =   SEED VALUE FOR THE RANDOM NUMBER GENERATOR, 
!                      SHOULD BE STARTED WITH A NEGATIVE INTEGER 
!                      VALUE, INTEGER
!   
!  OUTPUTS:
!    FITNESS_MIN   =   ARRAY OF MINIMUM FITNESS VALUES FOR EACH 
!                      GENERATION, DOUBLE PRECISION(N_GEN)
!    FITNESS_AVG   =   ARRAY OF THE AVERAGE FITNESS VALUES FOR EACH
!                      GENERATION, DOUBLE PRECISION(N_GEN)
!    INTEGER_MIN   =   INTEGER CHROMOSOME CORRESPONDING TO THE MINIMUM
!                      SOLUTION FOR EACH GENERATION, 
!                      INTEGER(N_GEN, N_INT)
!    DOUBLE_MIN    =   REAL VALUES CHROMOSOME CORRESPONDING TO THE 
!                      MINIMUM SOLUTION FOR EACH GENERATION, 
!                      DOUBLE PRECISION(N_GEN,N_DOUBLE)
!
!********************************************************************!
SUBROUTINE GENETIC_DRIVER(IPRINT, N_POP, N_GEN, N_INT, N_DOUBLE, N1, &
    N2, ITER_MAX_NLP, N_CON, INTEGER_UPPER, INTEGER_LOWER, P_CROSS, &
    P_REP,P_MUT, DOUBLE_UPPER, DOUBLE_LOWER, INPUT_ARRAY, CROSS_TYPE,&
    MUT_TYPE, SEL_TYPE, OPT_TYPE, SEED, FITNESS_MIN, FITNESS_AVG, &
    INTEGER_MIN, DOUBLE_MIN)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N_POP, N_GEN, N_INT, N_DOUBLE, N1, N2, N_CON, &
    INTEGER_UPPER(N_INT), INTEGER_LOWER(N_INT), IPRINT

INTEGER, INTENT(INOUT) :: SEED, ITER_MAX_NLP

DOUBLE PRECISION, INTENT(IN) :: P_CROSS, P_REP,P_MUT

DOUBLE PRECISION, INTENT(INOUT) :: DOUBLE_LOWER(N_DOUBLE), &
    DOUBLE_UPPER(N_DOUBLE), INPUT_ARRAY(N1,N2)

CHARACTER(LEN=30), INTENT(IN) :: CROSS_TYPE, MUT_TYPE, SEL_TYPE, &
    OPT_TYPE

DOUBLE PRECISION, INTENT(INOUT) :: FITNESS_MIN(N_GEN), &
    FITNESS_AVG(N_GEN), DOUBLE_MIN(N_GEN,N_DOUBLE)

INTEGER, INTENT(INOUT) :: INTEGER_MIN(N_GEN,N_INT)

DOUBLE PRECISION :: POP_DOUBLE(N_POP,N_DOUBLE), &
    POP_NEW_DOUBLE(N_POP,N_DOUBLE), FITNESS_INDV, FITNESS(N_POP), &
    FITNESS_NEW(N_POP), CHROM_DOUBLE(N_DOUBLE), TIME, &
    FITNESS_INDV_NLP, RAN, G_CON(N_CON)

INTEGER :: POP_INT(N_POP,N_INT), POP_NEW_INT(N_POP,N_INT), NAN_COUNT,&
    TEST, CHROM_INT(N_INT), P, Q, I, MIN_LOC, COUNT1, COUNT2, RATE, &
    ncon, count

DOUBLE PRECISION :: AVG, X0(N_DOUBLE), X(N_DOUBLE)

INTEGER :: NGEN_CONVERGE, INFO, NACMX1, N1_C, N2_C, N3_C, N4_C, N5_C

NCON=N_CON

CALL SYSTEM_CLOCK(COUNT1,RATE)    

NAN_COUNT=0
NGEN_CONVERGE=50

!GENERATE INITIAL POPULATION AND FIND THEIR FITNESS VALUES
CALL POPULATION_GENERATOR(N_POP, N_DOUBLE, N_INT, POP_DOUBLE, &
    POP_INT, SEED, INTEGER_UPPER, INTEGER_LOWER, DOUBLE_UPPER, &
    DOUBLE_LOWER)
   
!$OMP PARALLEL DO private(CHROM_INT, CHROM_DOUBLE, FITNESS_INDV, X0)&
!$OMP& PRIVATE(X, G_CON, INPUT_ARRAY)   

DO P=1,N_POP,1
    CHROM_DOUBLE=POP_DOUBLE(P,1:N_DOUBLE)
    CHROM_INT=POP_INT(P,1:N_INT)
    CALL COST(N_DOUBLE, N_INT, N1, N2, CHROM_DOUBLE, CHROM_INT, &
              FITNESS_INDV, INPUT_ARRAY, G_CON, NCON)
    ! THE INITIAL POPULATION SHOULD BE CHECKED TO ENSURE THAT THE COST
    ! FUNCTION DOESN'T GIVE ANY NANS OR VALUES THAT ARE TOO LARGER IF 
    ! IT DOES A NEW CHROMOSOME IS RANDOMLY GENERATED UNTIL YOU FIND 
    ! ONE THAT DOESN'T GIVE A NAN/VERY LARGE VALUES.
    IF (FITNESS_INDV .GE. 1.D24 .OR. ISNAN(FITNESS_INDV)) THEN
        COUNT=0
        DO WHILE(FITNESS_INDV .GE. 1.D24 .AND. COUNT.LT.50000)
            COUNT=COUNT+1
            CALL CHROMOSOME_GENERATOR(SEED, N_INT, N_DOUBLE, &
                CHROM_INT, CHROM_DOUBLE, INTEGER_UPPER, &
                INTEGER_LOWER, DOUBLE_LOWER, DOUBLE_UPPER)
        
            CALL COST(N_DOUBLE, N_INT, N1, N2, CHROM_DOUBLE, &
                CHROM_INT, FITNESS_INDV, INPUT_ARRAY, G_CON, NCON)           
        END DO
        IF (COUNT.EQ.50000) THEN
            WRITE(*,*) "UNABLE TO FIND CHROMOSOME" 
            WRITE(*,*) "THE PROBLEM MAY NOT BE WELL POSED"
        END IF
        POP_DOUBLE(P,1:N_DOUBLE)=CHROM_DOUBLE
        POP_INT(P,1:N_INT)=CHROM_INT
    END IF

    FITNESS(P)=FITNESS_INDV

END DO

!$OMP END PARALLEL DO
MIN_LOC=MINLOC(FITNESS,1)
FITNESS_MIN(1)=FITNESS(MIN_LOC)
INTEGER_MIN(1,1:N_INT)=POP_INT(MIN_LOC,1:N_INT)
DOUBLE_MIN(1,1:N_DOUBLE)=POP_DOUBLE(MIN_LOC,1:N_DOUBLE)
FITNESS_AVG(1)=SUM(FITNESS)/DBLE(N_POP)


CALL SYSTEM_CLOCK(COUNT2)

TIME=DBLE(COUNT2-COUNT1)/DBLE(RATE)
IF (IPRINT.GT.0) THEN
WRITE(*,*) ' GENERATION','   MINIMUM COST', &
    '              AVERAGE COST', '            NAN-COUNT     ',  &
    'RUN TIME'
WRITE(*,*) 1, FITNESS_MIN(1), FITNESS_AVG(1), NAN_COUNT, TIME
END IF

DO Q=2,N_GEN,1
		    CALL GENETIC_OPERATIONS(FITNESS,FITNESS_NEW, POP_INT,&
                POP_NEW_INT, POP_DOUBLE, POP_NEW_DOUBLE, P_REP, &
                P_CROSS, P_MUT, N_POP, CROSS_TYPE, MUT_TYPE, &
                SEL_TYPE, N_INT, N_DOUBLE, SEED, INTEGER_UPPER, &
                INTEGER_LOWER, DOUBLE_UPPER, DOUBLE_LOWER)
  
	FITNESS=FITNESS_NEW
    ! IF THE CONMIN SOVLER IS BEING USES, EVEN THE TOP TWO SOLUTIONS 
    ! SHOULD BE CHECKED TO ENSURE THAT CONSTRAINTS ARE MINIMIZED AS  
    ! MUCH AS POSSIBLE
    IF(TRIM(OPT_TYPE).EQ."HYB_CONMIN") THEN
        FITNESS(1)=1.D8
        FITNESS(2)=1.D8
    END IF

    !$OMP PARALLEL DO PRIVATE(P, CHROM_INT, CHROM_DOUBLE)&
    !$OMP& PRIVATE(FITNESS_INDV, X0, X, FITNESS_INDV_NLP, G_CON)&
    !$OMP& PRIVATE(INPUT_ARRAY)
    DO P=1, N_POP,1

	    !THE COST FUNCTION DOES NOT NEED TO BE CALCULATED FOR 
        !INDIVIDUALS THAT WENT THROUGH REPRODUCTION, BUT NOT MUTATION.
        !IF AN INDIVIDUAL IN THE POPULATION IS FROM A  CROSSOVER  OR 
        !MUTATION THE FITNESS VALUE IS SET TO 1.D8 AND A NEW FITNESS 
        !VALUE SHOULD BE COMPUTED.  THIS ENSURE THAT NO UNNECESSARY 
        !COST FUNCTION EVALUATIONS ARE PERFORMED.  THIS IS ALSO THE 
        !PORTION  OF  THE  CODE  THAT SHOULD BE EXECUTED IN PARALLEL.  
        !ALL OTHER PARTS  OF THE GENETIC  ALGORITHM TYPICALLY  REQUIRE 
        !VERY LITTLE TIME TO EXECUTE COMPARED TO THE EVALUATION OF THE
        !FITNESS FUNCTION.
		IF (ABS(FITNESS(P)-1.D8).LT.1.D0)THEN
            !MAKE SURE THAT THE NEW SOLUTION REGION IS FEASIBLE.  THIS 
            !SECTION ALSO PERFORMS ALL THE NECESSARY CALCULATIONS FOR  
            !THE PURE GENETIC ALGORITHM.   IF  ONE  OF  THE  HYBRID  
            !METHODS IS UTILIZED FURTHER CALCULATIONS WILL THEN BE 
            !PERFORMED.
            
            CHROM_INT=POP_NEW_INT(P,1:N_INT)
			CHROM_DOUBLE=POP_NEW_DOUBLE(P,1:N_DOUBLE)
            CALL COST(N_DOUBLE, N_INT, N1, N2, CHROM_DOUBLE, &
                CHROM_INT, FITNESS_INDV, INPUT_ARRAY, G_CON, NCON)
	            
            !NOW TEST TO MAKE SURE THE SOLUTION REGION HAS IS FEASIBLE  
            !SOLUTION AT THE POINT IN THE POPULATION.  IF NOT THE  
            !ORIGINAL SOLUTION FROM THE PRECIOUS POPULATION MEMBER P 
            !WILL BE USED INSTEAD.
            IF (FITNESS_INDV.GE.1.D24 .OR. ISNAN(FITNESS_INDV)) THEN
                CHROM_INT=POP_INT(P,1:N_INT)
                CHROM_DOUBLE=POP_DOUBLE(P,1:N_DOUBLE)
 
                CALL COST(N_DOUBLE, N_INT, N1, N2, CHROM_DOUBLE, &
                    CHROM_INT, FITNESS_INDV, INPUT_ARRAY, G_CON, NCON)
                POP_NEW_INT(P,1:N_INT)=CHROM_INT
                POP_NEW_DOUBLE(P,1:N_INT)=CHROM_DOUBLE
            END IF
            
            IF(TRIM(OPT_TYPE).EQ."GEN") THEN
                !DO NOTHING, THE PURE GENETIC OPERATIONS HAVE ALREADY  
                !BEEN PERFORMED
                X=CHROM_DOUBLE
          
            ELSE IF(TRIM(OPT_TYPE).EQ."HYB_UNCMIN") THEN
                !HYBRID ALGORITHM THAT USES AN UNCONSTRAINED  
                !MINIMIZATION NLP SOVLER TO ITERATE ON THE REAL VALUED
                !CHROMOSOME
            
                X0=CHROM_DOUBLE
                CALL UNCMIN_WRAPPER(N_DOUBLE, N_INT, N1, N2, &
                    ITER_MAX_NLP, X0, X, CHROM_INT, FITNESS_INDV_NLP,&
                    INPUT_ARRAY, G_CON, NCON)
                !MAKE SURE THAT  THE  NEW  SOLUTION  IS FEASIBLE.                
                !THIS SECTION CHECK  IF  THE  UNCMIN  CONVERGED  ON  
                !A  SOLUTION (HOPEFULLY OPTIMAL, BUT  IT MAY  NOT  BE,   
                !DON'T  WORRY  THE IDEA IS THE GENETIC  ALGORITHM  
                !WILL  EVENTUALLY  CONVERGE ON AND OPTIMAL SOLUTION,  
                !NOT  INDIVIDUAL  RUNS OF  THE NLP  SOLVER).  IF IT
                !DIDN'T THE COST FUNCTION SHOULD BE DESIGNED TO OUTPUT
                !A VALUELARGER THAN 1.D8.  THE ALGORITHM  ALSO CHECKS 
                !FOR  NANS,  BUT THE  USER  DESIGNED  COST  FUNCTION  
                !SHOULD BE DESIGNED NOT TO ALLOW ANY NANS (THIS WILL  
                !LIKELY BREAK NEARLY ALL NLP SOLVERS WITH EVEN 1 NAN).

                DO I=1,N_DOUBLE,1
                    IF(X(I) .GT. DOUBLE_UPPER(I) .OR. &
                        X(I).LT.DOUBLE_LOWER(I)) THEN
                        FITNESS_INDV_NLP=1.D30    
                    END IF
                END DO   
                
                IF(FITNESS_INDV_NLP .GE. 1.D24 .OR. &
                    ISNAN(FITNESS_INDV_NLP))THEN 
                    !UNCMIN FAILED, SO THE SOLUTION FROM THE PREVIOUS 
                    !GENERATION WILL BE USED INSTEAD.
                    CHROM_INT=POP_INT(P,1:N_INT)
                    CHROM_DOUBLE=POP_DOUBLE(P,1:N_DOUBLE)
                    POP_NEW_DOUBLE(P,1:N_DOUBLE)=CHROM_DOUBLE
                    POP_NEW_INT(P,1:N_INT)=CHROM_INT
                ELSE
        
                    !THE SOLUTION FOUND BY UNCMIN IS VALID, SO THE 
                    !NEW POPULATION  WITH  THE  DETERMINED SOLUTION 
                    !AND COST FUNCTION IS UPDATED.
                    POP_NEW_DOUBLE(P,1:N_DOUBLE)=X
                    FITNESS_INDV=FITNESS_INDV_NLP
                END IF    
            ELSE IF(TRIM(OPT_TYPE).EQ."HYB_CONMIN") THEN
                !HYBRID ALGORITHM THAT USES A CONMIN
               
                X0=CHROM_DOUBLE

                CALL CONMIN_WRAPPER(N_DOUBLE, N_INT, N1, N2, &
                    ITER_MAX_NLP, X0, X, CHROM_INT, FITNESS_INDV_NLP,&
                    INPUT_ARRAY, N_CON, DOUBLE_UPPER, DOUBLE_LOWER)  

                IF(FITNESS_INDV_NLP .GE. 1.D24 .OR. &
                    ISNAN(FITNESS_INDv_NLP))THEN 
                    !CONMIN FAILED, SO THE PREVIOUS GENERATION WILL 
                    !BE USED INSTEAD.
                    CHROM_INT=POP_INT(P,1:N_INT)
                    CHROM_DOUBLE=POP_DOUBLE(P,1:N_DOUBLE)
                    POP_NEW_DOUBLE(P,1:N_DOUBLE)=CHROM_DOUBLE
                    POP_NEW_INT(P,1:N_INT)=CHROM_INT
                ELSE
          
                    !THE SOLUTION FOUND BY CONMIN IS VALID, SO THE NEW 
                    !POPULATION  WITH  THE  DETERMINED SOLUTION AND  
                    !COST FUNCTION IS UPDATED.
                    POP_NEW_DOUBLE(P,1:N_DOUBLE)=X
                    FITNESS_INDV=FITNESS_INDV_NLP
                END IF
            ELSE IF(TRIM(OPT_TYPE).EQ."HYB_COBYLA") THEN
                !HYBRID ALGORITHM THAT USES A CONSTRAINED MINIMIZATION
                !WRITE(*,*) "CALLING CONMIN WRAPPER"

                X=CHROM_DOUBLE
                CALL COBYLA_DRIVER(N_DOUBLE, N_INT, N1, N2, X, &
                    CHROM_INT, FITNESS_INDV_NLP, INPUT_ARRAY, N_CON,&
                    ITER_MAX_NLP )
                
                DO I=1,N_DOUBLE,1
                    IF(X(I) .GT. DOUBLE_UPPER(I) .OR. &
                        X(I).LT.DOUBLE_LOWER(I)) THEN
                        FITNESS_INDV_NLP=1.D30    
                    END IF
                END DO  
                
                IF(FITNESS_INDV_NLP .GE. 1.D24 .OR. &
                    ISNAN(FITNESS_INDv_NLP))THEN 
                    !COBYLA.
                    CHROM_INT=POP_INT(P,1:N_INT)
                    CHROM_DOUBLE=POP_DOUBLE(P,1:N_DOUBLE)
                    POP_NEW_DOUBLE(P,1:N_DOUBLE)=CHROM_DOUBLE
                    POP_NEW_INT(P,1:N_INT)=CHROM_INT
                ELSE
       
                    !THE SOLUTION FOUND BY CONMIN IS VALID, SO THE NEW 
                    !POPULATION  WITH  THE  DETERMINED SOLUTION AND  
                    !COST FUNCTION IS UPDATED.
                    POP_NEW_DOUBLE(P,1:N_DOUBLE)=X
                    FITNESS_INDV=FITNESS_INDV_NLP
                END IF
            END IF


            
			FITNESS(P)=FITNESS_INDV
        END IF
	END DO
    !$OMP END PARALLEL DO


	POP_INT=POP_NEW_INT
	POP_DOUBLE=POP_NEW_DOUBLE
	MIN_LOC=MINLOC(FITNESS,1)
	FITNESS_MIN(Q)=FITNESS(MIN_LOC)

    !CHECK TO MAKE SURE THE SOLUTION HASN'T STAGNATED FOR MORE THAN 
    !25 GENERATIONS
    AVG=1.D0
    IF (Q>NGEN_CONVERGE) THEN
        AVG=ABS(SUM(FITNESS_MIN(Q-NGEN_CONVERGE:Q))-&
            SUM(FITNESS_MIN(Q-(NGEN_CONVERGE+1):Q-1)))
        IF (AVG.LE.1.D-1) THEN
            DO I=Q,N_GEN,1
               INTEGER_MIN(I,1:N_INT)=POP_INT(MIN_LOC,1:N_INT)
               DOUBLE_MIN(I,1:N_DOUBLE)=POP_DOUBLE(MIN_LOC,1:N_DOUBLE)
               FITNESS_MIN(I)=FITNESS(MIN_LOC)
            END DO
            EXIT
        END IF
        
    END IF
    
	INTEGER_MIN(Q,1:N_INT)=POP_INT(MIN_LOC,1:N_INT)
	DOUBLE_MIN(Q,1:N_DOUBLE)=POP_DOUBLE(MIN_LOC,1:N_DOUBLE)
	FITNESS_AVG(Q)=SUM(FITNESS)/DBLE(N_POP)
    CALL SYSTEM_CLOCK(COUNT2)

    TIME=DBLE(COUNT2-COUNT1)/DBLE(RATE)
    IF (IPRINT.GT.0) THEN
	    WRITE(*,*) Q, FITNESS_MIN(Q), FITNESS_AVG(Q), NAN_COUNT, TIME
    END IF
END DO

END SUBROUTINE GENETIC_DRIVER

!********************************************************************! 
!********************************************************************! 
!********************************************************************! 
!********************************************************************! 
! THIS SUBROUTINE PERFORMS ALL THE GENETIC OPERATIONS ON THE INITIAL 
!  POPULATION FOR EACH GENERATION
SUBROUTINE GENETIC_OPERATIONS(FITNESS,FITNESS_NEW, POP_INT,&
    POP_NEW_INT, POP_DOUBLE, POP_NEW_DOUBLE, P_REP, P_CROSS, P_MUT, &
    N_POP, CROSS_TYPE,  MUT_TYPE, SEL_TYPE, N_INT, N_DOUBLE, SEED, &
    INTEGER_UPPER, INTEGER_LOWER, DOUBLE_UPPER, DOUBLE_LOWER)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N_POP, N_INT, N_DOUBLE, POP_INT(N_POP,N_INT), &
    INTEGER_UPPER(N_INT), INTEGER_LOWER(N_INT)
INTEGER, INTENT(INOUT) :: SEED, POP_NEW_INT(N_POP,N_INT)
DOUBLE PRECISION, INTENT(IN) :: FITNESS(N_POP), P_REP, P_CROSS, &
    POP_DOUBLE(N_POP,N_DOUBLE), P_MUT, DOUBLE_UPPER(N_DOUBLE), &
    DOUBLE_LOWER(N_DOUBLE)
DOUBLE PRECISION, INTENT(INOUT) :: FITNESS_NEW(N_POP), &
    POP_NEW_DOUBLE(N_POP,N_DOUBLE)
CHARACTER(LEN=30), INTENT(IN) :: CROSS_TYPE, MUT_TYPE, SEL_TYPE
DOUBLE PRECISION :: RAN, CHROM_DOUBLE(N_DOUBLE)
INTEGER :: P, CHROM_INT(N_INT)

IF(TRIM(SEL_TYPE).EQ."ROULETTE")THEN
	CALL ROULETTE(SEED, N_POP, N_INT, N_DOUBLE, FITNESS, FITNESS_NEW,&
        POP_INT, POP_NEW_INT, POP_DOUBLE, POP_NEW_DOUBLE, P_CROSS, &
        P_REP, CROSS_TYPE, INTEGER_UPPER, INTEGER_LOWER, &
        DOUBLE_UPPER, DOUBLE_LOWER)
ELSE IF(TRIM(SEL_TYPE).EQ."TOURNAMENT")THEN
	CALL TOURNAMENT(SEED, N_POP, N_INT, N_DOUBLE, FITNESS, &
        FITNESS_NEW, POP_INT, POP_NEW_INT, POP_DOUBLE, &
        POP_NEW_DOUBLE, P_CROSS, P_REP, CROSS_TYPE, INTEGER_UPPER, &
        INTEGER_LOWER, DOUBLE_UPPER, DOUBLE_LOWER)
ELSE
	!WRITE(*,*) "INVALID SELECTION TYPE"
END IF

! AFTER NEW POPULATION HAS BEEN GENERATED THE MUTATIONS ARE DONE. THE  
! MUTATION PROBABILITY SHOULD BE USED RELATIVELY LOW TO AVOID A 
! COMPLETELY RANDOM SEARCH
DO P=3,N_POP,1
	RAN=RANDOM(SEED)
	IF (RAN.LE.P_MUT) THEN
		CHROM_INT=POP_INT(P,1:N_INT)
		CHROM_DOUBLE=POP_DOUBLE(P,1:N_DOUBLE)
		CALL MUTATION(CHROM_INT, CHROM_DOUBLE, N_INT, N_DOUBLE, SEED, &
            MUT_TYPE, INTEGER_LOWER, INTEGER_UPPER, DOUBLE_LOWER, &
            DOUBLE_UPPER)
		POP_NEW_INT(P,1:N_INT)=CHROM_INT
		POP_NEW_DOUBLE(P,1:N_DOUBLE)=CHROM_DOUBLE
		FITNESS_NEW(P)=1.D8
	END IF
END DO

END SUBROUTINE GENETIC_OPERATIONS

!********************************************************************! 
!********************************************************************! 
!********************************************************************! 
!********************************************************************! 
! THIS SUBROUTINE PERFORMES THE ROULETTE SELECTION METHOD OPERATOR
!
SUBROUTINE ROULETTE(SEED, N_POP, N_INT, N_DOUBLE, FITNESS, &
    FITNESS_NEW, POP_INT, POP_NEW_INT, POP_DOUBLE, POP_NEW_DOUBLE, &
    P_CROSS, P_REP, CROSS_TYPE, INTEGER_UPPER, INTEGER_LOWER, &
    DOUBLE_UPPER, DOUBLE_LOWER)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N_POP, N_INT, N_DOUBLE, POP_INT(N_POP, N_INT)

INTEGER, INTENT(INOUT) :: SEED, POP_NEW_INT(N_POP, N_INT)

DOUBLE PRECISION, INTENT(IN) :: P_CROSS, P_REP, FITNESS(N_POP), &
    POP_DOUBLE(N_POP, N_DOUBLE)

DOUBLE PRECISION, INTENT(INOUT) :: FITNESS_NEW(N_POP), &
    POP_NEW_DOUBLE(N_POP,N_DOUBLE)

CHARACTER(LEN=30), INTENT(IN) :: CROSS_TYPE

INTEGER, INTENT(IN) :: INTEGER_UPPER(N_INT), INTEGER_LOWER(N_INT)

DOUBLE PRECISION, INTENT(IN) :: DOUBLE_UPPER(N_DOUBLE), &
    DOUBLE_LOWER(N_DOUBLE)

DOUBLE PRECISION :: RAN, RAN1, RAN2, SUM, NORM_SWAP(N_POP), &
    FITNESS_SWAP(N_POP), NORMAL(N_POP), STD(N_POP), ADJ(N_POP), &
    SUM_ADJ, CHROM1_DOUBLE(N_DOUBLE), CHROM2_DOUBLE(N_DOUBLE), &
    POP_SWAP_DOUBLE(N_POP,N_DOUBLE), FIT1, FIT2
INTEGER :: TEST1, TEST2, ITER, CHROM1_INT(N_INT), CHROM2_INT(N_INT),&
    INDEX_SORTED(N_POP), I, POP_SWAP_INT(N_POP,N_INT), N, INDEX1, &
    INDEX2
! COMPUTE THE NORMALIZED DISTRIBUTION HERE
SUM_ADJ=0.D0
DO I=1,N_POP,1
    STD(I)=FITNESS(I)
    ADJ(I)=1.D0/(1.D0+STD(I))
    SUM_ADJ=SUM_ADJ+ADJ(I)
END DO

DO I=1,N_POP,1
    NORMAL(I)=ADJ(I)/SUM_ADJ
	INDEX_SORTED(I)=I
END DO

CALL HEAP_SORT(NORMAL,INDEX_SORTED,N_POP)

! THE NORM VECTOR IS RETURNED FROM SMALLEST TO LARGEST, BUT NORM NEED 
! TO BE SORTED FROMLARGEST TO SMALLEST, SO IT IS REVERSED HERE, ALONG 
!  WITH THE POPULATION INTO POP_SWAP
N=N_POP
DO I=1,N_POP,1
  NORM_SWAP(I)=NORMAL(N)
  POP_SWAP_DOUBLE(I,1:N_DOUBLE)=POP_DOUBLE(INDEX_SORTED(N),1:N_DOUBLE)
  POP_SWAP_INT(I,1:N_INT)=POP_INT(INDEX_SORTED(N),1:N_INT)
  FITNESS_SWAP(I)=FITNESS(INDEX_SORTED(N))
  N=N-1
END DO
!WRITE(*,*) "SWAPPING DONE"
NORMAL=NORM_SWAP
POP_NEW_INT(1,1:N_INT)=POP_SWAP_INT(1,1:N_INT)
POP_NEW_INT(2,1:N_INT)=POP_SWAP_INT(2,1:N_INT)
POP_NEW_DOUBLE(1,1:N_DOUBLE)=POP_SWAP_DOUBLE(1,1:N_DOUBLE)
POP_NEW_DOUBLE(2,1:N_DOUBLE)=POP_SWAP_DOUBLE(2,1:N_DOUBLE)


FITNESS_NEW(1)=FITNESS_SWAP(1)
FITNESS_NEW(2)=FITNESS_SWAP(2)

!WRITE(*,*) "STARTING ROULETTE SELECTION"

DO I=3,N_POP,2
    !WRITE(*,*) "I=", I

	!ROULETTE SELECTION IS PERFORMED HERE
	RAN1=RANDOM(SEED)
	RAN2=RANDOM(SEED)
	SUM=0.D0
	ITER=0
	TEST1=0
	TEST2=0
	DO WHILE (TEST1.EQ.0 .OR. TEST2.EQ.0)
		ITER=ITER+1
		SUM=SUM+NORMAL(ITER)
		IF (SUM.GT.RAN1 .AND. TEST1.EQ.0) THEN
			INDEX1=ITER
			TEST1=1
		END IF
		IF (SUM.GT.RAN2 .AND. TEST2.EQ.0) THEN
			INDEX2=ITER
			TEST2=1
		END IF
    END DO
   
    !WRITE(*,*) " REPRODUCTION AND CROSSOVER STARTED"
	!REPRODUCTION AND CROSSOVER CAN BE PERFORMED NOW THAT 2 PARENTS 
    !HAVE BEEN CHOSEN
	RAN=RANDOM(SEED)
	IF (RAN>P_CROSS)THEN
		POP_NEW_INT(I,1:N_INT)=POP_SWAP_INT(INDEX1,1:N_INT)
		POP_NEW_INT(I+1,1:N_INT)=POP_SWAP_INT(INDEX2,1:N_INT)
		POP_NEW_DOUBLE(I,1:N_DOUBLE)=&
            POP_SWAP_DOUBLE(INDEX1,1:N_DOUBLE)
		POP_NEW_DOUBLE(I+1,1:N_DOUBLE)=&
            POP_SWAP_DOUBLE(INDEX2,1:N_DOUBLE)
        FITNESS_NEW(I)=FITNESS_SWAP(INDEX1)
        FITNESS_NEW(I+1)=FITNESS_SWAP(INDEX2)
	ELSE
		CHROM1_INT=POP_SWAP_INT(INDEX1,1:N_INT)
		CHROM2_INT=POP_SWAP_INT(INDEX2,1:N_INT)
		CHROM1_DOUBLE=POP_SWAP_DOUBLE(INDEX1,1:N_DOUBLE)
		CHROM2_DOUBLE=POP_SWAP_DOUBLE(INDEX2,1:N_DOUBLE)
		FIT1=NORMAL(INDEX1)
		FIT2=NORMAL(INDEX2)

        CALL CROSSOVER(CHROM1_INT, CHROM2_INT, CHROM1_DOUBLE, &
            CHROM2_DOUBLE, N_INT, N_DOUBLE, SEED, CROSS_TYPE, &
            INTEGER_UPPER, INTEGER_LOWER, DOUBLE_UPPER, DOUBLE_LOWER,&
            FIT1, FIT2)
            
		POP_NEW_INT(I,1:N_INT)=CHROM1_INT
		POP_NEW_INT(I+1,1:N_INT)=CHROM2_INT
		POP_NEW_DOUBLE(I,1:N_DOUBLE)=CHROM1_DOUBLE
		POP_NEW_DOUBLE(I+1,1:N_DOUBLE)=CHROM2_DOUBLE
		FITNESS_NEW(I)=1.D8
		FITNESS_NEW(I+1)=1.D8
	END IF

END DO

END SUBROUTINE ROULETTE

!********************************************************************!  
!********************************************************************!  
!********************************************************************!  
!********************************************************************!  
! THIS SUBROUTINE EXECUTES THE TOURNAMENT SELECTION OPERATOR
!
SUBROUTINE TOURNAMENT(SEED, N_POP, N_INT, N_DOUBLE, FITNESS, &
    FITNESS_NEW, POP_INT, POP_NEW_INT, POP_DOUBLE, POP_NEW_DOUBLE, &
    P_CROSS, P_REP, CROSS_TYPE,  INTEGER_UPPER, INTEGER_LOWER, &
    DOUBLE_UPPER, DOUBLE_LOWER)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N_POP, N_INT, N_DOUBLE, POP_INT(N_POP, N_INT),&
    INTEGER_UPPER(N_INT), INTEGER_LOWER(N_INT)

INTEGER, INTENT(INOUT) :: SEED, POP_NEW_INT(N_POP, N_INT)

DOUBLE PRECISION, INTENT(IN) :: P_CROSS, P_REP, FITNESS(N_POP), &
    POP_DOUBLE(N_POP, N_DOUBLE), DOUBLE_UPPER(N_DOUBLE), &
    DOUBLE_LOWER(N_DOUBLE)

DOUBLE PRECISION, INTENT(INOUT) :: FITNESS_NEW(N_POP), &
    POP_NEW_DOUBLE(N_POP,N_DOUBLE)

CHARACTER(LEN=30), INTENT(IN) :: CROSS_TYPE

DOUBLE PRECISION :: NORMAL(N_POP), STD(N_POP), ADJ(N_POP), SUM_ADJ, &
    RAN, POP_DOUBLE_POOL(N_POP,N_DOUBLE), FITNESS_POOL(N_POP)

INTEGER :: I, N, INDEX, ORDER1(N_POP), ORDER2(N_POP), &
    LOCATION(N_POP), LOCATION_SWAP(N_POP), POP_INT_POOL(N_POP,N_INT)

 

! COMPUTE THE NORMALIZED DISTRIBUTION HERE
	SUM_ADJ=0.D0
	DO I=1,N_POP,1
		STD(I)=FITNESS(I)
		ADJ(I)=1.D0/(1.D0+STD(I))
		SUM_ADJ=SUM_ADJ+ADJ(I)
	END DO

	DO I=1,N_POP,1
	    NORMAL(I)=ADJ(I)/SUM_ADJ
	END DO

! FIRST RANDOM ORDER LOCATION
	DO I=1,N_POP,1
		LOCATION(I)=I
	END DO

	LOCATION_SWAP=LOCATION
	N=N_POP

	DO I=1,N_POP,1
		RAN=RANDOM(SEED)
		INDEX=FLOOR(RAN*DBLE(N)+1.D0)
		ORDER1(I)=LOCATION(INDEX)
		LOCATION_SWAP(1:INDEX-1)=LOCATION(1:INDEX-1)
		LOCATION_SWAP(INDEX:N-1)=LOCATION(INDEX+1:N)
		LOCATION=LOCATION_SWAP
		N=N-1
	END DO

! SECOND RANDOM ORDER LOCATION
	DO I=1,N_POP,1
		LOCATION(I)=I
	END DO

	LOCATION_SWAP=LOCATION
	N=N_POP

	DO I=1,N_POP,1
		RAN=RANDOM(SEED)
		INDEX=FLOOR(RAN*DBLE(N)+1.D0)
		ORDER2(I)=LOCATION(INDEX)
		LOCATION_SWAP(1:INDEX-1)=LOCATION(1:INDEX-1)
		LOCATION_SWAP(INDEX:N-1)=LOCATION(INDEX+1:N)
		LOCATION=LOCATION_SWAP
		N=N-1
	END DO
 
!PERFORM TOURNAMENT PRIOR TO ROULETTE SELECTION
	DO I=1,N_POP,2
		IF(NORMAL(ORDER1(I)).GT.NORMAL(ORDER1(I+1)))THEN
			POP_INT_POOL(I,1:N_INT)=POP_INT(ORDER1(I),1:N_INT)
			POP_DOUBLE_POOL(I,1:N_DOUBLE)=&
                POP_DOUBLE(ORDER1(I),1:N_DOUBLE)
			FITNESS_POOL(I)=FITNESS(ORDER1(I))
		ELSE
			POP_INT_POOL(I,1:N_INT)=POP_INT(ORDER1(I+1),1:N_INT)
			POP_DOUBLE_POOL(I,1:N_DOUBLE)=&
                POP_DOUBLE(ORDER1(I+1),1:N_DOUBLE)
			FITNESS_POOL(I)=FITNESS(ORDER1(I+1))
		END IF

		IF(NORMAL(ORDER2(I)).GT.NORMAL(ORDER2(I+1)))THEN
			POP_INT_POOL(I+1,1:N_INT)=POP_INT(ORDER2(I),1:N_INT)
			POP_DOUBLE_POOL(I+1,1:N_DOUBLE)=&
                POP_DOUBLE(ORDER2(I),1:N_DOUBLE)
			FITNESS_POOL(I+1)=FITNESS(ORDER2(I))
		ELSE
			POP_INT_POOL(I+1,1:N_INT)=POP_INT(ORDER2(I+1),1:N_INT)
			POP_DOUBLE_POOL(I+1,1:N_DOUBLE)=&
                POP_DOUBLE(ORDER2(I+1),1:N_DOUBLE)
			FITNESS_POOL(I+1)=FITNESS(ORDER2(I+1))
		END IF
	END DO

CALL ROULETTE(SEED, N_POP, N_INT, N_DOUBLE, FITNESS_POOL, &
    FITNESS_NEW, POP_INT_POOL, POP_NEW_INT, POP_DOUBLE_POOL, &
    POP_NEW_DOUBLE, P_CROSS, P_REP, CROSS_TYPE, INTEGER_UPPER, &
    INTEGER_LOWER, DOUBLE_UPPER, DOUBLE_LOWER)

END SUBROUTINE TOURNAMENT



!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS SUBROUTINE EXECUTES THE CROSSOVER OPERATORS
!
SUBROUTINE CROSSOVER(CHROM1_INT, CHROM2_INT, CHROM1_DOUBLE, &
    CHROM2_DOUBLE, N_INT, N_DOUBLE, SEED, CROSS_TYPE, INTEGER_UPPER,&
    INTEGER_LOWER, DOUBLE_UPPER, DOUBLE_LOWER, FIT1, FIT2)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N_INT, N_DOUBLE, INTEGER_UPPER(N_INT), &
    INTEGER_LOWER(N_INT)

INTEGER, INTENT(INOUT) :: SEED, CHROM1_INT(N_INT), CHROM2_INT(N_INT)

CHARACTER(LEN=30), INTENT(IN) :: CROSS_TYPE

DOUBLE PRECISION, INTENT(INOUT) :: CHROM1_DOUBLE(N_DOUBLE), &
    CHROM2_DOUBLE(N_DOUBLE)

DOUBLE PRECISION, INTENT(IN) :: DOUBLE_UPPER(N_DOUBLE), &
    DOUBLE_LOWER(N_DOUBLE), FIT1, FIT2

INTEGER :: I_DBLE, I, N, CHROM1_INT_NEW(N_INT), CHROM2_INT_NEW(N_INT),&
    LOC1, LOC2, LOC_SWAP, CHROM_INT_BEST(N_INT), &
    CHROM_INT_WORST(N_INT), TEST1, TEST2, COUNT

DOUBLE PRECISION :: RAN, CHROM1_DOUBLE_NEW(N_DOUBLE), &
    CHROM2_DOUBLE_NEW(N_DOUBLE), CHROM_DOUBLE_BEST(N_DOUBLE), &
    CHROM_DOUBLE_WORST(N_DOUBLE)

 

N=N_INT+N_DOUBLE
I_DBLE=0

IF(TRIM(CROSS_TYPE).EQ."UNIFORM")THEN
	DO I=1,N_INT,1
		IF(CHROM1_INT(I).EQ.CHROM2_INT(I))THEN 
            CHROM1_INT_NEW(I)=CHROM1_INT(I)
			CHROM2_INT_NEW(I)=CHROM2_INT(I)
        ELSE
			RAN=RANDOM(SEED)
			IF(RAN.LE.0.5D0)THEN
				CHROM1_INT_NEW(I)=CHROM1_INT(I)
				CHROM2_INT_NEW(I)=CHROM2_INT(I)
			ELSE
				CHROM1_INT_NEW(I)=CHROM2_INT(I)
				CHROM2_INT_NEW(I)=CHROM1_INT(I)
			END IF
		END IF
    END DO
    
    DO I=1,N_DOUBLE,1
		RAN=RANDOM(SEED)
        IF(RAN.LE.0.5D0)THEN
			CHROM1_DOUBLE_NEW(I)=CHROM1_DOUBLE(I)
			CHROM2_DOUBLE_NEW(I)=CHROM2_DOUBLE(I)				
		ELSE
			CHROM1_DOUBLE_NEW(I)=CHROM2_DOUBLE(I)
			CHROM2_DOUBLE_NEW(I)=CHROM1_DOUBLE(I)		
		END IF
	END DO

	CHROM1_INT=CHROM1_INT_NEW 
	CHROM2_INT=CHROM2_INT_NEW
	CHROM1_DOUBLE=CHROM1_DOUBLE_NEW
	CHROM2_DOUBLE=CHROM2_DOUBLE_NEW

ELSE IF(TRIM(CROSS_TYPE).EQ."SINGLE_POINT")THEN

	LOC1=RANDOM_INTEGER(SEED,N,0)

    IF (LOC1.EQ.0) LOC1=1

	DO I=1,N,1
		IF(I.LE.LOC1)THEN

			IF(I.LE.N_INT)THEN
				CHROM1_INT_NEW(I)=CHROM1_INT(I)
				CHROM2_INT_NEW(I)=CHROM2_INT(I)
			ELSE
				CHROM1_DOUBLE_NEW(I-N_INT)=CHROM1_DOUBLE(I-N_INT)
				CHROM2_DOUBLE_NEW(I-N_INT)=CHROM2_DOUBLE(I-N_INT)
			END IF

		ELSE

			IF(I.LE.N_INT)THEN
				CHROM1_INT_NEW(I)=CHROM2_INT(I)
				CHROM2_INT_NEW(I)=CHROM1_INT(I)
			ELSE
				CHROM1_DOUBLE_NEW(I-N_INT)=CHROM2_DOUBLE(I-N_INT)
				CHROM2_DOUBLE_NEW(I-N_INT)=CHROM1_DOUBLE(I-N_INT)
			END IF

		END IF
	END DO

	CHROM1_INT=CHROM1_INT_NEW
	CHROM2_INT=CHROM2_INT_NEW
	CHROM1_DOUBLE=CHROM1_DOUBLE_NEW
	CHROM2_DOUBLE=CHROM2_DOUBLE_NEW


ELSE IF(TRIM(CROSS_TYPE).EQ."DOUBLE_POINT")THEN

	LOC1=RANDOM_INTEGER(SEED,N,0)
    IF (LOC1.EQ.0) LOC1=1
	LOC2=RANDOM_INTEGER(SEED,N,0)
    IF (LOC2.EQ.0) LOC2=1
	!MAKE SURE THE TWO LOCATIONS AREN'T THE SAME
	
	DO  WHILE(LOC1.EQ.LOC2)
		LOC2=RANDOM_INTEGER(SEED,N,0)
        IF (LOC2.EQ.0) LOC2=1
	END DO
	
	IF(LOC1.GT.LOC2)THEN
		LOC_SWAP=LOC2
		LOC2=LOC1
		LOC1=LOC_SWAP
	END IF

	DO I=1,N,1
		IF(I.LE.LOC1)THEN

			IF(I.LE.N_INT)THEN
				CHROM1_INT_NEW(I)=CHROM1_INT(I)
				CHROM2_INT_NEW(I)=CHROM2_INT(I)
			ELSE
				CHROM1_DOUBLE_NEW(I-N_INT)=CHROM1_DOUBLE(I-N_INT)
				CHROM2_DOUBLE_NEW(I-N_INT)=CHROM2_DOUBLE(I-N_INT)
			END IF

		ELSE IF(I.GT.LOC2)THEN

			IF(I.LE.N_INT)THEN
				CHROM1_INT_NEW(I)=CHROM1_INT(I)
				CHROM2_INT_NEW(I)=CHROM2_INT(I)
			ELSE
				CHROM1_DOUBLE_NEW(I-N_INT)=CHROM1_DOUBLE(I-N_INT)
				CHROM2_DOUBLE_NEW(I-N_INT)=CHROM2_DOUBLE(I-N_INT)
			END IF

		ELSE

			IF(I.LE.N_INT)THEN
				CHROM1_INT_NEW(I)=CHROM2_INT(I)
				CHROM2_INT_NEW(I)=CHROM1_INT(I)
			ELSE
				CHROM1_DOUBLE_NEW(I-N_INT)=CHROM2_DOUBLE(I-N_INT)
				CHROM2_DOUBLE_NEW(I-N_INT)=CHROM1_DOUBLE(I-N_INT)
			END IF

		END IF
	END DO

	CHROM1_INT=CHROM1_INT_NEW
	CHROM2_INT=CHROM2_INT_NEW
	CHROM1_DOUBLE=CHROM1_DOUBLE_NEW
	CHROM2_DOUBLE=CHROM2_DOUBLE_NEW

ELSE IF(TRIM(CROSS_TYPE).EQ."ARITHMETIC")THEN

	RAN=RANDOM(SEED)
	CHROM1_INT_NEW=INT(RAN*DBLE(CHROM1_INT)+&
        (1.D0-RAN)*DBLE(CHROM2_INT))
	CHROM2_INT_NEW=INT((1.D0-RAN)*DBLE(CHROM1_INT)+&
        RAN*DBLE(CHROM2_INT))
	CHROM1_DOUBLE_NEW=RAN*CHROM1_DOUBLE+(1.D0-RAN)*CHROM2_DOUBLE
	CHROM2_DOUBLE_NEW=(1.D0-RAN)*CHROM1_DOUBLE+RAN*CHROM2_DOUBLE

	!TEST TO MAKE SURE THE RESULTING CHROMOSOMES AREN'T OUTSIDE THE 
    !ALLOWED BOUNDS
	DO I=1,N_INT,1
		IF(CHROM1_INT_NEW(I).LT.INTEGER_LOWER(I)) THEN
			CHROM1_INT_NEW(I)=INTEGER_LOWER(I)
		ELSE IF(CHROM1_INT_NEW(I).GT.INTEGER_UPPER(I))THEN
			CHROM1_INT_NEW(I)=INTEGER_UPPER(I)
		END IF

		IF(CHROM2_INT_NEW(I).LT.INTEGER_LOWER(I)) THEN
			CHROM2_INT_NEW(I)=INTEGER_LOWER(I)
		ELSE IF(CHROM2_INT_NEW(I).GT.INTEGER_UPPER(I))THEN
			CHROM2_INT_NEW(I)=INTEGER_UPPER(I)
		END IF
	END DO

	CHROM1_INT=CHROM1_INT_NEW
	CHROM2_INT=CHROM2_INT_NEW
	CHROM1_DOUBLE=CHROM1_DOUBLE_NEW
	CHROM2_DOUBLE=CHROM2_DOUBLE_NEW

ELSE IF(TRIM(CROSS_TYPE).EQ."HEURISTIC")THEN

	IF(FIT1.LT.FIT2)THEN
		CHROM_INT_BEST=CHROM1_INT
		CHROM_INT_WORST=CHROM2_INT
		CHROM_DOUBLE_BEST=CHROM1_DOUBLE
		CHROM_DOUBLE_WORST=CHROM2_DOUBLE
	ELSE
		CHROM_INT_BEST=CHROM2_INT
		CHROM_INT_WORST=CHROM1_INT
		CHROM_DOUBLE_BEST=CHROM2_DOUBLE
		CHROM_DOUBLE_WORST=CHROM1_DOUBLE
	END IF

	TEST1=0
	COUNT=1

	DO WHILE (TEST1.EQ.0 .AND. COUNT.LT.50)
		RAN=RANDOM(SEED)

		CHROM1_INT_NEW=CHROM_INT_BEST+&
            INT(RAN*DBLE(CHROM_INT_BEST-CHROM_INT_WORST))
		CHROM1_DOUBLE_NEW=CHROM_DOUBLE_BEST+&
            RAN*(CHROM_DOUBLE_BEST-CHROM_DOUBLE_WORST)

		CHROM2_INT_NEW=CHROM_INT_BEST		
		CHROM2_DOUBLE_NEW=CHROM_DOUBLE_BEST

		TEST2=0

		DO I=I,N_INT,1

			IF(CHROM1_INT(I).LT.INTEGER_LOWER(I))THEN
				TEST2=1
			ELSE IF(CHROM1_INT(I).GT.INTEGER_UPPER(I))THEN
				TEST2=1
			END IF

			IF(CHROM2_INT(I).LT.INTEGER_LOWER(I))THEN
				TEST2=1
			ELSE IF(CHROM2_INT(I).GT.INTEGER_UPPER(I))THEN
				TEST2=1
			END IF

		END DO

		DO I=1,N_DOUBLE,1

			IF(CHROM1_DOUBLE(I).LT.DOUBLE_LOWER(I))THEN
				TEST2=1
			ELSE IF(CHROM1_DOUBLE(I).GT.DOUBLE_UPPER(I))THEN
				TEST2=1
			END IF

			IF(CHROM2_DOUBLE(I).LT.DOUBLE_LOWER(I))THEN
				TEST2=1
			ELSE IF(CHROM2_DOUBLE(I).GT.DOUBLE_UPPER(I))THEN
				TEST2=1
			END IF

		END DO

		IF(TEST2.EQ.0)THEN
			TEST1=1
		END IF

		COUNT=COUNT+1

	END DO
	
	IF(COUNT.GE.50) THEN
		CHROM1_INT_NEW=CHROM_INT_WORST
		CHROM1_DOUBLE_NEW=CHROM_DOUBLE_WORST
	END IF

	CHROM1_INT=CHROM1_INT_NEW
	CHROM2_INT=CHROM2_INT_NEW
	CHROM1_DOUBLE=CHROM1_DOUBLE_NEW
	CHROM2_DOUBLE=CHROM2_DOUBLE_NEW

ELSE
	WRITE(*,*) "INVALID CROSSOVER TYPE: NO CROSSOVER PERFORMED"
END IF

END SUBROUTINE CROSSOVER
    
    
    
    
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS SUBROUTINE EXECUTES THE MUTATION OPERATOR
SUBROUTINE MUTATION(CHROM_INT, CHROM_DOUBLE, N_INT, N_DOUBLE, SEED, &
  MUT_TYPE, INTEGER_LOWER, INTEGER_UPPER, DOUBLE_LOWER, DOUBLE_UPPER)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N_INT, N_DOUBLE, INTEGER_UPPER(N_INT), &
    INTEGER_LOWER(N_INT)

INTEGER, INTENT(INOUT) :: CHROM_INT(N_INT), SEED

DOUBLE PRECISION, INTENT(IN) :: DOUBLE_LOWER(N_DOUBLE), &
    DOUBLE_UPPER(N_DOUBLE)

DOUBLE PRECISION, INTENT(INOUT) :: CHROM_DOUBLE(N_DOUBLE)

CHARACTER(LEN=30), INTENT(IN) :: MUT_TYPE

DOUBLE PRECISION :: CHROM_DOUBLE_NEW(N_DOUBLE), D_UP, D_LOW, RAN
INTEGER :: CHROM_INT_NEW(N_INT), LOC1, LOC2, N, N_UP, N_LOW

N=N_INT+N_DOUBLE
IF(TRIM(MUT_TYPE).EQ."UNIFORM")THEN
    LOC1=RANDOM_INTEGER(SEED,N,1)
    
    CHROM_INT_NEW=CHROM_INT
    CHROM_DOUBLE_NEW=CHROM_DOUBLE
    
    IF(LOC1.LE.N_INT) THEN
        N_UP=INTEGER_UPPER(LOC1)
        N_LOW=INTEGER_LOWER(LOC1)
        CHROM_INT_NEW(LOC1)=RANDOM_INTEGER(SEED,N_UP, N_LOW)
    ELSE
        LOC2=LOC1-N_INT
        D_UP=DOUBLE_UPPER(LOC2)
        D_LOW=DOUBLE_LOWER(LOC2)
        CHROM_DOUBLE_NEW(LOC2)=RANDOM_DOUBLE(SEED, D_UP, D_LOW)
    END IF
    
    CHROM_INT=CHROM_INT_NEW
    CHROM_DOUBLE=CHROM_DOUBLE_NEW
    
ELSE IF(TRIM(MUT_TYPE).EQ."SLIDING")THEN

    LOC1=RANDOM_INTEGER(SEED,N,1)

    CHROM_INT_NEW=CHROM_INT
    CHROM_DOUBLE_NEW=CHROM_DOUBLE
    
    RAN=RANDOM(SEED)
    
    IF(RAN.LE.0.5D0)THEN
        IF(LOC1.LE.N_INT) THEN
            N_UP=CHROM_INT(LOC1)
            N_LOW=INTEGER_LOWER(LOC1)
            CHROM_INT_NEW(LOC1)=RANDOM_INTEGER(SEED,N_UP, N_LOW)
        ELSE
            LOC2=LOC1-N_INT
            D_UP=CHROM_DOUBLE(LOC2)
            D_LOW=DOUBLE_LOWER(LOC2)
            CHROM_DOUBLE_NEW(LOC2)=RANDOM_DOUBLE(SEED,D_UP, D_LOW)
        END IF
    ELSE
         IF(LOC1.LE.N_INT) THEN
            N_UP=INTEGER_UPPER(LOC1)
            N_LOW=CHROM_INT(LOC1)
            CHROM_INT_NEW(LOC1)=RANDOM_INTEGER(SEED,N_UP, N_LOW)
        ELSE
            LOC2=LOC1-N_INT
            D_UP=DOUBLE_UPPER(LOC2)
            D_LOW=CHROM_DOUBLE(LOC2)
            CHROM_DOUBLE_NEW(LOC2)=RANDOM_DOUBLE(SEED,D_UP, D_LOW)
        END IF  
    END IF
     
    CHROM_INT=CHROM_INT_NEW
    CHROM_DOUBLE=CHROM_DOUBLE_NEW
    
ELSE IF(TRIM(MUT_TYPE).EQ."BOUNDARY")THEN    

    LOC1=RANDOM_INTEGER(SEED,N,1)

    CHROM_INT_NEW=CHROM_INT
    CHROM_DOUBLE_NEW=CHROM_DOUBLE
    
    RAN=RANDOM(SEED)
    
    IF (RAN.LE.0.5D0)THEN
        IF(LOC1.LE.N_INT)THEN
            CHROM_INT_NEW(LOC1)=INTEGER_UPPER(LOC1)
        ELSE
            LOC2=LOC1-N_INT
            CHROM_DOUBLE_NEW(LOC2)=DOUBLE_UPPER(LOC2)
        END IF    
    ELSE
        IF(LOC1.LE.N_INT)THEN
            CHROM_INT_NEW(LOC1)=INTEGER_LOWER(LOC1)
        ELSE
            LOC2=LOC1-N_INT
            CHROM_DOUBLE_NEW(LOC2)=DOUBLE_LOWER(LOC2)
        END IF     
    END IF

    CHROM_INT=CHROM_INT_NEW
    CHROM_DOUBLE=CHROM_DOUBLE_NEW
    
ELSE
    WRITE(*,*)"INVALID MUTATION TYPE: NO MUTATION PERFORMED"
END IF


END SUBROUTINE MUTATION

  
  

!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS SUBROUTINE RANDOMLY GENERATES A CHROMOSOME 
!
SUBROUTINE CHROMOSOME_GENERATOR(SEED, N_INT,N_DOUBLE, CHROM_INT, &
    CHROM_DOUBLE, INT_UPPER, INT_LOWER, DOUBLE_LOWER, DOUBLE_UPPER)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: SEED, CHROM_INT(N_INT)

INTEGER, INTENT(IN) :: N_INT, N_DOUBLE, INT_UPPER(N_INT), &
    INT_LOWER(N_INT)

DOUBLE PRECISION, INTENT(INOUT) :: CHROM_DOUBLE(N_DOUBLE)

DOUBLE PRECISION, INTENT(IN) ::DOUBLE_LOWER(N_DOUBLE), &
    DOUBLE_UPPER(N_DOUBLE)
INTEGER :: I

DO I=1,N_INT,1
    CHROM_INT(I)=RANDOM_INTEGER(SEED, INT_UPPER(I), INT_LOWER(I))
END DO

DO I=1,N_DOUBLE,1
    CHROM_DOUBLE(I)=RANDOM_DOUBLE(SEED, DOUBLE_UPPER(I), &
        DOUBLE_LOWER(I))
END DO

END SUBROUTINE CHROMOSOME_GENERATOR
    
    
    
    
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS SUBROUTINE GENERATES THE INITIAL POPULATION
!
SUBROUTINE POPULATION_GENERATOR(N_POP, N_DOUBLE, N_INT, POP_DOUBLE, &
    POP_INT, SEED, INT_UPPER,INT_LOWER, DOUBLE_UPPER, DOUBLE_LOWER)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N_POP, N_DOUBLE, N_INT, INT_UPPER(N_INT), &
    INT_LOWER(N_INT)
INTEGER, INTENT(INOUT) :: POP_INT(N_POP,N_INT), SEED
DOUBLE PRECISION, INTENT(INOUT) :: POP_DOUBLE(N_POP, N_DOUBLE)
DOUBLE PRECISION, INTENT(IN) :: DOUBLE_UPPER(N_DOUBLE), &
    DOUBLE_LOWER(N_DOUBLE)

INTEGER :: I, J

DO I=1,N_POP,1

    DO J=1,N_DOUBLE,1
        POP_DOUBLE(I,J)=RANDOM_DOUBLE(SEED, DOUBLE_UPPER(J), &
            DOUBLE_LOWER(J))
    END DO

    DO J=1,N_INT,1
        POP_INT(I,J)=RANDOM_INTEGER(SEED, INT_UPPER(J), INT_LOWER(J))
    END DO
END DO

END SUBROUTINE POPULATION_GENERATOR

    
    
    
!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS SUBROUTINE GENERATES A RANDOM INTEGER WITHIN THE SPECIFIED 
!  BOUNDS
!
FUNCTION RANDOM_INTEGER(SEED,UPPER,LOWER)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: SEED
INTEGER, INTENT(IN) :: UPPER, LOWER
DOUBLE PRECISION :: RAN
INTEGER :: RANDOM_INTEGER, DUM

IF(UPPER.EQ.LOWER) THEN
    RANDOM_INTEGER=UPPER
ELSE
    RAN=RANDOM(SEED)
    RANDOM_INTEGER=FLOOR(RAN*DBLE(UPPER+1-LOWER))+LOWER
END IF
!WRITE(*,*) RANDOM_INTEGER, LOWER, UPPER, SEED, RAN
END FUNCTION RANDOM_INTEGER




!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS SUBROUTINE GENERATES A RANDOM DOUBLE PRECISION NUMBER WITHIN  
!  THE SPECIFIED BOUNDS
!
FUNCTION RANDOM_DOUBLE(SEED, UPPER,LOWER)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: SEED
DOUBLE PRECISION, INTENT(IN) :: UPPER, LOWER
DOUBLE PRECISION :: RAN, RANDOM_DOUBLE

    RAN=RANDOM(SEED)
    RANDOM_DOUBLE=RAN*(UPPER-LOWER)+LOWER

END FUNCTION RANDOM_DOUBLE


!********************************************************************!
!********************************************************************!
!********************************************************************!
!********************************************************************!
! THIS  FUNCTION IS THE RAN FUNCTION FROM NUMERICAL RECIPES FOR 
! FORTRAN 90.  IT USES THE RANDOM NUMBER GENERATOR OF PARK AND MILLER 
!  COMBINED WITH A  MARSAGLIA SHIFT SEQUENCE.  THE PERIOD OF THIS 
!  GENERATOR HAS A PERIOD OF ABOUT 3.1x10^18.  TO  INITIALIZE IT  
!  IDUM SHOULD BE SET  TO  A NEGATIVE INTEGER VALUE.   AFTER THAT 
!  IT'S VALUE SHOULDN'T BE CHANGED,EXCEPT TO REINITIALIZE.  THIS 
!  FUNCTION IS TAKEN FROM THE NUMERICAL RECIPES FOR FORTRAN BOOK
FUNCTION RANDOM(IDUM)                                         
IMPLICIT NONE                                                       
INTEGER, PARAMETER :: K4B=SELECTED_INT_KIND(9)                      
INTEGER(K4B), INTENT(INOUT) :: IDUM                                 
DOUBLE PRECISION:: RANDOM                              
INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647,IQ=127773,IR=2836
DOUBLE PRECISION, SAVE :: AAM                                        
INTEGER(K4B), SAVE :: IIX=-1, IIY=-1,KK                                
                                                                    
    IF(IDUM<=0 .OR. IIY<0) THEN                                      
        AAM=NEAREST(1.D0,-1.D0)/IM                                   
        IIY=IOR(IEOR(888889999,ABS(IDUM)),1)                         
        IIX=IEOR(777755555,ABS(IDUM))                                
        IDUM=ABS(IDUM)+1                                            
    END IF                                                          
                                                                    
    IIX=IEOR(IIX,ISHFT(IIX,13))                                        
    IIX=IEOR(IIX,ISHFT(IIX,-17))                                       
    IIX=IEOR(IIX,ISHFT(IIX,5))                                         
    KK=IIY/IQ                                                         
    IIY=IA*(IIY-KK*IQ)-IR*KK                                            
                                                                    
    IF (IIY < 0) IIY=IIY+IM                                            
    RANDOM=AAM*IOR(IAND(IM,IEOR(IIX,IIY)),1)                              
                                                                    
END FUNCTION RANDOM

!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!THIS  FUNCTION IS THE RAN FUNCTION FROM NUMERICAL RECIPES FOR FORTRAN 90.  IT USES THE RANDOM
!NUMBER GENERATOR OF PARK AND MILLER COMBINED WITH A  MARSAGLIA SHIFT SEQUENCE.  THE PERIOD OF 
!THIS GENERATOR HAS A PERIOD OF ABOUT 3.1x10^18.  TO  INITIALIZE IT  IDUM SHOULD BE SET  TO  A 
!NEGATIVE INTEGER VALUE.   AFTER THAT IT'S VALUE SHOULDN'T  BE CHANGED,EXCEPT TO REINITIALIZE.
!THIS FUNCTION HAS BEEN MODIFIED TO ALLOW  MULTIPLE  PARALLEL  STREAMS  OF THE FUNCTION  TO BE 
!UTILIZED. INSTEAD OF GLOBALLY SAVING EACH AX, IX, AND IY THEY ARE PASSED IN AND OUT AS PARAM-
!ETERS.  EACH STREAM SHOULD BE  INITIIALIZED WITH A NEGATIVE INTEGER IDUM.  FROM THEN ON  NONE 
!OF INPUT VARIABLES SHOULD BE MODIFIED.

FUNCTION RANDOM_STREAM(IDUM, IX, IY, K, AM)                                         
IMPLICIT NONE                                                       
INTEGER, PARAMETER :: K4B=SELECTED_INT_KIND(9)                      
INTEGER(K4B), INTENT(INOUT) :: IDUM
DOUBLE PRECISION:: RANDOM_STREAM                              
INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647,IQ=127773,IR=2836
DOUBLE PRECISION, INTENT(INOUT) :: AM
INTEGER, INTENT(INOUT) :: IX, IY, K
!DOUBLE PRECISION, SAVE :: AM                                        
!INTEGER(K4B), SAVE :: IX=-1, IY=-1,K                                
                                                                    
    IF(IDUM<=0 .OR. IY<0) THEN                                      
        AM=NEAREST(1.D0,-1.D0)/IM                                   
        IY=IOR(IEOR(888889999,ABS(IDUM)),1)                         
        IX=IEOR(777755555,ABS(IDUM))                                
        IDUM=ABS(IDUM)+1                                            
    END IF                                                          
                                                                    
    IX=IEOR(IX,ISHFT(IX,13))                                        
    IX=IEOR(IX,ISHFT(IX,-17))                                       
    IX=IEOR(IX,ISHFT(IX,5))                                         
    K=IY/IQ                                                         
    IY=IA*(IY-K*IQ)-IR*K                                            
                                                                    
    IF (IY < 0) IY=IY+IM                                            
    RANDOM_STREAM=AM*IOR(IAND(IM,IEOR(IX,IY)),1)                              
                                                                    
END FUNCTION RANDOM_STREAM

!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
SUBROUTINE BUBBLE_SORT(ARRAY, INDEX,N)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N
INTEGER, INTENT(INOUT) :: INDEX(N)
DOUBLE PRECISION, INTENT(INOUT) :: ARRAY(N)
INTEGER :: p,q

DO p=1,N,1
    INDEX(p)=p
END DO    

DO p=1,N,1
    DO q=N,p+1,-1
        CALL ORDER_SORT(ARRAY(q-1),ARRAY(q), INDEX(q-1), INDEX(q))
    END DO
END DO    
END SUBROUTINE BUBBLE_SORT

SUBROUTINE ORDER_SORT(A1,A2,I1,I2)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: I1, I2
DOUBLE PRECISION, INTENT(INOUT) :: A1, A2
DOUBLE PRECISION :: TEMPA
INTEGER :: TEMPI
    IF (A1>A2) THEN
        TEMPA=A1
        TEMPI=I1
        A1=A2
        I1=I2
        A2=TEMPA
        I2=TEMPI
     END IF   
END SUBROUTINE ORDER_SORT





!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!********************************************************************************************!
!*********************THIS ALGORITHM WAS OBTAINED FROM THE FOLLOWING WIKI********************!
!**http://en.wikibooks.org/wiki/Algorithm_Implementation/Sorting/Quicksort#FORTRAN_90.2F95***!
  Subroutine Qsort(X, Ipt)
! ***********************************
! * Sort Array X(:) in ascendent order 
! * If present Ipt, a pointer with the 
! * changes is returned in Ipt.
! ***********************************
 
    Type Limits
       Integer :: Ileft, Iright
    End Type Limits
 
    ! For a list with Isw number of elements or
    ! less use Insrt
    Integer, Parameter :: Isw = 10
 
    Real (kind=8), Intent (inout) :: X(:)
    Integer, Intent (out), Optional :: Ipt(:)
 
    Integer :: I, Ipvn, Ileft, Iright, ISpos, ISmax
    Integer, Allocatable :: IIpt(:)
    Type (Limits), Allocatable :: Stack(:)
 
 
    Allocate(Stack(Size(X)))
 
    Stack(:)%Ileft = 0
    If (Present(Ipt)) Then
       Forall (I=1:Size(Ipt)) Ipt(I) = I
 
       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)
 
       Do While (Stack(ISpos)%Ileft /= 0)
 
          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL InsrtLC(X, Ipt, Ileft,Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition(X, Ileft, Iright, Ipvn, Ipt)
 
             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do
 
    Else
 
       ! Iniitialize the stack
       Ispos = 1
       Ismax = 1
       Stack(ISpos)%Ileft  = 1
       Stack(ISpos)%Iright = Size(X)
 
       Allocate(IIpt(10))
       Do While (Stack(ISpos)%Ileft /= 0)
!          !WRITE(*,*)Ispos, ISmax
 
          Ileft = Stack(ISPos)%Ileft
          Iright = Stack(ISPos)%Iright
          If (Iright-Ileft <= Isw) Then
             CALL InsrtLC(X, IIpt, Ileft, Iright)
             ISpos = ISPos + 1
          Else
             Ipvn = ChoosePiv(X, Ileft, Iright)
             Ipvn = Partition(X, Ileft, Iright, Ipvn)
 
             Stack(ISmax+1)%Ileft = Ileft
             Stack(ISmax+1) %Iright = Ipvn-1
             Stack(ISmax+2)%Ileft = Ipvn + 1
             Stack(ISmax+2)%Iright = Iright
             ISpos = ISpos + 1
             ISmax = ISmax + 2
          End If
       End Do
       Deallocate(IIpt)
 
    End If
 
    Deallocate(Stack)
 
    Return
 
  CONTAINS
 
    ! ***********************************
    Integer Function ChoosePiv(XX, IIleft, IIright) Result (IIpv)
    ! ***********************************
    ! * Choose a Pivot element from XX(Ileft:Iright)
    ! * for Qsort. This routine chooses the median
    ! * of the first, last and mid element of the 
    ! * list.
    ! ***********************************
 
      Real (kind=8), Intent (in) :: XX(:)
      Integer, Intent (in) :: IIleft, IIright
 
      Real (kind=8) :: XXcp(3)
      Integer :: IIpt(3), IImd
 
      IImd = Int((IIleft+IIright)/2)
      XXcp(1) = XX(IIleft)
      XXcp(2) = XX(IImd)
      XXcp(3) = XX(IIright)
      IIpt = (/1,2,3/)
 
      CALL InsrtLC(XXcp, IIpt, 1, 3)
 
      Select Case (IIpt(2))
      Case (1)
         IIpv = IIleft
      Case (2)
         IIpv = IImd
      Case (3)
         IIpv = IIright
      End Select
 
      Return
    End Function ChoosePiv
 
    ! ***********************************
    Subroutine InsrtLC(XX, IIpt, IIl, IIr)
    ! ***********************************
    ! * Perform an insertion sort of the list 
    ! * XX(:) between index values IIl and IIr.
    ! * IIpt(:) returns the permutations
    ! * made to sort.
    ! ***********************************
 
      Real (kind=8), Intent (inout) :: XX(:)
      Integer, Intent (inout) :: IIpt(:)
      Integer, Intent (in) :: IIl, IIr
 
      Real (kind=8) :: RRtmp
      Integer :: II, JJ, IItmp
 
 
      Do II = IIl+1, IIr
         RRtmp = XX(II)
         Do JJ = II-1, 1, -1
            If (RRtmp < XX(JJ)) Then
               XX(JJ+1) = XX(JJ)
               CALL Swap_IN(IIpt, JJ, JJ+1)
            Else
               Exit
            End If
         End Do
         XX(JJ+1) = RRtmp
      End Do
 
      Return
    End Subroutine InsrtLC
 
  End Subroutine Qsort
 
! ***********************************
! *
  Integer Function Partition(X, Ileft, Iright, Ipv, Ipt) Result (Ipvfn)
! *
! ***********************************
! * This routine arranges the array X
! * between the index values Ileft and Iright
! * positioning elements smallers than
! * X(Ipv) at the left and the others 
! * at the right.
! * Internal routine used by Qsort.
! ***********************************
 
    Real (kind=8), Intent (inout) :: X(:)
    Integer, Intent (in) :: Ileft, Iright, Ipv
    Integer, Intent (inout), Optional :: Ipt(:)
 
    Real (kind=8) :: Rpv
    Integer :: I
 
    Rpv = X(Ipv)
    CALL Swap(X, Ipv, Iright)
    If (Present(Ipt)) CALL Swap_IN(Ipt, Ipv, Iright)
    Ipvfn = Ileft
 
    If (Present(Ipt))  Then
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap(X, I, Ipvfn)
             CALL Swap_IN(Ipt, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    Else
       Do I = Ileft, Iright-1
          If (X(I) <= Rpv) Then
             CALL Swap(X, I, Ipvfn)
             Ipvfn = Ipvfn + 1
          End If
       End Do
    End If
 
    CALL Swap(X, Ipvfn, Iright)
    If (Present(Ipt)) CALL Swap_IN(Ipt, Ipvfn, Iright)
 
    Return
  End Function Partition
 
! ***********************************
! *
  Subroutine Swap(X, I, J)
! *
! ***********************************
! * Swaps elements I and J of array X(:). 
! ***********************************
 
    Real (kind=8), Intent (inout) :: X(:)
    Integer, Intent (in) :: I, J
 
    Real (kind=8) :: Itmp
 
    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp
 
    Return
  End Subroutine Swap
 
! ***********************************
! *
  Subroutine Swap_IN(X, I, J)
! *
! ***********************************
! * Swaps elements I and J of array X(:). 
! ***********************************
 
    Integer, Intent (inout) :: X(:)
    Integer, Intent (in) :: I, J
 
    Integer :: Itmp
 
    Itmp = X(I)
    X(I) = X(J)
    X(J) = Itmp
 
    Return
  End Subroutine Swap_IN

!********************************************************************!  
!********************************************************************!  
!********************************************************************!  
!********************************************************************! 
! THIS SUBROUTINE SORTS ARR FROM SMALLEST TO LARGEST, THE LOCATION OF
!  THE SORTED ORDER IS RETURN IN THE INDEX ARRAY.  HEAP SORT HAS AT
!  MOST AND ORDER OF O(NLOG(N)).  THIS ALGORITHM IS ADAPTED FROM THE
!  NUMERICAL RECIPES FOR FORTRAN BOOK.
!
SUBROUTINE HEAP_SORT(ARR, INDEX, N)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N
INTEGER, INTENT(INOUT) :: INDEX(N)
DOUBLE PRECISION, INTENT(INOUT) :: arr(N)

INTEGER :: i, DUM_INT
DOUBLE PRECISION :: DUM_REAL

!n=size(arr)
do i=n/2,1,-1
	call sift_down(i,n)
end do
do i=n,2,-1
	!call swap(arr(1),arr(i))
    DUM_REAL=ARR(1)
    ARR(1)=ARR(I)
    ARR(I)=DUM_REAL
    DUM_INT=INDEX(1)
    INDEX(1)=INDEX(I)
    INDEX(I)=DUM_INT
	call sift_down(1,i-1)
end do
CONTAINS

SUBROUTINE sift_down(l,r)
INTEGER, INTENT(IN) :: l,r
INTEGER :: j,jold
DOUBLE PRECISION :: A
INTEGER :: B
a=arr(l)
B=INDEX(L)
jold=l
j=l+l
do
	if (j > r) exit
	if (j < r) then
		if (arr(j) < arr(j+1)) j=j+1
	end if
	if (a >= arr(j)) exit
    arr(jold)=arr(j)
    INDEX(JOLD)=INDEX(J)
	jold=j
	j=j+j
end do
arr(jold)=a
INDEX(JOLD)=B
END SUBROUTINE sift_down
END SUBROUTINE HEAP_SORT

!********************************************************************!  
!********************************************************************!  
!********************************************************************!  
!********************************************************************! 
! THIS SUBROUTINE IS THE WRAPPER FOR THE UNCMIN PACKAGE. THE UNCMIN 
!  PACKAGES WAS MODIFIED TO PASS THE INTEGER CHROMOSOME/LENGTH, 
!  INPUT ARRAY/SIZES, AND THE CONSTRAINT ARRAY G_CON/LENGTH.  THE 
!  UNCMIN PACKAGE DOESN'T ACTUALLY USE THE CONSTRAINTS, BUT IT WAS 
!  MODIFIED SO THE SAME COST FUNCTION CAN BE USED FOR ALL 3 NLP 
!  SOLVERS.
!
SUBROUTINE UNCMIN_WRAPPER(N_DOUBLE, N_INT, N1, N2, ITMAX_NLP,  X0, &
    X, CHROM_INT, FITNESS, INPUT_ARRAY, G_CON, NCON)
IMPLICIT NONE

INTEGER, INTENT(IN) :: N_DOUBLE, N_INT, N1, N2, CHROM_INT(N_INT), NCON

INTEGER, INTENT(INOUT) :: ITMAX_NLP

DOUBLE PRECISION, INTENT(INOUT) ::  X(N_DOUBLE), X0(N_DOUBLE), &
    FITNESS, G_CON(NCON),INPUT_ARRAY(N1,N2)

INTEGER :: INFO

CALL UNCMIN(N_DOUBLE, N_INT, N1, N2, ITMAX_NLP, INFO, X0, X, &
    CHROM_INT, FITNESS, INPUT_ARRAY, G_CON, NCON)

    END SUBROUTINE UNCMIN_WRAPPER
    
    
    
    
!********************************************************************!  
!********************************************************************!  
!********************************************************************!  
!********************************************************************! 
! THIS SUBROUTINE IS THE WRAPPER FOR THE CONMIN PACKAGE.  VERY LITTLE 
!  MODIFICATION IS REQUIRED OF THE ORIGINAL CONMIN PACKAGE.  WITH
!  THIS VERSION OF CONMIN NO MORE THAN 2200 VARIABLES CAN BE 
!  OPTIMIZED. THIS NLP SOLVER SHOULDN'T BE USED WITH PROBLEMS 
!  CONTAINING MORE THAN 100-200 VARIABLES.  IF YOU HAVE THAT MANY 
!  VARIABLES RUN TIMES WILL BE HIGH.
!
SUBROUTINE CONMIN_WRAPPER(NDV, N_INT, N1, N2, ITMAX, X0, X_OUT, &
                            CHROM_INT, FITNESS, INPUT_ARRAY, N_CON, &
                            DOUBLE_UPPER, DOUBLE_LOWER)
IMPLICIT NONE

INTEGER, INTENT (IN) :: NDV, N_INT, N1, N2, ITMAX, CHROM_INT(N_INT), &
    N_CON

DOUBLE PRECISION, INTENT(IN) :: DOUBLE_UPPER(NDV), DOUBLE_LOWER(NDV),&
    X0(NDV)

DOUBLE PRECISION, INTENT(INOUT) :: FITNESS, INPUT_ARRAY(N1,N2), &
    X_OUT(NDV)

DOUBLE PRECISION :: X(22*100), VUB(22*100), VLB(22*100)

INTEGER :: ISC(50*100), IER, IPRINT, NSIDE

X(1:NDV)=X0
VUB(1:NDV)=DOUBLE_UPPER
VLB(1:NDV)=DOUBLE_LOWER
! 1 = ENFORCE BOUNDARY CONDITION
! 0 = NO BOUNDARY CONDITIONS
NSIDE=1
! 0-5 GOES FROM PRINTING NOTHING TO PROGRESSIVELY MORE AND MORE
IPRINT=0
CALL CMINEX (X, VLB, VUB, ISC, NDV, N_CON, NSIDE, IPRINT, &
             ITMAX, IER, N1, N2, N_INT, CHROM_INT, FITNESS, &
             INPUT_ARRAY)

X_OUT=X(1:NDV)
                            
END SUBROUTINE CONMIN_WRAPPER
                          
                          

!********************************************************************************************!
!********************************************************************************************!
!***********************                                               **********************!
!***********************       MONOTONIC BASIN HOPPING ALGORITHMS      **********************!
!***********************                                               **********************!
!********************************************************************************************!
!********************************************************************************************!
SUBROUTINE MBH(ND, NI, N1, N2, ITMAX_MBH, ITMAX_NLP, SEED, X0, X, X_LB, X_UB, CHROM_INT, &
               FITNESS, INPUT_ARRAY, NCON, RADIUS)
IMPLICIT NONE
INTEGER, INTENT(IN) :: ND, NI, N1, N2, ITMAX_MBH, CHROM_INT(NI), NCON
INTEGER, INTENT(INOUT) :: SEED, ITMAX_NLP
DOUBLE PRECISION, INTENT(IN) :: RADIUS
DOUBLE PRECISION, INTENT(INOUT) :: INPUT_ARRAY(N1,N2)
DOUBLE PRECISION, INTENT(INOUT) :: X0(ND), X(ND), FITNESS, X_LB(ND), X_UB(ND)

INTEGER :: ITER, I,J, ITDUM

DOUBLE PRECISION :: X_TEMP(ND), X_PRIME(ND), FITNESS_TEMP, G_CON(NCON)


ITDUM=150
!CALL UNCMIN_WRAPPER(ND, NI, N1, N2, ITMAX_NLP, X0, X_TEMP, CHROM_INT, FITNESS, INPUT_ARRAY, G_CON, NCON)
!CALL COST(ND, NI, N1, N2, X0, CHROM_INT, FITNESS, INPUT_ARRAY, G_CON, ncon )
CALL CONMIN_WRAPPER(ND, NI, N1, N2, ITDUM, X0, X, CHROM_INT, FITNESS, INPUT_ARRAY, NCON, &
                            X_UB, X_LB)
X=X0
write(*,*) fitness
!CHECK TO MAKE SURE  A  FEASIBLE  SOLUTION HAS  BEEN FOUND  AND THAT THE UPPER AND LOWER BOUND
!LIMITS HAVE NOT BEEN EXCEEDED.  IF EVERYTHING IS OK THE MONOTONIC BASIN HOPPING ALGORITHM CAN
!PROCEED.

!DO I=1,ND,1
!    IF(X_TEMP(I) .LT. X_LB(I) .OR. X_TEMP(I).GT.X_UB(I))THEN
!        FITNESS=1.1D30
!    END IF
!END DO
!WRITE(*,*) ISNAN(FITNESS)
IF( FITNESS .GE. 1.D24) THEN! .OR. ISNAN(FITNESS)) THEN
    !THE  STARTING  GUESS  WASN'T  A  FEASIBLE  SOLUTION.   THE ALGORITHM  EXITS AND RETURNS A 
    !FITNESS VALUE OF 1.D30
    !FITNESS=1.D30
    !WRITE(*,*) FITNESS
ELSE
    ITER=0
    !X_TEMP=
    !WRITE(*,*) FITNESS
    DO WHILE (ITER.LE.ITMAX_MBH)
        ITER=ITER+1
        
        CALL PERTURBING_FUNCTION(ND, SEED, X, X_PRIME, X_LB, X_UB, RADIUS)
        CALL CONMIN_WRAPPER(ND, NI, N1, N2, ITMAX_NLP, X_PRIME, X_TEMP, &
                            CHROM_INT, FITNESS_TEMP, INPUT_ARRAY, NCON, &
                            X_UB, X_LB)

        !write(*,*) "conmin wrapper exited"
        
        ! MAKE SURE THE BOUNDS AREN'T VIOLATED.  THIS SHOULDN'T HAPPEN WITH THE CONMIN SOLVER
        !DO I=1,ND,1
        !    IF(X_TEMP(I) .LT. X_LB(I) .OR. X_TEMP(I).GT.X_UB(I)) THEN
        !        WRITE(*,*) I, X_PRIME(I), X_TEMP(I), X_LB(I), X_UB(I)
        !        FITNESS_TEMP=1.1D30
        !    END IF
        !END DO
        
        IF(FITNESS_TEMP.LT.FITNESS)THEN
            X=X_TEMP
            ITER=0
            FITNESS=FITNESS_TEMP
            WRITE(*,*) "BEST FITNESS", FITNESS
        END IF
    END DO
END IF

!RANDOM_DOUBLE(SEED, UPPER,LOWER)

END SUBROUTINE MBH
    
SUBROUTINE PERTURBING_FUNCTION(ND, SEED, X, X_PRIME, X_LB, X_UB, R)    
IMPLICIT NONE
INTEGER, INTENT (IN) :: ND
INTEGER, INTENT (INOUT) :: SEED
DOUBLE PRECISION, INTENT(IN) :: X(ND), X_LB(ND), X_UB(ND), R
DOUBLE PRECISION, INTENT(INOUT) :: X_PRIME(ND)

INTEGER :: I, COUNTER
DOUBLE PRECISION :: RAN, RAN_CAUCHY, xb1(nd), xb2(nd)


X_PRIME=1.D8

!R=0.05D0

xb1=-r*(x_ub-x_lb)
xb2=r*(x_ub-x_lb)
    
!WRITE(*,*) "XB1", XB1
!WRITE(*,*) "XB2", XB2
DO I=1,ND,1
    !THE WHILE LOOP IS NEED SO ENSURE A PERTURBED SOLUTION IS FOUND THAT DOESN'T VIOLATE THE
    !VARIABLE BOUNDS
    COUNTER=0
    DO WHILE (X_PRIME(I).GT.X_UB(I) .OR.  X_PRIME(I) .LT. X_LB(I)) 
        COUNTER=COUNTER+1
        RAN=RANDOM_DOUBLE(SEED, Xb2(I), Xb1(I))
        !WRITE(*,*) RAN, XB1(I), XB2(I)
        !CONVERT  TO  CAUCHY  RANDON N UMBER  (WILL ENSURE MOST OF THE JUMPS ARE SMALL,  
        !BUT STILL ALLOWS FOR SOME LARGE JUMPS  LESS OFTEN THAN THE UNIFORM DISTRIBUTION)
        !RAN_CAUCHY=1.D0/(PI*(1.D0+RAN**2))
        
        X_PRIME(I)=X(I)+RAN
        !write(*,*) x_prime(i), x_LB(i), x_UB(i), ran, counter
        
    END DO
        IF (X_PRIME(I).GT.X_UB(I)) THEN
            WRITE(*,*) X_PRIME(I), X_UB(I)
        END IF
        IF (X_PRIME(I).LT.X_LB(I)) THEN
            WRITE(*,*) X_PRIME(I), X_LB(I)
        END IF
END DO

END SUBROUTINE

END MODULE OPTIMIZATION_MODULE

