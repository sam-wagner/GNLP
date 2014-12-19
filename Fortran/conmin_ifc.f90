MODULE CNMN1
IMPLICIT NONE
DOUBLE PRECISION, SAVE :: DELFUN, DABFUN, FDCH, FDCHM, CT, CTMIN, CTL, &
          CTLMIN, ALPHAX, ABOBJ1, THETA, OBJ
INTEGER, SAVE :: IDUM1, IDUM2, &
          IDUM3, IDUM4, NFDG, NSCAL, LINOBJ, IDUM5, ITRM, &
          ICNDIR, IGOTO, NAC, INFO, INFOG, ITER, IPRINT, NCON, NDV, NSIDE, ITMAX

END MODULE CNMN1

MODULE CONSAV
IMPLICIT NONE 

DOUBLE PRECISION, SAVE :: DM1,DM2,DM3,DM4,DM5,DM6,DM7,DM8,DM9,&
    DM10,DM11,DM12, DCT,DCTL,PHI,ABOBJ,CTA,CTAM,CTBM,OBJ1,SLOPE, &
    DX,DX1,FI,XI,DFTDF1, ALP,FFF,A1,A2,A3,A4,F1,F2,F3,F4,CV1,CV2, &
    CV3,CV4,APP,ALPCA,ALPFES,ALPLN, ALPMIN,ALPNC,ALPSAV,ALPSID, &
    ALPTOT,RSPACE

INTEGER, SAVE :: IOBJ,KOBJ,KCOUNT,NCAL(2),NFEAS,MSCAL,NCOBJ,NVC, &
    KOUNT,ICOUNT,IGOOD1, IGOOD2,IGOOD3,IGOOD4,IBEST,III,NLNC, &
    JGOTO,ISPACE(2),IDM1,IDM2,IDM3,JDIR   
END MODULE CONSAV
          


    
    
MODULE CONSTRAINED_MINIMIZATION
USE COST_MODULE  
IMPLICIT NONE

CONTAINS
	  

      
! LAST MODIFIED BY SAM WAGNER AT THE ASTEROID DEFLECTION RESEARCH
! CENTER-IOWA STATE UNIVERSITY
!
!  THIS ALGORITHM WAS MODIFIED BY PUTTING IT INTO A MODULE AND REQUIRE A 
!    MODULE CALLED COST MODULE THAT CONTAINS THE OBJECT FUNCTION,
!    WHICH MUST BE NAMED COST.  IT WAS ALSO MODIFIED TO BE COMPATABLE
!    WITH THE COST FUNCTIONS FOR THE GA-NLP ALGORITHM
!
!  ADDED VARIABLES:
!    
!

!********************************************************************!
!********************************************************************!
!                                                            
!  Murphy Laboratory   ---   Apollo Network   ---   Math Library    
!
!  Title:       CMINEX                   Date Created:   9-26-91    
!  Programmer:  Nick Thorp               Last Modified:  9-27-91    
!
!  DESCRIPTION:                                                     
!    This subroutine may be used to determine the set of independent
!    (design) variables, x, that will minimize a function, f(x),    **
!    subject to a set of constraints, g(x) <= 0, and a set of side  **
!    constraints, lower bound <= x <= upper bound.  In terms of the **
!    code's variables; determine the X(I), I=1,NDV that will        **
!    minimize OBJ subject to the set of constraints G(J) .LE. 0,    **
!    J=1,NCON and a set of side constraints VLB(I) .LE. X(I) .LE.   **
!    VUB(I), I=1,NDV.                                               **
!                                                                   **
!    CMINEX is just an interface between the user and the actual    **
!    optimization code CONMIN.                                      **
!                                                                   **
!  ARGUMENTS:                                                       **
!  Input:                                                           **
!    ANALYZ  = The name of the subroutine that will be used to      **
!             evaluate the problem's objective function and         **
!             constraints.                                          **
!    X      = A 1-D array (length: NDV) containing the initial      **
!             approximations to the optimal values of the set of    **
!             independent (design) variables.                       **
!    VLB    = A 1-D array (length: NDV) containing the lower limits **
!             on the values of the independent variables.  VLB need **
!             only be defined if NSIDE = 1.                         **
!    VUB    = A 1-D array (length: NDV) containing the upper limits **
!             on the values of the independent variables.  VUB need **
!             only be defined if NSIDE = 1.                         **
!    ISC    = A 1-D array (length: NCON) used to identify the type  **
!             of a constraint.                                      **
!       (i) = 0 -- The ith constraint, G(I), is non-linear.         **
!       (i) = 1 -- The ith constraint, G(I), is     linear.         **
!             If you are in doubt, set ISC(I) to 0.                 **
!    NDV    = The number of independent (design) variables          **
!             contained in X.                                       **
!    NCON   = The number of constraints contained in G.  NCON may   **
!             be 0.                                                 **
!    NSIDE  = The side constraint identifier.                       **
!           .NE. 0 -- Upper and lower bounds on X are imposed by    **
!                     VUB and VLB.                                  **
!           .EQ. 0 -- No bounds on X are imposed.  VUB and VLB are  **
!                     ignored.                                      **
!    IPRINT = Print control identifier.  IPRINT may take any        **
!             integer value from 0 to 5. IPRINT = 0 will suppress   **
!             all printing during the optimization process.  The    **
!             larger the value of IPRINT, the more output that is   **
!             generated during the optimization.                    **
!**    ITMAX  = The maximum number of allowed iterations.             **
!**                                                                   **
!**  Output:                                                          **
!**    X      = A 1-D array (length: NDV) containing the final        **
!**             approximations to the optimal values of the set of    **
!**             independent (design) variables.                       **
!**    IER    = An error code used to indicate the success or failure **
!**             of the optimization operation.                        **
!**           =  0 -- A successful optimization attempt.              **
!**           = -1 -- NDV is out of range.                            **
!**           = -2 -- NCON is out of range.                           **
!**                                                                   **
!**  NOTES:                                                           **
!**    1. Argument variable type declarations:                        **
!**       a) REAL:     X, VLB, VUB                                    **
!**       b) INTEGER:  ISC, NDV, NCON, NSIDE, IPRINT, ITMAX, IER      **
!**    2. Single precision arithmetic is used throught the entire     **
!**       subroutine.                                                 **
!**    3. The program module (main program, subroutine, etc.) which   **
!**       calls CMINEX should contain an EXTERNAL statement for the   **
!**       name of the subroutine which will evaluate the objective    **
!**       function, OBJ, and the constraints, G(I), I=1,NCON.         **
!**    4. All of the output from CMINEX and CONMIN will be sent to    **
!**       the screen.  If you want to have this output in a file, try **
!**       piping the output to a file, i.e., if your executable file  **
!**       is called example.x use the following UNIX command:         **
!**       UNIX:  example.x >example.out                               **
!**       This will cause everything that was sent to the screen to   **
!**       be sent to the file example.out instead. Thus, example.x    **
!**       could not be an interactive program since all of the screen **
!**       prompts would end up in example.out and not on the screen.  **
!**    5. The named COMMON blocks CNMN1 and CONSAV are used by CMINEX **
!**       and the CONMIN subroutines.  Thus the user should not have  **
!**       any named COMMON blocks with the names CNMN1 or CONSAV.     **
!**    6. The CONMIN subroutines are named CNMN01, CNMN02, CNMN03,    **
!**       CNMN04, CNMN05, CNMN06, CNMN07, and CNMN08.  Thus, the user **
!**       should not have any subroutine, function, or block data     **
!**       subprograms with these names.                               **
!**    7. The dimensions in this subroutine are currently set for a   **
!**       maximum of 20 independent (design) variables and a maximum  **
!**       of 10 constraints.  In the program module which calls       **
!**       CMINEX, the arrays X, VLB, VUB, and ISC should have the     **
!**       following dimensions:                                       **
!**       X(MAXN1), VLB(MAXN1), VUB(MAXN1) where MAXN1 = 22.          **
!**       ISC(MAXN2)                       where MAXN2 = 50.          **
!**    8. Of the argument variables, only the array X will have its   **
!**       values modified (AS EXPLAINED ABOVE) on exit from CMINEX.   **
!**                                                                   **
!**  REFERENCE:                                                       **
!**    "CONMIN - A FORTRAN Program for Constrained Function           **
!**     Minimization, User's Manual", G. N. Vanderplaats,             **
!**     NASA TMX-62,282, Aug., 1973.                                  **
!**                                                                   **
!**  REQUIRED SUBPROGRAMS:                                            **
!**  Non-User Supplied:                                               **
!**    CONMIN                                                         **
!**                                                                   **
!**  User Supplied:                                                   **
!**  If NCON >= 0:                                                     **
!**    ANALYZ (X, NDV, OBJ, G, NCON)                                   **
!**      A subroutine which will evaluate the objective function, OBJ **
!**      and the constraints, G(I), I=1,NCON, at the given set of     **
!**      independent (design) variables, X(J), J=1,NDV.               **
!**      REAL:      X, OBJ, G                                         **
!**      INTEGER:   NDV, NCON                                         **
!**                                                                   **
!********************************************************************!
!********************************************************************!
SUBROUTINE CMINEX (X, VLB, VUB, ISC, NDV_DUM, N_CON, NSIDE_DUM, IIPRINT, &
    ITMAX_DUM, IER, N1, N2, N_INT, ICHROM_INT, FITNESS, ARRAY_INPUT, &
    MAXN1, MAXN2, MAXN3, MAXN4, MAXN5)
USE CNMN1
IMPLICIT NONE 

!PARAMETER (MAXN1 = 22*100)
!PARAMETER (MAXN2 = 50*100)
!PARAMETER (MAXN3 = 21*100)
!PARAMETER (MAXN4 = 21*100)
!PARAMETER (MAXN5 = 42*100)

INTEGER, INTENT(IN) :: NDV_DUM, N_CON, NSIDE_DUM, IIPRINT, ITMAX_DUM, &
    N1, N2, N_INT,  ICHROM_INT(N_INT)
DOUBLE PRECISION, INTENT(INOUT) :: X(MAXN1), ARRAY_INPUT(N1,N2), &
    VLB(MAXN1), VUB(MAXN1), FITNESS
INTEGER, INTENT(INOUT) :: ISC(MAXN2), IER, MAXN1, MAXN2, MAXN3, MAXN4, MAXN5

!LOCAL VARIABLE
DOUBLE PRECISION :: G(MAXN2), SCAL(MAXN1), DF(MAXN1), &
    A(MAXN1,MAXN3), S(MAXN1), G1(MAXN2), G2(MAXN2), &
    B(MAXN3,MAXN3), C(MAXN4) 
INTEGER :: IC(MAXN3), MS1(MAXN5)

!INTEGER ::MAXN1, MAXN2, MAXN3, MAXN4, MAXN5
!MAXN1 = 22*100
!MAXN2 = 50*100
!MAXN3 = 21*100
!MAXN4 = 21*100
!MAXN5 = 42*100




!
!      COMMON /CNMN1/ DELFUN, DABFUN, FDCH, FDCHM, CT, CTMIN, CTL, &
!          CTLMIN, ALPHAX, ABOBJ1, THETA, OBJ, IDUM1, IDUM2, &
!          IDUM3, IDUM4, NFDG, NSCAL, LINOBJ, IDUM5, ITRM, &
!          ICNDIR, IGOTO, NAC, INFO, INFOG, ITER




!-----Making the argument list and the COMMON block compatible:

      
      IDUM1 = NDV_DUM
      IDUM2 = N_CON
      IDUM3 = NSIDE_DUM
      IDUM4 = IIPRINT
      IDUM5 = ITMAX_DUM
      
      IPRINT=IIPRINT
      NCON=N_CON
      NDV=NDV_DUM
      NSIDE=NSIDE_DUM
      ITMAX=ITMAX_DUM
!-----Check to be sure that the problem is not too big
!-----for the current set of parameters:

      IF (NDV .LE. 1 .OR. NDV .GT. (MAXN1 - 2)) THEN
        IER = -1
        RETURN
      END IF

      IF (NCON .LT. 0 .OR. NCON .GT. (MAXN2 - 2 * NDV)) THEN
        IER = -2
        RETURN
      END IF

      IF (IPRINT .LT. 0 .OR. IPRINT .GT. 5) IPRINT = 1

      ITMAX = ABS(ITMAX)

!-----Set the default parameters:
      DELFUN = 0.0
      DABFUN = 0.0
      FDCH   = 0.0
      FDCHM  = 0.0
      !CT     = -1.D-3
      CT     = 0.0
      CTMIN  = 0.0
      CTL    = 0.0
      CTLMIN = 0.0
      ALPHAX = 0.0
      ABOBJ1 = 0.0
      THETA  = 0.0
!	ASH: IF NFDG = 0	THEN GRADIENT IS CALCULATED BY FINITE DIFF METHOD
!	ASH: IF NSCAL =0	THEN INDEPENDANT VARIBALES ARE NOT SCALED 
      NFDG   = 0
      NSCAL  = 0
      LINOBJ = 0
      ITRM   = 5
      ICNDIR = 0
!
!-----Optimize:
!
      IGOTO = 0
!
10 CALL CONMIN (X, VLB, VUB, G, SCAL, DF, A, S, G1, G2, B, C, ISC, &
   IC, MS1, MAXN1, MAXN2, MAXN3, MAXN4, MAXN5)
!
!-----Check to see if we are done:
!
      IF (IGOTO .EQ. 0) THEN
        IER = 0
        RETURN
      END IF
!
!-----Evaluate the objective function and the constraints:
!
      CALL COST(NDV, N_INT, N1, N2, X(1:NDV), ICHROM_INT, OBJ,  &
          ARRAY_INPUT, G(1:NDV), NCON)
         FITNESS=OBJ

      GO TO 10
!
!-----Exit from the subroutine:
!
!      RETURN
END SUBROUTINE CMINEX
        
        
        
        
        
        

!**********************************************************************!
!**********************************************************************!
!**********************************************************************!
!     ROUTINE TO SOLVE CONSTRAINED OR UNCONSTRAINED FUNCTION
!     MINIMIZATION.
!     BY G. N. VANDERPLAATS                          APRIL, 1972.
!     * * * * * * * * * * *   MAY, 1978 VERSION    * * * * * * * * * * *
!     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
!     REFERENCE;  CONMIN - A FORTRAN PROGRAM FOR CONSTRAINED FUNCTION
!         MINIMIZATION:  USER'S MANUAL,  BY G. N. VANDERPLAATS,
!         NASA TM X-62,282, AUGUST, 1973.
!     STORAGE REQUIREMENTS:
!         PROGRAM - 7000 DECIMAL WORDS (CDC COMPUTER)
!         ARRAYS  - APPROX. 2*(NDV**2)+26*NDV+4*NCON,
!               WHERE N3 = NDV+2.
!     RE-SCALE VARIABLES IF REQUIRED.
!**********************************************************************!
!**********************************************************************!
!**********************************************************************!
!*****                                                            *****!
!*****               This version of CONMIN was                   *****!
!*****               modified to execute on the                   *****!
!*****               APOLLO SERIES 3000.  The                     *****!
!*****               source code was compiled                     *****!
!*****               using the FORTRAN 77                         *****!
!*****               compiler.                                    *****!
!*****                                                            *****!
!*****               Jerald M. Vogel                              *****!
!*****               Iowa State University                        *****!
!*****               September 26, 1983                           *****!
!*****                                                            *****!
!*****                                                            *****!
!*****                                                            *****!
!**********************************************************************!
!**********************************************************************!
!**********************************************************************!
SUBROUTINE CONMIN (X,VLB,VUB,G,SCAL,DF,A,S,G1,G2,B,C,ISC,IC,MS1,N1, &
N2,N3,N4,N5)
USE CNMN1
USE CONSAV
IMPLICIT NONE !DOUBLE (A-H,O-Z)
DOUBLE PRECISION, INTENT(INOUT) :: X(N1), VLB(N1), VUB(N1), G(N2), &
    SCAL(N1), DF(N1), A(N1,N3), S(N1), G1(N2), G2(N2), B(N3,N3), C(N4)
INTEGER, INTENT (INOUT) :: ISC(N2), IC(N3), MS1(N5), N1, N2, N3, N4, N5


!COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,ALPHAX, &
!ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,IT, &
!RM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER
!DIMENSION X(N1),VLB(N1),VUB(N1),G(N2),SCAL(N1),DF(N1),A(N1,N3),S(N1),G1(N2),G2(N2),B(N3,N3),C(N4),ISC(N2),IC(N3),MS1(N5)
!COMMON /CONSAV/ DM1,DM2,DM3,DM4,DM5,DM6,DM7,DM8,DM9,DM10,DM11,DM12, &
!    DCT,DCTL,PHI,ABOBJ,CTA,CTAM,CTBM,OBJ1,SLOPE,DX,DX1,FI,XI,DFTDF1, &
!    ALP,FFF,A1,A2,A3,A4,F1,F2,F3,F4,CV1,CV2,CV3,CV4,APP,ALPCA,ALPFES,ALPLN, &
!    ALPMIN,ALPNC,ALPSAV,ALPSID,ALPTOT,RSPACE,IDM1,IDM2,IDM3,JDIR, &
!    IOBJ,KOBJ,KCOUNT,NCAL(2),NFEAS,MSCAL,NCOBJ,NVC,KOUNT,ICOUNT,IGOOD1, &
!    IGOOD2,IGOOD3,IGOOD4,IBEST,III,NLNC,JGOTO,ISPACE(2)    

!LOCAL VARIABLES
DOUBLE PRECISION :: ALP1, ALP11, ALP12, C1, CT1, CTC, FF1, GI, OBJB, OBJD, SCJ, SI, &
    SIB, X1, X12, XID, XX
INTEGER :: I, II, IP1, J, K, M1, M2, M3, MCN1, NCI, NDV1, NDV2, NIC, NNAC


      IF (NSCAL.EQ.0.OR.IGOTO.EQ.0) GO TO 20
      DO 10 I=1,NDV
10    X(I)=C(I)
20    CONTINUE
!     CONSTANTS.
      NDV1=NDV+1
      NDV2=NDV+2
      IF (IGOTO.EQ.0) GO TO 30
      GO TO (150,370,360,650,670), IGOTO
!     ------------------------------------------------------------------
!                      SAVE INPUT CONTROL PARAMETERS
!     ------------------------------------------------------------------
30    CONTINUE
      IF (IPRINT.GT.0) WRITE (6,1230)
      IF (LINOBJ.EQ.0.OR.(NCON.GT.0.OR.NSIDE.GT.0)) GO TO 40
!     TOTALLY UNCONSTRAINED FUNCTION WITH LINEAR OBJECTIVE.
!     SOLUTION IS UNBOUNDED.
      WRITE (6,970) LINOBJ,NCON,NSIDE
      RETURN
40    CONTINUE
      IDM1=ITRM
      IDM2=ITMAX
      IDM3=ICNDIR
      DM1=DELFUN
      DM2=DABFUN
      DM3=CT
      DM4=CTMIN
      DM5=CTL
      DM6=CTLMIN
      DM7=THETA
      DM8=PHI
      DM9=FDCH
      DM10=FDCHM
      DM11=ABOBJ1
      DM12=ALPHAX
!     ------------------------------------------------------------------
!                                DEFAULTS
!     ------------------------------------------------------------------
      IF (ITRM.LE.0) ITRM=30
      IF (ITMAX.LE.0) ITMAX=200
      NDV1=NDV+1
      IF (ICNDIR.EQ.0) ICNDIR=NDV1
!HINDMAN 3/4/2002      IF (DELFUN.LE.0.) DELFUN=.0001
      IF (DELFUN.LE.0.D0) DELFUN=.00000000001
      CT=-ABS(CT)
      IF (CT.GE.0.D0) CT=-0.1D0
      CTMIN=ABS(CTMIN)
!	ASH: CTMIN ORGINAL = 0.004
      IF (CTMIN.LE.0.D0) CTMIN=.004D0
      CTL=-ABS(CTL)
      IF (CTL.GE.0.D0) CTL=-0.01D0
      CTLMIN=ABS(CTLMIN)
      IF (CTLMIN.LE.0.D0) CTLMIN=.001D0
!	ASH: 	THIS DROPPED THE OBJ WHEN THETA = 1 TO 10
!		ALSO, THE RESULTS BETWEEN THE LINEAR AND COSINE DISTR MATCHED			
!		WHEN ALL THESE VALUES WERE DROPPED THE ITER WAS MUCH LESS
!		AND THE OBJ WENT DOWN TO 10-6 MAGNITUDE 
!      IF (THETA.LE.0.) THETA=1000.0
!      IF (ABOBJ1.LE.0.) ABOBJ1=.000001
!      IF (ALPHAX.LE.0.) ALPHAX=.000001
!      IF (FDCH.LE.0.) FDCH=.000001
!      IF (FDCHM.LE.0.) FDCHM=.000001
     
      IF (THETA.LE.0.D0) THETA=1000.0D0
      IF (ABOBJ1.LE.0.D0) ABOBJ1=.000001D0
      IF (ALPHAX.LE.0.D0) ALPHAX=.000001D0
      IF (FDCH.LE.0.D0) FDCH=.000001D0
      IF (FDCHM.LE.0.D0) FDCHM=.000001D0
!     ------------------------------------------------------------------
!                     INITIALIZE INTERNAL PARAMETERS
!     ------------------------------------------------------------------
      INFOG=0
      ITER=0
      JDIR=0
      IOBJ=0
      KOBJ=0
      NDV2=NDV+2
      KCOUNT=0
      NCAL(1)=0
      NCAL(2)=0
      NAC=0
      NFEAS=0
      MSCAL=NSCAL
      CT1=ITRM
      CT1=1./CT1
      DCT=(CTMIN/ABS(CT))**CT1
      DCTL=(CTLMIN/ABS(CTL))**CT1
!	ASH: PHI ORIGINAL = 5.0
      PHI=1000.
      ABOBJ=ABOBJ1
      NCOBJ=0
      CTAM=ABS(CTMIN)
      CTBM=ABS(CTLMIN)
!     CALCULATE NUMBER OF LINEAR CONSTRAINTS, NLNC.
      NLNC=0
      IF (NCON.EQ.0) GO TO 60
      DO 50 I=1,NCON
      IF (ISC(I).GT.0) NLNC=NLNC+1
50    CONTINUE
60    CONTINUE
!     ------------------------------------------------------------------
!          CHECK TO BE SURE THAT SIDE CONSTRAINTS ARE SATISFIED
!     ------------------------------------------------------------------
      IF (NSIDE.EQ.0) GO TO 100
      DO 90 I=1,NDV
      IF (VLB(I).LE.VUB(I)) GO TO 70
      XX=.5*(VLB(I)+VUB(I))
      X(I)=XX
      VLB(I)=XX
      VUB(I)=XX
      !WRITE (6,1120) I
70    CONTINUE
      XX=X(I)-VLB(I)
      IF (XX.GE.0.) GO TO 80
!     LOWER BOUND VIOLATED.
      !WRITE (6,1130) X(I),VLB(I),I
      X(I)=VLB(I)
      GO TO 90
80    CONTINUE
      XX=VUB(I)-X(I)
      IF (XX.GE.0.) GO TO 90
      !WRITE (6,1140) X(I),VUB(I),I
      X(I)=VUB(I)
90    CONTINUE
100   CONTINUE
!     ------------------------------------------------------------------
!                        INITIALIZE SCALING VECTOR, SCAL
!     ------------------------------------------------------------------
      IF (NSCAL.EQ.0) GO TO 140
      IF (NSCAL.LT.0) GO TO 120
      DO 110 I=1,NDV
110   SCAL(I)=1.
      GO TO 140
120   CONTINUE
      DO 130 I=1,NDV
      SI=ABS(SCAL(I))
      IF (SI.LT.1.0E-20) SI=1.0E-5
      SCAL(I)=SI
      SI=1./SI
      X(I)=X(I)*SI
      IF (NSIDE.EQ.0) GO TO 130
      VLB(I)=VLB(I)*SI
      VUB(I)=VUB(I)*SI
130   CONTINUE
140   CONTINUE
!     ------------------------------------------------------------------
!     ***** CALCULATE INITIAL FUNCTION AND CONSTRAINT VALUES  *****
!     ------------------------------------------------------------------
      INFO=1
      NCAL(1)=1
      IGOTO=1
      GO TO 950
150   CONTINUE
      OBJ1=OBJ
!HINDMAN 3/4/2002      IF (DABFUN.LE.0.) DABFUN=.001*ABS(OBJ)
      IF (DABFUN.LE.0.) DABFUN=.00000000001*ABS(OBJ)
      IF (DABFUN.LT.1.0E-10) DABFUN=1.0E-10
      IF (IPRINT.LE.0) GO TO 260
!     ------------------------------------------------------------------
!                    PRINT INITIAL DESIGN INFORMATION
!     ------------------------------------------------------------------
      IF (IPRINT.LE.1) GO TO 220
      IF (NSIDE.EQ.0.AND.NCON.EQ.0) WRITE (6,1300)
      IF (NSIDE.NE.0.OR.NCON.GT.0) WRITE (6,1240)
      WRITE (6,1250) IPRINT,NDV,ITMAX,NCON,NSIDE,ICNDIR,NSCAL,NFDG,LINOBJ,ITRM,N1,N2,N3,N4,N5
      WRITE (6,1270) CT,CTMIN,CTL,CTLMIN,THETA,PHI,DELFUN,DABFUN
      WRITE (6,1260) FDCH,FDCHM,ALPHAX,ABOBJ1
      IF (NSIDE.EQ.0) GO TO 180
      WRITE (6,1280)
      DO 160 I=1,NDV,6
      M1=MIN0(NDV,I+5)
160   WRITE (6,1010) I,(VLB(J),J=I,M1)
      WRITE (6,1290)
      DO 170 I=1,NDV,6
      M1=MIN0(NDV,I+5)
170   WRITE (6,1010) I,(VUB(J),J=I,M1)
180   CONTINUE
      IF (NSCAL.GE.0) GO TO 190
      WRITE (6,1310)
      WRITE (6,1470) (SCAL(I),I=1,NDV)
190   CONTINUE
      IF (NCON.EQ.0) GO TO 220
      IF (NLNC.EQ.0.OR.NLNC.EQ.NCON) GO TO 210
      WRITE (6,1020)
      DO 200 I=1,NCON,15
      M1=MIN0(NCON,I+14)
200   WRITE (6,1030) I,(ISC(J),J=I,M1)
      GO TO 220
210   IF (NLNC.EQ.NCON) WRITE (6,1040)
      IF (NLNC.EQ.0) WRITE (6,1050)
220   CONTINUE
      WRITE (6,1450) OBJ
      WRITE (6,1460)
      DO 230 I=1,NDV
      X1=1.D0
      IF (NSCAL.NE.0) X1=SCAL(I)
230   G1(I)=X(I)*X1
      DO 240 I=1,NDV,6
      M1=MIN0(NDV,I+5)
240   WRITE (6,1010) I,(G1(J),J=I,M1)
      IF (NCON.EQ.0) GO TO 260
      WRITE (6,1480)
      DO 250 I=1,NCON,6
      M1=MIN0(NCON,I+5)
250   WRITE (6,1010) I,(G(J),J=I,M1)
260   CONTINUE
      IF (IPRINT.GT.1) WRITE (6,1370)
!     ------------------------------------------------------------------
!     ********************  BEGIN MINIMIZATION  ************************
!     ------------------------------------------------------------------
270   CONTINUE
      ITER=ITER+1
!	ASH IF THIS IS DROPPED THAT ITER GOES WAY UP !!!
      IF (ABOBJ1.LT..0001) ABOBJ1=.0001
      IF (ABOBJ1.GT..2) ABOBJ1=.2
      IF (ALPHAX.GT.1.) ALPHAX=1.
      IF (ALPHAX.LT..0001) ALPHAX=.0001
!	ASH
      NFEAS=NFEAS+1
!	ASH 10/22/02      IF (NFEAS.GT.10) GO TO 790
      IF (NFEAS.GT.30) GO TO 790
      IF (IPRINT.GT.2) WRITE (6,1320) ITER
      IF (IPRINT.GT.3.AND.NCON.GT.0) WRITE (6,1330) CT,CTL,PHI
      CTA=ABS(CT)
      IF (NCOBJ.EQ.0) GO TO 310
!     ------------------------------------------------------------------
!     NO MOVE ON LAST ITERATION.  DELETE CONSTRAINTS THAT ARE NO
!     LONGER ACTIVE.
!     ------------------------------------------------------------------
      NNAC=NAC
      DO 300 I=1,NNAC
      NIC=IC(I)
      IF (NIC.GT.NCON) NAC=NAC-1
      IF (NIC.GT.NCON) GO TO 300
      CT1=CT
      IF (ISC(NIC).GT.0) CT1=CTL
      IF (G(NIC).GT.CT1) GO TO 300
      NAC=NAC-1
      IF (I.EQ.NNAC) GO TO 300
      IP1=I+1
      DO 290 K=IP1,NNAC
      II=K-1
      DO 280 J=1,NDV2
280   A(J,II)=A(J,K)
290   IC(II)=IC(K)
300   CONTINUE
      GO TO 400
310   CONTINUE
      IF (MSCAL.LT.NSCAL.OR.NSCAL.EQ.0) GO TO 330
      IF (NSCAL.LT.0.AND.KCOUNT.LT.ICNDIR) GO TO 330
      MSCAL=0
      KCOUNT=0
!     ------------------------------------------------------------------
!                          SCALE VARIABLES
!     ------------------------------------------------------------------
      DO 320 I=1,NDV
      SI=SCAL(I)
      XI=SI*X(I)
      SIB=SI
      IF (NSCAL.GT.0) SI=ABS(XI)
      IF (SI.LT.1.0E-10) GO TO 320
      SCAL(I)=SI
      SI=1./SI
      X(I)=XI*SI
      IF (NSIDE.EQ.0) GO TO 320
      VLB(I)=SIB*SI*VLB(I)
      VUB(I)=SIB*SI*VUB(I)
320   CONTINUE
      IF (IPRINT.LT.4.OR.(NSCAL.LT.0.AND.ITER.GT.1)) GO TO 330
      WRITE (6,1340)
      WRITE (6,1470) (SCAL(I),I=1,NDV)
330   CONTINUE
      MSCAL=MSCAL+1
      NAC=0
!     ------------------------------------------------------------------
!          OBTAIN GRADIENTS OF OBJECTIVE AND ACTIVE CONSTRAINTS
!     ------------------------------------------------------------------
      INFO=2
      NCAL(2)=NCAL(2)+1
      IF (NFDG.NE.1) GO TO 350
      IGOTO=2
      GO TO 950
350   CONTINUE
      JGOTO=0
360   CONTINUE
      CALL CNMN01 (JGOTO,X,DF,G,ISC,IC,A,G1,VLB,VUB,SCAL,C,NCAL,DX,DX1,FI,XI,III,N1,N2,N3,N4)
      IGOTO=3
      IF (JGOTO.GT.0) GO TO 950
370   CONTINUE
      INFO=1
      IF (NAC.GE.N3) GO TO 790
      IF (NSCAL.EQ.0.OR.NFDG.EQ.0) GO TO 400
!     ------------------------------------------------------------------
!                              SCALE GRADIENTS
!     ------------------------------------------------------------------
!     SCALE GRADIENT OF OBJECTIVE FUNCTION.
      DO 380 I=1,NDV
380   DF(I)=DF(I)*SCAL(I)
      IF (NFDG.EQ.2.OR.NAC.EQ.0) GO TO 400
!     SCALE GRADIENTS OF ACTIVE CONSTRAINTS.
      DO 390 J=1,NDV
      SCJ=SCAL(J)
      DO 390 I=1,NAC
390   A(J,I)=A(J,I)*SCJ
400   CONTINUE
      IF (IPRINT.LT.3.OR.NCON.EQ.0) GO TO 450
!     ------------------------------------------------------------------
!                                   PRINT
!     ------------------------------------------------------------------
!     PRINT ACTIVE AND VIOLATED CONSTRAINT NUMBERS.
      M1=0
      M2=N3
      IF (NAC.EQ.0) GO TO 430
      DO 420 I=1,NAC
      J=IC(I)
      IF (J.GT.NCON) GO TO 420
      GI=G(J)
      C1=CTAM
      IF (ISC(J).GT.0) C1=CTBM
      GI=GI-C1
      IF (GI.GT.0.) GO TO 410
!     ACTIVE CONSTRAINT.
      M1=M1+1
      MS1(M1)=J
      GO TO 420
410   M2=M2+1
!     VIOLATED CONSTRAINT.
      MS1(M2)=J
420   CONTINUE
430   M3=M2-N3
      WRITE (6,1060) M1
      IF (M1.EQ.0) GO TO 440
      WRITE (6,1070)
      WRITE (6,1490) (MS1(I),I=1,M1)
440   WRITE (6,1080) M3
      IF (M3.EQ.0) GO TO 450
      WRITE (6,1070)
      M3=N3+1
      WRITE (6,1490) (MS1(I),I=M3,M2)
450   CONTINUE
!     ------------------------------------------------------------------
!            CALCULATE GRADIENTS OF ACTIVE SIDE CONSTRAINTS
!     ------------------------------------------------------------------
      IF (NSIDE.EQ.0) GO TO 510
      MCN1=NCON
      M1=0
      DO 490 I=1,NDV
!     LOWER BOUND.
      XI=X(I)
      XID=VLB(I)
      X12=ABS(XID)
      IF (X12.LT.1.) X12=1.
      GI=(XID-XI)/X12
      IF (GI.LT.-1.0E-6) GO TO 470
      M1=M1+1
      MS1(M1)=-I
      NAC=NAC+1
      IF (NAC.GE.N3) GO TO 790
      MCN1=MCN1+1
      DO 460 J=1,NDV
460   A(J,NAC)=0.
      A(I,NAC)=-1.
      IC(NAC)=MCN1
      G(MCN1)=GI
      ISC(MCN1)=1
!     UPPER BOUND.
470   XID=VUB(I)
      X12=ABS(XID)
      IF (X12.LT.1.) X12=1.
      GI=(XI-XID)/X12
      IF (GI.LT.-1.0E-6) GO TO 490
      M1=M1+1
      MS1(M1)=I
      NAC=NAC+1
      IF (NAC.GE.N3) GO TO 790
      MCN1=MCN1+1
      DO 480 J=1,NDV
480   A(J,NAC)=0.
      A(I,NAC)=1.
      IC(NAC)=MCN1
      G(MCN1)=GI
      ISC(MCN1)=1
490   CONTINUE
!     ------------------------------------------------------------------
!                                  PRINT
!     ------------------------------------------------------------------
!     PRINT ACTIVE SIDE CONSTRAINT NUMBERS.
      IF (IPRINT.LT.3) GO TO 510
      WRITE (6,1090) M1
      IF (M1.EQ.0) GO TO 510
      WRITE (6,1100)
      DO 500 I=1,M1,15
      M2=MIN0(M1,I+14)
500   WRITE (6,1490) (MS1(J),J=1,M2)
510   CONTINUE
!     PRINT GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS.
      IF (IPRINT.LT.4) GO TO 550
      WRITE (6,1350)
      DO 520 I=1,NDV,6
      M1=MIN0(NDV,I+5)
520   WRITE (6,1010) I,(DF(J),J=I,M1)
      IF (NAC.EQ.0) GO TO 550
      WRITE (6,1360)
      DO 540 I=1,NAC
      M1=IC(I)
      M2=M1-NCON
      M3=0
      IF (M2.GT.0) M3=IABS(MS1(M2))
      IF (M2.LE.0) WRITE (6,990) M1
      IF (M2.GT.0) WRITE (6,1000) M3
      DO 530 K=1,NDV,6
      M1=MIN0(NDV,K+5)
530   WRITE (6,1010) K,(A(J,I),J=K,M1)
540   WRITE (6,1370)
550   CONTINUE
!     ------------------------------------------------------------------
!     ******************  DETERMINE SEARCH DIRECTION *******************
!     ------------------------------------------------------------------
      ALP=1.0E+20
      IF (NAC.GT.0) GO TO 560
!     ------------------------------------------------------------------
!                        UNCONSTRAINED FUNCTION
!     ------------------------------------------------------------------
!     FIND DIRECTION OF STEEPEST DESCENT OR CONJUGATE DIRECTION.
      NVC=0
      NFEAS=0
      KCOUNT=KCOUNT+1
!     IF KCOUNT.GT.ICNDIR  RESTART CONJUGATE DIRECTION ALGORITHM.
      IF (KCOUNT.GT.ICNDIR.OR.IOBJ.EQ.2) KCOUNT=1
      IF (KCOUNT.EQ.1) JDIR=0
!     IF JDIR = 0 FIND DIRECTION OF STEEPEST DESCENT.
      CALL CNMN02 (JDIR,SLOPE,DFTDF1,DF,S,N1)
      GO TO 610
560   CONTINUE
!     ------------------------------------------------------------------
!                          CONSTRAINED FUNCTION
!     ------------------------------------------------------------------
!     FIND USABLE-FEASIBLE DIRECTION.
      KCOUNT=0
      JDIR=0
      PHI=10.*PHI
      IF (PHI.GT.1000.) PHI=1000.D0
      IF (NFEAS.EQ.1) PHI=5.
!     CALCULATE DIRECTION, S.
      CALL CNMN05 (G,DF,A,S,B,C,SLOPE,PHI,ISC,IC,MS1,NVC,N1,N2,N3,N4,N5)
!     IF THIS DESIGN IS FEASIBLE AND LAST ITERATION WAS INFEASIBLE,
!     SET ABOBJ1=.05 (5 PERCENT).
      IF (NVC.EQ.0.AND.NFEAS.GT.1) ABOBJ1=.05
      IF (NVC.EQ.0) NFEAS=0
      IF (IPRINT.LT.3) GO TO 580
      WRITE (6,1380)
      DO 570 I=1,NAC,6
      M1=MIN0(NAC,I+5)
570   WRITE (6,1010) I,(A(NDV1,J),J=I,M1)
      WRITE (6,1220) S(NDV1)
580   CONTINUE
!     ------------------------------------------------------------------
!     ****************** ONE-DIMENSIONAL SEARCH ************************
!     ------------------------------------------------------------------
      IF (S(NDV1).LT.1.0E-6.AND.NVC.EQ.0) GO TO 690
!     ------------------------------------------------------------------
!                 FIND ALPHA TO OBTAIN A FEASIBLE DESIGN
!     ------------------------------------------------------------------
      IF (NVC.EQ.0) GO TO 610
      ALP=-1.
      DO 600 I=1,NAC
      NCI=IC(I)
      C1=G(NCI)
      CTC=CTAM
      IF (ISC(NCI).GT.0) CTC=CTBM
      IF (C1.LE.CTC) GO TO 600
      ALP1=0.
      DO 590 J=1,NDV
590   ALP1=ALP1+S(J)*A(J,I)
      ALP1=ALP1*A(NDV2,I)
      IF (ABS(ALP1).LT.1.0E-20) GO TO 600
      ALP1=-C1/ALP1
      IF (ALP1.GT.ALP) ALP=ALP1
600   CONTINUE
610   CONTINUE
!     ------------------------------------------------------------------
!                       LIMIT CHANCE TO ABOBJ1*OBJ
!     ------------------------------------------------------------------
      ALP1=1.0E+20
      SI=ABS(OBJ)
      IF (SI.LT..01) SI=.01
      IF (ABS(SLOPE).GT.1.0E-20) ALP1=ABOBJ1*SI/SLOPE
      ALP1=ABS(ALP1)
      IF (NVC.GT.0) ALP1=10.*ALP1
      IF (ALP1.LT.ALP) ALP=ALP1
!     ------------------------------------------------------------------
!                   LIMIT CHANGE IN VARIABLE TO ALPHAX
!     ------------------------------------------------------------------
!	ASH
      ALP11=1.0E+20
      DO 620 I=1,NDV
      SI=ABS(S(I))
      XI=ABS(X(I))
      IF (SI.LT.1.0E-10.OR.XI.LT.0.1) GO TO 620
      ALP1=ALPHAX*XI/SI
      IF (ALP1.LT.ALP11) ALP11=ALP1
620   CONTINUE
      IF (NVC.GT.0) ALP11=10.*ALP11
      IF (ALP11.LT.ALP) ALP=ALP11
      IF (ALP.GT.1.0E+20) ALP=1.0E+20
      IF (ALP.LE.1.0E-20) ALP=1.0E-20
      IF (IPRINT.LT.3) GO TO 640
      WRITE (6,1390)
      DO 630 I=1,NDV,6
      M1=MIN0(NDV,I+5)
630   WRITE (6,1010) I,(S(J),J=I,M1)
      WRITE (6,1110) SLOPE,ALP
640   CONTINUE
      IF (NCON.GT.0.OR.NSIDE.GT.0) GO TO 660
!     ------------------------------------------------------------------
!           DO ONE-DIMENSIONAL SEARCH FOR UNCONSTRAINED FUNCTION
!     ------------------------------------------------------------------
      JGOTO=0
650   CONTINUE
      CALL CNMN03 (X,S,SLOPE,XI,ALP,FFF,A1,A2,A3,A4,F1,F2,F3,F4,APP,N1,NCAL,KOUNT,JGOTO)
      IGOTO=4
      IF (JGOTO.GT.0) GO TO 950
      JDIR=1
!     PROCEED TO CONVERGENCE CHECK.
      GO TO 680
!     ------------------------------------------------------------------
!       SOLVE ONE-DIMENSIONAL SEARCH PROBLEM FOR CONSTRAINED FUNCTION
!     ------------------------------------------------------------------
660   CONTINUE
      JGOTO=0
670   CONTINUE
      CALL CNMN06 (X,VLB,VUB,G,SCAL,DF,S,G1,G2,CTAM,CTBM,SLOPE,XI,ALP, &
          A2,A3,A4,F1,F2,F3,CV1,CV2,CV3,CV4,ALPCA,ALPFES,ALPLN,ALPMIN, &
          ALPNC,ALPSAV,ALPSID,ALPTOT,ISC,N1,N2,NCAL,NVC,ICOUNT,IGOOD1, &
          IGOOD2,IGOOD3, IGOOD4,IBEST,III,NLNC,JGOTO)
      IGOTO=5
      IF (JGOTO.GT.0) GO TO 950
      IF (NAC.EQ.0) JDIR=1
!     ------------------------------------------------------------------
!     *******************     UPDATE ALPHAX   **************************
!     ------------------------------------------------------------------
680   CONTINUE
690   CONTINUE
      IF (ALP.GT.1.0E+19) ALP=0.
!     UPDATE ALPHAX TO BE AVERAGE OF MAXIMUM CHANGE IN X(I)
!     AND ALHPAX.
      ALP11=0.
      DO 700 I=1,NDV
      SI=ABS(S(I))
      XI=ABS(X(I))
      IF (XI.LT.1.0E-10) GO TO 700
      ALP1=ALP*SI/XI
      IF (ALP1.GT.ALP11) ALP11=ALP1
700   CONTINUE
      ALP11=.5*(ALP11+ALPHAX)
      ALP12=5.*ALPHAX
      IF (ALP11.GT.ALP12) ALP11=ALP12
      ALPHAX=ALP11
      NCOBJ=NCOBJ+1
!	ASH
!     ABSOLUTE CHANGE IN OBJECTIVE.
      OBJD=OBJ1-OBJ
      OBJB=ABS(OBJD)
!	ASH
!	10/21/02 IF (OBJB.LT.1.0E-10) OBJB=0.
      IF (OBJB.LT.1.0E-10) OBJB=0.
      IF (NAC.EQ.0.OR.OBJB.GT.0.) NCOBJ=0
      IF (NCOBJ.GT.1) NCOBJ=0
!     ------------------------------------------------------------------
!                                  PRINT
!     ------------------------------------------------------------------
!     PRINT MOVE PARAMETER, NEW X-VECTOR AND CONSTRAINTS.
      IF (IPRINT.LT.3) GO TO 710
      WRITE (6,1400) ALP
710   IF (IPRINT.LT.2) GO TO 780
      IF (OBJB.GT.0.) GO TO 720
      IF (IPRINT.EQ.2) WRITE (6,1410) ITER,OBJ
      IF (IPRINT.GT.2) WRITE (6,1420) OBJ
      GO TO 740
720   IF (IPRINT.EQ.2) GO TO 730
      WRITE (6,1430) OBJ
      GO TO 740
730   WRITE (6,1440) ITER,OBJ
740   WRITE (6,1460)
      DO 750 I=1,NDV
      FF1=1.
      IF (NSCAL.NE.0) FF1=SCAL(I)
750   G1(I)=FF1*X(I)
      DO 760 I=1,NDV,6
      M1=MIN0(NDV,I+5)
760   WRITE (6,1010) I,(G1(J),J=I,M1)
      IF (NCON.EQ.0) GO TO 780
      WRITE (6,1480)
      DO 770 I=1,NCON,6
      M1=MIN0(NCON,I+5)
770   WRITE (6,1010) I,(G(J),J=I,M1)
780   CONTINUE
!     ------------------------------------------------------------------
!                          CHECK CONVERGENCE
!     ------------------------------------------------------------------
!     STOP IF ITER EQUALS ITMAX.
      IF (ITER.GE.ITMAX) GO TO 790
!     ------------------------------------------------------------------
!                     ABSOLUTE CHANGE IN OBJECTIVE
!     ------------------------------------------------------------------
!	ASH 
      OBJB=ABS(OBJD)
      KOBJ=KOBJ+1
      IF (OBJB.GE.DABFUN.OR.NFEAS.GT.0) KOBJ=0
!     ------------------------------------------------------------------
!                     RELATIVE CHANGE IN OBJECTIVE
!     ------------------------------------------------------------------
!	ASH
!	10/21/02 IF (ABS(OBJ1).GT.1.0E-10) OBJD=OBJD/ABS(OBJ1)
      IF (ABS(OBJ1).GT.1.0E-10 ) OBJD=OBJD/ABS(OBJ1)
      ABOBJ1=.5*(ABS(ABOBJ)+ABS(OBJD))
      ABOBJ=ABS(OBJD)
      IOBJ=IOBJ+1
      IF (NVC.GT.0.OR.OBJD.GE.DELFUN) IOBJ=0
      IF (IOBJ.GE.ITRM.OR.KOBJ.GE.ITRM) GO TO 790
      OBJ1=OBJ
!     ------------------------------------------------------------------
!           REDUCE CT IF OBJECTIVE FUNCTION IS CHANGING SLOWLY
!     ------------------------------------------------------------------
      IF (IOBJ.LT.1.OR.NAC.EQ.0) GO TO 270
      CT=DCT*CT
      CTL=CTL*DCTL
      IF(ABS(CT).LT.CTMIN) CT=-CTMIN
      IF(ABS(CTL).LT.CTLMIN) CTL=-CTLMIN
!     ------------------------------------------------------------------
!                     CHECK FOR UNBOUNDED SOLUTION
!     ------------------------------------------------------------------
!     STOP IF OBJ IS LESS THAN -1.0E+38.
      IF (OBJ.GT.-1.0E+38) GO TO 270
      WRITE (6,980)
790   CONTINUE
      IF (NAC.GE.N3) WRITE (6,1500)
!     ------------------------------------------------------------------
!     ****************  FINAL FUNCTION INFORMATION  ********************
!     ------------------------------------------------------------------
      IF (NSCAL.EQ.0) GO TO 820
!     UN-SCALE THE DESIGN VARIABLES.
      DO 810 I=1,NDV
      XI=SCAL(I)
      IF (NSIDE.EQ.0) GO TO 810
      VLB(I)=XI*VLB(I)
      VUB(I)=XI*VUB(I)
810   X(I)=XI*X(I)
!     ------------------------------------------------------------------
!                           PRINT FINAL RESULTS
!     ------------------------------------------------------------------
820   IF (IPRINT.EQ.0.OR.NAC.GE.N3) GO TO 940
      WRITE (6,1510)
      WRITE (6,1430) OBJ
      WRITE (6,1460)
      DO 830 I=1,NDV,6
      M1=MIN0(NDV,I+5)
830   WRITE (6,1010) I,(X(J),J=I,M1)
      IF (NCON.EQ.0) GO TO 890
      WRITE (6,1480)
      DO 840 I=1,NCON,6
      M1=MIN0(NCON,I+5)
840   WRITE (6,1010) I,(G(J),J=I,M1)
!     DETERMINE WHICH CONSTRAINTS ARE ACTIVE AND PRINT.
      NAC=0
      NVC=0
      DO 860 I=1,NCON
      CTA=CTAM
      IF (ISC(I).GT.0) CTA=CTBM
      GI=G(I)
      IF (GI.GT.CTA) GO TO 850
      IF (GI.LT.CT.AND.ISC(I).EQ.0) GO TO 860
      IF (GI.LT.CTL.AND.ISC(I).GT.0) GO TO 860
      NAC=NAC+1
      IC(NAC)=I
      GO TO 860
850   NVC=NVC+1
      MS1(NVC)=I
860   CONTINUE
      WRITE (6,1060) NAC
      IF (NAC.EQ.0) GO TO 870
      WRITE (6,1070)
      WRITE (6,1490) (IC(J),J=1,NAC)
870   WRITE (6,1080) NVC
      IF (NVC.EQ.0) GO TO 880
      WRITE (6,1070)
      WRITE (6,1490) (MS1(J),J=1,NVC)
880   CONTINUE
890   CONTINUE
      IF (NSIDE.EQ.0) GO TO 920
!     DETERMINE WHICH SIDE CONSTRAINTS ARE ACTIVE AND PRINT.
      NAC=0
      DO 910 I=1,NDV
      XI=X(I)
      XID=VLB(I)
      X12=ABS(XID)
      IF (X12.LT.1.) X12=1.
      GI=(XID-XI)/X12
      IF (GI.LT.-1.0E-6) GO TO 900
      NAC=NAC+1
      MS1(NAC)=-I
900   XID=VUB(I)
      X12=ABS(XID)
      IF (X12.LT.1.) X12=1.
      GI=(XI-XID)/X12
      IF (GI.LT.-1.0E-6) GO TO 910
      NAC=NAC+1
      MS1(NAC)=I
910   CONTINUE
      WRITE (6,1090) NAC
      IF (NAC.EQ.0) GO TO 920
      WRITE (6,1100)
      WRITE (6,1490) (MS1(J),J=1,NAC)
920   CONTINUE
      WRITE (6,1150)
!	ASH
      IF (ITER.GE.ITMAX) WRITE (6,1160)
      IF (NFEAS.GE.10) WRITE (6,1170)
      IF (IOBJ.GE.ITRM) WRITE (6,1190) ITRM
      IF (KOBJ.GE.ITRM) WRITE (6,1200) ITRM
      WRITE (6,1210) ITER
      WRITE (6,1520) NCAL(1)
      IF (NCON.GT.0) WRITE (6,1530) NCAL(1)
      IF (NFDG.NE.0) WRITE (6,1540) NCAL(2)
      IF (NCON.GT.0.AND.NFDG.EQ.1) WRITE (6,1550) NCAL(2)
!     ------------------------------------------------------------------
!                   RE-SET BASIC PARAMETERS TO INPUT VALUES
!     ------------------------------------------------------------------
940   ITRM=IDM1
      ITMAX=IDM2
      ICNDIR=IDM3
      DELFUN=DM1
      DABFUN=DM2
      CT=DM3
      CTMIN=DM4
      CTL=DM5
      CTLMIN=DM6
      THETA=DM7
      PHI=DM8
      FDCH=DM9
      FDCHM=DM10
      ABOBJ1=DM11
      ALPHAX=DM12
      IGOTO=0
950   CONTINUE
      IF (NSCAL.EQ.0.OR.IGOTO.EQ.0) RETURN
!     UN-SCALE VARIABLES.
      DO 960 I=1,NDV
      C(I)=X(I)
960   X(I)=X(I)*SCAL(I)
      RETURN
!     ------------------------------------------------------------------
!                                FORMATS
!     ------------------------------------------------------------------

      
      
      
      
970 FORMAT (//t6, 'A COMPLETELY UNCONSTRAINED FUNCTION WITH A LINEAR OBJECTIVE IS SPECIFIED'// &
      t11, 'LINOBJ =', i5/ t11, 'NCON   =', i5/  t11, 'NSIDE  =',i5// t6,&
      'CONTROL RETURNED TO CALLING PROGRAM')
980 FORMAT (//t6, 'CONMIN HAS ACHIEVED A SOLUTION OF OBJ LESS THAN -1.0E+40'/ &
         t6, 'SOLUTION APPEARS TO BE UNBOUNDED'/ t6, 'OPTIMIZATION IS TERMINATED')
990 FORMAT (t6, 'CONSTRAINT NUMBER', i5)   
1000 FORMAT (t6, 'SIDE CONSTRAINT ON VARIABLE', i5)
1010 FORMAT (t4, i5, ')  ', 6E13.5)    
1020 FORMAT (/t6, 'LINEAR CONSTRAINT IDENTIFIERS (ISC)'/  &
    t6, 'NON-ZERO INDICATES LINEAR CONSTRAINT')
1030 FORMAT (t4, i5, ')  ', 15I5)
1040 FORMAT (/t6, 'ALL CONSTRAINTS ARE LINEAR')
1050 FORMAT (/t6, 'ALL CONSTRAINTS ARE NON-LINEAR')
1060 FORMAT (/t6, 'THERE ARE',i5,' ACTIVE CONSTRAINTS') 
1070 FORMAT (t6, 'CONSTRAINT NUMBERS ARE')
1080 FORMAT (/t6, 'THERE ARE', i5, ' VIOLATED CONSTRAINTS')
1090 FORMAT (/t6, 'THERE ARE', i5, ' ACTIVE SIDE CONSTRAINTS')
1100 FORMAT (t6, 'DECISION VARIABLES AT LOWER OR UPPER BOUNDS',  &
    ' (MINUS INDICATES LOWER BOUND)')
1110 FORMAT (/t6, 'ONE-DIMENSIONAL SEARCH'/ t6, 'INITIAL SLOPE =', e12.4,  &
    '  PROPOSED ALPHA =', e12.4)     
1120 FORMAT (//t6, '* * CONMIN DETECTS VLB(I).GT.VUB(I)'/  &
    t6, 'FIX IS SET X(I)=VLB(I)=VUB(I) = .5*(VLB(I)+VUB(I) FOR I =', i5)
1130 FORMAT (//t6, '* * CONMIN DETECTS INITIAL X(I).LT.VLB(I)'/ t6,   &
    'X(I) =', e12.4, '  VLB(I) =', e12.4/ t6,   &
    'X(I) IS SET EQUAL TO VLB(I) FOR I =',i5)
1140 FORMAT (//t6, '* * CONMIN DETECTS INITIAL X(I).GT.VUB(I)'/t6,   &
    'X(I) =', e12.4, '  VUB(I) =', e12.4/ t6,   &
    'X(I) IS SET EQUAL TO VUB(I) FOR I =',i5)
1150 FORMAT (/t6, 'TERMINATION CRITERION')
1160 FORMAT (t11, 'ITER EQUALS ITMAX')
1170 FORMAT (t11,  &
    'NFEASCT CONSECUTIVE ITERATIONS FAILED TO PRODUCE A FEASIBLE DESIGN')     
1190 FORMAT (t11, 'ABS(1-OBJ(I-1)/OBJ(I)) LESS THAN DELFUN FOR', i3,  &
    ' ITERATIONS')
1200 FORMAT (t11,'ABS(OBJ(I)-OBJ(I-1))   LESS THAN DABFUN FOR',i3,  &
    ' ITERATIONS')
1210 FORMAT (/t6, 'NUMBER OF ITERATIONS =', i5)
1220 FORMAT (/t6, 'CONSTRAINT PARAMETER, BETA =', e14.5)
1230 FORMAT (// t13, 27('* ')/ t13, '*', t65, '*'/ t13,'*',t34,  &
    'C O N M I N', t65, '*'/ t13,'*', t65, '*'/ t13,'*',t29,  &
    ' FORTRAN PROGRAM FOR ', t65, '*'/ t13,'*',t65,'*'/ t13,'*',t23,  &
    'CONSTRAINED FUNCTION MINIMIZATION', t65, '*'/ t13, '*',t65, '*'/  &
    t13, 27('* '))
1240 FORMAT (//t6, 'CONSTRAINED FUNCTION MINIMIZATION'//   &
    t6, 'CONTROL PARAMETERS')
1250 FORMAT (/t6, 'IPRINT  NDV    ITMAX    NCON    NSIDE  ICNDIR   NSCAL NFDG'/ &
    8I8//t6, 'LINOBJ  ITRM     N1      N2      N3      N4      N5'/ 8I8)
1260 FORMAT (/t10, 'FDCH', t26, 'FDCHM', t42, 'ALPHAX', t58, 'ABOBJ1'/  &
    ' ', 4('  ', e14.5))
1270 FORMAT (/t10, 'CT', t26, 'CTMIN', t42, 'CTL', t58, 'CTLMIN'/  &
    ' ', 4('  ', e14.5)//  &
    t10, 'THETA', t26, 'PHI', t42, 'DELFUN', t58, 'DABFUN'/  &
    ' ',4('  ', e14.5))
1280 FORMAT (/t6, 'LOWER BOUNDS ON DECISION VARIABLES (VLB)')
1290 FORMAT (/t6, 'UPPER BOUNDS ON DECISION VARIABLES (VUB)')
1300 FORMAT (//t6, 'UNCONSTRAINED FUNCTION MINIMIZATION'//t6,   &
    'CONTROL PARAMETERS')
1310 FORMAT (/t6, 'SCALING VECTOR (SCAL)')
1320 FORMAT (//t6, 'BEGIN ITERATION NUMBER',i5)
1330 FORMAT (/t6, 'CT =', e14.5, '     CTL =', e14.5, '     PHI =', e14.5)
1340 FORMAT (/t6, 'NEW SCALING VECTOR (SCAL)')
1350 FORMAT (/t6, 'GRADIENT OF OBJ')
1360 FORMAT (/t6, 'GRADIENTS OF ACTIVE AND VIOLATED CONSTRAINTS')
1370 FORMAT (' ')
1380 FORMAT (/t6, 'PUSH-OFF FACTORS, (THETA(I), I=1,NAC)')
1390 FORMAT (/t6, 'SEARCH DIRECTION (S-VECTOR)')
1400 FORMAT (/t6, 'CALCULATED ALPHA =', e14.5)
1410 FORMAT (//t6, 'ITER =', i5, '     OBJ =', e14.5, '     NO CHANGE IN OBJ')
1420 FORMAT (/t6, 'OBJ =', e15.6, '     NO CHANGE ON OBJ')
1430 FORMAT (/t6, 'OBJ =', e15.6)
1440 FORMAT (//t6, 'ITER =', i5, '     OBJ =',e14.5)
1450 FORMAT (//t6, 'INITIAL FUNCTION INFORMATION'// t6, 'OBJ =', e15.6)
1460 FORMAT (/t6, 'DECISION VARIABLES (X-VECTOR)')
1470 FORMAT (t4, 7E13.4)
1480 FORMAT (/t6, 'CONSTRAINT VALUES (G-VECTOR)')
1490 FORMAT (t6,15I5)
1500 FORMAT (/t6,  'THE NUMBER OF ACTIVE AND VIOLATED CONSTRAINTS EXCEEDS N3-1.'/  &
    t6, 'DIMENSIONED SIZE OF MATRICES A AND B AND VECTOR IC IS INSUFFICIENT'/  &
    t6, 'OPTIMIZATION TERMINATED AND CONTROL RETURNED TO MAIN PROGRAM.')
1510 FORMAT (///'    FINAL OPTIMIZATION INFORMATION')
1520 FORMAT (/t6, 'OBJECTIVE FUNCTION WAS EVALUATED        ', i5, '  TIMES')
1530 FORMAT (/t6, 'CONSTRAINT FUNCTIONS WERE EVALUATED', i10, '  TIMES')
1540 FORMAT (/t6, 'GRADIENT OF OBJECTIVE WAS CALCULATED', i9, '  TIMES')
1550 FORMAT (/t6, 'GRADIENTS OF CONSTRAINTS WERE CALCULATED', i5, '  TIMES')

     

END SUBROUTINE CONMIN




SUBROUTINE CNMN01 (JGOTO,X,DF,G,ISC,IC,A,G1,VLB,VUB,SCAL,C,NCAL,DX, &
    DX1,FI,XI,III,N1,N2,N3,N4)
USE CNMN1
IMPLICIT NONE !DOUBLE (A-H,O-Z)
!      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,ALPHAX
!     *,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,IT
!     *RM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER
INTEGER, INTENT(INOUT) :: JGOTO, ISC(N2), IC(N3), III, N1, N2, N3, N4, NCAL(2)
DOUBLE PRECISION, INTENT(INOUT) :: X(N1), DF(N1), G(N2), VLB(N1), VUB(N1), &
    A(N1,N3), G1(N2), SCAL(N1), C(N4), DX, DX1, FI, XI

!LOCAL VARIABLES
DOUBLE PRECISION :: FDCH1, X1
INTEGER :: I, I1, INF, J


!     ROUTINE TO CALCULATE GRADIENT INFORMATION BY FINITE DIFFERENCE.
!     BY G. N. VANDERPLAATS                         JUNE, 1972.
!     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
      IF (JGOTO.EQ.1) GO TO 10
      IF (JGOTO.EQ.2) GO TO 70
      INFOG=0
      INF=INFO
      NAC=0
      IF (LINOBJ.NE.0.AND.ITER.GT.1) GO TO 10
!     ------------------------------------------------------------------
!                    GRADIENT OF LINEAR OBJECTIVE
!     ------------------------------------------------------------------
      IF (NFDG.EQ.2) JGOTO=1
      IF (NFDG.EQ.2) RETURN
10    CONTINUE
      JGOTO=0
      IF (NFDG.EQ.2.AND.NCON.EQ.0) RETURN
      IF (NCON.EQ.0) GO TO 40
!     ------------------------------------------------------------------
!       * * * DETERMINE WHICH CONSTRAINTS ARE ACTIVE OR VIOLATED * * *
!     ------------------------------------------------------------------
      DO 20 I=1,NCON
      IF (G(I).LT.CT) GO TO 20
      IF (ISC(I).GT.0.AND.G(I).LT.CTL) GO TO 20
      NAC=NAC+1
      IF (NAC.GE.N3) RETURN
      IC(NAC)=I
20    CONTINUE
      IF (NFDG.EQ.2.AND.NAC.EQ.0) RETURN
      IF ((LINOBJ.GT.0.AND.ITER.GT.1).AND.NAC.EQ.0) RETURN
!     ------------------------------------------------------------------
!                  STORE VALUES OF CONSTRAINTS IN G1
!     ------------------------------------------------------------------
      DO 30 I=1,NCON
30    G1(I)=G(I)
40    CONTINUE
      JGOTO=0
      IF (NAC.EQ.0.AND.NFDG.EQ.2) RETURN
!     ------------------------------------------------------------------
!                            CALCULATE GRADIENTS
!     ------------------------------------------------------------------
      INFOG=1
      INFO=1
      FI=OBJ
      III=0
50    III=III+1
      XI=X(III)
      DX=FDCH*XI
      DX=ABS(DX)
      FDCH1=FDCHM
      IF (NSCAL.NE.0) FDCH1=FDCHM/SCAL(III)
      IF (DX.LT.FDCH1) DX=FDCH1
      X1=XI+DX
      IF (NSIDE.EQ.0) GO TO 60
      IF (X1.GT.VUB(III)) DX=-DX
60    DX1=1./DX
      X(III)=XI+DX
      NCAL(1)=NCAL(1)+1
!     ------------------------------------------------------------------
!                         FUNCTION EVALUATION
!     ------------------------------------------------------------------
      JGOTO=2
      RETURN
70    CONTINUE
      X(III)=XI
!	ASH: NFDG = 0
      IF (NFDG.EQ.0) DF(III)=DX1*(OBJ-FI)
      IF (NAC.EQ.0) GO TO 90
!     ------------------------------------------------------------------
!             DETERMINE GRADIENT COMPONENTS OF ACTIVE CONSTRAINTS
!     ------------------------------------------------------------------
      DO 80 J=1,NAC
      I1=IC(J)
80    A(III,J)=DX1*(G(I1)-G1(I1))
90    CONTINUE
      IF (III.LT.NDV) GO TO 50
      INFOG=0
      INFO=INF
      JGOTO=0
      OBJ=FI
      IF (NCON.EQ.0) RETURN
!     ------------------------------------------------------------------
!             STORE CURRENT CONSTRAINT VALUES BACK IN G-VECTOR
!     ------------------------------------------------------------------
      DO 100 I=1,NCON
100   G(I)=G1(I)
      RETURN
    END SUBROUTINE CNMN01
    
    
    

!     ROUTINE TO DETERMINE CONJUGATE DIRECTION VECTOR OR DIRECTION
!     OF STEEPEST DESCENT FOR UNCONSTRAINED FUNCTION MINIMIZATION.
!     BY G. N. VANDERPLAATS                       APRIL, 1972.
!     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
!     NCALC = CALCULATION CONTROL.
!         NCALC = 0,     S = STEEPEST DESCENT.
!         NCALC = 1,     S = CONJUGATE DIRECTION.
!     CONJUGATE DIRECTION IS FOUND BY FLETCHER-REEVES ALGORITHM.
SUBROUTINE CNMN02 (NCALC,SLOPE,DFTDF1,DF,S,N1)
!IMPLICIT DOUBLE (A-H,O-Z)
!      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,ALPHAX, & ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,IT
!     *RM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER
!      DIMENSION DF(N1),S(N1)    
USE CNMN1
IMPLICIT NONE    
INTEGER, INTENT(INOUT) ::  NCALC   , N1
DOUBLE PRECISION, INTENT(INOUT) :: SLOPE, DFTDF1, S(N1), DF(N1)

!LOCAL VARIABLES
DOUBLE PRECISION :: BETA, DFI, DFTDF, S1, S2, SI
INTEGER :: I

!     ------------------------------------------------------------------
!                   CALCULATE NORM OF GRADIENT VECTOR
!     ------------------------------------------------------------------
      DFTDF=0.D0
      DO 10 I=1,NDV
      DFI=DF(I)
10    DFTDF=DFTDF+DFI*DFI
!     ------------------------------------------------------------------
!     **********                FIND DIRECTION S              **********
!     ------------------------------------------------------------------
      IF (NCALC.NE.1) GO TO 30
      IF (DFTDF1.LT.1.0E-20) GO TO 30
!     ------------------------------------------------------------------
!                 FIND FLETCHER-REEVES CONJUGATE DIRECTION
!     ------------------------------------------------------------------
      BETA=DFTDF/DFTDF1
      SLOPE=0.
      DO 20 I=1,NDV
      DFI=DF(I)
      SI=BETA*S(I)-DFI
      SLOPE=SLOPE+SI*DFI
20    S(I)=SI
      GO TO 50
30    CONTINUE
      NCALC=0
!     ------------------------------------------------------------------
!                  CALCULATE DIRECTION OF STEEPEST DESCENT
!     ------------------------------------------------------------------
      DO 40 I=1,NDV
40    S(I)=-DF(I)
      SLOPE=-DFTDF
50    CONTINUE
!     ------------------------------------------------------------------
!                  NORMALIZE S TO MAX ABS VALUE OF UNITY
!     ------------------------------------------------------------------
      S1=0.
      DO 60 I=1,NDV
      S2=ABS(S(I))
      IF (S2.GT.S1) S1=S2
60    CONTINUE
      IF (S1.LT.1.0E-20) S1=1.0E-20
      S1=1./S1
      DFTDF1=DFTDF*S1
      DO 70 I=1,NDV
70    S(I)=S1*S(I)
      SLOPE=S1*SLOPE
      RETURN
END SUBROUTINE CNMN02



!IMPLICIT DOUBLE (A-H,O-Z)
!COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,ALPHAX
!*,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,IT
!*RM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER
!DIMENSION X(N1),S(N1),NCAL(2)
!     ROUTINE TO SOLVE ONE-DIMENSIONAL SEARCH IN UNCONSTRAINED
!     MINIMIZATION USING 2-POINT QUADRATIC INTERPOLATION, 3-POINT
!     CUBIC INTERPOLATION AND 4-POINT CUBIC INTERPOLATION.
!     BY G. N. VANDERPLAATS                         APRIL, 1972.
!     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
!     ALP = PROPOSED MOVE PARAMETER.
!     SLOPE = INITIAL FUNCTION SLOPE = S-TRANSPOSE TIMES DF.
!     SLOPE MUST BE NEGATIVE.
!     OBJ = INITIAL FUNCTION VALUE.
SUBROUTINE CNMN03 (X,S,SLOPE,XI,ALP,FFF,A1,A2,A3,A4,F1,F2,F3,F4,APP,N1,NCAL,KOUNT,JGOTO)
USE CNMN1
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: N1, NCAL(2), JGOTO, KOUNT
DOUBLE PRECISION, INTENT(INOUT) :: X(N1), S(N1), XI, ALP, FFF, A1, A2, A3, A4, F1, &
    F2, F3, F4, APP, SLOPE


!LOCAL VARIABLES
DOUBLE PRECISION :: AA, AB, AB2, AB3, AP, AP1, FF, ZRO
INTEGER :: I, II, JJ

      ZRO=0.D0
      IF (JGOTO.EQ.0) GO TO 10
      GO TO (50,80,110,140,180,220,270), JGOTO
!     ------------------------------------------------------------------
!                     INITIAL INFORMATION  (ALPHA=0)
!     ------------------------------------------------------------------
10    IF (SLOPE.LT.0.) GO TO 20
      ALP=0.D0
      RETURN
20    CONTINUE
      IF (IPRINT.GT.4) WRITE (6,360)
      FFF=OBJ
      AP1=0.D0
      A1=0.D0
      F1=OBJ
      A2=ALP
      A3=0.D0
      F3=0.D0
      AP=A2
      KOUNT=0
!     ------------------------------------------------------------------
!            MOVE A DISTANCE AP*S AND UPDATE FUNCTION VALUE
!     ------------------------------------------------------------------
30    CONTINUE
      KOUNT=KOUNT+1
      DO 40 I=1,NDV
40    X(I)=X(I)+AP*S(I)
      IF (IPRINT.GT.4) WRITE (6,370) AP
      IF (IPRINT.GT.4) WRITE (6,380) (X(I),I=1,NDV)
      NCAL(1)=NCAL(1)+1
      JGOTO=1
      RETURN
50    CONTINUE
      F2=OBJ
      IF (IPRINT.GT.4) WRITE (6,390) F2
      IF (F2.LT.F1) GO TO 120
!     ------------------------------------------------------------------
!                     CHECK FOR ILL-CONDITIONING
!     ------------------------------------------------------------------
      IF (KOUNT.GT.5) GO TO 60
      FF=2.*ABS(F1)
      IF (F2.LT.FF) GO TO 90
      FF=5.*ABS(F1)
      IF (F2.LT.FF) GO TO 60
      A2=.5*A2
      AP=-A2
      ALP=A2
      GO TO 30
60    F3=F2
      A3=A2
      A2=.5*A2
!     ------------------------------------------------------------------
!                 UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
      AP=A2-ALP
      ALP=A2
      DO 70 I=1,NDV
70    X(I)=X(I)+AP*S(I)
      IF (IPRINT.GT.4) WRITE (6,370) A2
      IF (IPRINT.GT.4) WRITE (6,380) (X(I),I=1,NDV)
      NCAL(1)=NCAL(1)+1
      JGOTO=2
      RETURN
80    CONTINUE
      F2=OBJ
      IF (IPRINT.GT.4) WRITE (6,390) F2
!     PROCEED TO CUBIC INTERPOLATION.
      GO TO 160
90    CONTINUE
!     ------------------------------------------------------------------
!     **********        2-POINT QUADRATIC INTERPOLATION       **********
!     ------------------------------------------------------------------
      JJ=1
      II=1
      CALL CNMN04 (II,APP,ZRO,A1,F1,SLOPE,A2,F2,ZRO,ZRO,ZRO,ZRO)
      IF (APP.LT.ZRO.OR.APP.GT.A2) GO TO 120
      F3=F2
      A3=A2
      A2=APP
      JJ=0
!     ------------------------------------------------------------------
!                  UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
      AP=A2-ALP
      ALP=A2
      DO 100 I=1,NDV
100   X(I)=X(I)+AP*S(I)
      IF (IPRINT.GT.4) WRITE (6,370) A2
      IF (IPRINT.GT.4) WRITE (6,380) (X(I),I=1,NDV)
      NCAL(1)=NCAL(1)+1
      JGOTO=3
      RETURN
110   CONTINUE
      F2=OBJ
      IF (IPRINT.GT.4) WRITE (6,390) F2
      GO TO 150
120   A3=2.*A2
!     ------------------------------------------------------------------
!                  UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
      AP=A3-ALP
      ALP=A3
      DO 130 I=1,NDV
130   X(I)=X(I)+AP*S(I)
      IF (IPRINT.GT.4) WRITE (6,370) A3
      IF (IPRINT.GT.4) WRITE (6,380) (X(I),I=1,NDV)
      NCAL(1)=NCAL(1)+1
      JGOTO=4
      RETURN
140   CONTINUE
      F3=OBJ
      IF (IPRINT.GT.4) WRITE (6,390) F3
150   CONTINUE
      IF (F3.LT.F2) GO TO 190
160   CONTINUE
!     ------------------------------------------------------------------
!     **********       3-POINT CUBIC INTERPOLATION      **********
!     ------------------------------------------------------------------
      II=3
      CALL CNMN04 (II,APP,ZRO,A1,F1,SLOPE,A2,F2,A3,F3,ZRO,ZRO)
      IF (APP.LT.ZRO.OR.APP.GT.A3) GO TO 190
!     ------------------------------------------------------------------
!     UPDATE DESIGN VECTOR AND FUNCTION VALUE.
!     ------------------------------------------------------------------
      AP1=APP
      AP=APP-ALP
      ALP=APP
      DO 170 I=1,NDV
170   X(I)=X(I)+AP*S(I)
      IF (IPRINT.GT.4) WRITE (6,370) ALP
      IF (IPRINT.GT.4) WRITE (6,380) (X(I),I=1,NDV)
      NCAL(1)=NCAL(1)+1
      JGOTO=5
      RETURN
180   CONTINUE
      IF (IPRINT.GT.4) WRITE (6,390) OBJ
!     ------------------------------------------------------------------
!                         CHECK CONVERGENCE
!     ------------------------------------------------------------------
      AA=1.-APP/A2
      AB2=ABS(F2)
      AB3=ABS(OBJ)
      AB=AB2
      IF (AB3.GT.AB) AB=AB3
      IF (AB.LT.1.0E-15) AB=1.0E-15
      AB=(AB2-AB3)/AB
      IF (ABS(AB).LT.1.0E-15.AND.ABS(AA).LT..001) GO TO 330
      A4=A3
      F4=F3
      A3=APP
      F3=OBJ
      IF (A3.GT.A2) GO TO 230
      A3=A2
      F3=F2
      A2=APP
      F2=OBJ
      GO TO 230
190   CONTINUE
!     ------------------------------------------------------------------
!     **********        4-POINT CUBIC INTERPOLATION       **********
!     ------------------------------------------------------------------
200   CONTINUE
      A4=2.*A3
!     UPDATE DESIGN VECTOR AND FUNCTION VALUE.
      AP=A4-ALP
      ALP=A4
      DO 210 I=1,NDV
210   X(I)=X(I)+AP*S(I)
      IF (IPRINT.GT.4) WRITE (6,370) ALP
      IF (IPRINT.GT.4) WRITE (6,380) (X(I),I=1,NDV)
      NCAL(1)=NCAL(1)+1
      JGOTO=6
      RETURN
220   CONTINUE
      F4=OBJ
      IF (IPRINT.GT.4) WRITE (6,390) F4
      IF (F4.GT.F3) GO TO 230
      A1=A2
      F1=F2
      A2=A3
      F2=F3
      A3=A4
      F3=F4
      GO TO 200
230   CONTINUE
      II=4
      CALL CNMN04 (II,APP,A1,A1,F1,SLOPE,A2,F2,A3,F3,A4,F4)
      IF (APP.GT.A1) GO TO 250
      AP=A1-ALP
      ALP=A1
      OBJ=F1
      DO 240 I=1,NDV
240   X(I)=X(I)+AP*S(I)
      GO TO 280
250   CONTINUE
!     ------------------------------------------------------------------
!                 UPDATE DESIGN VECTOR AND FUNCTION VALUE
!     ------------------------------------------------------------------
      AP=APP-ALP
      ALP=APP
      DO 260 I=1,NDV
260   X(I)=X(I)+AP*S(I)
      IF (IPRINT.GT.4) WRITE (6,370) ALP
      IF (IPRINT.GT.4) WRITE (6,380) (X(I),I=1,NDV)
      NCAL(1)=NCAL(1)+1
      JGOTO=7
      RETURN
270   CONTINUE
      IF (IPRINT.GT.4) WRITE (6,390) OBJ
280   CONTINUE
!     ------------------------------------------------------------------
!                    CHECK FOR ILL-CONDITIONING
!     ------------------------------------------------------------------
      IF (OBJ.GT.F2.OR.OBJ.GT.F3) GO TO 290
      IF (OBJ.LE.F1) GO TO 330
      AP=A1-ALP
      ALP=A1
      OBJ=F1
      GO TO 310
290   CONTINUE
      IF (F2.LT.F3) GO TO 300
      OBJ=F3
      AP=A3-ALP
      ALP=A3
      GO TO 310
300   OBJ=F2
      AP=A2-ALP
      ALP=A2
310   CONTINUE
!     ------------------------------------------------------------------
!                       UPDATE DESIGN VECTOR
!     ------------------------------------------------------------------
      DO 320 I=1,NDV
320   X(I)=X(I)+AP*S(I)
330   CONTINUE
!     ------------------------------------------------------------------
!                     CHECK FOR MULTIPLE MINIMA
!     ------------------------------------------------------------------
      IF (OBJ.LE.FFF) GO TO 350
!     INITIAL FUNCTION IS MINIMUM.
      DO 340 I=1,NDV
340   X(I)=X(I)-ALP*S(I)
      ALP=0.
      OBJ=FFF
350   CONTINUE
      JGOTO=0
      RETURN
!     ------------------------------------------------------------------
!                                 FORMATS
!     ------------------------------------------------------------------

     
360 FORMAT (///t6,  &
     '* * * UNCONSTRAINED ONE-DIMENSIONAL SEARCH INFORMATION * * *')
370 FORMAT (/t6, 'ALPHA =', e14.5/ t6, 'X-VECTOR')
380 FORMAT (t6, 6E13.5)
390 FORMAT (/t6, 'OBJ =', e14.5)      
      
      
END SUBROUTINE CNMN03



!      SUBROUTINE CNMN04 (II,XBAR,EPS,X1,Y1,SLOPE,X2,Y2,X3,Y3,X4,Y4)
!      IMPLICIT DOUBLE (A-H,O-Z)
!     ROUTINE TO FIND FIRST XBAR.GE.EPS CORRESPONDING TO A MINIMUM
!     OF A ONE-DIMENSIONAL REAL FUNCTION BY POLYNOMIEL INTERPOLATION.
!     BY G. N. VANDERPLAATS                          APRIL, 1972.
!     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
!
!     II = CALCULATION CONTROL.
!          1:  2-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, SLOPE,
!              X2 AND Y2.
!          2:  3-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, X2, Y2,
!              X3 AND Y3.
!          3:  3-POINT CUBIC INTERPOLATION, GIVEN X1, Y1, SLOPE, X2, Y2,
!              X3 AND Y3.
!          4:  4-POINT CUBIC INTERPOLATION, GIVEN X1, Y1, X2, Y2, X3,
!              Y3, X4 AND Y4.
!     EPS MAY BE NEGATIVE.
!     IF REQUIRED MINIMUM ON Y DOES NOT EXITS, OR THE FUNCTION IS
!     ILL-CONDITIONED, XBAR = EPS-1.0 WILL BE RETURNED AS AN ERROR
!     INDICATOR.
!     IF DESIRED INTERPOLATION IS ILL-CONDITIONED, A LOWER ORDER
!     INTERPOLATION, CONSISTANT WITH INPUT DATA, WILL BE ATTEMPTED,
!     AND II WILL BE CHANGED ACCORDINGLY.
SUBROUTINE CNMN04 (II,XBAR,EPS,X1,Y1,SLOPE,X2,Y2,X3,Y3,X4,Y4)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: II
DOUBLE PRECISION, INTENT(INOUT) :: XBAR, EPS, X1, Y1, SLOPE, X2, Y2, X3, Y3, X4, Y4

!LOCAL VARIABLES
DOUBLE PRECISION :: AA, BAC, BB, CC, DNOM, DX, Q1, Q2, Q3, Q4, Q5, Q6, QQ, X11, &
    X111, X21, X22, X222, X31, X32, X33, X41, X42, X44, XBAR1
INTEGER :: NSLOP

      XBAR1=EPS-1.D0
      XBAR=XBAR1
      X21=X2-X1
      IF(ABS(X21).LT.1.0E-20) RETURN
      NSLOP=MOD(II,2)
      GO TO (10,20,40,50), II
10    CONTINUE
!     ------------------------------------------------------------------
!                 II=1: 2-POINT QUADRATIC INTERPOLATION
!     ------------------------------------------------------------------
      II=1
      DX=X1-X2
      IF (ABS(DX).LT.1.0E-20) RETURN
      AA=(SLOPE+(Y2-Y1)/DX)/DX
      IF (AA.LT.1.0E-20) RETURN
      BB=SLOPE-2.*AA*X1
      XBAR=-.5*BB/AA
      IF (XBAR.LT.EPS) XBAR=XBAR1
      RETURN
20    CONTINUE
!     ------------------------------------------------------------------
!                 II=2: 3-POINT QUADRATIC INTERPOLATION
!     ------------------------------------------------------------------
      II=2
      X21=X2-X1
      X31=X3-X1
      X32=X3-X2
      QQ=X21*X31*X32
      IF (ABS(QQ).LT.1.0E-20) RETURN
      AA=(Y1*X32-Y2*X31+Y3*X21)/QQ
      IF (AA.LT.1.0E-20) GO TO 30
      BB=(Y2-Y1)/X21-AA*(X1+X2)
      XBAR=-0.5D0*BB/AA
      IF (XBAR.LT.EPS) XBAR=XBAR1
      RETURN
30    CONTINUE
      IF (NSLOP.EQ.0) RETURN
      GO TO 10
40    CONTINUE
!     ------------------------------------------------------------------
!                   II=3: 3-POINT CUBIC INTERPOLATION
!     ------------------------------------------------------------------
      II=3
      X21=X2-X1
      X31=X3-X1
      X32=X3-X2
      QQ=X21*X31*X32
      IF (ABS(QQ).LT.1.0E-20) RETURN
      X11=X1*X1
      DNOM=X2*X2*X31-X11*X32-X3*X3*X21
      IF (ABS(DNOM).LT.1.0E-20) GO TO 20
      AA=((X31*X31*(Y2-Y1)-X21*X21*(Y3-Y1))/(X31*X21)-SLOPE*X32)/DNOM
      IF (ABS(AA).LT.1.0E-20) GO TO 20
      BB=((Y2-Y1)/X21-SLOPE-AA*(X2*X2+X1*X2-2.*X11))/X21
      CC=SLOPE-3.*AA*X11-2.*BB*X1
      BAC=BB*BB-3.*AA*CC
      IF (BAC.LT.0.) GO TO 20
      BAC=SQRT(BAC)
      XBAR=(BAC-BB)/(3.*AA)
      IF (XBAR.LT.EPS) XBAR=EPS
      RETURN
50    CONTINUE
!     ------------------------------------------------------------------
!                    II=4: 4-POINT CUBIC INTERPOLATION
!     ------------------------------------------------------------------
      X21=X2-X1
      X31=X3-X1
      X41=X4-X1
      X32=X3-X2
      X42=X4-X2
      X11=X1*X1
      X22=X2*X2
      X33=X3*X3
      X44=X4*X4
      X111=X1*X11
      X222=X2*X22
      Q2=X31*X21*X32
      IF (ABS(Q2).LT.1.0E-30) RETURN
      Q1=X111*X32-X222*X31+X3*X33*X21
      Q4=X111*X42-X222*X41+X4*X44*X21
      Q5=X41*X21*X42
      DNOM=Q2*Q4-Q1*Q5
      IF (ABS(DNOM).LT.1.0E-30) GO TO 60
      Q3=Y3*X21-Y2*X31+Y1*X32
      Q6=Y4*X21-Y2*X41+Y1*X42
      AA=(Q2*Q6-Q3*Q5)/DNOM
      BB=(Q3-Q1*AA)/Q2
      CC=(Y2-Y1-AA*(X222-X111))/X21-BB*(X1+X2)
      BAC=BB*BB-3.*AA*CC
      IF (ABS(AA).LT.1.0E-20.OR.BAC.LT.0.) GO TO 60
      BAC=SQRT(BAC)
      XBAR=(BAC-BB)/(3.*AA)
      IF (XBAR.LT.EPS) XBAR=XBAR1
      RETURN
60    CONTINUE
      IF (NSLOP.EQ.1) GO TO 40
      GO TO 20
END SUBROUTINE CNMN04






!      IMPLICIT DOUBLE (A-H,O-Z)
!      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,ALPHAX
!     *,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,IT
!     *RM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER
!      DIMENSION DF(N1),G(N2),ISC(N2),IC(N3),A(N1,N3),S(N1),C(N4),MS1(N5),B(N3,N3)

!     ROUTINE TO SOLVE DIRECTION FINDING PROBLEM IN MODIFIED METHOD OF
!     FEASIBLE DIRECTIONS.
!     BY G. N. VANDERPLAATS                            MAY, 1972.
!     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
!     NORM OF S VECTOR USED HERE IS S-TRANSPOSE TIMES S.LE.1.
!     IF NVC = 0 FIND DIRECTION BY ZOUTENDIJK'S METHOD.  OTHERWISE
!     FIND MODIFIED DIRECTION.
!     ------------------------------------------------------------------
!     ***  NORMALIZE GRADIENTS, CALCULATE THETA'S AND DETERMINE NVC  ***
!     ------------------------------------------------------------------
SUBROUTINE CNMN05 (G,DF,A,S,B,C,SLOPE,PHI,ISC,IC,MS1,NVC,N1,N2,N3,N4,N5)
USE CNMN1
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: N1, N2, N3, N4, N5, ISC(N2), IC(N3), MS1(N5), NVC
DOUBLE PRECISION, INTENT(INOUT) :: G(N2), DF(N1), A(N1,N3), S(N1), B(N3,N3), C(N4), SLOPE, PHI

!LOCAL VARIABLES
DOUBLE PRECISION :: A1, C1, CT1, CT2, CTA, CTAM, CTB, CTBM, CTC, CTD, GG, S1, THMAX, THT
INTEGER :: I, J, K, NAC1, NCI, NCJ, NDB, NDV1, NDV2, NER

      NDV1=NDV+1
      NDV2=NDV+2
      NAC1=NAC+1
      NVC=0
      THMAX=0.D0
      CTA=ABS(CT)
      CT1=1.D0/CTA
      CTAM=ABS(CTMIN)
      CTB=ABS(CTL)
      CT2=1.D0/CTB
      CTBM=ABS(CTLMIN)
      A1=1.D0
      DO 40 I=1,NAC
!     ASH CALCULATE THETA
      NCI=IC(I)
      NCJ=1
      IF (NCI.LE.NCON) NCJ=ISC(NCI)
      C1=G(NCI)
      CTD=CT1
      CTC=CTAM
      IF (NCJ.LE.0) GO TO 10
      CTC=CTBM
      CTD=CT2
10    IF (C1.GT.CTC) NVC=NVC+1
      THT=0.D0
      GG=1.D0+CTD*C1
      IF (NCJ.EQ.0.OR.C1.GT.CTC) THT=THETA*GG*GG
      IF (NCJ.GT.0.AND.C1.GT.CTC) THT=THT-3.*THETA
      IF (THT.GT.50.) THT=50.
      IF (THT.GT.THMAX) THMAX=THT
      A(NDV1,I)=THT
!     ------------------------------------------------------------------
!                    NORMALIZE GRADIENTS OF CONSTRAINTS
!     ------------------------------------------------------------------
      A(NDV2,I)=1.
      IF (NCI.GT.NCON) GO TO 40
      A1=0.D0
      DO 20 J=1,NDV
      A1=A1+A(J,I)**2
20    CONTINUE
      IF (A1.LT.1.0E-20) A1=1.0E-20
      A1=SQRT(A1)
      A(NDV2,I)=A1
      A1=1.D0/A1
      DO 30 J=1,NDV
30    A(J,I)=A1*A(J,I)
40    CONTINUE
!     ------------------------------------------------------------------
!     NORMALIZE GRADIENT OF OBJECTIVE FUNCTION AND STORE IN NAC+1
!     ROW OF A
!     ------------------------------------------------------------------
      A1=0.
      DO 50 I=1,NDV
      A1=A1+DF(I)**2
50    CONTINUE
      IF (A1.LT.1.0E-20) A1=1.0E-20
      A1=SQRT(A1)
      A1=1./A1
      DO 60 I=1,NDV
60    A(I,NAC1)=A1*DF(I)
!     BUILD C VECTOR.
      IF (NVC.GT.0) GO TO 80
!     ------------------------------------------------------------------
!                 BUILD C FOR CLASSICAL METHOD
!     ------------------------------------------------------------------
      NDB=NAC1
      A(NDV1,NDB)=1.
      DO 70 I=1,NDB
70    C(I)=-A(NDV1,I)
      GO TO 110
80    CONTINUE
!     ------------------------------------------------------------------
!                   BUILD C FOR MODIFIED METHOD
!     ------------------------------------------------------------------
      NDB=NAC
      A(NDV1,NAC1)=-PHI
!     ------------------------------------------------------------------
!           SCALE THETA'S SO THAT MAXIMUM THETA IS UNITY
!     ------------------------------------------------------------------
!	ASH THETA    
      IF(THMAX.GT.0.00001) THMAX=1./THMAX
      DO 90 I=1,NDB
      NCI=IC(I)
      C1=CTA
      IF (ISC(NCI).GT.0) C1=CTB
      A(NDV1,I)=A(NDV1,I)*THMAX
90    CONTINUE
      DO 100 I=1,NDB
      C(I)=0.
      DO 100 J=1,NDV1
100   C(I)=C(I)+A(J,I)*A(J,NAC1)
110   CONTINUE
!     ------------------------------------------------------------------
!                      BUILD B MATRIX
!     ------------------------------------------------------------------
      DO 120 I=1,NDB
      DO 120 J=1,NDB
      B(I,J)=0.
      DO 120 K=1,NDV1
120   B(I,J)=B(I,J)-A(K,I)*A(K,J)
!     ------------------------------------------------------------------
!                    SOLVE SPECIAL L. P. PROBLEM
!     ------------------------------------------------------------------
      CALL CNMN08 (NDB,NER,C,MS1,B,N3,N4,N5)
      IF (IPRINT.GT.1.AND.NER.GT.0) WRITE (6,180)
!     CALCULATE RESULTING DIRECTION VECTOR, S.
      SLOPE=0.
!     ------------------------------------------------------------------
!                  USABLE-FEASIBLE DIRECTION
!     ------------------------------------------------------------------
      DO 140 I=1,NDV
      S1=0.
      IF (NVC.GT.0) S1=-A(I,NAC1)
      DO 130 J=1,NDB
130   S1=S1-A(I,J)*C(J)
      SLOPE=SLOPE+S1*DF(I)
140   S(I)=S1
      S(NDV1)=1.
      IF (NVC.GT.0) S(NDV1)=-A(NDV1,NAC1)
      DO 150 J=1,NDB
150   S(NDV1)=S(NDV1)-A(NDV1,J)*C(J)
!     ------------------------------------------------------------------
!                  NORMALIZE S TO MAX ABS OF UNITY
!     ------------------------------------------------------------------
      S1=0.
      DO 160 I=1,NDV
      A1=ABS(S(I))
      IF (A1.GT.S1) S1=A1
160   CONTINUE
      IF (S1.LT.1.0E-10) S1=1.0E-10
      S1=1./S1
      DO 170 I=1,NDV
170   S(I)=S1*S(I)
      SLOPE=S1*SLOPE
      S(NDV1)=S1*S(NDV1)
      RETURN

180 FORMAT (//t6, '* * DIRECTION FINDING PROCESS DID NOT CONVERGE'/t6,   &
    '* * S-VECTOR MAY NOT BE VALID')
      

END SUBROUTINE CNMN05


SUBROUTINE CNMN06 (X,VLB,VUB,G,SCAL,DF,S,G1,G2,CTAM,CTBM,SLOPE,XI, &
    ALP,A2,A3,A4,F1,F2,F3,CV1,CV2,CV3,CV4,ALPCA,ALPFES,ALPLN,ALPMIN,ALPNC,&
    ALPSAV,ALPSID,ALPTOT,ISC,N1,N2,NCAL,NVC,ICOUNT,IGOOD1,IGOOD2,IGOOD3,&
    IGOOD4,IBEST,III,NLNC,JGOTO)
USE CNMN1
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: ISC(N2), N1, N2, NCAL(2), NVC, ICOUNT, IGOOD1, IGOOD2, IGOOD3, &
    IGOOD4, IBEST, III, NLNC, JGOTO
DOUBLE PRECISION, INTENT(INOUT) :: X(N1), VLB(N1), VUB(N1), G(N2), SCAL(N1), DF(N1), &
    S(N1), G1(N2), G2(N2), CTAM, CTBM, SLOPE, XI, ALP, A2, A3, A4, F1, F2, F3, CV1, &
    CV2, CV3, CV4, ALPCA, ALPFES, ALPLN, ALPMIN, ALPNC, ALPSAV, ALPSID, ALPTOT

!LOCAL VARIABLES
DOUBLE PRECISION :: ALPA, ALPB, C1, C2, C3, CC, F4, FI, GI, SI, XI1, XI2, ZRO
INTEGER :: I, II, JBEST, KSID, NVC1

!      IMPLICIT DOUBLE (A-H,O-Z)
!      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,ALPHAX
!     *,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ,ITMAX,IT
!     *RM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER
!      DIMENSION X(N1),VLB(N1),VUB(N1),G(N2),SCAL(N1),DF(N1),S(N1),G1(N2)
!     *,G2(N2),ISC(N2),NCAL(2)
!     ROUTINE TO SOLVE ONE-DIMENSIONAL SEARCH PROBLEM FOR CONSTRAINED
!     FUNCTION MINIMIZATION.
!     BY G. N. VANDERPLAATS                           AUG., 1974.
!     NASA-AMES RESEARCH CENTER, MOFFETT FIELD, CALIF.
!     OBJ = INITIAL AND FINAL FUNCTION VALUE.
!     ALP = MOVE PARAMETER.
!     SLOPE = INITIAL SLOPE.
!
!     ALPSID = MOVE TO SIDE CONSTRAINT.
!     ALPFES = MOVE TO FEASIBLE REGION.
!     ALPNC = MOVE TO NEW NON-LINEAR CONSTRAINT.
!     ALPLN = MOVE TO LINEAR CONSTRAINT.
!     ALPCA = MOVE TO RE-ENCOUNTER CURRENTLY ACTIVE CONSTRAINT.
!     ALPMIN = MOVE TO MINIMIZE FUNCTION.
!     ALPTOT = TOTAL MOVE PARAMETER.
      ZRO=0.D0
      IF (JGOTO.EQ.0) GO TO 10
      GO TO (140,310,520), JGOTO
10    IF (IPRINT.GE.5) WRITE (6,730)
      ALPSAV=ALP
      ICOUNT=0
      ALPTOT=0.D0
!     TOLERANCES.
      CTAM=ABS(CTMIN)
      CTBM=ABS(CTLMIN)
!     PROPOSED MOVE.
20    CONTINUE
!     ------------------------------------------------------------------
!     *****  BEGIN SEARCH OR IMPOSE SIDE CONSTRAINT MODIFICATION  *****
!     ------------------------------------------------------------------
      A2=ALPSAV
      ICOUNT=ICOUNT+1
      ALPSID=1.0E+20
!     INITIAL ALPHA AND OBJ.
      ALP=0.
      F1=OBJ
      KSID=0
      IF (NSIDE.EQ.0) GO TO 70
!     ------------------------------------------------------------------
!     FIND MOVE TO SIDE CONSTRAINT AND INSURE AGAINST VIOLATION OF
!     SIDE CONSTRAINTS
!     ------------------------------------------------------------------
      DO 60 I=1,NDV
      SI=S(I)
      IF (ABS(SI).GT.1.0E-20) GO TO 30
!     ITH COMPONENT OF S IS SMALL.  SET TO ZERO.
      S(I)=0.
      SLOPE=SLOPE-SI*DF(I)
      GO TO 60
30    CONTINUE
      XI=X(I)
      SI=1./SI
      IF (SI.GT.0.) GO TO 40
!     LOWER BOUND.
      XI2=VLB(I)
      XI1=ABS(XI2)
      IF (XI1.LT.1.) XI1=1.
!     CONSTRAINT VALUE.
      GI=(XI2-XI)/XI1
      IF (GI.GT.-1.0E-6) GO TO 50
!     PROPOSED MOVE TO LOWER BOUND.
      ALPA=(XI2-XI)*SI
      IF (ALPA.LT.ALPSID) ALPSID=ALPA
      GO TO 60
40    CONTINUE
!     UPPER BOUND.
      XI2=VUB(I)
      XI1=ABS(XI2)
      IF (XI1.LT.1.) XI1=1.
!     CONSTRAINT VALUE.
      GI=(XI-XI2)/XI1
      IF (GI.GT.-1.0E-6) GO TO 50
!     PROPOSED MOVE TO UPPER BOUND.
      ALPA=(XI2-XI)*SI
      IF (ALPA.LT.ALPSID) ALPSID=ALPA
      GO TO 60
50    CONTINUE
!     MOVE WILL VIOLATE SIDE CONSTRAINT.  SET S(I)=0.
      SLOPE=SLOPE-S(I)*DF(I)
      S(I)=0.
      KSID=KSID+1
60    CONTINUE
!     ALPSID IS UPPER BOUND ON ALPHA.
      IF (A2.GT.ALPSID) A2=ALPSID
70    CONTINUE
!     ------------------------------------------------------------------
!               CHECK ILL-CONDITIONING
!     ------------------------------------------------------------------
      IF (KSID.EQ.NDV.OR.ICOUNT.GT.10) GO TO 710
      IF (NVC.EQ.0.AND.SLOPE.GT.0.) GO TO 710
      ALPFES=-1.
      ALPMIN=-1.
      ALPLN=1.1*ALPSID
      ALPNC=ALPSID
      ALPCA=ALPSID
      IF (NCON.EQ.0) GO TO 90
!     STORE CONSTRAINT VALUES IN G1.
      DO 80 I=1,NCON
      G1(I)=G(I)
80    CONTINUE
90    CONTINUE
!     ------------------------------------------------------------------
!                  MOVE A DISTANCE A2*S
!     ------------------------------------------------------------------
      ALPTOT=ALPTOT+A2
      DO 100 I=1,NDV
      X(I)=X(I)+A2*S(I)
100   CONTINUE
      IF (IPRINT.LT.5) GO TO 130
      WRITE (6,740) A2
      IF (NSCAL.EQ.0) GO TO 120
      DO 110 I=1,NDV
110   G(I)=SCAL(I)*X(I)
      WRITE (6,750) (G(I),I=1,NDV)
      GO TO 130
120   WRITE (6,750) (X(I),I=1,NDV)
!     ------------------------------------------------------------------
!                   UPDATE FUNCTION AND CONSTRAINT VALUES
!     ------------------------------------------------------------------
!                   UPDATE FUNCTION AND CONSTRAINT VALUES
!     ------------------------------------------------------------------
130   NCAL(1)=NCAL(1)+1
      JGOTO=1
      RETURN
140   CONTINUE
      F2=OBJ
      IF (IPRINT.GE.5) WRITE (6,760) F2
      IF (IPRINT.LT.5.OR.NCON.EQ.0) GO TO 150
      WRITE (6,770)
      WRITE (6,750) (G(I),I=1,NCON)
150   CONTINUE
!     ------------------------------------------------------------------
!               IDENTIFY ACCAPTABILITY OF DESIGNS F1 AND F2
!     ------------------------------------------------------------------
!     IGOOD = 0 IS ACCAPTABLE.
!     CV = MAXIMUM CONSTRAINT VIOLATION.
      IGOOD1=0
      IGOOD2=0
      CV1=0.
      CV2=0.
      NVC1=0
      IF (NCON.EQ.0) GO TO 170
      DO 160 I=1,NCON
      CC=CTAM
      IF (ISC(I).GT.0) CC=CTBM
      C1=G1(I)-CC
      C2=G(I)-CC
      IF (C2.GT.0.) NVC1=NVC1+1
      IF (C1.GT.CV1) CV1=C1
      IF (C2.GT.CV2) CV2=C2
160   CONTINUE
      IF (CV1.GT.0.) IGOOD1=1
      IF (CV2.GT.0.) IGOOD2=1
170   CONTINUE
      ALP=A2
      OBJ=F2
!     ------------------------------------------------------------------
!     IF F2 VIOLATES FEWER CONSTRAINTS THAN F1 BUT STILL HAS CONSTRAINT
!     VIOLATIONS RETURN
!     ------------------------------------------------------------------
      IF (NVC1.LT.NVC.AND.NVC1.GT.0) GO TO 710
!     ------------------------------------------------------------------
!             IDENTIFY BEST OF DESIGNS F1 ANF F2
!     ------------------------------------------------------------------
!     IBEST CORRESPONDS TO MINIMUM VALUE DESIGN.
!     IF CONSTRAINTS ARE VIOLATED, IBEST CORRESPONDS TO MINIMUM
!     CONSTRAINT VIOLATION.
      IF (IGOOD1.EQ.0.AND.IGOOD2.EQ.0) GO TO 180
!     VIOLATED CONSTRAINTS.  PICK MINIMUM VIOLATION.
      IBEST=1
      IF (CV1.GE.CV2) IBEST=2
      GO TO 190
180   CONTINUE
!     NO CONSTRAINT VIOLATION.  PICK MINIMUM F.
      IBEST=1
      IF (F2.LE.F1) IBEST=2
190   CONTINUE
      II=1
      IF (NCON.EQ.0) GO TO 230
!     ------------------------------------------------------------------
!     *****                 2 - POINT INTERPOLATION                *****
!     ------------------------------------------------------------------
      III=0
200   III=III+1
      C1=G1(III)
      C2=G(III)
      IF (ISC(III).EQ.0) GO TO 210
!     ------------------------------------------------------------------
!                        LINEAR CONSTRAINT
!     ------------------------------------------------------------------
      IF (C1.GE.1.0E-5.AND.C1.LE.CTBM) GO TO 220
      CALL CNMN07 (II,ALP,ZRO,ZRO,C1,A2,C2,ZRO,ZRO)
      IF (ALP.LE.0.) GO TO 220
      IF (C1.GT.CTBM.AND.ALP.GT.ALPFES) ALPFES=ALP
      IF (C1.LT.CTL.AND.ALP.LT.ALPLN) ALPLN=ALP
      GO TO 220
210   CONTINUE
!     ------------------------------------------------------------------
!                     NON-LINEAR CONSTRAINT
!     ------------------------------------------------------------------
      IF (C1.GE.1.0E-5.AND.C1.LE.CTAM) GO TO 220
      CALL CNMN07 (II,ALP,ZRO,ZRO,C1,A2,C2,ZRO,ZRO)
      IF (ALP.LE.0.) GO TO 220
      IF (C1.GT.CTAM.AND.ALP.GT.ALPFES) ALPFES=ALP
      IF (C1.LT.CT.AND.ALP.LT.ALPNC) ALPNC=ALP
220   CONTINUE
      IF (III.LT.NCON) GO TO 200
230   CONTINUE
      IF (LINOBJ.GT.0.OR.SLOPE.GE.0.) GO TO 240
!     CALCULATE ALPHA TO MINIMIZE FUNCTION.
      CALL CNMN04 (II,ALPMIN,ZRO,ZRO,F1,SLOPE,A2,F2,ZRO,ZRO,ZRO,ZRO)
240   CONTINUE
!     ------------------------------------------------------------------
!                         PROPOSED MOVE
!     ------------------------------------------------------------------
!     MOVE AT LEAST FAR ENOUGH TO OVERCOME CONSTRAINT VIOLATIONS.
      A3=ALPFES
!     MOVE TO MINIMIZE FUNCTION.
      IF (ALPMIN.GT.A3) A3=ALPMIN
!     IF A3.LE.0, SET A3 = ALPSID.
      IF (A3.LE.0.) A3=ALPSID
!     LIMIT MOVE TO NEW CONSTRAINT ENCOUNTER.
      IF (A3.GT.ALPNC) A3=ALPNC
      IF (A3.GT.ALPLN) A3=ALPLN
!     MAKE A3 NON-ZERO.
      IF (A3.LE.1.0E-20) A3=1.0E-20
!     IF A3=A2=ALPSID AND F2 IS BEST, GO INVOKE SIDE CONSTRAINT
!     MODIFICATION.
      ALPB=1.-A2/A3
      ALPA=1.-ALPSID/A3
      JBEST=0
      IF (ABS(ALPB).LT.1.0E-10.AND.ABS(ALPA).LT.1.0E-10) JBEST=1
      IF (JBEST.EQ.1.AND.IBEST.EQ.2) GO TO 20
!     SIDE CONSTRAINT CHECK NOT SATISFIED.
      IF (NCON.EQ.0) GO TO 260
!     STORE CONSTRAINT VALUES IN G2.
      DO 250 I=1,NCON
      G2(I)=G(I)
250   CONTINUE
260   CONTINUE
!     IF A3=A2, SET A3=.9*A2.
      IF (ABS(ALPB).LT.1.0E-10) A3=.9*A2
!     MOVE AT LEAST .01*A2.
      IF (A3.LT.(.01*A2)) A3=.01*A2
!     LIMIT MOVE TO 5.*A2.
      IF (A3.GT.(5.*A2)) A3=5.*A2
!     LIMIT MOVE TO ALPSID.
      IF (A3.GT.ALPSID) A3=ALPSID
!     MOVE A DISTANCE A3*S.
      ALP=A3-A2
      ALPTOT=ALPTOT+ALP
      DO 270 I=1,NDV
      X(I)=X(I)+ALP*S(I)
270   CONTINUE
      IF (IPRINT.LT.5) GO TO 300
      WRITE (6,780)
      WRITE (6,740) A3
      IF (NSCAL.EQ.0) GO TO 290
      DO 280 I=1,NDV
280   G(I)=SCAL(I)*X(I)
      WRITE (6,750) (G(I),I=1,NDV)
      GO TO 300
290   WRITE (6,750) (X(I),I=1,NDV)
300   CONTINUE
!     ------------------------------------------------------------------
!              UPDATE FUNCTION AND CONSTRAINT VALUES
!     ------------------------------------------------------------------
      NCAL(1)=NCAL(1)+1
      JGOTO=2
      RETURN
310   CONTINUE
      F3=OBJ
      IF (IPRINT.GE.5) WRITE (6,760) F3
      IF (IPRINT.LT.5.OR.NCON.EQ.0) GO TO 320
      WRITE (6,770)
      WRITE (6,750) (G(I),I=1,NCON)
320   CONTINUE
!     ------------------------------------------------------------------
!       CALCULATE MAXIMUM CONSTRAINT VIOLATION AND PICK BEST DESIGN
!     ---------------------------------------------------------R--------
     CV3=0.
      IGOOD3=0
      NVC1=0                
      IF (NCON.EQ.0) GO TO 340
      DO 330 I=1,NCON
      CC=CTAM
      IF (ISC(I).GT.0) CC=CTBM
      C1=G(I)-CC
      IF (C1.GT.CV3) CV3=C1
      IF (C1.GT.0.) NVC1=NVC1+1
330   CONTINUE
      IF (CV3.GT.0.) IGOOD3=1
340   CONTINUE
!     DETERMINE BEST DESIGN.
      IF (IBEST.EQ.2) GO TO 360
!     CHOOSE BETWEEN F1 AND F3.
      IF (IGOOD1.EQ.0.AND.IGOOD3.EQ.0) GO TO 350
      IF (CV1.GE.CV3) IBEST=3
      GO TO 380
350   IF (F3.LE.F1) IBEST=3
      GO TO 380
360   CONTINUE
!     CHOOSE BETWEEN F2 AND F3.
      IF (IGOOD2.EQ.0.AND.IGOOD3.EQ.0) GO TO 370
      IF (CV2.GE.CV3) IBEST=3
      GO TO 380
370   IF (F3.LE.F2) IBEST=3
380   CONTINUE
      ALP=A3
      OBJ=F3
!     IF F3 VIOLATES FEWER CONSTRAINTS THAN F1 RETURN.
      IF (NVC1.LT.NVC) GO TO 710
!     IF OBJECTIVE AND ALL CONSTRAINTS ARE LINEAR, RETURN.
      IF (LINOBJ.NE.0.AND.NLNC.EQ.NCON) GO TO 710
!     IF A3 = ALPLN AND F3 IS BOTH GOOD AND BEST RETURN.
      ALPB=1.-ALPLN/A3
      IF ((ABS(ALPB).LT.1.0E-20.AND.IBEST.EQ.3).AND.(IGOOD3.EQ.0)) GO TO 710
!     IF A3 = ALPSID AND F3 IS BEST, GO INVOKE SIDE CONSTRAINT
!     MODIFICATION.
      ALPA=1.-ALPSID/A3
      IF (ABS(ALPA).LT.1.0E-20.AND.IBEST.EQ.3) GO TO 20
!     ------------------------------------------------------------------
!     **********            3 - POINT INTERPOLATION            *********
!     ------------------------------------------------------------------
      ALPNC=ALPSID
      ALPCA=ALPSID
      ALPFES=-1.
      ALPMIN=-1.
      IF (NCON.EQ.0) GO TO 440
      III=0
390   III=III+1
      C1=G1(III)
      C2=G2(III)
      C3=G(III)
      IF (ISC(III).EQ.0) GO TO 400
!     ------------------------------------------------------------------
!     LINEAR CONSTRAINT.  FIND ALPFES ONLY.  ALPLN SAME AS BEFORE.
!     ------------------------------------------------------------------
      IF (C1.LE.CTBM) GO TO 430
      II=1
      CALL CNMN07 (II,ALP,ZRO,ZRO,C1,A3,C3,ZRO,ZRO)
      IF (ALP.GT.ALPFES) ALPFES=ALP
      GO TO 430
400   CONTINUE
!     ------------------------------------------------------------------
!                     NON-LINEAR CONSTRAINT
!     ------------------------------------------------------------------
      II=2
      CALL CNMN07 (II,ALP,ZRO,ZRO,C1,A2,C2,A3,C3)
      IF (ALP.LE.ZRO) GO TO 430
      IF (C1.GE.CT.AND.C1.LE.0.) GO TO 410
      IF (C1.GT.CTAM.OR.C1.LT.0.) GO TO 420
!     ALP IS MINIMUM MOVE.  UPDATE FOR NEXT CONSTRAINT ENCOUNTER.
410   ALPA=ALP
      CALL CNMN07 (II,ALP,ALPA,ZRO,C1,A2,C2,A3,C3)
      IF (ALP.LT.ALPCA.AND.ALP.GE.ALPA) ALPCA=ALP
      GO TO 430
420   CONTINUE
      IF (ALP.GT.ALPFES.AND.C1.GT.CTAM) ALPFES=ALP
      IF (ALP.LT.ALPNC.AND.C1.LT.0.) ALPNC=ALP
430   CONTINUE
      IF (III.LT.NCON) GO TO 390
440   CONTINUE
      IF (LINOBJ.GT.0.OR.SLOPE.GT.0.) GO TO 450
!     ------------------------------------------------------------------
!              CALCULATE ALPHA TO MINIMIZE FUNCTION
!     ------------------------------------------------------------------
      II=3
      IF (A2.GT.A3.AND.(IGOOD2.EQ.0.AND.IBEST.EQ.2)) II=2
      CALL CNMN04 (II,ALPMIN,ZRO,ZRO,F1,SLOPE,A2,F2,A3,F3,ZRO,ZRO)
450   CONTINUE
!     ------------------------------------------------------------------
!                       PROPOSED MOVE
!     ------------------------------------------------------------------
!     MOVE AT LEAST ENOUGH TO OVERCOME CONSTRAINT VIOLATIONS.
      A4=ALPFES
!     MOVE TO MINIMIZE FUNCTION.
      IF (ALPMIN.GT.A4) A4=ALPMIN
!     IF A4.LE.0, SET A4 = ALPSID.
      IF (A4.LE.0.) A4=ALPSID
!     LIMIT MOVE TO NEW CONSTRAINT ENCOUNTER.
      IF (A4.GT.ALPLN) A4=ALPLN
      IF (A4.GT.ALPNC) A4=ALPNC
!     LIMIT MOVE TO RE-ENCOUNTER CURRENTLY ACTIVE CONSTRAINT.
      IF (A4.GT.ALPCA) A4=ALPCA
!     LIMIT A4 TO 5.*A3.
      IF (A4.GT.(5.*A3)) A4=5.*A3
!     UPDATE DESIGN.
      IF (IBEST.NE.3.OR.NCON.EQ.0) GO TO 470
!     STORE CONSTRAINT VALUES IN G2.  F3 IS BEST.  F2 IS NOT.
      DO 460 I=1,NCON
      G2(I)=G(I)
460   CONTINUE
470   CONTINUE
!     IF A4=A3 AND IGOOD1=0 AND IGOOD3=1, SET A4=.9*A3.
      ALP=A4-A3
      IF ((IGOOD1.EQ.0.AND.IGOOD3.EQ.1).AND.ABS(ALP).LT.1.0E-20) A4=.9*A3
!     ------------------------------------------------------------------
!                   MOVE A DISTANCE A4*S
!     ------------------------------------------------------------------
      ALP=A4-A3
      ALPTOT=ALPTOT+ALP
      DO 480 I=1,NDV
      X(I)=X(I)+ALP*S(I)
480   CONTINUE
      IF (IPRINT.LT.5) GO TO 510
      WRITE (6,720)
      WRITE (6,740) A4
      IF (NSCAL.EQ.0) GO TO 500
      DO 490 I=1,NDV
490   G(I)=SCAL(I)*X(I)
      WRITE (6,750) (G(I),I=1,NDV)
      GO TO 510
500   WRITE (6,750) (X(I),I=1,NDV)
510   CONTINUE
!     ------------------------------------------------------------------
!              UPDATE FUNCTION AND CONSTRAINT VALUES
!     ------------------------------------------------------------------
      NCAL(1)=NCAL(1)+1
      JGOTO=3
      RETURN
520   CONTINUE
      F4=OBJ
      IF (IPRINT.GE.5) WRITE (6,760) F4
      IF (IPRINT.LT.5.OR.NCON.EQ.0) GO TO 530
      WRITE (6,770)
      WRITE (6,750) (G(I),I=1,NCON)
530   CONTINUE
!     DETERMINE ACCAPTABILITY OF F4.
      IGOOD4=0
      CV4=0.
      IF (NCON.EQ.0) GO TO 550
      DO 540 I=1,NCON
      CC=CTAM
      IF (ISC(I).GT.0) CC=CTBM
      C1=G(I)-CC
      IF (C1.GT.CV4) CV4=C1
540   CONTINUE
      IF (CV4.GT.0.) IGOOD4=1
550   CONTINUE
      ALP=A4
      OBJ=F4
!     ------------------------------------------------------------------
!                     DETERMINE BEST DESIGN
!     ------------------------------------------------------------------
      GO TO (560,610,660), IBEST
560   CONTINUE
!     CHOOSE BETWEEN F1 AND F4.
      IF (IGOOD1.EQ.0.AND.IGOOD4.EQ.0) GO TO 570
      IF (CV1.GT.CV4) GO TO 710
      GO TO 580
570   CONTINUE
      IF (F4.LE.F1) GO TO 710
580   CONTINUE
!     F1 IS BEST.
      ALPTOT=ALPTOT-A4
      OBJ=F1
      DO 590 I=1,NDV
      X(I)=X(I)-A4*S(I)
590   CONTINUE
      IF (NCON.EQ.0) GO TO 710
      DO 600 I=1,NCON
      G(I)=G1(I)
600   CONTINUE
      GO TO 710
610   CONTINUE
!     CHOOSE BETWEEN F2 AND F4.
      IF (IGOOD2.EQ.0.AND.IGOOD4.EQ.0) GO TO 620
      IF (CV2.GT.CV4) GO TO 710
      GO TO 630
620   CONTINUE
      IF (F4.LE.F2) GO TO 710
630   CONTINUE
!     F2 IS BEST.
      OBJ=F2
      A2=A4-A2
      ALPTOT=ALPTOT-A2
      DO 640 I=1,NDV
      X(I)=X(I)-A2*S(I)
640   CONTINUE
      IF (NCON.EQ.0) GO TO 710
      DO 650 I=1,NCON
      G(I)=G2(I)
650   CONTINUE
      GO TO 710
660   CONTINUE
!     CHOOSE BETWEEN F3 AND F4.
      IF (IGOOD3.EQ.0.AND.IGOOD4.EQ.0) GO TO 670
      IF (CV3.GT.CV4) GO TO 710
      GO TO 680
670   CONTINUE
      IF (F4.LE.F3) GO TO 710
680   CONTINUE
!     F3 IS BEST.
      OBJ=F3
      A3=A4-A3
      ALPTOT=ALPTOT-A3
      DO 690 I=1,NDV
      X(I)=X(I)-A3*S(I)
690   CONTINUE
      IF (NCON.EQ.0) GO TO 710
      DO 700 I=1,NCON
      G(I)=G2(I)
700   CONTINUE
710   CONTINUE
      ALP=ALPTOT
      IF (IPRINT.GE.5) WRITE (6,790)
      JGOTO=0
      RETURN
!     ------------------------------------------------------------------
!                                  FORMATS
!     ------------------------------------------------------------------

720 FORMAT (/t6, 'THREE-POINT INTERPOLATION')
730 FORMAT (/// '* * * CONSTRAINED ONE-DIMENSIONAL SEARCH INFORMATION * * *')
740 FORMAT (//t6, 'PROPOSED DESIGN'/ t6, 'ALPHA =', e12.5/ t6, 'X-VECTOR')
750 FORMAT (' ', 8E12.4)
760 FORMAT (/t6, 'OBJ =', e13.5)
770 FORMAT (/t6, 'CONSTRAINT VALUES')
780 FORMAT (/t6, 'TWO-POINT INTERPOLATION')
790 FORMAT (/t6, '* * * END OF ONE-DIMENSIONAL SEARCH')

END SUBROUTINE CNMN06
      

!     ROUTINE TO FIND FIRST XBAR.GE.EPS CORRESPONDING TO A REAL ZERO
!     OF A ONE-DIMENSIONAL FUNCTION BY POLYNOMIEL INTERPOLATION.
!     BY G. N. VANDERPLAATS                          APRIL, 1972.
!     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
!     II = CALCULATION CONTROL.
!          1:  2-POINT LINEAR INTERPOLATION, GIVEN X1, Y1, X2 AND Y2.
!          2:  3-POINT QUADRATIC INTERPOLATION, GIVEN X1, Y1, X2, Y2,
!              X3 AND Y3.
!     EPS MAY BE NEGATIVE.
!     IF REQUIRED ZERO ON Y DOES NOT EXITS, OR THE FUNCTION IS
!     ILL-CONDITIONED, XBAR = EPS-1.0 WILL BE RETURNED AS AN ERROR
!     INDICATOR.
!     IF DESIRED INTERPOLATION IS ILL-CONDITIONED, A LOWER ORDER
!     INTERPOLATION, CONSISTANT WITH INPUT DATA, WILL BE ATTEMPTED AND
!     II WILL BE CHANGED ACCORDINGLY.
SUBROUTINE CNMN07 (II,XBAR,EPS,X1,Y1,X2,Y2,X3,Y3)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: II
DOUBLE PRECISION, INTENT(INOUT) :: XBAR, EPS, X1, Y1, X2, Y2, X3, Y3

!LOCAL VARIABLES
INTEGER :: JJ
DOUBLE PRECISION :: XBAR1, X21, YY, DY, X31, X32, QQ, AA, BB, CC, BAC, XB2
    
      XBAR1=EPS-1.
      XBAR=XBAR1
      JJ=0
      X21=X2-X1
      IF (ABS(X21).LT.1.0E-20) RETURN
      IF (II.EQ.2) GO TO 30

10    CONTINUE
!     ------------------------------------------------------------------
!                  II=1: 2-POINT LINEAR INTERPOLATION
!     ------------------------------------------------------------------
      II=1
      YY=Y1*Y2
      IF (JJ.EQ.0.OR.YY.LT.0.) GO TO 20
!     INTERPOLATE BETWEEN X2 AND X3.
      DY=Y3-Y2
      IF (ABS(DY).LT.1.0E-20) GO TO 20
      XBAR=X2+Y2*(X2-X3)/DY
      IF (XBAR.LT.EPS) XBAR=XBAR1
      RETURN
20    DY=Y2-Y1
!     INTERPOLATE BETWEEN X1 AND X2.
      IF (ABS(DY).LT.1.0E-20) RETURN
      XBAR=X1+Y1*(X1-X2)/DY
      IF (XBAR.LT.EPS) XBAR=XBAR1
      RETURN
30    CONTINUE
!     ------------------------------------------------------------------
!                 II=2: 3-POINT QUADRATIC INTERPOLATION
!     ------------------------------------------------------------------
      JJ=1
      X31=X3-X1
      X32=X3-X2
      QQ=X21*X31*X32
      IF (ABS(QQ).LT.1.0E-20) RETURN
      AA=(Y1*X32-Y2*X31+Y3*X21)/QQ
      IF (ABS(AA).LT.1.0E-20) GO TO 10
      BB=(Y2-Y1)/X21-AA*(X1+X2)
      CC=Y1-X1*(AA*X1+BB)
      BAC=BB*BB-4.*AA*CC
      IF (BAC.LT.0.) GO TO 10
      BAC=SQRT(BAC)
      AA=.5/AA
      XBAR=AA*(BAC-BB)
      XB2=-AA*(BAC+BB)
      IF (XBAR.LT.EPS) XBAR=XB2
      IF(XB2.LT.XBAR.AND.XB2.GT.EPS) XBAR=XB2
      IF (XBAR.LT.EPS) XBAR=XBAR1
      RETURN
END SUBROUTINE CNMN07



!     ROUTINE TO SOLVE SPECIAL LINEAR PROBLEM FOR IMPOSING S-TRANSPOSE
!     TIMES S.LE.1 BOUNDS IN THE MODIFIED METHOD OF FEASIBLE DIRECTIONS.
!     BY G. N. VANDERPLAATS                             APRIL, 1972.
!     NASA-AMES RESEARCH CENTER,  MOFFETT FIELD, CALIF.
!     REF.  'STRUCTURAL OPTIMIZATION BY METHODS OF FEASIBLE DIRECTIONS',
!     G. N. VANDERPLAATS AND F. MOSES, JOURNAL OF COMPUTERS
!     AND STRUCTURES, VOL 3, PP 739-755, 1973.
!     FORM OF L. P. IS BX=C WHERE 1ST NDB COMPONENTS OF X CONTAIN VECTOR
!     U AND LAST NDB COMPONENTS CONTAIN VECTOR V.  CONSTRAINTS ARE
!     U.GE.0, V.GE.0, AND U-TRANSPOSE TIMES V = 0.
!     NER = ERROR FLAG.  IF NER.NE.0 ON RETURN, PROCESS HAS NOT
!     CONVERGED IN 5*NDB ITERATIONS.
!     VECTOR MS1 IDENTIFIES THE SET OF BASIC VARIABLES.
SUBROUTINE CNMN08 (NDB,NER,C,MS1,B,N3,N4,N5)
IMPLICIT NONE
INTEGER, INTENT(INOUT) :: MS1(N5), N3, N4, N5, NDB, NER
DOUBLE PRECISION, INTENT(INOUT) :: C(N4), B(N3, N3)

!LOCAL VARIABLES
INTEGER :: I, J, M2, ITER1, NMAX, JJ, ICHK, KK
DOUBLE PRECISION :: EPS, CBMIN, CBMAX, C1, BI,  CB, BB, BB1
!IMPLICIT DOUBLE (A-H,O-Z)
!DIMENSION C(N4),B(N3,N3),MS1(N5)
!     ------------------------------------------------------------------
!     CHOOSE INITIAL BASIC VARIABLES AS V, AND INITIALIZE VECTOR MS1
!     ------------------------------------------------------------------
      NER=1
      M2=2*NDB
!     CALCULATE CBMIN AND EPS AND INITIALIZE MS1.
      EPS=-1.D10
      CBMIN=0.
      DO 10 I=1,NDB
      BI=B(I,I)
      CBMAX=0.
      IF (BI.LT.-1.0E-6) CBMAX=C(I)/BI
      IF (BI.GT.EPS) EPS=BI
      IF (CBMAX.GT.CBMIN) CBMIN=CBMAX
10    MS1(I)=0
      EPS=.0001*EPS
      IF (EPS.LT.-1.0E-10) EPS=-1.0E-10
      IF (EPS.GT.-.0001) EPS=-.0001
      CBMIN=CBMIN*1.0E-6
      IF (CBMIN.LT.1.0E-10) CBMIN=1.0E-10
      ITER1=0
      NMAX=5*NDB
!     ------------------------------------------------------------------
!     **********             BEGIN NEW ITERATION              **********
!     ------------------------------------------------------------------
20    ITER1=ITER1+1
      IF (ITER1.GT.NMAX) RETURN
!     FIND MAX. C(I)/B(I,I) FOR I=1,NDB.
      CBMAX=.9*CBMIN
      ICHK=0
      DO 30 I=1,NDB
      C1=C(I)
      BI=B(I,I)
      IF (BI.GT.EPS.OR.C1.GT.0.) GO TO 30
      CB=C1/BI
      IF (CB.LE.CBMAX) GO TO 30
      ICHK=I
      CBMAX=CB
30    CONTINUE
      IF (CBMAX.LT.CBMIN) GO TO 70
      IF (ICHK.EQ.0) GO TO 70
!     UPDATE VECTOR MS1.
      JJ=ICHK
      IF (MS1(JJ).EQ.0) JJ=ICHK+NDB
      KK=JJ+NDB
      IF (KK.GT.M2) KK=JJ-NDB
      MS1(KK)=ICHK
      MS1(JJ)=0
!     ------------------------------------------------------------------
!                     PIVOT OF B(ICHK,ICHK)
!     ------------------------------------------------------------------
      BB=1./B(ICHK,ICHK)
      DO 40 J=1,NDB
40    B(ICHK,J)=BB*B(ICHK,J)
      C(ICHK)=CBMAX
      B(ICHK,ICHK)=BB
!     ELIMINATE COEFICIENTS ON VARIABLE ENTERING BASIS AND STORE
!     COEFICIENTS ON VARIABLE LEAVING BASIS IN THEIR PLACE.
      DO 60 I=1,NDB
      IF (I.EQ.ICHK) GO TO 60
      BB1=B(I,ICHK)
      B(I,ICHK)=0.
      DO 50 J=1,NDB
50    B(I,J)=B(I,J)-BB1*B(ICHK,J)
      C(I)=C(I)-BB1*CBMAX
60    CONTINUE
      GO TO 20
70    CONTINUE
      NER=0
!     ------------------------------------------------------------------
!     STORE ONLY COMPONENTS OF U-VECTOR IN 'C'.  USE B(I,1) FOR
!     TEMPORARY STORAGE
!     ------------------------------------------------------------------
      DO 80 I=1,NDB
      B(I,1)=C(I)
80    CONTINUE
      DO 90 I=1,NDB
      C(I)=0.
      J=MS1(I)
      IF (J.GT.0) C(I)=B(J,1)
      IF (C(I).LT.0.) C(I)=0.
90    CONTINUE
      RETURN
END SUBROUTINE CNMN08
      
      
END MODULE CONSTRAINED_MINIMIZATION
