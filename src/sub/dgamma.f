#ifdef ALLOW_TAPENADE
C*ID* DGAMMA   C109.PTOLEMY.VAXSTUFF                        PER704  21:
      FUNCTION DGAMMA ( X )
C
C     Complete Gamma Function and related functions.
C
C
C     Entry Points (all are double precision):
C
C     DGAMMA:  Complete Gamma Function.
C     DGAMAP:  First derivative of the Complete Gamma Function.
C     DLGAMA, DLNGAM:  Natural logarithm of DABS( GAMMA(X) ).
C     DIGAMA:  Digamma, or Psi, function.
C
C
C          All of these functions must be declared DOUBLE PRECISION
C     (or REAL*8) in the calling program, and all require a DOUBLE
C     PRECISION argument.
C
C
C          Reference:  "A&S" below refers to Abramowitz and Stegun,
C     Handbook of Mathematical Functions, published by the National
C     Bureau of Standards (1964) and by Dover (1965).
C
C     METHODS:
C          DGAMMA and DLGAMA (or DLNGAM) are found using the same
C     methods (below).  Logs or exponentials are taken at the end.
C          For X < o.5 (including negative X), the function at 1-X
C     is computed (as below) and transformed back using A&S 6.1.17.
C     Hence we assume below that X >= 0.5.
C          For X < 12, Gamma(X-N) is computed (as below), where N is
C     the integer such that 0.5 <= X-N < 1.5.  The result is transformed
C     back up (if necesary) using A&S 6.1.15.
C          For 0.5 <= X < 1.5, we compute 1/Gamma(X) using the first
C     19 terms of the series expansion in A&S 6.1.34.
C          For 12 <= X < 164, we compute Log(Gamma) using 6 terms of
C     the asymptotic expansion in A&S 6.1.40.  For X >= 164, only 2
C     terms are necessary.
C
C          DIGAMA(X) is found using A&S 6.3.18, using 8 terms, for
C     ºXº >= 9.0.  For ºXº < 9.0, the argument is transformed using
C     A&S 6.3.6 to bring it above 10.0.  For X < 0.5, DIGAMA(1-X)
C     is found and transformed back using A&S 6.3.7.
C
C          ACCURACY:  In the following, "error" refers to absolute
C     or relative error, whichever is smaller.
C          DLGAMA/DLNGAM and DIGAMA:  For X > 0.5, errors in these
C     functions are less than 1.E-16 on a VAX.  Absolute errors for
C     X < 0.5 should be comparable to those for 1-X; hence, they
C     will be much worse when the function is much smaller at
C     X (<0.5) than at 1-X.
C          DGAMMA:  For -11 < X < 12, errors are less than 1.E-16.
C     For X >= 12 and X <= 11, accuracy is limited by the DEXP function,
C     because to first order in D we have
C            EXP[(1+D)*Y] = EXP(Y) * (1 + Y*D).
C     Hence, the relative error in DGAMMA(X) is expected to be DLGAMA(X)
C     times the relative error in DLGAMA(X).  This expectation is
C     verified by the tests.  Relative errors in DGAMMA are about 3.E-16
C     at X = 12, and climb to 8.E-16 at X = 19.
C          DGAMAP:  Errors are comparable to DGAMMA.
C
C          ERRORS:  If the argument is too big or too near a negative
C     integer or zero, an appropriate error message is printed, and then
C     execution is allowed to continue.  In most cases, the program will
C     subsequently bomb with a divide by zero or an overflow in the DEXP
C     function.  The maximum argument for DGAMMA or DGAMAP is about 34
C     for a VAX, 57 for IBM, and 300 for CDC.
C
C     1/26/80 - Robert P. Goddard
C
C
      IMPLICIT REAL*8  ( A-H, O-Z )
C
      LOGICAL  LNEG, LEXP
      CHARACTER*6 NAME(4)
#if 0
      DATA  NAME / 'DGAMMA', 'DGAMAP', 'DLGAMA', 'DIGAMA' /
      DATA  PI / 3.1415 92653 58979 324 0D0 /
      DATA  R2PILN / 0.91893 85332 04672 742 0D0 /
#else
      NAME(1) = 'DGAMMA'
      NAME(2) = 'DGAMAP'
      NAME(3) = 'DLGAMA'
      NAME(4) = 'DIGAMA'
      PI = 3.1415 92653 58979 324 0D0
      R2PILN = 0.91893 85332 04672 742 0D0
#endif
C
C     ELIM is the biggest number s.t. DEXP(ELIM) doesn't overflow.
C
c     very rought guess for alpha g-format:
c
#if 0
      DATA  ELIM / 695. 0D0 /
#else
      ELIM = 695.0D0
#endif
CVA      DATA  ELIM / 88.029 0D0 /
CIB      DATA  ELIM / 176.058 0D0 /
CUN      DATA  ELIM / 352.116 0D0 /
CDC      DATA  ELIM / 704.23 /
CRA      DATA  ELIM / 5633.8 /
C
C
C - - - - - - - - - - DGAMMA - - - - - - - - - -
C
      NE = 1
      GO TO 100
C
 10   DGAMMA = G
      RETURN
C
C
C     Log Gamma code common to DGAMMA, DLGAMA, and DLNGAM.
C
C
C     Check for argument < 0.5 and transform.
C
 100  Z = X
      LEXP = .FALSE.
      LNEG = Z .LT. 0.50D0
      IF ( LNEG )  Z = 1.D0 - Z
C
C     If Z < 12, we transform it down into the range 0.5 <= Z < 1.5
C     (A&S 6.1.15) and use a 19-term series expansion of 1/Gamma(1+Z)
C     (from A&S 6.1.34).
C
      IF ( Z .GE. 12.D0 )  GO TO 200
      H = 1.D0
      N = Z - 0.50D0
      IF ( N .EQ. 0 )  GO TO 120
      DO 110  I = 1, N
         Z = Z - 1.D0
         H = H * Z
 110  CONTINUE
C
 120  Y = Z - 1.D0
      G = -.00002 01348 547807 0D0
      IF ( DABS( Y ) .GT. 0.10D0 )  G = G
     D   + Y * ( -.00000 12504 934821 0D0
     E   + Y * ( 0.00000 11330 272320 0D0
     F   + Y * ( -.00000 02056 338417 0D0
     G   + Y * ( 0.00000 00061 160950 0D0
     H   + Y * ( 0.00000 00050 020075 0D0
     I   + Y * ( -.00000 00011 812746 0D0
     J   + Y * ( 0.00000 00001 043427 0D0
     X   ) ) ) ) ) ) )
      G = 1.0D0
     2   + Y * ( 0.57721 56649 015329 0D0
     3   + Y * ( -.65587 80715 202538 0D0
     4   + Y * ( -.04200 26350 340952 0D0
     5   + Y * ( 0.16653 86113 822915 0D0
     6   + Y * ( -.04219 77345 555443 0D0
     7   + Y * ( -.00962 19715 278770 0D0
     8   + Y * ( 0.00721 89432 466630 0D0
     9   + Y * ( -.00116 51675 918591 0D0
     A   + Y * ( -.00021 52416 741149 0D0
     B   + Y * ( 0.00012 80502 823882 0D0
     C   + Y * G ) ) ) ) ) ) ) ) ) )
      G = H / G
      LEXP = .TRUE.
      GO TO 250
C
C     Calculate log gamma using 6 terms of the Bernoulli expansion.
C     For ºXº > 164, we can get 17 decimal places with 2 terms.
C     The constant expressions for the Bernoulli numbers [divided
C     by (2*m)*(2*m-1)] are evaluated at compile time.  The Bernoulli
C     numbers come from A&S Table 23.2.
C
 200  ZI = 1.D0 / Z
      ZIS = ZI * ZI
      G = (-1.D0/360.D0)
      IF ( Z .LT. 164.D0 )  G = G
     3        + ZIS * ( (1.D0/1260.D0)
     4          + ZIS * ( (-1.D0/1680.D0)
     5            + ZIS * ( (1.D0/1188.D0)
     6              + ZIS * ( (-691.D0/312840.D0)
     X                ) ) ) )
      G = R2PILN - Z + ( Z - 0.50D0 ) * DLOG( Z )
     1    + ZI * ( (1.D0/12.D0)
     2      + ZIS * G )
C
C      If X < 0.5, transform the result using A&S 6.1.17.
C
 250  IF ( .NOT. LNEG )  GO TO 300
         T = DSIN( PI*Z )
         IF ( T .NE. 0.D0 )  GO TO 260
            write (6, 803)  NAME(NE), X
            write (6, 813)
 260     IF ( LEXP )  GO TO 270
            G = DLOG( PI / DABS( T ) ) - G
            LNEG = T .LT. 0
            GO TO 300
 270     CONTINUE
            G = PI / ( G * T )
            LNEG = .FALSE.
C
C     We now have log gamma (or gamma if LEXP = .TRUE.).
C     Branch depends on the entry point used.
C
 300  IF ( LEXP )  GO TO 10
      IF ( G .LE. ELIM )  GO TO 350
         write (6, 803)  NAME(NE), X
         write (6, 913)
 350  G = DEXP( G )
      IF ( LNEG )  G = -G
      IF ( NE .EQ. 1 )  GO TO 10
C
C     Error messages.
C
 803  FORMAT ( '0**** ERROR IN ', A6, ':  ARGUMENT =', G25.16, ' ****' )
 813  FORMAT ( ' **** ARGUMENT IS A NEGATIVE INTEGER OR ZERO ****' )
 913  FORMAT ( ' **** EXPONENT OVERFLOW ****' )
      END
#else
C*ID* DGAMMA   C109.PTOLEMY.VAXSTUFF                        PER704  21:
      FUNCTION DGAMMA ( X )
C
C     Complete Gamma Function and related functions.
C
C
C     Entry Points (all are double precision):
C
C     DGAMMA:  Complete Gamma Function.
C     DGAMAP:  First derivative of the Complete Gamma Function.
C     DLGAMA, DLNGAM:  Natural logarithm of DABS( GAMMA(X) ).
C     DIGAMA:  Digamma, or Psi, function.
C
C
C          All of these functions must be declared DOUBLE PRECISION
C     (or REAL*8) in the calling program, and all require a DOUBLE
C     PRECISION argument.
C
C
C          Reference:  "A&S" below refers to Abramowitz and Stegun,
C     Handbook of Mathematical Functions, published by the National
C     Bureau of Standards (1964) and by Dover (1965).
C
C     METHODS:
C          DGAMMA and DLGAMA (or DLNGAM) are found using the same
C     methods (below).  Logs or exponentials are taken at the end.
C          For X < o.5 (including negative X), the function at 1-X
C     is computed (as below) and transformed back using A&S 6.1.17.
C     Hence we assume below that X >= 0.5.
C          For X < 12, Gamma(X-N) is computed (as below), where N is
C     the integer such that 0.5 <= X-N < 1.5.  The result is transformed
C     back up (if necesary) using A&S 6.1.15.
C          For 0.5 <= X < 1.5, we compute 1/Gamma(X) using the first
C     19 terms of the series expansion in A&S 6.1.34.
C          For 12 <= X < 164, we compute Log(Gamma) using 6 terms of
C     the asymptotic expansion in A&S 6.1.40.  For X >= 164, only 2
C     terms are necessary.
C
C          DIGAMA(X) is found using A&S 6.3.18, using 8 terms, for
C     ºXº >= 9.0.  For ºXº < 9.0, the argument is transformed using
C     A&S 6.3.6 to bring it above 10.0.  For X < 0.5, DIGAMA(1-X)
C     is found and transformed back using A&S 6.3.7.
C
C          ACCURACY:  In the following, "error" refers to absolute
C     or relative error, whichever is smaller.
C          DLGAMA/DLNGAM and DIGAMA:  For X > 0.5, errors in these
C     functions are less than 1.E-16 on a VAX.  Absolute errors for
C     X < 0.5 should be comparable to those for 1-X; hence, they
C     will be much worse when the function is much smaller at
C     X (<0.5) than at 1-X.
C          DGAMMA:  For -11 < X < 12, errors are less than 1.E-16.
C     For X >= 12 and X <= 11, accuracy is limited by the DEXP function,
C     because to first order in D we have
C            EXP[(1+D)*Y] = EXP(Y) * (1 + Y*D).
C     Hence, the relative error in DGAMMA(X) is expected to be DLGAMA(X)
C     times the relative error in DLGAMA(X).  This expectation is
C     verified by the tests.  Relative errors in DGAMMA are about 3.E-16
C     at X = 12, and climb to 8.E-16 at X = 19.
C          DGAMAP:  Errors are comparable to DGAMMA.
C
C          ERRORS:  If the argument is too big or too near a negative
C     integer or zero, an appropriate error message is printed, and then
C     execution is allowed to continue.  In most cases, the program will
C     subsequently bomb with a divide by zero or an overflow in the DEXP
C     function.  The maximum argument for DGAMMA or DGAMAP is about 34
C     for a VAX, 57 for IBM, and 300 for CDC.
C
C     1/26/80 - Robert P. Goddard
C
C
      IMPLICIT REAL*8  ( A-H, O-Z )
C
      LOGICAL  LNEG, LEXP
      CHARACTER*6 NAME(4)
      DATA  NAME / 'DGAMMA', 'DGAMAP', 'DLGAMA', 'DIGAMA' /
      DATA  PI / 3.1415 92653 58979 324 0D0 /
      DATA  R2PILN / 0.91893 85332 04672 742 0D0 /
C
C     ELIM is the biggest number s.t. DEXP(ELIM) doesn't overflow.
C
c     very rought guess for alpha g-format:
c
      DATA  ELIM / 695. 0D0 /
CVA      DATA  ELIM / 88.029 0D0 /
CIB      DATA  ELIM / 176.058 0D0 /
CUN      DATA  ELIM / 352.116 0D0 /
CDC      DATA  ELIM / 704.23 /
CRA      DATA  ELIM / 5633.8 /
C
C
C - - - - - - - - - - DGAMMA - - - - - - - - - -
C
      NE = 1
      GO TO 100
C
 10   DGAMMA = G
      RETURN
C
C - - - - - - - - - - DGAMAP - - - - - - - - - -
C
      ENTRY DGAMAP (X)
C
      NE = 2
      GO TO 100
C
 20   DGAMAP = G * D
      RETURN
C
C - - - - - - - - DLGAMA, DLNGAM - - - - - - - -
C
      ENTRY DLGAMA ( X )
      ENTRY DLNGAM ( X )
C
      NE = 3
      GO TO 100
C
 30   IF ( LEXP )  G = DLOG( G )
      DLGAMA = G
      RETURN
C
C - - - - - - - - - - DIGAMA - - - - - - - - - -
C
      ENTRY DIGAMA ( X )
C
      NE = 4
      GO TO 500
C
 40   DIGAMA = D
      RETURN
C
C
C     Log Gamma code common to DGAMMA, DLGAMA, and DLNGAM.
C
C
C     Check for argument < 0.5 and transform.
C
 100  Z = X
      LEXP = .FALSE.
      LNEG = Z .LT. 0.50D0
      IF ( LNEG )  Z = 1.D0 - Z
C
C     If Z < 12, we transform it down into the range 0.5 <= Z < 1.5
C     (A&S 6.1.15) and use a 19-term series expansion of 1/Gamma(1+Z)
C     (from A&S 6.1.34).
C
      IF ( Z .GE. 12.D0 )  GO TO 200
      H = 1.D0
      N = Z - 0.50D0
      IF ( N .EQ. 0 )  GO TO 120
      DO 110  I = 1, N
         Z = Z - 1.D0
         H = H * Z
 110  CONTINUE
C
 120  Y = Z - 1.D0
      G = -.00002 01348 547807 0D0
      IF ( DABS( Y ) .GT. 0.10D0 )  G = G
     D   + Y * ( -.00000 12504 934821 0D0
     E   + Y * ( 0.00000 11330 272320 0D0
     F   + Y * ( -.00000 02056 338417 0D0
     G   + Y * ( 0.00000 00061 160950 0D0
     H   + Y * ( 0.00000 00050 020075 0D0
     I   + Y * ( -.00000 00011 812746 0D0
     J   + Y * ( 0.00000 00001 043427 0D0
     X   ) ) ) ) ) ) )
      G = 1.0D0
     2   + Y * ( 0.57721 56649 015329 0D0
     3   + Y * ( -.65587 80715 202538 0D0
     4   + Y * ( -.04200 26350 340952 0D0
     5   + Y * ( 0.16653 86113 822915 0D0
     6   + Y * ( -.04219 77345 555443 0D0
     7   + Y * ( -.00962 19715 278770 0D0
     8   + Y * ( 0.00721 89432 466630 0D0
     9   + Y * ( -.00116 51675 918591 0D0
     A   + Y * ( -.00021 52416 741149 0D0
     B   + Y * ( 0.00012 80502 823882 0D0
     C   + Y * G ) ) ) ) ) ) ) ) ) )
      G = H / G
      LEXP = .TRUE.
      GO TO 250
C
C     Calculate log gamma using 6 terms of the Bernoulli expansion.
C     For ºXº > 164, we can get 17 decimal places with 2 terms.
C     The constant expressions for the Bernoulli numbers [divided
C     by (2*m)*(2*m-1)] are evaluated at compile time.  The Bernoulli
C     numbers come from A&S Table 23.2.
C
 200  ZI = 1.D0 / Z
      ZIS = ZI * ZI
      G = (-1.D0/360.D0)
      IF ( Z .LT. 164.D0 )  G = G
     3        + ZIS * ( (1.D0/1260.D0)
     4          + ZIS * ( (-1.D0/1680.D0)
     5            + ZIS * ( (1.D0/1188.D0)
     6              + ZIS * ( (-691.D0/312840.D0)
     X                ) ) ) )
      G = R2PILN - Z + ( Z - 0.50D0 ) * DLOG( Z )
     1    + ZI * ( (1.D0/12.D0)
     2      + ZIS * G )
C
C      If X < 0.5, transform the result using A&S 6.1.17.
C
 250  IF ( .NOT. LNEG )  GO TO 300
         T = DSIN( PI*Z )
         IF ( T .NE. 0.D0 )  GO TO 260
            write (6, 803)  NAME(NE), X
            write (6, 813)
 260     IF ( LEXP )  GO TO 270
            G = DLOG( PI / DABS( T ) ) - G
            LNEG = T .LT. 0
            GO TO 300
 270     CONTINUE
            G = PI / ( G * T )
            LNEG = .FALSE.
C
C     We now have log gamma (or gamma if LEXP = .TRUE.).
C     Branch depends on the entry point used.
C
 300  IF ( NE .EQ. 3 )  GO TO 30
      IF ( LEXP )  GO TO ( 10, 500 ), NE
      IF ( G .LE. ELIM )  GO TO 350
         write (6, 803)  NAME(NE), X
         write (6, 913)
 350  G = DEXP( G )
      IF ( LNEG )  G = -G
      IF ( NE .EQ. 1 )  GO TO 10
C
C
C     Digamma code common to DGAMAP, DIGAMMA.
C
C     Check for argument < 0.5 and transform.
C
 500  Z = X
      LNEG = Z .LT. 0.50D0
      IF ( LNEG )  Z = 1.D0 - Z
C
C     Transform to Z >= 9 if necessary.
C
 550  H = 0.D0
      IF ( Z .GE. 9.D0 )  GO TO 600
      N = 10.D0 - Z
      DO 580  I = 1, N
         H = H + 1.D0/Z
         Z = Z + 1.D0
 580  CONTINUE
C
C     Calculate Digamma using 8 terms of the Bernoulli expansion.
C     The constant expressions for the Bernoulli numbers (divided
C     by 2*n) are evaluated at compile time.
C
 600  ZI = 1.D0/Z
      ZIS = ZI*ZI
      D = DLOG(Z) - 0.50D0*ZI - H
     1    - ZIS * ( (1.D0/12.D0)
     2      + ZIS * ( (-1.D0/120.D0)
     3        + ZIS * ( (1.D0/252.D0)
     4          + ZIS * ( (-1.D0/240.D0)
     5            + ZIS * ( (1.D0/132.D0)
     6              + ZIS * ( (-691.D0/32760.D0)
     7                + ZIS * ( (1.D0/12.D0)
     8                  + ZIS * ( (-3617.D0/8160.D0)
     A            ) ) ) ) ) ) ) )
C
C     Modify the result if the argument was negative.
C
      IF ( .NOT. LNEG )  GO TO 700
         T = DSIN( PI*Z )
         IF ( T .NE. 0.D0 )  GO TO 690
            write (6, 803)  NAME(NE), X
            write (6, 813)
 690     D = D + PI*DCOS( PI*Z )/T
 700  IF ( NE .EQ. 4 )  GO TO 40
      GO TO 20
C
C     Error messages.
C
 803  FORMAT ( '0**** ERROR IN ', A6, ':  ARGUMENT =', G25.16, ' ****' )
 813  FORMAT ( ' **** ARGUMENT IS A NEGATIVE INTEGER OR ZERO ****' )
 913  FORMAT ( ' **** EXPONENT OVERFLOW ****' )
      END
#endif
