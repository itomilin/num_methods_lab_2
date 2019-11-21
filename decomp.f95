      SUBROUTINE DECOMP(NDIM,N,A,COND,IPVT,WORK)


      INTEGER NDIM,N

      REAL A(NDIM,N),COND,WORK(N)

      INTEGER IPVT(N)

      REAL EK,T,ANORM,YNORM,ZNORM

      INTEGER NM1,I,J,K,KP1,KB,KM1,M

      IPVT(N)=1

      IF(N.EQ.1)GO TO 80

      NM1=N-1

      ANORM=0.0

      DO 10 J=1,N

        T=0.0

        DO 5 I=1,N

          T=T+ABS(A(I,J))

    5   CONTINUE

        IF(T.GT.ANORM) ANORM=T

   10 CONTINUE


      DO 35 K=1,NM1

        KP1=K+1

        M=K

        DO 15 I=KP1,N

          IF(ABS(A(I,K)).GT.ABS(A(M,K))) M=I

   15   CONTINUE

        IPVT(K)=M

        IF(M.NE.K)IPVT(N)=-IPVT(N)

        T=A(M,K)

        A(M,K)=A(K,K)

        A(K,K)=T

        IF(T.EQ.0.0)GO TO 35

        DO 20 I=KP1,N

          A(I,K)=-A(I,K)/T

   20   CONTINUE

        DO 30 J=KP1,N

          T=A(M,J)

          A(M,J)=A(K,J)

          A(K,J)=T

          IF(T.EQ.0.0)GO TO 30

          DO 25 I=KP1,N

            A(I,J)=A(I,J)+A(I,K)*T

   25     CONTINUE

   30   CONTINUE

   35 CONTINUE

      DO 50 K=1,N

        T=0.0

        IF(K.EQ.1)GO TO 45

        KM1=K-1

        DO 40 I=1,KM1

          T=T+A(I,K)*WORK(I)

   40   CONTINUE

   45   EK=1.0

        IF(T.LT.0.0)EK=-1.0

        IF(A(K,K).EQ.0.0)GO TO 90

        WORK(K)=-(EK+T)/A(K,K)

   50 CONTINUE

      DO 60 KB=1,NM1

        K=N-KB

        T=WORK(K)

        KP1=K+1

        DO 55 I=KP1,N

          T=T+A(I,K)*WORK(I)

   55   CONTINUE

        WORK(K)=T

        M=IPVT(K)

        IF(M.EQ.K)GO TO 60

        T=WORK(M)

        WORK(M)=WORK(K)

        WORK(K)=T

   60 CONTINUE

      YNORM=0.0

      DO 65 I=1,N

        YNORM=YNORM+ABS(WORK(I))

   65 CONTINUE

      CALL SOLVE(NDIM,N,A,WORK,IPVT)

      ZNORM=0.0

      DO 70 I=1,N

        ZNORM=ZNORM+ABS(WORK(I))

   70 CONTINUE

      COND=ANORM*ZNORM/YNORM

      IF(COND.LT.1.0)COND=1.0

      RETURN

   80 COND=1.0

      IF(A(1,1).NE.0.0)RETURN

   90 CONTINUE

      COND=1.0E+32

      RETURN

      END
