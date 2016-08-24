!
!=================================================================================
!         G=-1/(4*PI)*(1/r+1/r1)
!         
!=================================================================================
!

        SUBROUTINE SOLIDANGLE(INODE,SANGLE)
	  USE MVAR_MOD
	  USE MFUNC_mod
	   
        IMPLICIT   NONE 
	  
	  INTEGER  INODE,IELEM,I,J,LK,LJ,LI,IEL
	  INTEGER  MCELE,MCNT(100),MEND(100),MNOD(100)
	  INTEGER  NCONT(100),NSIDE(100,2)
	  REAL(4)   SANGLE,DR,SI,ETA,DUM
	  REAL(4)   SIQ(12),ETAQ(12),SIT(9),ETAT(9),SF(8),DSF(2,8)
	  REAL(4)   DSIQ(12,2),DSIT(9,2)
	  REAL(4)   DXYZN(100,3),TXYZN(100,2,3),XJ(2,3)

!
!       XITSI, XITETA are the coordinates of nodes for the 
!       quadrilateral element                   
!
        DATA  SIQ/-1.0D0, 0.0D0, 0.0D0, 1.0D0, 1.0D0, 1.0D0,
	1		     1.0D0, 0.0D0, 0.0D0,-1.0D0,-1.0D0,-1.0D0/ 

        DATA ETAQ/-1.0D0,-1.0D0,-1.0D0,-1.0D0, 0.0D0, 0.0D0, 
	1		     1.0D0, 1.0D0, 1.0D0, 1.0D0, 0.0D0, 0.0D0/ 

	  DATA  DSIQ / 1.0d0,-1.0d0, 1.0d0,-1.0d0,-1.0d0,-1.0d0,
     1	          -1.0d0, 1.0d0,-1.0d0, 1.0d0, 1.0d0, 1.0d0,      
     2               1.0d0, 1.0d0, 1.0d0, 1.0d0,-1.0d0, 1.0d0,
     3	          -1.0d0,-1.0d0,-1.0d0,-1.0d0, 1.0d0,-1.0d0/ 
!
!
!       XITSI, XITETA are the coordinates of nodes for the 
!       triangle element                   
!

        DATA SIT /0.0D0,1.0D0,0.0D0,0.5D0,0.5D0,0.5D0,0.5D0,0.0D0,0.0D0/
        DATA ETAT/0.0D0,0.0D0,1.0D0,0.0D0,0.0D0,0.5D0,0.5D0,0.5D0,0.5D0/

	  DATA  DSIT /1.0d0,-1.0d0, 1.0d0, 1.0d0,-1.0d0, 
	1		      1.0d0, 1.0d0, 1.0d0, 1.0d0,
     2              1.0d0, 1.0d0,-1.0d0, 1.0d0, 1.0d0,
     3			  1.0d0,-1.0d0,-1.0d0, 1.0d0/ 

!
!
! ============================================================
! 

	  MCELE=0
	  DO IELEM=1, NELEM
          DO I=1, NCN(IELEM)
            IF(INODE .EQ. NCON(IELEM,I)) THEN
!
	        IF(MCELE .GT. 99) THEN
	          Print *,' Number of elements at the node is ',MCELE
	        ENDIF
!
	        IF(NCN(IELEM) .EQ. 8)  THEN
	          IF(I .EQ. 1)  THEN
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=1
	          ELSE IF(I .EQ. 2)  THEN
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=2
!
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=3
	          ELSE IF(I .EQ. 3)  THEN
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=4
	          ELSE IF(I .EQ. 4)  THEN
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=5
!
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=6
	          ELSE IF(I .EQ. 5)  THEN
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=7
	          ELSE IF(I .EQ. 6)  THEN
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=8
!
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=9
	          ELSE IF(I .EQ. 7)  THEN
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
		        MNOD(MCELE)=10
	          ELSE IF(I .EQ. 8)  THEN
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=11
		      		 
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
   	            MEND(MCELE)=I
	            MNOD(MCELE)=12 		 
		      END IF
!
	        ELSE IF(NCN(IELEM) .EQ. 6)  THEN
	          IF(I .EQ. 1)  THEN
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=1
	          ELSE IF(I .EQ. 2)  THEN
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=2
	          ELSE IF(I .EQ. 3)  THEN
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=3
	          ELSE IF(I .EQ. 4)  THEN
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=4

	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=5
	          ELSE IF(I .EQ. 5)  THEN
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=6

	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=7
	          ELSE IF(I .EQ. 6)  THEN
	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=8

	            MCELE=MCELE+1
	            MCNT(MCELE)=IELEM
	            MEND(MCELE)=I
	            MNOD(MCELE)=9
	          END IF
!
              END IF             ! NCN(IELEM) .EQ. 8 or 6
!
	      END IF

          END DO
        END DO
!
!
	  DXYZN=0.0d0
	  TXYZN=0.0d0


	  DO IEL=1, MCELE
!	    MCELE=MCELE+1
	    IELEM=MCNT(IEL)

		DXYZN(IEL,1)=DXYZE(1,MEND(IEL),IELEM)
		DXYZN(IEL,2)=DXYZE(2,MEND(IEL),IELEM)
		DXYZN(IEL,3)=DXYZE(3,MEND(IEL),IELEM)
!
          IF(NCN(IELEM).EQ.8) THEN
!
!           Quadrilateral element
!
            SI =SIQ(MNOD(IEL))
            ETA=ETAQ(MNOD(IEL))
        
            CALL SPFUNC8(SI,ETA,SF,DSF)
!
          ELSE IF(NCN(IELEM).EQ.6) THEN
!
!           Trianglar element
!
            SI =SIT(MNOD(IEL))
            ETA=ETAT(MNOD(IEL))
       
            CALL SPFUNC6(SI,ETA,SF,DSF)
	    END IF


!         evaluate the Jacobian matrix at (SI,ETA),  XJ(2,3)
!        
!         LI: 1--SI,  2--ETA
!         LJ: 1--X,   2--Y,   3--Z
!

          DO LI=1,2
            DO LJ=1,3
              DUM=0.0D0
              DO LK=1, NCN(IELEM)
                DUM=DUM+DSF(LI,LK)*XYZ(LJ,NCON(IELEM,LK))
              END DO
	        XJ(LI,LJ)=DUM
            END DO
!
	      
		  DR=SQRT(XJ(LI,1)*XJ(LI,1)+XJ(LI,2)*XJ(LI,2)
	1	                           +XJ(LI,3)*XJ(LI,3))

	      IF(NCN(IELEM).EQ.8) THEN
	        TXYZN(IEL,LI,1)=DSIQ(MNOD(IEL),LI)*XJ(LI,1)/DR
	        TXYZN(IEL,LI,2)=DSIQ(MNOD(IEL),LI)*XJ(LI,2)/DR
	        TXYZN(IEL,LI,3)=DSIQ(MNOD(IEL),LI)*XJ(LI,3)/DR
            ELSE IF(NCN(IELEM).EQ.6) THEN
	        TXYZN(IEL,LI,1)=DSIT(MNOD(IEL),LI)*XJ(LI,1)/DR
	        TXYZN(IEL,LI,2)=DSIT(MNOD(IEL),LI)*XJ(LI,2)/DR
	        TXYZN(IEL,LI,3)=DSIT(MNOD(IEL),LI)*XJ(LI,3)/DR
	      END IF
!
!	      print *, TXYZN(IEL,LI,1),TXYZN(IEL,LI,2),TXYZN(IEL,LI,3)
          END DO

        END DO

! 

	 
	  IF(ABS(XYZ(3,INODE)+H) .LT. 1.0D-6)  THEN

!         The Node is at the waterline, on the sea bed, or at the symmetric plans
!
	    CALL XYSYM(DXYZN,TXYZN,MCELE)

	  END IF 

!
!	  Print *,'NODQUA=',NODQUA(INODE)
	   
	  IF(NODQUA(INODE) .EQ. 0) THEN
	    CALL DOMAINFAS(DXYZN,TXYZN,MCELE,NCONT,NSIDE)
	    CALL ANGLE0(DXYZN,TXYZN,MCELE,NCONT,NSIDE,SANGLE)
	 
	  ELSE IF(NODQUA(INODE) .EQ. 2) THEN
	    CALL ZXSYM(DXYZN,TXYZN,MCELE)
	    CALL DOMAINFAS(DXYZN,TXYZN,MCELE,NCONT,NSIDE)
	    CALL ANGLE0(DXYZN,TXYZN,MCELE,NCONT,NSIDE,SANGLE) 
	  ELSE IF(NODQUA(INODE) .EQ. 4) THEN
	    CALL YZSYM(DXYZN,TXYZN,MCELE)
	    CALL DOMAINFAS(DXYZN,TXYZN,MCELE,NCONT,NSIDE)
	    CALL ANGLE0(DXYZN,TXYZN,MCELE,NCONT,NSIDE,SANGLE)
	 
	  ELSE IF(NODQUA(INODE) .EQ. 5) THEN
	    CALL YZZXSYM(DXYZN,TXYZN,MCELE)
	    CALL DOMAINFAS(DXYZN,TXYZN,MCELE,NCONT,NSIDE)
	    CALL ANGLE0(DXYZN,TXYZN,MCELE,NCONT,NSIDE,SANGLE)
	 
	  END IF

	  SANGLE=1.0d0-SANGLE/(4.0d0*PI)

!	  Print *,'     Solid_ANGLE=',SANGLE
!
!	  pause
!
        RETURN
        END SUBROUTINE SOLIDANGLE         

!
!
!    
        SUBROUTINE  XYSYM(DXYZN,TXYZN,MCELE)
        IMPLICIT   NONE 
	  
	  INTEGER  I,J
	  INTEGER  MCELE

	  REAL(4)   DXYZN(100,3),TXYZN(100,2,3)
!
	  DO I=1, MCELE
	    J=I+MCELE
	    DXYZN(J,1)= DXYZN(I,1)
	    DXYZN(J,2)= DXYZN(I,2)
	    DXYZN(J,3)=-DXYZN(I,3)

	    TXYZN(J,1,1)= TXYZN(I,1,1)
	    TXYZN(J,1,2)= TXYZN(I,1,2)
	    TXYZN(J,1,3)=-TXYZN(I,1,3)

	    TXYZN(J,2,1)= TXYZN(I,2,1)
	    TXYZN(J,2,2)= TXYZN(I,2,2)
	    TXYZN(J,2,3)=-TXYZN(I,2,3)
	  ENDDO
	  MCELE=MCELE+MCELE

	  END SUBROUTINE  XYSYM

!
!
!  
        SUBROUTINE  YZSYM(DXYZN,TXYZN,MCELE)
        IMPLICIT   NONE 
	  
	  INTEGER  I,J
	  INTEGER  MCELE

	  REAL(4)   DXYZN(100,3),TXYZN(100,2,3)
!
	  DO I=1, MCELE
	    J=I+MCELE
	    DXYZN(J,1)=-DXYZN(I,1)
	    DXYZN(J,2)= DXYZN(I,2)
	    DXYZN(J,3)= DXYZN(I,3)

	    TXYZN(J,1,1)=-TXYZN(I,1,1)
	    TXYZN(J,1,2)= TXYZN(I,1,2)
	    TXYZN(J,1,3)= TXYZN(I,1,3)

	    TXYZN(J,2,1)=-TXYZN(I,2,1)
	    TXYZN(J,2,2)= TXYZN(I,2,2)
	    TXYZN(J,2,3)= TXYZN(I,2,3)
	  ENDDO

	  MCELE=MCELE+MCELE

	  END SUBROUTINE  YZSYM
!
!
        SUBROUTINE  ZXSYM(DXYZN,TXYZN,MCELE)
        IMPLICIT   NONE 
	  
	  INTEGER  I,J
	  INTEGER  MCELE

	  REAL(4)   DXYZN(100,3),TXYZN(100,2,3)
!
	  DO I=1, MCELE
	    J=I+MCELE
	    DXYZN(J,1)= DXYZN(I,1)
	    DXYZN(J,2)=-DXYZN(I,2)
	    DXYZN(J,3)= DXYZN(I,3)

	    TXYZN(J,1,1)= TXYZN(I,1,1)
	    TXYZN(J,1,2)=-TXYZN(I,1,2)
	    TXYZN(J,1,3)= TXYZN(I,1,3)

	    TXYZN(J,2,1)= TXYZN(I,2,1)
	    TXYZN(J,2,2)=-TXYZN(I,2,2)
	    TXYZN(J,2,3)= TXYZN(I,2,3)
	  ENDDO

	  MCELE=MCELE+MCELE

	  END SUBROUTINE  ZXSYM

!
!  
        SUBROUTINE  YZZXSYM(DXYZN,TXYZN,MCELE)
        IMPLICIT   NONE 
	  
	  INTEGER  I,J
	  INTEGER  MCELE

	  REAL(4)   DXYZN(100,3),TXYZN(100,2,3)
!
	  DO I=1, MCELE
	    J=I+MCELE
	    DXYZN(J,1)= DXYZN(I,1)
	    DXYZN(J,2)=-DXYZN(I,2)
	    DXYZN(J,3)= DXYZN(I,3)

	    TXYZN(J,1,1)= TXYZN(I,1,1)
	    TXYZN(J,1,2)=-TXYZN(I,1,2)
	    TXYZN(J,1,3)= TXYZN(I,1,3)

	    TXYZN(J,2,1)= TXYZN(I,2,1)
	    TXYZN(J,2,2)=-TXYZN(I,2,2)
	    TXYZN(J,2,3)= TXYZN(I,2,3)
	  ENDDO
!
	  DO I=1, MCELE
	    J=I+2*MCELE
	    DXYZN(J,1)=-DXYZN(I,1)
	    DXYZN(J,2)=-DXYZN(I,2)
	    DXYZN(J,3)= DXYZN(I,3)

	    TXYZN(J,1,1)=-TXYZN(I,1,1)
	    TXYZN(J,1,2)=-TXYZN(I,1,2)
	    TXYZN(J,1,3)= TXYZN(I,1,3)

	    TXYZN(J,2,1)=-TXYZN(I,2,1)
	    TXYZN(J,2,2)=-TXYZN(I,2,2)
	    TXYZN(J,2,3)= TXYZN(I,2,3)
	  ENDDO
!
	  DO I=1, MCELE
	    J=I+3*MCELE
	    DXYZN(J,1)=-DXYZN(I,1)
	    DXYZN(J,2)= DXYZN(I,2)
	    DXYZN(J,3)= DXYZN(I,3)

	    TXYZN(J,1,1)=-TXYZN(I,1,1)
	    TXYZN(J,1,2)= TXYZN(I,1,2)
	    TXYZN(J,1,3)= TXYZN(I,1,3)

	    TXYZN(J,2,1)=-TXYZN(I,2,1)
	    TXYZN(J,2,2)= TXYZN(I,2,2)
	    TXYZN(J,2,3)= TXYZN(I,2,3)
	  ENDDO

	  MCELE=4*MCELE

	  END SUBROUTINE  YZZXSYM
!
!
!  
        SUBROUTINE  DOMAINFAS(DXYZN,TXYZN,MCELE,NCONT,NSIDE)
        IMPLICIT   NONE 
	  
	  INTEGER  I,J,K,N,NA,NAXIS
	  INTEGER  MCELE,NCONT(100),NSIDE(100,2)

	  INTEGER  LCT(100),LSD(100),MIN

	  REAL(4)   DXYZN(100,3),TXYZN(100,2,3)
	  REAL(4)   VEC(3),DN(100),DMIN,DOT
!
!
	  I=1  
	  NCONT(I)=1

	  VEC(1)=TXYZN(I,1,2)*TXYZN(I,2,3)-TXYZN(I,1,3)*TXYZN(I,2,2)
	  VEC(2)=TXYZN(I,1,3)*TXYZN(I,2,1)-TXYZN(I,1,1)*TXYZN(I,2,3)
	  VEC(3)=TXYZN(I,1,1)*TXYZN(I,2,2)-TXYZN(I,1,2)*TXYZN(I,2,1)

	  DOT=VEC(1)*DXYZN(1,1)+VEC(2)*DXYZN(1,2)+VEC(3)*DXYZN(1,3)

	  IF(DOT .LT. 0.0D0) THEN
	    NSIDE(I,1)=1
	    NSIDE(I,2)=2
	  ELSE
	    NSIDE(I,1)=2
	    NSIDE(I,2)=1
	  ENDIF


	  DO 100 I=2, MCELE
!	    PRINT *,' I=',I
!
!	   
          NA=0
	    DO  40 J=2, MCELE
!	      PRINT *,'     J=',J
!
!	  
	      DO N=1, I-1
	        IF(NCONT(N).EQ.J)  GOTO  40
	      END DO

!	      PRINT *,'    J=',J,' NCONT(I-1)=',NCONT(I-1)
!	      Print *,' NSIDE(I-1,2)=',NSIDE(I-1,2)
	      DO 30 K=1, 2

	        NA=NA+1
!
	        DN(NA)=SQRT((TXYZN(J,K,1)-
	1			   TXYZN(NCONT(I-1),NSIDE(I-1,2),1))**2+
     2	           (TXYZN(J,K,2)-
     2			   TXYZN(NCONT(I-1),NSIDE(I-1,2),2))**2+
     3	           (TXYZN(J,K,3)-
     3			   TXYZN(NCONT(I-1),NSIDE(I-1,2),3))**2 )	         
!
	        LCT(NA)=J
	        LSD(NA)=K
30	      CONTINUE

40	    CONTINUE

!	    PRINT *,' NA=',NA

	    DMIN=DN(1)+1.0d0
	    DO 80 N=1, NA
	      IF(DN(N) .LT. DMIN) THEN
		    DMIN=DN(N)
	        MIN=N
	      ENDIF
80	    CONTINUE

	    NCONT(I)=LCT(MIN)

	    NSIDE(I,1)=LSD(MIN) 
	    IF(LSD(MIN) .EQ. 1) THEN
	      NSIDE(I,2)=2
	    ELSE IF(LSD(MIN) .EQ. 2) THEN
	      NSIDE(I,2)=1
	    ENDIF

100	  CONTINUE

!	  PAUSE
!
	  END SUBROUTINE  DOMAINFAS

!
!  

        SUBROUTINE  ANGLE0(DXYZN,TXYZN,MCELE,NCONT,NSIDE,ANGLE)
	  IMPLICIT NONE

	  INTEGER  I,N,M
	  INTEGER  MCELE,MCNT(100),MNOD(100)
	  INTEGER  NCONT(100),NSIDE(100,2)
	  REAL(4)   ANGLE,DXYZN(100,3),TXYZN(100,2,3)
	  REAL(4)   ANG(100),VEC(3),DOT1,DOT2,SIGN,PI

	  PI=4.0d0*DATAN(1.0d0)

	  DO I=1, MCELE-1
	    N=NCONT(I)
	    M=NCONT(I+1)
	    DOT1=DXYZN(N,1)*DXYZN(M,1)+DXYZN(N,2)*DXYZN(M,2)+
	1	     DXYZN(N,3)*DXYZN(M,3)
	   
	    VEC(1)=DXYZN(N,2)*DXYZN(M,3)-DXYZN(N,3)*DXYZN(M,2)
	    VEC(2)=DXYZN(N,3)*DXYZN(M,1)-DXYZN(N,1)*DXYZN(M,3)
	    VEC(3)=DXYZN(N,1)*DXYZN(M,2)-DXYZN(N,2)*DXYZN(M,1)
	    DOT2=VEC(1)*TXYZN(N,NSIDE(I,2),1)
	1	    +VEC(2)*TXYZN(N,NSIDE(I,2),2)
	2	    +VEC(3)*TXYZN(N,NSIDE(I,2),3)

    
	    IF(ABS(DOT2) .LT. 1.0D-10) THEN
	      SIGN=1.0d0
	    ELSE
	      SIGN=DOT2/ABS(DOT2)
	    ENDIF

	    IF(DOT1  .GT. 1.0d0)   DOT1= 1.0d0
	    IF(DOT1  .LT.-1.0d0)   DOT1=-1.0d0

	    ANG(I)=SIGN*ACOS(DOT1)

	  ENDDO
	  
	  I=MCELE
	  N=NCONT(I)
	  M=NCONT(1)
	  DOT1=DXYZN(N,1)*DXYZN(M,1)+DXYZN(N,2)*DXYZN(M,2)+
	1	   DXYZN(N,3)*DXYZN(M,3)
	  VEC(1)=DXYZN(N,2)*DXYZN(M,3)-DXYZN(N,3)*DXYZN(M,2)
	  VEC(2)=DXYZN(N,3)*DXYZN(M,1)-DXYZN(N,1)*DXYZN(M,3)
	  VEC(3)=DXYZN(N,1)*DXYZN(M,2)-DXYZN(N,2)*DXYZN(M,1)   
	  DOT2=VEC(1)*TXYZN(N,NSIDE(I,2),1)
	1      +VEC(2)*TXYZN(N,NSIDE(I,2),2)
	2      +VEC(3)*TXYZN(N,NSIDE(I,2),3)


	  IF(ABS(DOT2) .LT. 1.0D-10) THEN
	    SIGN=1.0d0
	  ELSE
	    SIGN=DOT2/ABS(DOT2)
	  ENDIF
!
	  IF(DOT1  .GT. 1.0d0)   DOT1= 1.0d0
	  IF(DOT1  .LT.-1.0d0)   DOT1=-1.0d0
!
	  ANG(I)=SIGN*ACOS(DOT1)
!
	  ANGLE=2.0d0*PI

	  DO I=1, MCELE
	    ANGLE=ANGLE-ANG(I)
	  ENDDO

	 END SUBROUTINE  ANGLE0
