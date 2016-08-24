C
C *************************************************
C
C  Compute acceleration of body motion
C  NUM :
C DISP : Displcacement,  Input data
C DSDT : Velocity of body motion, Input data
C COMP :                 Output data
C *************************************************
C
	 SUBROUTINE ACCEL(NUM,COMP)

       USE MVar_MOD
	 USE FCVar_MOD
       USE PVar_MOD
       IMPLICIT NONE
C
       INTEGER I,J,K,NUM,IPrint,IMOD,ICTIM,IDTIM,NPRI
	 REAL*4  CABLEF(6),FenderF(6)

	 REAL*4  FORCES(6),FORCEV(6),FORCEK(6),
	1		 FORCEC(6),FORCED(6)
	 REAL*4  COMP(6)
	 REAL*4  FCABL(50),SCABL(50),FFender(50),SFender(50)
!
!  -------------------------------------------
!
	 DO 10 K=1, 6
	  FORCE(K) =FORCE0(K)+FORCEW(K)
10	 CONTINUE
c
c  Restoring force
c
	  FORCES=0.0d0
	 DO 22 K=1, 6
	 DO 20 J=1, 6
	  FORCE(K) =FORCE(K) -CRS(K,J)*DISP(J)
	  FORCES(K)=FORCES(K)-CRS(K,J)*DISP(J)
20	 CONTINUE
22	 CONTINUE
!
c
c  Force for elastic springs
c
	  FORCEK=0.0d0
	 DO  K=1, 6
	  DO  J=1, 6
	   FORCE(K) =FORCE(K) -STKM(K,J)*DISP(J)
	   FORCEK(K)=FORCEK(K)-STKM(K,J)*DISP(J)
	  ENDDO
	 ENDDO
!
!

!  Viscous damping force
!
	  FORCEV=0.0d0
	 DO 28 K=1, 6
	 DO 26 J=1, 6
	  FORCE(K) =FORCE(K) -RAMPV*VISC(K,J)*DSDT(J)
	  FORCEV(K)=FORCEV(K)-RAMPV*VISC(K,J)*DSDT(J)

!	  FORCE(K) =FORCE(K) -RAMPV*VISC(K,J)*DXDT(J)*DABS(DXDT(J))
!	  FORCEV(K)=FORCEV(K)-RAMPV*VISC(K,J)*DXDT(J)*DABS(DXDT(J))

26	 CONTINUE
28	 CONTINUE
C
C  Fasting force from cables
C
!	 CALL CFCABLE(CABLEF,DISP,FCABL,SCABL)
	
	CABLEF(:)=0.0
	
	 DO 30 K=1, 6
	 FORCE(K)=FORCE(K)+CABLEF(K)
	 FORCEC(K)=CABLEF(K)
30	 CONTINUE
C
C  Acting force from fenders
C
!	 CALL CFFender(FenderF,FFender,SFender)

       FenderF(:)=0.0
	 DO 40 K=1, 6
	 FORCE(K)=FORCE(K)+FenderF(K)
	 FORCED(K)=FenderF(K)
40	 CONTINUE
C
C ===========================================
C
       COMP=0.0D0
	DO 101 K=1, 6
	 DO 100 J=1, 6
        COMP(K)=COMP(K)+TRMAS(K,J)*FORCE(J)
100	 CONTINUE
101	 CONTINUE
C
!	 Print *,' FORCE(1)=', FORCE(1)
!	 Print *,' FORCE(3)=', FORCE(3)
!	 Print *,' FORCE(5)=', FORCE(5)
!	 Print *,' Tstep=', Tstep

        COMP(:)=COMP(:)*Tstep

	 
!	 Print *,' COMP(1)=',COMP(1)
!	 Print *,' COMP(3)=',COMP(3)
!	 Print *,' COMP(5)=',COMP(5)


C
C Time history of displacement, velocity and different kind forces	 
C

!	  IF(TimeRK .GT. 2.0d0*TPER) THEN

	  NPRI=2.0*PI/W1/Tstep/16.0
	  IF(NPRI .LE. 0) NPRI=1
	  IMOD=MOD(ITIME, NPRI)
! 每IMOD步输出一次
!
	 IF(IMOD. EQ. 0 .and. NUM .EQ. 1) THEN

	  WRITE(31, 1010)  TIMERK, DISP(1), DSDT(1), FORCEW(1), FORCES(1), 
	1					    FORCEK(1),FORCEV(1), FORCEC(1), FORCED(1)
	  WRITE(32, 1010)  TIMERK, DISP(2), DSDT(2), FORCEW(2), FORCES(2),
	1					    FORCEK(2),FORCEV(2), FORCEC(2), FORCED(2)
	  WRITE(33, 1010)  TIMERK, DISP(3), DSDT(3), FORCEW(3), FORCES(3), 
	1					    FORCEK(3),FORCEV(3), FORCEC(3), FORCED(3)
	  WRITE(34, 1010)  TIMERK, DISP(4), DSDT(4), FORCEW(4), FORCES(4),
	1					    FORCEK(4),FORCEV(4), FORCEC(4), FORCED(4)
	  WRITE(35, 1010)  TIMERK, DISP(5), DSDT(5), FORCEW(5), FORCES(5), 
	1					    FORCEK(5),FORCEV(5), FORCEC(5), FORCED(5)
	  WRITE(36, 1010)  TIMERK, DISP(6), DSDT(6), FORCEW(6), FORCES(6),
	1					    FORCEK(6),FORCEV(6), FORCEC(6), FORCED(6)
C
C time history of fender force
C  	  
	 IDTIM=(NOFender+4)/5
	 IF(IDTIM .GT. 5) IDTIM=5
	 DO I=1, IDTIM
	  J=(I-1)*5
	  IPRINT=40+I
	  WRITE(IPRINT,   1010)  TIMERK, FFender(J+1), FFender(J+2),
	1					     FFender(J+3), FFender(J+4), FFender(J+5)
	  WRITE(IPRINT+5, 1010)  TIMERK, SFender(J+1), SFender(J+2),
	1					     SFender(J+3), SFender(J+4), SFender(J+5)

	  END DO
C
C Time history of cable tensile force 
C
	 ICTIM=(NOCABLE+4)/5
	 IF(ICTIM .GT. 5) ICTIM=5
	 DO I=1, ICTIM
	  J=(I-1)*5
	  IPRINT=50+I
	  WRITE(IPRINT,   1010)  TIMERK, FCABL(J+1), FCABL(J+2),
	1					     FCABL(J+3), FCABL(J+4), FCABL(J+5)
	  WRITE(IPRINT+5, 1010)  TIMERK, SCABL(J+1), SCABL(J+2),
	1					     SCABL(J+3), SCABL(J+4), SCABL(J+5)

	  END DO

	  END IF
C
!	 END IF
!
1010    FORMAT(F7.3,1x,8E12.4) 
C
	 END 
