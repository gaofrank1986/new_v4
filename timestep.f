!
! **********************************************
!     	
!     4th order Runge_Kutta method
!
! **********************************************
!
	 SUBROUTINE Time_intg_RK4
	 USE MVAR_MOD
	 USE PVAR_MOD
!
	 IMPLICIT NONE
!	 
	 Integer INODE,JNODE,N0,NELE,NLOC,MOD,IP,IPLOT
	 Integer K
	 REAL*4  VECT(3)
	 REAL*4  AMPF(6),RAL
!
!  ============================================
!  RK1

	 RAL=180.0d0/PI

	 TimeRK=TIME
!  Computing velocity on the body surface, and normal velocity of free surface by BEM
! 
	 CALL Runge_Kutta(1)

!	 WRITE(*,*) 'RK1 COMPLETED'
!
!  ============================================
!  RK2
!
	 TimeRK=TIME+Tstep/2.0d0
! 
!  Compute wave force, and body response
!

	 CALL Runge_Kutta(2)

!	 WRITE(*,*) 'RK2 COMPLETED' 
!
!  =============================================
!  R-K3

	  TimeRK=TIME+Tstep/2.0d0
 
!
	  CALL Runge_Kutta(3)

!	  WRITE(*,*) 'RK3 COMPLETED' 
!
!  =============================================
!  R-K4
!
	  TimeRK=TIME+Tstep
! 
	 CALL Runge_Kutta(4)

!	 WRITE(*,*) 'RK4 COMPLETED'
!
!
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
       DO IP=1, NSYS
	 DO INODE=1, NNF
	    ET(INODE,IP)=ET_O(INODE,IP)+Tstep/6.0*
	1			(DH(1,INODE,IP)+2.0*DH(2,INODE,IP)+
     2		     2.0*DH(3,INODE,IP)+DH(4,INODE,IP) )

	    BKN(INODE,IP)=BKN_O(INODE,IP)+TSTEP/6.0D0*
	1		     (DP(1,INODE,IP)+2.0D0*DP(2,INODE,IP)+
     2		      2.0d0*DP(3,INODE,IP)+DP(4,INODE,IP) )

       ENDDO
       ENDDO


	 IF (NPLOUT.EQ.1) THEN
	   IPLOT=MOD(ITIME, 1)
	   IF (IPLOT .EQ. 0  )  CALL PLOTOUT8
	 ELSE IF (NPLOUT.EQ.2) THEN
	   IPLOT=MOD(ITIME, 2)
	   IF (IPLOT .EQ. 0  )  CALL PLOTOUT8
	 ELSE IF (NPLOUT.EQ.3) THEN
	   IPLOT=MOD(ITIME, 4)
	   IF (IPLOT .EQ. 0  )  CALL PLOTOUT8
	 ELSE IF (NPLOUT.EQ.4) THEN
	   IPLOT=MOD(ITIME, 8)
	   IF (IPLOT .EQ. 0  )  CALL PLOTOUT8
	 ELSE IF (NPLOUT.EQ.5) THEN
	   IPLOT=MOD(ITIME, 16)
	   IF (IPLOT .EQ. 0  )  CALL PLOTOUT8
	 ENDIF

	  FORCE=FORCEW
	  WRITE(21, 1010)  TIME, FORCE(1), FORCE(2), FORCE(3),
	1				   FORCE(4), FORCE(5), FORCE(6)

C

C Displacement and velocity of the body
C
	 DO 120 K=1, 6
 	   DSDT(K)=DSDT_O(K)+(Dposi(1,K)+2.0D0*Dposi(2,K)+
	1		        2.0D0*Dposi(3,K)+Dposi(4,K))/6.0D0
120	 CONTINUE
C

	 DO 140 K=1, 6
         DISP(K)=DISP_O(K)+
	1		 Tstep*DSDT_O(K)+Tstep*(Dposi(1,K)+Dposi(2,K)
     2		     +Dposi(3,K))/6.0D0
140	 CONTINUE

C
!

!	  IF(TIME .GT. 2.0d0*TPER) THEN

	  WRITE(22, 1010)  TIME, DISP(1), DISP(2), DISP(3),
	1             RAL*DISP(4), RAL*DISP(5), RAL*DISP(6)
C        
!
	  WRITE(23, 1010)  TIME, DSDT(1),     DSDT(2),
	2						 DSDT(3),     RAL*DSDT(4),
	3						 RAL*DSDT(5), RAL*DSDT(6)
!	  ENDIF
C
	 DO 150 K=1, 6
        IF( DISP(K). GT. 1000.0d0) DISP(K)= 1000.0d0
        IF( DISP(K). LT.-1000.0d0) DISP(K)=-1000.0d0
150	 CONTINUE


!
1010    FORMAT(F7.3,1x,7E12.4) 
!
C *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
!
10    FORMAT(1X,I6,2X,3(E15.9,3X))
21    FORMAT(1X,5(F10.5,3X))
22    FORMAT(1X,9(I7,2X))
	RETURN
	END



!
!  ==========================================
!
!
!  ==========================================
!
	 SUBROUTINE Runge_Kutta(N)
	 USE MVAR_MOD
	 USE MFUNC_mod
	 USE PVar_mod

!
	 IMPLICIT REAL*4 (A-H,O-Z)
! 
	 INTEGER INODE,N
	 REAL*4 BMAT(4),RSN(4,4)
	 REAL*4 COMP(6),AMPF(6)


       DATA RSN /1.,  1.,  1.,  1., 
     1           1., -1.,  1., -1.,
     2           1.,  1., -1., -1., 
     3           1., -1., -1.,  1./
!
!
! RAMPF: ramp function for incient potential
! RAMPV: ramp function for damping 
!
	  IF(TimeRK .LT. 2.0d0*TPER) THEN
         RAMPF=0.5d0*(1.0d0-COS(PI*TimeRK/2.0d0/TPER))
	  ELSE
	   RAMPF=1.0D0
	  END IF

	  IF(TimeRk .LT. 6.0d0*TPER) THEN
         RAMPV=0.5d0*(3.0D0+Cos(PI*TimeRK/6.0d0/TPER))
	  ELSE
	   RAMPV=1.0D0
	  END IF
!
!  ================================================ at RK1 step
!
	  IF(N  .EQ. 1)   THEN
!
!  Computing wave height, potential at water surface 
!
!
	    DISP(:)= DISP_O(:)
	    DSDT(:)= DSDT_O(:)
!
	  CALL TASSBT
!

	  DO IP=1,    NSYS
	  DO INODE=1, NNF
	    DH(1,INODE,IP)=  UNKN(INODE,IP) -DAMPF(INODE)*ET(INODE,IP)
          DP(1,INODE,IP)= -G*ET(INODE,IP) -DAMPF(INODE)*BKN(INODE,IP)
	  ENDDO  
	  ENDDO
! 
!  Compute wave force, and body response
!
	   DPDT(1:NNF,:)=DP(1,1:NNF,:)
	   DPDT(NNF+1:NNODE,:)=(UNKN(NNF+1:NNODE,:)
     1	                 -UNKN_O(NNF+1:NNODE,:))/Tstep


!         DPDT(:,:)=(UNKN(:,:)-UNKN_O(-1,:,:))/Tstep


	  CALL TFORCE			  
            UNKN_O(:,:)=UNKN(:,:)

!
	  CALL TRMAX
!	  CALL MESHT
	  CALL ACCEL(1,COMP)
C	 
	  DO K=1, 6
	    Dposi(1,K)=COMP(K)
	  END DO
!	 Print *,' Rk1  DPosi(1)=',DPosi(1,1)

!
!  ================================================ at RK2 step
!
!
	  ELSE IF(N .EQ. 2)  THEN

	   ET(:,:)      = ET_O(:,:)     +DH(1,:,:)*Tstep/2.0d0
	   BKN(1:NNF,:) = BKN_O(1:NNF,:)+DP(1,:,:)*Tstep/2.0d0

	    DISP(:)=DISP_O(:) +Tstep*DSDT_O(:)/2.0d0
	    DSDT(:)= DSDT_O(:)+Dposi(1,:)/2.0d0
C
	  CALL TASSBT
!
	  DO IP=1,    NSYS
	  DO INODE=1, NNF
	    DH(2,INODE,IP)=  UNKN(INODE,IP) -DAMPF(INODE)*ET(INODE,IP)
          DP(2,INODE,IP)= -G*ET(INODE,IP) -DAMPF(INODE)*BKN(INODE,IP)
	  ENDDO  
	  ENDDO  


! 
!  Compute wave force, and body response
!
	      DPDT(1:NNF,:)=DP(2,1:NNF,:)
	      DPDT(NNF+1:NNODE,:)=(UNKN(NNF+1:NNODE,:)
     1		                  -UNKN_O(NNF+1:NNODE,:))/(0.5*Tstep)

	  CALL TFORCE
!
	  CALL TRMAX
!	  CALL MESHT
	  CALL ACCEL(2,COMP)
C	 
	  DO K=1, 6
	    Dposi(2,K)=COMP(K)
	  END DO

!	 Print *,' Rk2  DPosi(1)=',DPosi(2,1)
!
!  ================================================ at RK3 step
!
	  ELSE IF(N .EQ. 3)  THEN

	   ET(:,:)      = ET_O(:,:)     +DH(2,:,:)*Tstep/2.0d0
	   BKN(1:NNF,:) = BKN_O(1:NNF,:)+DP(2,:,:)*Tstep/2.0d0
	  
	    DISP(:)=DISP_O(:)	+
	1		 	 Tstep*DSDT_O(:)/2.0d0+Tstep*Dposi(1,:)/4.0d0
	    DSDT(:)= DSDT_O(:)+Dposi(2,:)/2.0d0

	  CALL TASSBT
!
	  DO IP=1,    NSYS
	  DO INODE=1, NNF
	    DH(3,INODE,IP)=  UNKN(INODE,IP) -DAMPF(INODE)*ET(INODE,IP)
          DP(3,INODE,IP)= -G*ET(INODE,IP) -DAMPF(INODE)*BKN(INODE,IP)
	  ENDDO  
	  ENDDO  
!
! 
!  Compute wave force, and body response
!
	      DPDT(1:NNF,:)=DP(3,1:NNF,:)
	      DPDT(NNF+1:NNODE,:)=(UNKN(NNF+1:NNODE,:)
     1		                 -UNKN_O(NNF+1:NNODE,:))/(0.5*Tstep)

!         DPDT(:,:)=(UNKN(:,:)-UNKN_O(-1,:,:))/(1.5*Tstep)

	  CALL TFORCE	  
C
	  CALL TRMAX
!	  CALL MESHT
	  CALL ACCEL(3,COMP)
C	 
	  DO K=1, 6
	    Dposi(3,K)=COMP(K)
	  END DO
C
!	 Print *,' Rk3  DPosi(1)=',DPosi(3,1)
!
!  ================================================ at RK4 step
!
	  ELSE IF(N .EQ. 4)  THEN

	   ET(:,:)      = ET_O(:,:)     +DH(3,:,:)*Tstep
	   BKN(1:NNF,:) = BKN_O(1:NNF,:)+DP(3,:,:)*Tstep

	    DISP(:)=DISP_O(:)	+
	1		 	 Tstep*DSDT_O(:)+Tstep*Dposi(2,:)/2.0d0
	    DSDT(:)= DSDT_O(:)+Dposi(3,:)
!
!	   HEIGHT(4,:,:)=HEIGHT(4,:,:)+DH(3,:,:)*Tstep
!	   PFREEN(4,:,:)=PFREEN(4,:,:)+DP(3,:,:)*Tstep
!	   BKN(:,:) = BKN_O(:,:)

	  CALL TASSBT
!
	  DO IP=1,    NSYS
	  DO INODE=1, NNF
	    DH(4,INODE,IP)=  UNKN(INODE,IP) -DAMPF(INODE)*ET(INODE,IP)
          DP(4,INODE,IP)= -G*ET(INODE,IP) -DAMPF(INODE)*BKN(INODE,IP)
	  ENDDO  
	  ENDDO  
!
! 
!  Compute wave force, and body response
!
	      DPDT(1:NNF,:)=DP(4,1:NNF,:)
	      DPDT(NNF+1:NNODE,:)=(UNKN(NNF+1:NNODE,:)
     1		                 -UNKN_O(NNF+1:NNODE,:))/Tstep

!        DPDT(:,:)=(UNKN(:,:)-UNKN_O(-1,:,:))/(2.0*Tstep)

	  CALL TFORCE

			  
C
	  CALL TRMAX
!	  CALL MESHT

	  CALL ACCEL(4,COMP)
C	 
	  DO K=1, 6
	    Dposi(4,K)=COMP(K)
	  END DO
!
!	 Print *,' Rk4  DPosi(1)=',DPosi(4,1)

	  ENDIF
!
!
!	 pause
	 END
!
! ====================================================
!
!
!
!
	 SUBROUTINE PLOTOUT8
	  USE MVAR_MOD
	  USE MFUNC_mod
      use io
        IMPLICIT   NONE  

	 !CHARACTER*16 NAME
       character(:),allocatable :: filename
	 CHARACTER*6  FIRST

	 INTEGER  I,INODE,NE

	 CALL OPENFILE(ITIME,FIRST)


	  Print *,' Itime=',itime,' First=',first
      

	 !NAME='OUTtime\WAVE'//FIRST//'.DAT'
	 !NAME='OUTtime\WAVE'//FIRST//'.DAT'
       filename = getfilename("./wave/wave_elev",itime,'.dat')
       OPEN(102,FILE=filename,STATUS='UNKNOWN')

!        OPEN(9,  FILE='OUTPUT\OUTPUT1.txt',    STATUS='UNKNOWN')
        
        

        WRITE(102,*) 'TITLE = "3D Mesh Grid Data for Element Boundary"'
         WRITE(102,*) 'VARIABLES = "X", "Y", "Z"'

        !WRITE(102,*) 'ZONE N=',NNF,',','E=',NELEMF,',','F=FEPOINT,ET=QUADRILATERAL'
      
      !write(102,'(a,i0,a,i0,a)') "ZONE N=",NNF,",E=",NELEMF,",F=FEPOINT,ET=QUADRILATERAL"

      write(102,*) "ZONE N=",NNF,",E=",NELEMF
      write(102,*) ",F=FEPOINT,ET=QUADRILATERAL"
      write(102,*) 'SolutionTime=',itime

      write(102,*) 'StrandID=1'


       DO INODE=1, NNF
	  WRITE(102,21)     (XYZ(I,INODE), I=1, 2), 
	1		        ET(INODE,1)+ETI(XYZ(1,INODE),XYZ(2,INODE))
       ENDDO

	 DO NE=1,  NELEMF
        WRITE(102,22) NCON(NE,1),NCON(NE,3),NCON(NE,5),NCON(NE,7)
	 ENDDO
!
	 Close(102)
!
10     FORMAT(1X,I6,2X,3(E15.9,3X))
21     FORMAT(1X,5(F10.5,3X))
22     FORMAT(1X,9(I7,2X))
!
	 END

!
! ====================================================
!

	 SUBROUTINE PLOTOUT9
	  USE MVAR_MOD

        IMPLICIT   NONE  

	 CHARACTER*14 NAME,NAME1
	 CHARACTER*6  FIRST

	 INTEGER  I,INODE,NE

	  CALL OPENFILE(ITIME,FIRST)

	 NAME='WAVE'//FIRST//'.DAT'
       OPEN(102,FILE=NAME,STATUS='UNKNOWN')

	 WRITE(102,*) 'TITLE = "3D Mesh Grid Data for Element Boundary"'
	 WRITE(102,*) 'VARIABLES = "X", "Y", "Z"'
	 WRITE(102,*) 'ZONE N=',NNTCH,',','E=',NELEMF*4,',',
     &	'F=FEPOINT,ET=QUADRILATERAL'

       DO INODE=1, NNTCH
	  WRITE(102,21)     (XYZ(I,INODE), I=1, 3)
       ENDDO

	 DO NE=1,  NELEMF
        WRITE(102,22) NCON(NE,1),NCON(NE,2),NCON(NE,9),NCON(NE,8)
        WRITE(102,22) NCON(NE,2),NCON(NE,3),NCON(NE,4),NCON(NE,9)
        WRITE(102,22) NCON(NE,4),NCON(NE,5),NCON(NE,6),NCON(NE,9)
        WRITE(102,22) NCON(NE,8),NCON(NE,9),NCON(NE,6),NCON(NE,7)
	 ENDDO

	 Close(102)
!
10     FORMAT(1X,I6,2X,3(E15.9,3X))
21     FORMAT(1X,5(F10.5,3X))
22     FORMAT(1X,9(I7,2X))
!
	 END


C =========================================

	 SUBROUTINE OPENFILE(MTIME,FIRST)
	!INTEGER MTIME
	 CHARACTER*6 FIRST
	 
	 IF(MTIME.LT.10) THEN
	  MN0=MTIME
	  FIRST=CHAR(MN0+48)
	 ELSE IF(MTIME.LT.100) THEN
	  MN0=MOD(MTIME,10)	
	  MN1=(MTIME-MN0)/10
	  FIRST=CHAR(MN1+48)//CHAR(MN0+48)
	 ELSE IF(MTIME.LT.1000) THEN
	  MN01=MOD(MTIME,100)
	  MN0=MOD(MN01,10)	
	  MN1=(MN01-MN0)/10
	  MN2=(MTIME-10*MN1-MN0)/100
	  FIRST=CHAR(MN2+48)//CHAR(MN1+48)//CHAR(MN0+48)
	 ELSE IF(MTIME.LT.10000) THEN 
	  MN012=MOD(MTIME,1000)
	  MN01=MOD(MN012,100)
	  MN0=MOD(MN01,10)
	  MN1=(MN01-MN0)/10
	  MN2=(MN012-10*MN1-MN0)/100
	  MN3=(MTIME-MN012)/1000
	  FIRST=CHAR(MN3+48)//CHAR(MN2+48)//CHAR(MN1+48)//CHAR(MN0+48)
	 ENDIF
!
! *-*-*--*
!
	 RETURN
	 END

