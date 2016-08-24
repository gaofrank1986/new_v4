C  
C **************************************************************** 
C *                                                              * 
C *  Programme for the first order diffraction and radiation     * 
C *  wave analysis of 3-dimensional bodies                       * 
C *															   *
C * Revised 											           *
C *  4 into NSYS	 Sept.16, 2000                                 *
C *                                                              * 
C **************************************************************** 
C 
        PROGRAM WVTDUT1
	  USE MVAR_MOD
	  USE PVAR_MOD
      use io
 
!
        IMPLICIT   NONE  
	  INTEGER I,IFWKO,IS,IND,M,IP,IPOL,K,NDMAX
	  INTEGER NTnum,IFLAG_T
!
	  REAL*4  WL,AMFJ,R,EX(4),Alpha,WAVES  
	  REAL*4  FAMPR(6),FAMPI(6),FORCER(6),FORCEI(6)
	  REAL*4  PL_AMP(6),FORAMP
	  REAL*4  FCD_AMR,  FCD_AMI
!
! IFLAG_T,FCD_AMR, FCD_AMI: not used in this program
!
        DATA EX /  1.,  1., -1., -1./ 
!
! FAMPR, FAMPI:   real and imaginary parts of body motion 
!                 from frequency domain calculation
! FORCER, FORCEI: real and imaginary parts of wave force 
!                 from frequency domain calculation
! PL_AMP:         amplitude of plotting variable
C
C ----------------------------------------
C Input data files
        fd = create_folder()
        OPEN(1, FILE='INPUT/DATIN.txt',      STATUS='OLD') 
        OPEN(2, FILE='INPUT/DATBDMS.txt',    STATUS='OLD') 
        OPEN(3, FILE='INPUT/DATFSMS.txt',    STATUS='OLD') 
        OPEN(7, FILE='INPUT/DATMASS.txt',    STATUS='OLD')  
!
!
! -----------------------------------
!  Output data files
!
        OPEN(9,  FILE='OUTPUT/OUTPUT1.txt',    STATUS='UNKNOWN')
	  OPEN(10, FILE='OUTPUT/OUTPUT.txt' ,    STATUS='UNKNOWN')
	  OPEN(101,FILE='OUTPUT/OUTAMT.txt' ,    STATUS='UNKNOWN')
	  OPEN(102,FILE='OUTPUT/OUTBMT.txt' ,    STATUS='UNKNOWN')
	  OPEN(103,FILE='OUTPUT/OUTELA.txt' ,    STATUS='UNKNOWN')
C                                                                              

C
C ---------------------------------------- 
C
C  IFWKO=0 : input wave number, otherwise input wave frequency
C  H<0: infinite water depth
C  Nwave: Number of waves to simulate
C
        READ(1,*)      IFWKO

        IF (IFWKO .EQ. 0)  THEN
          READ(1,*)      H, AMP, WK, BETA
          IF (H .LE. 0.0D0) THEN
            W1=SQRT(G*WK)
          ELSE
            W1=SQRT(G*WK*TANH(WK*H))
          END IF
        ELSE
          READ(1,*)      H, AMP, W1, BETA
          IF(H .LE. 0.0D0) THEN
            WK=W1**2/G
          ELSE
            CALL WAVECK(W1,H,WK)				 ! Compute wave number
          END IF
        END IF
!
	  READ(1,*)  FCD_AMR,  FCD_AMI
	  READ(1,*)  IFLAG_T,WAVES,  NTnum
	  READ(1,*)  NPLOUT    
!
!  Waves: number of waves to simulate (may be 10 or 0.5)
!  NTnum: number of time steps in one wave
!  NPLOUT=1-5, Output the wave profile with step intervals 

	  PI4=4.0*PI
        IORDER=1
        WRITE(6,1113)
!
! 
        TPER=2.*PI/W1 
        BETA=BETA*PI/180.0D0 
C 
        V=W1*W1/G 
	  WL=2.0D0*PI/WK
C 
	  TStep=TPer/NTnum
	  NTIME=INT(WAVES*TPER/TStep)
        
        Print *,' Ntime=', Ntime


        WRITE(6,*) 
        WRITE(9,*)
	  WRITE(6,*) '                   ================='
        WRITE(9,*) '                   ================='
C
        WRITE(6,1111)  H,AMP,WK,V,WL,W1,TPER,BETA*180./PI 
        WRITE(9,1111)  H,AMP,WK,V,WL,W1,TPER,BETA*180./PI 


!
!  ------------------------------------
!
        CALL BODMASS 				            ! Read in data of body mass, etc.
        WRITE(6,*),'  After BODMASS' 
C 
        READ(2,*) ISYS 
        READ(2,*)   NELEMB, NNB, NNBD, IPOL
!
   	 ALLOCATE (NCONB(NELEMB,8),NCONDB(NELEMB,8))
       ALLOCATE (XYZB(3,NNB),DXYZB(3,NNBD))

!	  INTEGER NELEM,NELEMB,NELEMF,NNODE,NNODED,NNB,NNF,NNTCH

        IF(ISYS.EQ.0) NSYS=1
        IF(ISYS.EQ.1) NSYS=2
        IF(ISYS.EQ.2) NSYS=4
C

        READ(3,*)   NELEMF
!
	  WRITE(6,*) ' ISYS=',ISYS,' NSYS=',NSYS
	  WRITE(6,*) ' NELEMB=',NELEMB,' NELEMF=',NELEMF

	  WRITE(10,*) ' ISYS=',ISYS,' NSYS=',NSYS
	  WRITE(10,*) ' NELEMB=',NELEMB,' NELEMF=',NELEMF
        WRITE(10,*)
        
	  NELEM=NELEMB+NELEMF
!
	  WRITE(6,*) ' NELEM=',NELEM,'  IOPL=',IPOL
	  WRITE(10,*) ' NELEM=',NELEM,'  IOPL=',IPOL
        WRITE(10,*)
        
        
        ALLOCATE(SAMB(NELEM,16,0:8),SAMBXY(NELEM,16,3),
	1		   DSAMB(NELEM,16,6),NCN(NELEM),NCON(NELEM,8),
     1		   NCOND(NELEM,8),IETYPE(NELEM),NNORMN(8*NELEM) )
        ALLOCATE( XYZE(3,8,NELEM),DXYZE(3,8,NELEM),DAMPE(8,NELEM))
        ALLOCATE( TXYZE(3,8,NELEM))
	  ALLOCATE( XYZTP(3,8*NELEM),DXYZTP(3,8*NELEM),DAMPTP(8*NELEM))
!
! --------------------------------------------
!
        CALL MESHFS4   		        ! Read in data on free surface mesh
        WRITE(6,*),'  After MESHFS4' 
        WRITE(10,*),'  After MESHFS4' 
        WRITE(10,*)
                
        CALL MESHBD(IPOL) 		    ! Read in data on body mesh
        WRITE(6,*),'  After MESHBD4' 
        WRITE(10,*),'  After MESHBD4'!
        WRITE(10,*)	  
	  
	  CLOSE(2)
        !OPEN(50, FILE='OUTPUT\DATBDMS.txt',    STATUS='UNKNOWN')

	  CALL CONVSB
	  
        WRITE(6,*),'  After CONVSB          NNODE=',NNODE 
        WRITE(10,*),'  After CONVSB          NNODE=',NNODE !
        WRITE(10,*)

	  NDMAX=MAX(NNODE,NNODED)
       ALLOCATE(NODELE(NNODE,64),NODNOE(NNODE),NODELJ(NNODE,64),
     1         NODQUA(NNODE),NNORMC(NNODED),
     2		 XYZ(3,NDMAX),DXYZ(6,NDMAX),ANGLE(NNODE),
     2		 DAMPF(NNF))
!
! ---------------------------------------------------
!
       WRITE(10,*)
	 WRITE(10,*) '  IND       X          Y          Z        DAMPing'
	 DO IND=1, NNF
	  DAMPF(IND)=DAMPTP(IND)
	  XYZ(1,IND)=XYZTP(1,IND)
	  XYZ(2,IND)=XYZTP(2,IND)
	  XYZ(3,IND)=XYZTP(3,IND)
	  WRITE(10,111) IND, XYZ(1,IND),XYZ(2,IND),XYZ(3,IND),DAMPF(IND)
	 END DO
!	  	 
	 DO IND=NNF+1, NDMAX
	  XYZ(1,IND)=XYZTP(1,IND)
	  XYZ(2,IND)=XYZTP(2,IND)
	  XYZ(3,IND)=XYZTP(3,IND)
	  WRITE(10,111) IND, XYZ(1,IND),XYZ(2,IND),XYZ(3,IND)
	 END DO
!

       WRITE(10,*)
	 WRITE(10,*) '  IND       DX         DY         DZ'

	 DO IND=1, NNODED
	  DXYZ(1,IND)=DXYZTP(1,IND)
	  DXYZ(2,IND)=DXYZTP(2,IND)
	  DXYZ(3,IND)=DXYZTP(3,IND)
	  NNORMC(IND)=NNORMN(IND)
	  WRITE(10,111) IND, DXYZ(1,IND),DXYZ(2,IND),DXYZ(3,IND)
	 END DO
!
111	 FORMAT(I4,5(2x,F12.5))
!
	 DEALLOCATE(XYZTP,DXYZTP,NNORMN)
!
	 DAMPF(:)=W1*DAMPF(:)
!
!
!  --------------------------------------------
!
       ALLOCATE(AMATA(NNODE,NNODE,NSYS),CMATA(NNODE,NNODED,NSYS),
	1		  BMATA(NNODE,NSYS), INDX(NNODE,NSYS))
       ALLOCATE(UNKN(NNODE,NSYS),  BKN(NNODED,NSYS),
	1		  UNKN_O(NNODE,NSYS),BKN_O(NNODED,NSYS),
	1	      ET(NNF,NSYS),    ET_O(NNF,NSYS),DPDT(NNODE,NSYS))
!	 ALLOCATE(HEIGHT(4,NNF,NSYS),PFREEN(4,NNF,NSYS))
	 ALLOCATE(DH(4,NNF,NSYS),DP(4,NNF,NSYS),Dposi(4,6))
!
!      
! Identify the Gaussian sampling points and evaluate the    
! Corresponding values for the shape functions and Jacobian matrices       
!
        CALL BODYFD                  
        WRITE(6,*)  'AFTER BODYFD' 
!   
! Assembling matrix and computing diffraction and radiation potentials
!
       Print *,' Ntime=', Ntime
 
        CALL TASSB0        
        WRITE(6,*) 'AFTER ASSEMB0' 
! 
	  ITIME=0
        TIME=0.0d0
!	  
	  ET(:,:)  =0.0
	  BKN(:,:) =0.0
	  UNKN(:,:)=0.0
	  UNKN_O(:,:)=0.0
	  FORCE(:) =0.0
	  DISP(:)  =0.0
	  DSDT(:)  =0.0
!
!	  CALL PLOTOUT8
!
!  ==================================================
!
	  WRITE(21, *)  ' TIME   FORCE(1)   FORCE(2)   FORCE(3)',
     1				     '   FORCE(4)   FORCE(5)   FORCE(6)'
	  WRITE(22, *)  ' TIME    Dis(1)     Dis(2)     Dis(3)',
     1		              '   Dis(4)     Dis(5)     Dis(6)'
	  WRITE(23, *)  ' TIME    Vel(1)     Vel(2)     Vel(3)',
     1		              '   Vel(4)     Vel(5)     Vel(6)'
!
	  WRITE(31, *)  'TIMERK   DISP   DSDT    FORCEW    FORCES', 
	1		'	   FORCEK   FORCEV   FORCEC   FORCED'
	  WRITE(32, *)  'TIMERK   DISP   DSDT    FORCEW    FORCES', 
	1		'	   FORCEK   FORCEV   FORCEC   FORCED'
	  WRITE(33, *)  'TIMERK   DISP   DSDT    FORCEW    FORCES', 
	1		'	   FORCEK   FORCEV   FORCEC   FORCED'
	  WRITE(34, *)  'TIMERK   DISP   DSDT    FORCEW    FORCES', 
	1		'	   FORCEK   FORCEV   FORCEC   FORCED'
	  WRITE(35, *)  'TIMERK   DISP   DSDT    FORCEW    FORCES', 
	1		'	   FORCEK   FORCEV   FORCEC   FORCED'
	  WRITE(36, *)  'TIMERK   DISP   DSDT    FORCEW    FORCES', 
	1		'	   FORCEK   FORCEV   FORCEC   FORCED'
!
! ======================================
!

	  DO 500 ITIME=0, NTIME
       WRITE(10,*)
	 WRITE(10,*) '  Itime=',Itime

	   TIME=ITIME*Tstep
	   IF(MOD(ITIME,10) .EQ. 0)  THEN
           WRITE(6,1115) ITIME,Time     
	   ENDIF
!
	   BKN_O(:,:) = BKN(:,:)
	   ET_O(:,:)  = ET(:,:)  

	   DISP_O(:) =DISP(:)
	   DSDT_O(:)= DSDT(:)

	   FORCE_O=FORCE
!	   
	   CALL Time_intg_RK4
!

C 
500	 CONTINUE      
C 
C =================
C 
       DEALLOCATE(AMATA,CMATA,BMATA,INDX)

C                              
1010    FORMAT(F7.3,1x,F7.3,1x,6E14.5) 
C 
1111    FORMAT(//,'  WATER DEPTH=',F9.3,'    WAVE AMPLITUDE=', F6.2,/,
     1    '  WAVE NUMBER=',F9.5,'  K0=',F9.5,'  WAVE LENGTH=',F9.4,/, 
     3    '  ANGULAR FREQU.=',F9.5,'   WAVE PERIOD=',F7.3,/,      
     2    '  WAVE DIRECTION:',F7.3,'  Degree',/)
C 
1113    FORMAT(/,15x,' FIRST  ORDER PROBLEM')  
1114    FORMAT(/,15x,' SECOND ORDER PROBLEM') 
1115	  FORMAT(' I_time=',I5,'    at the time:',F10.5,' s')
1200	FORMAT(2x,I3,3F12.5,1x,2F13.5)
C 
        STOP 
        END      
 
 

C 
C
C
C  *************************************************
C  *    The subroutine computes wave number        *
C  *  by wave frequency and water depth.           *
C  *    For infinite water depth, give negative H. *
C  *************************************************

        SUBROUTINE WAVECK(SIGMA,H,WK)
        IMPLICIT  NONE

	  REAL*4 SIGMA,H,WK,B,G,A,Y,C
C
C  H: WATER DEPTH    SIGMA: WAVE FREQUENCY
C  C: WAVE CELERITY  Wk: wave number
C
        DATA G/9.807D0/
C
        IF( SIGMA .LE. 0.0D0 )  THEN
	    Print *,' IN WAVECK:  W=',SIGMA
		STOP
	
        ELSE 

          IF( H .GT. 0.0D0) THEN
           B = G * H
           Y = SIGMA * SIGMA *H/G
  12       A=1.0D0/(1.0D0+Y*(0.66667D0+Y*(0.35550D0+
     1       Y*(0.16084D0+Y*(0.063201D0+Y*(0.02174D0+
     2       Y*(0.00654D0+Y*(0.00171D0+
     2       Y*(0.00039D0+Y*0.00011D0)))))))))
  22       C=SQRT(B/(Y+A))
           WK=SIGMA/C
           ELSE IF( H .LE. 0.0D0) THEN
           WK=SIGMA**2/G
        END IF
        RETURN
C

	 END IF


        RETURN
        END

