C **************************************************
C 												
C   This is a module for declaring variables  
C                                                 
C       Nov. 3, 2000           by   Bin Teng		
C **************************************************
C
C  MODULES: MVAR_MOD + PVAR_MOD + TRVAR_MOD
C
C  MVAR_MOD : for constants, wave parameter and mesh
C  BVAR_MOD : for B-spline function
C  PVAR_MOD : for potentials, forces, and body mass, etc
C  TRVAR_MOD: for finer mesh from tri-pole transformation
C
C ===================================================
C   Constants and Variables for body mesh
C
        MODULE MVar_mod
C
        INTEGER NTIME,ITIME,ISYS,NSYS,IORDER
	  INTEGER NELEM,NELEMB,NELEMF,NNODE,NNODED,NNB,NNBD,NNF,NNTCH

	  INTEGER NPLOUT
!
!  NPLOUT=1, Output the wave profile at every time step; =2, every 2 time steps;
!        =3, every 4 time step;  =4, every 8 time steps; =5, every 16 time steps.
!
! NTIME:  total time steps for simulation
! ITIME:  present time step 
! ISYS:   number of symmetric planes
! IORDER: perturbation order of the problem 
! NELEM:  number of total elements 
! NELEMB: number of elements on body surface
! NELEMF: number of elements on the free surface
! NNODE:  total number of nodes according to the coordinate
! NNODED: total number of nodes according to the normals
! NNB:    number of nodes on the body surface according to coordinate
! NNBD:   number of nodes on the body surface according to directives
! NNF:    number of nodes on the free surface
! NNTCH:  number of nodes and the centers of element (For plotting by Techplot)
! 
        character(:),allocatable :: fd !output folder

        REAL*4 G,RHO,PI,PI4     
        REAL*4 H,AMP,BETA,W1,V 
        REAL*4 WK,TPER
	  REAL*4 Tstep,TIME,TimeRK,RAMPF,RAMPV

! H   : water depth
! Amp : amplitude of incident waves
! BETA: incident angle
! W1  : angular frequency of incident waves
! V   : wave number in deep water
! WK  : wave number
! TPER: wave period
! Tstep: time step
! Time : simulation time at each time step
! TimeRK: simulation time at each RK step
! RampF: ramp function for incident potential
! RampV: ramp function for damping
!    
   	 INTEGER, ALLOCATABLE:: NCN(:),IETYPE(:),NCON(:,:),NCOND(:,:),
	1                        NNORMC(:)
       INTEGER, ALLOCATABLE:: NODELE(:,:),NODNOE(:),
	1			            NODELJ(:,:),NODQUA(:)
   	 INTEGER, ALLOCATABLE:: NCONB(:,:),NCONDB(:,:)

! NCN: number of nodes in the element
! IETYPE: type of the element; =1, on body surface; =2, on free surface
! NCON: 
! NCOND:
! NNORMC:
! NODELE:
! NODNOE:
! NODELJ:
! NODQUA:
!
       REAL*4,ALLOCATABLE:: XYZE(:,:,:),DXYZE(:,:,:),XYZ(:,:),DXYZ(:,:),
	1                      DAMPE(:,:),DAMPF(:)
       REAL*4,ALLOCATABLE:: XYZB(:,:),DXYZB(:,:)

       REAL*4,ALLOCATABLE:: TXYZE(:,:,:)
!
! XYZE  : Initial Coordinates of nodes of body mesh
! DXYZE:
! XYZ:
! DXYZ:
! DAMPE:
! DAMPF:
! TXYZE : Coordinates of nodes of body mesh at the simulation time 
!
	 REAL*4,ALLOCATABLE:: BKN(:,:),BKN_O(:,:),UNKN(:,:),UNKN_O(:,:)
	 REAL*4,ALLOCATABLE:: ET(:,:),ET_O(:,:),DPDT(:,:)
!	 REAL*4,ALLOCATABLE:: BKN(:,:),UNKN(:,:),BKN_O(:,:),UNKN_O(:,:,:),
!	1					  ET(:,:),ET_O(:,:),DPDT(:,:)
!
! BKN:  1:NNF, Potential; NNF+1:NNODE, normal derivative
! UNKN: 1:NNF, normal derivative; NNF+1:NNODE, potential 
! BKN_O: store one-step data
! UNKN_O: store two-step data
! ET: Wave profile
! ET_O: store one-step data
! DPDT: time difference of potential

!
       REAL*4,ALLOCATABLE:: SAMB(:,:,:),SAMBXY(:,:,:),
	1	                  DSAMB(:,:,:),ANGLE(:)
!
! SAMB: 
! SAMBXY: Coordinates of Gaussin points
! DSAMB:  Normal direvatives at Gaussian points
! ANGLE:  Solid angle
!
	 REAL*4,ALLOCATABLE:: DH(:,:,:),DP(:,:,:),Dposi(:,:)
!
! DH:
! DP:
! Dposi:
!
!  For linear equations
!
	 INTEGER, ALLOCATABLE:: INDX(:,:)
       REAL*8,ALLOCATABLE::   AMATA(:,:,:),CMATA(:,:,:),BMATA(:,:)


       DATA G,PI,RHO/9.807,3.14159265359,
	1	               1023.0/  

! ------------------------------------
! Temporary arrayes
! 
   	 INTEGER, ALLOCATABLE:: NNORMN(:)
       REAL*4,  ALLOCATABLE:: XYZTP(:,:),DXYZTP(:,:),DAMPTP(:)
!
       END MODULE MVar_mod


C
C Variable for potential, force and body motion
C =============================================
C
        MODULE PVar_mod
C
	  INTEGER  NNT
	  PARAMETER (NNT=2000)    

        REAL*4 XC,YC,ZC,XTC,YTC,ZTC
        REAL*4 ARE,XF,YF,XK2,YK2,XCF
        REAL*4 VOLM,XB,YB,ZB

        REAL*4 AMAS(6,6),BDMP(6,6)  
        REAL*4 RMAS(6,6),CDMP(6,6),CRS(6,6),STKM(6,6),XIA(3,3)
C
        REAL*4 FORCEW(6),FORCE0(6),FORSCD(6),AMPJ(6)
	  REAL*4 FORCE(6),FORCE_O(6)
C
        REAL*4 TRMAS(6,6),VISC(6,6)

        REAL*4  DISP(6),DSDT(6),DISP_O(6),
	2		  DSDT_O(6),DSDDTL(6),TRXYZ(3,3)   

        REAL*4  RESPR(6),RESPI(6)  


        END MODULE PVar_mod

C
C =========================================
C  Varibles used in the Tri-pole transform 
C							  
       MODULE TRVar_mod

	 INTEGER NOSAMP
	 REAL*4 XYNOD(3,50),DXYNOD(6,50),SAMNOD(50,0:8)

       END MODULE TRVar_mod

C
C  **************************************************
C  *												  *
C  *  This is a module for declaring variables      *
C  *      in B-spline expansion                     *
C  *                                                *
C  *      Aug. 10, 2002          by   Bin Teng      *
C  **************************************************				
C
       MODULE  BVAR_mod
C
       INTEGER KBSPL,NB,NT,NA,NAA,NBU,NA1
       INTEGER MUSTA,MUEND	 
C
C --------------------------------------------------------------------
C KBSPL : DEGREE (order-1) of the B-splines
C
C NB    : Number of given points 
C NT    : Maximum number of given points
C
C NBU   : Number of intervals in the u-direction 
C NA1	  : Number of control factors    (NBU+KBSPL)
C NAA	  : Number of control factors ??   (NBU+2*KBSPL+1)
C NA    : Maximum number of total control factors  
C -------------------------------------------------------------------
C
	 PARAMETER (NT=500, NA=100, KBSPL=3)
	 REAL*4 BXYZ(2,NT),BXYZ1(2,NT),UXYZ(NT)
C
C XYZ     : Catersian coordinates of given points
C UXYZ    : U coordinates of given points on the body surface
C
       REAL*4 BJU(NA),B(NA,4),DBJU(NA)
	 REAL*4 XYZCON(NA,2),UTJ(NA),UXYZNE(NA)
C
C UTJ    : U coordinates of nodes and co-located points
C XYZCON : factors for body geometry expansion
C UXYZNE : U coordinates of controlling points 
C
       END MODULE BVAR_mod
C


!
C  **************************************************
C  *												  *
C  *  variables for fenders and cables              *
C  *                                                *
C  *      April. 20, 2005        by   Bin Teng      *
C  **************************************************				
C
       MODULE FCVAR_mod
C
	  INTEGER NOCABLE,NCKIND,NOFender,NFender(50),JCPRO(50)
	  REAL*4  WTL
	  REAL*4  DLS(50),CLNT(50),CXYZ0(3,50),CXYZ1(3,50)
	  REAL*4  DLNT(50),DXYZ0(3,50),DXYZ1(3,50)
C 
C NCABLE: Number of cables;        
C NCKIND: Number of cable kinds;
C NOFender: Number of Fenders;
C NFender(KD): Number of data for tension curve for each Fender

!  WTL: water level
C    CLNT: length of each cable;     DLNT: length of each Fender
C  CXYZ0, CXYZ1: the coordinates of the two ends of each cable
C  DXYZ0, DXYZ1: the coordinates of the two ends of each Fender

        INTEGER NCABLEU,NFenderU,JDRMXY(50),JDRM(50),KindDrm
	  REAL*4 CABX1(50),CABX2(50),CABXYCON(200,2,50),CABUTJ(200,50)
	  REAL*4 DRMX1(50),DRMX2(50),DRMXYCON(200,2,50),DRMUTJ(200,50)
	  REAL*4 CEPSMAX(50), FCRMAX(50)
	  REAL*4 FDRMAX(6,50),EPSMAX(50)
C
C  JDRM    : The code of Fender character
C JDRMXY(N): the code for Fender direction, 
C          =1, force in x-direction; =2, force in y-direction
C KindDrm  : Number of Fender's kinds

       END MODULE FCVAR_mod

