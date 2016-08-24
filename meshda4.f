!       MESHFS4 + MESHBD +MESHT
!
C *******************************************************************
C *                                                                 *
C *  Read in the data file of the free surface mesh                 *
C *                                                                 *
C *                                                                 *
C *******************************************************************
C 
        SUBROUTINE MESHFS4

	  USE MVAR_MOD
        IMPLICIT   NONE  

	  INTEGER IE,J,M
C
	  XYZE(3,1:8, 1:NELEMF)=0.0d0
!
        DO 100 IE=1, NELEMF
	  IETYPE(IE)=2
	  READ(3, *)    M, NCN(IE)
        READ(3, *) (XYZE(1,J,IE), J=1, NCN(IE))
        READ(3, *) (XYZE(2,J,IE), J=1, NCN(IE))
        READ(3, *) (DAMPE(J,IE), J=1, NCN(IE))
C
100    CONTINUE
!
!
!  数据文件中物面法方向指入流体为正
!
	  DXYZE(1, 1:8, 1:NELEMF)= 0.0d0
	  DXYZE(2, 1:8, 1:NELEMF)= 0.0d0
	  DXYZE(3, 1:8, 1:NELEMF)=-1.0d0 
C

      RETURN
      END



C ***************************************************************
C *                                                             *
C *  Generate nodal and element data on the body surface of     *
C *  an arbitrary body                                          *
C *                                                             *
C ***************************************************************
C 
        SUBROUTINE MESHBD(NCOR)

	  USE MVAR_MOD
        IMPLICIT   NONE  

!	  INTEGER NELEM,NELEMB,NELEMF,NNODE,NNODED,NNB,NNBD,NNF,NNTCH


	  INTEGER NCOR,IPL,IPOLAR(50)
	  INTEGER I,IE,M,INODE,NCNN,K,KK
	  REAL*8  A1,A2,R1,R2,Z1
        REAL*8  XOFSET(50),YOFSET(50),ZOFSET(50)
!
! --------------------------------------------
!
        DO 5 I=1, NCOR
        READ(2,*)  M, IPOLAR(I), XOFSET(I),YOFSET(I),ZOFSET(I)
5       CONTINUE
!
        DO 10 INODE=1, NNB
        READ(2,*) M, IPL, R1, A1, Z1
        IF(IPOLAR(IPL).EQ.0) THEN
         XYZB(1,INODE)=R1+XOFSET(IPL)
         XYZB(2,INODE)=A1+YOFSET(IPL)
        ELSE IF(IPOLAR(IPL).EQ.1) THEN
         XYZB(1,INODE)=R1*DCOS(A1*PI/180.0D0)+XOFSET(IPL)
         XYZB(2,INODE)=R1*DSIN(A1*PI/180.0D0)+YOFSET(IPL)
        ENDIF
	   XYZB(3,INODE)=Z1+ZOFSET(IPL)
10      CONTINUE
!       
        DO 11 INODE=1, NNBD
        READ(2,*) M,IPL,R2,A2,DXYZB(3,INODE)
        IF(IPOLAR(IPL).EQ.0) THEN
          DXYZB(1,INODE)= R2
          DXYZB(2,INODE)= A2
        ELSE IF(IPOLAR(IPL).NE.0) THEN
          DXYZB(1,INODE)=R2*DCOS( A2*PI/180.0D0 )
          DXYZB(2,INODE)=R2*DSIN( A2*PI/180.0D0 )
        ENDIF
11      CONTINUE
C
C      
        DO 20 I=1, NELEMB
	  IE=I+NELEMF
	  IETYPE(IE)=1
        READ(2,*) M, NCN(IE)    
        READ(2,*) (NCONB(I,K), K=1, NCN(IE))
20    CONTINUE
!
	 DO 30 I=1,  NELEMB
	  IE=I+NELEMF
        READ(2,*)  M, NCNN
        READ(2,*) (NCONDB(I,K), K=1, NCN(IE))
30     CONTINUE
!
! ===============================================================
!    转换网格数据格式，将水线上节点从物面上消除，归入到水面上
!
       DO 100 I=1, NELEMB
	  IE=I+NELEMF
	  DO 80 K=1, NCN(IE)
         XYZE(1,K,IE)=XYZB(1, NCONB(I,K))
         XYZE(2,K,IE)=XYZB(2, NCONB(I,K))
         XYZE(3,K,IE)=XYZB(3, NCONB(I,K))
80      CONTINUE
!
100    CONTINUE
!
       DO 200 I=1, NELEMB
	  IE=I+NELEMF
	  DO 180 K=1, NCN(IE)
         DXYZE(1,K,IE)=DXYZB(1,NCONDB(I,K))
         DXYZE(2,K,IE)=DXYZB(2,NCONDB(I,K))
         DXYZE(3,K,IE)=DXYZB(3,NCONDB(I,K))
180     CONTINUE
!
200    CONTINUE
!
       RETURN
       END


C  
C *******************************************************
C *                                                     *
C *  Compute mesh position at T-time                    *
C *                                                     *
C *******************************************************
C 
        SUBROUTINE MESHT

	  USE MVAR_MOD
	  USE PVar_mod
        IMPLICIT   NONE  

	  INTEGER IE,J
C
C    Rotation 
C
        DO 100 IE=1, NELEM
	   DO 80 J=1, NCN(IE)
          TXYZE(1,J,IE)=(TXYZE(1,J,IE)-XC)*TRXYZ(1,1)+
	1			      (TXYZE(2,J,IE)-YC)*TRXYZ(1,2)+
     2                  (TXYZE(3,J,IE)-ZC)*TRXYZ(1,3)   
          TXYZE(2,J,IE)=(TXYZE(1,J,IE)-XC)*TRXYZ(2,1)+
	1	              (TXYZE(2,J,IE)-YC)*TRXYZ(2,2)+
     2                  (TXYZE(3,J,IE)-ZC)*TRXYZ(2,3)   
          TXYZE(3,J,IE)=(TXYZE(1,J,IE)-XC)*TRXYZ(3,1)+
	1	              (TXYZE(2,J,IE)-YC)*TRXYZ(3,2)+
     2                  (TXYZE(3,J,IE)-ZC)*TRXYZ(3,3)        
80      CONTINUE
C
100    CONTINUE
C
C  Translation
C
        DO 200 IE=1, NELEM

	  DO 180 J=1, NCN(IE)
         TXYZE(1,J,IE)=TXYZE(1,J,IE)+XC+DISP(1)
         TXYZE(2,J,IE)=TXYZE(2,J,IE)+YC+DISP(2)  
         TXYZE(3,J,IE)=TXYZE(3,J,IE)+ZC+DISP(3)  
180      CONTINUE
C      
200    CONTINUE
C
      RETURN
      END


!C         
!C ******************************************************
!C *  Attention: the order of rotations                 *
!C *                                                    *
!C ******************************************************
!C
        SUBROUTINE TRMAX
        USE MVar_mod
	  USE PVar_mod
	  IMPLICIT NONE

	  REAL*4 CSAX,SNAX,CSAY,SNAY,CSAZ,SNAZ
        REAL*4 TRX(3,3),TRY(3,3),TRZ(3,3),TRXY(3,3)
!C
        CSAX=COS(DISP(4))
        SNAX=SIN(DISP(4))
        CSAY=COS(DISP(5))
        SNAY=SIN(DISP(5))
        CSAZ=COS(DISP(6))
        SNAZ=SIN(DISP(6))

!                           
        TRXYZ(1,1)= CSAY*CSAZ
        TRXYZ(2,1)= CSAX*SNAZ+SNAX*SNAY*CSAZ
        TRXYZ(3,1)= SNAX*SNAZ-CSAX*SNAY*CSAZ

        TRXYZ(1,2)=-CSAY*SNAZ
        TRXYZ(2,2)= CSAX*CSAZ-SNAX*SNAY*SNAZ
        TRXYZ(3,2)= SNAX*CSAZ+CSAX*SNAY*SNAZ

        TRXYZ(1,3)= SNAY
        TRXYZ(2,3)=-SNAX*CSAY
        TRXYZ(3,3)= CSAX*CSAY
!
        RETURN
        END

