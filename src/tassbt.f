!  TASSBT+ TINBODT+ NSBDT
!
C *********************************************************
C *                                                       *
C * Calculate the element contribution and assembly the   *
C * coefficients of the corresponding system of equationn *
C *                                                       *
C *********************************************************
C
        SUBROUTINE TASSBT
C  
	  USE MVAR_MOD
	  USE PVAR_MOD
	  USE MFUNC_mod

        IMPLICIT NONE  
	  INTEGER  INODE,IND,IS,IP,J
	  REAL*4 XP,YP,ZP,DPOX,DPOY,DPOZ,DPDN
	  REAL*4 RSN(4,4),EX(4),EY(4)
        REAL*4 BMAT(NSYS),CMAT(NNODED,NSYS)
C
        DATA RSN /1.,  1.,  1.,  1., 
     1            1., -1.,  1., -1.,
     2            1.,  1., -1., -1., 
     3            1., -1., -1.,  1./ 
	  DATA EX / 1.,  1., -1., -1./                                            
        DATA EY / 1., -1., -1.,  1./
C
C ***************************************
C
	  CMAT(:,:)=0.0

       DO  20   INODE=1,  NNF  
	  DO IS=1, NSYS
	   DO IP=1, NSYS
	    CMAT(INODE,IS)=CMAT(INODE,IS)+RSN(IS,IP)*BKN(INODE,IP)
	   ENDDO
	  ENDDO
20	 CONTINUE
!
! ==============================================
!
	 DO 40 INODE=NNF+1, NNODED
!	   
	   DO 40 IP=1, NSYS 
	    XP=EX(IP)*XYZ(1,INODE)
	    YP=EY(IP)*XYZ(2,INODE)
	    ZP=       XYZ(3,INODE)
	    CALL DINP(XP,YP,ZP,DPOX,DPOY,DPOZ)
	    DPDN=(DPOX*EX(IP)*DXYZ(1,INODE)+
	1		  DPOY*EY(IP)*DXYZ(2,INODE)+
	2		  DPOZ       *DXYZ(3,INODE) )
!		
		!DPDN=DPDN-DSDT(1)*EX(IP)*DXYZ(1,INODE)
		!DPDN=DPDN-DSDT(2)*EY(IP)*DXYZ(2,INODE)
		!DPDN=DPDN-DSDT(3)*       DXYZ(3,INODE)
		!DPDN=DPDN-DSDT(4)*EY(IP)*DXYZ(4,INODE)
		!DPDN=DPDN-DSDT(5)*EX(IP)*DXYZ(5,INODE)
		!DPDN=DPDN-DSDT(6)*EX(IP)*EY(IP)*DXYZ(6,INODE)
        dpdn=-dpdn
!
	  DO  IS=1, NSYS 
	   CMAT(INODE,IS)=CMAT(INODE,IS)+RSN(IS,IP)*DPDN
	  ENDDO
C
40      Continue   
!
! ==============================================
!
	  BMATA(:,:)=0.0d0
	  DO 200 IS=1, NSYS
	  DO 200 IND=1,   NNODE
	  DO 200 INODE=1, NNODED
	   BMATA(IND,IS)=BMATA(IND,IS)+CMATA(IND,INODE,IS)*CMAT(INODE,IS)
200     continue   
! 
!        WRITE(102, 610)  TimeRK                 
!       DO 600   INODE=1,  NNODE
!        XP=XYZ(1,INODE)
!        YP=XYZ(2,INODE)
!        ZP=XYZ(3,INODE) 
!        WRITE(102, 620)  INODE,XP,YP,ZP,BMATA(INODE,1)
! 600   CONTINUE 
!
 610   FORMAT(10x,'NNODE     BMATA(IND,1)     T=',F14.6)
 620   FORMAT(1X,I4,4(2X,F13.6))                            
C      
          DO 300 IS=1, NSYS   
            CALL RLUBKSB(IS,AMATA,NNODE,NNODE, 1,NSYS, 1,INDX,BMATA)
300       CONTINUE
C                 
C ** output the results
C      
        DO 500 IP=1, NSYS 
        DO 360 IND=1, NNODE 
        BMAT(IP)=(0.0D0, 0.0D0)
        DO 350 IS=1, NSYS
350     BMAT(IP)=BMAT(IP)+BMATA(IND,IS)*RSN(IP,IS)
        BMAT(IP)=BMAT(IP)/NSYS
        UNKN(IND,IP)=BMAT(IP)
360     CONTINUE
C
500     CONTINUE
C                 
       RETURN
       END
