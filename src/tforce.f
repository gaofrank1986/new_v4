C
C ******************************************************
C *                                                    *
C *  Evaluate the wave force on a 3-D body             *
C *                                                    *
C ******************************************************
C
        SUBROUTINE TFORCE
	  USE MVAR_MOD
	  USE PVAR_MOD
	  USE MFUNC_mod
!
        IMPLICIT   NONE
!
	  INTEGER IP,IELEM,NSAMB,N,K
        REAL*4  EXY(4,2),XP,YP,ZP    
        REAL*4  DDUM
C
        DATA EXY /1., 1.,-1.,-1.,  1.,-1.,-1., 1./
C              
        FORCEW=0.0D0
C
        DO 100  IELEM=NELEMF+1,  NELEM 
        NSAMB=16
        IF(NCN(IELEM) .EQ.6 ) NSAMB=4
C                    
        DO  100  IP=1,  NSYS   

        DO  80  N =1, NSAMB  
C
        DDUM=0.0D0
        DO 40   K=1,  NCN(IELEM) 
          DDUM=DDUM+ DPDT(NCON(IELEM,K),IP)*SAMB(IELEM,N,K) 
40      CONTINUE 
C

        XP=EXY(IP,1)*SAMBXY(IELEM,N,1)
        YP=EXY(IP,2)*SAMBXY(IELEM,N,2)  
        ZP=          SAMBXY(IELEM,N,3)      
C 
        DDUM=DDUM+DPOT(XP,YP,ZP)*SAMB(IELEM,N,0) 
C                
        FORCEW(1)=FORCEW(1)+DDUM*
     1                        EXY(IP,1)*DSAMB(IELEM,N,1)
        FORCEW(2)=FORCEW(2)+ DDUM*
     1                        EXY(IP,2)*DSAMB(IELEM,N,2)

        FORCEW(3)=FORCEW(3) + DDUM* DSAMB(IELEM,N,3)       
C
        FORCEW(4)=FORCEW(4)+DDUM*
     2               EXY(IP,2)* DSAMB(IELEM,N,4) 
        FORCEW(5)=FORCEW(5)+DDUM*
     2               EXY(IP,1)* DSAMB(IELEM,N,5)
        FORCEW(6)=FORCEW(6)+DDUM*
     2               EXY(IP,2)* EXY(IP,1)*DSAMB(IELEM,N,6)

!	  Print *,'  DDUM=',DDUM,' DX=',DSAMB(IELEM,N,1),' EX=',EXY(IP,1)

80      CONTINUE
!
!	  Print *,' IELEM=',IELEM,' F1=',FORCEW(1),' F3=',FORCEW(3)

100     CONTINUE
!
C
	  FORCEW(:)=RHO*FORCEW(:)
	  
!	  Print *,  ' F1=',FORCEW(1),' F3=',FORCEW(3)
	  
!	  pause

        RETURN
        END                      


