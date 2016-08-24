!  
!  G=-1/4/pi*(1/r+1/r1) 
!
!
        SUBROUTINE TGRN(H,X,X0,Y,Y0,Z,Z0,XHF) 
! 
        IMPLICIT NONE
		INTEGER I

	    REAL*4,INTENT(IN)::  H,X,Y,Z,X0,Y0,Z0
		REAL*4,INTENT(OUT):: XHF(4)
	    REAL*4 DX,DY,DZ,DZ1
		REAL*4 RXY2,DZ02,DZ12,R02,R12,SR,SR1,PI,PI4
!
	    DATA PI/3.14159265358979/ 
!C 
!C H: water depth, negtive for infinity water depth 
!C
        PI4 = 4.0D0*PI

        DX=X-X0
        DY=Y-Y0
        DZ  =Z- Z0
        DZ1 =Z+ Z0 +2.0*H


		RXY2=DX*DX+DY*DY
	    DZ02=DZ*DZ
		DZ12=DZ1*DZ1
!
	   IF(H .GT. 0) THEN
        R02=RXY2+DZ02                      
        R12=RXY2+DZ12           
        SR =SQRT(R02) 
        SR1=SQRT(R12)   
! 
        XHF(1)= 1.D0/SR +  1.D0/SR1 
        XHF(2)=-DX/SR**3 -  DX /SR1**3 
        XHF(3)=-DY/SR**3 -  DY /SR1**3 
        XHF(4)=-DZ/SR**3 -  DZ1/SR1**3    
! 
	   ELSE
        R02=RXY2+DZ02                      
        SR =SQRT(R02) 
! 
        XHF(1)= 1.D0/SR  
        XHF(2)=-DX/SR**3  
        XHF(3)=-DY/SR**3  
        XHF(4)=-DZ/SR**3     
! 
	   ENDIF
!
       DO 120   I=1,  4
         XHF(I)=-XHF(I)/PI4 
120	   CONTINUE
!
        RETURN 
        END 
		           

