MODULE NODE_CLASS
USE GlobalDataFun
IMPLICIT NONE
!->+++++++++++++++++++++++++++++++++++++++++++++
TYPE,PUBLIC :: NODE

    INTEGER :: ID
    REAL(RDT) :: RLUD,BTA,TA0,BSC

    REAL(RDT) :: EA,EAETA,EAZTA,EIZZ,EIEZ,EIEE,EAC0,EAC1,EAC2,EAC3,&
                        GEA,GZA,GEEA,GZZA,GEZA,GEAEB,GEAZB,GZAEC,GZAZC,& 
                        GJ,GEJ,GZJ,&
            
                        EAD0,EAD1,EAD2,EAD3,EAD4,EAD5,EAD6,EAD7,&
                        EAD0P,EAD1P,EAD2P,EAD3P,EAD4P,EAD6P,EAD7P,&
            
                        EAB0,EAB1,EAB2,EAB3,EAB4,EAB5,EAB6,EAB7,EAB8,&
                        EAB9,EAB10,EAB11,EAB12,EAB13,EAB14,EAB15,&
                        EAB3P,EAB8P

    REAL(RDT) :: MSL,MEM,MZM,IMEE,IMZZ,IMEZ,MD0,MD1,MD2,MD3

    REAL(RDT) :: CTFGRL

CONTAINS

    PROCEDURE,PUBLIC :: CONSTRUCT_NODE

    PROCEDURE,PRIVATE :: INITIAL_NODE

END TYPE NODE
!->+++++++++++++++++++++++++++++++++++++++++++++
CONTAINS
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PUBLIC++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE CONSTRUCT_NODE(THIS,ID,BTA,TAO,RL,STIFF,MASS,CHO)
USE GlobalDataFun
IMPLICIT NONE
CLASS(NODE) :: THIS
INTEGER,INTENT(IN) :: ID
INTEGER,INTENT(IN) :: CHO
REAL(RDT),INTENT(IN) :: STIFF(6,6),MASS(6,6),BTA,TAO,RL

    CALL THIS%INITIAL_NODE()

    THIS%ID=ID
    THIS%RLUD=RL
    THIS%BTA=BTA   
    THIS%TA0=TAO
	
    THIS%MSL=MASS(1,1)
    THIS%MEM=-MASS(1,6)
    THIS%MZM=MASS(1,5)
    THIS%IMEZ=MASS(1,1)*MASS(5,6)
    THIS%IMZZ=MASS(1,1)*MASS(6,6)
    THIS%IMEE=MASS(1,1)*MASS(5,5)

    IF(CHO.EQ.0) THEN!14自由度

        THIS%EA=STIFF(1,1)
        THIS%EAETA=-STIFF(1,4)
        THIS%EAZTA=STIFF(1,3)
        THIS%EIZZ=STIFF(4,4)
        THIS%EIEE=STIFF(3,3)
        THIS%EIEZ=-STIFF(3,4) 
        THIS%EAC0=THIS%EIZZ+THIS%EIEE
        THIS%EAC1=THIS%EAC0*THIS%EAETA/THIS%EA
        THIS%EAC2=THIS%EAC0*THIS%EAZTA/THIS%EA
        THIS%EAC3=THIS%EAC0*THIS%EAC0/THIS%EA
        THIS%GJ=STIFF(2,2)
        THIS%EAB10=THIS%GJ*THIS%EA/THIS%EAC0*THIS%EAETA/THIS%EA      	
        THIS%GEEA=THIS%GJ*THIS%EA/THIS%EAC0
        THIS%GZZA=THIS%GJ*THIS%EA/THIS%EAC0

    ELSE IF(CHO.EQ.1) THEN!23自由度

        THIS%EA=STIFF(1,1)
        THIS%EAETA=-STIFF(1,6)
        THIS%EAZTA=STIFF(1,5)
        THIS%EIZZ=STIFF(6,6)
        THIS%EIEE=STIFF(5,5)
        THIS%EIEZ=-STIFF(5,6) 
        THIS%EAC0=THIS%EIZZ+THIS%EIEE
        THIS%EAC1=THIS%EAC0*THIS%EAETA/THIS%EA
        THIS%EAC2=THIS%EAC0*THIS%EAZTA/THIS%EA
        THIS%EAC3=THIS%EAC0*THIS%EAC0/THIS%EA

        THIS%GEA=STIFF(1,2)   
        THIS%GZA=STIFF(1,3)  
        THIS%GEEA=STIFF(2,2)  
        THIS%GZZA=STIFF(3,3)  
        THIS%GEZA=STIFF(2,3)   
        THIS%GEAEB=-STIFF(2,6)  
        THIS%GEAZB=STIFF(2,5)  
        THIS%GZAEC=-STIFF(3,6) 
        THIS%GZAZC=STIFF(3,5) 
!        THIS%GEJ=(-STIFF(4,6)/STIFF(1,3)-STIFF(4,5)/STIFF(1,2))*STIFF(1,3)
!        THIS%GZJ=(-STIFF(4,6)/STIFF(1,3)-STIFF(4,5)/STIFF(1,2))*STIFF(1,2)
        THIS%GEJ=-STIFF(4,6)
        THIS%GZJ=-STIFF(4,5)
        THIS%GJ=STIFF(4,4)

        THIS%EAB0=STIFF(1,4)
        THIS%EAB1=-STIFF(4,6)
        THIS%EAB2=STIFF(4,5)
        THIS%EAB10=STIFF(3,4)   
        THIS%EAB12=STIFF(2,4)    	

    END IF


END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++PRIVATE++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
!->+++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE INITIAL_NODE(THIS)
IMPLICIT NONE
CLASS(NODE) :: THIS
    
    THIS%RLUD=0.0D0    
    THIS%BTA=0.0D0    
    THIS%TA0=0.0D0    

    THIS%EA=0.0D0     
    THIS%EAETA=0.0D0  
    THIS%EAZTA=0.0D0  
    THIS%EIZZ=0.0D0   
    THIS%EIEZ=0.0D0   
    THIS%EIEE=0.0D0   
    THIS%EAB0=0.0D0   
    THIS%EAB1=0.0D0   
    THIS%EAB2=0.0D0   
    THIS%EAB3=0.0D0   
    THIS%EAB4=0.0D0   
    THIS%EAB5=0.0D0   
    THIS%EAB6=0.0D0   
    THIS%EAB7=0.0D0   
    THIS%EAB8=0.0D0   
    THIS%EAB9=0.0D0   
    THIS%EAB10=0.0D0  
    THIS%EAB11=0.0D0  
    THIS%EAB12=0.0D0  
    THIS%EAB13=0.0D0  
    THIS%EAB14=0.0D0  
    THIS%EAB15=0.0D0  
    THIS%EAB3P=0.0D0  
    THIS%EAB8P=0.0D0  
    THIS%EAC1=0.0D0   
    THIS%EAC2=0.0D0   
    THIS%EAC3=0.0D0   
    THIS%EAD0=0.0D0   
    THIS%EAD1=0.0D0   
    THIS%EAD2=0.0D0   
    THIS%EAD3=0.0D0   
    THIS%EAD4=0.0D0   
    THIS%EAD5=0.0D0   
    THIS%EAD6=0.0D0   
    THIS%EAD7=0.0D0   
    THIS%EAD0P=0.0D0  
    THIS%EAD1P=0.0D0  
    THIS%EAD2P=0.0D0  
    THIS%EAD3P=0.0D0  
    THIS%EAD4P=0.0D0  
    THIS%EAD6P=0.0D0  
    THIS%EAD7P=0.0D0  
    THIS%GJ=0.0D0     
    THIS%GEJ=0.0D0    
    THIS%GZJ=0.0D0    
    THIS%GEA=0.0D0   
    THIS%GZA=0.0D0   
    THIS%GEEA=0.0D0   
    THIS%GZZA=0.0D0   
    THIS%GEZA=0.0D0   
    THIS%GEAEB=0.0D0  
    THIS%GEAZB=0.0D0  
    THIS%GZAEC=0.0D0  
    THIS%GZAZC=0.0D0  
    THIS%MSL=0.0D0    
    THIS%MEM=0.0D0    
    THIS%MZM=0.0D0    
    THIS%MD0=0.0D0    
    THIS%MD1=0.0D0    
    THIS%MD2=0.0D0    
    THIS%MD3=0.0D0    
    THIS%IMEE=0.0D0   
    THIS%IMZZ=0.0D0	 
    THIS%IMEZ=0.0D0   

    THIS%CTFGRL=0.0D0

END SUBROUTINE
!->+++++++++++++++++++++++++++++++++++++++++++++
END MODULE
