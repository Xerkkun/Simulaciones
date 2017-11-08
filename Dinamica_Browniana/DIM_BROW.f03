PROGRAM DBGAYLOR
  USE INICIAL
  USE FRZAS
  USE FUNC_RAD
!-------------------------------------------------------------------------------
!PROGRAMA MAESTRO DE SIMULACION
!-------------------------------------------------------------------------------
IMPLICIT REAL (A-H,O-Z)

!DECLARACION DE VARIABLES
CHARACTER(len=20)              ::FINAL
INTEGER,DIMENSION(1000)        ::nseed
REAL,DIMENSION(3)              ::R,S
REAL,ALLOCATABLE,DIMENSION(:,:)::cx,cy,cz,cxr,cyr,czr,xr,x,f
REAL,PARAMETER                 ::PI = 3.14159
OPEN(20,FILE="terma.dat",STATUS="REPLACE",ACTION="WRITE")
WRITE(FINAL , '(A,F6.3,A)') "fin",DENS,"dat"
OPEN(30,FILE=FINAL,STATUS="REPLACE",ACTION="WRITE")

!===============================================================================
WRITE(*,*)"CODIGO DE DINAMICA BROWNIANA V0"
!ADQUIRIR PARAMETROS

nCONF  = 1 !1) ALEATORIA, 2) REGULAR
nd     = 3       !2) 2D , 3) 3D
N      = 400     !NUMERO DE PARTICULAS
SIGMA  = 1.0     !DIAMETRO DE LAS PARTICULAS
NSTEP  = 50000  !NUMERO DE CONFIGURACIONES
NENER  = 20000   !CONFIGURACIONES FUERA DE EQUILIBRIO
NFREC  = 100
NFREC2 = 100
DT     = 0.0004

!CALCULOS PRELIMINARES
A      = 1.0/REAL(nd)
!AA Y ZK PARAMETROS DEL POTENCIAL
AA     = 556.0
ZK     = 0.149

PHI    = 4.4E-4
DENS   = 6.0*PHI/PI

!DIMENSIONES DE LA CAJA
BOXL   = ((REAL(N))/DENS)**(A)
RCUT   = BOXL/2.0
VAR    = SQRT(2.0*DT)

!DELIMITAR QUE LAS PARTICULAS ESTEN COMPLETAMENTE DENTRO DE LA CELDA
BOXL1 = BOXL-1.0
NN2 = (NSTEP-NENER)/NFREC
ALLOCATE(cx(N,NN2),cy(N,NN2),cz(N,NN2),xr(N,nd),cxr(N,NN2),cyr(N,NN2),czr(N,NN2),f(N,3))

WRITE(*,*) "SEMILLA ENTERA PARA EL MOV. DE LAS PARTICULAS"
READ*, iseedn

  SELECT CASE(nCONF)

    CASE(1)
      CALL ALEATORIA(N,nd,BOXL1,SIGMA,DENS,x)

    CASE(2)
      CALL REGULAR(N,nd,BOXL1,DENS,x)

  END SELECT


!===============================================================================
ISTEP = 1
!SUBRUTINA PARA LEER LAS FUERZAS SOBRE CADA PARTICULA EN LA CONF INICIAL
CALL FUERZAS(AA,ZK,RCUT,BOXL,N,ISTEP,ENPOT,x,f)
!ESCRIBIR LA ENERGIA DE LA CONFIGURACION IMPORTANTE
!PARA VERIFICAR COMO TERMALIZA EL SISTEMA
WRITE(20,*) ISTEP,ENPOT/REAL(N)
!===============================================================================
!OBTENER NUMEROS ALEATORIOS
nseed = iseedn

!CALL RANDOM_SEED(SIZE = ms)
CALL RANDOM_SEED(PUT = nseed)
!===============================================================================
KI  = 0
KI2 = 0
!CICLO DE CONSTRUCCION DE CONFIGURACIONES
DO ISTEP=1,NSTEP
  !CICLO PARA MOVER A LAS PARTICULAS
  DO i=1,N
    !GENERACION DE 3 NUMEROS ALEATORIOS
    CALL RANDOM_NUMBER(R(:))
    CALL RANDOM_NUMBER(S(:))

    AX=SQRT(-2.0*LOG(R(1)))*COS(2.0*PI*S(1))
    AY=SQRT(-2.0*LOG(R(2)))*COS(2.0*PI*S(2))
    AZ=SQRT(-2.0*LOG(R(3)))*COS(2.0*PI*S(3))

    x(i,1)=x(i,1)+f(i,1)*DT+(VAR*AX)
    x(i,2)=x(i,2)+f(i,2)*DT+(VAR*AY)
    x(i,3)=x(i,3)+f(i,3)*DT+(VAR*AZ)

    xr(i,1)=xr(i,1)+f(i,1)*DT+(VAR*AX)
    xr(i,2)=xr(i,2)+f(i,2)*DT+(VAR*AY)
    xr(i,3)=xr(i,3)+f(i,3)*DT+(VAR*AZ)

    !INCLUYENDO CONDICIONES PERIODICAS
    x(i,1)=x(i,1)-BOXL*ANINT(x(i,1)/BOXL)
    x(i,2)=x(i,2)-BOXL*ANINT(x(i,2)/BOXL)
    x(i,3)=x(i,3)-BOXL*ANINT(x(i,3)/BOXL)
    IF (ISTEP == NSTEP)THEN
      WRITE(30,*) x(i,1), x(i,2), x(i,3)
    END IF
  END DO
!===============================================================================
  !VERIFICANDO SI DEBE ALMACENAR CONFIGURACIONES DE EQUILIBRIO
  !EN LAS MATRICES cx,cy,cz PARA EL CALCULO DE LA G(R)
  IF (MOD(ISTEP,NFREC) == 0.AND.ISTEP > NENER ) THEN
    KI = KI + 1
    DO i=1,N
      cx(i,KI) = x(i,1)
      cy(i,KI) = x(i,2)
      cz(i,KI) = x(i,3)
    END DO

  END IF
    !DECIDIENDO SI ALMACENAMOS LAS CONGIFURACIONES DE EQUILIBRIO
    !EN LAS MATRICES cxr,cyr,czr PARA EL CALCULO DE LAS PROPIEDADES DE
    !AUTODIFUSION W(t) Y D(t)
  IF (MOD(ISTEP,NFREC2) == 0.AND.ISTEP > NENER ) THEN
    KI2 = KI2 + 1
    DO i=1,N
      cxr(i,KI2) = xr(i,1)
      cyr(i,KI2) = xr(i,2)
      czr(i,KI2) = xr(i,3)

    END DO

  END IF
CALL FUERZAS(AA,ZK,RCUT,BOXL,N,ISTEP,ENPOT,x,f)
EN = ENPOT/REAL(N)
WRITE(20,*) ISTEP,EN

END DO
CLOSE(20)
CLOSE(30)
!CALCULO DE LA FUNCION DE DISTRIBUCION RADIAL
CALL GDR(cx,cy,cz,KI,BOXL,DENS,N,GRC,PRESION)


!===============================================================================
  PRINT*, "NUMERO DE ATOMOS:          ", N
  PRINT*, "DIMENSIONES:               ", nd
  PRINT*, "DIAMETRO DE LAS PARTS.     ", SIGMA
  PRINT*, "CONCENTRACION REDUCIDA:    ", DENS
  PRINT*, "LONGITUD DE CELDA:         ", BOXL
  PRINT*, "NUM. TOTAL DE CONFIGS.:    ", NSTEP
!===============================================================================
DEALLOCATE(cx,cy,cz,xr,cxr,cyr,czr,f,x)

END PROGRAM DBGAYLOR

!===============================================================================
