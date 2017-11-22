PROGRAM SIMULACION
  USE INICIAL
  USE ENERGIA
  USE FUNC_RAD
!-------------------------------------------------------------------------------
!PROGRAMA MAESTRO DE SIMULACION
!-------------------------------------------------------------------------------
IMPLICIT REAL*4 (A-H,O-Z)

!DECLARACION DE VARIABLES
CHARACTER(len=20)       ::ARCHIVO,FINAL,VIS
INTEGER                 ::nCONF,N,nd,NSTEP,KI2,l,ms,ISTEP,iseedn,NENER,NN2,i
INTEGER,DIMENSION(1000) ::nseed
REAL*4,ALLOCATABLE      ::xseed(:)
REAL*4,ALLOCATABLE      ::x(:,:),rnum(:,:),cx(:,:),cy(:,:),cz(:,:)
REAL,PARAMETER          ::PI = 3.14159
!===============================================================================
WRITE(*,*)"CODIGO DE MONTECARLO V0"
!ADQUIRIR PARAMETROS
!WRITE(*,*) "SELECCIONE EL TIPO DE CONFIGURACION QUE DESEA:"
!WRITE(*,*) "1)ALEATORIA"
!WRITE(*,*) "2)REGULAR"
!READ*, nCONF
nCONF = 1
!WRITE(*,*) "INDIQUE EL NUMERO DE DIMENSIONES DESEADA:"
!WRITE(*,*) "2) BIDIMENSIONAL"
!WRITE(*,*) "3) TRIDIMENSIONAL"
!READ*,nd
nd = 3

!WRITE(*,*)"NUMERO DE PARTICULAS:"
!READ*,N
N = 400

!WRITE(*,*)"CONCENTRACION REDUCIDA:"
!READ*,DENS
PHI    = 4.4E-4
DENS   = 6.0*PHI/PI

!WRITE(*,*)"RAZON DE CRITERIO DE ACEPTACION"
!READ*,ACEPT
ACEPT = 0.5

!WRITE(*,*) "DIAMETRO DE LAS PARTICULAS:"
!READ*,SIGMA
SIGMA = 1.0

!PRINT*, "NUMERO DE CONFIGURACIONES:"
!READ*, NSTEP
NSTEP = 100000
!PRINT*, "NUMERO DE CONFIGURACIONES FUERA DE EQUILIBRIO"
!READ*, NENER
NENER = 15000
WRITE(*,*) "SEMILLA ENTERA PARA EL MOV. DE LAS PARTICULAS Y PARA EL ALGORITMO DE MC:"
READ*, iseedn

!PRINT*, "FRECUENCIA DE IMPRESION EN PANTALLA:"
!READ*, IPRINT
IPRINT = 100
!PRINT*, "FRECUENCIA DE GUARDADO DE CONFIGURACIONES:"
!READ*, ISAVE
ISAVE = 100
!PRINT*, "FRECUENCIA PARA CORREGIR DESPLAZAMIENTO MAXIMO:"
!READ*, IRATIO
IRATIO = 100
!OPEN(60,FILE="PvsDENS.dat",STATUS="REPLACE",ACTION="WRITE")
!===============================================================================
!CALCULOS PRELIMINARES
A=1.0/(nd*1.0)
!DO i=0,18
!  DENS = 0.05 + 0.05*i
!  PRINT*, "DENSIDAD:",DENS
!  NSTEP = 12000
!  NENER = 2000

!  IF (DENS <= 0.2)THEN
!    nCONF = 1
!    ACEPT = 0.9

!  ELSE IF (0.2 < DENS.AND.DENS <= 0.6)THEN
!    nCONF = 1
!    ACEPT = 0.5

!  ELSE IF (0.6 < DENS.AND.DENS <= 0.85)THEN
!    nCONF = 2
!    ACEPT = 0.3

!  ELSE IF (DENS > 0.85)THEN
!    N     = 400
!    nCONF = 2
!    ACEPT = 0.2
!    NSTEP = 12000
!    NENER = 2000
!  END IF

  !DIMENSIONES DE LA CAJA
  BOXL = ((REAL(N))/DENS)**(A)
  !DESPLAZAMIENTO MAXMO
  DRMAX = 0.1
  KI2 = 0
  ACM = 0
  ACATMA = 0.0
  !RADIO DE CORTE
  RCUT = SIGMA
  !DELIMITAR QUELAS PARTICULAS ESTEN COMPLETAMENTE DENTRO DE LA CELDA
  BOXL1 = BOXL-1.0
  NN2 = (NSTEP-NENER)/ISAVE

  SELECT CASE(nCONF)

    CASE(1)
      CALL ALEATORIA(N,nd,BOXL1,SIGMA,DENS,x)

    CASE(2)
      CALL REGULAR(N,nd,BOXL1,DENS,x)

    END SELECT
!===============================================================================
  !CORRECCION DE LARGO ALCANCE
  VLRC = 0.0
!===============================================================================
  !SUBRUTINA PARA LEER LA ENERGIA DE LA CONFIGURACION INICIAL
  CALL ENER_CONF(N,x,BOXL,V)
  VI = V + VLRC
!===============================================================================
  !OBTENER NUMEROS ALEATORIOS
  nseed = iseedn

  CALL RANDOM_SEED(SIZE = ms)
  CALL RANDOM_SEED(PUT = nseed)

!===============================================================================
  WRITE(ARCHIVO , '(F6.3,A)') DENS,".dat"
  OPEN(20,FILE=ARCHIVO,STATUS="REPLACE",ACTION="WRITE")
  WRITE(FINAL, '(A,F6.3,A)') "FIN", DENS, ".dat"
  OPEN(70,FILE=FINAL,STATUS="REPLACE",ACTION="WRITE")

  ALLOCATE(rnum(N,nd),xseed(N),cx(N,NN2),cy(N,NN2),cz(N,NN2)) !lyr

    DO 100 ISTEP=1,NSTEP

      DO 90 l=1,N
        xo=x(l,1) !lyr
        yo=x(l,2) !lyr
        zo=x(l,3)
        !ENERGIA DE LA L-ESIMA PARTICULA EN LA CONFIGURACION VIEJA
        CALL ENER_PART(N,l,xo,yo,zo,x,BOXL,Vold) !lyr

        !MOVER A LAS PARTICULAS (MOVIMIENTO TENTATIVO)
        CALL RANDOM_NUMBER(rnum(l,:))
        xn=xo+((2.0*rnum(l,1)-1.0)*DRMAX) !lyr
        yn=yo+((2.0*rnum(l,2)-1.0)*DRMAX) !lyr
        zn=zo+((2.0*rnum(l,3)-1.0)*DRMAX)

        !INCLUYENDO CONDICIONES PERIODICAS
        xn=xn-BOXL*ANINT(xn/BOXL) !lyr
        yn=yn-BOXL*ANINT(yn/BOXL) !lyr
        zn=zn-BOXL*ANINT(zn/BOXL)

        !ENERGIA DE LA L-ESIMA PARTICULA EN LA CONFIGURACION NUEVA
        CALL ENER_PART(N,l,xn,yn,zn,x,BOXL,Vnew) !lyr

        !CRITERIOS DE ACEPTACION O RECHAZO
        DELTAV = Vnew - Vold

        IF (DELTAV < 75.0)THEN
          CALL RANDOM_NUMBER(rmc)
          IF (DELTAV.LE.0.0.OR.EXP(-DELTAV).GT.rmc)THEN
            V = V + DELTAV
            x(l,1) = xn
            x(l,2) = yn
            x(l,3) = zn

            ACATMA = ACATMA + 1.0
          END IF

        END IF
        ACM = ACM + 1.0
        90 END DO

        !GUARDANDO ENERGIA/PARTICULA DE LA CONFIGURACION PARA VERIFICAR TERMALIZACION
        VN = (V + VLRC)/ REAL(N) !lyr

        !VERIFICANDO SI AJUSTA EL DESPLAZAMIENTO MAXIMO
        IF (MOD(ISTEP,IRATIO) == 0)THEN
          RATIO = ACATMA / REAL(N*IRATIO)

          IF (RATIO > ACEPT) THEN !lyr
            DRMAX = DRMAX*1.05
          ELSE
            DRMAX = DRMAX*0.95
          END IF
          ACATMA = 0.0
        END IF

        !VERIFICANDO SI REQUIERE ESCRIBIR INFORMACION DE EJECUCION
        IF (MOD (ISTEP, IPRINT) == 0) THEN
          WRITE(*,*) ISTEP, RATIO, DRMAX, VN !lyr
        END IF

!===============================================================================
!        IF (MOD(ISTEP,ISAVE) == 0)THEN
!          IF (ISTEP < 10)THEN
!            WRITE(VIS,'(I1,A)')ISTEP,".dat"
!          ELSE IF (ISTEP < 100)THEN
!            WRITE(VIS,'(I2,A)')ISTEP,".dat"
!          ELSE IF (ISTEP < 1000)THEN
!            WRITE(VIS,'(I3,A)')ISTEP,".dat"
!          ELSE
!            WRITE(VIS,'(I4,A)')ISTEP,".dat"
!          END IF
!        END IF
!        OPEN(15,FILE=VIS,STATUS="REPLACE",ACTION="WRITE")

!        DO kk=1,N
!          write(15,*) x(kk,1),x(kk,2),x(kk,3)
!        END DO
!===============================================================================

        !VERIFICANDO SI DEBE ALMACENAR CONFIGURACIONES DE EQUILIBRIO
        IF (MOD(ISTEP,ISAVE) == 0.AND.ISTEP > NENER ) THEN
          WRITE(20,*) ISTEP,VN
          KI2 = KI2 + 1

          !CONSTRUIR ARREGLOS PARA LAS CONFIGURACIONES SELECCIONADAS
          DO k=1,N
            cx(k,KI2)=x(k,1)
            cy(k,KI2)=x(k,2)
            cz(k,KI2)=x(k,3)
          END DO

        END IF
        CLOSE(15)
  100 END DO

  !GUARDAR CONFIGURACION FINAL
  DO i=1,N
      WRITE(70,*) x(i,1), x(i,2), x(i,3)
  END DO
      !CALCULO DE LA FUNCION DE DISTRIBUCION RADIAL
      CALL GDR(cx,cy,cz,KI2,BOXL,DENS,N,GRC,PRESION)
      WRITE(60,*)GRC,PRESION,DENS


!===============================================================================
  PRINT*, "NUMERO DE ATOMOS:          ", N
  PRINT*, "DIMENSIONES:               ", nd
  PRINT*, "DIAMETRO DE LAS PARTS.     ", SIGMA
  PRINT*, "CONCENTRACION REDUCIDA:    ", DENS
  PRINT*, "LONGITUD DE CELDA:         ", BOXL
  PRINT*, "ENERGIA DE LA CONF. INICIAL", VI
  PRINT*, "NUM. TOTAL DE CONFIGS.:    ", NSTEP
!===============================================================================
  DEALLOCATE(x,rnum,cx,cy,cz,xseed)
  CLOSE(20)
  CLOSE(70)
!END DO
!  CLOSE(60)
END PROGRAM SIMULACION
