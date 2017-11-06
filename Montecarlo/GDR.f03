MODULE FUNC_RAD
CONTAINS
SUBROUTINE GDR(cx,cy,NTMAX,BOXL,DENS,N,GRC,PRESION)
  IMPLICIT REAL*4(A-H,O-Z)
!===============================================================================
!DECLARACION DE VARIABLES
  CHARACTER(len=20)   :: ARCHIVO
  INTEGER,ALLOCATABLE :: NHIST(:)
  INTEGER             :: i,j,k,NTMAX,MAXBIN,N
  REAL,PARAMETER      :: PI = 3.14159
  REAL,ALLOCATABLE    :: cx(:,:),cy(:,:)
!===============================================================================
!CALCULOS PRELIMINARES
  DELTAR = 0.01
  MAXBIN = INT(BOXL/(2*DELTAR))
  ALLOCATE(NHIST(MAXBIN))
!===============================================================================

  NHIST(:) = 0
  DO i=1,N
    DO j=1,N
      IF(i.NE.j)THEN
        DO k=1,NTMAX
          xi0 = cx(i,k)
          xit = cx(j,k)
          xi0t = xi0 - xit

          yi0 = cy(i,k)
          yit = cy(j,k)
          yi0t = yi0 - yit

          !CONDICION DE IMAGEN MINIMA
          xi0t = xi0t-BOXL*ANINT(xi0t/BOXL)
          yi0t = yi0t-BOXL*ANINT(yi0t/BOXL)

          !CALCULO DE LA DISTANCIA ENTRE LA I-ESIMA Y LA J-ESIMA
          r0t = SQRT(xi0t**2 + yi0t**2)

          !CUANTAS VECES CABE DELTAR EN R0T
          NBIN = INT(r0t/DELTAR)+1

          !HISTOGRAMA
          IF(NBIN.LE.MAXBIN)THEN
            NHIST(NBIN) = NHIST(NBIN)+1
          END IF

        END DO
      END IF
    END DO
  END DO
!===============================================================================
  WRITE(ARCHIVO , '(A,F6.3,A)') "gr",DENS,".dat"
  OPEN(50,FILE=ARCHIVO,STATUS="REPLACE",ACTION="WRITE")
!===============================================================================
!CALCULO DE LA G(R)
C1 = PI*DENS
  DO NBIN=1,MAXBIN
    RL = REAL(NBIN-1)*DELTAR
    RU = RL + DELTAR
    RT = RL + DELTAR/2
    C2 = C1*(RU**2-RL**2)
    GDRTA = REAL(NHIST(NBIN))/(REAL(NTMAX)*REAL(N)*C2)
    WRITE(50,*)RT,GDRTA
!===============================================================================
!CALCULO DE LA G DE CONTACTO
    IF ((RT-1.0 <= 0.005.AND.GDRTA /= 0.0))THEN
      PRESION = 1.0+(0.5*PI*DENS*GDRTA)
      GRC = GDRTA
      PRINT*, "G DE CONTACTO:             ", GRC
      PRINT*, "PRESION REDUCIDA:          ", PRESION
    END IF
  END DO
!===============================================================================
CLOSE(50)
DEALLOCATE(NHIST)
END SUBROUTINE GDR
END MODULE FUNC_RAD
