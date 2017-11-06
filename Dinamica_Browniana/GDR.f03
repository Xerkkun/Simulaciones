MODULE FUNC_RAD
CONTAINS
SUBROUTINE GDR(cx,cy,cz,NTMAX,BOXL,DENS,N,GRC,PRESION)
  IMPLICIT REAL(A-H,O-Z)

  CHARACTER(len=20)   :: ARCHIVO
  INTEGER,ALLOCATABLE :: NHIST(:)
  INTEGER             :: i,j,k,NTMAX,MAXBIN,N
  REAL,PARAMETER      :: PI = 3.14159
  REAL,ALLOCATABLE    :: cx(:,:),cy(:,:),cz(:,:)

  DELTAR = 0.01
  MAXBIN = INT(BOXL/(2*DELTAR))

  ALLOCATE(NHIST(MAXBIN+50))
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

          zi0 = cz(i,k)
          zit = cz(j,k)
          zi0t = zi0 - zit

          !CONDICION DE IMAGEN MINIMA
          xi0t = xi0t-BOXL*ANINT(xi0t/BOXL)
          yi0t = yi0t-BOXL*ANINT(yi0t/BOXL)
          zi0t = zi0t-BOXL*ANINT(zi0t/BOXL)

          r0t = SQRT(xi0t**2 + yi0t**2 + zi0t**2)
          NBIN = INT(r0t/DELTAR)+1
          print*, NBIN
          IF(NBIN.LE.MAXBIN)THEN
            NHIST(NBIN) = NHIST(NBIN)+1
          END IF

        END DO
      END IF
    END DO
  END DO
  C1 = PI*DENS

  WRITE(ARCHIVO , '(A,F6.3,A)') "gr",DENS,".dat"
  OPEN(50,FILE=ARCHIVO,STATUS="REPLACE",ACTION="WRITE")

  DO NBIN=1,MAXBIN
    RL = REAL(NBIN-1)*DELTAR
    RU = RL + DELTAR
    RT = RL + DELTAR/2
    C2 = C1*(RU**2-RL**2)
    GDRTA = REAL(NHIST(NBIN))/(REAL(NTMAX)*REAL(N)*C2)

    WRITE(50,*)RT,GDRTA

    IF ((RT-1.0 <= 0.005.AND.GDRTA /= 0.0))THEN
      PRESION = 1.0+(0.5*PI*DENS*GDRTA)
      GRC = GDRTA
      PRINT*, "G DE CONTACTO:             ", GRC
      PRINT*, "PRESION REDUCIDA:          ", PRESION
    END IF
  END DO

CLOSE(50)
DEALLOCATE(NHIST)
END SUBROUTINE GDR
END MODULE FUNC_RAD
