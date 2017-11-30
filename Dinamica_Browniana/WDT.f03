MODULE W
CONTAINS
  SUBROUTINE WDT(cxr,cyr,czr,KI2,DT,NFREC2,N)
    IMPLICIT REAL (A-H,O-Z)
    REAL,ALLOCATABLE,DIMENSION(:,:)::cxr,cyr,czr
    OPEN(40,FILE="wdt.dat",STATUS="REPLACE",ACTION="WRITE")
    !TIEMPO ENTRE CONFIGURACIONES ALMACENADAS
    TIM=REAL(NFREC2)*DT
    !CICLO PARA DAR LA CADENCIA EN EL BARRIMIENTO TEMPORAL
    DO i=1,KI2-1
      NTMAX=KI2-i
      wtx = 0.0
      wty = 0.0
      wtz = 0.0
      wt = 0.0

      !CICLO PARA INICIAR EL BARRIMIENTO DE PARTICULAS
      DO l=1,N
        !CICLO PARA AVANCE TEMPORAL
        DO j= 1,NTMAX
          wtx = wtx + (cxr(l,i+j)-cxr(l,j))**2
          wty = wty + (cyr(l,i+j)-cyr(l,j))**2
          wtz = wtz + (czr(l,i+j)-czr(l,j))**2
        END DO
      END DO
      TIME = TIM*REAL(i)
      wt = (wtx+wty+wtz)/(REAL(NTMAX)*REAL(N)*6.0)
      DIF = wt/TIME
      WRITE(40,*)TIME,wt,DIF

    END DO
    CLOSE(40)
  END SUBROUTINE WDT
END MODULE W
