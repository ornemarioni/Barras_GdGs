! hola_mundo_sub.f90
SUBROUTINE hola_mundo(msg)
    CHARACTER(LEN=12), INTENT(OUT) :: msg
    msg = 'Hola, mundo!'
END SUBROUTINE
