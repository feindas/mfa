! Aca ire poniendo las ultimas modmificaciones del programa.
!01-12-2011 CM_V8: Se arreglo el problema del seed. La funcion INR250
! cambia el valor de iseed a -1. Esto mejora la secuencia de numeros random.
! se implementa tambien la fucnion wall_time. que mide el tiempo de ejecucion
! del programa.

!11-7-11 CM_2.8: Se arreglaron los problemas asociados a la configuracion inicial
!( habia que pensar en la funcion refold() ). La termalizacion erronea estaba asociada
! a la implementacion incorrecta de las condiciones iniciales sobre los anillos.
! (faltaba setear las velocidades iniciales). Aunque detectado el problema aun no esta resuelto.
! Por ahora funcionan los anillos fijos, (de forma tal que no aportan a la energia cinetica )
! En el caso de anillos moviles hay que pensar como re-calcular la energ.cinetica.
!
!Para la version 2.9 es necesario:
! A- Mejorar el output, para que los datos indispensables puedan leerse con los scripts de medicion.
! B- Independizar los scripts de medicion del arbol de directorios (error del primer dia)

! 3-6-11 CM_2.7: La idea es limpiar el codigo de comentarios espureos y mejorar
! el output para que sea facil de awkear (EJ: que los comentarios que
! especifican un numero correspondan a un solo campo.
!          numero de beads: 14     ---->   num_beads: 14 )
! se tocaron los archivos:
! chain_fftw.f90 (se mejora la salida del archivo parametros_fft)
! obser_out.f90 (se imprimen mas variables en parametros_fft)
!
!




! Hasta 2.6 se logro una cadea que funciona con STORE =1 y STORE = 0
! esta implementado correctamente:
! fftw --> difrac_promedio
! las condiciones de contorno(periódicas, bead N interactuando con bead 1) --> film_xmol.xyz
! los potenciales de FENE y Lennar-Jones se computan correctamente (habia un bug
! que sobre escribia el epsilon de LJ para los beads de la cadeja. (esta linea
! se borro y se dejo un comentario de WARNING)
