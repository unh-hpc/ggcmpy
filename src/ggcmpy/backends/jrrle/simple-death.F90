! ----------------------------------------
subroutine death()
!----------------------------------------
!...... causes the sudden death syndrome
   implicit none

   write(0,*) 'ERROR occurred -- entering death.'
   stop
end
