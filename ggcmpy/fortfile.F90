
#define FSEEKABLE

#ifdef FSEEKABLE
#define HAVE_FSEEK 1
#define HAVE_FTELL 1
#endif

module c_string
   implicit none
   private

   public :: c_strlen, c_strdup, c_strdup_f

   interface
      function c_strlen(cstr) bind(c, name='strlen')
         use iso_c_binding
         type(c_ptr), intent(in), value :: cstr
         integer(c_size_t) :: c_strlen
      end function c_strlen

      function c_strdup(cstr) bind(c, name='strdup')
         use iso_c_binding
         type(c_ptr), intent(in), value :: cstr
         type(c_ptr) :: c_strdup
      end function c_strdup
   end interface

contains
   function c_strdup_f(fstr)
      use iso_c_binding
      character(kind=C_CHAR, len=*), intent(in) :: fstr
      type(c_ptr) :: c_strdup_f

      character(kind=C_CHAR, len=:), allocatable, target :: tmp

      tmp = trim(fstr) // c_null_char
      c_strdup_f = c_strdup(c_loc(tmp))
   end function c_strdup_f
end module c_string

subroutine seek(unit,offset,whence,status)
   implicit none
   integer unit
   INTEGER*8 offset
   integer whence
   integer status
   !f2py intent(in) unit
   !f2py intent(in) offset
   !f2py integer optional, intent(in) :: whence=0
   !f2py intent(out) status
   !Non-zero exit status on exit if this routine fails
   status = -1

#ifdef HAVE_FSEEK
   call fseek(unit,offset,whence,status)
#else
   read(unit,'()',ADVANCE='NO',POS=offset,IOSTAT=status)
#endif
   return

end subroutine seek

function tell(unit) result(offset) bind(c)
   use iso_c_binding
   implicit none
   integer(kind=C_INT), value, intent(in) :: unit
   integer(kind=C_LONG) :: offset
   offset = -1

#ifdef HAVE_FTELL
   call ftell(unit,offset)
#else
   inquire(UNIT=unit,POS=offset)
#endif
end function tell

subroutine freefileunit(uu,funit)
   IMPLICIT NONE
   integer uu
   !f2py integer optional uu=-1
   integer funit
   !f2py intent(out) funit

   logical isopen
   integer i

   ! Specifically check for this unit
   if(uu.gt.0)then
      inquire(unit=uu,opened=isopen)
      if(.not.isopen)then
         funit=uu
         return
      endif
   endif

   !look for a free file unit
   do i=10,10000
      if((i.ge.100) .and. (i.le.102))then
         !some implementations reserve these units
         cycle
      endif
      inquire(unit=i,opened=isopen)
      if(.not.isopen)then
         funit=i
         exit
      endif
   enddo

end subroutine freefileunit

subroutine frewind(funit, debug) bind(c)
   use iso_c_binding
   implicit none
   integer(kind=C_INT), value :: funit, debug
   if(debug.gt.0) print*,"Rewinding unit:",funit

   rewind(funit)
end subroutine frewind

subroutine fbackspace(funit,debug)
   IMPLICIT NONE
   integer :: funit
   !f2py intent(in) funit
   integer debug
   !f2py integer optional debug=0

   if(debug.gt.1) print*,"backspacing unit:",funit

   backspace(funit)
end subroutine fbackspace

function fopen(fname, debug) result(funit) bind(c)
   !TODO: open as binary or for appending
   ! open a fortran file
   use iso_c_binding
   use c_string
   implicit none
   type(c_ptr), value :: fname
   integer, value, intent(in) :: debug
   integer :: funit

   integer :: uu = -1
   integer :: len, openstat
   logical isopen
   character*10 access_method
   character(len=:), pointer :: fname_f

   len = c_strlen(fname)
   block
      character(len=len), pointer :: fname_t
      call c_f_pointer(fname, fname_t)
      fname_f => fname_t
   end block

#ifdef HAVE_STREAM
   access_method = 'STREAM'
#else
   access_method = 'SEQUENTIAL'
#endif

   inquire(file=fname_f, opened=isopen, number=funit)
   if(isopen)then
      if(debug.gt.0) print*,"file already opened:",fname
      if(debug.gt.0) print*,"associated unit:",funit
      return
   endif

   call freefileunit(uu,funit)
   open(unit=funit,file=fname_F, status='UNKNOWN', form='FORMATTED', access=access_method, IOSTAT=openstat)

   if (openstat.ne.0) then
      funit = -1 * openstat
   endif

   if(debug.gt.0) print*,"opened file:",fname_f
   if(debug.gt.0) print*,"associated unit:",funit
   if(debug.gt.0) print*,"access method:",access_method

end function fopen

function fadvance_one_line(funit,debug) result(success) bind(c)
   !advance one line in the file
   use iso_c_binding
   implicit none
   integer(C_INT), value, intent(in) :: funit, debug
   integer(C_INT) :: success
   character(len=1) :: a

   if(debug.gt.1) print*,"advancing unit by 1 line:",funit

   success=0
   read(funit,'(A)',err=777,end=777) a
   success=1
777 continue ! continue on error

end function fadvance_one_line

function fisopen(funit, debug) bind(c)
   ! Check if a fortran file unit is open
   use iso_c_binding
   implicit none
   integer(c_int), value :: funit, debug
   integer(c_int) :: fisopen

   logical :: o

   inquire(unit=funit,opened=o)
   if(o)then
      fisopen=1
   else
      fisopen=0
   endif
end function fisopen

subroutine fclose(funit, debug) bind(c)
   ! close a fortran file
   use iso_c_binding
   implicit none
   integer(c_int), value :: funit, debug

   logical isopen
   if(debug.gt.0) print*,"closing unit:",funit

   inquire(unit=funit,opened=isopen)
   if(isopen)then
      close(unit=funit)
   else
      print*,"Cannot close unit that isn't open",funit
   endif

end subroutine fclose

 !TODO: fseek a binary unit
