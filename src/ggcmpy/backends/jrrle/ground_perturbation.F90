module ground_perturbation_m
   implicit none
   private

   public :: iopar, gradpt, ground

   integer :: np, nt
   real, allocatable, dimension(:) :: hh1, hh2, hh3, f1, f2, awork

   integer, parameter :: NIOWK_ = 19000000
contains

! -------------------------------------------------------------------
   subroutine iopar_initialize(np_, nt_)
! -------------------------------------------------------------------
      integer, intent(in) :: np_, nt_

      integer :: npt
      logical, save :: inited = .false.

      if (.not. inited) then
         np = np_
         nt = nt_
         npt = np * nt
         allocate(hh1(npt+1), hh2(npt+1), hh3(npt+1), f1(npt+1), f2(npt+1))
         allocate(awork(NIOWK_))
         inited = .true.
      end if
      if (np /= np_ .or. nt /= nt_) then
         write(0,*) 'Error in iopar_initialize: np,nt do not match'
         stop
      end if
   end subroutine iopar_initialize

! -------------------------------------------------------------------
   subroutine iopar_finalize()
! -------------------------------------------------------------------

      deallocate(hh1, hh2, hh3, f1, f2)
      deallocate(awork)
   end subroutine iopar_finalize


! -------------------------------------------------------------------
   subroutine iopar(np,nt,fac,pot,sigp,sigh,ep,et, &
         ctiot,ctiop,tau,ctaup,ctaut,cpolp,cpolt,delbr,delbp,delbt,xjh)
! -------------------------------------------------------------------
      integer, intent(in) :: np, nt
      !f2py intent(hide) :: np, nt
      real, dimension(np, nt), intent(in) :: fac, pot, sigp, sigh
      real, dimension(np, nt), intent(out) :: tau
      real, dimension(np, nt), intent(out) :: ctiot, ctiop, ctaup, ctaut, cpolp, cpolt
      real, dimension(np, nt), intent(out) :: ep, et, delbr, delbp, delbt, xjh

      real t, cost, sint, sini, dp, dt, f, p, pi, re
      real spp, spt, stt
      integer i, ip, it, iv, niogo, niopy, niowk, niox, nioy
      integer nniowk

      call iopar_initialize(np, nt)

      ! ..... electric field
      call gradpt(pot,np,nt,et,ep)

      pi=4.0*atan(1.0)
      dp=2.0*pi/float(np-1)
      dt=pi/float(nt-1)
      re=6372.0e3
      f=1.0/re

      iv=0
      i=0
      do it=1,nt
         t=float(it-1)*dt
         cost=cos(t)
         sint=sin(t)
         sini=amax1(0.05,2.0*abs(cost)/sqrt(1.0+3.0*cost*cost))
         do ip=1,np
            p=float(ip-1)*dp
            i=i+1
            f1(i)=f*f
            f2(i)=0.0
            !..... total ionospheric current
            stt=sigp(ip,it)/(sini*sini)
            spt=sigh(ip,it)/sini
            spp=sigp(ip,it)
            ctiot(ip,it)=stt*et(ip,it)+spt*ep(ip,it)
            ctiop(ip,it)=spp*ep(ip,it)-spt*et(ip,it)
            !..... joule heat rate
            xjh(ip,it)=ctiot(ip,it)*et(ip,it)+ctiop(ip,it)*ep(ip,it)
            delbt(ip,it)=fac(ip,it)
         end do
      end do

      !...... solve poisson eq for the potential of the irrotational
      !       (poloidal) current
      if(iv.ne.0)write(0,*)' solving poisson equation '
      niox=20
      nioy=10
      niogo=4
      niopy=2*nioy
      nniowk =(size(awork,1)/2)-10
      niowk=nniowk/2
      !.... set up solver
      call io4(0,niox,nioy,niogo,niopy,f1,f2,delbt,tau,np,nt,awork,niowk,hh1,hh2,hh3)
      !.... factor matrix
      call io4(1,niox,nioy,niogo,niopy,f1,f2,delbt,tau,np,nt,awork,niowk,hh1,hh2,hh3)
      !.... solve matrix and interpolate
      call io4(2,niox,nioy,niogo,niopy,f1,f2,delbt,tau,np,nt,awork,niowk,hh1,hh2,hh3)

      !...... irrotational (poloidal)  current
      call gradpt(tau,np,nt,ctaut,ctaup)

      !...... toroidal (div free) part
      do it=1,nt
         do ip=1,np
            cpolt(ip,it)=ctiot(ip,it)-ctaut(ip,it)
            cpolp(ip,it)=ctiop(ip,it)-ctaup(ip,it)
         end do
      end do
      !..... ground magnetic perturbation
      call ground(np,nt,cpolp,cpolt,delbp,delbt,delbr)

      return
   end subroutine iopar


! -------------------------------------------------------------------
   subroutine gradpt(pot,np,nt,et,ep)
! -------------------------------------------------------------------
! calculate meridional and azimuthal electric field from
! the potential, i.e.  (ep,et) = -grad pot

      integer, intent(in) :: np, nt
      !f2py intent(hide) :: np, nt
      real, dimension(np, nt), intent(in) :: pot
      real, dimension(np, nt), intent(out) :: ep, et

      integer ip,it,it0,it1,ip0,ip1
      real pi,dp,dt,f,st,t

      pi=4.0*atan(1.0)
      dp=2.0*pi/float(np-1)
      dt=pi/float(nt-1)
      f=1.0/6372.0e3
      do ip=1,np
         do it=1,nt
            it0=max0(1,it-1)
            it1=min0(nt,it+1)
            et(ip,it)=-f*(pot(ip,it1)-pot(ip,it0))/(float(it1-it0)*dt)
            ep(ip,it)=0.0
            if( (it.gt.1) .and. (it.lt.nt) ) then
               t=float(it-1)*dt
               st=sin(t)
               ip0=ip-1
               if(ip0.lt.1)ip0=np-1
               ip1=ip+1
               if(ip1.gt.np)ip1=2
               ep(ip,it)=-f*(pot(ip1,it)-pot(ip0,it))/(2.0*dp*st)
            endif
         end do
      end do
   end subroutine gradpt

! ----------------------------------------------------------------
   subroutine ground(np,nt,cpolp,cpolt,delbp,delbt,delbr)
!----------------------------------------------------------------
      integer, intent(in) :: np, nt
      !f2py intent(hide) :: np, nt
      real, dimension(np, nt), intent(in) :: cpolp, cpolt
      real, dimension(np, nt), intent(out) :: delbp, delbt, delbr

      real, dimension(500) :: sp,cp,st,ct,st1,st2,ddf
      integer ip,it,jp,jt,nt1,nt2
      real pi,dp,dt,f
      real xmue0,re,reh,h,dd,d2
      real x0,y0,z0,x1,y1,z1,r3i
      real xjx,xjy,xjz,rx,ry,rz,x,y,z
      real dbx,dby,dbz,dbr,dbp,dbt
      real df,p,t,tdeg,tmp_1,xjp,xjr,xjt

      ! write(0,*)'start ground',np,nt
      pi=4.0*atan(1.0)
      xmue0=4.0e-7*pi
      dp=2.0*pi/float(np-1)
      dt=pi/float(nt-1)
      re=6372.0e3
      h=90.0e3    !  ionosphere height
      reh=re+h
      dd=600.0e3   !  radius to which to integrate
      d2=dd*dd

      do ip=1,np
         p=float(ip-1)*dp
         sp(ip)=sin(p)
         cp(ip)=cos(p)
      end do
      do it=1,nt
         t=float(it-1)*dt
         st(it)=sin(t)
         ct(it)=cos(t)
         st1(it)=sin(t-0.5*dt)
         st2(it)=sin(t+0.5*dt)
         ddf(it)=0.5*re*re*dp*dt*(st1(it)+st2(it))*xmue0/(4.0*pi)
      end do

      do it=1,nt
         nt1=max0(2,it-14)
         nt2=min0(nt-1,it+14)
         t=float(it-1)*dt
         tdeg=t*180.0/pi
         do ip=1,np
            delbr(ip,it)=0.0
            delbp(ip,it)=0.0
            delbt(ip,it)=0.0
            if(tdeg.gt.45.0.and.tdeg.lt.135.0) cycle
            p=float(ip-1)*dp
            x0=cp(ip)*st(it)*re
            y0=sp(ip)*st(it)*re
            z0=ct(it)*re
            dbx=0.0
            dby=0.0
            dbz=0.0

            !..... integration loop
            do jt=nt1,nt2
               do jp=1,np-1
                  x1=cp(jp)*st(jt)*reh
                  y1=sp(jp)*st(jt)*reh
                  z1=ct(jt)*reh
                  dd=(x1-x0)**2+(y1-y0)**2+(z1-z0)**2
                  if(dd.gt.d2) cycle
                  rx=x0-x1
                  ry=y0-y1
                  rz=z0-z1
                  r3i=1.0/(dd*sqrt(dd))
                  df=ddf(jt)*r3i
                  xjr=0.0
                  xjp=cpolp(ip,it)
                  xjt=cpolt(ip,it)
                  call vec_pol_cart(xjr,xjp,xjt,xjx,xjy,xjz,sp(jp),cp(jp),st(jt),ct(jt))
                  call cross11(xjx,xjy,xjz,rx,ry,rz,x,y,z)
                  dbx=dbx+df*x
                  dby=dby+df*y
                  dbz=dbz+df*z
               end do
            end do ! integration loop

            call vec_cart_pol(dbx,dby,dbz,dbr,dbp,dbt,sp(ip),cp(ip),st(it),ct(it))
            delbr(ip,it)=dbr
            delbp(ip,it)=dbp
            delbt(ip,it)=dbt
         end do
      end do

   contains

      subroutine vec_cart_pol(AX,AY,AZ,AR,AP,AT,SP,CP,ST,CT)
         real AX,AY,AZ,SP,CP,ST,CT
         real AR,AP,AT

         real tmp_1

         tmp_1=AX*CP+AY*SP
         AR=tmp_1*ST+AZ*CT
         AT=tmp_1*CT-AZ*ST
         AP=AY*CP-AX*SP

      end subroutine vec_cart_pol


      subroutine vec_pol_cart(AR,AP,AT,AX,AY,AZ,SP,CP,ST,CT)
         real AX,AY,AZ,SP,CP,ST,CT
         real AR,AP,AT

         real tmp_1

         ! tmp_1=AR*ST+AT*CT
         tmp_1=AT*CT
         AX=tmp_1*CP-AP*SP
         AY=tmp_1*SP+AP*CP
         ! AZ=AR*CT-AT*ST
         AZ=AT*ST

      end subroutine vec_pol_cart

      subroutine cross11(AX,AY,AZ,BX,BY,BZ,CX,CY,CZ)
         real AX,AY,AZ,BX,BY,BZ,CX,CY,CZ

         CX=AY*BZ-AZ*BY
         CY=AZ*BX-AX*BZ
         CZ=AX*BY-AY*BX

      end subroutine cross11

   end subroutine ground

end module ground_perturbation_m
