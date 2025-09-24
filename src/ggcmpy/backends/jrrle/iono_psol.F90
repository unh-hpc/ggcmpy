!......  This module contains the ionosphere potential solver
!        and some helper routines.
!        Completely self contained, no side effects through common
!        blocks, but takes hash-constant parameters.
!        Temp storage must be maintained between calls.
!
!        JR  3/7/2000:  from ionos08.for
!
! ---------------------------------------------------------------------
subroutine io4(ijob,nx,ny,ngord,npany,sigp,sigh,xjpa,pot,nphi,nthe,wk,nwk, t1,t2,t3)
! ---------------------------------------------------------------------
   real sigp(*),sigh(*),xjpa(*),pot(*),wk(*)
   real t1(nthe/2,nphi)
   real t2(nthe/2,nphi)
   real t3(nthe/2,nphi)
   parameter (pi=3.1415926)

   mx=nthe/2
   my=nphi
   xx2=pi*(89.95)/180.0

!    call second1(cpu1)

!.... initialize
   if(ijob.eq.0) then
      call iosol2(0,nx,ny,xx2,ngord,npany,t1,t2,t3,mx,my,wk(1),nwk)
      call iosol2(0,nx,ny,xx2,ngord,npany,t1,t2,t3,mx,my,wk(1+nwk),nwk)
   endif

!.... set up matrix and factor it
   if(ijob.eq.1) then
      call trio1n(sigp,t1,nphi,nthe,mx,1.0)
      call trio1n(sigh,t2,nphi,nthe,mx,1.0)
      call iosol2(1,nx,ny,xx2,ngord,npany,t1,t2,t3,mx,my,wk(1),nwk)
      call trio1s(sigp,t1,nphi,nthe,mx,1.0)
      call trio1s(sigh,t2,nphi,nthe,mx,1.0)
      call iosol2(1,nx,ny,xx2,ngord,npany,t1,t2,t3,mx,my,wk(1+nwk),nwk)
   endif

!.... solve equation with matrix factored
   if(ijob.eq.2) then
      call trio1n(xjpa,t3,nphi,nthe,mx,-1.0)
      call iosol2(2,nx,ny,xx2,ngord,npany,t1,t2,t3,mx,my,wk(1),nwk)
      call iosol2(4,nx,ny,xx2,ngord,npany,t1,t2,t3,mx,my,wk(1),nwk)
      call trio2n(t3,pot,nthe,nphi,mx)
      call trio1s(xjpa,t3,nphi,nthe,mx,-1.0)
      call iosol2(2,nx,ny,xx2,ngord,npany,t1,t2,t3,mx,my,wk(1+nwk),nwk)
      call iosol2(4,nx,ny,xx2,ngord,npany,t1,t2,t3,mx,my,wk(1+nwk),nwk)
      call trio2s(t3,pot,nthe,nphi,mx)
   endif

!    call second1(cpu2)
!      write(0,*)'io4: ',ijob,cpu2-cpu1

   return
end
! ------------------------------------------------------------------
subroutine trio1n(a,b,np,nt,nt2,s)
! ------------------------------------------------------------------
!.... transpose array for the northern hemisphere
   real a(np,nt),b(nt2,np)
   it2=1
   do it=1,nt2
      do ip=1,np
         b(it,ip)=a(ip,it2)*s
      end do
      it2=it2+1
   end do
   return
end
! ------------------------------------------------------------------
subroutine trio2n(a,b,nt,np,nt2)
! ------------------------------------------------------------------
!.... copy back and flip into total potential array, northern part
   real a(nt2,np),b(np,nt)
   it2=1
   do it=1,nt2
      do ip=1,np
         b(ip,it2)=a(it,ip)
      end do
      it2=it2+1
   end do
   return
end
! ------------------------------------------------------------------
subroutine trio1s(a,b,np,nt,nt2,s)
! ------------------------------------------------------------------
!.... transpose array for the southern hemisphere
   real a(np,nt),b(nt2,np)
   it2=nt
   do it=1,nt2
      do ip=1,np
         b(it,ip)=a(ip,it2)*s
      end do
      it2=it2-1
   end do
   return
end
! ------------------------------------------------------------------
subroutine trio2s(a,b,nt,np,nt2)
! ------------------------------------------------------------------
!.... copy back and flip into total potential array, southern part
   real a(nt2,np),b(np,nt)
   it2=nt
   do it=1,nt2
      do ip=1,np
         b(ip,it2)=a(it,ip)
      end do
      it2=it2-1
   end do
   return
end
! ------------------------------------------------------------------
subroutine iosol2(ijob,nx,ny,xx2,mgord,mpany,sigp,sigh,xjpa,mx,my,wk,mwk)
! ---------------------------------------------------------------------
!   galerkin pseudospectral routine to solve the ionosphere
!   potential equation
!
!     d/dx( a(x,y) * df/dx + b(x,y) * df/dy )
!   + d/dy( c(x,y) * df/dx + d(x,y) * df/dy )  =  rhs
!
!   with boundary conditions:  f(0,.)=0, df(xx2,.)/dx=0, periodic y
!
!   basis functions are:
!
!   B(i,j,x,y) = H(mod(i-1,4),x) * F(j,y) , i=1,..,4*nx;  j=1,...,2*ny
!
!   with:   H(i,x) : Hermite polynominal with support t1,..,t2
!                    cubics with continuous first derivative
!           F(j,y) : sin((((j-1)/2)+1)*y), j odd
!                    cos( ((j-1)/2   )*y), j even
!
!   nx:     no of panels in x (theta)
!   ny:     no of modes in y (phi)
!   xx2:    upper x boundary (radians)
!   mgord:  order of gauss integration
!   mpany:  no of panels in y for integrations
!   sigp:   array holding pedersen conductivity
!   sigh:   array holding hall conductivity
!   xjpa:   array holding parallel current
!           for ijob=4: potential
!   mx,my:  dimensions of sigp,sigh,xjpa
!   wk:     workspace, must not be touched between calls
!   mwk:    size of workspace, will stop if not enough
!
!... ijob=0   initialize
!... ijob=1   set up global matrix and factor it
!... ijob=2   set up rhs and solve
!... ijob=4   evaluate solution on reg grid xjpa
!
!    written by: Jimmy Raeder, IGPP/UCLA, Nov 1992
!
   use constants_m, only: tty
   real wk(*),sigp(mx,my),sigh(mx,my),xjpa(mx,my)

   common /iosola/neq,nb,ndof,ndof1,ny2,ny4,na,ngord,npany, nwa, &
         nwrhs,nwem,nwleq,nwemat,ncof,ncofr,nwcof,nwcofr,nwcfx, nwcfy, &
         nwtx,nwsol
   character*80 rec
   save
   parameter (pi=3.1415926)


!... initialize
   if(ijob.eq.0) then
      ndof=4*2*ny
      ndof1=ndof+1
      ny2=2*ny
      ny4=2*ny2
      npany=mpany
      ngord=mgord
!.... do a mock assembly to determine neq and bandwidth
      neq=0
      nb=0
      do kx=1,nx
         do jy=1,ndof
            if(kx.eq.1) then
               iy=jy
            else
               iy=ny4*(kx-1)+jy
            endif
            do jx=1,ndof
               if(kx.eq.1) then
                  ix=jx
               else
                  ix=ny4*(kx-1)+jx
               endif
               nb=max0(nb,iabs(ix-iy))
               neq=max0(neq,ix)
               neq=max0(neq,iy)
            end do
         end do
      end do
!      write(0,*)' iosol2 neq,nb,na ',neq,nb,neq*(nb+nb+1)
!..... workspace assignments
      nw = 1
      nwtx = nw
      nw = nw + nx+1
      nwa = nw
      na = neq*(nb+nb+1)
      nw = nw + na
      nwrhs = nw
      nw = nw + neq
      nwem = nw
      nw = nw + (ndof+1)*ndof
      nwleq = nw
      nw = nw + neq*(nb+1)
      nwemat = nw
      nw = nw + 4*npany*ngord*ny
      nwcof = nw
      ncof = 8*ngord*ngord*nx*npany
      nw = nw + ncof
      ncofr = 3*ngord*npany
      nwcofr = nw
      nw = nw + ncofr
      nwsol = nw
      nw = nw + neq
      nwcfx = nw
      nw = nw + 5*mx
      nwcfy = nw
      nw = nw + 2*ny*my

!      write(tty,*)'iosol2 nw mwk neq nb ',nw,mwk,neq,nb
      if(nw.gt.mwk) call death()
!.... create panels in x (theta)
      x1=0.0
      x2=xx2
      if (.false.) then
         do 120 kx=1,nx+1
            if (.false.) then
               if(kx.le.60) then
                  wk(nwtx+kx-1)=(pi/180.0)*0.5*float(kx-1)
                  xxl=wk(nwtx+kx-1)
               else
                  wk(nwtx+kx-1)=wk(nwtx+60-1)+(xx2-xxl)*float(kx-60)/float(nx-60+1)
               endif
!      write(0,*)' th-elements ',kx,s,wk(nwtx+kx-1)*57.29
            endif
            if (.false.) then
               s=float(kx-1)/float(nx)
               s1=s*s
               s2=0.5*s
               ss=amax1(s1,s2)
               wk(nwtx+kx-1)=x1+(x2-x1)*ss
!      write(0,*)' th-elements ',kx,s,wk(nwtx+kx-1)*57.29
            endif
            if (.false.) then
               s=float(kx-1)/float(nx)
               s1=s*s
               s2=0.5*s
               ss=amax1(s1,s2)
               wk(nwtx+kx-1)=x1+(x2-x1)*ss
!      write(0,*)' th-elements ',kx,s,wk(nwtx+kx-1)*57.29
            endif
120      continue
      endif
!      write(0,*)'call mkel ',nx,nwtx
      call mkel(nx+1,wk(nwtx),x1,x2)
!      write(0,*)'call mkgaus ',nx,ny,npany,ngord,nwtx,nwcof,nwcofr
      call mkgaus(nx,ny,npany,ngord,wk(nwtx),wk(nwcof),wk(nwcofr))
!      write(0,*)'call mkcoe1',nx,ny,npany,ngord,nwtx,nwcof,nwcofr
      call mkcoe1(nx,ny,npany,ngord,wk(nwtx),wk(nwcof),wk(nwcofr),xx2,sigp,sigh,xjpa,mx,my,0)
!      write(0,*)'calli eval2 ',nwtx,nwsol,neq,nx,ny,mx,my,nwcfx,nwcfy
      call eval2(0,wk(nwtx),wk(nwsol),neq,nx,ny,xjpa,mx,my,xx2,wk(nwcfx),wk(nwcfy))
!      write(0,*)'eval2 back'
      call emat(wk(nwem),wk(nwtx),1,nx,ny,ndof,ngord,npany, wk(nwemat), 0 ,wk(nwcof),wk(nwcofr))
      return
   endif



   if(ijob.eq.1.or.ijob.eq.2) then

      call xzero(wk(nwrhs),neq)
      if(ijob.eq.1)call xzero(wk(nwa),na)

      call mkcoe1(nx,ny,npany,ngord,wk(nwtx),wk(nwcof),wk(nwcofr),xx2, sigp,sigh,xjpa,mx,my,ijob)

!..... assemble global matrix or rhs
      do k=1,nx
         call xzero(wk(nwem),ndof*ndof1)
         call emat(wk(nwem),wk(nwtx),k,nx,ny,ndof,ngord,npany, wk(nwemat),ijob,wk(nwcof),wk(nwcofr))
         call addmat(wk(nwa),wk(nwrhs),wk(nwem),neq,nb,ndof,ny4,k,ijob)
      end do

!.... factor matrix
      if(ijob.eq.1) then
         call leqt1b(wk(nwa),neq,nb,nb,neq,wk(nwsol),1,neq,1,wk(nwleq),ier)
         if(ier.ne.0) then
            write(tty,*)' leqt1b ier= ',ier
         endif
!      if(ier.eq.129) call death()
      endif

!.... solve equations
      if(ijob.eq.2) then
         call xcopy(wk(nwrhs),wk(nwsol),neq)
         call leqt1b(wk(nwa),neq,nb,nb,neq,wk(nwsol),1,neq,2,wk(nwleq),ier)
         if(ier.ne.0) then
            write(tty,*)' leqt1b ier= ',ier
         endif
         if(ier.eq.129) call death()
      endif

   endif

   if(ijob.eq.4) then
      call eval2(4,wk(nwtx),wk(nwsol),neq,nx,ny,xjpa,mx,my,xx2,wk(nwcfx),wk(nwcfy))
   endif


   return
end
! -----------------------------------------------------------------------
subroutine mkel(nn,w,x1,x2)
! -----------------------------------------------------------------------
   real w(nn)
   deg=45.0/atan(1.0)
   do i=1,nn
      w(i)=deg*(x1+float(i-1)*(x2-x1)/float(nn-1))
   end do
   if (.true.) then
      s=0.0
      do i=2,nn
         t=deg*(x1+float(i-1)*(x2-x1)/float(nn-1))
         s=s+(1.0/xmkel(t))
      end do
      w(1)=x1
      do i=2,nn
         w(i)=w(i-1)+(s/xmkel(w(i)))
      end do
      s=(x2-x1)/(w(nn)-w(1))
      do i=2,nn
         w(i)=x1+s*(w(i)-x1)
      end do
      if (.false.) then
         dd=0.0
         do i=1,nn
            if(i.gt.1)dd=w(i)-w(i-1)
            write(0,'(a,i3,1x,3(1x,g12.5))')'mkel ',i,deg*w(i),deg*dd,deg*xmkel(w(i))*dd
         end do
      endif
   endif
   if (.false.) then
      f=1.0
      do it=1,10000
         do i=2,nn-1
            s1=0.5*(w(i-1)+w(i))
            s2=0.5*(w(i)+w(i+1))
            d1=w(i)-w(i-1)
            d2=w(i+1)-w(i)
            df1=d1*xmkel(s1)
            df2=d2*xmkel(s2)
!      if(mod(it,500).eq.0)write(0,*)'mkel: ',i,df1,df2,d1,d2
            if(df1.gt.df2) then
               w(i)=w(i)-0.2*f*d1
            else
               w(i)=w(i)+0.2*f*d2
            endif
         end do
         if(mod(it,1000).eq.0)f=0.8*f
      end do
      dd=0.0
      do i=1,nn
         if(i.gt.1)dd=w(i)-w(i-1)
         write(0,'(a,i3,1x,3(1x,g12.5))')'mkel ',i,w(i),dd,xmkel(w(i))*dd
      end do
      do i=1,nn
         w(i)=w(i)/deg
      end do
   endif
   return
end
!-----------------------------------------------------------------------
real function xmkel(x)
!-----------------------------------------------------------------------
   d=-1.0*( ((x-40.0)/13.0)**2 )
   xmkel = 1.0 + 2.0*exp( d ) - 0.85*x/90.0
   return
end
!-----------------------------------------------------------------------
subroutine eval2(ijob,tx,sol,neq,nx,ny,pot,mx,my,xend,cfx,cfy)
!-----------------------------------------------------------------------
   real tx(*),sol(*),pot(mx,my),cfx(5,mx),cfy(2*ny,my)
   real bb(4),bbx(4)
   parameter (pi=3.1415926)

!      write(0,*)'evale: ',neq,nx,ny,mx,my,xend
   ny2=2*ny
   ny4=2*ny2
   if(ijob.eq.0) then
      do ix=1,mx
         x=xend*float(ix-1)/float(mx-1)
         kx=nx
!... find x interval
         do j=1,nx
            if(tx(j).le.x.and.tx(j+1).ge.x) kx=j
         end do
         call hermite(bb,bbx,x,tx(kx),tx(kx+1))
         cfx(1,ix)=bb(1)
         cfx(2,ix)=bb(2)
         cfx(3,ix)=bb(3)
         cfx(4,ix)=bb(4)
         cfx(5,ix)=float(kx)+0.1
      end do
      do iy=1,my
         y=(-pi)+2.0*pi*float(iy-1)/float(my-1)
         do jy=1,ny2
            itmp1=(jy-1)/2
            itmp2=mod(jy,2)
            tmp1=float(itmp1+itmp2)
            if(itmp2.eq.0) fy=cos(y*tmp1)
            if(itmp2.eq.1) fy=sin(y*tmp1)
            cfy(jy,iy)=fy
         end do
      end do
      return
   endif


   if(ijob.eq.4) then
      do iy=1,my
         do ix=1,mx
            kx=cfx(5,ix)
            pot(ix,iy)=0.0
            do j1=1,ny2
               do i1=1,4
                  jy=ny2*(i1-1)+j1
                  if(kx.eq.1) then
                     ky=jy
                  else
                     ky=ny4*(kx-1)+jy
                  endif
                  pot(ix,iy)=pot(ix,iy)+sol(ky)*cfx(i1,ix)*cfy(j1,iy)
               end do
            end do
         end do
      end do
   endif

   return
end
! -----------------------------------------------------------------------
subroutine mkcoe1(nx,ny,npany,ngord,tx,coef,coefr,xx2,sigp,sigh,xjpa,mx,my,ijob)
! -----------------------------------------------------------------------
! .... create gauss legendre integration points and weights
   use constants_m, only: tty
   real coef(ngord,ngord,npany,nx,8)
   real coefr(ngord,npany,3)
   real sigp(mx,my),sigh(mx,my),xjpa(mx,my)
   real tx(*)
   parameter (NIOGO=4) ! max order of gauss legendre integration
   parameter (NIOX=40) ! max no of panels in theta for one hemisphere
   common /mkcoea/sico(3,NIOGO,NIOX)
   common /mkcoeb/xa,xb,ya,yb,hx,hy,hxi,hyi
   common /mkcoec/mxm1,mym1
   parameter (pi=3.1415926)

   if (ngord .gt. NIOGO) then
      write(tty,*)' mkcoe1: ngord too large ',ngord,NIOGO
      call death()
   endif
   if (nx .gt. NIOX) then
      write(tty,*)' mkcoe1: nx too large ',nx,NIOX
      call death()
   endif

!      call minmax(sigp,mx*my,zmin,zmax)
!      write(tty,*)' minmax sigp  ',zmin,zmax
!      call minmax(sigh,mx*my,zmin,zmax)
!      write(tty,*)' minmax sigh  ',zmin,zmax
!      call minmax(xjpa,mx*my,zmin,zmax)
!      write(tty,*)' minmax xjpa  ',zmin,zmax
!
!      write(0,*)'mkcoe:',ijob,mx,my,nx,ny,xx2

   if(ijob.eq.0) then
      xa=0.0
      xb=xx2
      ya=(-pi)
      yb=pi
      mxm1=mx-1
      hx=(xb-xa)/float(mxm1)
      hxi=1.0/hx
      mym1=my-1
      hy=(yb-ya)/float(mym1)
      hyi=1.0/hy
!.... make sin/cos/sini table
      do kx=1,nx
         do igx=1,ngord
            x=coef(igx,1,1,kx,2)
            st=sin(x)
            ct=cos(x)
            si=2.0*ct/sqrt(1.0+3.0*ct*ct)
            sico(1,igx,kx)=st
            sico(2,igx,kx)=ct
            sico(3,igx,kx)=si
         end do
      end do
      return
   endif



   if(ijob.eq.1) then
      do kx=1,nx
         do igx=1,ngord
            x=coef(igx,1,1,kx,2)
            st=sico(1,igx,kx)
            ct=sico(2,igx,kx)
            si=sico(3,igx,kx)
            tmpx=hxi*(x-xa)
            ix=tmpx
            sx2=tmpx-float(ix)
            sx1=1.0-sx2
            ix=max0(1,ix+1)
            ix=min0(mxm1,ix)
            do iyp=1,npany
               do igy=1,ngord
                  y=coef(igx,igy,iyp,kx,3)
                  tmpy=hyi*(y-ya)
                  iy=tmpy
                  sy2=tmpy-float(iy)
                  sy1=1.0-sy2
                  iy=max0(1,iy+1)
                  iy=min0(mym1,iy)
                  sp=sy1*(sx1*sigp(ix,iy)+sx2*sigp(ix+1,iy)) +sy2*(sx1*sigp(ix,iy+1)+sx2*sigp(ix+1,iy+1))
                  sh=sy1*(sx1*sigh(ix,iy)+sx2*sigh(ix+1,iy)) +sy2*(sx1*sigh(ix,iy+1)+sx2*sigh(ix+1,iy+1))
                  coef(igx,igy,iyp,kx,5)=st*sp/(si*si)
                  coef(igx,igy,iyp,kx,6)=sh/si
                  coef(igx,igy,iyp,kx,7)=(-coef(igx,igy,iyp,kx,6))
                  coef(igx,igy,iyp,kx,8)=sp/st
               end do
            end do
         end do
      end do
      x=tx(nx+1)
      st=sin(x)
      ct=cos(x)
      si=2.0*ct/sqrt(1.0+3.0*ct*ct)
      kx=1
      igx=1
      tmpx=hxi*(x-xa)
      ix=tmpx
      sx2=tmpx-float(ix)
      sx1=1.0-sx2
      ix=max0(1,ix+1)
      ix=min0(mxm1,ix)
      do iyp=1,npany
         do igy=1,ngord
            y=coef(igx,igy,iyp,kx,3)
            tmpy=hyi*(y-ya)
            iy=tmpy
            sy2=tmpy-float(iy)
            sy1=1.0-sy2
            iy=max0(1,iy+1)
            iy=min0(mym1,iy)
            sh=sy1*(sx1*sigh(ix,iy)+sx2*sigh(ix+1,iy)) +sy2*(sx1*sigh(ix,iy+1)+sx2*sigh(ix+1,iy+1))
            coefr(igy,iyp,1)=sh/si
         end do
      end do
   endif

   if(ijob.eq.2) then
      do kx=1,nx
         do igx=1,ngord
            x=coef(igx,1,1,kx,2)
            st=sico(1,igx,kx)
            ct=sico(2,igx,kx)
            si=sico(3,igx,kx)
            tmpx=hxi*(x-xa)
            ix=tmpx
            sx2=tmpx-float(ix)
            sx1=1.0-sx2
            ix=max0(1,ix+1)
            ix=min0(mxm1,ix)
            do iyp=1,npany
               do igy=1,ngord
                  y=coef(igx,igy,iyp,kx,3)
                  tmpy=hyi*(y-ya)
                  iy=tmpy
                  sy2=tmpy-float(iy)
                  sy1=1.0-sy2
                  iy=max0(1,iy+1)
                  iy=min0(mym1,iy)
                  xj=sy1*(sx1*xjpa(ix,iy)+sx2*xjpa(ix+1,iy)) +sy2*(sx1*xjpa(ix,iy+1)+sx2*xjpa(ix+1,iy+1))
                  coef(igx,igy,iyp,kx,1)=st*si*xj
               end do
            end do
         end do
      end do
   endif

   return
end
! ---------------------------------------------------------------------
subroutine mkgaus(nx,ny,npany,ngord,tx,coef,coefr)
! ---------------------------------------------------------------------
!.... create gauss legendre integration points and weights
   real coef(ngord,ngord,npany,nx,8)
   real coefr(ngord,npany,3)
   real tx(*)
   real gpxx(20),gpyy(20),gpwx(20),gpwy(20)
   parameter (pi=3.1415926)

   do kx=1,nx
      call gauleg(tx(kx),tx(kx+1),gpxx,gpwx,ngord)
      do iyp=1,npany
         ty1=(-pi)+pi*2.0*float(iyp-1)/float(npany)
         ty2=(-pi)+pi*2.0*float(iyp )/float(npany)
         call gauleg(ty1,ty2,gpyy,gpwy,ngord)
         do igy=1,ngord
            coefr(igy,iyp,2)=gpyy(igy)
            coefr(igy,iyp,3)=gpwy(igy)
            do igx=1,ngord
               coef(igx,igy,iyp,kx,2)=gpxx(igx)
               coef(igx,igy,iyp,kx,3)=gpyy(igy)
               coef(igx,igy,iyp,kx,4)=gpwx(igx)*gpwy(igy)
            end do
         end do
      end do
   end do

   return
end

! -------------------------------------------------------------------
subroutine emat(em,tx,kk,nx,ny,ndof,ngord,npany, wk,ijob,coef,coefr)
! -------------------------------------------------------------------
! ..... coefficients etc at gausspoints
! ..... 1: rhs 2: gpx 3: gpy  4: ww 5: stt 6: stp  7: spt  8: spp
   real coef(ngord,ngord,npany,nx,8)
   real coefr(ngord,npany,3)
   real em(ndof+1,ndof)
   real tx(*)
   real wk(2,npany,ngord,2*ny)
   real bb(4),bx(4)
   parameter (pi=3.1415926)

   ny2=2*ny
!.... make table for fourier coeff
   if(ijob.eq.0) then
      do iyp=1,npany
         do igy=1,ngord
            yy=coef(1,igy,iyp,1,3)
            do j=1,ny2
               itmp1=(j-1)/2
               itmp2=mod(j,2)
               tmp1=float(itmp1+itmp2)
               if(itmp2.eq.0) wk(1,iyp,igy,j)=cos(yy*tmp1)
               if(itmp2.eq.1) wk(1,iyp,igy,j)=sin(yy*tmp1)
               itmp1=(j-1)/2
               itmp2=mod(j,2)
               tmp1=float(itmp1+itmp2)
               if(itmp2.eq.0) wk(2,iyp,igy,j)=(-tmp1*sin(yy*tmp1))
               if(itmp2.eq.1) wk(2,iyp,igy,j)= tmp1*cos(yy*tmp1)
            end do
         end do
      end do
   endif

!.... numerically integrate element matrix and rhs



   if(ijob.eq.1) then

      tx1=tx(kk)
      tx2=tx(kk+1)
      do 94711 iyp=1,npany
         do 94712 igy=1,ngord
            do 94713 igx=1,ngord
!..... 1: rhs 2: gpx 3: gpy  4: stt 5: stp  6: spt  7: spp
               ww =coef(igx,igy,iyp,kk,4)
               xx =coef(igx,igy,iyp,kk,2)
               yy =coef(igx,igy,iyp,kk,3)
               acof=coef(igx,igy,iyp,kk,5)*ww
               bcof=coef(igx,igy,iyp,kk,6)*ww
               ccof=coef(igx,igy,iyp,kk,7)*ww
               dcof=coef(igx,igy,iyp,kk,8)*ww
               call hermite(bb,bx,xx,tx1,tx2)


               do 94714 j2=1,ny2
                  f20=wk(1,iyp,igy,j2)
                  f21=wk(2,iyp,igy,j2)
                  do 94714 i2=1,4
                     k2=ny2*(i2-1)+j2
                     tmp1=acof*bx(i2)*f20
                     tmp2=bcof*bx(i2)*f20
                     tmp3=ccof*bb(i2)*f21
                     tmp4=dcof*bb(i2)*f21
                     do 94715 j1=1,ny2
                        tmp5=wk(1,iyp,igy,j1)*(tmp1+tmp3)
                        tmp6=wk(2,iyp,igy,j1)*(tmp2+tmp4)
                        do 94715 i1=1,4
                           k1=ny2*(i1-1)+j1
                           em(k1,k2)=em(k1,k2)+bx(i1)*tmp5+bb(i1)*tmp6
94715                continue
94714          continue


94713       continue
94712    continue
94711 continue

!.... boundary conditions

      if(kk.eq.1) then
         do 94682 i2=2,2
            do 94682 j2=1,2*ny
               k2=ny2*(i2-1)+j2
               do 94682 i1=1,4
                  do 94682 j1=1,2*ny
                     k1=ny2*(i1-1)+j1
                     em(k1,k2)=0.0
                     if(k1.eq.k2)em(k1,k2)=1.0
94682    continue
      endif

      if(kk.eq.nx) then
         do 94683 i2=3,3
            do 94683 j2=1,2*ny
               k2=ny2*(i2-1)+j2
               do 94683 i1=1,4
                  do 94683 j1=1,2*ny
                     k1=ny2*(i1-1)+j1
                     em(k1,k2)=0.0
                     if(k1.eq.k2)em(k1,k2)=1.0
94683    continue
      endif

   endif
   if(ijob.eq.2) then

      tx1=tx(kk)
      tx2=tx(kk+1)
      do 94591 iyp=1,npany
         do 94592 igy=1,ngord
            do 94593 igx=1,ngord
!..... 1: rhs 2: gpx 3: gpy  4: stt 5: stp  6: spt  7: spp
               ww =coef(igx,igy,iyp,kk,4)
               xx =coef(igx,igy,iyp,kk,2)
               yy =coef(igx,igy,iyp,kk,3)
               rrhs=coef(igx,igy,iyp,kk,1)*ww
               call hermite(bb,bx,xx,tx1,tx2)


               do 94594 j2=1,ny2
                  f20=wk(1,iyp,igy,j2)
                  f21=wk(2,iyp,igy,j2)
                  do 94594 i2=1,4
                     k2=ny2*(i2-1)+j2
                     em(ndof+1,k2)=em(ndof+1,k2)-rrhs*bb(i2)*f20
94594          continue


94593       continue
94592    continue
94591 continue

!.... boundary conditions

      if(kk.eq.1) then
         do 94562 i2=2,2
            do 94562 j2=1,2*ny
               k2=ny2*(i2-1)+j2
               em(ndof+1,k2)=0.0
94562    continue
      endif

      if(kk.eq.nx) then
         do 94563 i2=3,3
            do 94563 j2=1,2*ny
               k2=ny2*(i2-1)+j2
               em(ndof+1,k2)=0.0
94563    continue
      endif

   endif

   return
end
! ---------------------------------------------------------------------
subroutine addmat(a,rhs,em,neq,nb,ndof,ny4,k,ijob)
! ---------------------------------------------------------------------
   real a(*),rhs(*),em(ndof+1,ndof)

!.... add element matrix to global array
   if(ijob.eq.1) then
      do j=1,ndof
         if(k.eq.1) then
            iy=j
         else
            iy=ny4*(k-1)+j
         endif
         do i=1,ndof
            if(k.eq.1) then
               ix=i
            else
               ix=ny4*(k-1)+i
            endif
            l=(ix-iy+nb)*neq+iy
            a(l)=a(l)+em(i,j)
         end do
      end do
   endif

!.... compile rhs only
   if(ijob.eq.2) then
      do j=1,ndof
         if(k.eq.1) then
            iy=j
         else
            iy=ny4*(k-1)+j
         endif
         rhs(iy)=rhs(iy)+em(ndof+1,j)
      end do
   endif

   return
end
! ---------------------------------------------------------------------
subroutine xzero(a,n)
! ---------------------------------------------------------------------
   real a(n)
   do i=1,n
      a(i)=0.0
   end do
   return
end
! ---------------------------------------------------------------------
subroutine xcopy(a,b,n)
! ---------------------------------------------------------------------
   real a(n),b(n)
   do i=1,n
      b(i)=a(i)
   end do
   return
end

! -----------------------------------------------------
subroutine hermite(bb,bx,x,x1,x2)
! -----------------------------------------------------
   real bb(4),bx(4)
!
!.... hermite polys have:
!                H(1,x1) = 1
!           d/dx/H(2,x1) = 1
!                H(3,x2) = 1
!           d/dx H(4,x2) = 1
!... makes global numbering a lot easier
   z=x-x1
   h=x2-x1
   rh=1.0/h
   s=z*rh
   bb(1) = s*s*(2.0*s-3.0)+1.0
   bb(3) = s*s*(3.0-2.0*s)
   bb(2) = s*(1.0-s)*(h-z)
   bb(4) = s*z*(s-1.0)
   bx(1) = 6.0*(s-1.0)*rh*s
   bx(3) = 6.0*(1.0-s)*rh*s
   bx(2) = (1.0-s)*(1.0-3.0*s)
   bx(4) = s*(3.0*s-2.0)

   return
end
!-------------------------------------------------
subroutine gauleg(x1,x2,x,w,n)
!-------------------------------------------------

   real x1,x2,x(n),w(n)

   parameter (pi=3.1415926)

   eps=1.0e-6

   m=(n+1)/2
   xm=0.5*(x2+x1)
   xl=0.5*(x2-x1)
   do 12 i=1,m
      z=cos(pi*(i-.25)/(n+.5))
      nrep=0
1     continue
      nrep=nrep+1
      p1=1.0
      p2=0.0
      do 11 j=1,n
         p3=p2
         p2=p1
         p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
11    continue
      pp=n*(z*p1-p2)/(z*z-1.0)
      z1=z
      z=z1-p1/pp
      if(abs(z-z1).gt.eps.and.nrep.le.100)go to 1
      x(i)=xm-xl*z
      x(n+1-i)=xm+xl*z
      w(i)=2.0*xl/((1.0-z*z)*pp*pp)
      w(n+1-i)=w(i)
12 continue
   return
end
! ---------------------------------------------------------------
subroutine leqt1b (a,n,nlc,nuc,ia,b,m,ib,ijob,xl,ier)
! --------------------------------------------------------------
   real a(ia,*),xl(n,*),b(ib,*)
   data zero/0./,one/1.0/
   ier = 0
   jbeg = nlc+1
   nlc1 = jbeg
   if (ijob .eq. 2) go to 80
   rn = n
   i = 1
   nc = jbeg+nuc
   nn = nc
   jend = nc
   if (n .eq. 1 .or. nlc .eq. 0) go to 25
5  k = 1
   p = zero
   do j = jbeg,jend
      a(i,k) = a(i,j)
      q = abs(a(i,k))
      if (q .gt. p) p = q
      k = k+1
   end do
   if (p .eq. zero) go to 135
   xl(i,nlc1) = one/p
   if (k .gt. nc) go to 20
   do j = k,nc
      a(i,j) = zero
   end do
20 i = i+1
   jbeg = jbeg-1
   if (jend-jbeg .eq. n) jend = jend-1
   if (i .le. nlc) go to 5
   jbeg = i
   nn = jend
25 jend = n-nuc
   do i = jbeg,n
      p = zero
      do j = 1,nn
         q = abs(a(i,j))
         if (q .gt. p) p = q
      end do
      if (p .eq. zero) go to 135
      xl(i,nlc1) = one/p
      if (i .eq. jend) go to 37
      if (i .lt. jend) cycle
      k = nn+1
      do j = k,nc
         a(i,j) = zero
      end do
37    nn = nn-1
   end do
   l = nlc
   do 75 k = 1,n
      p = abs(a(k,1))*xl(k,nlc1)
      i = k
      if (l .lt. n) l = l+1
      k1 = k+1
      if (k1 .gt. l) go to 50
      do j = k1,l
         q = abs(a(j,1))*xl(j,nlc1)
         if (q .le. p) cycle
         p = q
         i = j
      end do
50    xl(i,nlc1) = xl(k,nlc1)
      xl(k,nlc1) = i
      q = rn+p
      if (q .eq. rn) go to 135
      if (k .eq. i) go to 60
      do j = 1,nc
         p = a(k,j)
         a(k,j) = a(i,j)
         a(i,j) = p
      end do
60    if (k1 .gt. l) go to 75
      do 70 i = k1,l
         p = a(i,1)/a(k,1)
         ik = i-k
         xl(k1,ik) = p
         if(p.eq.0.0) then
            do j = 2,nc
               a(i,j-1) = a(i,j)
            end do
         else
            do j = 2,nc
               a(i,j-1) = a(i,j)-p*a(k,j)
            end do
         endif
         a(i,nc) = zero
70    continue
75 continue
   if (ijob .eq. 1) go to 9005
80 l = nlc
   do 105 k = 1,n
      i = xl(k,nlc1)
      if (i .eq. k) go to 90
      do j = 1,m
         p = b(k,j)
         b(k,j) = b(i,j)
         b(i,j) = p
      end do
90    if (l .lt. n) l = l+1
      k1 = k+1
      if (k1 .gt. l) go to 105
      do 100 i = k1,l
         ik = i-k
         p = xl(k1,ik)
         if(p.ne.0.0) then
            do j = 1,m
               b(i,j) = b(i,j)-p*b(k,j)
            end do
         endif
100   continue
105 continue
   jbeg = nuc+nlc
   do j = 1,m
      l = 1
      k1 = n+1
      do i = 1,n
         k = k1-i
         p = b(k,j)
         if (l .eq. 1) go to 115
         do kk = 2,l
            ik = kk+k
            p = p-a(k,kk)*b(ik-1,j)
         end do
115      b(k,j) = p/a(k,1)
         if (l .le. jbeg) l = l+1
      end do
   end do
   go to 9005
135 ier = 129
9000 continue
9005 return
end
