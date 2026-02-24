# Calculating Ground Magnetic Perturbations

OpenGGCM calculates the $\theta$-component (the southward component) of the
ground magnetic perturbation caused by substorm electrojets, $\delta B_\theta$,
as follows:

1. Inputs are the ionospheric electric potential, $\Phi$, Pedersen conductance,
   $\Sigma_P$, and Hall conductance, $\Sigma_H$.
2. Ionospheric electric field components perpendicular to the local magnetic
   field, $\left(E_\perp\right)_\phi$ and $\left(E_\perp\right)_\theta$, are
   calculated as the negative gradient of the potential,
   $\mathbf{E}_\perp=-\nabla_\perp\Phi$.
3. The total current densities are calculated by combining Pedersen and Hall
   current densities: $\mathbf{j}_\perp=\mathbf{j_P}+\mathbf{j_H}$.
4. Toroidal current density components perpendicular to the local magnetic
   field, $\left(j_{tor}\right)_\phi$ and $\left(j_{tor}\right)_\theta$,
   responsible for the ground magnetic perturbations according to the Fukushima
   theorem, are isolated from total current densities,
   $\left(j_\perp\right)_\phi$ and $\left(j_\perp\right)_\theta$.
5. A Biot-Savart integration is performed over the ionosphere to get
   $\delta B_\theta$.

The variables in the code are,
|Symbols|Variables|
|-------|---------|
|$\delta B_\theta$|`delbt`|
|$\Phi$|`pot`|
|$\Sigma_P$|`sigp`|
|$\Sigma_H$|`sigh`|
|$\left(E_\perp\right)_\phi$|`ep`|
|$\left(E_\perp\right)_\theta$|`et`|
|$\left(j_{pol}\right)_\phi$|`ctaup`|
|$\left(j_{pol}\right)_\theta$|`ctaut`|
|$\left(j_{tor}\right)_\phi$|`cpolp`|
|$\left(j_{tor}\right)_\theta$|`cpolt`|
|$\left(j_\perp\right)_\phi$|`ctiop`|
|$\left(j_\perp\right)_\theta$|`ctiot`|


## 1. Inputs

The process begins in the `iono_diag_post` subroutine in the `iono.F90` file:

```fortran
! -------------------------------------------------------------------
   subroutine iono_diag_post(mhd_tim)
! -------------------------------------------------------------------
!  ionosphere post-processing
      use constants_m, only: sig0, phi0, j0
      use ggcm_params_m
      use iono_m, only: pot, delphi, fac_tot, sigp, sigh, ppio, rrio, &
            hh1, hh2, hh3, hh4, hh6, awork
      real :: mhd_tim

      integer :: niowk

      call iono_diag_post_m_initialize()
      niowk = (size(awork,1)/2)-10

      ! Copy into temporary buffers.  This way we don't run into any problems
      ! with normalizing/unnormalizing.
      pot_tmp = pot * phi0
      delphi_tmp = delphi * phi0
      fac_tmp = fac_tot * j0
      sigp_tmp = sigp * sig0
      sigh_tmp = sigh * sig0

      call iopar(nphi,ntheta,fac_tmp,pot_tmp,sigp_tmp,sigh_tmp,delphi_tmp,rrio,ppio, &
            epio,etio,ctiot,ctiop,tau,ctaup,ctaut,cpolp,cpolt,delbr,delbp,delbt,xjh, &
            hh1,hh2,hh3,hh4,hh6,awork,niowk)

      ! write out
      call wrpot(ctiot, nphi, ntheta, 'ctiot ', 1.0, mhd_tim)
      call wrpot(ctiop, nphi, ntheta, 'ctiop ', 1.0, mhd_tim)
      call wrpot(tau  , nphi, ntheta, 'tau   ', 1.0, mhd_tim)
      call wrpot(ctaup, nphi, ntheta, 'ctaup ', 1.0, mhd_tim)
      call wrpot(ctaut, nphi, ntheta, 'ctaut ', 1.0, mhd_tim)
      call wrpot(cpolp, nphi, ntheta, 'cpolp ', 1.0, mhd_tim)
      call wrpot(cpolt, nphi, ntheta, 'cpolt ', 1.0, mhd_tim)
      call wrpot(delbr, nphi, ntheta, 'delbr ', 1.0, mhd_tim)
      call wrpot(delbp, nphi, ntheta, 'delbp ', 1.0, mhd_tim)
      call wrpot(delbt, nphi, ntheta, 'delbt ', 1.0, mhd_tim)
      call wrpot(epio , nphi, ntheta, 'epio  ', 1.0, mhd_tim)
      call wrpot(etio , nphi, ntheta, 'etio  ', 1.0, mhd_tim)
      call wrpot(xjh  , nphi, ntheta, 'xjh   ', 1.0, mhd_tim)
   end subroutine iono_diag_post
```

- The raw simulation data, such as `pot`, `sigp`, and `sigh`, are scaled by
  physical constants (e.g., `phi0`) and copied into temporary arrays (e.g.,
  `pot_temp`).
- The code calls the `iopar` subroutine in the `io-post.for` file, passing these
  scaled arrays. The output variable for the $\theta$-component of the magnetic
  perturbation is passed as an argument, `delbt`.

```fortran
c---------------------------------------------------------------------
      subroutine iopar(np,nt,fac,pot,sigp,sigh,delphi,rrio,ppio,ep,et,
     *ctiot,ctiop,tau,ctaup,ctaut,cpolp,cpolt,delbr,delbp,delbt,xjh,
     *hh1,hh2,hh3,f1,f2,wk,nniowk)
c---------------------------------------------------------------------
.for A in (fac,pot,sigp,sigh,delphi,rrio,ppio){real A(np,nt) }
.for A in (ctiot,ctiop,tau,ctaup,ctaut,cpolp,cpolt){real A(np,nt) }
.for A in (ep,et,delbr,delbp,delbt,xjh){real A(np,nt) }
.for A in (hh1,hh2,hh3){real A(np,nt) }
      real wk(nniowk),f1(*),f2(*)
c
c..... electric field
      call gradpt(pot,np,nt,et,ep)
c
      pi=4.0*atan(1.0)
      dp=2.0*pi/float(np-1)
      dt=pi/float(nt-1)
      re=6372.0e3
      f=1.0/re
c
      iv=0
      i=0
      do 100 it=1,nt
      t=float(it-1)*dt
      cost=cos(t)
      sint=sin(t)
      sini=amax1(0.05,2.0*abs(cost)/sqrt(1.0+3.0*cost*cost))
      do 100 ip=1,np
      p=float(ip-1)*dp
      i=i+1
      f1(i)=f*f
      f2(i)=0.0
c..... total ionospheric current
      stt=sigp(ip,it)/(sini*sini)
      spt=sigh(ip,it)/sini
      spp=sigp(ip,it)
      ctiot(ip,it)=stt*et(ip,it)+spt*ep(ip,it)
      ctiop(ip,it)=spp*ep(ip,it)-spt*et(ip,it)
c..... joule heat rate
      xjh(ip,it)=ctiot(ip,it)*et(ip,it)+ctiop(ip,it)*ep(ip,it)
      delbt(ip,it)=fac(ip,it)
100   continue
c
c...... solve poisson eq for the potential of the irrotational
c       (poloidal) current
      if(iv.ne.0)write(0,*)' solving poisson equation '
      niox=20
      nioy=10
      niogo=4
      niopy=2*nioy
      niowk=nniowk/2
c.... set up solver
      call io4(0,niox,nioy,niogo,niopy,f1,f2,delbt,tau,np,nt,wk,niowk,hh1,hh2,hh3)
c.... factor matrix
      call io4(1,niox,nioy,niogo,niopy,f1,f2,delbt,tau,np,nt,wk,niowk,hh1,hh2,hh3)
c.... solve matrix and interpolate
      call io4(2,niox,nioy,niogo,niopy,f1,f2,delbt,tau,np,nt,wk,niowk,hh1,hh2,hh3)
c
c...... irrotational (poloidal)  current
      call gradpt(tau,np,nt,ctaut,ctaup)
c
c...... toroidal (div free) part
      do 120 it=1,nt
      do 120 ip=1,np
      cpolt(ip,it)=ctiot(ip,it)-ctaut(ip,it)
      cpolp(ip,it)=ctiop(ip,it)-ctaup(ip,it)
120   continue
c
c..... ground magnetic perturbation
      call ground(np,nt,cpolp,cpolt,delbp,delbt,delbr)
c
      return
      end
```

## 2. Electric Field Calculation

Inside `iopar`, the first physics step is deriving the electric field from the
potential.

- The subroutine `gradpt` from the `io-post.for` file is called to calculate
  `et` and `ep` from `pot` using the finite difference method on the spherical
  grid.

```fortran
c----------------------------------------------------------------
      subroutine gradpt(pot,np,nt,et,ep)
c----------------------------------------------------------------
c.... calculate meridional and azimuthal electric field from
c     the potential, i.e.  (ep,et) = -grad pot
c
      real pot(np,nt),ep(np,nt),et(np,nt)
      pi=4.0*atan(1.0)
      dp=2.0*pi/float(np-1)
      dt=pi/float(nt-1)
      f=1.0/6372.0e3
      do 100 ip=1,np
      do 100 it=1,nt
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
100   continue
      return
      end
```

## 3. Total Ionospheric Current Density Calculation

Next, `iopar` calculates the total height-integrated ionospheric current
densities using Ohm's law.

- The code adjusts the conductances based on the magnetic inclination (dip)
  angle ($I$) to account for the geometry of the magnetic field lines,

$$
\sin I=\frac{2|\cos(\theta)|}{\sqrt{1+3\cos^2(\theta)}},
$$

- `ctiop` and `ctiot` and are computed by multiplying the `sigp` and `sigh` terms with the electric fields, i.e.,

$$
\begin{pmatrix}
\left(j_\perp\right)_\phi \\
\left(j_\perp\right)_\theta
\end{pmatrix}

=

\begin{pmatrix}
\Sigma_{\phi\phi} & -\Sigma_{\phi\theta} \\
\Sigma_{\theta\phi} & \Sigma_{\theta\theta}
\end{pmatrix}
\begin{pmatrix}
E_\phi \\
E_\theta
\end{pmatrix},
$$

where $\Sigma_{\phi\phi}\equiv\Sigma_P$,
$\Sigma_{\phi\theta}\equiv-\Sigma_H/\sin I$, and
$\Sigma_{\theta\theta}\equiv\Sigma_P/\sin^2 I$.

## 4. Current Density Decomposition

The next step in calculating ground magnetic perturbations is separating the
total current density into poloidal (curl-free) and toroidal (divergence-free)
parts. Ground magnetometers primarily detect the **toroidal part**.

- `iopar` solves a Poisson equation using the `io4` subroutine in the
  `io-psol.for` file to obtain the poloidal current density potential, `tau`,

$$
\nabla^2\tau=\nabla\cdot\mathbf{j_\perp}.
$$

- `ctaup` and `ctaut` are calculated from `tau` using `gradpt`.
- `cpolp` and `cpolt` are calculated by subtracting `ctaup` and `ctaut` from
  `ctiop` and `ctiot`, respectively.
  - `cpolp` = `ctiop` - `ctaup`
  - `cpolt` = `ctiot` - `ctaut`

## 5. Biot-Savart Integration
In theory, the ground magnetic perturbation in the $i$-direction at $\mathbf{r}$ due to a variation of a current density in the $j$-direction at $\mathbf{r}'$ is,

$$
\begin{align}
\delta B_i(\mathbf{r},\mathbf{r}') &= \frac{\mu_0}{4\pi}\frac{d^3 r'}{|\mathbf{r}-\mathbf{r}'|^3}\epsilon_{ijk}\delta j\left(\mathbf{r'}\right)^j \left(\mathbf{r}-\mathbf{r}'\right)^k, \\
&\equiv \frac{\mu_0}{4\pi}\frac{d^3 r'}{|\mathfrak{r}|^3}\epsilon_{ijk}\delta j^j\mathfrak{r}^k.
\end{align}
$$

Integrating over space, the total ground magnetic perturbation at $\mathbf{r}$ is then,
$$
\begin{equation}
\delta B_i(\mathbf{r})=\frac{\mu_0}{4\pi}\int\frac{d^3 r'}{|\mathfrak{r}|^3}\epsilon_{ijk}\delta j^j\mathfrak{r}^k.
\end{equation}
$$

In practice, `iopar` calls the `ground` subroutine in the `io-post.for` file to
calculate `delbt` using `cpolp` and `cpolt`.

```fortran
c----------------------------------------------------------------
      subroutine ground(np,nt,cpolp,cpolt,delbp,delbt,delbr)
c----------------------------------------------------------------
.for A in (cpolp,cpolt,delbp,delbt,delbr){real A(np,nt)}
.for A in (sp,cp,st,ct,st1,st2,ddf){real A(500)}
c
c      write(0,*)'start ground',np,nt
      pi=4.0*atan(1.0)
      xmue0=4.0e-7*pi
      dp=2.0*pi/float(np-1)
      dt=pi/float(nt-1)
      re=6372.0e3
      h=90.0e3    !  ionosphere height
      reh=re+h
      dd=600.0e3   !  radius to which to integrate
      d2=dd*dd
c
      do 400 ip=1,np
      p=float(ip-1)*dp
      sp(ip)=sin(p)
      cp(ip)=cos(p)
400   continue
      do 410 it=1,nt
      t=float(it-1)*dt
      st(it)=sin(t)
      ct(it)=cos(t)
      st1(it)=sin(t-0.5*dt)
      st2(it)=sin(t+0.5*dt)
      ddf(it)=0.5*re*re*dp*dt*(st1(it)+st2(it))*xmue0/(4.0*pi)
410   continue
c
      do 100 it=1,nt
      nt1=max0(2,it-14)
      nt2=min0(nt-1,it+14)
      t=float(it-1)*dt
      tdeg=t*180.0/pi
      do 100 ip=1,np
.for A in (r,p,t){delbA(ip,it)=0.0}
      if(tdeg.gt.45.0.and.tdeg.lt.135.0) goto 101
      p=float(ip-1)*dp
      x0=cp(ip)*st(it)*re
      y0=sp(ip)*st(it)*re
      z0=ct(it)*re
.for A in (x,y,z){dbA=0.0}
c
c..... integration loop
      do 110 jt=nt1,nt2
      do 110 jp=1,np-1
      x1=cp(jp)*st(jt)*reh
      y1=sp(jp)*st(jt)*reh
      z1=ct(jt)*reh
      dd=(x1-x0)**2+(y1-y0)**2+(z1-z0)**2
      if(dd.gt.d2) goto 111
.for A in (x,y,z){rA=A0-A1}
      r3i=1.0/(dd*sqrt(dd))
      df=ddf(jt)*r3i
      xjr=0.0
      xjp=cpolp(ip,it)
      xjt=cpolt(ip,it)
.call vec_pol_cart(xjr,xjp,xjt,xjx,xjy,xjz,sp(jp),cp(jp),st(jt),ct(jt))
.call cross11(xjx,xjy,xjz,rx,ry,rz,x,y,z)
.for A in (x,y,z){dbA=dbA+df*A}
111   continue
110   continue
c
.call vec_cart_pol(dbx,dby,dbz,dbr,dbp,dbt,sp(ip),cp(ip),st(it),ct(it))
.for A in (r,p,t){delbA(ip,it)=dbA}
c
101   continue
100   continue
c
      return
      end
```

The subroutine...

- Defines the ionosphere at a height `h = 90.0e3` (90.0 km)
- Defines an integration length `dd = 600e3` and area `d2 = dd*dd`
- Loops over every ground point, `ip` and `it`
- Loops over every source point in the ionosphere within the cutoff distance,
  `jp` and `jt`
- Computes the cross-product of the current density vector (`xjp`=`cpolp`,
  `xjt`=`cpolt`) and the position vector to get the magnetic field in Cartesian
  coordinates (`dbx`, `dby`, `dbz`)
- Transforms the Cartesian coordinates into the local spherical coordinates
  using the `vec_cart_pol` macro:

```fortran
.macro vec_cart_pol(AX,AY,AZ,AR,AP,AT,SP,CP,ST,CT)
      tmp_1=AX*CP+AY*SP
      AR=tmp_1*ST+AZ*CT
      AT=tmp_1*CT-AZ*ST
      AP=AY*CP-AX*SP
.endmacro
```

- Extracts the desired $\theta$-component via `AT = tmp_1*CT - AZ*ST`, i.e.,

$$
\delta B_\theta=\left(\delta B_x\cos{\phi}+\delta B_y\sin{\phi}\right)\cos{\theta}-\delta B_z\sin{\theta}
$$

- Accumulates the value into the array, `delbt(ip,it)`, i.e.,

$$
\boxed{
\delta B_\theta(\mathbf{r},\mathbf{r}')=\frac{\mu_0}{4\pi}\int\frac{d^3 r'}{|\mathfrak{r}|^3}\left[\left(\delta j^y\mathfrak{r}^z-\delta j^z\mathfrak{r}^y\right)\cos{\theta}\cos{\phi}+\left(\delta j^z\mathfrak{r}^x-\delta j^x\mathfrak{r}^z\right)\cos{\theta}\sin{\phi}+\left(\delta j^y\mathfrak{r}^x-\delta j^x\mathfrak{r}^y\right)\sin{\theta}\right]
}
$$
