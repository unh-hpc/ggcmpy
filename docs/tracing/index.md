# Particle Tracing

```{note}

    These docs aren't really even a good start yet :(.

    What's currently here are pieces that do the work, in particular various
    implementations of interpolating fields to particle positions, and
    implementations of the Boris pusher algorithm,
    which takes care of actually computing the particle trajectories.

    Eventually, there will be a higher level interface so one doesn't
    have to deal with implementation details.
```

Most functionality is currently more or less under development and implemented
in Python, where it's slow but easier to debug.

Some parts have been implemented in Fortran and interfaced into Python using
fwpy: Interpolation (cell-centered and Yee grid) and Boris integrator (on Yee
grid).

```{toctree}
:maxdepth: 1

simple-gyration
dipole
openggcm
```
