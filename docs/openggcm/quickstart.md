# Installing and Running OpenGGCM

## System Requirements

OpenGGCM has fairly minimal dependencies:

**Required dependencies:**

- Fortran90 compiler (recommended: `gfortran`)
- MPI
- perl
- autotools

**Optional dependencies:**

- HDF5 (if using XDMF/HDF5 output)
- ADIOS2 (for code coupling)
- cmake

```{todo}
A devcontainer can be a convenient way of getting everything in place to build / run / develop OpenGGCM, but needs documentation and unfortunately is not usually a feasible way of running production.
```

### Marvin

On UNH's Cray "Marvin", I set up the environment for OpenGGCM using conda. (This
should work on other machines, too, if one opts to use conda. Most machines are
going to have a way to set up environments, it's generally probably best to
stick with that.)

```{todo}
Add how to installing (mini-)conda.
```

Create a `.condarc` in your home directory with the following content:

```sh
.condarc
channels:                                                                                                                      - conda-forge
 - defaults
channel_priority: strict
```

Create an environment for OpenGGCM:

```sh
conda create -n ggcm

conda install openmpi-mpifort
conda install "adios2=*=mpi_openmpi_*"
conda install cmake
conda install hdf5
conda install jupyter
conda install xarray
conda install dask
conda install matplotlib
pip install ggcmpy
```

You also want to set the environment variable `OPENGGCMDIR` to the location of
the OpenGGCM source you cloned from github above, e.g. by adding the following
to your `.bashrc` or `.zshrc`.

```sh
export OPENGGCMDIR=$HOME/src/openggcm # adjust for your OpenGGCM source location
```

Activate your environment.

```sh
[kaig@ln-0001 openggcm]$ conda activate ggcm
(ggcm) [kaig@ln-0001 openggcm]$
```

## Obtaining OpenGGCM

### Repository Access

OpenGGCM is maintained in a private github repository. Please contact Kai
<kai.germaschewski@unh.edu> or Doug <william.cramer@unh.edu> to be given access.

### Cloning the repository

```{todo}
* Setting up ssh or alternative github authentication
* forking the repo if planning to do development
```

```sh
(ggcm) [kaig@ln-0001 ~]$ mkdir src
(ggcm) [kaig@ln-0001 ~]$ cd src
(ggcm) [kaig@ln-0001 src]$ git clone git@github.com:unh-hpc/openggcm.git
Cloning into 'openggcm'...
remote: Enumerating objects: 70227, done.
remote: Counting objects: 100% (940/940), done.
remote: Compressing objects: 100% (316/316), done.
remote: Total 70227 (delta 682), reused 756 (delta 621), pack-reused 69287 (from 2)
Receiving objects: 100% (70227/70227), 86.08 MiB | 42.67 MiB/s, done.
Resolving deltas: 100% (54561/54561), done.
Updating files: 100% (1615/1615), done.
```

This should get you a copy of the OpenGGCM source in `src/openggcm`.

Have the autotools prepare the code for building:

```sh
(ggcm) [kaig@ln-0001 src]$ cd openggcm
(ggcm) [kaig@ln-0001 openggcm]$ ./autogen.sh
```

## Building and Installing OpenGGCM Utilities

Compile the utilities:

```sh
(ggcm) [kaig@ln-0001 openggcm]$ cd util
(ggcm) [kaig@ln-0001 util]$ ./configure --prefix=$HOME
[...]
(ggcm) [kaig@ln-0001 util]$ make
[...]
(ggcm) [kaig@ln-0001 util]$ make install
[...]
```

This compiles and install the utilities into `$PREFIX/bin`, so with the prefix
set to `$HOME` into the `bin/` folder in your home directory. Add that folder to
your path (again, something you may want to do in `.bashrc` or similar to
persist it).

```sh
(ggcm) [kaig@ln-0001 util]$ export PATH=$HOME/bin:$PATH
(ggcm) [kaig@ln-0001 util]$ fppn # make sure the fppn tool is now available
usage: fppn [-q] [-c] [-v1] [-v2] [-v3] [-v4] [-longout] [-fold] file
```

## Setting up a run

Let's set up a sample run, and run it.

```sh
(ggcm) [kaig@ln-0001 work]$ mkdir test0003
(ggcm) [kaig@ln-0001 work]$ cp -a $OPENGGCMDIR/run-template/test0003/{runme,swdata} test0003/
(ggcm) [kaig@ln-0001 work]$ cd test0003
(ggcm) [kaig@ln-0001 test0003]$ ./runme
[...]
```

Get a coffee and cross your fingers...

Once the build successfully finishes, one should make a batch script and submit
a job.

```{todo}
Really, the `runme` should be adapted to make a batch script for Marvin here.
```

But for now, we'll just run the code by hand.

```sh
(ggcm) [kaig@ln-0001 test0003]$ cd target
(ggcm) [kaig@ln-0001 target]$ srun -n 6 ./openggcm
srun: job 27370 queued and waiting for resources
```

Well, since the queue is full, I just ran on the head node:

```sh
(ggcm) [kaig@ln-0001 target]$ mpirun -n 6 ./openggcm
[...]
```

This run is set up to be a tiny test case, so it only runs for a minute or so.
But it should produce some MHD and ionosphere output files that can be plotted
using ggcmpy as shown in the
[sample notebook](../getting-started-guide/quick-overview.ipynb).
