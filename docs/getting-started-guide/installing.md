# Installation

To install the package from PyPI, you can use `pip`. Open your terminal and run
the following command:

```sh
pip install ggcmpy
```

Alternatively, you can do the same through

```sh
python -mpip install ggcmpy
```

This will download and install the latest version of the package from the Python
Package Index (PyPI).

ggcmpy uses an extension written in Fortran to read OpenGGCM jrlle custom binary
format, which means it contains binary code. We currently provide binary
"wheels" for Linux x86_64-based systems, but not yet for other architectures and
operating systems. If no binary wheel is available for your use case, the
command above will download the source and compile the code on your system --
however, this means that you'll need to having a working set of compilers
(gfortran) and libraries installed.

To verify the installation, you can run:

```sh
pip show ggcmpy
```

This will display information about the installed package.

## Required dependencies

Required dependencies will be automatically if you use `pip` to install ggcmpy
as shown above. Currently, `ggcmpy` has these required dependencies:

- xarray
- dask

The "dask" library helps to distribute array operations, but that should usually
happen under the hood and be visible to the user. At the current time, however,
`xarray.open_mfdataset()` needs dask to be installed to work.
