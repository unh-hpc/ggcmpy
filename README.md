# ggcmpy -- Python utilities for the OpenGGCM Global Magnetosphere Code

[![Actions Status][actions-badge]][actions-link]
[![Documentation Status][rtd-badge]][rtd-link]

[![PyPI version][pypi-version]][pypi-link]
[![PyPI platforms][pypi-platforms]][pypi-link]

<!-- [![Conda-Forge][conda-badge]][conda-link] -->

[![GitHub Discussion][github-discussions-badge]][github-discussions-link]

<!-- prettier-ignore-start -->
[actions-badge]:            https://github.com/unh-hpc/ggcmpy/workflows/CI/badge.svg
[actions-link]:             https://github.com/unh-hpc/ggcmpy/actions
[conda-badge]:              https://img.shields.io/conda/vn/conda-forge/ggcmpy
[conda-link]:               https://github.com/conda-forge/ggcmpy-feedstock
[github-discussions-badge]: https://img.shields.io/static/v1?label=Discussions&message=Ask&color=blue&logo=github
[github-discussions-link]:  https://github.com/unh-hpc/ggcmpy/discussions
[pypi-link]:                https://pypi.org/project/ggcmpy/
[pypi-platforms]:           https://img.shields.io/pypi/pyversions/ggcmpy
[pypi-version]:             https://img.shields.io/pypi/v/ggcmpy
[rtd-badge]:                https://readthedocs.org/projects/ggcmpy/badge/?version=latest
[rtd-link]:                 https://ggcmpy.readthedocs.io/en/latest/?badge=latest

<!-- prettier-ignore-end -->

[Documentation](https://ggcmpy.readthedocs.io/)

<!-- SPHINX-START -->

The ggcmpy package aims to be the place to collect various OpenGGCM related
tools.

- Support for reading OpenGGCM data files as Xarray datasets.

  Xarray is essentially an extension of numpy arrays, adding additional
  information like dimension names and coordinates, and making many data
  analysis tasks much easier. This package allows to read OpenGGCM binary files
  (`.iof`, `.3df`, `.p[xyz]_N`) simply by `ds = xr.open_dataset(filename)`.

- OpenGGCM-specific Xarray extensions

  - for now, this is limited to providing the `mlts` and `colats` coordinates in
    addition to the standard `longs`, `lats`.

- (TBD): OpenGGCM specific plotting support

- (TBD): Setting up an OpenGGCM run

  - Generating a runme
  - Generating a non-uniform grid
  - Preparing event solar wind data

- (external): Support for reading OpenGGCM XDMF/HDF5 data is available through
  [xarray-pschdf5](https://github.com/psc-code/xarray-pschdf5)

- (external): Support for reading OpenGGCM ADIOS2 data is available through
  [xarray-adios2](https://github.com/unh-hpc/xarray-adios2)

Needless to say, there is a lot of work the remains to be done, feedback /
requests and help are always appreciated!
