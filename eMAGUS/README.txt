Here I (Neilen Marais) is describing the steps neccesary to get eMAGUS
built on Ubuntu 10.10.

Installs
========

gfortran: wajig install   gfortran-multilib gfortran-doc gfortran-4.4-multilib gfortran-4.4-doc libgfortran3-dbg lib32mudflap0

umfpack: wajig install libsuitesparse-dev

minpack: wajig install minpack-dev

arpack: wajig install libarpack2-dev

libatlas-base-dev: wajig install libatlas-base-dev

Code modifications
==================

- Shortened some lines by removing whitespace in geometry.f90
- removed -lg2c from linking; seems to have been a leftover from g77
- Changed FC ang F90 in makefile to use gfortran and added
  --ffixed-line-length-none option
- Had to modify format strings in output_mesh.f90. E.g. I -> I0, to
   use defaultformat

