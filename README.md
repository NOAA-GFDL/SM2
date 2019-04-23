# GFDL FMS Slab Ocean Model

This repository includes the public release of the GFDL FMS Slab Ocean
Model.  The technical documentation of this model is available on
[GFDL's web
pages](https://www.gfdl.noaa.gov/fms-slab-ocean-model-technical-documentation/)

The layout of this package includes the following directories:

* src - The source code for the FMS Slab Ocean Model
* exec - Build directory with Makefiles for building the executable
* run - Sample run script

## Code Version

GFDL tracks code using the git version control system.  This package
includes a single version of the GFDL Slab Ocean model code.  The git
commit hash of the code used in this package is included in the table
below.

Component | Commit Hash
--------- | -----------
atmos_drivers | 3be6ed406de2db29766746a69115fd6a47048692
atmos_fv_dynamics | 75be30e6977582f011d0460e75cc7c59d56a3afc
atmos_param_am3 | b59c0450b5d1219a7d10c029354d16917884b26e
atmos_shared_am3 | 357aed1a76a7d7fec8cdcf80b62e9820974a7e45
coupler | 424267bf43adf01b6bee9d1f13befe807f9d18bc
ice_param | 154bd2b4bf523f3e699de5017679b156242ec13f
ice_sis | 0d15043e03bcaaca76bf41b76fbd2a383bd0e0b8
land_lad2 | a2a4cff2fab3ae040bd8ce24c31a4518fe4a49ae
ocean_slab | b7e4a54179d967074630e1f34ce7515f56f58cd0
fms | e8940fe90d68c3dc7c0d6bf1b8f552a577251754

While this repository contains only a single version of the code, some
of the source code components (e.g. [fms](src/fms)) are available on
the [NOAA-GFDL](https://github.com/NOAA-GFDL) github organization.

## Building the FMS Slab Ocean Model

The [exec](exec) directory contains Makefiles that can be used to
build the executable.  These Makefiles were generated using the [Make
Makefile (mkmf)](https://github.com/NOAA-GFDL/mkmf) program.  Included
in the exec directory is a sample make template file for the Intel
compilers ([intel.mk](exec/templates/intel.mk)).  This make template
can be used on any system with a relatively recent version of the
Intel compilers, the netCDF 4 library, and the MPICH2 MPI library.
Included in the [intel.mk](exec/templates/intel.mk) file are additional
settings that can be modified during the build.

To build the executable, modify the [intel.mk](exec/templates/intel.mk) to
match your system, and then run `make`.

## Obtaining the Input data

The input data required for running the FMS Slab Ocean model can be
obtained on the GFDL FTP site directory: ftp://ftp.gfdl.noaa.gov/perm/GFDL_pubrelease/SM2

The file
[SM2_Control_inputData.tar.gz](ftp://ftp.gfdl.noaa.gov/pub/perm/GFDL_pubrelease/SM2_Control_inputData.tar.gz)
contains the input data, and the runtime configuration required to run
the SM2 Control experiment.  The file
[SM2_Control_initCond_07610101.tar.gz](ftp://ftp.gfdl.noaa.gov/pub/perm/GFDL_pubrelease/SM2_Control_initCond_07610101.tar.gz)
contains the initial condition (restart) data to start the SM2 Control
experiment on January 01, 0761.

## Running the SM2_Control experiment

Included in the run directory is a sample run script.  To run the
SM2_Control experiment, first download the two data files mentioned in
the [Obtaining the Input data](obtaining-the-input-data) section.
Modify the variables in the configuration section in the sample run
script, and then run the script.

The sample data and run script are configured to run on 30 processors.
To run on a different number of processors, some settings in the
`input.nml` file (in the `SM2_Control_inputData.tar.gz` file) and run
script will need to be modified.  For example, to run the SM2 Control
experiment on 15 processors, do the following:

* Change the variable `total_npes` in the run script to `15`
* In `coupler_nml`
  * Change `atmos_npes` to `15`
  * Change `ocean_npes` to `15`
* In `fv_core_nml`, change `layout` to `1,15`
* In `land_model_nml`, change `layout` to `1,15`

Note: The `input.nml` file contains Fortran namelist and namelist
variables that modify, at run time, the model.  To learn more about
the settings in the `input.nml` file, please refer to source code
where the namelist/variable are defined.

## Disclaimer

The United States Department of Commerce (DOC) GitHub project code is
provided on an 'as is' basis and the user assumes responsibility for
its use. DOC has relinquished control of the information and no longer
has responsibility to protect the integrity, confidentiality, or
availability of the information. Any claims against the Department of
Commerce stemming from the use of its GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial
products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of
Commerce. The Department of Commerce seal and logo, or the seal and
logo of a DOC bureau, shall not be used in any manner to imply
endorsement of any commercial product or activity by DOC or the United
States Government.

This project code is made available through GitHub but is managed by
NOAA-GFDL at https://gitlab.gfdl.noaa.gov.
