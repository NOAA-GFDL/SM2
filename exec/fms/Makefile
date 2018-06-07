# Makefile created by mkmf 19.3.0

CPPDEFS = -Duse_libMPI -Duse_netCDF -DINTERNAL_FILE_NML

include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libfms.a
MersenneTwister.o: $(SRCROOT)fms/random_numbers/MersenneTwister.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/random_numbers/MersenneTwister.F90
affinity.o: $(SRCROOT)fms/mpp/affinity.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)fms/mpp/affinity.c
amip_interp.o: $(SRCROOT)fms/amip_interp/amip_interp.F90 $(SRCROOT)fms/include/file_version.h time_interp.o time_manager.o get_cal_time.o mpp_io.o horiz_interp.o fms.o fms_io.o constants.o platform.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/amip_interp/amip_interp.F90
astronomy.o: $(SRCROOT)fms/astronomy/astronomy.F90 $(SRCROOT)fms/include/file_version.h fms.o time_manager.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/astronomy/astronomy.F90
atmos_ocean_fluxes.o: $(SRCROOT)fms/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)fms/include/file_version.h mpp.o fms.o coupler_types.o field_manager.o fm_util.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/coupler/atmos_ocean_fluxes.F90
axis_utils.o: $(SRCROOT)fms/axis_utils/axis_utils.F90 $(SRCROOT)fms/include/file_version.h mpp_io.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/axis_utils/axis_utils.F90
block_control.o: $(SRCROOT)fms/block_control/block_control.F90 $(SRCROOT)fms/include/fms_platform.h mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/block_control/block_control.F90
cloud_interpolator.o: $(SRCROOT)fms/drifters/cloud_interpolator.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/drifters/cloud_interpolator.F90
column_diagnostics.o: $(SRCROOT)fms/column_diagnostics/column_diagnostics.F90 $(SRCROOT)fms/include/file_version.h mpp_io.o fms.o time_manager.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/column_diagnostics/column_diagnostics.F90
constants.o: $(SRCROOT)fms/constants/constants.F90 $(SRCROOT)fms/include/file_version.h platform.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/constants/constants.F90
coupler_types.o: $(SRCROOT)fms/coupler/coupler_types.F90 $(SRCROOT)fms/include/file_version.h fms.o fms_io.o time_manager.o diag_manager.o data_override.o mpp_domains.o mpp.o fm_util.o field_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/coupler/coupler_types.F90
create_xgrid.o: $(SRCROOT)fms/mosaic/create_xgrid.c $(SRCROOT)fms/mosaic/mosaic_util.h $(SRCROOT)fms/mosaic/create_xgrid.h $(SRCROOT)fms/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)fms/mosaic	$(SRCROOT)fms/mosaic/create_xgrid.c
data_override.o: $(SRCROOT)fms/data_override/data_override.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h platform.o constants.o mpp_io.o mpp.o horiz_interp.o time_interp_external.o fms_io.o fms.o axis_utils.o mpp_domains.o time_manager.o diag_manager.o mpp_memutils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/data_override/data_override.F90
diag_axis.o: $(SRCROOT)fms/diag_manager/diag_axis.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h mpp_domains.o fms.o diag_data.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/diag_manager/diag_axis.F90
diag_data.o: $(SRCROOT)fms/diag_manager/diag_data.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h time_manager.o mpp_domains.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/diag_manager/diag_data.F90
diag_grid.o: $(SRCROOT)fms/diag_manager/diag_grid.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h constants.o fms.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/diag_manager/diag_grid.F90
diag_manager.o: $(SRCROOT)fms/diag_manager/diag_manager.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h time_manager.o mpp_io.o mpp.o fms.o fms_io.o diag_axis.o diag_util.o diag_data.o diag_table.o diag_output.o diag_grid.o diag_manifest.o constants.o mpp_domains.o mpp_parameter.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/diag_manager/diag_manager.F90
diag_manifest.o: $(SRCROOT)fms/diag_manager/diag_manifest.F90 diag_data.o mpp.o fms.o fms_io.o time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/diag_manager/diag_manifest.F90
diag_output.o: $(SRCROOT)fms/diag_manager/diag_output.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h mpp_io.o mpp_domains.o mpp.o diag_axis.o diag_data.o time_manager.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/diag_manager/diag_output.F90
diag_table.o: $(SRCROOT)fms/diag_manager/diag_table.F90 mpp_io.o mpp.o fms.o time_manager.o constants.o diag_data.o diag_util.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/diag_manager/diag_table.F90
diag_util.o: $(SRCROOT)fms/diag_manager/diag_util.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h diag_data.o diag_axis.o diag_output.o diag_grid.o fms.o fms_io.o mpp_domains.o time_manager.o mpp_io.o mpp.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/diag_manager/diag_util.F90
drifters.o: $(SRCROOT)fms/drifters/drifters.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/drifters/fms_switches.h $(SRCROOT)fms/include/file_version.h $(SRCROOT)fms/drifters/drifters_push.h $(SRCROOT)fms/drifters/drifters_set_field.h $(SRCROOT)fms/drifters/drifters_compute_k.h mpp.o mpp_domains.o drifters_core.o drifters_input.o drifters_io.o drifters_comm.o cloud_interpolator.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include -I$(SRCROOT)fms/drifters	$(SRCROOT)fms/drifters/drifters.F90
drifters_comm.o: $(SRCROOT)fms/drifters/drifters_comm.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/drifters/fms_switches.h mpp.o mpp_domains.o drifters_core.o drifters_input.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include -I$(SRCROOT)fms/drifters	$(SRCROOT)fms/drifters/drifters_comm.F90
drifters_core.o: $(SRCROOT)fms/drifters/drifters_core.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/drifters/drifters_core.F90
drifters_input.o: $(SRCROOT)fms/drifters/drifters_input.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/drifters/drifters_input.F90
drifters_io.o: $(SRCROOT)fms/drifters/drifters_io.F90 $(SRCROOT)fms/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/drifters/drifters_io.F90
ensemble_manager.o: $(SRCROOT)fms/coupler/ensemble_manager.F90 fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/coupler/ensemble_manager.F90
fft.o: $(SRCROOT)fms/fft/fft.F90 $(SRCROOT)fms/include/file_version.h platform.o fms.o fft99.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/fft/fft.F90
fft99.o: $(SRCROOT)fms/fft/fft99.F90 constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/fft/fft99.F90
field_manager.o: $(SRCROOT)fms/field_manager/field_manager.F90 $(SRCROOT)fms/include/file_version.h $(SRCROOT)fms/field_manager/parse.inc mpp.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include -I$(SRCROOT)fms/field_manager	$(SRCROOT)fms/field_manager/field_manager.F90
fm_util.o: $(SRCROOT)fms/field_manager/fm_util.F90 $(SRCROOT)fms/include/file_version.h field_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/field_manager/fm_util.F90
fms.o: $(SRCROOT)fms/fms/fms.F90 $(SRCROOT)fms/include/file_version.h mpp.o mpp_domains.o mpp_io.o fms_io.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/fms/fms.F90
fms_io.o: $(SRCROOT)fms/fms/fms_io.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h $(SRCROOT)fms/fms/read_data_2d.inc $(SRCROOT)fms/fms/read_data_3d.inc $(SRCROOT)fms/fms/read_data_4d.inc $(SRCROOT)fms/fms/write_data.inc $(SRCROOT)fms/fms/fms_io_unstructured_register_restart_axis.inc $(SRCROOT)fms/fms/fms_io_unstructured_setup_one_field.inc $(SRCROOT)fms/fms/fms_io_unstructured_register_restart_field.inc $(SRCROOT)fms/fms/fms_io_unstructured_save_restart.inc $(SRCROOT)fms/fms/fms_io_unstructured_read.inc $(SRCROOT)fms/fms/fms_io_unstructured_get_file_name.inc $(SRCROOT)fms/fms/fms_io_unstructured_get_file_unit.inc $(SRCROOT)fms/fms/fms_io_unstructured_file_unit.inc $(SRCROOT)fms/fms/fms_io_unstructured_get_field_size.inc $(SRCROOT)fms/fms/fms_io_unstructured_field_exist.inc mpp_io.o mpp_domains.o mpp.o platform.o mpp_parameter.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include -I$(SRCROOT)fms/fms	$(SRCROOT)fms/fms/fms_io.F90
gaussian_topog.o: $(SRCROOT)fms/topography/gaussian_topog.F90 $(SRCROOT)fms/include/file_version.h fms.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/topography/gaussian_topog.F90
get_cal_time.o: $(SRCROOT)fms/time_manager/get_cal_time.F90 $(SRCROOT)fms/include/file_version.h fms.o time_manager.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/time_manager/get_cal_time.F90
gradient.o: $(SRCROOT)fms/mosaic/gradient.F90 $(SRCROOT)fms/include/file_version.h mpp.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/mosaic/gradient.F90
gradient_c2l.o: $(SRCROOT)fms/mosaic/gradient_c2l.c $(SRCROOT)fms/mosaic/constant.h $(SRCROOT)fms/mosaic/mosaic_util.h $(SRCROOT)fms/mosaic/gradient_c2l.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)fms/mosaic	$(SRCROOT)fms/mosaic/gradient_c2l.c
grid.o: $(SRCROOT)fms/mosaic/grid.F90 $(SRCROOT)fms/include/file_version.h mpp.o constants.o fms.o fms_io.o mosaic.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/mosaic/grid.F90
horiz_interp.o: $(SRCROOT)fms/horiz_interp/horiz_interp.F90 $(SRCROOT)fms/include/file_version.h fms.o mpp.o constants.o horiz_interp_type.o horiz_interp_conserve.o horiz_interp_bilinear.o horiz_interp_bicubic.o horiz_interp_spherical.o mpp_io.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/horiz_interp/horiz_interp.F90
horiz_interp_bicubic.o: $(SRCROOT)fms/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)fms/include/file_version.h mpp.o fms.o horiz_interp_type.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/horiz_interp/horiz_interp_bicubic.F90
horiz_interp_bilinear.o: $(SRCROOT)fms/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)fms/include/file_version.h mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/horiz_interp/horiz_interp_bilinear.F90
horiz_interp_conserve.o: $(SRCROOT)fms/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)fms/include/file_version.h mpp.o fms.o fms_io.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/horiz_interp/horiz_interp_conserve.F90
horiz_interp_spherical.o: $(SRCROOT)fms/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)fms/include/file_version.h mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/horiz_interp/horiz_interp_spherical.F90
horiz_interp_type.o: $(SRCROOT)fms/horiz_interp/horiz_interp_type.F90 mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/horiz_interp/horiz_interp_type.F90
interp.o: $(SRCROOT)fms/mosaic/interp.c $(SRCROOT)fms/mosaic/mosaic_util.h $(SRCROOT)fms/mosaic/interp.h $(SRCROOT)fms/mosaic/create_xgrid.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)fms/mosaic	$(SRCROOT)fms/mosaic/interp.c
interpolator.o: $(SRCROOT)fms/interpolator/interpolator.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h mpp.o mpp_io.o mpp_domains.o diag_manager.o fms.o horiz_interp.o time_manager.o time_interp.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/interpolator/interpolator.F90
memuse.o: $(SRCROOT)fms/memutils/memuse.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)fms/memutils/memuse.c
memutils.o: $(SRCROOT)fms/memutils/memutils.F90 mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/memutils/memutils.F90
mosaic.o: $(SRCROOT)fms/mosaic/mosaic.F90 $(SRCROOT)fms/include/file_version.h fms.o mpp.o mpp_io.o fms_io.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/mosaic/mosaic.F90
mosaic_util.o: $(SRCROOT)fms/mosaic/mosaic_util.c $(SRCROOT)fms/mosaic/mosaic_util.h $(SRCROOT)fms/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)fms/mosaic	$(SRCROOT)fms/mosaic/mosaic_util.c
mpp.o: $(SRCROOT)fms/mpp/mpp.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h $(SRCROOT)fms/mpp/include/system_clock.h $(SRCROOT)fms/mpp/include/mpp_util.inc $(SRCROOT)fms/mpp/include/mpp_util_sma.inc $(SRCROOT)fms/mpp/include/mpp_util_mpi.inc $(SRCROOT)fms/mpp/include/mpp_util_nocomm.inc $(SRCROOT)fms/mpp/include/mpp_error_a_a.h $(SRCROOT)fms/mpp/include/mpp_error_a_s.h $(SRCROOT)fms/mpp/include/mpp_error_s_a.h $(SRCROOT)fms/mpp/include/mpp_error_s_s.h $(SRCROOT)fms/mpp/include/mpp_comm.inc $(SRCROOT)fms/mpp/include/mpp_comm_sma.inc $(SRCROOT)fms/mpp/include/mpp_transmit_sma.h $(SRCROOT)fms/mpp/include/mpp_transmit.inc $(SRCROOT)fms/mpp/include/mpp_reduce_sma.h $(SRCROOT)fms/mpp/include/mpp_sum_sma.h $(SRCROOT)fms/mpp/include/mpp_sum.inc $(SRCROOT)fms/mpp/include/mpp_alltoall_sma.h $(SRCROOT)fms/mpp/include/mpp_comm_mpi.inc $(SRCROOT)fms/mpp/include/mpp_transmit_mpi.h $(SRCROOT)fms/mpp/include/mpp_reduce_mpi.h $(SRCROOT)fms/mpp/include/mpp_sum_mpi.h $(SRCROOT)fms/mpp/include/mpp_alltoall_mpi.h $(SRCROOT)fms/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)fms/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)fms/mpp/include/mpp_reduce_nocomm.h $(SRCROOT)fms/mpp/include/mpp_sum_nocomm.h $(SRCROOT)fms/mpp/include/mpp_alltoall_nocomm.h $(SRCROOT)fms/mpp/include/mpp_chksum_int.h $(SRCROOT)fms/mpp/include/mpp_chksum_scalar.h $(SRCROOT)fms/mpp/include/mpp_chksum.h $(SRCROOT)fms/mpp/include/mpp_gather.h $(SRCROOT)fms/mpp/include/mpp_scatter.h mpp_parameter.o mpp_data.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include -I$(SRCROOT)fms/mpp/include	$(SRCROOT)fms/mpp/mpp.F90
mpp_data.o: $(SRCROOT)fms/mpp/mpp_data.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h $(SRCROOT)fms/mpp/include/mpp_data_sma.inc $(SRCROOT)fms/mpp/include/mpp_data_mpi.inc $(SRCROOT)fms/mpp/include/mpp_data_nocomm.inc mpp_parameter.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include -I$(SRCROOT)fms/mpp/include	$(SRCROOT)fms/mpp/mpp_data.F90
mpp_domains.o: $(SRCROOT)fms/mpp/mpp_domains.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h $(SRCROOT)fms/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)fms/mpp/include/mpp_domains_util.inc $(SRCROOT)fms/mpp/include/mpp_domains_comm.inc $(SRCROOT)fms/mpp/include/mpp_domains_define.inc $(SRCROOT)fms/mpp/include/mpp_domains_misc.inc $(SRCROOT)fms/mpp/include/mpp_update_domains2D.h $(SRCROOT)fms/mpp/include/mpp_update_domains2D_nonblock.h $(SRCROOT)fms/mpp/include/mpp_do_update_nonblock.h $(SRCROOT)fms/mpp/include/mpp_do_updateV_nonblock.h $(SRCROOT)fms/mpp/include/mpp_do_update.h $(SRCROOT)fms/mpp/include/mpp_do_updateV.h $(SRCROOT)fms/mpp/include/mpp_do_check.h $(SRCROOT)fms/mpp/include/mpp_do_checkV.h $(SRCROOT)fms/mpp/include/mpp_update_nest_domains.h $(SRCROOT)fms/mpp/include/mpp_do_update_nest.h $(SRCROOT)fms/mpp/include/mpp_do_redistribute.h $(SRCROOT)fms/mpp/include/mpp_get_boundary.h $(SRCROOT)fms/mpp/include/mpp_do_get_boundary.h $(SRCROOT)fms/mpp/include/mpp_group_update.h $(SRCROOT)fms/mpp/include/group_update_pack.inc $(SRCROOT)fms/mpp/include/group_update_unpack.inc $(SRCROOT)fms/mpp/include/mpp_domains_reduce.inc $(SRCROOT)fms/mpp/include/mpp_global_reduce.h $(SRCROOT)fms/mpp/include/mpp_global_sum.h $(SRCROOT)fms/mpp/include/mpp_global_sum_tl.h $(SRCROOT)fms/mpp/include/mpp_global_field.h $(SRCROOT)fms/mpp/include/mpp_do_global_field.h $(SRCROOT)fms/mpp/include/mpp_unstruct_domain.inc $(SRCROOT)fms/mpp/include/mpp_unstruct_pass_data.h $(SRCROOT)fms/mpp/include/mpp_global_field_ug.h mpp_parameter.o mpp_data.o mpp.o mpp_memutils.o mpp_pset.o mpp_efp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include -I$(SRCROOT)fms/mpp/include	$(SRCROOT)fms/mpp/mpp_domains.F90
mpp_efp.o: $(SRCROOT)fms/mpp/mpp_efp.F90 $(SRCROOT)fms/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/mpp/mpp_efp.F90
mpp_io.o: $(SRCROOT)fms/mpp/mpp_io.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h $(SRCROOT)fms/mpp/include/mpp_io_util.inc $(SRCROOT)fms/mpp/include/mpp_io_misc.inc $(SRCROOT)fms/mpp/include/mpp_io_connect.inc $(SRCROOT)fms/mpp/include/mpp_io_read.inc $(SRCROOT)fms/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)fms/mpp/include/mpp_read_compressed.h $(SRCROOT)fms/mpp/include/mpp_read_distributed_ascii.inc $(SRCROOT)fms/mpp/include/mpp_read_distributed_ascii.h $(SRCROOT)fms/mpp/include/mpp_io_write.inc $(SRCROOT)fms/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)fms/mpp/include/mpp_write_compressed.h $(SRCROOT)fms/mpp/include/mpp_write_unlimited_axis.h $(SRCROOT)fms/mpp/include/mpp_write.h $(SRCROOT)fms/mpp/include/mpp_io_unstructured_write.inc $(SRCROOT)fms/mpp/include/mpp_io_unstructured_read.inc mpp_parameter.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include -I$(SRCROOT)fms/mpp/include	$(SRCROOT)fms/mpp/mpp_io.F90
mpp_memutils.o: $(SRCROOT)fms/mpp/mpp_memutils.F90 mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/mpp/mpp_memutils.F90
mpp_parameter.o: $(SRCROOT)fms/mpp/mpp_parameter.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/mpp/mpp_parameter.F90
mpp_pset.o: $(SRCROOT)fms/mpp/mpp_pset.F90 $(SRCROOT)fms/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/mpp/mpp_pset.F90
mpp_utilities.o: $(SRCROOT)fms/mpp/mpp_utilities.F90 $(SRCROOT)fms/include/file_version.h mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/mpp/mpp_utilities.F90
nsclock.o: $(SRCROOT)fms/mpp/nsclock.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)fms/mpp/nsclock.c
oda_core.o: $(SRCROOT)fms/oda_tools/oda_core.F90 fms.o mpp.o mpp_domains.o time_manager.o get_cal_time.o axis_utils.o constants.o oda_types.o write_ocean_data.o mpp_io.o field_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/oda_tools/oda_core.F90
oda_core_ecda.o: $(SRCROOT)fms/oda_tools/oda_core_ecda.F90 fms.o mpp.o mpp_io.o mpp_domains.o mpp_memutils.o time_manager.o get_cal_time.o axis_utils.o horiz_interp_type.o horiz_interp_bilinear.o constants.o oda_types.o xbt_drop_rate_adjust.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/oda_tools/oda_core_ecda.F90
oda_types.o: $(SRCROOT)fms/oda_tools/oda_types.F90 time_manager.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/oda_tools/oda_types.F90
platform.o: $(SRCROOT)fms/platform/platform.F90 $(SRCROOT)fms/include/fms_platform.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/platform/platform.F90
quicksort.o: $(SRCROOT)fms/drifters/quicksort.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/drifters/quicksort.F90
random_numbers.o: $(SRCROOT)fms/random_numbers/random_numbers.F90 MersenneTwister.o time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/random_numbers/random_numbers.F90
read_mosaic.o: $(SRCROOT)fms/mosaic/read_mosaic.c $(SRCROOT)fms/mosaic/read_mosaic.h $(SRCROOT)fms/mosaic/constant.h $(SRCROOT)fms/mosaic/mosaic_util.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)fms/mosaic	$(SRCROOT)fms/mosaic/read_mosaic.c
sat_vapor_pres.o: $(SRCROOT)fms/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)fms/include/file_version.h constants.o fms.o mpp_io.o mpp.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/sat_vapor_pres/sat_vapor_pres.F90
sat_vapor_pres_k.o: $(SRCROOT)fms/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)fms/include/file_version.h
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/sat_vapor_pres/sat_vapor_pres_k.F90
station_data.o: $(SRCROOT)fms/station_data/station_data.F90 $(SRCROOT)fms/include/file_version.h axis_utils.o mpp_io.o fms.o mpp.o mpp_domains.o diag_axis.o diag_output.o diag_manager.o diag_util.o time_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/station_data/station_data.F90
stock_constants.o: $(SRCROOT)fms/exchange/stock_constants.F90 $(SRCROOT)fms/include/file_version.h mpp.o fms.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/exchange/stock_constants.F90
test_fms_io.o: $(SRCROOT)fms/fms/test_fms_io.F90 $(SRCROOT)fms/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o fms_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/fms/test_fms_io.F90
test_horiz_interp.o: $(SRCROOT)fms/horiz_interp/test_horiz_interp.F90 mpp.o mpp_io.o mpp_domains.o mosaic.o horiz_interp.o axis_utils.o fms_io.o fms.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/horiz_interp/test_horiz_interp.F90
test_mpp.o: $(SRCROOT)fms/mpp/test_mpp.F90 $(SRCROOT)fms/include/fms_platform.h mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/mpp/test_mpp.F90
test_mpp_domains.o: $(SRCROOT)fms/mpp/test_mpp_domains.F90 $(SRCROOT)fms/include/fms_platform.h mpp.o mpp_domains.o mpp_memutils.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/mpp/test_mpp_domains.F90
test_mpp_io.o: $(SRCROOT)fms/mpp/test_mpp_io.F90 $(SRCROOT)fms/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/mpp/test_mpp_io.F90
test_mpp_pset.o: $(SRCROOT)fms/mpp/test_mpp_pset.F90 $(SRCROOT)fms/include/fms_platform.h mpp.o mpp_pset.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/mpp/test_mpp_pset.F90
test_unstructured_fms_io.o: $(SRCROOT)fms/fms/test_unstructured_fms_io.F90 $(SRCROOT)fms/include/fms_platform.h mpp_parameter.o mpp.o mpp_domains.o mpp_io.o fms_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/fms/test_unstructured_fms_io.F90
test_xgrid.o: $(SRCROOT)fms/exchange/test_xgrid.F90 $(SRCROOT)fms/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o fms.o fms_io.o constants.o xgrid.o mosaic.o grid.o gradient.o ensemble_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/exchange/test_xgrid.F90
threadloc.o: $(SRCROOT)fms/mpp/threadloc.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)fms/mpp/threadloc.c
time_interp.o: $(SRCROOT)fms/time_interp/time_interp.F90 $(SRCROOT)fms/include/file_version.h time_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/time_interp/time_interp.F90
time_interp_external.o: $(SRCROOT)fms/time_interp/time_interp_external.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h fms.o mpp.o mpp_io.o time_manager.o get_cal_time.o mpp_domains.o time_interp.o axis_utils.o platform.o horiz_interp.o constants.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/time_interp/time_interp_external.F90
time_manager.o: $(SRCROOT)fms/time_manager/time_manager.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h constants.o fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/time_manager/time_manager.F90
topography.o: $(SRCROOT)fms/topography/topography.F90 $(SRCROOT)fms/include/file_version.h gaussian_topog.o horiz_interp.o fms.o fms_io.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/topography/topography.F90
tracer_manager.o: $(SRCROOT)fms/tracer_manager/tracer_manager.F90 $(SRCROOT)fms/include/file_version.h mpp.o mpp_io.o fms.o field_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/tracer_manager/tracer_manager.F90
tridiagonal.o: $(SRCROOT)fms/tridiagonal/tridiagonal.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/tridiagonal/tridiagonal.F90
write_ocean_data.o: $(SRCROOT)fms/oda_tools/write_ocean_data.F90 mpp_io.o mpp.o oda_types.o time_manager.o mpp_domains.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/oda_tools/write_ocean_data.F90
xbt_drop_rate_adjust.o: $(SRCROOT)fms/oda_tools/xbt_drop_rate_adjust.f90 oda_types.o
	$(FC) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)fms/oda_tools/xbt_drop_rate_adjust.f90
xgrid.o: $(SRCROOT)fms/exchange/xgrid.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/include/file_version.h fms.o fms_io.o mpp.o mpp_domains.o mpp_io.o constants.o mosaic.o stock_constants.o gradient.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include	$(SRCROOT)fms/exchange/xgrid.F90

SRC = $(SRCROOT)fms/random_numbers/random_numbers.F90 $(SRCROOT)fms/random_numbers/MersenneTwister.F90 $(SRCROOT)fms/exchange/xgrid.F90 $(SRCROOT)fms/oda_tools/write_ocean_data.F90 $(SRCROOT)fms/station_data/station_data.F90 $(SRCROOT)fms/axis_utils/axis_utils.F90 $(SRCROOT)fms/mosaic/grid.F90 $(SRCROOT)fms/mpp/nsclock.c $(SRCROOT)fms/horiz_interp/horiz_interp_type.F90 $(SRCROOT)fms/mpp/mpp_efp.F90 $(SRCROOT)fms/constants/constants.F90 $(SRCROOT)fms/time_manager/get_cal_time.F90 $(SRCROOT)fms/fft/fft.F90 $(SRCROOT)fms/field_manager/fm_util.F90 $(SRCROOT)fms/topography/gaussian_topog.F90 $(SRCROOT)fms/topography/topography.F90 $(SRCROOT)fms/horiz_interp/horiz_interp.F90 $(SRCROOT)fms/time_interp/time_interp.F90 $(SRCROOT)fms/mpp/test_mpp_io.F90 $(SRCROOT)fms/fms/fms.F90 $(SRCROOT)fms/diag_manager/diag_grid.F90 $(SRCROOT)fms/drifters/drifters_io.F90 $(SRCROOT)fms/fms/fms_io.F90 $(SRCROOT)fms/mpp/mpp_io.F90 $(SRCROOT)fms/oda_tools/oda_core_ecda.F90 $(SRCROOT)fms/field_manager/field_manager.F90 $(SRCROOT)fms/column_diagnostics/column_diagnostics.F90 $(SRCROOT)fms/mosaic/gradient.F90 $(SRCROOT)fms/interpolator/interpolator.F90 $(SRCROOT)fms/mpp/test_mpp_domains.F90 $(SRCROOT)fms/time_manager/time_manager.F90 $(SRCROOT)fms/diag_manager/diag_data.F90 $(SRCROOT)fms/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)fms/drifters/quicksort.F90 $(SRCROOT)fms/diag_manager/diag_table.F90 $(SRCROOT)fms/mosaic/mosaic.F90 $(SRCROOT)fms/astronomy/astronomy.F90 $(SRCROOT)fms/diag_manager/diag_manager.F90 $(SRCROOT)fms/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)fms/mosaic/create_xgrid.c $(SRCROOT)fms/mpp/affinity.c $(SRCROOT)fms/time_interp/time_interp_external.F90 $(SRCROOT)fms/data_override/data_override.F90 $(SRCROOT)fms/amip_interp/amip_interp.F90 $(SRCROOT)fms/exchange/test_xgrid.F90 $(SRCROOT)fms/drifters/drifters.F90 $(SRCROOT)fms/mpp/test_mpp_pset.F90 $(SRCROOT)fms/mpp/mpp_memutils.F90 $(SRCROOT)fms/mpp/mpp.F90 $(SRCROOT)fms/fms/test_fms_io.F90 $(SRCROOT)fms/diag_manager/diag_manifest.F90 $(SRCROOT)fms/drifters/cloud_interpolator.F90 $(SRCROOT)fms/mpp/threadloc.c $(SRCROOT)fms/coupler/coupler_types.F90 $(SRCROOT)fms/mpp/mpp_utilities.F90 $(SRCROOT)fms/mpp/mpp_parameter.F90 $(SRCROOT)fms/oda_tools/oda_core.F90 $(SRCROOT)fms/mosaic/read_mosaic.c $(SRCROOT)fms/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)fms/mpp/test_mpp.F90 $(SRCROOT)fms/oda_tools/xbt_drop_rate_adjust.f90 $(SRCROOT)fms/coupler/ensemble_manager.F90 $(SRCROOT)fms/fms/test_unstructured_fms_io.F90 $(SRCROOT)fms/drifters/drifters_input.F90 $(SRCROOT)fms/diag_manager/diag_output.F90 $(SRCROOT)fms/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)fms/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)fms/horiz_interp/test_horiz_interp.F90 $(SRCROOT)fms/mosaic/gradient_c2l.c $(SRCROOT)fms/tracer_manager/tracer_manager.F90 $(SRCROOT)fms/memutils/memutils.F90 $(SRCROOT)fms/block_control/block_control.F90 $(SRCROOT)fms/mpp/mpp_domains.F90 $(SRCROOT)fms/mpp/mpp_pset.F90 $(SRCROOT)fms/mpp/mpp_data.F90 $(SRCROOT)fms/drifters/drifters_comm.F90 $(SRCROOT)fms/platform/platform.F90 $(SRCROOT)fms/exchange/stock_constants.F90 $(SRCROOT)fms/fft/fft99.F90 $(SRCROOT)fms/diag_manager/diag_axis.F90 $(SRCROOT)fms/oda_tools/oda_types.F90 $(SRCROOT)fms/tridiagonal/tridiagonal.F90 $(SRCROOT)fms/drifters/drifters_core.F90 $(SRCROOT)fms/mosaic/interp.c $(SRCROOT)fms/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)fms/diag_manager/diag_util.F90 $(SRCROOT)fms/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)fms/memutils/memuse.c $(SRCROOT)fms/mosaic/mosaic_util.c $(SRCROOT)fms/mosaic/interp.h $(SRCROOT)fms/drifters/drifters_push.h $(SRCROOT)fms/fms/fms_io_unstructured_save_restart.inc $(SRCROOT)fms/mpp/include/mpp_transmit_sma.h $(SRCROOT)fms/drifters/fms_switches.h $(SRCROOT)fms/mpp/include/mpp_do_checkV.h $(SRCROOT)fms/mpp/include/mpp_write_compressed.h $(SRCROOT)fms/mpp/include/mpp_data_sma.inc $(SRCROOT)fms/fms/fms_io_unstructured_get_file_unit.inc $(SRCROOT)fms/mpp/include/mpp_read_distributed_ascii.inc $(SRCROOT)fms/mpp/include/mpp_alltoall_nocomm.h $(SRCROOT)fms/mosaic/constant.h $(SRCROOT)fms/mpp/include/mpp_do_check.h $(SRCROOT)fms/mpp/include/mpp_do_update_nonblock.h $(SRCROOT)fms/drifters/drifters_compute_k.h $(SRCROOT)fms/mosaic/read_mosaic.h $(SRCROOT)fms/fms/fms_io_unstructured_register_restart_field.inc $(SRCROOT)fms/mpp/include/mpp_reduce_nocomm.h $(SRCROOT)fms/mpp/include/mpp_domains_comm.inc $(SRCROOT)fms/mpp/include/mpp_get_boundary.h $(SRCROOT)fms/mpp/include/group_update_pack.inc $(SRCROOT)fms/mpp/include/mpp_sum_sma.h $(SRCROOT)fms/mpp/include/mpp_transmit_mpi.h $(SRCROOT)fms/mpp/include/mpp_reduce_mpi.h $(SRCROOT)fms/mpp/include/mpp_sum_nocomm.h $(SRCROOT)fms/drifters/drifters_set_field.h $(SRCROOT)fms/mpp/include/mpp_error_a_s.h $(SRCROOT)fms/mpp/include/mpp_alltoall_mpi.h $(SRCROOT)fms/mpp/include/mpp_util_sma.inc $(SRCROOT)fms/mpp/include/mpp_domains_reduce.inc $(SRCROOT)fms/mpp/include/mpp_io_read.inc $(SRCROOT)fms/mpp/include/mpp_chksum.h $(SRCROOT)fms/mpp/include/mpp_io_unstructured_read.inc $(SRCROOT)fms/field_manager/parse.inc $(SRCROOT)fms/mpp/include/mpp_comm_sma.inc $(SRCROOT)fms/mpp/include/mpp_comm.inc $(SRCROOT)fms/fms/read_data_2d.inc $(SRCROOT)fms/mpp/include/mpp_io_write.inc $(SRCROOT)fms/mpp/include/system_clock.h $(SRCROOT)fms/fms/fms_io_unstructured_field_exist.inc $(SRCROOT)fms/mpp/include/mpp_global_reduce.h $(SRCROOT)fms/mpp/include/mpp_global_field.h $(SRCROOT)fms/fms/fms_io_unstructured_setup_one_field.inc $(SRCROOT)fms/mpp/include/mpp_do_updateV_nonblock.h $(SRCROOT)fms/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)fms/mpp/include/mpp_do_global_field.h $(SRCROOT)fms/mpp/include/mpp_do_update.h $(SRCROOT)fms/mpp/include/mpp_update_domains2D_nonblock.h $(SRCROOT)fms/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)fms/fms/read_data_3d.inc $(SRCROOT)fms/mpp/include/mpp_global_sum.h $(SRCROOT)fms/include/file_version.h $(SRCROOT)fms/mpp/include/mpp_global_field_ug.h $(SRCROOT)fms/mpp/include/mpp_global_sum_tl.h $(SRCROOT)fms/mpp/include/mpp_domains_define.inc $(SRCROOT)fms/mpp/include/mpp_util_nocomm.inc $(SRCROOT)fms/mpp/include/mpp_write_unlimited_axis.h $(SRCROOT)fms/mpp/include/mpp_util.inc $(SRCROOT)fms/mpp/include/mpp_io_util.inc $(SRCROOT)fms/mpp/include/mpp_do_updateV.h $(SRCROOT)fms/mpp/include/mpp_group_update.h $(SRCROOT)fms/mpp/include/mpp_data_nocomm.inc $(SRCROOT)fms/mosaic/gradient_c2l.h $(SRCROOT)fms/mpp/include/mpp_domains_misc.inc $(SRCROOT)fms/mpp/include/mpp_sum_mpi.h $(SRCROOT)fms/mpp/include/mpp_error_s_s.h $(SRCROOT)fms/mpp/include/mpp_reduce_sma.h $(SRCROOT)fms/mpp/include/group_update_unpack.inc $(SRCROOT)fms/mpp/include/mpp_gather.h $(SRCROOT)fms/mpp/include/mpp_update_nest_domains.h $(SRCROOT)fms/mpp/include/mpp_domains_util.inc $(SRCROOT)fms/mpp/include/mpp_read_distributed_ascii.h $(SRCROOT)fms/mpp/include/mpp_comm_mpi.inc $(SRCROOT)fms/mpp/include/mpp_unstruct_pass_data.h $(SRCROOT)fms/mpp/include/mpp_do_get_boundary.h $(SRCROOT)fms/mpp/include/mpp_error_s_a.h $(SRCROOT)fms/mpp/include/mpp_util_mpi.inc $(SRCROOT)fms/mpp/include/mpp_data_mpi.inc $(SRCROOT)fms/mpp/include/mpp_sum.inc $(SRCROOT)fms/fms/write_data.inc $(SRCROOT)fms/mpp/include/mpp_error_a_a.h $(SRCROOT)fms/mpp/include/mpp_chksum_int.h $(SRCROOT)fms/fms/fms_io_unstructured_file_unit.inc $(SRCROOT)fms/fms/fms_io_unstructured_get_file_name.inc $(SRCROOT)fms/mpp/include/mpp_scatter.h $(SRCROOT)fms/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)fms/mpp/include/mpp_do_update_nest.h $(SRCROOT)fms/fms/read_data_4d.inc $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)fms/mpp/include/mpp_io_connect.inc $(SRCROOT)fms/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)fms/fms/fms_io_unstructured_register_restart_axis.inc $(SRCROOT)fms/mpp/include/mpp_write.h $(SRCROOT)fms/fms/fms_io_unstructured_read.inc $(SRCROOT)fms/mosaic/create_xgrid.h $(SRCROOT)fms/mpp/include/mpp_unstruct_domain.inc $(SRCROOT)fms/mpp/include/mpp_transmit.inc $(SRCROOT)fms/mpp/include/mpp_io_unstructured_write.inc $(SRCROOT)fms/mpp/include/mpp_read_compressed.h $(SRCROOT)fms/mpp/include/mpp_update_domains2D.h $(SRCROOT)fms/mpp/include/mpp_io_misc.inc $(SRCROOT)fms/mpp/include/mpp_alltoall_sma.h $(SRCROOT)fms/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)fms/mpp/include/mpp_do_redistribute.h $(SRCROOT)fms/mpp/include/mpp_chksum_scalar.h $(SRCROOT)fms/fms/fms_io_unstructured_get_field_size.inc $(SRCROOT)fms/mosaic/mosaic_util.h
OBJ = random_numbers.o MersenneTwister.o xgrid.o write_ocean_data.o station_data.o axis_utils.o grid.o nsclock.o horiz_interp_type.o mpp_efp.o constants.o get_cal_time.o fft.o fm_util.o gaussian_topog.o topography.o horiz_interp.o time_interp.o test_mpp_io.o fms.o diag_grid.o drifters_io.o fms_io.o mpp_io.o oda_core_ecda.o field_manager.o column_diagnostics.o gradient.o interpolator.o test_mpp_domains.o time_manager.o diag_data.o sat_vapor_pres_k.o quicksort.o diag_table.o mosaic.o astronomy.o diag_manager.o atmos_ocean_fluxes.o create_xgrid.o affinity.o time_interp_external.o data_override.o amip_interp.o test_xgrid.o drifters.o test_mpp_pset.o mpp_memutils.o mpp.o test_fms_io.o diag_manifest.o cloud_interpolator.o threadloc.o coupler_types.o mpp_utilities.o mpp_parameter.o oda_core.o read_mosaic.o horiz_interp_spherical.o test_mpp.o xbt_drop_rate_adjust.o ensemble_manager.o test_unstructured_fms_io.o drifters_input.o diag_output.o sat_vapor_pres.o horiz_interp_bilinear.o test_horiz_interp.o gradient_c2l.o tracer_manager.o memutils.o block_control.o mpp_domains.o mpp_pset.o mpp_data.o drifters_comm.o platform.o stock_constants.o fft99.o diag_axis.o oda_types.o tridiagonal.o drifters_core.o interp.o horiz_interp_conserve.o diag_util.o horiz_interp_bicubic.o memuse.o mosaic_util.o

clean: neat
	$(RM) $(OBJ) libfms.a
	$(RM) *.mod
neat:
	$(RM) $(TMPFILES)

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libfms.a: $(OBJ)
	$(AR) $(ARFLAGS) libfms.a $(OBJ)
