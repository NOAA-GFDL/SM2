# Makefile created by mkmf 19.3.0

CPPDEFS = -DSPMD -Duse_shared_pointers -Duse_AM3_physics -DINTERNAL_FILE_NML

OTHERFLAGS = -I$(BUILDROOT)atmos_phys -I$(BUILDROOT)fms

include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libatmos_dyn.a
age_of_air.o: $(SRCROOT)atmos_fv_dynamics/tools/age_of_air.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/tools/age_of_air.F90
atmos_model.o: $(SRCROOT)atmos_drivers/coupled/atmos_model.F90 atmosphere.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_drivers/coupled/atmos_model.F90
atmosphere.o: $(SRCROOT)atmos_fv_dynamics/driver/coupled/atmosphere.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)atmos_fv_dynamics/model/fv_arrays.h $(SRCROOT)atmos_fv_dynamics/model/fv_point.inc fv_pack.o fv_diagnostics.o timingModule.o fv_restart.o fv_dynamics.o fv_arrays.o update_fv_phys.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include -I$(SRCROOT)atmos_fv_dynamics/model	$(SRCROOT)atmos_fv_dynamics/driver/coupled/atmosphere.F90
dyn_core.o: $(SRCROOT)atmos_fv_dynamics/model/dyn_core.F90 fv_pack.o shr_kind_mod.o pft_module.o sw_core.o timingModule.o mod_comm.o fv_arrays.o pmaxmin.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/model/dyn_core.F90
ecmfft.o: $(SRCROOT)atmos_fv_dynamics/model/ecmfft.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/model/ecmfft.F90
fill_module.o: $(SRCROOT)atmos_fv_dynamics/model/fill_module.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/model/fill_module.F90
fv_arrays.o: $(SRCROOT)atmos_fv_dynamics/model/fv_arrays.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)atmos_fv_dynamics/model/fv_point.inc mod_comm.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)fms/include -I$(SRCROOT)atmos_fv_dynamics/model	$(SRCROOT)atmos_fv_dynamics/model/fv_arrays.F90
fv_diagnostics.o: $(SRCROOT)atmos_fv_dynamics/tools/fv_diagnostics.F90 $(SRCROOT)atmos_fv_dynamics/model/fv_arrays.h $(SRCROOT)atmos_fv_dynamics/model/fv_point.inc fv_pack.o fv_arrays.o pmaxmin.o pv_module.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_fv_dynamics/model	$(SRCROOT)atmos_fv_dynamics/tools/fv_diagnostics.F90
fv_dynamics.o: $(SRCROOT)atmos_fv_dynamics/model/fv_dynamics.F90 fv_pack.o mapz_module.o dyn_core.o pft_module.o timingModule.o mod_comm.o fv_arrays.o shr_kind_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/model/fv_dynamics.F90
fv_pack.o: $(SRCROOT)atmos_fv_dynamics/model/fv_pack.F90 $(SRCROOT)atmos_fv_dynamics/model/fv_arrays.h $(SRCROOT)atmos_fv_dynamics/model/fv_point.inc timingModule.o shr_kind_mod.o mod_comm.o fv_arrays.o pft_module.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_fv_dynamics/model	$(SRCROOT)atmos_fv_dynamics/model/fv_pack.F90
fv_restart.o: $(SRCROOT)atmos_fv_dynamics/tools/fv_restart.F90 $(SRCROOT)atmos_fv_dynamics/model/fv_arrays.h $(SRCROOT)atmos_fv_dynamics/model/fv_point.inc fv_pack.o mod_comm.o fv_arrays.o pmaxmin.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_fv_dynamics/model	$(SRCROOT)atmos_fv_dynamics/tools/fv_restart.F90
getmax.o: $(SRCROOT)atmos_fv_dynamics/tools/getmax.F90 mod_comm.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/tools/getmax.F90
gmean.o: $(SRCROOT)atmos_fv_dynamics/tools/gmean.F90 fv_pack.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/tools/gmean.F90
init_dry_atm.o: $(SRCROOT)atmos_fv_dynamics/tools/init_dry_atm.F90 $(SRCROOT)atmos_fv_dynamics/model/fv_arrays.h $(SRCROOT)atmos_fv_dynamics/model/fv_point.inc fv_pack.o mod_comm.o fv_arrays.o pmaxmin.o shr_kind_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_fv_dynamics/model	$(SRCROOT)atmos_fv_dynamics/tools/init_dry_atm.F90
init_sw_ic.o: $(SRCROOT)atmos_fv_dynamics/tools/init_sw_ic.F90 fv_pack.o pv_module.o mod_comm.o shr_kind_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/tools/init_sw_ic.F90
mapz_module.o: $(SRCROOT)atmos_fv_dynamics/model/mapz_module.F90 fv_pack.o fill_module.o fv_diagnostics.o mod_comm.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/model/mapz_module.F90
mod_comm.o: $(SRCROOT)atmos_fv_dynamics/tools/mod_comm.F90 timingModule.o shr_kind_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/tools/mod_comm.F90
par_vecsum.o: $(SRCROOT)atmos_fv_dynamics/tools/par_vecsum.F90 mod_comm.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/tools/par_vecsum.F90
pft_module.o: $(SRCROOT)atmos_fv_dynamics/model/pft_module.F90 ecmfft.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/model/pft_module.F90
pmaxmin.o: $(SRCROOT)atmos_fv_dynamics/tools/pmaxmin.F90 mod_comm.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/tools/pmaxmin.F90
pv_module.o: $(SRCROOT)atmos_fv_dynamics/tools/pv_module.F90 fv_pack.o mod_comm.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/tools/pv_module.F90
set_eta.o: $(SRCROOT)atmos_fv_dynamics/tools/set_eta.F90 fv_pack.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/tools/set_eta.F90
shr_kind_mod.o: $(SRCROOT)atmos_fv_dynamics/model/shr_kind_mod.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/model/shr_kind_mod.F90
sw_core.o: $(SRCROOT)atmos_fv_dynamics/model/sw_core.F90 tp_core.o shr_kind_mod.o pft_module.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/model/sw_core.F90
timingModule.o: $(SRCROOT)atmos_fv_dynamics/tools/timingModule.F90 shr_kind_mod.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/tools/timingModule.F90
tp_core.o: $(SRCROOT)atmos_fv_dynamics/model/tp_core.F90 shr_kind_mod.o fill_module.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/model/tp_core.F90
tracer_2d.o: $(SRCROOT)atmos_fv_dynamics/model/tracer_2d.F90 fv_pack.o tp_core.o fv_diagnostics.o shr_kind_mod.o fill_module.o mod_comm.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/model/tracer_2d.F90
update_fv_phys.o: $(SRCROOT)atmos_fv_dynamics/model/update_fv_phys.F90 $(SRCROOT)atmos_fv_dynamics/model/fv_arrays.h $(SRCROOT)atmos_fv_dynamics/model/fv_point.inc fv_pack.o pft_module.o mod_comm.o fv_arrays.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)atmos_fv_dynamics/model	$(SRCROOT)atmos_fv_dynamics/model/update_fv_phys.F90
upper.o: $(SRCROOT)atmos_fv_dynamics/tools/upper.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)atmos_fv_dynamics/tools/upper.F90

SRC = $(SRCROOT)atmos_fv_dynamics/tools/getmax.F90 $(SRCROOT)atmos_fv_dynamics/tools/pmaxmin.F90 $(SRCROOT)atmos_fv_dynamics/model/sw_core.F90 $(SRCROOT)atmos_fv_dynamics/tools/age_of_air.F90 $(SRCROOT)atmos_fv_dynamics/model/fv_dynamics.F90 $(SRCROOT)atmos_fv_dynamics/model/fv_arrays.F90 $(SRCROOT)atmos_fv_dynamics/tools/par_vecsum.F90 $(SRCROOT)atmos_fv_dynamics/model/shr_kind_mod.F90 $(SRCROOT)atmos_fv_dynamics/driver/coupled/atmosphere.F90 $(SRCROOT)atmos_fv_dynamics/model/tracer_2d.F90 $(SRCROOT)atmos_fv_dynamics/model/update_fv_phys.F90 $(SRCROOT)atmos_fv_dynamics/tools/fv_restart.F90 $(SRCROOT)atmos_fv_dynamics/model/mapz_module.F90 $(SRCROOT)atmos_fv_dynamics/tools/init_dry_atm.F90 $(SRCROOT)atmos_fv_dynamics/tools/gmean.F90 $(SRCROOT)atmos_fv_dynamics/model/tp_core.F90 $(SRCROOT)atmos_fv_dynamics/tools/pv_module.F90 $(SRCROOT)atmos_fv_dynamics/tools/init_sw_ic.F90 $(SRCROOT)atmos_fv_dynamics/model/dyn_core.F90 $(SRCROOT)atmos_fv_dynamics/model/fill_module.F90 $(SRCROOT)atmos_fv_dynamics/tools/upper.F90 $(SRCROOT)atmos_drivers/coupled/atmos_model.F90 $(SRCROOT)atmos_fv_dynamics/tools/mod_comm.F90 $(SRCROOT)atmos_fv_dynamics/tools/timingModule.F90 $(SRCROOT)atmos_fv_dynamics/tools/fv_diagnostics.F90 $(SRCROOT)atmos_fv_dynamics/model/ecmfft.F90 $(SRCROOT)atmos_fv_dynamics/tools/set_eta.F90 $(SRCROOT)atmos_fv_dynamics/model/pft_module.F90 $(SRCROOT)atmos_fv_dynamics/model/fv_pack.F90 $(SRCROOT)fms/include/fms_platform.h $(SRCROOT)atmos_fv_dynamics/model/fv_arrays.h $(SRCROOT)atmos_fv_dynamics/model/fv_point.inc
OBJ = getmax.o pmaxmin.o sw_core.o age_of_air.o fv_dynamics.o fv_arrays.o par_vecsum.o shr_kind_mod.o atmosphere.o tracer_2d.o update_fv_phys.o fv_restart.o mapz_module.o init_dry_atm.o gmean.o tp_core.o pv_module.o init_sw_ic.o dyn_core.o fill_module.o upper.o atmos_model.o mod_comm.o timingModule.o fv_diagnostics.o ecmfft.o set_eta.o pft_module.o fv_pack.o

clean: neat
	$(RM) $(OBJ) libatmos_dyn.a
	$(RM) *.mod

neat:
	$(RM) $(TMPFILES)

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libatmos_dyn.a: $(OBJ)
	$(AR) $(ARFLAGS) libatmos_dyn.a $(OBJ)
