# Makefile for Experiment 'SM2_exec'


BUILDROOT = $(dir $(abspath $(firstword $(MAKEFILE_LIST))))
SRCROOT = $(abspath $(BUILDROOT)../src)/ # Final '/' added to denote directory

MK_TEMPLATE = $(BUILDROOT)templates/intel.mk

# Override default compile options. See templates/intel.mk for more compile options.
# BLD_TYPE configures the type of build, PROD (production), REPRO (regression) or
# DEBUG (debug options)
BLD_TYPE = PROD
OPENMP = t # Build with OpenMP enabled

include $(MK_TEMPLATE)

SUBMAKEFLAGS = BUILDROOT=$(BUILDROOT) SRCROOT=$(SRCROOT) MK_TEMPLATE=$(MK_TEMPLATE) OPENMP=$(OPENMP)

fms_SM2_exec.x: coupler/libcoupler.a atmos_dyn/libatmos_dyn.a ice/libice.a land/libland.a atmos_phys/libatmos_phys.a ocean/libocean.a fms/libfms.a
	$(LD) $(addprefix $(BUILDROOT),$^) $(LDFLAGS) -o $(BUILDROOT)$@ $(STATIC_LIBS)

fms/libfms.a:  FORCE
	@echo $(MAKE) $(SUBMAKEFLAGS) --directory=$(BUILDROOT)fms $(@F)
	$(MAKE) $(SUBMAKEFLAGS) --directory=$(BUILDROOT)fms $(@F)

ocean/libocean.a: fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) --directory=$(BUILDROOT)ocean $(@F)

atmos_phys/libatmos_phys.a: fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) --directory=$(BUILDROOT)atmos_phys $(@F)

atmos_dyn/libatmos_dyn.a: atmos_phys/libatmos_phys.a fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) --directory=$(BUILDROOT)atmos_dyn $(@F)

ice/libice.a: ocean/libocean.a fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) --directory=$(BUILDROOT)ice $(@F)

land/libland.a: fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) --directory=$(BUILDROOT)land $(@F)

coupler/libcoupler.a: ice/libice.a atmos_dyn/libatmos_dyn.a atmos_phys/libatmos_phys.a land/libland.a ocean/libocean.a fms/libfms.a FORCE
	$(MAKE) $(SUBMAKEFLAGS) --directory=$(BUILDROOT)coupler $(@F)

FORCE:

clean:
	$(MAKE) --directory=fms clean
	$(MAKE) --directory=ocean clean
	$(MAKE) --directory=atmos_phys clean
	$(MAKE) --directory=atmos_dyn clean
	$(MAKE) --directory=ice clean
	$(MAKE) --directory=land clean
	$(MAKE) --directory=coupler clean

clean_all: clean
	$(RM) *.x
