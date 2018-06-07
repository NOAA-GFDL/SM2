# Makefile created by mkmf 19.3.0

CPPDEFS = -DINTERNAL_FILE_NML

OTHERFLAGS = -I$(BUILDROOT)fms

include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libocean.a
ocean_grids.o: $(SRCROOT)ocean_slab/ocean_grids.F90
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ocean_slab/ocean_grids.F90
ocean_model.o: $(SRCROOT)ocean_slab/ocean_model.F90 ocean_grids.o
	$(FC) $(CPPDEFS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)ocean_slab/ocean_model.F90

SRC = $(SRCROOT)ocean_slab/ocean_model.F90 $(SRCROOT)ocean_slab/ocean_grids.F90
OBJ = ocean_model.o ocean_grids.o

clean: neat
	$(RM) $(OBJ) libocean.a
	$(RM) *.mod

neat:
	$(RM) $(TMPFILES)

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libocean.a: $(OBJ)
	$(AR) $(ARFLAGS) libocean.a $(OBJ)
