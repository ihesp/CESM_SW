CC=sw5cc -slave
CFLAGS=-msimd -O2 -ffast-math -OPT:IEEE_arith=2
FC=sw5f90 -slave
FFLAGS=-msimd -O2 -mp -ffast-math -OPT:IEEE_arith=2 -J$(SHAREDPATH)/$(lib) -I$(SHAREDPATH)/$(lib)
CPPFLAGS=-cpp
DEP=cpp
DEPFLAGS=-MM
lib=cpe

C_SRC=$(wildcard *.c)
F90_SRC=$(wildcard *.f90)
F90PP_SRC=$(wildcard *.F90)
#MKDIR:=$(shell if [ ! -d build ]; then mkdir build; fi; if [ ! -d dep ]; then mkdir dep; fi)
OBJS=$(patsubst %.c, $(SHAREDPATH)/$(lib)/%.o, $(C_SRC)) $(patsubst %.F90, $(SHAREDPATH)/$(lib)/%.o, $(F90PP_SRC)) $(patsubst %.f90, $(SHAREDPATH)/$(lib)/%.o, $(F90_SRC))
DEPS=$(patsubst %.c, $(SHAREDPATH)/$(lib)/%.d, $(C_SRC)) $(patsubst %.F90, $(SHAREDPATH)/$(lib)/%.d, $(F90PP_SRC))
all: $(SHAREDPATH)/$(lib)/libcpe.a
$(SHAREDPATH)/$(lib)/libcpe.a: $(OBJS)
	if [ -e $@ ]; then rm $@ ; fi
	ar cr $@ $^
$(SHAREDPATH)/$(lib)/%.o: %.c $(SHAREDPATH)/$(lib)/%.d
	$(CC) $(CFLAGS) -c $< -o $@
$(SHAREDPATH)/$(lib)/%.o: %.F90 $(SHAREDPATH)/$(lib)/%.d
	echo $(DEPS)
	$(FC) $(FFLAGS) $(CPPFLAGS) -c $< -o $@
$(SHAREDPATH)/$(lib)/%.o: %.f90 $(SHAREDPATH)/$(lib)/%.d
	$(FC) $(FFLAGS) -c $< -o $@

$(SHAREDPATH)/$(lib)/%.d: %.c
	$(DEP) $(DEPFLAGS) -MT $(patsubst $(SHAREDPATH)/$(lib)/%.d, $(SHAREDPATH)/$(lib)/%.o, $@) $< -o $@
$(SHAREDPATH)/$(lib)/%.d: %.F90
	$(DEP) $(DEPFLAGS) -MT $(patsubst $(SHAREDPATH)/$(lib)/%.d, $(SHAREDPATH)/$(lib)/%.o, $@) $< -o $@
install: $(SHAREDPATH)/$(lib)/libcpe.a
	cp $(SHAREDPATH)/$(lib)/libcpe.a $(SHAREDPATH)/lib
-include $(DEPS)
.PHONY: clean
clean:
	rm $(SHAREDPATH)/$(lib)/*.o $(SHAREDPATH)/$(lib)/*.mod dep/*
