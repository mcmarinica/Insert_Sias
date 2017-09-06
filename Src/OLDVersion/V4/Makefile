BIN = bin
MODULES = ./Modules
F90 = ifort
F77 = ifort
F90FLAGS = -O0 -fpe0 -CB -check all -traceback -warn all -module $(MODULES) -I $(MODULES)
LIBLOC =  -Vaxlib


##################################

LD = $(F90)
LDFLAGS = $(F90FLAGS)

.PRECIOUS: libbabel.a

.SUFFIXES: .f .f90 .o .a

.f.o:
	$(F77) $(F77FLAGS) -c $<

.f.a:
	$(F77) $(F77FLAGS) -c $<
	ar rv $@ $%
	rm -f $%

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.f90.a:
	$(F90) $(F90FLAGS) -c $<
	ar rv $@ $%
	rm -f $%

all:	in_gaos  put_v4.x

in_gaos:	$(BIN)/in_gaos

put_v4.x:	$(BIN)/put_v4.x


clean:
	@rm -f *.o *.il $(MODULES)/*.mod libbabel.a

libbabel.a: libbabel.a(module_octahedre.o) \
            libbabel.a(some_math.o) \
            libbabel.a(sub_octa.o) \
            libbabel.a(nr.o) 

$(BIN)/in_gaos: libbabel.a gaos_v4.o
	$(LD) $(LDFLAGS) -o $(BIN)/in_gaos gaos_v4.o libbabel.a $(LIBLOC)

$(BIN)/put_v4.x: libbabel.a put_v4.o
	$(LD) $(LDFLAGS) -o $(BIN)/put_v4.x put_v4.o libbabel.a $(LIBLOC)

