BIN = bin
MODULES = ./Modules
F90 = ifort
F77 = ifort
#F90FLAGS = -O0 -fpe0 -CB -check all -traceback -warn all -module $(MODULES) -I $(MODULES)
F90FLAGS = -O3  -module $(MODULES) -I $(MODULES)
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

#all:	cut_vct.x put_vct.x put_vat.x put.x
all:	cut_vct.x put_vct.x  put.x cube.x cube_shell.x read_write_n00.x

cut_vct.x:	$(BIN)/cut_vct.x
put_vct.x:	$(BIN)/put_vct.x
put_vat.x:	$(BIN)/put_vat.x
put.x:	        $(BIN)/put.x
cube.x:		$(BIN)/cube.x
cube_shell.x:	$(BIN)/cube_shell.x
read_write_n00.x:	$(BIN)/read_write_n00.x

clean:
	@rm -f *.o *.il *.mod $(MODULES)/*.mod libbabel.a $(BIN)/*.x

libbabel.a: libbabel.a(module_octahedre.o) \
            libbabel.a(module_migration_sia.o) \
            libbabel.a(module_111_smooth.o) \
            libbabel.a(module_100_smooth.o) \
            libbabel.a(math.o) \
            libbabel.a(read_tsia.o) \
            libbabel.a(module_pbc.o) \
            libbabel.a(some_math.o) \
            libbabel.a(sub_octa.o) \
            libbabel.a(common_diverse.o) \
	    libbabel.a(read_input_write_output.o) \
	    libbabel.a(kazuto_style.o) \
            libbabel.a(nr.o) \
            libbabel.a(lapack.o) 
 


$(BIN)/cut_vct.x: libbabel.a cut_vct.o
	$(LD) $(LDFLAGS) -o $(BIN)/cut_vct.x cut_vct.o libbabel.a $(LIBLOC)

$(BIN)/put_vct.x: libbabel.a put_vct.o
	$(LD) $(LDFLAGS) -o $(BIN)/put_vct.x put_vct.o libbabel.a $(LIBLOC)

$(BIN)/put_vat.x: libbabel.a put_vat.o
	$(LD) $(LDFLAGS) -o $(BIN)/put_vat.x put_vat.o libbabel.a $(LIBLOC)
$(BIN)/put.x:  
	$(F90)  module_pbc.f90 math.f90 lapack.f putSIAinBCC_V4.f90 -o $(BIN)/put.x  $(LIBLOC)
$(BIN)/cube.x:  
	$(F90)  cubeFCCandBCC_v1.f -o $(BIN)/cube.x  $(LIBLOC)
$(BIN)/cube_shell.x:  
	$(F90)  cubeFCCandBCC_v3.f90 -o $(BIN)/cube_shell.x  $(LIBLOC)
$(BIN)/read_write_n00.x:  
	$(F90)  read_write_n00.f90 -o $(BIN)/read_write_n00.x  $(LIBLOC)

