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

all:	cut_vct.x put_vct.x put_vat.x put.x

cut_vct.x:	$(BIN)/cut_vct.x
put_vct.x:	$(BIN)/put_vct.x
put_vat.x:	$(BIN)/put_vat.x
put.x:	        $(BIN)/put.x


clean:
	@rm -f *.o *.il $(MODULES)/*.mod libbabel.a

libbabel.a: libbabel.a(module_octahedre.o) \
            libbabel.a(some_math.o) \
            libbabel.a(sub_octa.o) \
            libbabel.a(common_diverse.o) \
	    libbabel.a(read_input_write_output.o) \
            libbabel.a(nr.o) 


$(BIN)/cut_vct.x: libbabel.a cut_vct.o
	$(LD) $(LDFLAGS) -o $(BIN)/cut_vct.x cut_vct.o libbabel.a $(LIBLOC)

$(BIN)/put_vct.x: libbabel.a put_vct.o
	$(LD) $(LDFLAGS) -o $(BIN)/put_vct.x put_vct.o libbabel.a $(LIBLOC)

$(BIN)/put_vat.x: libbabel.a put_vat.o
	$(LD) $(LDFLAGS) -o $(BIN)/put_vat.x put_vat.o libbabel.a $(LIBLOC)
$(BIN)/put.x: 
	$(F90) putSIAinBCC_V4.f90 -o $(BIN)/put.x  $(LIBLOC)

