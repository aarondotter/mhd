FC = gfortran
FLAGS= -g -fbacktrace -Wall -O3


LIB=libMHD.a
LIB_OBJS = cray_emul.o mod_mhd.o fedi.o dfermi200.o \
	 lineqsolv_sun.o free_energy_rel_relf4.o saha_neweta_rel.o lir.o \
	 vinta.o utils.o matgen.o lambda.o window.o output.o fscan.o \
	 setphi.o fjk.o finddn.o

$(LIB): $(LIB_OBJS)
	ar crs $(LIB) $(LIB_OBJS)

test: test_MHD.o $(LIB)
	$(FC) -o ../mhd.exe $^

all: test

test_MHD.o: $(LIB)

%.o: %.f
	$(FC) $(FLAGS) -c $<

%.o: %.f90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f test_MHD.o ${LIB_OBJS} $(LIB) core *~
