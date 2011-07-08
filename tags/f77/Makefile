PROG = kb3d	

SRCS =	kb3d.f90 

OBJS =	kb3d.f90.o 

MATHKEI =  /SX/opt/mathkeisan/inst

MATHKEI_LIB = $(MATHKEI)/lib0
MATHKEI_INC = $(MATHKEI)/include

#on cross-compile machines
#32bit
#F90 = sxf90 -dw
#64bit
F90 = sxf90 -ew

CFLAGS = -I$(MATHKEI_INC) -g
FFTFLAG = -lfft 
LDFLAGS = $(FFTFLAG) -L$(MATHKEI_LIB)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(OBJS) $(LDFLAGS) -o $@

clean:
	rm -f $(PROG) $(OBJS) *.mod

%.SUFFIXES: $(SUFFIXES) .f90

%.f90.o: %.f90
	$(F90) $(CFLAGS) -c $< -o $@

%.mod.o:
	$(F90) $(CFLAGS) -c $< -o $@
	
