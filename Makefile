
COMMON = coords.o swapbytes.o common.o traj.o
OBJECTS = $(COMMON) rdc_funcs.o 
#TARGETS = rdc boxsize qcalc_gmx noe_gmx
TARGETS = noe_gmx rdc s2 \
	  noe_corfun 

# CFLAGS = -I$(HOME)/include -O3
CFLAGS = -I/usr/include -O3
#LDFLAGS = -L$(HOME)/lib #-static
LDFLAGS =  -L/usr/lib -lgsl -lgslcblas #$(HOME)/lib #-static
#LDFLAGS = /usr/lib/libgsl.so.0 /usr/lib/libgslcblas.so.0
GSLFLAGS = -lgsl -lgslcblas -lm 
CXX = g++ $(CFLAGS)

all: $(TARGETS)

coords.o: coords.cc coords.h

traj.o: traj.cc traj.h

rdc_funcs.o: rdc_funcs.cc rdc_funcs.h

swapbytes.o: swapbytes.cc swapbytes.h

common.o: common.cc common.h

rdc: rdc.cc $(COMMON) rdc_funcs.o 
	$(CXX) $(GSLFLAGS) -o rdc rdc.cc rdc_funcs.o  $(COMMON) $(LDFLAGS)

noe: noe.cc $(COMMON) rdc_funcs.o 
	$(CXX) $(GSLFLAGS) -o noe noe.cc rdc_funcs.o  $(COMMON) $(LDFLAGS)

noe_corfun: noe_corfun.cc $(COMMON) rdc_funcs.o 
	$(CXX) $(GSLFLAGS) -o noe_corfun noe_corfun.cc rdc_funcs.o  $(COMMON) $(LDFLAGS)

noe_gmx: noe_gmx.cc coords.o swapbytes.o common.o traj.o
	$(CXX) -o noe_gmx noe_gmx.cc coords.o swapbytes.o common.o traj.o $(LDFLAGS)

s2: s2.cc $(OBJECTS)
	$(CXX) -o s2 s2.cc $(OBJECTS) $(LDFLAGS)

install: $(TARGETS)
	cp $(TARGETS) $(HOME)/bin

clean: 
	-rm $(OBJECTS) 

distclean: clean
	-rm $(TARGETS)
