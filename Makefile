CXX=g++
CFLAGS=-O6 -g -Wall
#-I/usr/recon/lib -I/usr/image/lib -L/usr/recon/lib -L/usr/image/lib 

LFLAGS=-lm
#-lim

SRCS=$(addsuffix .cpp,brain_phantom dncat_main dncatsubs global_vars nurbs)
OBJS=$(SRCS:%.cpp=%.o)

dncat_im: $(OBJS) dncat_output_im.o
	$(CXX) $(CFLAGS) -o $@ $^ $(LFLAGS)

dncat_bin: $(OBJS) dncat_output_bin.o
	$(CXX) $(CFLAGS) -o $@ $^ $(LFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS)
	rm -f dncat_output_im.o dncat_output_bin.o
	rm -f dncat_im dncat_bin
