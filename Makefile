CXX=g++
CFLAGS=-O6 -g -Wall -std=c++11
#-I/usr/recon/lib -I/usr/image/lib -L/usr/recon/lib -L/usr/image/lib 

LFLAGS=-lm
#-lim

SRCS=$(addsuffix .cpp,dncat_main dncatsubs global_vars nurbs)
OBJS=$(SRCS:%.cpp=%.o)

run: dncat_bin test.par
	./$^ test
	python -BO showslice.py test_1.bin

all: dncat_bin
# dncat_im

dncat_%: $(OBJS) dncat_output_%.o
	$(CXX) $(CFLAGS) -o $@ $^ $(LFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(OBJS)
	rm -f dncat_output_im.o dncat_output_bin.o
	rm -f dncat_im dncat_bin dncat_bin.exe dncat_im.exe
