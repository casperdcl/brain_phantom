CFLAGS= -O6 -I/usr/recon/lib -I/usr/image/lib -L/usr/recon/lib -L/usr/image/lib 
LIBS= -lm #-lim 


#dncat_im 	: dncat_main.o dncat_output_im.o 
#		cc $(CFLAGS) -o dncat_im dncat_main.o  \
#			dncat_output_im.o $(LIBS)

dncat_bin	: dncat_main.o
		cc $(CFLAGS) -o dncat_bin dncat_main.o -lm 
clean:
	rm -f *.o dncat_im 
