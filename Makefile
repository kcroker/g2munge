all : g2munge.o tkgadget2.o burns.o
	gcc -o g2munge tkgadget2.o g2munge.o burns.o -g -lgsl -lm -lgslcblas
g2munge.o : g2munge.c	
	gcc -DDEBUG -c g2munge.c -o g2munge.o -g
burns.o : burns.c
	gcc -c burns.c -o burns.o -g
tkgadget2.o : tkgadget2.c
	gcc -DDEBUG -c tkgadget2.c -o tkgadget2.o -g 
clean :
	rm -f *.o
