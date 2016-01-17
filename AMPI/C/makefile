CFLAGS  = -O
CC       = cc
LIBS     = -lm -llapack -lblas ./quark-0.9.0/libquark.a -lpthread -lhwloc


run.out: qr.o
	$(CC) -o $@ qr.o ./tools.c $(LDFLAGS) $(LIBS)

debug.out: qr.o
	$(CC) -o $@ qr.o ./tools.c $(LDFLAGS) $(LIBS) -g
	
clean: 
	rm *.out
	rm *.o
