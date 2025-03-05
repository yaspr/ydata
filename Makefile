CC=gcc
CFLAGS=-Wall -Wextra -g
LFLAGS=-shared -fPIC
OFLAGS=-Ofast
MFLAGS=-lm

all: libydata.so test

libydata.so: ydata.c
	$(CC) $(CFLAGS) $(OFLAGS) $(LFLAGS) $< ../ynotif/ynotif.c -o $@ $(MFLAGS)

test: test.c
	$(CC) $(CFLAGS) $(OFLAGS) $< ydata.c ../ynotif/ynotif.c -o $@ $(MFLAGS)

runtest: test
	@./test

clean:
	rm -Rf *.so *.csv test
