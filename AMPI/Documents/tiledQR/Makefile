include Makefile.opts

BIN_FILES=tiledQR

INCLUDES = -I$(INSTALL_PATH)/include -I./include
DEFS = -L$(INSTALL_PATH)/lib/
LDADD = MC_scheduler.o TT-FT_scheduler.o TS-FT_scheduler.o TT-PT_scheduler.o TS-PT_scheduler.o ASAP_scheduler.o GREEDY_scheduler.o main_scheduler.o task.o worker.o mapping_policy.o main.o

all: $(BIN_FILES)

$(BIN_FILES): $(LDADD)
	$(CC) $(INCLUDES) $(DEFS) $(CFLAGS) $^ $(LIBS) -o $@ 

%.o: ./src/%.c
	$(CC) $(INCLUDES) $(DEFS) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.o *~
cleaner:clean
	rm -f $(BIN_FILES)

.PHONY: clean
.SUFFIXES:
