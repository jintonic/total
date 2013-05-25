PREFIX=$(HOME)
LIBNAME=TOTAL

CC     =gfortran
CFLAGS =-fPIC -w
SOFLAGS=-shared

SRCS=test.f
EXES=$(SRCS:.f=.exe)

SOURCES=$(filter-out $(SRCS), $(wildcard *.f))
OBJECTS=$(SOURCES:.f=.o)
LIBRARY=lib$(LIBNAME).so

all:$(EXES)

$(EXES):$(LIBRARY)
	$(CC) $(CFLAGS) -L./ -l$(LIBNAME) -o $@ $(SRCS)

$(LIBRARY):$(OBJECTS)
	$(CC) $(CFLAGS) $(SOFLAGS) $^ -o $@

$(OBJECTS):%.o:%.f
	$(CC) $(CFLAGS) -c $^
clean:
	rm -f $(EXES) *.o $(LIBRARY)

install:
	@echo "PREFIX=$(PREFIX)"
	@echo -n "checking if $(PREFIX) exists..."
	@if [ -d $(PREFIX) ]; then \
	  echo "yes."; \
	else \
	  echo "no."; \
	  echo "mkdir $(PREFIX)..."; \
	  mkdir $(PREFIX); \
	fi
	@echo -n "copying $(LIBRARY) to $(PREFIX)/lib..."
	@if [ -d $(PREFIX)/lib ]; then \
	  cp $(LIBRARY) $(PREFIX)/lib; \
	else \
	  mkdir $(PREFIX)/lib; \
	  cp $(LIBRARY) $(PREFIX)/lib; \
	fi; 
	@echo "done."; 

uninstall:
	rm -f $(PREFIX)/lib/$(LIBRARY)

.PHONY: all clean install uninstall
