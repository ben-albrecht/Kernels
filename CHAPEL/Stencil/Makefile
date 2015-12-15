# Chapel Makefile

# TODO Find a way to make this fit into the existing Makefile system

OPTFLAGS = --fast

CHPL = chpl

SOURCES = stencil-serial.chpl stencil-datablock.chpl
EXECUTABLE = $(SOURCES:.chpl=)

all: $(EXECUTABLE)

%: %.chpl
	$(CHPL) $(OPTFLAGS) $*.chpl -o $*

clean:
	rm $(EXECUTABLE)
