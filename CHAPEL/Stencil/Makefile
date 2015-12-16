# Chapel Makefile

# TODO Find a way to make this fit into the existing Makefile system

OPTFLAGS = --fast -sassertNoSlicing

CHPL = chpl

SOURCES = stencil-serial.chpl stencil-data.chpl stencil-task.chpl
EXECUTABLE = $(SOURCES:.chpl=)

all: $(EXECUTABLE)

%: %.chpl
	$(CHPL) $(OPTFLAGS) $*.chpl -o $*

clean:
	rm -f $(EXECUTABLE)
