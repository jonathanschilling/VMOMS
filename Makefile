
.PHONY: all clean

all: vmoms
clean:
	rm -f vmoms

vmoms: vmoms.f
	gfortran --std=legacy vmoms.f -o vmoms
