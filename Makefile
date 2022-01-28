
VMOMS=aux bcaux curfun d02agf d02agy d02agz f03aff \
      f04ajf geotrn gpasma i1mach mhdeq p01aaf     \
      pltmeq prsfun prsol prtmeq raaux sdot        \
      seval sincos speval splaan spmpar vmoms      \
      x03aaf x03aaz

OBJECTS=$(VMOMS:%=obj/%.o)

F77:=gfortran --std=legacy

.PHONY: all clean

all: vmoms
clean:
	rm -rf obj
	rm -f vmoms

obj/%.o: src/%.f
	@mkdir -p obj
	$(F77) -c $^ -o $@

vmoms: $(OBJECTS)
	$(F77) $(OBJECTS) -o vmoms
