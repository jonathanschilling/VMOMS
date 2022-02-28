
VMOMS=aux bcaux curfun \
      geotrn gpasma mhdeq \
      pltmeq prsfun prsol prtmeq raaux \
      sincos vmoms

THIRD_PARTY=d02agf d02agy d02agz f03aff \
      f04ajf i1mach p01aaf sdot \
      seval  speval splaan spmpar x03aaf x03aaz

OBJECTS=$(VMOMS:%=obj/%.o)
THIRD_OBJECTS=$(THIRD_PARTY:%=obj_3rd/%.o)

F77:=gfortran --std=legacy

.PHONY: all clean

all: vmoms
clean:
	rm -rf obj
	rm -f vmoms

obj/%.o: src/%.f
	@mkdir -p obj
	$(F77) -c $^ -o $@

obj_3rd/%.o: src/3rd_party/%.f
	@mkdir -p obj_3rd
	$(F77) -c $^ -o $@

vmoms: $(OBJECTS) $(THIRD_OBJECTS)
	$(F77) $(OBJECTS) $(THIRD_OBJECTS) -o vmoms
