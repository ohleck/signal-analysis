all: oct/gauss_mod2_cc.oct oct/pll_cc.oct oct/early_late_cc.oct

oct/%.oct: src/%.cc
	cd oct; mkoctfile -I../include ../$< && rm $*.o

clean:
	rm -f oct/*
