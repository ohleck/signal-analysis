all: oct/s4285_cc.oct oct/gauss_mod2_cc.oct oct/pll_cc.oct oct/early_late_cc.oct oct/lfsr_gen_cc.oct oct/demod_msk_helper_cc.oct

oct/%.oct: src/%.cc include/*.hpp
	cd oct; mkoctfile -I../include ../$< && rm $*.o

clean:
	rm -f oct/*
