all: analyzer univmake chi_square_cc0pi

univmake: univmake.C stv_root_dict.o
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

analyzer: analyzer.C stv_root_dict.o
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

chi_square_cc0pi: chi_square_cc0pi.cpp stv_root_dict.o
	$(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^

stv_root_dict.o:
	$(RM) stv_root_dict*.*
	rootcling -f stv_root_dict.cc -c LinkDef.h
	$(CXX) $(shell root-config --cflags --libs) -O3 \
	-fPIC -o stv_root_dict.o -c stv_root_dict.cc
	$(RM) stv_root_dict.cc

.PHONY: clean

#.INTERMEDIATE: stv_root_dict.o

clean:
	$(RM) univmake analyzer chi_square_cc0pi stv_root_dict.o stv_root_dict_rdict.pcm
