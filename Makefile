CXX = g++
CXXFLAGS = -g -Wall -fPIC -Wno-unused-variable
ROOTFLAGS = `root-config --cflags --glibs --libs` -lTreePlayer -lEG -lMinuit


# make a binary for every .cxx file
# all : $(patsubst %.cpp, %.o, $(wildcard *.cpp)) chi_square_cc0pi_christian univmake
 all : chi_square_cc0pi_christian univmake
# cc0pi_analyzer_org
# cc0pi_analyzer
# # rule for each one
# 
%: %.cpp
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o  $@ $< includes/*.o 
	
%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $*.o  -c $*.cpp 
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -o $* $*.o includes/*.o 


stv_root_dict.o:
	$(RM) stv_root_dict*.*
	rootcling -f stv_root_dict.cc -c LinkDef.h
	$(CXX) $(shell root-config --cflags --libs) -O3 \
	-fPIC -o stv_root_dict.o -c stv_root_dict.cc
	$(RM) stv_root_dict.cc
	
chi_square_cc0pi_christian: chi_square_cc0pi_christian.cpp
	$(CXX) $(CXXFLAGS) $(ROOTFLAGS) -O3 -o $@ $^ includes/*.o

univmake: univmake.C
	 $(CXX) $(shell root-config --cflags --libs) -O3 -o $@ $^
	
.PHONY: clean

.INTERMEDIATE: stv_root_dict.o

clean:
	rm -f $(wildcard *.o) $(patsubst %.cpp, %, $(wildcard *.cpp)) chi_square_cc0pi_christian univmake

