OBJECTS  = match.o graph.o tools.o cal.o 
CXXFLAGS = -O3 -g -std=c++0x
LDFLAGS  = -lemon

matching : $(OBJECTS)
	g++ -o matching $(OBJECTS) $(CXXFLAGS) $(LDFLAGS)

match.o : match.cc
	g++ -c match.cc $(CXXFLAGS)

graph.o : graph.cc
	g++ -c graph.cc $(CXXFLAGS)

cal.o : cal.cc
	g++ -c cal.cc $(CXXFLAGS)

tools.o : tools.cc
	g++ -c tools.cc $(CXXFLAGS)

.PHONY : clean

clean : 
	rm -f matching *.o
