SOURCES = $(wildcard /home/savan/projects/petscInit/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = main

CXX = /usr/bin/mpic++
CXXFLAGS = -fdiagnostics-color=always -g -I/usr/include/petsc -I/usr/include
LDFLAGS = -L/usr/lib -lpetsc -lstdc++ -lgmsh

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)
