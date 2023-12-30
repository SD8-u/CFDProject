SRCDIR = src
OUTDIR = build

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OUTDIR)/%.o,$(SOURCES))
TARGET = $(OUTDIR)/main


CXX = /usr/bin/mpic++
CXXFLAGS = -fdiagnostics-color=always -g -I/usr/include/petsc -I/usr/include -I/usr/include/python3.10 -I$(CURDIR)/venv/lib/python3.10/site-packages/pybind11/include -I$(CURDIR)/src
LDFLAGS = -L/usr/lib -lpetsc -lstdc++ -lgmsh -lpython3.10

all: directories $(TARGET)

directories:
	@mkdir -p $(OUTDIR)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OUTDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OUTDIR)
