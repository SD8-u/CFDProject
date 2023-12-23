SRCDIR = src
OUTDIR = build

SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OUTDIR)/%.o,$(SOURCES))
TARGET = $(OUTDIR)/main


CXX = /usr/bin/mpic++
CXXFLAGS = -fdiagnostics-color=always -g -I/usr/include/petsc -I/usr/include
LDFLAGS = -L/usr/lib -lpetsc -lstdc++ -lgmsh

all: directories $(TARGET)

directories:
	@mkdir -p $(OUTDIR)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OUTDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OUTDIR)
