INCLUDE=-I"src"
OPTIMIZE=-O3
CXXFLAGS=$(OPTIMIZE) $(INCLUDE) -c -std=c++11 -Wall
LDFLAGS=
LDLIBS=

.PHONY: all clean

SOURCEPATH=src
SOURCES=$(SOURCEPATH)/aabb.cpp \
		$(SOURCEPATH)/fluidsim.cpp \
		$(SOURCEPATH)/interpolation.cpp \
		$(SOURCEPATH)/levelsetutils.cpp \
		$(SOURCEPATH)/macvelocityfield.cpp \
		$(SOURCEPATH)/main.cpp \
		$(SOURCEPATH)/meshlevelset.cpp \
		$(SOURCEPATH)/particlelevelset.cpp \
		$(SOURCEPATH)/pressuresolver.cpp \
		$(SOURCEPATH)/trianglemesh.cpp \
		$(SOURCEPATH)/viscositysolver.cpp \
		$(SOURCEPATH)/vmath.cpp

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=fluidsim

ifeq ($(OS),Windows_NT)
	RM = del /F /Q 
	SEP=\\
	ERRIGNORE = 2>NUL || true
else
	RM = rm -rf 
	SEP=/
	ERRIGNORE = 2>/dev/null
endif

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(LDFLAGS) $(OBJECTS) $(LDLIBS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	@$(RM) $(SOURCEPATH)$(SEP)*.o $(EXECUTABLE)* $(ERRIGNORE)