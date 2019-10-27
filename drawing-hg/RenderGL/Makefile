# Makefile for Json scene file viewer and rasterizer
#
# Conrad Steenberg <conrad.steenberg@caltech.edu>
# Nov 6, 2008

INCLUDES= -I/usr/X11R6/include -Ijsoncpp/include
LIB_DIR= -L/usr/X11R6/lib -Llib
LIBS= -ljsoncpp -lgd -lgle -lGL -lGLU -lglut -lm
OSMESALIB=-lOSMesa
JSONLIB=
CXXFLAGS=-O3 -g $(INCLUDES)

# Objetcs
BASEOBJECTS=RenderGL.o Objects.o
DOBJECTS=DemoScene.o
JOBJECTS=JsonViewer.o JsonScene.o
ROBJECTS=JsonRenderer.o JsonScene.o

# Headers
BASEHEADERS=RenderGL.hpp Objects.hpp 
JHEADERS=JsonScene.hpp

# Sources
BASESOURCES=RenderGL.cpp Objects.cpp 
DSOURCES=DemoScene.cpp
JSOURCES=JsonScene.cpp JsonViewer.cpp

# Executables
DEXE = DemoScene
JEXE = JsonViewer
REXE = JsonRenderer

#Targets
all : $(DEXE) $(JEXE) $(REXE)

$(DEXE) : $(JSONLIB) $(BASEOBJECTS) $(DOBJECTS)
	$(CXX) $(BASEOBJECTS) $(DOBJECTS) $(LIB_DIR) $(LIBS) -o $@
	strip $@

$(JEXE) : $(JSONLIB) $(BASEOBJECTS) $(JOBJECTS) 
	$(CXX) $(BASEOBJECTS) $(JOBJECTS) $(LIB_DIR) $(LIBS) -o $@
	strip $@

$(REXE) : $(JSONLIB) $(BASEOBJECTS) $(ROBJECTS) 
	$(CXX) -DUSE_OSMESA $(BASEOBJECTS) $(ROBJECTS) $(LIB_DIR) $(OSMESALIB) $(LIBS) -o $@
	strip $@

clean:
	rm -f $(BASEOBJECTS) \
	$(DEXE) $(DOBJECTS) \
	$(JEXE) $(JOBJECTS) \
	$(REXE) JsonRenderer.o

$(DOBJECTS): $(BASEOBJECTS)

JsonRenderer.o: JsonViewer.cpp
	$(CXX) -DUSE_OSMESA -c $(CXXFLAGS) $< -o $@

#$(JOBJECTS): $(JHEADERS)



$(JSONLIB):
	$(MAKE) -C jsoncpp
