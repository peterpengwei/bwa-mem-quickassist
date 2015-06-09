TARGET  := bwa

SPATH := .
OPATH := .

SOURCES := $(notdir $(wildcard $(SPATH)/*.c))
OBJECTS := $(patsubst %.c, %.o, $(SOURCES))
LIBPATH := .

CXX = gcc
CXXFLAGS 	  = -g -Wall -Wno-unused-function -O2
WRAP_MALLOC = -DUSE_MALLOC_WRAPPERS
LDFLAGS     = -DHAVE_PTHREAD $(WRAP_MALLOC) -lm -lz -lpthread

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CFLAGS) *.o -o $(TARGET) $(LDFLAGS)

$(OBJECTS): %.o: $(SPATH)/%.c
	$(CXX) -c $(CXXFLAGS) -o $(OPATH)/$@ $<

clean:
	rm -f $(OPATH)/*.o $(TARGET)

