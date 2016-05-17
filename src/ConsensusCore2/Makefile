CFLAGS=-fPIC -Wall -I../../inst/include/ -I./include/ --std=c++11 -Wno-deprecated-register -Wno-unused-local-typedefs
LDFLAGS=-force_load
#MODULES := align matrix models poa
#SRC_DIR   := $(addprefix src/,$(MODULES))


# Should be equivalent to your list of C files, if you don't build selectively
#SRC=$(wildcard *.cpp)
#CPP_FILES := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cpp))
#CPP_FILES := $(wildcard **/*.cpp)
#CPP_FILES := $(shell find . -name '.cpp')

# http://stackoverflow.com/questions/2483182/recursive-wildcards-in-gnu-make
rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *,%,$2),$d))
CPP_FILES := $(call rwildcard, , *.cpp)
OBJ_FILES := $(CPP_FILES:.cpp=.o)




# $@ - The name of the target
# $^ - All dependencies
# %.o - wildcard matching,
# See http://www.cprogramming.com/tutorial/makefiles_continued.html

libpbconsensus.a: $(OBJ_FILES)
	echo $(OBJ_FILES)
	echo $(CPP_FILES)
	$(AR) -r $@ $(OBJ_FILES)


%.o: %.cpp
	g++ $(CFLAGS) $(LIBS) -c -o $@ $<


clean:
	rm *.o
	rm *.a
