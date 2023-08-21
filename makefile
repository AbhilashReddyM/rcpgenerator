CXX = g++
FLAGS = -O3 -march=native

INCLUDES = include
bin = rcpgenerator

src = $(wildcard ./src/*.cpp)
obj = $(src:.cpp=.o)

$(bin): $(obj)
	$(CXX) $(FLAGS) $^ -o $@

%.o: %.cpp
	$(CXX) $(FLAGS) -I$(INCLUDES) -c $< -o $@

# clean for windows and linux
ifeq ($(OS),Windows_NT)
  clean_cmd = del /Q /S *.o *.exe *.a
else
  clean_cmd = rm -f ./src/*.o $(bin)
endif

.PHONY: clean
clean:
	@$(clean_cmd)
