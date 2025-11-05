CXX       := clang++
CXX_FLAGS := -std=c++17 -v -I/usr/src/myapp/include/ -L/usr/local/lib/ -O3

BIN     := ./bin
SRC     := ./src/pulsatility
INCLUDE := include

LIBRARIES   := -liomp5 -Xclang -fopenmp  
EXECUTABLE  := main


all: $(BIN)/$(EXECUTABLE)

run: clean all
	clear
	$(BIN)/$(EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp
	$(CXX) $(CXX_FLAGS) $^ -o $@ $(LIBRARIES)

clean:
	-rm $(BIN)/*