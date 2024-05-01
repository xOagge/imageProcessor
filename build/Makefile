# Directory structure
SRC_DIR := ../src
BUILD_DIR := .
BIN_DIR := ../bin
INC_DIR := ../inc
RES_DIR := ../res  # New res folder

# Find all source files
SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC))

# Compiler and flags
CXX := g++
ROOTCONFIG := root-config
ROOTCXXFLAGS := $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS := $(shell $(ROOTCONFIG) --libs)
CXXFLAGS := -std=c++20 -I$(INC_DIR) $(ROOTCXXFLAGS)
LDLIBS := $(ROOTLIBS)

# Executable name
TARGET := main.exe

# Rules to build object files from source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to link object files into executable
$(BIN_DIR)/$(TARGET): $(OBJ)
	@echo "Linking $@"
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDLIBS)
	@echo "Assimilating existence of res directory"  # New message for assimilation
	test -d $(RES_DIR) || mkdir $(RES_DIR)

# Phony target to clean up object files and executables
.PHONY: clean
clean:
	@echo "Cleaning..."
	rm -f $(OBJ) $(BIN_DIR)/$(TARGET)
