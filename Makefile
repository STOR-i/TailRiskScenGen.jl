###################################
# Include and Library Directories #
###################################

INC_DIR:=-I include
LIB_DIR:=

#############################
# Libraries and other flags #
#############################

LIB:= -lcblas
CXX_FLAGS:=-std=c++11 -g -Wall -fPIC
C_FLAGS :=-std=c99 -fPIC -Wall

#########################
## Library source code ##
#########################

VPATH = include src bin

LIB_SOURCES := LCP_Solver.c LCP_julia.c

LDFLAGS+=-ldl

LIB_OBJ := $(patsubst %.c, bin/%.o, $(LIB_SOURCES))

bin/libScenGen.so: $(LIB_OBJ)
	g++ -shared $^ $(LIB) -o $@

bin/libScenGen.a: $(LIB_OBJ)
	ar rv $@ $?

bin/%.o: %.c %.h
	gcc -c $< $(C_FLAGS) $(INC_DIR) -o $@

bin/%.o: %.c
	gcc -c $< $(C_FLAGS) $(INC_DIR) -o $@

###########
## Tests ##
###########

# TEST_DIR := test
# TEST_SOURCES := $(shell find $(TEST_DIR) -type f -name "*.cpp")
# TEST_OBJ := $(patsubst $(TEST_DIR)/test_%.cpp, bin/test_%.o, $(TEST_SOURCES))
# TEST_EXE := $(patsubst $(TEST_DIR)/test_%.cpp, bin/test_%, $(TEST_SOURCES))

# $(TEST_OBJ): bin/test_%.o: $(TEST_DIR)/test_%.cpp bin/libScenGen.a
# 	g++ -c $(INC_DIR) $(CXX_FLAGS) $< -o $@

# $(TEST_EXE): bin/test_%: bin/test_%.o libScenGen.a
# 	g++ $^ $(LIB) $(LDFLAGS) -o $@

# .PHONY: test
# test: $(TEST_EXE)
# 	$(foreach var, $^, ./$(var);)

###############
## Utilities ##
###############

.PHONY: clean
clean:
	rm -rf *.o bin/* $(TEST_EXE)  *~ \#*\# *.pyc include/*~ src/*~ $(TEST_DIR)/*~
