FC = gfortran

TARGET = 2dU1.exe

SRC = src
BIN = bin

SOURCE = statistics.f90 pbc.f90 arrays.f90 parameters.f90 dynamics.f90 main.f90

OBJECT = $(patsubst %, $(BIN)/%, $(SOURCE:.f90=.o ) )

#FFLAGS = -Wall -Wextra -fcheck=all -O0 -J$(BIN) -I$(BIN)
FFLAGS = -O3 -J$(BIN) -I$(BIN)


$(BIN)/$(TARGET): $(OBJECT)
	$(FC) -o $@ $^

$(BIN)/%.o: $(SRC)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

.PHONY: help run clean

run:
	@echo input_parameters.par | time $(BIN)/$(TARGET)

clean:
	rm -f $(OBJECT) $(BIN)/$(TARGET)

help:
	@echo "src: $(SOURCE)"
