CFLAGS     = -Wall
CPPFLAGS   =
LIBFILES   = 
TARGET     = hype
SRC        = $(wildcard *.c)
OBJ        = $(SRC:.c=.o)
CLEANFILES = $(TARGET)
LOCDIR     = $(CURDIR)

#PETSC_DIR=/home/magu/library/petsc

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

all: $(TARGET)

$(TARGET): $(OBJ) 
	-${CLINKER} -o $(TARGET) $(OBJ) ${PETSC_LIB}
	$(RM) $(OBJ)

distclean :
	$(RM) $(TARGET)

.PHONY: clean_kirchhoff

clean_kirchhoff:
	@rm -r kirchhoff/input/top/* kirchhoff/input/bottom/* kirchhoff/input/curved/* kirchhoff/p_weno.dat kirchhoff/input/t.dat

