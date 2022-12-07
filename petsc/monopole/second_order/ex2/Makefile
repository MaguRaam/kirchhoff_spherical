CFLAGS     = -Wall
CPPFLAGS   =
LIBFILES   = 
TARGET     = hype
SRC        = $(wildcard *.c)
OBJ        = $(SRC:.c=.o)
CLEANFILES = $(TARGET)
LOCDIR     = $(CURDIR)

#PETSC_DIR=/scratch/mecmagur/library/petsc

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test

all: $(TARGET)

$(TARGET): $(OBJ) 
	-${CLINKER} -o $(TARGET) $(OBJ) ${PETSC_LIB}
	$(RM) $(OBJ)

distclean :
	$(RM) $(TARGET)
