# Author: Ryan Houlihan

SHELL = /bin/sh

# ------ FILES ------
SRC_MAT =	LatUtil.cpp LatPro.cpp GranPro.cpp GrassPro.cpp
INC_MAT =       LatPro.h GranPro.h LatUtil.h

SRCDIR=SRC
INCDIR=INCLUDE

SRC=SRC_MAT
INC=$(patsubst %,$(INCDIR)/%,$(INC_MAT))

FILES = $(SRC) $(INCDIR) $(LIBDIR) Makefile

# ------ DEFINITIONS ------
ODIR=../OBJ

_OBJ = $(SRC_MAT:.cpp=.o)
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


# ------ SETTINGS ------

CC =	        g++
CCFLAGS =       -g
ARCHIVE =	ar
ARCHFLAG =	-rc
DEPFLAGS =      -M
LINK =         	g++
LINKFLAGS =	-O
USRLIB =	-L$(LIBDIR)/BIN -lmath	
INCLUDE =	-I$(LIBDIR)/MATH
SYSLIB =
LIBDIR=../LIB

# ------ MAKE PROCEDURE ------

all: Math LatPro	

LatticeProducer: Math LatPro

# ------ COMPILE RULES ------

$(ODIR)/%.o:%.cpp $(INC)
	$(CC) $< -o $@ -c $(INCLUDE)

LatPro: $(OBJ)
	$(CC) -o $@ $^ $(USRLIB) 		

Math:
	cd $(LIBDIR)/MATH; make; cd ../../$(SRCDIR)

# ------ DEPENDENCIES ------

#DEPENDS = $(OBJ:.o=.d)

# ------ CLEAN ------

clean:	
	rm LatPro; 
	cd ../OBJ; rm *.o; cd ../SRC
	cd ../LIB/MATH; make clean; cd ../../SRC
tar:
	tar -cvf LatPro.tar $(FILES)
