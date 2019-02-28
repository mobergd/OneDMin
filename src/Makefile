#Makefile for dint

#set POT here or use 'gmake POT=potname'
TINKER=0
POT=molpro
#POT=tbplusexp6all
#POT=rohM
#TINKER=1
METH=auto1dmin

EXEDIR = ./
TARG = $(EXEDIR)/auto1dmin.x
OBJDIR = obj
OBJDIR_OPT = obj

  FC = gfortran
  CC = gfortran
  LD = gfortran
  OPT_FFLAGS = -O3
#  FTPPFLAGS = -I../../sprng/SRC -I.. -c
  FTPPFLAGS = -I.. -c
  FFLAGS = -c
  FFLAGS = -c
  CFLAGS = -c
  COMPLIB = -llapack -lblas
#  RNGLIB = -lmlfg
  LIBS = $(COMPLIB) $(RNGLIB)
#  LDFLAGS = -O3 -g -L../../sprng/lib -L/usr/lib64 -L/home/ajasper/lib
  LDFLAGS = -O3 -g -L/usr/lib64 -L/home/ajasper/lib

VPATH = ..

SRC =	$(POT).f \
	spin.f \
	mass.f \
	$(METH).F 

OBJ = $(patsubst %.f,%.o,$(SRC))
OBJ := $(patsubst %.F,%.o,$(OBJ))
OBJ := $(patsubst %.cc,%.o,$(OBJ))

EXEDIR:=../$(EXEDIR)
MAKEOPT = -f ../Makefile
all: yes_opt

no_opt:
	$(MAKE) $(MAKEOPT) -C obj targ
yes_opt:
	$(MAKE) $(MAKEOPT) -C obj EXESUFFIX=.opt OPT="$(OPT_FFLAGS)" targ

targ: $(TARG)

# dependencies

COMPILE = $(FC) $(FFLAGS) $(OPT)
COMPILE = $(CC) $(CFLAGS) $(OPT)

$(TARG): $(OBJ)
	$(LD) $(LDFLAGS) $(OBJ) $(LIBS) -o $(TARG)

%.o : %.f
	$(COMPILE) $< -o $@

clean:
	-rm -f obj*/* 

realclean:
	-rm -f obj*/* ./auto1dmin.x

