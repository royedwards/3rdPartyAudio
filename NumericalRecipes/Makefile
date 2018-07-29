BASEDIR = ../..
LDIR = ${BASEDIR}/lib
OPTIM = -O 
INC = -I${BASEDIR}/include
CFLAGS = ${OPTIM} ${INC}
A = recipes
all:	lib$A.a
lib$A.a:	used
	build
	ar r lib$A.a *.o
	ranlib lib$A.a
mlint:	lib$A.a
	mlint $A ${INC} `cat used`
lint:
	lint ${INC} -bhux `cat used`
