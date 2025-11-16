#
#	Makefile for all programs in this directory
#
#	implements:
#		all (default), create, analyze, graphlet, clean
#

# Global configuration for SNAP makefiles
GLIB = glib-core
SNAP = snap-core
GLIBADV = glib-adv
SNAPADV = snap-adv
SNAPEXP = snap-exp

CGLIB = ../snap-master/$(GLIB)
CSNAP = ../snap-master/$(SNAP)

EXGLIB = ../snap-master/$(GLIB)
EXSNAP = ../snap-master/$(SNAP)
EXGLIBADV = ../snap-master/$(GLIBADV)
EXSNAPADV = ../snap-master/$(SNAPADV)
EXSNAPEXP = ../snap-master/$(SNAPEXP)

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	CC = g++
	CXXFLAGS += -std=c++0x -Wall
	CXXFLAGS += -O3 -DNDEBUG -fopenmp
	LDFLAGS +=
	LIBS += -lrt -lm
endif

ANALYZE = analyze
GRAPHLET = graphlet
RES = resilience

all: $(ANALYZE) $(GRAPHLET) $(RES)

# COMPILE

$(ANALYZE): $(ANALYZE).cpp $(EXSNAP)/Snap.o 
	$(CC) $(CXXFLAGS) -o $(ANALYZE) $(ANALYZE).cpp $(EXSNAP)/Snap.o -I$(EXSNAP) -I$(EXSNAPADV) -I$(EXGLIB) -I$(EXSNAPEXP) $(LDFLAGS) $(LIBS)

$(RES): $(RES).cpp $(EXSNAP)/Snap.o 
	$(CC) $(CXXFLAGS) -o $(RES) $(RES).cpp $(EXSNAP)/Snap.o -I$(EXSNAP) -I$(EXSNAPADV) -I$(EXGLIB) -I$(EXSNAPEXP) $(LDFLAGS) $(LIBS)

$(GRAPHLET): $(GRAPHLET).cpp $(EXSNAP)/Snap.o 
	$(CC) $(CXXFLAGS) -o $(GRAPHLET) $(GRAPHLET).cpp $(EXSNAP)/Snap.o -I$(EXSNAP) -I$(EXSNAPADV) -I$(EXGLIB) -I$(EXSNAPEXP) $(LDFLAGS) $(LIBS)

$(EXSNAP)/Snap.o: 
	make -C $(EXSNAP)

clean:
	rm -f *.o  $(ANALYZE) $(GRAPHLET)
