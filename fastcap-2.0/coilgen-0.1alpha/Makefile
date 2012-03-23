# Makefile
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.


SHELL = /bin/sh

CC      = gcc
CFLAGS  = -O2 -g -Wall -pedantic # -DDEBUG 

LD = gcc
LDFLAGS = -lm

# This is where executables, man-pages and support files are installed
BINDIR = /usr/local/bin
MANDIR = /usr/local/man/man1

############### No configuration parameters below this line ###############

# fasthenry input file generation for coils systems

SRCS = coilgen.c
SOURCES = Makefile $(SRCS) # plus any .h files

OBJ = coilgen.o

DOCFILS = COPYING README TODO coilrc

all: coilgen

coilgen: coilgen.o

.PHONY : all depend clean dist

depend:
	makedepend -- $(CFLAGS) -- $(SRCS)

clean: 
	-rm coilgen *~ core $(OBJ)

dist: $(SOURCES) $(DOCFILS)
	mkdir coilgen-0.1alpha
	ln $(SOURCES) $(DOCFILS) coilgen-0.1alpha
	tar cf - coilgen-0.1alpha | gzip -best > coilgen-0.1alpha.tar.gz
	rm -rf coilgen-0.1alpha

# DO NOT DELETE THIS LINE -- make  depend  depends  on  it.
