# Makefile for creating our standalone Cython program
PYTHON := python
PYVERSION := $(shell $(PYTHON) -c "import sys; print(sys.version[:3])")
PYPREFIX := $(shell $(PYTHON) -c "import sys; print(sys.prefix)")

INCDIR := $(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_python_inc())")
PLATINCDIR := $(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_python_inc(plat_specific=True))")
LIBDIR1 := $(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
LIBDIR2 := $(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBPL'))")
PYLIB := $(shell $(PYTHON) -c "from distutils import sysconfig; print(sysconfig.get_config_var('LIBRARY')[3:-2])")

CC := $(shell $(PYTHON) -c "import distutils.sysconfig; print(distutils.sysconfig.get_config_var('CC'))")
LINKCC := $(shell $(PYTHON) -c "import distutils.sysconfig; print(distutils.sysconfig.get_config_var('LINKCC'))")
LINKFORSHARED := $(shell $(PYTHON) -c "import distutils.sysconfig; print(distutils.sysconfig.get_config_var('LINKFORSHARED'))")
LIBS := $(shell $(PYTHON) -c "import distutils.sysconfig; print(distutils.sysconfig.get_config_var('LIBS'))")
SYSLIBS :=  $(shell $(PYTHON) -c "import distutils.sysconfig; print(distutils.sysconfig.get_config_var('SYSLIBS'))")

paths:
	@echo "PYTHON=$(PYTHON)"
	@echo "PYVERSION=$(PYVERSION)"
	@echo "PYPREFIX=$(PYPREFIX)"
	@echo "INCDIR=$(INCDIR)"
	@echo "PLATINCDIR=$(PLATINCDIR)"
	@echo "LIBDIR1=$(LIBDIR1)"
	@echo "LIBDIR2=$(LIBDIR2)"
	@echo "PYLIB=$(PYLIB)"
	@echo "CC=$(CC)"
	@echo "LINKCC=$(LINKCC)"
	@echo "LINKFORSHARED=$(LINKFORSHARED)"
	@echo "LIBS=$(LIBS)"
	@echo "SYSLIBS=$(SYSLIBS)"


build:
	mkdir build

embed-cython: build
	cython fastqe/fastqe.py -o build/fastqe.cpp --embed --cplus --3 -Wextra --no-docstrings

compile: embed-cython
	$(CC) -c build/fastqe.cpp -I$(INCDIR) -I$(PLATINCDIR) -fPIC -o build/fastqe.o

embed: compile
	$(LINKCC) -o build/fastqe build/fastqe.o -L$(LIBDIR1) -L$(LIBDIR2) -l$(PYLIB) $(LIBS) $(SYSLIBS) $(LINKFORSHARED) -rpath $(LIBDIR1)

clean:
	rm build/*