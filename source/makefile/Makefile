SHELL := /bin/bash

hn ?= $(shell hostname)

ifneq "$(findstring kraken,$(shell hostname))" ""
CPPFLAGS += -D KRAKEN
include make.kraken
endif

ifneq "$(findstring mpp,$(shell hostname))" ""
CPPFLAGS += -D LRZ_CLUSTER
include make.lrz
endif

ifneq "$(findstring anslab,$(shell hostname))" ""
CPPFLAGS += -D ANSLAB
include make.anslab
endif


ifeq "$(helmholtz)" "hypre"
	CPPFLAGS += -D_MULTIGRID_
endif

ifeq "$(config)" "production"
CPPFLAGS += -DNDEBUG
CPPFLAGS += -D_FMMSILENT
endif

#set parameters
dim ?= 3 
bpd ?= 8
bs ?= 14 
#bs 18, bpd = 8 for lower resolution, bs=16 
CPPFLAGS += -D_DIM=$(dim) -D_BLOCKSIZE_=$(bs) -D_BPD_=$(bpd) 
CPPFLAGS += -D HGG

#CPPFLAGS += -D Propagation
#CPPFLAGS += -D HGG_UQ_MALA
#CPPFLAGS += -D Visualisation
#CPPFLAGS += -D Necrosis
#CPPFLAGS += -D VP


ifeq "$(dim)" "2"
	bsz ?= 1
	CPPFLAGS += -D_BLOCKSIZE_Z_=$(bsz)
endif

ifeq "$(bpd)" "16"
	ml ?= 4
	CPPFLAGS += -D _MAXLEVEL_=$(ml)
endif

ifeq "$(bpd)" "8"
        ml ?= 3
        CPPFLAGS += -D _MAXLEVEL_=$(ml)
endif

	
LIBS += \
	-L$(TBB_LIB_DIR) \
	-ltbb \
	-ltbbmalloc \
	-L$(VTK_LIB_DIR) \
	-lvtkViews \
	-lvtkInfovis \
	-lvtkWidgets \
	-lvtkHybrid \
	-lvtkRendering \
	-lvtkGraphics \
	-lvtkverdict \
	-lvtkImaging \
	-lvtkftgl \
	-lvtkfreetype \
	-lvtkIO \
	-lvtkFiltering \
	-lvtkCommon \
	-lm \
	-lvtkDICOMParser \
	-lvtkmetaio \
	-lvtksqlite \
	-lvtkpng \
	-lvtktiff \
	-lvtkjpeg \
	-lvtkexpat \
	-lvtksys \
	-lvtkexoIIc \
	-lvtkNetCDF \
	-lvtklibxml2 \
	-lvtkzlib \
	-lpthread \
	-ldl \
	$(OPENMP_FLAG)
	
	
##################

ifeq "$(helmholtz)" "hypre"
	CC=mpic++
	LD=$(CC)
	CPPFLAGS += -I$(HYPRE_INC_DIR)
	LIBS += -L$(HYPRE_LIB_DIR) -lHYPRE
endif	


VPATH := ../Glioma/ ../MRAGcore/
.DEFAULT_GOAL := brain

########################################################################
# List of all programs to be compiled
########################################################################
	


#HGG PART
OBJECTS = \
Glioma_main.o \
Glioma_dat2VP.o\
Glioma_RAT_UQ.o\
Glioma_RAT_preprocessing.o\
MRAGBoundaryBlockInfo.o \
MRAGProfiler.o \
MRAGWavelets_StaticData.o \
Test.o\

#all: brain

brain: $(OBJECTS)
	$(CC) $^ -o $@ $(LIBS)
	@echo done

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o 
	rm -f brain 
	echo $(OBJECTS)

