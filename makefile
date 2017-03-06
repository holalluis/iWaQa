#
#  Makefile to build an executable from all source 
#  files in SDIR and its 1st level subdirectories
#  specified in SUBDIRS. 
#  All .cpp files are compiled, .o files are put 
#  into the ODIR directory (situated on the same
#  level as SDIR). 
#  The executable is placed in the project root.  


#######################################################
#            PROJECT SPECIFIC SETTINGS
#######################################################

# main target name
SERVEROUT = server
CLIENTOUT = client
LIBRARYOUT = libmodel

TXMLFILES = tinystr tinyxml tinyxmlerror tinyxmlparser
SERVERFILES = setup datatable complink modelfactory solver evaluator evaluatormethod particleswarm server main sampleutils biasmatrices seriesinterface filter script $(TXMLFILES)
CLIENTFILES = client
LIBRARYFILES = model mathutils lsodaintegrator

# compiler 
CC = g++
#for debug
# -O2
#-msse3
CFLAGS = -O2 -std=c++0x

# source directory
SDIR = src
ODIR = build
LDIR = lib
IDIR = include
MDIR = models
PDIR = modelsrc

# additional 1st level subdirectories in SDIR to compile
SUBDIRS = tinyxml

#######################################################
#              DYNAMIC DEBUG BUILDS
#######################################################

ifneq ($(DEBUG), )
	CFLAGS = -g
endif

#######################################################
#            FILE LISTS, OS DETECTION
#######################################################

# specify all source paths
DIRPATHS = $(SDIR) $(foreach dir, $(SUBDIRS), $(SDIR)/$(dir))
FILES = $(foreach dir,$(DIRPATHS),$(wildcard $(dir)/*.cpp)) 

# specific lists of .o files in the build directory
SERVEROBJS = $(foreach file,$(SERVERFILES),$(ODIR)/$(file).o) 
CLIENTOBJS = $(foreach file,$(CLIENTFILES),$(ODIR)/$(file).o)
LIBRARYOBJS = $(foreach file,$(LIBRARYFILES),$(ODIR)/$(file).o)

SERVERLFLAGS = -L"$(LDIR)" -lmodel
SOCKLFLAGS = 
LIBCFLAGS = 
PLATFORMLFLAGS = 

#######################################################
#                   OS DETECTION
#######################################################

# flag indicating a Windows environment
ISWIN =

# OS detection process
OSNAME = $(shell uname -s)
	
ifneq ($(OSNAME), )
	ifneq ($(findstring MINGW,$(OSNAME)), ) 
		ISWIN = 1
	endif
	ifneq ($(findstring CYGWIN,$(OSNAME)), )
		ISWIN = 1
	endif
else
	ISWIN = 1
	OSNAME = Windows
endif

OSID = $(OSNAME)
DLLEXT = so

#Windows specific general settings
ifeq ($(ISWIN), 1)
	DLLEXT = dll
	OSID = win
endif

# OSX specific general settings
ifneq ($(findstring Darwin,$(OSNAME)), )
	OSID = mac
endif

#Linux specific general settings
ifneq ($(findstring Linux,$(OSNAME)), )
	OSID = linux
endif

#######################################################
#    PRESERVE 32 BIT COMPILATION ON 64 BIT SYSTEMS
#######################################################

ifneq ($(ISWIN), 1)
	
	LBITS := $(shell getconf LONG_BIT)
	
	ifeq ($(LBITS),64)
   		# do 64 bit stuff here, like set some CFLAGS
		
		#CFLAGS = $(CGLAGS) -m32
	
	endif

endif
#######################################################
#                OS SPECIFIC COMMANDS
#######################################################

# defaults for unix-like OS
DELCMDO = rm -f $(ODIR)/*.o
DELCMDMODELS = rm -f $(MDIR)/*_$(OSID).$(DLLEXT)
DELCMDEXE = rm -f $(SERVEROUT)
DELCMDCLIENT = rm -f $(CLIENTOUT)
DELCMDLIB = rm -f $(LDIR)/$(LIBRARYOUT).a
DELCMDHDR = rm -f $(IDIR)/model.h $(IDIR)/lsodaintegrator.h
DELCMDPLUGINO = rm -f $(PLUGINOBJS)
CPCMDHEADER = cp $(SDIR)/model.h $(IDIR)/model.h
CPCMDHEADER2 = cp $(SDIR)/lsodaintegrator.h $(IDIR)/lsodaintegrator.h
DEFFILENAME = $(PDIR)/$*.def
DLLFILENAME = $(MDIR)/$*_$(OSID).$(DLLEXT)
LISTPLUGINSCMD = ls -dl $(PDIR)/*/ | awk '{print $$9}'

STRIPCMD = 

# Windows specific settings
ifeq ($(ISWIN), 1)
	DELCMDO = del $(ODIR)\*.o
	DELCMDEXE = del $(SERVEROUT).exe
	DELCMDMODELS = del $(MDIR)\*_$(OSID).$(DLLEXT)
	DELCMDCLIENT = del $(CLIENTOUT).exe
	SOCKLFLAGS = -lws2_32
	PLATFORMLFLAGS = -static-libgcc -static-libstdc++
	DELCMDLIB = del $(LDIR)\$(LIBRARYOUT).a
	DELCMDHDR = del $(IDIR)\model.h $(IDIR)\lsodaintegrator.h
	DELCMDPLUGINO = $(subst /,\,del $(PLUGINOBJS))
	CPCMDHEADER = copy $(SDIR)\model.h $(IDIR)\model.h
	CPCMDHEADER2 = copy $(SDIR)\lsodaintegrator.h $(IDIR)\lsodaintegrator.h
	PLUGINCFLAGS = -shared $(CFLAGS) -Wl,--kill-at,--output-def,$(DEFFILENAME)
	LISTPLUGINSCMD = dir $(PDIR) /AD /B
endif

# OSX specific dynamic library creation
ifneq ($(findstring Darwin,$(OSNAME)), )
	PLUGINCFLAGS = $(CFLAGS) -dynamiclib -exported_symbols_list interface_protocol -Wno-return-type-c-linkage
	#was above: -arch i386 
	#mimic def file creation + filter out all exported but undefined symbols
	DEFCREATECMD = nm -gm $(DLLFILENAME) | grep -v "undefined" > $(DEFFILENAME)
	#PLATFORMLFLAGS = -m32
	CC = g++
endif

# Linux specific dynamic library creation
ifneq ($(findstring Linux,$(OSNAME)), )
	PLUGINCFLAGS = $(CFLAGS) -fPIC -shared -Wl,-soname,$(DLLFILENAME) -o $(DLLFILENAME) -lc -s
	#mimic def file creation + filter out all exported but undefined symbols
	STRIPCMD = 
	DEFCREATECMD = sed '/\#/d' interface_protocol | sed -e 's/_//g' > $(DEFFILENAME)
	SOCKLFLAGS = -ldl
	LIBCFLAGS = -fPIC
endif

#######################################################
#          DYNAMIC TARGETS FOR ALL PLUGINS
#######################################################

ifeq ($(ISWIN), 1)
	PLUGINNAMES = $(shell $(LISTPLUGINSCMD))
else
# need some post-processing
	PLUGINFULLPATHS = $(shell $(LISTPLUGINSCMD))
	PLUGINNAMES = $(foreach dir,$(PLUGINFULLPATHS),$(subst /,,$(subst $(PDIR)/,,$(dir))))
endif

# a virtual filename for each plugin target
PLUGINDLLNAMES  = $(foreach plugin,$(PLUGINNAMES),$(plugin).plugin)

#######################################################
#                      RULES
#######################################################

.PHONY : printosid clean cleanbuild all plugin

# build all
all: printosid $(SERVEROUT) $(CLIENTOUT) $(PLUGINDLLNAMES)

#identify operating system
osid: 
		@echo OS detection: $(OSID) '('$(OSNAME)')'.

# clean everything
clean: 
		$(DELCMDO) 
		$(DELCMDEXE)
		$(DELCMDCLIENT)
		$(DELCMDLIB)
		$(DELCMDHDR)
		$(DELCMDMODELS)

# clean intermediate objects
cleanbuild:
		$(DELCMDO)
		
# compile any .o files
$(ODIR)/%.o : $(filter %.cpp, $(FILES))
		@echo Compiling $@
		@$(CC) -c $(filter %/$(addsuffix .cpp, $(basename $(notdir $@))), $(FILES)) -o $@ $(CFLAGS) $(if $(findstring $(basename $(notdir $@)),$(LIBRARYOBJS)),$(LIBCFLAGS))

# static library of iWQModel
$(LIBRARYOUT): $(LIBRARYOBJS)
		@echo Making $(LIBRARYOUT)
		@ar rcs $(LDIR)/$(LIBRARYOUT).a $(LIBRARYOBJS)
		@$(CPCMDHEADER)
		@$(CPCMDHEADER2)

# model server		
$(SERVEROUT): $(LIBRARYOUT) $(SERVEROBJS) 
		@echo Making $(SERVEROUT)
		@$(CC) $(SERVEROBJS) -o $(SERVEROUT) $(SERVERLFLAGS) $(SOCKLFLAGS) $(PLATFORMLFLAGS)
		
# model client
$(CLIENTOUT): $(CLIENTOBJS)
		@echo Making $(CLIENTOUT)
		@$(CC) $(CLIENTOBJS) -o $(CLIENTOUT) $(SOCKLFLAGS) $(PLATFORMLFLAGS)

# model plugins 		
plugin: $(NAME).plugin

%.plugin: $(LDIR)/$(LIBRARYOUT).a
		@echo Making model: $*
		@$(CC) $(PDIR)/$*/*.cpp -o $(MDIR)/$*_$(OSID).$(DLLEXT) -I"$(IDIR)" -DIWQ_MODEL_NAME=$* $(PLUGINCFLAGS) $(SERVERLFLAGS) $(PLATFORMLFLAGS)
		@$(STRIPCMD)
		@$(DEFCREATECMD)
		
