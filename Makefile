# Winter 2017
#
BANG    =       $(shell expr match `hostname` ccom-bang)
BANG-COMPUTE   =  $(shell expr match `hostname` compute)
BANG-BANG = $(shell echo $(BANG)\&$(BANG-COMPUTE))


include $(PUB)/Arch/arch.gnu-c++11.generic

OPTIMIZATION += -O3
CMD        = apf
THE_DATE    := $(shell date "+%y.%m.%d-%H.%M.%S-%Z")
THE_ARCHIVE = $(CMD)--$(THE_DATE)

# set vec=1 if you want to vectorize by hand
ifeq ($(vec),1)
C++FLAGS += -DSSE_VEC
C++FLAGS += -msse -msse2
endif

app:		apf

OBJECTS = apf.o solve.o Plotting.o cmdLine.o Report.o utils.o helper.o Timer.o

apf:	        $(OBJECTS) 
		$(C++LINK) $(LDFLAGS) -o $@ $(OBJECTS)  $(LDLIBS)

tgz:
	-tar cvfz $(THE_ARCHIVE).tgz Makefile *.cpp *.h
	-chmod 640 $(THE_ARCHIVE).tgz

.PHONY: clean
clean:	
	$(RM) *.o apf;
	$(RM) core;
