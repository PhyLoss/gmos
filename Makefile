EXECFILE = gmos
EXECFILE64 = gmos64

DIRECTORY = Gmos
################################################

VERSION = 1.0
.PHONY : all
all : $(EXECFILE) $(EXECFILE64)

# 32 bit version
$(EXECFILE):
	cd Src/Gmos; make; cp gmos ../../
# 64 bit version
$(EXECFILE64):
	cd Src/Gmos; make; cp gmos64 ../../

###############################################
#
# Other Standard make rules
#
clean:
	rm $(EXECFILE) $(EXECFILE64)
	cd ExternSrc/DeepShallow/Src/; make clean;
	cd ExternSrc/DeepShallow64/Src; make clean;
	cd ExternSrc/LibDivSufSort; make clean;
	cd Src/Common/Interval; make clean;
	cd Src/Common/Sequence/; make clean;
	cd Src/Common/Util/; make clean;
	cd Src/Gmos; make clean;
tarfile:
	mkdir $(DIRECTORY)_$(VERSION)
	#cd Doc; make pdf 
	cp -r Doc $(DIRECTORY)_$(VERSION)
	cp -r Data $(DIRECTORY)_$(VERSION)
	cp -r Auxiliary $(DIRECTORY)_$(VERSION)
	cp ACKNOWLEDGMENTS COPYRIGHT.GPL.txt COPYRIGHT ChangeLog INSTALL README Makefile $(DIRECTORY)_$(VERSION)
	mkdir $(DIRECTORY)_$(VERSION)/ExternSrc $(DIRECTORY)_$(VERSION)/Lib \
	$(DIRECTORY)_$(VERSION)/Lib/DeepShallow $(DIRECTORY)_$(VERSION)/Lib/DeepShallow64 $(DIRECTORY)_$(VERSION)/Lib/LibDivSufSort
	mkdir $(DIRECTORY)_$(VERSION)/ExternSrc/DeepShallow $(DIRECTORY)_$(VERSION)/ExternSrc/DeepShallow64 $(DIRECTORY)_$(VERSION)/ExternSrc/LibDivSufSort 
	mkdir $(DIRECTORY)_$(VERSION)/ExternSrc/DeepShallow/Src $(DIRECTORY)_$(VERSION)/ExternSrc/DeepShallow64/Src \
	$(DIRECTORY)_$(VERSION)/ExternSrc/LibDivSufSort/examples $(DIRECTORY)_$(VERSION)/ExternSrc/LibDivSufSort/include $(DIRECTORY)_$(VERSION)/ExternSrc/LibDivSufSort/lib	
	cp ExternSrc/DeepShallow/COPYRIGHT.* ExternSrc/DeepShallow/README $(DIRECTORY)_$(VERSION)/ExternSrc/DeepShallow/
	cp -r ExternSrc/DeepShallow/Doc $(DIRECTORY)_$(VERSION)/ExternSrc/DeepShallow/
	cp ExternSrc/DeepShallow64/COPYRIGHT.* ExternSrc/DeepShallow64/README $(DIRECTORY)_$(VERSION)/ExternSrc/DeepShallow64/
	cp -r ExternSrc/DeepShallow64/Doc $(DIRECTORY)_$(VERSION)/ExternSrc/DeepShallow64/
	cp ExternSrc/DeepShallow/Src/Makefile ExternSrc/DeepShallow/Src/*.c ExternSrc/DeepShallow/Src/*.h $(DIRECTORY)_$(VERSION)/ExternSrc/DeepShallow/Src/
	cp ExternSrc/DeepShallow64/Src/Makefile ExternSrc/DeepShallow64/Src/*.c ExternSrc/DeepShallow64/Src/*.h $(DIRECTORY)_$(VERSION)/ExternSrc/DeepShallow64/Src/
	cp -r ExternSrc/LibDivSufSort/* $(DIRECTORY)_$(VERSION)/ExternSrc/LibDivSufSort	
	mkdir $(DIRECTORY)_$(VERSION)/Src $(DIRECTORY)_$(VERSION)/Src/Gmos \
	$(DIRECTORY)_$(VERSION)/Src/Common $(DIRECTORY)_$(VERSION)/Src/Common/Interval \
	$(DIRECTORY)_$(VERSION)/Src/Common/Sequence $(DIRECTORY)_$(VERSION)/Src/Common/Util
	cp Src/Common/*.h $(DIRECTORY)_$(VERSION)/Src/Common/
	cp Src/Common/Interval/Makefile Src/Common/Interval/*.c Src/Common/Interval/*.h $(DIRECTORY)_$(VERSION)/Src/Common/Interval
	cp Src/Common/Sequence/Makefile Src/Common/Sequence/*.c Src/Common/Sequence/*.h $(DIRECTORY)_$(VERSION)/Src/Common/Sequence
	cp Src/Common/Util/Makefile Src/Common/Util/*.c Src/Common/Util/*.h $(DIRECTORY)_$(VERSION)/Src/Common/Util
	cp Src/Gmos/Makefile Src/Gmos/*.c Src/Gmos/*.h $(DIRECTORY)_$(VERSION)/Src/Gmos
	tar cvzfh $(EXECFILE)_$(VERSION).tgz $(DIRECTORY)_$(VERSION)
	mv $(EXECFILE)_$(VERSION).tgz ../
	rm -r $(DIRECTORY)_$(VERSION)