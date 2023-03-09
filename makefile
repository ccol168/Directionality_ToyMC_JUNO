LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

fai: Directionality_ToyMC.cxx
	g++ Directionality_ToyMC.cxx -o Directionality_ToyMC -lMinuit ${INCS} ${LIBS}
clean:
	rm *.o output
