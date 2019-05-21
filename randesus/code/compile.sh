#
# NAUTY
#
echo Compiling Nauty
cd nauty
gcc -c -O3 *.c
cd ..

#
# MAIN PROGRAM
#
echo Compiling Randesu

g++ -O3 randesu.cpp output.cpp init.cpp graph64.cpp random.cpp maingraph.cpp nauty/gtools.o nauty/naugraph.o nauty/nautil.o  nauty/nautinv.o nauty/nauty.o  -o randesu
