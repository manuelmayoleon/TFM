OBJECTS = gates_dani.o rutinas.o 
#MODULES = gates.mod 

FFLAGS = -g 

.PHONY: main.exe  clean

main: main.exe
	  ./main.exe

main.exe: $(MODULES) $(OBJECTS) 
		gfortran $(FFLAGS) $(OBJECTS) int_2d.f90 prueba_3.f90  -o  main.exe  
		
		

%.o : %.f90 
		gfortran $(FFLAGS)  -c   $< 



