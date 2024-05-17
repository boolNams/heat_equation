heat: main.o func.o scheme.o
	g++ -o heat main.o scheme.o func.o

main.o: main.cpp HEAT.h
	g++ -c main.cpp
	
func.o: func.cpp HEAT.h
	g++ -c func.cpp
	
scheme.o: scheme.cpp
	g++ -c scheme.cpp