CC = g++
LDFLAGS = `pkg-config --libs ibsimu-1.0.6dev`
CXXFLAGS = -Wall -g `pkg-config --cflags ibsimu-1.0.6dev`

simulation_3d: simulation_3d.o
	$(CC) -o simulation_3d simulation_3d.o $(LDFLAGS)

simulation_3d.o: simulation_3d.cpp
	$(CC) -c -o simulation_3d.o simulation_3d.cpp $(CXXFLAGS)

clean:
	$(RM) *~ *.o simulation_3d
