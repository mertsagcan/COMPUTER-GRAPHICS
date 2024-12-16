all: main

main: main.cpp
	g++ main.cpp -o main -g -lglfw -lpthread -lX11 -ldl -lXrandr -lGLEW -lGL -lfreetype -DGL_SILENCE_DEPRECATION -DGLM_ENABLE_EXPERIMENTAL -I. -I/usr/include/freetype2

clean:
	rm -f main

