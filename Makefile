all: 
	g++ main.cpp rbfcolor.cpp utils.cpp png/lodepng.cpp -o main -std=c++17

script:
	g++ main.cpp rbfcolor.cpp utils.cpp png/lodepng.cpp -o main -std=c++17 -Dscript