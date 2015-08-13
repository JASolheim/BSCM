del *.o
del *.exe
cls
H:\JASolheim\MinGW\bin\g++.exe Spline.cpp ^
-Wall -c -O2 ^
-o Spline.o ^
-I"H:\JASolheim\EIGEN-~1\EIGEN-~1"

H:\JASolheim\MinGW\bin\g++.exe main.cpp ^
-Wall -c -O2 ^
-o main.o ^
-I"H:\JASolheim\EIGEN-~1\EIGEN-~1"

H:\JASolheim\MinGW\bin\g++.exe -o main.exe Spline.o main.o
