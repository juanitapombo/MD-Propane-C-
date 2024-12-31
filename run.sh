#!/bin/bash 

g++ -o program forces2.cpp iniconfig.cpp vverlet.cpp main_moldyn.cpp -stdlib=libc++
./program
