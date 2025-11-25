#!/bin/bash
gcc -o 2d_test 2d_test_main.c ../contact_full.c -lm -Wall
./2d_test
echo "2D Test Complete!"