mkdir -p bin/intermediate && g++ -std=c++23 -I. tests/demo.cpp -o bin/demo
mkdir -p bin/intermediate && g++ -std=c++23 -I. tests/exact_tests.cpp -o bin/exact_tests
mkdir -p bin/intermediate && g++ -std=c++23 -I. tests/ergonomics_tests.cpp -o bin/ergonomics_tests
