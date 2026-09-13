mkdir -p bin/intermediate
for file in tests/*.cpp
do
    filename="${file##*/}"
    cleanname="${filename%.*}"
    g++ -std=c++23 -I. $file -o bin/$cleanname
done
