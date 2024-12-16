make clean
make
for file in ../test_scenes/inputs/*
do
    echo $file
    ./raytracer $file
done
