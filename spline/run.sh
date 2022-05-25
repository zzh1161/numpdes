#!/bin/bash

./bin/exercise_B > ./output/exercise_B.txt
echo "Exercise B is finished!"

./bin/exercise_CD > ./output/exercise_CD.txt
echo "Exercise C & D are finished!"

./bin/exercise_E > ./output/exercise_E.txt
echo "Exercise E is finished!"

./bin/exercise_F > ./output/exercise_F.txt
echo "Exercise F is finished!"

if [[ ! -f "./bin/test" ]]; then
    echo "Please \"make test\" first!"
else
    ./bin/test > ./output/test.txt
    echo "test is finished!"
fi