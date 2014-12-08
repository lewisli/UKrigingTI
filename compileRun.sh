cd Release
make clean
make -j4
cd ../python

python runKrigingNonStationary.py
python runKrigingNonStationaryAll.py
python runKrigingWalkerLake.py
python runKrigingWalkerLakeAll3.py
python runKrigingWalkerLakeAll5.py
python runKrigingWalkerLakeRaw.py
python runKrigingWalkerLakeRawAll3.py
python runKrigingWalkerLakeRawAll5.py
