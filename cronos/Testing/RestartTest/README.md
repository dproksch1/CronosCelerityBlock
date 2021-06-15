# Simple Restart Test

## Quickstart

* Compile the respective code version (serial, parallel)
* Run `runTestSerial.sh` or `runTestMpi.sh`, respectively
* Analyse the results with `evalTest.sh`

## Under the hood

* Write random data to hdf5 (`RestartTest_first.cat`)
* Restart from it and write it again (`RestartTest_second.cat`)
* Dump the produced float files and compare them (`evalTest.sh`)
