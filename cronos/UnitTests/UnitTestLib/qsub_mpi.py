def get_qsub_mpi(catDir, catName, binaryPath, nproc):
    return '''#!/bin/bash
    
NAME={1}
DATADIR={0}

source /etc/profile.d/modules.sh
module load openmpi/intel/gcc-4.8.2/1.8.2
module load szip/intel/gcc-4.8.2/2.1
module load hdf5/intel/gcc-4.8.2/1.8.13

MPIRUN=mpirun
PARALLEL_PROGRAM={2}
PARALLEL_PROGRAM_OPTIONS="$DATADIR $NAME"

#$ -cwd
#$ -N {1}
#$ -pe openmpi {3}
#$ -m ea
#$ -o go
#$ -e go
#$ -q all.q

##-bind-to core:overload-allowed
##-bind-to numa
# -report-bindings

echo "Using $NSLOTS slots to run job" > $DATADIR/$NAME.out

# setenv pname
echo "$MPIRUN -np $NSLOTS -v $PARALLEL_PROGRAM $PARALLEL_PROGRAM_OPTIONS"  >>$DATADIR/$NAME.out

$MPIRUN -np $NSLOTS -v $PARALLEL_PROGRAM $PARALLEL_PROGRAM_OPTIONS >> $DATADIR/$NAME.out 2>> $DATADIR/$NAME.out
'''.format(catDir, catName, binaryPath, nproc)


if __name__ == '__main__':
    print(get_qsub_mpi('testDir', 'testName', '.test/testBin', 3))