import os, sys
from UnitTestLib.UnitTest import UnitTest
import multiprocessing

if __name__ == "__main__":
	rootDir = os.getcwd()

	if len(sys.argv) > 1:
		unitTests = sys.argv[1:]
	else:
		testRootDir = os.path.join(rootDir, 'tests')
		unitTests = ['tests/' + thing for thing in os.listdir(testRootDir)
						if os.path.isdir(os.path.join(testRootDir, thing))
						and not thing.startswith('.')
					]

	for test in unitTests:
		print('================================================================')
		print('Test: ', test)
		print('................................................................')

		testPath = os.path.join(rootDir, test)

		unitTest = UnitTest(test, testPath)
		unitTest.get_setupFiles()
		build_ok = unitTest.build_all(multiprocessing.cpu_count())
		if not build_ok:
			continue

		unitTest.find_catFiles()

		# run the job
		#unitTest.qsub_serial()
		#unitTest.qsub_mpi()
		unitTest.qsub_auto()
		#unitTest.run_serial()
		#unitTest.run_mpi()
		#unitTest.run_auto()
