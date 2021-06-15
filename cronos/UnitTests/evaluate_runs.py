import h5py
import os
import numpy as np
import matplotlib.pyplot as plt


class EvalTests(object):
	def __init__(self):

		self.rootDir = os.getcwd()

		# check which reference files are available
		self.refDir = os.path.join(self.rootDir, 'References')
		self.refFilesNames = [thing for thing in os.listdir(self.refDir)
							  if os.path.isfile(os.path.join(self.refDir, thing))
							  and not thing.startswith('.')
							  and thing.endswith('.h5')
							  ]
		self.refNames = ['_'.join(refFile.split('_')[:-2]) for refFile in self.refFilesNames]
		self.refFiles = [os.path.join(self.refDir, file) for file in self.refFilesNames]

		# check which test are available
		self.testRootDir = os.path.join(self.rootDir, 'tests')
		self.testDirs = [os.path.join(self.testRootDir, dir) for dir in os.listdir(self.testRootDir)
						 if os.path.isdir(os.path.join(self.testRootDir, dir))
						 and not dir.startswith('.')
						 ]
		self.testNames = []
		self.testResultDirs = []
		for testDir in self.testDirs:
			resultDir = [dir for dir in os.listdir(testDir)
						 if os.path.isdir(os.path.join(testDir, dir))
						 and dir.endswith('_float')
						 ]

			testNames = [dir.replace('_float', '').replace('.modded', '') for dir in resultDir]
			resultDir = [os.path.join(testDir, dir) for dir in resultDir]

			for name, dir in zip(testNames, resultDir):
				self.testNames.append(name)
				self.testResultDirs.append(os.path.join(testDir, dir))

		self.comp_dset_name = {}

		self.comp_dset_name["BrioWuHll"] = "v_x"
		self.comp_dset_name["BrioWuHlld"] = "v_x"
		self.comp_dset_name["MagRotor2DCartHll"] = "v_x"
		self.comp_dset_name["MagRotor2DCartHlld"] = "v_x"
		self.comp_dset_name["MagRotor2DPolarHlld"] = "v_x"
		self.comp_dset_name["OrszagTangVortexHll"] = "rho"
		self.comp_dset_name["OrszagTangVortexHlld"] = "rho"
		self.comp_dset_name["SedovCartHll"] = "v_x"
		self.comp_dset_name["SedovCartHllc"] = "v_x"
		self.comp_dset_name["SedovSphericalHll"] = "v_x"
		self.comp_dset_name["SedovSphericalHllc"] = "v_x"
		self.comp_dset_name["SodToro5hll"] = "rho"
		self.comp_dset_name["SodToro5hllc"] = "rho"

		# Dictionary used to store the result (passed/failed) for every test
		# for later printing
		self.results = {}

	def get_TimeCat(self, test_name):
		# read output time from cat file
		my_cat = test_name + ".cat"
		catFile = open(my_cat, "r")
		while True:
			line = catFile.readline()
			if not line:
				break
			if ((line.find("t_end") > -1)):
				posSign = line.find("=")
				posDot = line.find(".", posSign)
				posEnd = line.find(" ", posDot)
				t_end = float(line[posSign + 1:posEnd])
				break
		# print(" t_end:",t_end)
		return t_end

	def get_TimeFile(self, filePath):
		time = None
		if os.path.exists(filePath):
			h5in = h5py.File(filePath, 'r')
			h5group = h5in['Data']
			time = h5group.attrs['time']

		return time

	def get_KeysFile(self, filePath):
		# get all dataset names from given reference file
		h5in = h5py.File(filePath, 'r')
		h5group = h5in['fluid00']
		datasetNames = list(h5group.keys())
		# print("Names:",datasetNames)

		return datasetNames

	def get_TestFile(self, testResultsDir, tRef):

		# check if directory exists
		has_dir = os.path.isdir(testResultsDir)

		# print("Entering dir: ", testResultsDir)
		resultsFile = None

		if has_dir:
			for fileName in os.listdir(testResultsDir):
				file = os.path.join(testResultsDir, fileName)

				# get timestamp of the result
				time = self.get_TimeFile(file)

				if time is None:
					continue

				if abs(time - tRef) < 1.e-6:
					print(" The results file is ", os.path.relpath(file, self.rootDir))
					print(" -> Time:", time, "- Delta:", time - tRef)
					resultsFile = file
					break

		return resultsFile

	def comp_dset(self, refFile, testFile, dset_name):
		# Compare given dataset and write error

		# read data from reference file
		h5f_ref = h5py.File(refFile, 'r')
		h5g_ref = h5f_ref["fluid00"]
		dset_ref = h5g_ref[dset_name]

		if dset_ref.shape[0] == 1:
			if dset_ref.shape[1] == 1:
				data_ref = dset_ref[0, 0, :]
			else:
				data_ref = dset_ref[0, :, :]
		else:
			data_ref = dset_ref[:, :, :]

		rms_ref = np.sqrt(np.mean(np.square(data_ref)))
		# print("Data(ref): ",data_ref.shape)
		# print("RMS: ",rms_ref)

		# read data from test file
		h5f_test = h5py.File(testFile, 'r')
		h5g_test = h5f_test["fluid00"]
		dset_test = h5g_test[dset_name]
		#    print(dset_ref.shape)

		if dset_test.shape[0] == 1:
			if dset_test.shape[1] == 1:
				data_test = dset_test[0, 0, :]
			else:
				data_test = dset_test[0, :, :]
		else:
			data_test = dset_test[:, :, :]

		# Compute difference of data sets
		diff = data_ref - data_test
		rms_diff = np.sqrt(np.mean(np.square(diff)))
		# print("Data(test): ",data_test.shape)
		# print("Data(diff): ",diff.shape)

		if (rms_ref < 1.e-9 and rms_diff < 1.e-9):
			l2_err = 0.
			linf_err = 0.
		else:
			l2_err = rms_diff / rms_ref
			# print("l2 deviation: ",l2err)

			max_diff = np.amax(np.absolute(diff))
			linf_err = max_diff / rms_ref
			# print("linf deviation: ",linf_err)

		print(" ", dset_name, "->", "l2 deviation:", "%g" % l2_err, "linf deviation:", "%g" % linf_err)
		return l2_err, linf_err

	def plot_dset(self, ref_name, test_name, store_name, dset_name):
		# Compare given dataset and write error
		# read data from reference file
		h5f_ref = h5py.File(ref_name, 'r')

		h5g_ref = h5f_ref["fluid00"]
		dset_ref = h5g_ref[dset_name]

		# get size
		shape = dset_ref.shape
		Nx = shape[2]
		Ny = shape[1]
		Nz = shape[0]

		data_ref = dset_ref[Nz / 2, Ny / 2, :]

		# get grid
		h5g_grid = h5f_ref["Data"]
		dx = 0.
		dx = h5g_grid.attrs["dx"]
		xb = 0.
		xb = h5g_grid.attrs["xmin"]

		# check whether a nonlinear grid has been used:
		nonLin = "NonlinGrid" in h5g_grid
		if (nonLin):
			h5GridGroup = h5group["NonlinGrid"]

		xPosC = np.array(range(Nx), dtype='f')

		if (nonLin):
			dataset = h5GridGroup["xCentres"]
			self.xPosC = dataset[:]
		else:
			for ix in range(Nx):
				xPosC[ix] = xb[0] + dx[0] * ix

		# read data from test file
		h5f_test = h5py.File(test_name, 'r')
		h5g_test = h5f_test["fluid00"]
		dset_test = h5g_test[dset_name]

		data_test = dset_test[Nz / 2, Ny / 2, :]
		# print("Data(test): ",data_test.shape)


		xMin = xPosC.min()
		xMax = xPosC.max()
		yMin = 0.9 * data_ref.min()
		yMax = 1.1 * data_ref.max()

		# Plot both datasets
		fig = plt.figure()
		plt.axis([xMin, xMax, yMin, yMax])
		ax = fig.add_subplot(111)

		ax.plot(xPosC, data_ref, linestyle="-", linewidth=1.5, color="k")
		ax.plot(xPosC, data_test, linestyle="", marker='o', color="r")

		figName = store_name + ".pdf"
		print(" Saving: ", figName)

		fig.savefig(figName, transparent=True)


	def executeTest(self, testName, testResultDir):
		# find the corresponding reference file
		try:
			refIndex = self.refNames.index(testName)
		except:
			print(" Reference file not found for %s" % testName)
			status = "N/A"
			return status

		refFile = self.refFiles[refIndex]
		refTime = self.get_TimeFile(refFile)

		print(" The reference file is ", os.path.relpath(refFile, self.rootDir))

		# print(" Now we look at:",self.test_names[iTest])

		if refTime is not None:
			# get correct test file to compare to
			testResultFile = self.get_TestFile(testResultDir, refTime)

			if testResultFile is not None:
				dset_names = self.get_KeysFile(refFile)

				l2_err = 0.
				linf_err = 0.

				print("")
				for dset_name in dset_names:
					# print(dset_name)
					# compare given dataste
					errs = self.comp_dset(refFile, testResultFile, dset_name)
					l2_err = max(l2_err, errs[0])
					linf_err = max(linf_err, errs[1])

				# store pdf for comparison
				print("")
				self.plot_dset(refFile, testResultFile, testName, self.comp_dset_name.get(testName, 'rho'))

				if l2_err > 1.e-6 or linf_err > 1.e-6:
					status = "FAILED"
				else:
					status = "PASSED"

				print("")
				print(" l2 deviation:", "%g" % l2_err, "linf deviation:", "%g" % linf_err)
			else:
				print(" Test file not found ")
				status = "FAILED"
		else:
			print(" Ref file not found ")
			status = "FAILED"

		return status


	def execute(self):
		for testName, testResultDir in zip(self.testNames, self.testResultDirs):
			print("----------------------------------------------------------------")
			print("")
			print(" >>>> Comparing results for test ", testName)
			print("")

			status = self.executeTest(testName, testResultDir)				
			self.results[testName] = status

			print("    **** Test", testName, status, "****")
			print("")
			print("----------------------------------------------------------------")

	def print_summary(self):
		print("----------------------------------------------------------------")
		print("")
		print("Summary:")
		print("")
		for testName, status in self.results.items():
			print("  {:<35} {:<10}".format(testName, status) )

		print("----------------------------------------------------------------")


if __name__ == "__main__":
	Tester = EvalTests()
	Tester.execute()
	Tester.print_summary()
