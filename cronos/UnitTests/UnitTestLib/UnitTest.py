import os, sys, platform
from shutil import copyfile
import subprocess
from .util import *
from .qsub_mpi import *
from .qsub_serial import *


class UnitTest:
	def __init__(self, testName, testPath):
		self.testName = testName
		self.testPath = testPath
		self.rootPath = os.getcwd()

		relativeSetupPath = ''

		with open(os.path.join(self.testPath, 'setup_path'), 'r') as file:
			relativeSetupPath = file.readlines()[0]
			relativeSetupPath = relativeSetupPath.replace('\n', '')

		self.setupPath = os.path.join(self.testPath, relativeSetupPath)


	def modify_makefile(self):
		with open(os.path.join(self.testPath, "Makefile"),"r") as makefile:
			content = makefile.read()

		content = content.split("\n")

		for idx, line in enumerate(content):
			if "path_include" in line:
				break

		theline = content[idx].split(" ")
		theline[-1] = self.setupPath + '/' + theline[-1]
		content[idx] = " ".join(theline)
		content = "\n".join(content)

		with open(os.path.join(self.testPath, "Makefile"), "w") as makefile:
			content = makefile.write(content)

	def get_setupFiles(self):
		setupFiles = [file for file in os.listdir(self.setupPath) if os.path.isfile(os.path.join(self.setupPath, file))]

		# filter some unwanted files
		setupFiles = [name for name in setupFiles if name != 'constants.H' and not name.endswith('.cat')]

		print("Copying ", len(setupFiles), "setup files...")

		# copy base setup
		for file in setupFiles:
			src = os.path.join(self.setupPath, file)
			dst = os.path.join(self.testPath, file)
			copyfile(src, dst)

		self.modify_makefile()


	def get_compileDir(self):
		system = platform.system()
		machine = platform.machine()

		# make it compatible with the os_type of NumLib
		# no clue if that works in any case...
		machine = machine.replace('x86_64', 'amd64')
		machine = machine.replace('.86', '386')
		machine = machine.replace('/800', '')

		dirName = system + '-' + machine
		return os.path.join(self.testPath, dirName)

	def build_all(self, nCores=1, clean=True):
		# build stuff
		success = True
		success &= self.build_serial(nCores, clean)
		success &= self.build_mpi(nCores, False)
		return success

	def build_serial(self, nCores=1, clean=True):
		# build stuff
		print('Building Serial Version:', self.testName, '...')
		success = True
		with open(os.path.join(self.testPath, 'compile.out'), "w") as out:
			try:
				if clean:
					p = subprocess.check_call(['make', 'clean'], cwd=self.testPath, stdout=out, stderr=out)
				p = subprocess.check_call(['make', '-j' + str(nCores), 'proj'], cwd=self.testPath, stdout=out, stderr=out)
				print('  done!')
			except subprocess.CalledProcessError as e:
				print('  failed!')
				success = False

		if success:
			self.compilePath_serial = os.path.join(self.get_compileDir(), 'proj')

		return success

	def build_mpi(self, nCores=1, clean=True):
		# build stuff
		print('Building MPI Version:', self.testName, '...')
		success = True
		with open(os.path.join(self.testPath, 'compile_mpi.out'), "w") as out:
			try:
				if clean:
					p = subprocess.check_call(['make', 'clean'], cwd=self.testPath, stdout=out, stderr=out)
				p = subprocess.check_call(['make', '-j' + str(nCores), 'proj_MPI'], cwd=self.testPath, stdout=out, stderr=out)
				print('  done!')
			except subprocess.CalledProcessError as e:
				print('  failed!')
				success = False

		if success:
			self.compilePath_mpi = os.path.join(self.get_compileDir(), 'proj_MPI')

		return success


	def find_catFiles(self):
		# find cat files
		catFiles = [file for file in os.listdir(self.testPath)
					if file.endswith('.cat') and not file.endswith('.modded.cat')
					]
		self.catNames = [file[:-4] for file in catFiles]
		self.catFiles = [os.path.join(self.testPath, file) for file in catFiles]
		print('Cat Files found: ', len(self.catFiles), self.catNames)


	def get_nproc(self, catFile):
		# check out how many cores are needed
		nProcTotal = 1

		cat_content = None
		with open(catFile, 'r') as file:
			cat_content = file.readlines()

		nproclines = [line for line in cat_content if line.startswith('nproc')]
		nprocs = [line.split('=')[1].replace('\n', '').split(' ') for line in nproclines]

		# filter empties and pick first
		nprocs = [int(list(filter(None, line))[0]) for line in nprocs]

		# Multiply
		for n in nprocs:
			nProcTotal *= n

		return nProcTotal

	def factorize(self, n):
		primfac = []
		d = 2
		while d * d <= n:
			while (n % d) == 0:
				primfac.append(d)  # supposing you want multiple factors repeated
				n //= d
			d += 1
		if n > 1:
			primfac.append(n)
		return primfac

	def product(self, matrix):
		x = 1
		for row in matrix:
			for val in row:
				x*=val
		return x

	def mod_cat_nprocMax(self, catName, catFile, nProcMax):
		assert(nProcMax > 0)

		cat_content = None
		with open(catFile, 'r') as file:
			cat_content = file.readlines()

		nprocLines = [(idx, line) for idx, line in enumerate(cat_content) if line.startswith('nproc')]
		nprocIdx, nprocLines = zip(*nprocLines)

		nprocs = [line.split('=')[1].replace('\n', '').split(' ') for line in nprocLines]

		# filter empties and pick first
		nprocs = [int(list(filter(None, line))[0]) for line in nprocs]

		nProcTotal = 1
		for n in nprocs:
			nProcTotal *= n

		if nProcTotal < nProcMax:
			return catFile

		#
		#	here starts the fun
		#

		# factorize everything we have
		nprocs_fac = [self.factorize(n) for n in nprocs]
		max_fac = self.factorize(nProcMax)
		nprocs_fac_new=[[],[],[]]

		accumulated = 1


		# go through the factors of nproc_max
		for x in max_fac:

			# recompute the iteration order after every iteration
			# we want to go from the largest nproc to the smallest
			for i, facs in enumerate(nprocs_fac):
				nprocs[i] = self.product([facs])

			order = [None]*3
			order[0] = argmax_index(nprocs)
			order[2] = argmin_index(nprocs)
			for i in range(3):
				if i not in order:
					order[1] = i
					break


			associated = False

			for i in order:
				for j, y in enumerate(nprocs_fac[i]):
					# if we have the same factor available in one of our nprocs its easy
					# we just move it then to the new container
					if y == x:
						nprocs_fac_new[i].append(y)
						del nprocs_fac[i][j]
						associated = True
						break
				if associated:
					break

			# if we dont find anything accumulate it
			if not associated:
				accumulated *= x
				for i in order:
					for j, y in enumerate(nprocs_fac[i]):
						# we shall now trade it to the next smaller available factor
						if y < accumulated:
							nprocs_fac_new[i].append(y)
							del nprocs_fac[i][j]
							associated = True
							accumulated = 1
							break
					if associated:
						break

		for i, facs in enumerate(nprocs_fac_new):
			nprocs[i] = self.product([facs])

		nprocMod = self.product([nprocs])
		print('Modded nprocs', nprocs, nprocMod)

		# Modify the cat content
		cat_content[nprocIdx[0]] = "nprocx = " + str(nprocs[0]) + '\n'
		cat_content[nprocIdx[1]] = "nprocy = " + str(nprocs[1]) + '\n'
		cat_content[nprocIdx[2]] = "nprocz = " + str(nprocs[2]) + '\n'

		catNameMod = catName + '.modded'
		catFileMod = catFile[:-4] + '.modded.cat'

		with open(catFileMod, 'w') as file:
			file.write(''.join(cat_content))

		return catNameMod, catFileMod, nprocMod


	def __qsub_serial(self, catName, catFile):
		# Individual qscript
		qscript_name = 'qscript_' + catName + '_serial.sh'
		qscript_path = os.path.join(self.testPath, qscript_name)
		qscript_content = get_qsub_serial(self.testPath, catName, self.compilePath_mpi)

		# write the script
		with open(qscript_path, 'w') as file:
			file.write("".join(qscript_content))

		# submit script
		try:
			print('  submitting serial test %s' % catName)
			p = subprocess.check_call(['qsub', qscript_name], cwd=self.testPath)
		except subprocess.CalledProcessError as e:
			print('  failed!', e.message)

	def qsub_serial(self):
		for catName, catFile in zip(self.catNames, self.catFiles):
			self.__qsub_serial(catName, catFile)


	def __qsub_mpi(self, catName, catFile, nProcMax):
		nproc = self.get_nproc(catFile)

		# reduce amound of processes if necessary
		if nProcMax is not None and nproc > nProcMax:
			catName, catFile, nproc = self.mod_cat_nprocMax(catName, catFile, nProcMax)

		# Individual qscript
		qscript_name = 'qscript_' + catName + '.sh'
		qscript_path = os.path.join(self.testPath, qscript_name)
		qscript_content = get_qsub_mpi(self.testPath, catName, self.compilePath_mpi, nproc)

		# write the script
		with open(qscript_path, 'w') as file:
			file.write("".join(qscript_content))

		# submit script
		try:
			print('  submitting mpi test %s' % catName)
			p = subprocess.check_call(['qsub', qscript_name], cwd=self.testPath)
		except subprocess.CalledProcessError as e:
			print('  failed!', e.message)

	def qsub_mpi(self, nProcMax=None):
		for catName, catFile in zip(self.catNames, self.catFiles):
			self.__qsub_mpi(catName, catFile, nProcMax)

	def qsub_auto(self, nProcMax=None):
		for catName, catFile in zip(self.catNames, self.catFiles):
			nproc = self.get_nproc(catFile)
			if nproc == 1:
				self.__qsub_serial(catName, catFile)
			else:
				self.__qsub_mpi(catName, catFile, nProcMax)
		
	def __run_serial(self, catName, catFile):
		# run test
		try:
			print('  running serial test %s' % catName)
			p = subprocess.check_call(
					[self.compilePath_serial, self.testPath, catName]
				, cwd=self.testPath)
		except subprocess.CalledProcessError as e:
			print('  failed!')
			print(e)

	def run_serial(self):
		for catName, catFile in zip(self.catNames, self.catFiles):
			self.__run_serial(catName, catFile)

	def __run_mpi(self, catName, catFile, nProcMax):
		nproc = self.get_nproc(catFile)

		# reduce amound of processes if necessary
		if nProcMax is not None and nproc > nProcMax:
			catName, catFile, nproc = self.mod_cat_nprocMax(catName, catFile, nProcMax)

		# run test
		try:
			print('  running mpi test %s' % catName)
			p = subprocess.check_call(
				['mpirun', '-n', str(nproc), '-v', self.compilePath_mpi, self.testPath, catName]
				, cwd=self.testPath)
		except subprocess.CalledProcessError as e:
			print('  failed!')
			print(e)

	def run_mpi(self, nProcMax=None):
		for catName, catFile in zip(self.catNames, self.catFiles):
			self.__run_mpi(catName, catFile, nProcMax)

	def run_auto(self, nProcMax=None):
		for catName, catFile in zip(self.catNames, self.catFiles):
			nproc = self.get_nproc(catFile)
			if nproc == 1:
				self.__run_serial(catName, catFile)
			else:
				self.__run_mpi(catName, catFile, nProcMax)


