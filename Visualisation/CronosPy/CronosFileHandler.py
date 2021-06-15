import os
import sys
import re
import h5py
import numpy as np

class CronosFileHandler(object) :
    def set_dataDir(self, DataDir) :
        self.datadir = DataDir
        print " Data dir set to: '" + self.datadir +"'"

    def scan_Projs(self) :
        
        filenames = []
        fullnames = []

        for file in os.listdir(self.datadir):
            if(self.string_has_patt(file, ".*_float\Z")) :
                filenames.append(file)
                fullnames.append(self.datadir + "/" + file)
#                print file


        numDirs = len(filenames)

        self.numProjs = 00

        self.projs = []
        self.nFrames = []

        for iDir in range(numDirs) :
            dirloc = filenames[iDir]
            nameLength = len(dirloc)

            pname = dirloc[0:nameLength-6]
            dirname = self.datadir + "/" + pname +"_float";

            numLoc = self.count_files_serial(dirname, ".h5", True)

            if(numLoc > 0) :
                self.numProjs += 1
                self.projs.append(pname)
                self.nFrames.append(numLoc)

        print " Reduced num of projects: " + str(self.numProjs)

        for iProj in range(self.numProjs) :
            message = " Project: " + self.projs[iProj] + " with " + str(self.nFrames[iProj])
            if(self.nFrames[iProj] > 1) :
                message += " files"
            else :
                message += " file"
            print message

        return self.numProjs


    def scan_tsteps(self, i_proj) :
        # determine available timesteps
	float_dir = self.datadir + "/";
	float_dir += self.projs[i_proj];
	float_dir += "_float/";
        print " Reading files in dir: " + float_dir

	# Now we get all available files in the float-data dir -- files
	# are saved in vector "files"
	self.numFiles = self.count_files_serial(float_dir, ".h5", True);

        # clear lists
	self.files_steps = [];
	self.files_times = [];
	self.filenames = [];
        times_loc = [];
        pos = np.array(range(self.numFiles), dtype='i')

        self.files.sort()

	# Array 'files' holds all filenames:
	for ifile in range(self.numFiles) :
            filename = float_dir + self.files[ifile];

            # Open the file
            h5in = h5py.File(filename, 'r');

            # open data group
            h5group = h5in['Data']

            time = 0.;
            time = h5group.attrs['time']
            times_loc.append(time)

        # Now sort the stuff...
	pos[0] = 0;

        for ifile in range(self.numFiles) :
            pos[ifile] = 0

        for ifile in range(1, self.numFiles) :

            jfile=0;

            while(times_loc[ifile] > times_loc[pos[jfile]] and jfile<ifile) :
                print " other file " + str(jfile) + " " + str(ifile)
                jfile += 1

            print " Found jfile: " + str(jfile);
            for ipos in range(ifile, jfile, -1) :
                pos[ipos] = pos[ipos-1];
                print " Copy " + str(ipos-1) + " " + str(ipos)

            pos[jfile] = ifile;



        # Now store filenames, timesteps and times in correct order
        for ifile in range(self.numFiles) :

            filenum = pos[ifile];
            print " pos " + str(pos[ifile])
                
            # Open file:
            filename = float_dir + self.files[filenum];

            self.filenames.append(self.files[filenum]);

            # Open the file
            h5in = h5py.File(filename, 'r');

            # open data group
            h5group = h5in['Data']
            time = 0.;
            time = h5group.attrs['time']
            
            self.files_times.append(time);
#            print " time: " + str(ifile) + " " + str(time) + " " + filename

            tstep = 0
            tstep = h5group.attrs["timestep"];
            self.files_steps.append(tstep);

        # Return number of timesteps
	return self.numFiles;



    def scan_DataFields(self, i_proj, i_time) :
	#! Generate the list of fields available in hdf5 file
	#!
	#  We begin by obtaining all names of the data fields
	#
	

	# Get the directory
        float_dir = self.datadir + "/";
	float_dir += self.projs[i_proj] + "_float/";

	# Get filename:
        filename = float_dir + self.filenames[i_time];
	print " My filename " + filename

	# Now open the hdf5-file
        h5in = h5py.File(filename,'r')
        # Open the group
        h5group = h5in['Data']

	# Get number of datasets
        self.numDatasets = 0
        self.numDatasets = h5group.attrs['Entries']

	# Get all field names:
	self.fieldNames = [];
	
        for i_field in range(0, self.numDatasets) :
            
            fieldName = self.GetDatasetName(h5group, i_field);
            self.fieldNames.append(fieldName);

            
        # Write field names:
        for i_field in range(0, self.numDatasets) :
            print " Field: " + str(i_field) + " --> " + self.fieldNames[i_field];

	return self.numDatasets;

        
    def get_ArraySizes(self, i_proj, i_time, i_field) :
	#! Determine extent of arrays
	#! Here we load all properties of the arrays - like, e.g., the
	#  gridpoint-size and the spatial extent.

	# Make filename
	filename = self.make_filename(i_proj, i_time);

	# Now open the hdf5-file
        h5in = h5py.File(filename, 'r');

        # open data group
	#h5group = h5in['Data']
	h5group = h5in['fluid00']

	# Pick array in the file:
        dset = h5group[self.fieldNames[i_field]]

	# Get attributes of dataset
	# Get extent of array
        Nx_inv = dset.shape
        xb = dset.attrs["origin"]
        dx = dset.attrs["delta"]

        # Save corrected array extent
        self.Nx = [1]*3
        self.Nx[0] = Nx_inv[2]
        self.Nx[1] = Nx_inv[1]
        self.Nx[2] = Nx_inv[0]

	# Use extent info to build position arrays:
        self.xPos = np.array(range(self.Nx[0]), dtype='f')
        for ix in range(self.Nx[0]) :
            self.xPos[ix] = xb[0] + dx[0]*ix

        self.yPos = np.array(range(self.Nx[1]), dtype='f')
        for iy in range(self.Nx[1]) :
            self.yPos[iy] = xb[1] + dx[1]*iy

        self.zPos = np.array(range(self.Nx[2]), dtype='f')
        for iz in range(self.Nx[2]) :
            self.zPos[iz] = xb[2] + dx[2]*iz
	
        # Now we return a vector just holding the number of gridpoints:
	return self.Nx;


    def get_pname(self, iproj) :
	if(iproj >= 0 and iproj < self.numProjs) :
            return self.projs[iproj];
        else :
            print " no such iproj "
            sys.exit();

            
    def get_nFrames(self, iproj) :
	if(iproj >= 0 and iproj < self.numProjs) :
            print " my num " + str(self.nFrames[iproj])
            return self.nFrames[iproj];
        else :
            print " no such iproj "
            sys.exit();

    def get_numFields(self) :
        return self.numDatasets
	
    def get_FieldName(self, i_field) :
	#! Return name of chosen dataset.
	return self.fieldNames[i_field];

    def get_time(self, i_time) :
	#! Return time saved in hdf5-file
	return self.files_times[i_time]

    def get_step(self, i_time) :
	#! Return timestep saved in hdf5-file
	return self.files_steps[i_time];


    def get_pos(self, i_plane, i_pos) :
	if(i_plane == 0) : # x-y plane:
            return self.zPos[i_pos];
        elif (i_plane == 1) : # x-z plane
            return self.yPos[i_pos];
	else : # y-z plane
            return self.xPos[i_pos];

    def get_PosArray(self, dir) :
        if(dir==0) :
            PosArr = self.xPos;
        elif (dir==1) :
            PosArr = self.yPos
        elif (dir==2) :
            PosArr = self.zPos
        else :
            print " No such direction " + str(dir)
            sys.exit();
        return PosArr



    def Load_2DData(self, i_proj, i_time, i_dataset, i_plotPlane, i_slice) :

	# Open the file
	# Set up filename:
	pname = self.projs[i_proj]; # Get name of project

	# directory of project
	proj_dir = self.datadir + "/" + pname + "_float/";
	
	# filename
        filename = proj_dir + self.filenames[i_time]

	# Open the file
        h5in = h5py.File(filename, 'r');

        # open data group
        #h5group = h5in['Data']
	h5group = h5in['fluid00']	

	# Get name of dataset:
	dset_name = self.fieldNames[i_dataset]

        # Get the cut-data
        data = self.Read2DFrom3DMatrix(h5group, dset_name,
                                       i_plotPlane, i_slice);

#        print " More of data " + " " + str(data.shape) +  " " + str(data[2,3])
           
        return data   


    def Load_1DData(self, i_proj, i_time, i_dataset, i_dir, i_posPerp) :

	# Open the file
	# Set up filename:
	pname = self.projs[i_proj]; # Get name of project

	# directory of project
	proj_dir = self.datadir + "/" + pname + "_float/";
	
	# filename
        filename = proj_dir + self.filenames[i_time]

	# Open the file
        h5in = h5py.File(filename, 'r');

        # open data group
        #h5group = h5in['Data']
	h5group = h5in['fluid00']

	# Get name of dataset:
	dset_name = self.fieldNames[i_dataset]

        # Get the cut-data
        data = self.Read1DFrom3DMatrix(h5group, dset_name,
                                       i_dir, i_posPerp);

#        print " More of data " + " " + str(data.shape) +  " " + str(data[2,3])
           
        return data   



    def Read2DFrom3DMatrix(self, h5group, ArrayName, dir, posPerp) :

	# Open the dataset:
        dataset = h5group[ArrayName]
#        data = dset[0,2:10,1:9:3]

	# Get dimensions of dataset:
	DIM = 3;

        dims_out = dataset.shape
        dimhdf = len(dims_out)
        
	if not (dimhdf == DIM) :
            print " Wrong dimensionality of input data: "
            sys.exit();
	 
        # Beware - dimensions are swapped
        # 2D array holding size
        mx = [1]*2

	if(dir==0) :
            mx[0] = dims_out[1];
            mx[1] = dims_out[2];
	elif (dir==1) : # x,z
            mx[0] = dims_out[0];
            mx[1] = dims_out[2];
        elif (dir==2) : # y,z
            mx[0] = dims_out[0];
            mx[1] = dims_out[1];
        else :
            print " Error no such direction "
            sys.exit();
	
        data = np.zeros((mx[0],mx[1]), dtype=np.float32)
        if(dir==0) :#/ x,y plane

            data = dataset[posPerp,:,:]

	elif (dir==1) : # x,z plane

#            data = np.zeros((mx[0],mx[1]), dtype=np.float32)
            data = dataset[:,posPerp,:]

        else :

 #           data = np.zeros((mx[0],mx[1]), dtype=np.float32)
            data = dataset[:,:,posPerp]


#        print " Size of data " + str(dir) + " " + str(data.shape) + " " + str(mx[0]) + " " + str(mx[1]) + " " + str(data[2,3])

        return data




    def Read1DFrom3DMatrix(self, h5group, ArrayName, dir, posPerp) :

	# Open the dataset:
        dataset = h5group[ArrayName]
#        data = dset[0,2:10,1:9:3]

	# Get dimensions of dataset:
	DIM = 3;

        dims_out = dataset.shape
        dimhdf = len(dims_out)
        
	if not (dimhdf == DIM) :
            print " Wrong dimensionality of input data: "
            sys.exit();
	 
        # Beware - dimensions are swapped
        # 2D array holding size

	if(dir==0) :
            mx = dims_out[2];
	elif (dir==1) : # x,z
            mx = dims_out[1];
        elif (dir==2) : # y,z
            mx = dims_out[0];
        else :
            print " Error no such direction "
            sys.exit();
	
        data = np.zeros((mx), dtype=np.float32)
        if(dir==0) :#/ x,y plane

            data = dataset[posPerp[1],posPerp[0],:]

	elif (dir==1) : # x,z plane

#            data = np.zeros((mx[0],mx[1]), dtype=np.float32)
            data = dataset[posPerp[1],:,posPerp[0]]

        else :

 #           data = np.zeros((mx[0],mx[1]), dtype=np.float32)
            data = dataset[:,posPerp[1],posPerp[0]]


#        print " Size of data " + str(dir) + " " + str(data.shape) + " " + str(mx[0]) + " " + str(mx[1]) + " " + str(data[2,3])

        return data




    def count_files_serial(self, DirName, patt, from_end):
        
        if from_end == True:
            my_string = ".*" + patt + "\Z"
        else:
#            my_string = "\A" + patt + ".*"
            my_string = ".*" + patt + ".*"

        my_regexp = re.compile(my_string)

        coord_regexp = re.compile(".*_coord.*")

        numFiles = 0
        myNum = 0
        self.files = []
        
        # get all files first:
        for file in os.listdir(DirName):
            myNum += 1

            if self.string_has_regexp(file, my_regexp):
                if not self.string_has_regexp(file, coord_regexp):
                    numFiles += 1
                    self.files.append(file)

        return numFiles

                        
    def string_has_regexp(self, string, regexp):

        number = len(regexp.findall(string))
        if number == 1:
            return True
        else:
            return False
    def string_has_patt(self, string, patt):

        regexp = re.compile(patt)
        number = len(regexp.findall(string))
        if number == 1:
            return True
        else:
            return False


    def GetDatasetName(self, h5group, num) :

        AttrName = "Name_om" + "%2.2d" % num
	
	DatasetName = h5group.attrs[AttrName]

	return DatasetName

    def make_filename(self, i_proj, i_time) :
	#! Make the filename

	# Make the directory
	float_dir = self.datadir;
	float_dir += "/";
	float_dir += self.projs[i_proj];
	float_dir += "_float/";

	# Make filename:
	filename = float_dir + self.filenames[i_time];
	return filename;
