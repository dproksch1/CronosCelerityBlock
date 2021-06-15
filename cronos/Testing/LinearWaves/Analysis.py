import math
import os
import re
import h5py
import numpy as np
import matplotlib.pyplot as plt

class WaveAnalysis(object) :
    def __init__(self) :
        homeDir = os.environ['HOME']
        self.fileDir = homeDir + "/data/Cronos/StandardTests/ShearAlfven/"
        #self.waveName = "DecayingShearAlfWave_64"
        #self.waveName = "DecayingFastWave_128"
        #self.waveName = "DecayingSlowWave_128"
        #        self.waveName = "DecayingSlowWave_128_hll_minMod"
        #        self.waveName = "DecayingShearAlfWave_16_hll_vL"
        #self.waveName = "DecayingShearAlfWave_64_hll_vL"
        #self.waveName = "DecayingShearAlfWave_64_hll_long"
        #self.waveName = "DecayingShearAlfWave_64_hll_100"
        #self.waveName = "DecayingShearAlfWave_64_d1440_100"
        #self.waveName = "DecayingSoundWave_64"
        #self.waveName = "DecayingShearFlow_32_Dens"
        #self.waveName = "DecayingShearFlow_8_L25_cs150"
        self.waveName = "DecayingShearAlfWaveHll_64"
        #self.waveName = "DecayingShearFlow_32_L50_cs150"
        #self.waveName = "DecayingFastWave_64_hll_var"
        #self.waveName = "DecayingFastWave_128_hll_vL"
        #self.waveName = "DecayingSlowWave_8_hll_vL"
        #self.waveName = "DecayingSlowWave_128"
        #self.waveName = "DecayingShearAlfWave_16_hlld"
        #self.waveName = "DecayingFastWave_128_hlld"
        #self.waveName = "DecayingSlowWave_8_hlld"
        self.waveDir = ""

        self.vz_rms = []
        self.times = []
        self.Gamma = 1.
        self.v0_rms = 0.

        self.cAlfven = 0.
        self.cWave_Alfven = 0.
        self.cWave_Fast = 0.
        self.cWave_Slow = 0.
        self.vamp = 0.
        self.wave_type = 0

        self.version = 1

    def setup(self, waveName) :
        self.waveName = waveName

    def listFiles(self) :
        dirName = self.fileDir + self.waveName + "_float/"
#        print dirName
        self.waveDir = dirName

        self.filenames = []
        for file in os.listdir(dirName):
            if(self.string_has_patt(file, ".h5\Z")) :
               self.filenames.append(file)
#            print file

    def run(self) :
        numFiles = len(self.filenames)
        print str(numFiles) + " to process "

        self.times = []
        self.vz_rms = []

        # get properties:
        theName = self.waveDir + "/" + self.filenames[0]
        self.getPropertiesFromFile(theName)

        # get reference data
        self.getRefDataFromFile(theName)

        if(self.version==0) :
            for iFile in range(numFiles) :
                theName = self.waveDir + "/" + self.filenames[iFile]
                self.getDataFromFile(theName)
        else :
            # get last file
            numFiles = len(self.filenames)
            # New version, where output files store full data array
            fname = self.filenames[numFiles-1]
            print " Reading " + fname
            self.times = self.Load_1DData(fname, "times");
            self.vz_rms = self.Load_1DData(fname, "rmsVelo");
            self.v0_rms = self.vz_rms[0]


        if(self.wave_type<4) :
            self.set_Gamma()
        else :
            self.set_GammaDirect()
        


        #       print self.times
        print " Decay rate: " + str(self.Gamma)
        print " Decay time: " + str(1./self.Gamma)
        if(self.wave_type==1) :
            cA = self.cWave_Alfven
            R_A = 8.*math.pi**2*cA/(self.Gamma*1.)
            print " R_A: " + str(R_A)
            print "Gamma: {:.3g}, tau: {:.4g}, RA: {:.2f}".format(self.Gamma, 1./self.Gamma, R_A)
#            print("Gamma: %.3f tau: %.3f" % self.Gamma % 1./self.Gamma)
        elif(self.wave_type==2) :
            cF = self.cWave_Fast
            R_F = 8.*math.pi**2*cF/(self.Gamma*1.)
            print " R_F: " + str(R_F)
            print "Gamma: {:.3g}, tau: {:.4g}, RF: {:.2f}".format(self.Gamma, 1./self.Gamma, R_F)
        elif(self.wave_type==3) :
            cS = self.cWave_Slow
            R_S = 8.*math.pi**2*cS/(self.Gamma*1.)
            print " R_S: " + str(R_S)
            print "Gamma: {:.3g}, tau: {:.4g}, RS: {:.2f}".format(self.Gamma, 1./self.Gamma, R_S)
        elif(self.wave_type==4) :
            cS = self.cSound
            print "Gamma: {:.3g}, tau: {:.4g}".format(self.Gamma, 1./self.Gamma)
        elif(self.wave_type==5) :
            cS = self.cSound
            #rho0 = rho0
            eta = self.Gamma*self.rho0*self.Len**2/(8.*math.pi**2)
            RShear = cS*self.Len*self.rho0/(2.*eta)
            print " eta: " + str(eta)
            print " Gamma: " + str(self.Gamma)
            print " Len: " + str(self.Len)
            print " Reyn: " + str(RShear)
            print "Gamma: {:.3g}, tau: {:.4g}, Rshear: {:.2f}".format(self.Gamma, 1./self.Gamma,RShear)
        


        # make plot of decay rates
        times_gam = []
        val_gam = []
        i_t_num = 256
        numTimes = len(self.times)
        for i_t in range(i_t_num) :
            t_val = i_t/(i_t_num-1.)*self.times[numTimes-1]
            g_val = self.v0_rms*math.exp(-t_val*self.Gamma)

            times_gam.append(t_val)
            val_gam.append(g_val)


        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        ax.plot(self.times, self.vz_rms,
                color='k',linestyle='-');
        ax.plot(times_gam, val_gam,
                color='r',linestyle=':');
        ax.set_yscale('log')
        plt.ylim([self.AMax/20.,2.*self.AMax])

        #        print times_gam
        #        print val_gam

        if(self.wave_type<4) :
            plt.ylim(0.01,1.1)
        else :
            maxVal = max(plt.ylim(self.vz_rms[numTimes-1],self.vz_rms[0]))
            minVal = min(plt.ylim(self.vz_rms[numTimes-1],self.vz_rms[0]))
            delta = maxVal/minVal
            print " range: " + str(minVal) + " " + str(maxVal)
            plt.ylim(minVal/(delta),maxVal*delta)
        plt.xlim(0,self.times[numTimes-1]*1.05)
        fig.savefig("wave.pdf")

    def set_Gamma(self) :
        # find time of maxium:
        numTimes = len(self.times)

        oldVal = self.vz_rms[numTimes-2]
        vOldVal = self.vz_rms[numTimes-1]
        for i_t in range(numTimes-3,0,-1) :
            print " time: " + str(i_t) + " " + str(self.times[i_t]) + " " + str(oldVal) + " " + str(self.vz_rms[i_t])
            currVal = self.vz_rms[i_t]
            if(oldVal > currVal and oldVal > vOldVal) :
                maxTime = self.times[i_t+1]
                maxVal = self.vz_rms[i_t+1]
                break
            else :
                vOldVal = oldVal
                oldVal = currVal

        print " mehr " + str(self.v0_rms)

        print " Values: " + str(maxTime) + " " + str(maxVal)
        self.Gamma = 1./maxTime*math.log(self.v0_rms/maxVal)

 
    def set_GammaDirect(self) :
        # find time of maxium:
        numTimes = len(self.times)
        maxTime = self.times[numTimes-1]
        minTime = self.times[0]

        maxVal = self.vz_rms[numTimes-1]
        minVal = self.vz_rms[0]
        self.Gamma = 1./(maxTime-minTime)*math.log(minVal/maxVal)

        for i_t in range(numTimes) :
            print " time: " + str(self.times[i_t]) + " " + str(self.vz_rms[i_t])
        

    def getPropertiesFromFile(self, fname) :
        h5in = h5py.File(fname, 'r');
        # open data group
        h5group = h5in['Data/WavePars']

        self.wave_type = h5group.attrs["Wave_Type"]
        self.vamp = h5group.attrs["Initial_Amplitude"]
        self.cAlfven = h5group.attrs["cA"]
        self.cWave_Alfven = h5group.attrs["cWave_Alfven"]
        self.cWave_Fast = h5group.attrs["cWave_Fast"]
        self.cWave_Slow = h5group.attrs["cWave_Slow"]
        self.cSound = h5group.attrs["cs"]
        self.Len = h5group.attrs["Len"]
        if(self.wave_type==5) :
            self.rho0 = h5group.attrs["BackgroundDensity"]
        
        if(self.wave_type==1) :
            self.AMax = self.vamp*self.cWave_Alfven
        elif(self.wave_type==2) :
            self.AMax = self.vamp*self.cWave_Fast
            self.AMax *= 2.*self.cWave_Fast**2 - self.cAlfven**2
            self.AMax /= self.cWave_Fast**2 - self.cAlfven**2
        elif(self.wave_type==3) :
            self.AMax = self.vamp*self.cWave_Slow
            self.AMax *= 2.*self.cWave_Slow**4 - self.cWave_Slow**2*self.cAlfven**2 + self.cAlfven**4
            self.AMax /= (self.cWave_Slow**2 - self.cAlfven**2)**2
#            self.AMax *= 2.*self.cWave_Slow**2 - self.cAlfven**2
#            self.AMax /= self.cWave_Slow**2 - self.cAlfven**2
        elif(self.wave_type==4 or self.wave_type==5) :
            self.AMax = self.vamp*self.cSound


    def getRefDataFromFile(self, fname) :
        # Open the file
#        print fname
        h5in = h5py.File(fname, 'r');
        # open data group
        h5group = h5in['Data']
        
        time = 0.;
        time = h5group.attrs['time']


        # Load velocity
        h5group = h5in['fluid00']
        #        print " v_z " + str(v_z[12,3]) + " " + str(v_z[14,3])

        self.ref_vx = self.Load_2DData(h5group, "v_x", 0, 0)
        self.ref_vy = self.Load_2DData(h5group, "v_y", 0, 0)
        self.ref_vz = self.Load_2DData(h5group, "v_z", 0, 0)



    def getDataFromFile(self, fname) :
        # Open the file
#        print fname
        h5in = h5py.File(fname, 'r');
        # open data group
        h5group = h5in['Data']
        
        time = 0.;
        time = h5group.attrs['time']
        self.times.append(time)


        # Load velocity
        h5group = h5in['fluid00']
        #        print " v_z " + str(v_z[12,3]) + " " + str(v_z[14,3])

        if(self.wave_type==1) :
            v_z = self.Load_2DData(h5group, "v_z", 0, 0)
            rms_velo = self.rms(v_z)
        else :
            v_x = self.Load_2DData(h5group, "v_x", 0, 0)
            v_y = self.Load_2DData(h5group, "v_y", 0, 0)
            rms_velo = self.rms(v_x) + self.rms(v_y)
            rms_velo = self.rms(v_y)
            rms_velo = self.rms2(v_x, v_y)
            rms_velo = self.rms(v_x/self.ref_vx)
            #            rms_velo = self.get_amplitude(v_x, v_y)
            rms_velo = self.rms2(v_x, v_y)


        print " rms: " + str(rms_velo) + " " + str(time)
        self.vz_rms.append(rms_velo)

        if(time < 1.e-8) :
            self.v0_rms = rms_velo

        if(time > 10-1.e-6 and time < 10.+1.e-6) :
            # Compute decay rate
#            print " gamma at time " + str(time)
            self.Gamma = 0.1*math.log(self.v0_rms/rms_velo)

    def get_amplitude(self, field1, field2) :
        if(self.wave_type==3) :
            famp = field1/(self.ref_vx + 1.e-15)
            amp = math.sqrt(self.rms(famp)**2)
            print amp
        return amp
            


    def rms(self, field) :
        value = np.mean(np.square(field))
        value = math.sqrt(value)
        return value

    def rms2(self, field1, field2) :
        value1 = np.mean(np.square(field1))
        value2 = np.mean(np.square(field2))
        value = value1 + value2
        value = math.sqrt(value)
        if(self.wave_type==2) :
            value = self.cWave_Fast**2*value1 
            value += (self.cWave_Fast**2/(self.cWave_Fast**2 - self.cAlfven**2))**2*value2
            value = math.sqrt(value)
        if(self.wave_type==3) :
            value = self.cWave_Slow**2*value1 
            value += (self.cWave_Slow**2/(self.cWave_Slow**2 - self.cAlfven**2))**2*value2
            value = math.sqrt(value)
        return value


    def Load_1DData(self, fname, dset_name) :
        fname = self.waveDir + "/" + fname
        #print " Opening " + fname
        h5in = h5py.File(fname, 'r');
        h5group = h5in['Data/WavePars']
        #        h5group = h5in['Data']
        print " dataset: " + str(dset_name)
        #print h5group.items()
        dataset = h5group[dset_name]
        data = dataset[:]
        #data = [1,1]
        return data

    def Load_2DData(self, group, dset_name, dir, posPerp) :
        dataset = group[dset_name]
        DIM = 3;

        dims_out = dataset.shape
        dimhdf = len(dims_out)

        if not (dimhdf == DIM) :
            print " Wrong dimensionality of input data: "
            sys.exit();

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


        return data

#        print " Size of data " + str(dir) + " " + str(data.shape) + " " + str(mx[0]) + " " + str(mx[1]) + " " + str(data[2,3])



    def string_has_patt(self, string, patt):

        regexp = re.compile(patt)
        number = len(regexp.findall(string))
        if number == 1:
            return True
        else:
            return False

if __name__ == "__main__":
    print " What?? "
    Analysis = WaveAnalysis()
    Analysis.listFiles()
    Analysis.run()
