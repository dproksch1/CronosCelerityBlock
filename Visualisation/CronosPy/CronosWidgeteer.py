from PyQt4 import QtCore, QtGui
import CronosWidgetDesign
import CronosFileHandler
import re
import math
import sys
import PlotDataStorage
import numpy as np
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
matplotlib.rcParams['ps.useafm'] = 'True'


class CronosWidgeteer(CronosWidgetDesign.Ui_Form, QtGui.QWidget):
    def __init__(self, parent=None):
        CronosWidgetDesign.Ui_Form.__init__(self,parent)
        QtGui.QWidget.__init__(self,parent)

        # Display widget - inherited from CronosWidgetDesign
        self.setupUi(self)

    def setupWidget(self, Widget):

        # Make file handler object
        self.ProjHandler = CronosFileHandler.CronosFileHandler()
        self.ProjHandler.set_dataDir("/home/dave/work/Cronos/data")

        self.connect(self.butt_Exit, QtCore.SIGNAL('clicked()'),
                     QtGui.qApp, QtCore.SLOT('quit()'))
        self.connect(self.butt_ScanFolder, QtCore.SIGNAL('clicked()'),
                     self.scanProjs)

        # Scan for timesteps and fields
        self.connect(self.butt_ScanFrames, QtCore.SIGNAL('clicked()'),
                     self.scanTSteps);


        # Get chosen project
	self.connect(self.list_projs, QtCore.SIGNAL('itemSelectionChanged()'),
                     self.set_iproj)

        # Selection of field
        self.connect(self.select_fields, QtCore.SIGNAL('activated(int)'),
                     self.set_iQuant);
        # Selection of mode
        self.connect(self.select_mod, QtCore.SIGNAL('activated(int)'),
                     self.set_iMod);

        # get chosen timestep
        self.spinBox_timestep.setRange(0, 0);
	self.slider_timestep.setRange(0, 0);
        self.connect(self.spinBox_timestep, QtCore.SIGNAL('valueChanged(int)'),
                     self.slider_timestep, QtCore.SLOT('setValue(int)'));
	self.connect(self.slider_timestep, QtCore.SIGNAL('valueChanged(int)'),
                     self.spinBox_timestep, QtCore.SLOT('setValue(int)'));
        self.connect(self.slider_timestep, QtCore.SIGNAL('valueChanged(int)'),
                     self.set_iTime);

        # Choice of plane
        self.connect(self.butt_plane_xy, QtCore.SIGNAL('toggled(bool)'),
                     self.set_xyPlane);
	self.connect(self.butt_plane_xz, QtCore.SIGNAL('toggled(bool)'),
                     self.set_xzPlane);
	self.connect(self.butt_plane_yz, QtCore.SIGNAL('toggled(bool)'),
                     self.set_yzPlane);
        self.butt_plane_xy.click()


        # Choice of slice
        self.spinBox_slice.setRange(0, 0);
	self.slider_slice.setRange(0, 0);
        self.connect(self.spinBox_slice, QtCore.SIGNAL('valueChanged(int)'),
                     self.slider_slice, QtCore.SLOT('setValue(int)'));
	self.connect(self.slider_slice, QtCore.SIGNAL('valueChanged(int)'),
                     self.spinBox_slice, QtCore.SLOT('setValue(int)'));
        self.connect(self.slider_slice, QtCore.SIGNAL('valueChanged(int)'),
                     self.set_iSlice);
	self.spinBox_slice.setEnabled(False);
	self.slider_slice.setEnabled(False);


        # Check if log is used
	self.connect(self.butt_linear, QtCore.SIGNAL('toggled(bool)'),
                     self.set_linear);
	self.butt_linear.click();


        # Selection of output format:
	self.connect(self.select_OutputFormat, QtCore.SIGNAL('activated(int)'),
                     self.set_OutputFormat);
	self.set_FormatOptions();
	self.set_OutputFormat(0); # x-window as default


         # Selection of contour format
	self.connect(self.butt_showCells, QtCore.SIGNAL('toggled(bool)'),
                     self.set_contourFormat);
	self.butt_showCells.click();


        # Selection of geometry type
        self.connect(self.butt_CoordCart, QtCore.SIGNAL('toggled(bool)'),
                     self.set_setGeomCartesian);
        self.connect(self.butt_CoordCyl, QtCore.SIGNAL('toggled(bool)'),
                     self.set_setGeomCylindrical);
        self.connect(self.butt_CoordSph, QtCore.SIGNAL('toggled(bool)'),
                     self.set_setGeomSpherical);
        self.butt_CoordCart.click();

        # Selection of color table:
        self.spinBox_colorTable.setRange(0, 70);
	self.connect(self.spinBox_colorTable,
                     QtCore.SIGNAL('valueChanged(int)'), self.set_colorTable)
        self.spinBox_colorTable.setValue(68)


        # Do a countour plot:
	self.connect(self.butt_plotContour, QtCore.SIGNAL('clicked()'),
                     self.do_plotContour);

        # Initial scan:
        self.scanProjs();

        # default values
        self.i_slice = 0
        self.i_plotPlane = 0
        self.has_figure = False

        # switch off some options
        self.butt_plane_xy.setEnabled(False);
	self.butt_plane_xz.setEnabled(False);
	self.butt_plane_yz.setEnabled(False);


        # Buttons on line plot tab
        # Choice of direction
        self.connect(self.butt_Dir1D_x, QtCore.SIGNAL('toggled(bool)'),
                     self.set_Dir1D_x);
	self.connect(self.butt_Dir1D_y, QtCore.SIGNAL('toggled(bool)'),
                     self.set_Dir1D_y);
	self.connect(self.butt_Dir1D_z, QtCore.SIGNAL('toggled(bool)'),
                     self.set_Dir1D_z);
        self.butt_Dir1D_x.click()

        self.butt_Dir1D_x.setEnabled(False);
        self.butt_Dir1D_y.setEnabled(False);
        self.butt_Dir1D_z.setEnabled(False);

        # Choice of perpendicular position
        self.spinBox_slice1DPerp1.setRange(0, 0);
	self.slider_slice1DPerp1.setRange(0, 0);
        self.connect(self.spinBox_slice1DPerp1,
                     QtCore.SIGNAL('valueChanged(int)'),
                     self.slider_slice1DPerp1, QtCore.SLOT('setValue(int)'));
	self.connect(self.slider_slice1DPerp1,
                     QtCore.SIGNAL('valueChanged(int)'),
                     self.spinBox_slice1DPerp1, QtCore.SLOT('setValue(int)'));
        self.connect(self.slider_slice1DPerp1,
                     QtCore.SIGNAL('valueChanged(int)'),
                     self.set_perpPos1);
	self.spinBox_slice1DPerp1.setEnabled(False);
	self.slider_slice1DPerp1.setEnabled(False);

        self.spinBox_slice1DPerp2.setRange(0, 0);
	self.slider_slice1DPerp2.setRange(0, 0);
        self.connect(self.spinBox_slice1DPerp2,
                     QtCore.SIGNAL('valueChanged(int)'),
                     self.slider_slice1DPerp2, QtCore.SLOT('setValue(int)'));
	self.connect(self.slider_slice1DPerp2,
                     QtCore.SIGNAL('valueChanged(int)'),
                     self.spinBox_slice1DPerp2, QtCore.SLOT('setValue(int)'));
        self.connect(self.slider_slice1DPerp2,
                     QtCore.SIGNAL('valueChanged(int)'),
                     self.set_perpPos2);
	self.spinBox_slice1DPerp2.setEnabled(False);
	self.slider_slice1DPerp2.setEnabled(False);

        self.i_posPerp = [0,0]

        # Do a line plot:
	self.connect(self.butt_LinePlot, QtCore.SIGNAL('clicked()'),
                     self.do_plotLine);


    def scanProjs(self) :
        self.WriteInfo(" Scanning for projects...")
        self.n_proj = self.ProjHandler.scan_Projs();
        
        MaxFrames = 0
        for iproj in range(self.n_proj) :
            MaxFrames = max(MaxFrames, self.ProjHandler.get_nFrames(iproj));

        print " Max number of frames " + str(MaxFrames)

        proj_list = []
        for iproj in range(self.n_proj) :
            list_entry = str(self.ProjHandler.get_nFrames(iproj))
            list_entry += "   "
            list_entry += self.ProjHandler.get_pname(iproj)
            proj_list.append(list_entry)

        # Add projects to the list
        for iproj in range(self.n_proj) :
            self.list_projs.addItem(proj_list[iproj])

        self.butt_ScanFrames.setEnabled(False);
	self.WriteInfo("...done (scan projs)");
        
        
    def scanTSteps(self) :
	self.WriteInfo("Scanning for timesteps...");
	# Get chosen project:
	pname = self.ProjHandler.get_pname(self.i_proj);
	print " Using pname " + pname

	print " getting steps "

	self.n_times = self.ProjHandler.scan_tsteps(self.i_proj);

	print " Got steps "

	# Fixing max number of timesteps for slider:
	self.spinBox_timestep.setRange(0, self.n_times-1);
	self.slider_timestep.setRange(0, self.n_times-1);

	# Set step=0 as default:
	self.i_time = 0;
	self.set_iTime(0);

	# Now identify all available files (use timestep 0 as standard)
	self.numFields = self.ProjHandler.scan_DataFields(self.i_proj, 0);

	self.WriteInfo("...done (Scanning timesteps)");

        

    def set_iproj(self) :
	# Get selection of project when clicked
	self.i_proj = self.list_projs.currentRow();
        message = " Project selected: " + str(self.i_proj)
	self.WriteInfo(message);
	message = " chosen project: ";
	message += self.ProjHandler.get_pname(self.i_proj);
	self.butt_ScanFrames.setEnabled(True);
        # Do initial scan:
        self.scanTSteps()
	self.WriteInfo(message);


    def set_iTime(self, i_Time) :
        self.WriteInfo("Setting timestep to" + str(i_Time));
        self.i_time = i_Time

        # Obtain corresponding time:
        time = self.ProjHandler.get_time(self.i_time);
	print " time is: " + str(time);

        # Now write a string holding the time:
	message = " Time: " + str(time)
        message += "  at step: " + str(self.ProjHandler.get_step(self.i_time))

        print message
        self.info_timestep.setText(message)

        # Now identify all available files (use timestep 0 as standard)
        numFields = self.ProjHandler.scan_DataFields(self.i_proj, self.i_time);

        # Make list of all sensible fields:
	self.make_QuantList();
        
        # Set default quantity:
        self.set_iQuant(0)

    def set_iQuant(self, i_Quant) :
	#! Select quantity to plot (like, i.e., velocity
	self.i_quant = i_Quant;
	print " My Quantity: " + str(self.i_quant);
	print " corresponds to " + self.AllFields[self.i_quant]

	# Now determine possible sub-options:
	self.make_ModList();
	self.set_iMod(0);


    def set_iMod(self, i_Mod) :
	self.i_mod = i_Mod;
	print " my mod: " + str(self.i_mod) + " " + self.Mods[self.i_mod]

	# Compute arr-index
	# Start with vectorial quantities:
	
	if(self.string_has_patt(self.AllFields[self.i_quant],"Velocity")) :
            if(self.string_has_patt(self.Mods[self.i_mod],"v_y")) :
                self.i_dataset = self.NameIDs["v_y"];
            elif (self.string_has_patt(self.Mods[self.i_mod],"v_z")) :
                self.i_dataset = self.NameIDs["v_z"];
            else :
                self.i_dataset = self.NameIDs["v_x"];
        elif(self.string_has_patt(self.AllFields[self.i_quant],"Mag. field")) :
            if(self.string_has_patt(self.Mods[self.i_mod],"B_y")) :
                self.i_dataset = self.NameIDs["B_y"];
            elif (self.string_has_patt(self.Mods[self.i_mod],"B_z")) :
                self.i_dataset = self.NameIDs["B_z"];
            else :
                self.i_dataset = self.NameIDs["B_x"];
        else :
            my_name = self.AllFields[self.i_quant];
            self.i_dataset = self.NameIDs[my_name];
	

	# Test if id is correctly identified
	print " You chose (at least): " + self.ProjHandler.get_FieldName(self.i_dataset);

	self.Nx = self.ProjHandler.get_ArraySizes(self.i_proj, self.i_time,
                                                  self.i_dataset);

	# Enable plot choices
	self.butt_plane_xy.setEnabled(True);
	self.butt_plane_xz.setEnabled(True);
	self.butt_plane_yz.setEnabled(True);
	self.spinBox_slice.setEnabled(True);
	self.slider_slice.setEnabled(True);

        self.butt_Dir1D_x.setEnabled(True);
        self.butt_Dir1D_y.setEnabled(True);
        self.butt_Dir1D_z.setEnabled(True);
	self.spinBox_slice1DPerp1.setEnabled(True);
	self.slider_slice1DPerp1.setEnabled(True);
	self.spinBox_slice1DPerp2.setEnabled(True);
	self.slider_slice1DPerp2.setEnabled(True);

	# Set corresponding default values:
	self.set_iSlice(self.i_slice);
	self.spinBox_slice.setValue(self.i_slice);

	print " The number of gridpoints: " + str(self.Nx[0]) + " " + str(self.Nx[1]) + " " + str(self.Nx[2]);

	# Write data info into list widget:
	self.show_DataProperties();

	# Set standard extent for slider (standard is x-y plane):
	if(self.butt_plane_xy.isChecked()) :
            self.spinBox_slice.setRange(0,self.Nx[2]-1);
            self.slider_slice.setRange(0,self.Nx[2]-1);
        elif (self.butt_plane_xz.isChecked()) :
            self.spinBox_slice.setRange(0,self.Nx[1]-1);
            self.slider_slice.setRange(0,self.Nx[1]-1);
        else :
            self.spinBox_slice.setRange(0,self.Nx[0]-1);
            self.slider_slice.setRange(0,self.Nx[0]-1);
            
	
    def set_iSlice(self, i_Slice) :
	print " Setting slice to " + str(i_Slice);
	self.i_slice = i_Slice;

	# Obtain corresponding position:
        pos = self.ProjHandler.get_pos(self.i_plotPlane, self.i_slice);

	# Now write a string holding the position
        message = " Position: " + str(pos)
	self.info_slice.setText(message);

    def set_xyPlane(self, checked) :
	#! Do plot for x-y plane
	self.i_plotPlane = 0;
	if(self.butt_plane_xy.isEnabled()) :
            self.spinBox_slice.setRange(0,self.Nx[2]-1);
            self.slider_slice.setRange(0,self.Nx[2]-1);
            # Initial value for slice
            self.i_slice = self.Nx[2]/2;
            self.set_iSlice(self.i_slice);
            self.spinBox_slice.setValue(self.i_slice);

    def set_xzPlane(self, checked) :
	#! Do plot for x-z plane
	self.i_plotPlane = 1;
	if(self.butt_plane_xz.isEnabled()) :
            self.spinBox_slice.setRange(0,self.Nx[1]-1);
            self.slider_slice.setRange(0,self.Nx[1]-1);
            # Initial value for slice
            self.i_slice = self.Nx[1]/2;
            self.set_iSlice(self.i_slice);
            self.spinBox_slice.setValue(self.i_slice);


    def set_yzPlane(self, checked) :
	#! Do plot for y-z plane
	self.i_plotPlane = 2;
	if(self.butt_plane_yz.isEnabled()) :
            self.spinBox_slice.setRange(0,self.Nx[0]-1);
            self.slider_slice.setRange(0,self.Nx[0]-1);
            # Initial value for slice
            self.i_slice = self.Nx[0]/2;
            self.set_iSlice(self.i_slice);
            self.spinBox_slice.setValue(self.i_slice);

    def set_contourFormat(self, checked) :
        # Distinguish two types of contour plotting:
        self.cont_asCells = checked;



    def set_setGeomCartesian(self, checked) :
        self.geometry = 0

    def set_setGeomCylindrical(self, checked) :
        self.geometry = 1

    def set_setGeomSpherical(self, checked) :
        self.geometry = 2

    def make_QuantList(self) :
	#! Make a list of plotable quantities
	#!  From the list of fields obtained by the data-handler we
	#  extract a new list of modified quantities
	#

	# Get number of fields
#        numFields = self.ProjHandler.get_DataFields(self.i_proj, self.i_time);
        numFields = self.ProjHandler.get_numFields();
        
        # Now obtain a list of all available fields:
	name_fields = [];
	is_used = []
        
        for i_field in range(0, numFields) :
            name_fields.append(self.ProjHandler.get_FieldName(i_field))
            is_used.append(False);
	
        # Set all indicators to false:
        self.has_Density = False;
	self.has_Velocity = False;
	self.has_MagField = False;
	self.has_Etherm = False;

        self.AllFields = []
        self.NameIDs = {}


        for i_field in range(0, numFields) :

#            if(self.string_has_patt(name_fields[i_field], ".*rho.*")) :
            if(name_fields[i_field] == "rho") :
                print " Found density: " + name_fields[i_field]
                if not(self.has_Density) :
                    self.AllFields.append("Density")
                    
                self.has_Density = True
                self.NameIDs["Density"] = i_field

            elif(self.string_has_patt(name_fields[i_field], "^v_\w$")) :
#                elif(self.string_has_patt(name_fields[i_field], ".*v_.*")) :
                print " Found a velocity component: " + name_fields[i_field]
                if not(self.has_Velocity) :
                    self.AllFields.append("Velocity")
                    
                self.has_Velocity = True
                if(self.string_has_patt(name_fields[i_field], ".*v_x.*")) :
                    self.NameIDs["v_x"] = i_field
                elif(self.string_has_patt(name_fields[i_field], ".*v_y.*")) :
                    self.NameIDs["v_y"] = i_field
                elif(self.string_has_patt(name_fields[i_field], ".*v_z.*")) :
                    self.NameIDs["v_z"] = i_field

            elif(self.string_has_patt(name_fields[i_field], ".*B_.*")) :
                print " Found a mag component: " + name_fields[i_field]
                if not(self.has_MagField) :
                    self.AllFields.append("Mag. field")
                    
                self.has_MagField = True
                if(self.string_has_patt(name_fields[i_field], ".*B_x.*")) :
                    self.NameIDs["B_x"] = i_field
                elif(self.string_has_patt(name_fields[i_field], ".*B_y.*")) :
                    self.NameIDs["B_y"] = i_field
                elif(self.string_has_patt(name_fields[i_field], ".*B_z.*")) :
                    self.NameIDs["B_z"] = i_field

            elif(self.string_has_patt(name_fields[i_field], ".*Etherm.*")) :
                print " Found energy: " + name_fields[i_field]
                if not(self.has_Etherm) :
                    self.AllFields.append("Etherm")
                    
                self.has_Etherm = True
                self.NameIDs["Etherm"] = i_field

            else :
                self.AllFields.append(name_fields[i_field])
                self.NameIDs[name_fields[i_field]] = i_field

        # Write full list:
        list_entries = len(self.AllFields)
        for i_field in range(list_entries) :
            print self.AllFields[i_field]

        # Now use the list of fields to generate choosable entries
	self.select_fields.clear();

        for i_sel in range(list_entries) :
            field = self.AllFields[i_sel];
            self.select_fields.insertItem(i_sel, field);
	

    def make_ModList(self) :
	#! Generate list of possible modifiers
	#! Possible modifiers are e.g. abs values
	 

	self.Mods = []

	# Walk through options for vector fields -- the rest is standard:
	if(self.string_has_patt(self.AllFields[self.i_quant],".*Velocity.*")) :
            # Test if all components are available
            allComp = False;
            if(self.NameIDs.has_key("v_x") and self.NameIDs.has_key("v_y") and
               self.NameIDs.has_key("v_z")) :
                allComp = True;
		
		
            if(allComp) :
                self.Mods.append(" Abs ");

            if(self.NameIDs.has_key("v_x")) :
                self.Mods.append(" v_x ");

            if(self.NameIDs.has_key("v_y")) :
                self.Mods.append(" v_y ");
		
            if(self.NameIDs.has_key("v_z")) :
                self.Mods.append(" v_z ");
		
            if(allComp) :
                self.Mods.append(" Div ");
                self.Mods.append(" Curl ");
		

        elif (self.string_has_patt(self.AllFields[self.i_quant],".*Mag. field.*")) : # Mag. field
            # Test if all components are available
            allComp = False;
            if(self.NameIDs.has_key("B_x") and self.NameIDs.has_key("B_y") and
               self.NameIDs.has_key("B_z")) :
                allComp = True;
		
            if(allComp) :
                self.Mods.append(" Abs ");
		
            if(self.NameIDs.has_key("B_x")) :
                self.Mods.append(" B_x ");
		
            if(self.NameIDs.has_key("B_y")) :
                self.Mods.append(" B_y ");
		
            if(self.NameIDs.has_key("B_z")) :
                self.Mods.append(" B_z ");
		
            if(allComp) :
                self.Mods.append(" Div ");
                self.Mods.append(" Curl ");
		
        else : # others
            self.Mods.append(" None ");
            self.Mods.append(" Abs ");
	
	
	# Now use the list of options to generate choosable entries
	self.select_mod.clear();

	for i_sel in range(len(self.Mods)) :
            field = self.Mods[i_sel];
            self.select_mod.insertItem(i_sel, field);
	

    def set_linear(self, checked) :
        self.lin_plot = checked;
        print " Doing a linear plot: " + str(self.lin_plot)


    def set_FormatOptions(self) :
	#! Set available output format options
	self.select_OutputFormat.insertItem(0, "x-window");
	self.select_OutputFormat.insertItem(1, "Postscript");
	self.select_OutputFormat.insertItem(2, "PNG");
	self.select_OutputFormat.insertItem(3, "PDF");
	self.select_OutputFormat.insertItem(4, "SVG");

    def set_OutputFormat(self, outOption) :
	 
        if(outOption==0) :
            self.fileFormatContour = "display";
        elif(outOption==1) :
            self.fileFormatContour = "ps";
        elif(outOption==2) :
            self.fileFormatContour = "png";
        elif(outOption==3) :
            self.fileFormatContour = "pdf";
        elif(outOption==4) :
            self.fileFormatContour = "svg";
        self.WriteInfo(" Output format is: " + self.fileFormatContour)



    def set_colorTable(self, colorOption) :
	self.i_color = colorOption
#        print " My color " + self.select_colorTable.currentText()

        if(self.i_color==0) :
            self.name_colorTable = "afmhot"
        elif(self.i_color==1) :
            self.name_colorTable = "autumn";
        elif(self.i_color==2) :
            self.name_colorTable = "bone";
        elif(self.i_color==3) :
            self.name_colorTable = "binary";
        elif(self.i_color==4) :
            self.name_colorTable = "bwr";
        elif(self.i_color==5) :
            self.name_colorTable = "brg";
        elif(self.i_color==6) :
            self.name_colorTable = "cool";
        elif(self.i_color==7) :
            self.name_colorTable = "copper";
        elif(self.i_color==8) :
            self.name_colorTable = "cubehelix";
        elif(self.i_color==9) :
            self.name_colorTable = "flag";
        elif(self.i_color==10) :
            self.name_colorTable = "gnuplot";
        elif(self.i_color==11) :
            self.name_colorTable = "gnuplot2";
        elif(self.i_color==12) :
            self.name_colorTable = "gray";
        elif(self.i_color==13) :
            self.name_colorTable = "hot";
        elif(self.i_color==14) :
            self.name_colorTable = "hsv";
        elif(self.i_color==15) :
            self.name_colorTable = "jet";
        elif(self.i_color==16) :
            self.name_colorTable = "ocean";
        elif(self.i_color==17) :
            self.name_colorTable = "pink";
        elif(self.i_color==18) :
            self.name_colorTable = "prism";
        elif(self.i_color==19) :
            self.name_colorTable = "rainbow";
        elif(self.i_color==20) :
            self.name_colorTable = "seismic";
        elif(self.i_color==21) :
            self.name_colorTable = "spring";
        elif(self.i_color==22) :
            self.name_colorTable = "summer";
        elif(self.i_color==23) :
            self.name_colorTable = "terrain";
        elif(self.i_color==24) :
            self.name_colorTable = "winter";
        elif(self.i_color==25) :
            self.name_colorTable = "spectral";
        elif(self.i_color==26) :
            self.name_colorTable = "redtoblue";
        elif(self.i_color==27) :
            self.name_colorTable = "redtoblue2";
        elif(self.i_color==28) :
            self.name_colorTable = "Accent";
        elif(self.i_color==29) :
            self.name_colorTable = "Blues";
        elif(self.i_color==30) :
            self.name_colorTable = "BrBG";
        elif(self.i_color==31) :
            self.name_colorTable = "BuGn";
        elif(self.i_color==32) :
            self.name_colorTable = "BuPu";
        elif(self.i_color==33) :
            self.name_colorTable = "Dark2";
        elif(self.i_color==34) :
            self.name_colorTable = "GnBu";
        elif(self.i_color==35) :
            self.name_colorTable = "Greens";
        elif(self.i_color==36) :
            self.name_colorTable = "Greys";
        elif(self.i_color==37) :
            self.name_colorTable = "Oranges";
        elif(self.i_color==38) :
            self.name_colorTable = "OrRd";
        elif(self.i_color==39) :
            self.name_colorTable = "Paired";
        elif(self.i_color==40) :
            self.name_colorTable = "Pastel1";
        elif(self.i_color==41) :
            self.name_colorTable = "Pastel2";
        elif(self.i_color==42) :
            self.name_colorTable = "PiYG"
        elif(self.i_color==43) :
            self.name_colorTable = "PRGn"
        elif(self.i_color==44) :
            self.name_colorTable = "PuBu"
        elif(self.i_color==45) :
            self.name_colorTable = "PuBuGn"
        elif(self.i_color==46) :
            self.name_colorTable = "PuOr"
        elif(self.i_color==47) :
            self.name_colorTable = "PuRd"
        elif(self.i_color==48) :
            self.name_colorTable = "Purples"
        elif(self.i_color==49) :
            self.name_colorTable = "RdBu"
        elif(self.i_color==50) :
            self.name_colorTable = "RdGy"
        elif(self.i_color==51) :
            self.name_colorTable = 'RdPu'
        elif(self.i_color==52) :
            self.name_colorTable = 'RdYlBu'
        elif(self.i_color==53) :
            self.name_colorTable = 'RdYlGn'
        elif(self.i_color==54) :
            self.name_colorTable = 'Reds'
        elif(self.i_color==55) :
            self.name_colorTable = 'Set1'
        elif(self.i_color==56) :
            self.name_colorTable = 'Set2'
        elif(self.i_color==57) :
            self.name_colorTable = 'Set3'
        elif(self.i_color==58) :
            self.name_colorTable = 'Spectral'
        elif(self.i_color==59) :
            self.name_colorTable = 'YlGn'
        elif(self.i_color==60) :
            self.name_colorTable = 'YlGnBu'
        elif(self.i_color==61) :
            self.name_colorTable = 'YlOrBr'
        elif(self.i_color==62) :
            self.name_colorTable = 'YlOrRd'
        elif(self.i_color==63) :
            self.name_colorTable = 'gist_earth'
        elif(self.i_color==64) :
            self.name_colorTable = 'gist_gray'
        elif(self.i_color==65) :
            self.name_colorTable = 'gist_heat'
        elif(self.i_color==66) :
            self.name_colorTable = 'gist_ncar'
        elif(self.i_color==67) :
            self.name_colorTable = 'gist_rainbow'
        elif(self.i_color==68) :
            self.name_colorTable = 'gist_stern'
        elif(self.i_color==69) :
            self.name_colorTable = 'gist_yarg'
        elif(self.i_color==70) :
            self.name_colorTable = 'coolwarm'


        print " My color table " + self.name_colorTable

        message = self.name_colorTable
        self.lineEdit_colorTable.clear();
	self.lineEdit_colorTable.insert(message);



    def LoadData(self) :
        # Load chosen data files from data handler
	#! Here we load the data corresponding to the choices in the widget. In
	#  paricular for the vector-fields we need to load several arrays,
	#  e.g., for the abs value

        if(self.string_has_patt(self.AllFields[self.i_quant],"Velocity")) :
            if(self.string_has_patt(self.Mods[self.i_mod],"Abs") or
               self.string_has_patt(self.Mods[self.i_mod],"Div") or
               self.string_has_patt(self.Mods[self.i_mod],"Curl")) :

                vx = self.ProjHandler.Load_2DData(self.i_proj, self.i_time,
                                                  self.NameIDs["v_x"],
                                                  self.i_plotPlane,
                                                  self.i_slice);
                vy = self.ProjHandler.Load_2DData(self.i_proj, self.i_time,
                                                  self.NameIDs["v_y"],
                                                  self.i_plotPlane,
                                                  self.i_slice);
                vz = self.ProjHandler.Load_2DData(self.i_proj, self.i_time,
                                                  self.NameIDs["v_z"],
                                                  self.i_plotPlane,
                                                  self.i_slice);
                data = np.copy(vx);

                if(self.string_has_patt(self.Mods[self.i_mod],"Abs")) :
                    data = np.sqrt(vx*vx + vy*vy + vz*vz)
                    PlotName = "v_Abs";
                elif(self.string_has_patt(self.Mods[self.i_mod],"Div") or
                     self.string_has_patt(self.Mods[self.i_mod],"Curl")) :
                    print " Still working on this "
                    sys.exit()
            else :
                # Use just a single field
                data = self.ProjHandler.Load_2DData(self.i_proj, self.i_time,
                                                    self.i_dataset,
                                                    self.i_plotPlane,
                                                    self.i_slice);
                PlotName = self.ProjHandler.get_FieldName(self.i_dataset);

        elif(self.string_has_patt(self.AllFields[self.i_quant],"Mag. field")) :
            if(self.string_has_patt(self.Mods[self.i_mod],"Abs") or
               self.string_has_patt(self.Mods[self.i_mod],"Div") or
               self.string_has_patt(self.Mods[self.i_mod],"Curl")) :

                Bx = self.ProjHandler.Load_2DData(self.i_proj, self.i_time,
                                                  self.NameIDs["B_x"],
                                                  self.i_plotPlane,
                                                  self.i_slice);
                By = self.ProjHandler.Load_2DData(self.i_proj, self.i_time,
                                                  self.NameIDs["B_y"],
                                                  self.i_plotPlane,
                                                  self.i_slice);
                Bz = self.ProjHandler.Load_2DData(self.i_proj, self.i_time,
                                                  self.NameIDs["B_z"],
                                                  self.i_plotPlane,
                                                  self.i_slice);
                data = np.copy(Bx);

                if(self.string_has_patt(self.Mods[self.i_mod],"Abs")) :
                    data = np.sqrt(Bx*Bx + By*By + Bz*Bz)
                    PlotName = "B_Abs";
                elif(self.string_has_patt(self.Mods[self.i_mod],"Div") or
                     self.string_has_patt(self.Mods[self.i_mod],"Curl")) :
                    print " Still working on this "
                    sys.exit()

            else :
                # Use just a single field
                data = self.ProjHandler.Load_2DData(self.i_proj, self.i_time,
                                                    self.i_dataset,
                                                    self.i_plotPlane,
                                                    self.i_slice);
                PlotName = self.ProjHandler.get_FieldName(self.i_dataset);
            
        else :
            # Use just a single field
            print " Reading "
            data = self.ProjHandler.Load_2DData(self.i_proj, self.i_time,
                                                self.i_dataset,
                                                self.i_plotPlane,
                                                self.i_slice);
            PlotName = self.ProjHandler.get_FieldName(self.i_dataset);
    


        # Start by getting the position
        if(self.i_plotPlane==0) :
            xPos = self.ProjHandler.get_PosArray(0);
            yPos = self.ProjHandler.get_PosArray(1);
            xName = "x";
            yName = "y";
        elif(self.i_plotPlane==1) :
            xPos = self.ProjHandler.get_PosArray(0);
            yPos = self.ProjHandler.get_PosArray(2);
            xName = "x";
            yName = "z";
        else :
            xPos = self.ProjHandler.get_PosArray(1);
            yPos = self.ProjHandler.get_PosArray(2);
            xName = "y";
            yName = "z";



        # Now store stuff in plot thingy
        Stash2D = PlotDataStorage.PlotData2D(xPos, yPos, data, xName,
                                             yName, PlotName)

        return Stash2D



    def do_plotContour(self) :
        # Finally make a countour plot
        
        self.WriteInfo(" Reading data ");
        

        stash2D = self.LoadData()


#        print " data "
#        print stash2D.data

        # Now plot the stuff:
        self.Plot2DStash(stash2D, stash2D.xAxisName(), stash2D.yAxisName())


        self.set_colorMap()


    def Plot2DStash(self, stash, xName, yName) :
        self.WriteInfo(" Doing actual 2D plot ")

#  http://scipy-lectures.github.io/intro/matplotlib/matplotlib.html#contour-plots

        if(self.has_figure) :
            plt.close()
            self.has_figure = False
        self.fig = plt.figure()
        self.has_figure = True
        plt.xlabel(xName)
        plt.ylabel(yName)
        plt.title("data")

        self.set_colorMap()

        if(self.geometry == 0) :
            is_polar = False
        elif(self.geometry == 1) :
            if(self.i_plotPlane == 0) :
                is_polar = True
            else :
                is_polar = False
        elif(self.geometry == 2) :
            if(self.i_plotPlane == 0 or self.i_plotPlane==1) :
                is_polar = True
            else :
                is_polar = False
#        ax = self.fig.add_subplot(111)
#        self.fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        if not (is_polar) :
            self.fig, ax = plt.subplots()
#            ax = self.fig.add_subplot(111)
        else :
            self.fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
            rad, theta = np.meshgrid(stash.xAxis, stash.yAxis)
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)


        num_pixels = len(stash.xAxis)*len(stash.yAxis)

        if(num_pixels > 100000 and self.cont_asCells) :
            print " pixels: " + str(num_pixels)
            self.butt_showSmooth.click();

        if(self.cont_asCells) :
            if(self.lin_plot) :
                if not (is_polar) :
                    im = ax.pcolor(stash.xAxis, stash.yAxis, stash.data,
                                   cmap=self.cmap)
                else :
                    im = ax.pcolor(theta, rad, stash.data,
                                   cmap=self.cmap)
            else :
                LogData = np.log(stash.data)
                if not (is_polar) :
                    im = ax.pcolor(stash.xAxis, stash.yAxis, LogData,
                                   cmap=self.cmap)
                else :
                    im = ax.pcolor(theta, rad, LogData,
                                   cmap=self.cmap)
            if not (is_polar) :
                ax.axis([np.min(stash.xAxis), np.max(stash.xAxis),
                         np.min(stash.yAxis), np.max(stash.yAxis)])
        else :
            if(self.lin_plot) :
                if not (is_polar) :
                    im = plt.contourf(stash.xAxis,stash.yAxis,stash.data,
                                      100,interpolation='linear',cmap=self.cmap)
                else :

                    im = ax.contourf(theta,rad,stash.data,30,cmap=self.cmap)
#                ax.axis([-np.max(stash.xAxis), np.max(stash.xAxis),
#                         -np.max(stash.xAxis), np.max(stash.xAxis)])
            else :
                # set levels manually
#                zmax = np.ceil(np.log10(stash.data.max())+1)
#                zmin = np.floor(np.log10(stash.data.min())-1)
#                dz = (zmax - zmin)/100.
#                lev_exp = np.arange(zmin, zmax, dz)
#                levs = np.power(10, lev_exp)
#                cs = P.contourf(X, Y, z, levs, norm=colors.LogNorm())
#                im = plt.contourf(stash.xAxis,stash.yAxis,stash.data,
#                                  levs,interpolation='linear',cmap=self.cmap)
                LogData = np.log(stash.data)
                if not (is_polar) :
                    im = plt.contourf(stash.xAxis,stash.yAxis,LogData,
                                      100,interpolation='linear',cmap=self.cmap)
                else :
                    im = plt.contourf(theta,rad,LogData,
                                      100,interpolation='linear',cmap=self.cmap)
#                im = plt.contourf(stash.xAxis,stash.yAxis,stash.data,
#                                  100,interpolation='linear',cmap=self.cmap,
#                                  norm=colors.LogNorm())


        ax.set_aspect('equal')

        delx = np.max(stash.xAxis) - np.min(stash.xAxis)
        dely = np.max(stash.yAxis) - np.min(stash.yAxis)
        if(delx/dely > 2.) :
            self.fig.colorbar(im, orientation='horizontal')
        else :
            self.fig.colorbar(im)

        self.show_Fig("Plot2D");

        self.WriteInfo(" done(actual 2D plot) ")




    def show_Fig(self, filename) :

        print self.fileFormatContour

        if(self.fileFormatContour == 'display') :
            self.WriteInfo(" Plotting on screen ")
            self.fig.show()
        else :
            print " in other plot "
            if(self.fileFormatContour == 'ps') :
                filename += '.ps'
            elif (self.fileFormatContour == 'png') :
                filename += '.png'
            elif (self.fileFormatContour == 'pdf') :
                filename += '.pdf'
            elif (self.fileFormatContour == 'svg') :
                filename += '.svg'
                
            self.WriteInfo(" Saving the file " + filename)
            self.fig.savefig(filename)
        

    def set_colorMap(self) :
        # set color map (still ugly):
        if(self.i_color==0) :
            self.cmap = plt.cm.afmhot
        elif(self.i_color==1) :
            self.cmap = plt.cm.autumn
        elif(self.i_color==2) :
            self.cmap = plt.cm.bone
        elif(self.i_color==3) :
            self.cmap = plt.cm.binary
        elif(self.i_color==4) :
            self.cmap = plt.cm.bwr
        elif(self.i_color==5) :
            self.cmap = plt.cm.brg
        elif(self.i_color==6) :
            self.cmap = plt.cm.cool
        elif(self.i_color==7) :
            self.cmap = plt.cm.copper
        elif(self.i_color==8) :
            self.cmap = plt.cm.cubehelix
        elif(self.i_color==9) :
            self.cmap = plt.cm.flag
        elif(self.i_color==10) :
            self.cmap = plt.cm.gnuplot
        elif(self.i_color==11) :
            self.cmap = plt.cm.gnuplot2
        elif(self.i_color==12) :
            self.cmap = plt.cm.gray
        elif(self.i_color==13) :
            self.cmap = plt.cm.hot
        elif(self.i_color==14) :
            self.cmap = plt.cm.hsv
        elif(self.i_color==15) :
            self.cmap = plt.cm.jet
        elif(self.i_color==16) :
            self.cmap = plt.cm.ocean
        elif(self.i_color==17) :
            self.cmap = plt.cm.pink
        elif(self.i_color==18) :
            self.cmap = plt.cm.prism
        elif(self.i_color==19) :
            self.cmap = plt.cm.rainbow
        elif(self.i_color==20) :
            self.cmap = plt.cm.seismic
        elif(self.i_color==21) :
            self.cmap = plt.cm.spring
        elif(self.i_color==22) :
            self.cmap = plt.cm.summer
        elif(self.i_color==23) :
            self.cmap = plt.cm.terrain
        elif(self.i_color==24) :
            self.cmap = plt.cm.winter
        elif(self.i_color==25) :
            self.cmap = plt.cm.spectral
        elif(self.i_color==26) :
            self.cmap = plt.cm.redtoblue
        elif(self.i_color==27) :
            self.cmap = plt.cm.redtoblue2
        elif(self.i_color==28) :
            self.cmap = plt.cm.Accent
        elif(self.i_color==29) :
            self.cmap = plt.cm.Blues
        elif(self.i_color==30) :
            self.cmap = plt.cm.BrBG
        elif(self.i_color==31) :
            self.cmap = plt.cm.BuGn
        elif(self.i_color==32) :
            self.cmap = plt.cm.BuPu
        elif(self.i_color==33) :
            self.cmap = plt.cm.Dark2
        elif(self.i_color==34) :
            self.cmap = plt.cm.GnBu
        elif(self.i_color==35) :
            self.cmap = plt.cm.Greens
        elif(self.i_color==36) :
            self.cmap = plt.cm.Greys
        elif(self.i_color==37) :
            self.cmap = plt.cm.Oranges
        elif(self.i_color==38) :
            self.cmap = plt.cm.OrRd
        elif(self.i_color==39) :
            self.cmap = plt.cm.Paired
        elif(self.i_color==40) :
            self.cmap = plt.cm.Pastel1
        elif(self.i_color==41) :
            self.cmap = plt.cm.Pastel2
        elif(self.i_color==42) :
            self.cmap = plt.cm.PiYG
        elif(self.i_color==43) :
            self.cmap = plt.cm.PRGn
        elif(self.i_color==44) :
            self.cmap = plt.cm.PuBu
        elif(self.i_color==45) :
            self.cmap = plt.cm.PuBuGn
        elif(self.i_color==46) :
            self.cmap = plt.cm.PuOr
        elif(self.i_color==47) :
            self.cmap = plt.cm.PuRd
        elif(self.i_color==48) :
            self.cmap = plt.cm.Purples
        elif(self.i_color==49) :
            self.cmap = plt.cm.RdBu
        elif(self.i_color==50) :
            self.cmap = plt.cm.RdGy
        elif(self.i_color==51) :
            self.cmap = plt.cm.RdPu
        elif(self.i_color==52) :
            self.cmap = plt.cm.RdYlBu
        elif(self.i_color==53) :
            self.cmap = plt.cm.RdYlGn
        elif(self.i_color==54) :
            self.cmap = plt.cm.Reds
        elif(self.i_color==55) :
            self.cmap = plt.cm.Set1
        elif(self.i_color==56) :
            self.cmap = plt.cm.Set2
        elif(self.i_color==57) :
            self.cmap = plt.cm.Set3
        elif(self.i_color==58) :
            self.cmap = plt.cm.Spectral
        elif(self.i_color==59) :
            self.cmap = plt.cm.YlGn
        elif(self.i_color==60) :
            self.cmap = plt.cm.YlGnBu
        elif(self.i_color==61) :
            self.cmap = plt.cm.YlOrBr
        elif(self.i_color==62) :
            self.cmap = plt.cm.YlOrRd
        elif(self.i_color==63) :
            self.cmap = plt.cm.gist_earth
        elif(self.i_color==64) :
            self.cmap = plt.cm.gist_gray
        elif(self.i_color==65) :
            self.cmap = plt.cm.gist_heat
        elif(self.i_color==66) :
            self.cmap = plt.cm.gist_ncar
        elif(self.i_color==67) :
            self.cmap = plt.cm.gist_rainbow
        elif(self.i_color==68) :
            self.cmap = plt.cm.gist_stern
        elif(self.i_color==69) :
            self.cmap = plt.cm.gist_yarg
        elif(self.i_color==70) :
            self.cmap = plt.cm.coolwarm



    # now the stuff for the line tab:
    def set_Dir1D_x(self, checked) :
        self.i_lineDir = 0
        self.label_slice1D_1.setText("y:");
        self.label_slice1D_2.setText("z:");

        if(self.butt_Dir1D_x.isEnabled()) :

            self.spinBox_slice1DPerp1.setRange(0,self.Nx[1]-1);
            self.spinBox_slice1DPerp2.setRange(0,self.Nx[2]-1);
            self.slider_slice1DPerp1.setRange(0,self.Nx[1]-1);
            self.slider_slice1DPerp2.setRange(0,self.Nx[2]-1);

            self.dir_Perp1 = 1
            self.dir_Perp2 = 2

            self.set_perpPos1(self.Nx[1]/2)
            self.set_perpPos2(self.Nx[2]/2)
            self.spinBox_slice1DPerp1.setValue(self.Nx[1]/2)
            self.spinBox_slice1DPerp2.setValue(self.Nx[2]/2)


    def set_Dir1D_y(self, checked) :
        self.i_lineDir = 1
        self.label_slice1D_1.setText("x:");
        self.label_slice1D_2.setText("z:");

        if(self.butt_Dir1D_y.isEnabled()) :

            self.spinBox_slice1DPerp1.setRange(0,self.Nx[0]-1);
            self.spinBox_slice1DPerp2.setRange(0,self.Nx[2]-1);
            self.slider_slice1DPerp1.setRange(0,self.Nx[0]-1);
            self.slider_slice1DPerp2.setRange(0,self.Nx[2]-1);
            
            self.dir_Perp1 = 0
            self.dir_Perp2 = 2
            
            self.set_perpPos1(self.Nx[0]/2)
            self.set_perpPos2(self.Nx[2]/2)
            self.spinBox_slice1DPerp1.setValue(self.Nx[0]/2)
            self.spinBox_slice1DPerp2.setValue(self.Nx[2]/2)

    def set_Dir1D_z(self, checked) :
        self.i_lineDir = 2
        self.label_slice1D_1.setText("x:");
        self.label_slice1D_2.setText("y:");

        if(self.butt_Dir1D_z.isEnabled()) :

            self.spinBox_slice1DPerp1.setRange(0,self.Nx[0]-1);
            self.spinBox_slice1DPerp2.setRange(0,self.Nx[1]-1);
            self.slider_slice1DPerp1.setRange(0,self.Nx[0]-1);
            self.slider_slice1DPerp2.setRange(0,self.Nx[1]-1);

            self.dir_Perp1 = 0
            self.dir_Perp2 = 1

            self.set_perpPos1(self.Nx[0]/2)
            self.set_perpPos2(self.Nx[1]/2)
            self.spinBox_slice1DPerp1.setValue(self.Nx[0]/2)
            self.spinBox_slice1DPerp2.setValue(self.Nx[1]/2)

    def set_perpPos1(self, i_PosPerp1) :
        print " setting perp pos " + str(i_PosPerp1)
        self.i_posPerp[0] = i_PosPerp1

        # Obtain corresponding position:
        # Now write a string holding the position

        if(self.i_lineDir == 0) :
            message = "y = "
            pos = self.ProjHandler.get_pos(1, self.i_posPerp[0]);
        else :
            message = "x = "
            pos = self.ProjHandler.get_pos(2, self.i_posPerp[0]);

        message += str(pos)
        self.label_slice1D_1.setText(message);

    def set_perpPos2(self, i_PosPerp2) :
        self.i_posPerp[1] = i_PosPerp2


        # Obtain corresponding position:
        if(self.i_lineDir == 2) :
            message = "y = "
            pos = self.ProjHandler.get_pos(1, self.i_posPerp[1]);
        else :
            message = "z = "
            pos = self.ProjHandler.get_pos(0, self.i_posPerp[1]);

        message += str(pos)
        self.label_slice1D_2.setText(message);


    def LoadData1D(self) :

        if(self.string_has_patt(self.AllFields[self.i_quant],"Velocity")) :
            if(self.string_has_patt(self.Mods[self.i_mod],"Abs") or
               self.string_has_patt(self.Mods[self.i_mod],"Div") or
               self.string_has_patt(self.Mods[self.i_mod],"Curl")) :

                vx = self.ProjHandler.Load_1DData(self.i_proj, self.i_time,
                                                  self.NameIDs["v_x"],
                                                  self.i_lineDir,
                                                  self.i_posPerp);
                vy = self.ProjHandler.Load_1DData(self.i_proj, self.i_time,
                                                  self.NameIDs["v_y"],
                                                  self.i_lineDir,
                                                  self.i_posPerp);
                vz = self.ProjHandler.Load_1DData(self.i_proj, self.i_time,
                                                  self.NameIDs["v_z"],
                                                  self.i_lineDir,
                                                  self.i_posPerp);
                data = np.copy(vx);

                if(self.string_has_patt(self.Mods[self.i_mod],"Abs")) :
                    data = np.sqrt(vx*vx + vy*vy + vz*vz)
                    PlotName = "v_Abs";
                elif(self.string_has_patt(self.Mods[self.i_mod],"Div") or
                     self.string_has_patt(self.Mods[self.i_mod],"Curl")) :
                    print " Still working on this "
                    sys.exit()
            else :
                # Use just a single field
                data = self.ProjHandler.Load_1DData(self.i_proj, self.i_time,
                                                    self.i_dataset,
                                                    self.i_lineDir,
                                                    self.i_posPerp);
                PlotName = self.ProjHandler.get_FieldName(self.i_dataset);

        elif(self.string_has_patt(self.AllFields[self.i_quant],"Mag. field")) :
            if(self.string_has_patt(self.Mods[self.i_mod],"Abs") or
               self.string_has_patt(self.Mods[self.i_mod],"Div") or
               self.string_has_patt(self.Mods[self.i_mod],"Curl")) :

                Bx = self.ProjHandler.Load_1DData(self.i_proj, self.i_time,
                                                  self.NameIDs["B_x"],
                                                  self.i_lineDir,
                                                  self.i_posPerp);
                By = self.ProjHandler.Load_1DData(self.i_proj, self.i_time,
                                                  self.NameIDs["B_y"],
                                                  self.i_lineDir,
                                                  self.i_posPerp);
                Bz = self.ProjHandler.Load_1DData(self.i_proj, self.i_time,
                                                  self.NameIDs["B_z"],
                                                  self.i_lineDir,
                                                  self.i_posPerp);
                data = np.copy(Bx);

                if(self.string_has_patt(self.Mods[self.i_mod],"Abs")) :
                    data = np.sqrt(Bx*Bx + By*By + Bz*Bz)
                    PlotName = "B_Abs";
                elif(self.string_has_patt(self.Mods[self.i_mod],"Div") or
                     self.string_has_patt(self.Mods[self.i_mod],"Curl")) :
                    print " Still working on this "
                    sys.exit()

            else :
                # Use just a single field
                data = self.ProjHandler.Load_1DData(self.i_proj, self.i_time,
                                                    self.i_dataset,
                                                    self.i_lineDir,
                                                    self.i_posPerp); 
                PlotName = self.ProjHandler.get_FieldName(self.i_dataset);
            
        else :
            # Use just a single field
            print " Reading "
            data = self.ProjHandler.Load_1DData(self.i_proj, self.i_time,
                                                self.i_dataset,
                                                self.i_lineDir,
                                                self.i_posPerp); 
            PlotName = self.ProjHandler.get_FieldName(self.i_dataset);

        # Now obtain position array:
        if(self.i_lineDir==0) :
            xPos = self.ProjHandler.get_PosArray(0);
            xName = "x";
        elif(self.i_lineDir==0) :
            xPos = self.ProjHandler.get_PosArray(1);
            xName = "y";
        else :
            xPos = self.ProjHandler.get_PosArray(2);
            xName = "z";

        yName = " Data "
        PlotName = " Line Plot "


        # Now store stuff in stash
        Stash1D = PlotDataStorage.PlotData1D(xPos, data, xName,
                                             yName, PlotName)

        return Stash1D

    def do_plotLine(self) :
        self.WriteInfo(" Doing a line plot ");

        stash1D = self.LoadData1D()
        
        # Now plot the data
        self.Plot1DStash(stash1D, stash1D.xAxisName(), stash1D.yAxisName())

        self.WriteInfo(" done (line plot) ");


    def Plot1DStash(self, stash1D, xName, yName) :
        self.WriteInfo(" Doing actual line plot ")

        # get axis extreme values:
        xMin = stash1D.xAxis.min()
        xMax = stash1D.xAxis.max()

        yMin = stash1D.data.min()
        yMax = stash1D.data.max()

        print yMax
        dely = yMax - yMin

        yAxisMax = yMax + 0.1*dely
        yAxisMin = yMin - 0.1*dely


        if(self.has_figure) :
            plt.close()
            self.has_figure = False
        self.fig = plt.figure()
        self.has_figure = True
        plt.xlabel(xName)
        plt.ylabel(yName)
        plt.title("line plot")
        plt.axis([xMin, xMax, yAxisMin, yAxisMax])


        ax = self.fig.add_subplot(111)
        ax.plot(stash1D.xAxis, stash1D.data, 'ro') #  ro means - red circles
        
        self.show_Fig("LinePlot")


    def WriteInfo(self, message) :
        self.info_list.addItem(message)
        self.info_list.scrollToBottom()
            

    def show_DataProperties(self) :

	self.list_dataProperties.clear();

        infoLine = "Resolution: ";
	self.list_dataProperties.insertItem(0, infoLine);
	infoLine = "  [" + str(self.Nx[0]) + " x " + str(self.Nx[1])
        infoLine += " x " + str(self.Nx[2]) + "]";
	self.list_dataProperties.insertItem(1, infoLine);

        xPos = self.ProjHandler.get_PosArray(0);
	yPos = self.ProjHandler.get_PosArray(1);
	zPos = self.ProjHandler.get_PosArray(2);

	infoLine = "Extent: ";
	self.list_dataProperties.insertItem(3, infoLine);
	
        infoLine = " [" + str(xPos[0]) + "x" + str(xPos[self.Nx[0]-1]) + "]"
	self.list_dataProperties.insertItem(4, infoLine);
        infoLine = " [" + str(yPos[0]) + "x" + str(yPos[self.Nx[1]-1]) + "]"
	self.list_dataProperties.insertItem(5, infoLine);
        infoLine = " [" + str(zPos[0]) + "x" + str(zPos[self.Nx[2]-1]) + "]"
	self.list_dataProperties.insertItem(6, infoLine);


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

        
