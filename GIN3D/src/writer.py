import gi
import re
from utilities import *
from gi.repository import Gtk

gi.require_version('Gtk', '3.0')


class GIN3DConfigWriter:
  
    def __init__(self, builder):
        self.cnfg = {}
        self.win = builder
        self.required = {
	  'domain_lx' : 'Domain Size Lx', 
	  'domain_ly': 'Domain Size Ly', 
	  'domain_lz': 'Domain Size Lz', 
	  'mesh_nx': 'Mesh Size Nx', 
	  'mesh_ny': 'Mesh Size Ny', 
	  'mesh_nz': 'Mesh Size Nz',
	  'refLength' : 'Reference Length', 
	  'laminar_viscosity' : 'Laminar Viscosity', 
	  'sim_endTime' : 'Simulation End Time',
	  'inletProfile_u' : 'Inlet Profile Parameter u',
	  'inletProfile_v' : 'Inlet Profile Parameter v',
	  'inletProfile_w' : 'Inlet Profile Parameter w',
	  'ic_uniform_u' : 'Initial Condition u',
	  'ic_uniform_v' : 'Initial Condition v',
	  'ic_uniform_w' : 'Initial Condition w',
	  'fsp_text_momentum_hybrid' : 'Flow Solver Parameter Hybrid',
	  'fsp_text_temperature_hybrid' : 'Temperature Equation Hybrid'
	  }

    def to_file(self, filename):
        #print("Filename: ", filename)
        
        # initialize that all required values are present
        self.gotAllReqValues = True
        self.reqValuesLog = ''
        
        self.get_data()
        
        if not self.gotAllReqValues:
	  log(self.reqValuesLog)
	  return self.reqValuesLog
	
        self.writeToFile(filename)

    def text(self, obj_id, evaluate=False):
        obj = self.win.get_object(obj_id)
        if obj is not None:
            if isinstance(obj, Gtk.ComboBoxText):
                return obj.get_active_id()
            elif isinstance(obj, Gtk.FileChooserButton):
                return obj.get_filename()

            else:
                val = obj.get_text()
                if val != '' and evaluate:
                    return str(eval(val))
		if obj_id in self.required and val.strip() == "":
		    self.gotAllReqValues = False
		    self.reqValuesLog += 'Missing... ' + self.required[obj_id] + '\n'
                return val
        return ''

    def active(self, obj_id):
        return self.win.get_object(obj_id).get_active()

    def get_obj(self, obj_id):
        return self.win.get_object(obj_id)

    # Read Data from the GUI
    def get_data(self):
        
        self.cnfg['Main'] = self.get_main_parameters()
        self.cnfg['BC'] = self.get_bc_parameters()
        self.cnfg['IC'] = self.get_ic_parameters()
        self.cnfg['FS'] = self.get_fs_parameters()
        self.cnfg['SO'] = self.get_so_parameters()
        self.cnfg['TP'] = self.get_tp_parameters()
        self.cnfg['TM'] = self.get_tm_parameters()
        self.cnfg['GE'] = self.get_ge_parameters()

    def get_main_parameters(self):
        cnfg_main = {'Dimensions': "%s %s %s" % (self.text("domain_lx"), self.text("domain_ly"), self.text("domain_lz")),
                     'Grid': "%s %s %s" % (self.text("mesh_nx"), self.text("mesh_ny"), self.text("mesh_nz")),
                     'ReferenceLength': self.text("refLength", True), 'Nu': self.text("laminar_viscosity", True),
                     'StopTime': self.text("sim_endTime"), 'Temperature': self.active("simParam_switch_temperature"),
                     'Turbulence': self.active("simParam_switch_turbulence"),
                     'SolidGeometry': self.active("simParam_switch_geometry")}

        self.forcing = "off"  # variable set to check when writing to file
        forcing_off = self.active("simParam_rb_forcing_off")

        if not forcing_off:
            forcing_cpg = self.active("simParam_rb_forcing_constPresGrad")
            if forcing_cpg:
                self.forcing = "constPresGrad"
                cnfg_main['Forcing'] = "%s %s %s" % (
                    self.text("simParam_txtbox_fx", True), self.text("simParam_txtbox_fy", True),
                    self.text("simParam_txtbox_fz", True))
            else:
                self.forcing = "constMassFlowRate"
                cnfg_main['NumberForcingRegions'] = self.text("no_of_forcing_regions")
                numforcings = int(cnfg_main['NumberForcingRegions'])
                heights = [0.0]
                mdot = []
                for i in range(numforcings):
                    h = first_child("simParam_forcing_cmfr_box", "height_", i, self.win)
                    if h is not None:
                        heights.append(eval(h.get_text()))
                        tmp = []
                        x = first_child("simParam_forcing_cmfr_box", "frfr_x_", i, self.win)
                        if x is not None:
                            tmp.append(str(eval(x.get_text())))

                        y = first_child("simParam_forcing_cmfr_box", "frfr_y_", i, self.win)
                        if y is not None:
                            tmp.append(str(eval(y.get_text())))

                        z = first_child("simParam_forcing_cmfr_box", "frfr_z_", i, self.win)
                        if z is not None:
                            tmp.append(str(eval(z.get_text())))

                        mdot.append(tmp)

                cnfg_main['ForcingRegionBounds'] = heights
                cnfg_main['ForcingRegionFlowRate'] = mdot

        return cnfg_main

    def get_bc_parameters(self):
        cnfg_bc = {'Face_West': self.text("west_face"), 'Face_East': self.text("east_face"),
                   'Face_South': self.text("south_face"), 'Face_North': self.text("north_face"),
                   'Face_Bottom': self.text("bottom_face"), 'Face_Top': self.text("top_face")}
	
      
        if "Inlet" in [cnfg_bc['Face_South'], cnfg_bc['Face_West']]:                      
	    cnfg_bc['InletVelocities'] = ("%s %s %s") %(self.text("inletProfile_u"), self.text("inletProfile_v"), self.text("inletProfile_w"))
	    cnfg_bc['VariableProfile'] = self.active("inletProfile_switch_variable")
	    cnfg_bc['PerturbTurbInflow'] = self.active("inletProfile_switch_turbulent")
      
	    if cnfg_bc['PerturbTurbInflow']:
		cnfg_bc['PerturbZoneOffset'] = self.text('inletProfile_offset')
		cnfg_bc['PerturbZoneRichardson'] = self.text('inletProfile_richardson')
		cnfg_bc['PerturbZoneDepth'] = self.text('inletProfile_depth')
		cnfg_bc['PerturbBoxLength'] = self.text('inletProfile_length')
		cnfg_bc['PerturbBoxWidth'] = self.text('inletProfile_width')
		cnfg_bc['PerturbBoxHeight'] = self.text('inletProfile_height')

	    if cnfg_bc['VariableProfile']:
		val = "Inlet"
		if self.active("rb_inletProfile_roughWall"):
		    val = "RoughLogLawInlet"
		elif self.active("rb_inletProfile_smoothWall"):
		    val = "SmoothLogLawInlet"
		elif self.active("rb_inletProfile_powerLaw"):
		    val = "OneSeventhsInlet"
		elif self.active("rb_inletProfile_parabolic"):
		    val = "Parabolic"
		elif self.active("rb_inletProfile_dns"):
		    val = "DNS2000"

		if cnfg_bc['Face_South'] == "Inlet":
		    cnfg_bc['Face_South'] = val

		if cnfg_bc['Face_West'] == "Inlet":
		    cnfg_bc['Face_West'] = val

        return cnfg_bc

    def get_ic_parameters(self):
        cnfg_ic = {}

        if self.active('ic_rb_sameAsInletProfile'):
            # removed 'Parabolic' since it has not been implemented
            inlet_values = ['Inlet', 'RoughLogLawInlet', 'SmoothLogLawInlet', 'OneSeventhsInlet', 'DNS2000']

            cnfg_bc = self.cnfg['BC']
            faceW = cnfg_bc['Face_West']
            faceS = cnfg_bc['Face_South']

            if faceW in inlet_values or faceS in inlet_values:
                if faceW == "Inlet" or faceS == "Inlet":
                    cnfg_ic['InitialCondition'] = 'Uniform ' + self.cnfg['BC']['InletVelocities']
                elif faceW == "RoughLogLawInlet" or faceS == "RoughLogLawInlet":
                    cnfg_ic['InitialCondition'] = 'ComplexTerrain'
                elif faceW == "SmoothLogLawInlet" or faceS == "SmoothLogLawInlet":
                    cnfg_ic['InitialCondition'] = 'TurbulentChannelXMean'
                    # Don't have a power initial condition implemented yet, so we stick with smooth wall log law
                elif faceW == "OneSeventhsInlet" or faceS == "OneSeventhsInlet":
                    cnfg_ic['InitialCondition'] = 'TurbulentChannelXMean'
                    # Same with DNS2000
                elif faceW == "DNS2000" or faceS == "DNS2000":
                    cnfg_ic['InitialCondition'] = 'TurbulentChannelXMean'
            else:
                self.get_obj('ic_rb_sameAsInletProfile').set_active('False')
        elif self.active('ic_rb_zeroField'):
            cnfg_ic['InitialCondition'] = 'ZeroField'
        elif self.active('ic_rb_uniform'):
            cnfg_ic['InitialCondition'] = 'Uniform ' + self.text('ic_uniform_u') + " " + self.text(
                'ic_uniform_v') + " " + self.text('ic_uniform_w')

        elif self.active('ic_rb_specialCases'):
            special_case = self.text('ic_list_specialCases')
            # Smooth Wall Log Law w/Sinusoidal Perturbations
            if special_case == '1':
                cnfg_ic['InitialCondition'] = 'ComplexTerrainSinPerturb'
            # Rough Wall Log Law w/Sinusoidal Perturbations
            elif special_case == '2':
                cnfg_ic['InitialCondition'] = 'TurbulentChannelX'
            # ABC Flow
            elif special_case == '3':
                cnfg_ic['InitialCondition'] = 'ABCFlow'

        return cnfg_ic

    def get_fs_parameters(self):

        cnfg_fs = {'AdvectionScheme': '0.00'}

        # Advection Scheme

        # Momentum
        if self.active('fsp_rb_momentum_cds'):
            cnfg_fs['AdvectionScheme'] = '0.00'
        elif self.active('fsp_rb_momentum_fou'):
            cnfg_fs['AdvectionScheme'] = '1.00'
        elif self.active('fsp_rb_momentum_quick'):
            cnfg_fs['AdvectionScheme'] = '-1.00'
        elif self.active('fsp_rb_momentum_kappa'):
            cnfg_fs['AdvectionScheme'] = ''
        elif self.active('fsp_rb_momentum_hybrid'):
            cnfg_fs['AdvectionScheme'] = self.text('fsp_text_momentum_hybrid')

        # Time Marching
        cnfg_fs['TimeMethod'] = self.text('fsp_timeMarching')

        # Time Step Size
        cnfg_fs['TimeStepSize'] = self.text('fsp_timeStepSize')

        if cnfg_fs['TimeStepSize'] == 'Variable':
            cnfg_fs['ConstTimeStep'] = '-0.1'
            cnfg_fs['CFL'] = self.text('fsp_cfl')
            cnfg_fs['DTStability'] = self.text('fsp_viscousLimit')
        else:
            cnfg_fs['ConstTimeStep'] = self.text('fsp_dt')
            cnfg_fs['CFL'] = '0.4'
            cnfg_fs['DTStability'] = '0.4'

        # Poisson Solver
        cnfg_fs['SolverMethod'] = self.text('fsp_poissonSolver')

        if cnfg_fs['SolverMethod'] == 'Multigrid':

            cnfg_fs['MG_Loops'] = self.text('fsp_text_poissonSolver_outerLoops')
            cnfg_fs['MG_Levels'] = self.text('fsp_text_poissonSolver_mgLevels')
            cnfg_fs['MG_Smoother_Iterations'] = [self.text('fsp_text_poissonSolver_innerLoops_down'),
                                                 self.text('fsp_text_poissonSolver_innerLoops_up')]
            cnfg_fs['MG_Jacobi_Weight'] = self.text('fsp_text_poissonSolver_mg_jacobiWeight')

        else:
            cnfg_fs['Iterative_Loops'] = self.text('fsp_text_poissonSolver_iterations')
            cnfg_fs['Iterative_Weight'] = self.text('fsp_text_poissonSolver_pj_jacobiWeight')
        return cnfg_fs

    def get_so_parameters(self):
        cnfg_so = {'OutputPath': self.text('so_folder_outputPath'),
                   'OutputPrefix': self.text('so_text_outputPrefix')
                   }

        # Restart Solution
        if self.active('so_switch_restartSolution'):
            cnfg_so['RestartDirectory'] = self.text('so_folder_restartDirectory') + '/'

        outputFormat = []

        # Formatted Output
        if self.active('so_cb_fileFormat_hdf5'):
            outputFormat.append('HDF5_Multiple')

        # need to confirm this
        if self.active('so_cb_fileFormat_vti'):
            outputFormat.append('VTK') 

        if self.active('so_cb_fileFormat_hdf5') or self.active('so_cb_fileFormat_vti'):

            if self.active('so_cb_formatOutput_timeStep'):
                cnfg_so['OutputFrequencyVTK1'] = self.text('so_text_formatOutput_timeStep_interval1')
                if self.active('so_cb_formatOutput_timeStep_switch'):
                    cnfg_so['OutputFrequencyVTK2'] = self.text('so_text_formatOutput_timeStep_interval2')
                    cnfg_so['OutputFrequencyVTKSwitch'] = self.text('so_text_formatOutput_timeStep_switch')

            if self.active('so_cb_formatOutput_physicalTime'):
                cnfg_so['OutputTimeIntervalVTK1'] = self.text('so_text_formatOutput_physicalTime_interval1')
                if self.active('so_cb_formatOutput_physicalTime_switch'):
                    cnfg_so['OutputTimeIntervalVTK2'] = self.text('so_text_formatOutput_physicalTime_interval2')
                    cnfg_so['OutputTimeIntervalVTKSwitch'] = self.text('so_text_formatOutput_physicalTime_switch')

            outputVariablesVTK = []

            if self.active('so_cb_outputVTK_u'):
                outputVariablesVTK.append('u')
            if self.active('so_cb_outputVTK_p'):
                outputVariablesVTK.append('p')
            if self.active('so_cb_outputVTK_t'):
                outputVariablesVTK.append('t')
            if self.active('so_cb_outputVTK_n'):
                outputVariablesVTK.append('n')
            if self.active('so_cb_outputVTK_z'):
                outputVariablesVTK.append('z')
            if self.active('so_cb_outputVTK_q'):
                outputVariablesVTK.append('q')
            if self.active('so_cb_outputVTK_g'):
                outputVariablesVTK.append('g')
            if self.active('so_cb_outputVTK_c'):
                outputVariablesVTK.append('c')
            if self.active('so_cb_outputVTK_i'):
                outputVariablesVTK.append('i')

            cnfg_so['OutputVariablesVTK'] = outputVariablesVTK

        # Raw Data
        if self.active('so_cb_rawData_binary'):
            outputFormat.append('Binary')
        elif self.active('so_cb_rawData_parallelBinary'):
            outputFormat.append('PBinary')

        if self.active('so_cb_rawData_ASCII'):
            outputFormat.append('VTK')
        elif self.active('so_cb_rawData_parallelASCII'):
            outputFormat.append('PVTK')

        if len(set(outputFormat).intersection(set(['Binary', 'PBinary', 'VTK', 'PVTK']))) > 0:
            if self.active('so_cb_rawData_timeStep'):
                cnfg_so['OutputFrequencyText1'] = self.text('so_text_rawData_timeStep_interval1')
                if self.active('so_cb_rawData_timeStep_switch'):
                    cnfg_so['OutputFrequencyText2'] = self.text('so_text_rawData_timeStep_interval2')
                    cnfg_so['OutputFrequencyTextSwitch'] = self.text('so_text_rawData_timeStep_switch')

            if self.active('so_cb_rawData_physicalTime'):
                cnfg_so['OutputTimeIntervalText1'] = self.text('so_text_rawData_physicalTime_interval1')
                if self.active('so_cb_formatOutput_physicalTime_switch'):
                    cnfg_so['OutputTimeIntervalText2'] = self.text('so_text_rawData_physicalTime_interval2')
                    cnfg_so['OutputTimeIntervalTextSwitch'] = self.text('so_text_rawData_physicalTime_switch')

            outputVariablesTXT = []

            if self.active('so_cb_outputTXT_u'):
                outputVariablesTXT.append('u')
            if self.active('so_cb_outputTXT_p'):
                outputVariablesTXT.append('p')
            if self.active('so_cb_outputTXT_t'):
                outputVariablesTXT.append('t')
            if self.active('so_cb_outputTXT_n'):
                outputVariablesTXT.append('n')
            if self.active('so_cb_outputTXT_z'):
                outputVariablesTXT.append('z')
            if self.active('so_cb_outputTXT_q'):
                outputVariablesTXT.append('q')
            if self.active('so_cb_outputTXT_g'):
                outputVariablesTXT.append('g')
            if self.active('so_cb_outputTXT_c'):
                outputVariablesTXT.append('c')
            if self.active('so_cb_outputTXT_i'):
                outputVariablesTXT.append('i')

            cnfg_so['OutputVariablesTXT'] = outputVariablesTXT

        # Screen Dump
        if self.active('so_cb_screenDump'):
            outputFormat.append('Matrix')

        # Output Plane
        if self.active('so_cb_screenDump'):
	  if self.active('so_rb_screenDump_xy'):
	      cnfg_so['OutputPlane'] = 'XY'
	      cnfg_so['OutputSlice'] = self.text('so_text_screenDump_xy')
	  elif self.active('so_rb_screenDump_xz'):
	      cnfg_so['OutputPlane'] = 'XZ'
	      cnfg_so['OutputSlice'] = self.text('so_text_screenDump_xz')
	  elif self.active('so_rb_screenDump_yz'):
	      cnfg_so['OutputPlane'] = 'YZ'
	      cnfg_so['OutputSlice'] = self.text('so_text_screenDump_yz')

	  if self.active('so_rb_screenDump_u'):
	      cnfg_so['MatrixData'] = 'u'
	  elif self.active('so_rb_screenDump_v'):
	      cnfg_so['MatrixData'] = 'v'
	  elif self.active('so_rb_screenDump_w'):
	      cnfg_so['MatrixData'] = 'w'
	  elif self.active('so_rb_screenDump_p'):
	      cnfg_so['MatrixData'] = 'p'
	  elif self.active('so_rb_screenDump_temperature'):
	      cnfg_so['MatrixData'] = 'phi'

        cnfg_so['OutputFormat'] = outputFormat

        # Sampling
        cnfg_so['TrackInterval'] = self.text('so_text_physicalTime_interval')
        coord_file = self.text('so_file_coordFile')

        cnfg_so['TrackCoordFile'] = coord_file
        if coord_file is None:
            cnfg_so['TrackCoordinates'] = 0
        else:
            # parse the coord file
            count = 0;
            trackCoord = []
            with open(coord_file, "r") as fp:
                for line in fp:
                    l = re.split('[, ]', line.strip())
                    trackCoord.append(l)
                    count += 1

            cnfg_so['TrackCoordinates'] = count
            if count > 0:
                cnfg_so['TrackCoord'] = trackCoord

        # Time Series OPtions
        # Collect Statistics
        cnfg_so['TimeAvg'] = str(self.active('so_switch_collectStatistics'))

        # Sampling Time Begin
        cnfg_so['StartingTime'] = self.text('so_text_samplingTimeBegin')

        # Spatial Averaging
        spatialAvg = self.active('so_switch_spatialAveraging')
        # -- need to review this
        turbStats = []
        if self.active('so_rb_firstOrderStatistics'):
            turbStats.append('m')
        else:
            if self.active('so_cb_stats_s'):
                turbStats.append('s')
            if self.active('so_cb_stats_d'):
                turbStats.append('d')
            if self.active('so_cb_stats_i'):
                turbStats.append('i')

        if spatialAvg:
            turbStats = [e.upper() for e in turbStats]

        cnfg_so['TurbStats'] = turbStats

        # Mean quantity files for fluctuation calculations
        if self.active('so_cb_mean_velocity'):
            cnfg_so['MeanProfile_U'] = self.text('so_file_meanprofile_u')
            cnfg_so['MeanProfile_V'] = self.text('so_file_meanprofile_v')
            cnfg_so['MeanProfile_W'] = self.text('so_file_meanprofile_w')

        if self.active('so_cb_mean_temperature'):
            cnfg_so['MeanProfile_Phi'] = self.text('so_file_meanprofile_phi')

        if self.active('so_cb_mean_pressure'):
            cnfg_so['MeanProfile_P'] = self.text('so_file_meanprofile_p')

        return cnfg_so

    def get_tp_parameters(self):
        cnfg_tp = {}
        if self.cnfg['Main']['Temperature']:
            cnfg_tp = {'Temp_West': [self.text('temp_west_face'), self.text('temp_text_westFace')],
                       'Temp_East': [self.text('temp_east_face'), self.text('temp_text_eastFace')],
                       'Temp_North': [self.text('temp_north_face'), self.text('temp_text_northFace')],
                       'Temp_South': [self.text('temp_south_face'), self.text('temp_text_southFace')],
                       'Temp_Top': [self.text('temp_top_face'), self.text('temp_text_topFace')],
                       'Temp_Bottom': [self.text('temp_bottom_face'), self.text('temp_text_bottomFace')],
                       'Gravity': [self.text('temp_text_gravity_x'), self.text('temp_text_gravity_y'),
                                   self.text('temp_text_gravity_z')],
                       'TurbulentPrandtl': self.text('temp_text_prandtlNumber'),
                       'Temp_Infinity': self.text('temp_text_refTemperature'),
                       'Beta': self.text('temp_text_thermalExpansionCoefficient'),
                       'Gamma': self.text('temp_text_thermalDiffusivity'),
                       'Source_PHI': self.text('temp_text_sourceMagnitude'),
                       'TemperatureAdvectionScheme': '0.00'}

            # Temperature
            if self.active('fsp_rb_temperature_cds'):
                cnfg_tp['TemperatureAdvectionScheme'] = '0.00'
            elif self.active('fsp_rb_temperature_fou'):
                cnfg_tp['TemperatureAdvectionScheme'] = '1.00'
            elif self.active('fsp_rb_temperature_quick'):
                cnfg_tp['TemperatureAdvectionScheme'] = '-1.00'
            elif self.active('fsp_rb_temperature_kappa'):
                cnfg_tp['TemperatureAdvectionScheme'] = ''
            elif self.active('fsp_rb_temperature_hybrid'):
                cnfg_tp['TemperatureAdvectionScheme'] = self.text('fsp_text_temperature_hybrid')

        return cnfg_tp

    def get_tm_parameters(self):
        cnfg_tm = {}
        # Turbulence Model
        # LES
        if self.cnfg['Main']['Turbulence']:
            if self.active('tm_rb_originalSmagorinsky'):
                cnfg_tm['TurbulenceModel'] = 'Smagorinsky'
            elif self.active('tm_rb_lagrangianDynamicSmagorinsky'):
                cnfg_tm['TurbulenceModel'] = 'LagDynamic'
            elif self.active('tm_rb_prandtlSmagorinsky'):
                cnfg_tm['TurbulenceModel'] = 'SmagRANS'
            elif self.active('tm_rb_prandtlLagrangianDynamicSmagorinsky'):
                cnfg_tm['TurbulenceModel'] = 'LagDynRANS'

            # Subgrid Scale Parameters
            # Smagorinsky Coefficient
            cnfg_tm['TurbCS'] = self.text('tm_text_smagorinskyCoefficient')
            # RANS-LES Blending Height
            cnfg_tm['TransitionHeight'] = self.text('tm_text_blendingHeight')

            # Distance Field
            distance_field_opt = self.text('tm_list_distanceField')
            if distance_field_opt == 0:
                file = self.text('tm_file_distanceField')
                if file is not None:
                    cnfg_tm['DistanceField'] = file
                else:
                    cnfg_tm['DistanceField'] = ''
            elif distance_field_opt == 1:
                if self.active('tm_rb_channel_x'):
                    cnfg_tm['DistanceField'] = 'TPC_X_WallNormal'
                elif self.active('tm_rb_channel_y'):
                    cnfg_tm['DistanceField'] = 'TPC_Y_WallNormal'
                elif self.active('tm_rb_channel_z'):
                    cnfg_tm['DistanceField'] = 'TPC_Z_WallNormal'
        else:
            cnfg_tm['TurbulenceModel'] = 'Laminar'
        return cnfg_tm

    def get_ge_parameters(self):
        cnfg_ge = {}
        if self.cnfg['Main']['SolidGeometry']:
            ib_file = self.text('geom_file_ibfiles')

            count = 0;
            # parse the IB Files
            if ib_file is not None:
                with open(ib_file, "r") as fp:
                    for line in fp:
                        l = line.strip().replace('\t', '')
                        if not (l.startswith('#') or '' == l):
                            if '_ibnode_s' in l:
                                cnfg_ge['IBNodeFile_S'] = l
                                count += 1
                            elif '_ibnode_u' in l:
                                cnfg_ge['IBNodeFile_U'] = l
                                count += 1
                            elif '_ibnode_v' in l:
                                cnfg_ge['IBNodeFile_V'] = l
                                count += 1
                            elif '_ibnode_w' in l:
                                cnfg_ge['IBNodeFile_W'] = l
                                count += 1
                            elif '_flag_s' in l:
                                cnfg_ge['IBFlagFile_S'] = l
                                count += 1
                            elif '_flag_u' in l:
                                cnfg_ge['IBFlagFile_U'] = l
                                count += 1
                            elif '_flag_v' in l:
                                cnfg_ge['IBFlagFile_V'] = l
                                count += 1
                            elif '_flag_w' in l:
                                cnfg_ge['IBFlagFile_W'] = l
                                count += 1
                if count != 8:
                    log('Error parsing, please review:  ' + ib_file, 'e')
                    ib_file = None

            cnfg_ge['IBFiles'] = ib_file

            # IB Reconstruction Scheme
            if self.active('geom_rb_ib_linear'):
                cnfg_ge['IBReconstruction'] = 'Linear'
            elif self.active('geom_rb_ib_roughLogLaw'):
                cnfg_ge['IBReconstruction'] = 'RoughLogLaw'
                cnfg_ge['Roughness'] = self.text('geom_text_roughness')
            elif self.active('geom_rb_ib_powerLaw'):
                cnfg_ge['IBReconstruction'] = 'OneSevenths'

            # Distance Field
            distance_field_opt = self.text('geom_list_distanceField')
            if distance_field_opt == 0:
                file = self.text('geom_file_distanceField')
                if file is not None:
                    cnfg_ge['DistanceField'] = file
                else:
                    cnfg_ge['DistanceField'] = ''
            elif distance_field_opt == 1:
                if self.active('geom_rb_channel_x'):
                    cnfg_ge['DistanceField'] = 'TPC_X_WallNormal'
                elif self.active('geom_rb_channel_y'):
                    cnfg_ge['DistanceField'] = 'TPC_Y_WallNormal'
                elif self.active('geom_rb_channel_z'):
                    cnfg_ge['DistanceField'] = 'TPC_Z_WallNormal'
        return cnfg_ge
    
    def write(self, tag, value=None):
      
      newlineOnly = re.compile('^\n+$', re.DOTALL)
      if tag.startswith("#") or newlineOnly.match(tag) != None:
	self.fp.write(tag)
      elif value is not None:
	val = value.strip();
	if len(val) > 0:
	  self.fp.write(tag + value)
	else:
	  self.fp.write('#' + tag + value)
      
      
    # Write Data to File
    def writeToFile(self, filename):

        try:
            self.fp = open(filename.replace('.cfg', '') + ".cfg", "w")

            self.write("################################################################################\n")
            self.write("# Parameters to control the model input/output\n")
            self.write("# Generated from GIN3D Configuration GUI\n")
            self.write("################################################################################\n")
            self.write("\n")
            self.write("# The # comments out the line.  If a parameter is commented, its value is set to\n")
            self.write("# zero.\n")
            self.write("# NOTE: The density is considered 1.0 for all simulations.\n")
            self.write("\n")

            # Main Parameters
            cnfg_main = self.cnfg["Main"]

            self.write("#-------------------------------------------------------------------------------\n")
            self.write("#                          Main Parameters\n")
            self.write("#-------------------------------------------------------------------------------\n")
            self.write("\n")

            # Item 1
            self.write("# Physical domain size, requires three floating point arguments after\n")
            self.write("# Dimensions, separated by spaces.\n")
            self.write("Dimensions ", cnfg_main["Dimensions"])
            self.write("\n\n")

            # Item 2
            self.write("# Simple mesh size, requires three integer arguments after Grid, separated by\n")
            self.write("# spaces.\n")
            self.write("Grid ", cnfg_main["Grid"])
            self.write("\n\n")

            # Item 3
            self.write("# Reference Length\n")
            self.write("# This is used as the characteristic length 'L' in calculations of the\n")
            self.write("# Reynolds and Rayleigh numbers.  Leave unset if you do not know what it is.\n")
            self.write("#\n")
            self.write("# Enter either 'LX', 'LY', 'LZ', or a number.\n")
            self.write("ReferenceLength ", cnfg_main["ReferenceLength"])
            self.write("\n\n")

            # Item 4
            self.write("# Kinematic viscosity\n")
            self.write("Nu ", cnfg_main["Nu"])
            self.write("\n\n")

            # Item 5
            self.write("# Specify a specific physical time (in seconds) for stopping the simulation.\n")
            self.write("StopTime  ", cnfg_main['StopTime'])
            self.write("\n\n")

            # Item 6
	    # TurbulenceModel

            # Item 7
            self.write("# Whether or not to solve for and output temperature\n")
            self.write("SolveTemperature   ", str(cnfg_main['Temperature']) + "\n")
            self.write("\n")

            # Item 8
            self.write("# To solve for and output geometry the line below must read #<SolidGeometry True>\n")            
            self.write("#<SolidGeometry " + str(cnfg_main['SolidGeometry']) + ">\n")
            self.write("\n")

            # Item 9
            # Do nothing


            # Item 10
            if self.forcing == "constPresGrad":
                self.write("# Constant forcing, requires three floating point arguments after Forcing,\n")
                self.write("# separated by spaces.\n")
                self.write("Forcing   ", cnfg_main["Forcing"])
                self.write("\n\n")

            # Item 11
            if self.forcing == "constMassFlowRate":
                self.write("# To drive a periodic flow with a constant mass flow rate, we can define \n")
                self.write("# different regions that will be given independent forcing values. The regions\n")
                self.write("# are defined by the ForcingRegionBounds tag which requires two values, the\n")
                self.write("# first for the lower bound and the second for the upper bound. The bounds are\n")
                self.write("# determined by the distance field so ensure that is loaded. The forcing will \n")
                self.write("# be adjusted by the difference between a prescribed mass flow rate and the \n")
                self.write("# simulated mass flow rate. The prescribed mass flow rate is initially set by\n")
                self.write("# the three values given to the ForcingRegionFlowRate tag. The three values\n")
                self.write("# correspond to the x-, y- and z-directions, respectively.  The number of \n")
                self.write("# ForcingRegionBounds and ForcingRegionFlowRate tag sets need to be equal to\n")
                self.write("# the number of regions given to the NumberForcingRegions tag.  The first\n")
                self.write("# ForcingRegionFlowRate tag applies to the first ForcingRegionBounds tag.\n")
                self.write("ConstantMassFlowRate  True\n")

                numforcings = int(cnfg_main["NumberForcingRegions"])
                self.write("NumberForcingRegions   ", str(numforcings) + "\n")

                heights = cnfg_main["ForcingRegionBounds"]
                mdot = cnfg_main["ForcingRegionFlowRate"]
                for i in range(numforcings):
                    self.write("ForcingRegionBounds    ", str(heights[i]) + " " + str(heights[i + 1]) + "\n")
                    self.write("ForcingRegionFlowRate  ",str(mdot[i][0]) + " " + str(mdot[i][1]) + " " + str(mdot[i][2]) + "\n")
                self.write("\n")

            # Boundary Conditions
            cnfg_bc = self.cnfg["BC"]

            # Items
            self.write("#-------------------------------------------------------------------------------\n")
            self.write("#                          Domain Boundary Conditions\n")
            self.write("#-------------------------------------------------------------------------------\n")
            self.write("\n")
            self.write("#    Face type           Boundary setting\n")
            self.write("#    ------------------  -------------------------------------\n")
            self.write("#    NoSlip              velocity 0 at boundary\n")
            self.write("#    FreeSlip            velocity unchanged at boundary\n")
            self.write("#    Inlet               BC = inlet\n")
            self.write("#    Outlet              BC = interior velocity\n")
            self.write("#    ConvectiveOutlet    See Ferziger (2001)\n")
            self.write("#    Driven              BC = inlet in normal direction\n")
            self.write("#    Periodic            BC = opposite side\n")
            self.write("#    SmoothLogLawInlet   an inlet with a constant smooth-wall logarithmic profile\n")
            self.write("#                          in the z-direction.\n")
            self.write("#    RoughLogLawInlet    an inlet with a constant rough-wall logarithmic profile\n")
            self.write("#                          in the z-direction.\n")
            self.write("#    OneSeventhsInlet    an inlet with a constant 1/7 power law profile\n")
            self.write("#                          in the z-direction.\n")
            self.write("#    DataBase            inflow database BC, usage: DataBase path_to_database\n")
            self.write("#    WallModel           Set boundary condition to Schumann wall model \n")
            self.write("\n")
            self.write("# Only one each of an Inlet and Driven are allowed.\n")
            self.write("\n")
            self.write("# Bottom is Z=0, Top   is in the Z+ direction (w component of velocity).\n")
            self.write("# South  is Y=0, North is in the Y+ direction (v component of velocity).\n")
            self.write("# West   is X=0, East  is in the X+ direction (u component of velocity).\n")
            self.write("\n")
            self.write("Face_West     ", cnfg_bc['Face_West'] + "\n")
            self.write("Face_East     ", cnfg_bc['Face_East'] + "\n")
            self.write("Face_South    ", cnfg_bc['Face_South'] + "\n")
            self.write("Face_North    ", cnfg_bc['Face_North'] + "\n")
            self.write("Face_Bottom   ", cnfg_bc['Face_Bottom'] + "\n")
            self.write("Face_Top      ", cnfg_bc['Face_Top'] + "\n")
            self.write("\n")

            # Items
            try:
	      inletVel = cnfg_bc["InletVelocities"]
	      self.write("# InletVelocities requires three floating point values separated by spaces.\n")
	      self.write("# These values apply to Inlet and LogarithmicInlet boundary conditions. With \n")
	      self.write("# the LogarithmicInlet, it is only used to determine the angle of the flow so\n")
	      self.write("# magnitude does not matter.  The magnitude is found from the FrictionVelocity.\n")
	      self.write("# Note: These are signed vectors, negative is west/south\n")
	      self.write("\n")
	      self.write("InletVelocities  ", inletVel + "\n")
	      self.write("\n")
	    except Exception as err:
	      pass

            # Items
            try:
	      perturbInflow = cnfg_bc["PerturbTurbInflow"]
	      self.write("# The perturbation cell turbulent inflow conditions (first proposed by \n")
	      self.write("# Munoz-Esparza et al. 2014) applies random perturbations to groups of \n")
	      self.write("# temperature grid cells referred to as perturbation cells. These\n")
	      self.write("# perturbation cells need their length, width and height defined with\n")
	      self.write("# PerturbBoxLength, PerturbBoxWidth and PerturbBoxHeight, respectively.\n")
	      self.write("# The number of perturbation cells from the inflow is defined by the\n")
	      self.write("# PerturbZoneDepth tag as an integer value.\n")
	      self.write("# The perturbation amplitude is defined by PerturbZoneRichardson.\n")
	      self.write("# The perturbation cell method is activated by setting a boundary to an\n")
	      self.write("# inlet and setting PerturbTurbInflow to True.\n")
	      self.write("# Note: PerturbBoxLength, PerturbBoxWidth and PerturbBoxHeight can\n")
	      self.write("# be set to the strings dx, dy and dz, respectively as a shortcut to\n")
	      self.write("# set the perturbation dimensions equal to the grid spacing.\n")
	      self.write("\n")
	      self.write("PerturbTurbInflow      ", str(perturbInflow) + "\n")
	      if cnfg_bc["PerturbTurbInflow"]:
		self.write("PerturbZoneOffset      ", cnfg_bc["PerturbZoneOffset"] + "\n")
		self.write("PerturbZoneRichardson  ", cnfg_bc["PerturbZoneRichardson"] + "\n")
		self.write("PerturbZoneDepth       ", cnfg_bc["PerturbZoneDepth"] + "\n")
		self.write("PerturbBoxLength  ", cnfg_bc["PerturbBoxLength"] + "\n")
		self.write("PerturbBoxWidth   ", cnfg_bc["PerturbBoxWidth"] + "\n")
		self.write("PerturbBoxHeight  ", cnfg_bc["PerturbBoxHeight"] + "\n")
		self.write("\n")
	    except Exception as err:
	      pass
	    
            # Initial Conditions
            cnfg_ic = self.cnfg["IC"]

            self.write("#-------------------------------------------------------------------------------\n")
            self.write("#                             Initial Conditions\n")
            self.write("#-------------------------------------------------------------------------------\n")
            self.write("\n")
            self.write("# Set the initial condition. Choose from the following:\n")
            self.write("#    File                      Load from file. Fill in Initial_* tags below.\n")
            self.write("#    TaylorGreenVortex         u =  sin(2*PI*x)*cos(2*PI*y)\n")
            self.write("#                              v = -cos(2*PI*x)*sin(2*PI*y)\n")
            self.write("#                              w = 0.0\n")
            self.write("#    TurbulentChannelX         Note: Give FrictionVelocity and ReferenceVelocity tags\n")
            self.write(
                "#                              u = u_tau ( 1/kappa * log(z+) + 5.2) + 2.0 * C/wx * cos(wx*x)*sin(wy*y)sin(wz*z)\n")
            self.write("#                              v = -C/wy * sin(wx*x)*cos(wy*y)sin(wz*z)\n")
            self.write("#                              w = -C/wz * sin(wx*x)*sin(wy*y)cos(wz*z)\n")
            self.write("#                              wx = round(0.5*LX) * Pi / LX\n")
            self.write("#                              wy = round(0.5*LY) * Pi / LY\n")
            self.write("#                              wz = round(3.0*LZ) * Pi / LZ\n")
            self.write("#                              C = 0.15 * Vref\n")
            self.write("#    TurbulentChannelXMean     u = u_tau ( 1/kappa * log(z+) + 5.2)\n")
            self.write("#                              v = 0.0\n")
            self.write("#                              w = 0.0\n")
            self.write("#    TurbulentChannelXPoise    u = 1.5 * Uinlet * ( 1 - y^2 / h^2 ) + perturbation\n")
            self.write("#                              v = perturbation\n")
            self.write("#                              w = perturbation\n")
            self.write(
                "#    TurbulentChannelZ         u =  u_tau( 1/kappa * log(x+) + 5.5) + sin(Pi*x)*cos(2*z)*sin(2*y)\n")
            self.write("#                              v =  -( 1 + cos(Pi*x) ) * sin(z) * sin(4.1*y)\n")
            self.write("#                              w = -0.5*Pi*sin(z) * sin(Pi*x) * cos(1.25*Pi*y)\n")
            self.write("#    ComplexTerrain            Note: theta determined from InletVelocity u and v components\n")
            self.write("#                              u = u_tau/kappa * log(z/z_not) * cos(theta)\n")
            self.write("#                              v = u_tau/kappa * log(z/z_not) * sin(theta)\n")
            self.write("#                              w = 0.0\n")
            self.write("#    ComplexTerrainSinPerturb  Note: theta determined from InletVelocity u and v components\n")
            self.write(
                "#                              u = u_tau/kappa * log(z/z_not) * cos(theta) + 2.0*cos(x)*sin(y)*sin(z)\n")
            self.write(
                "#                              v = u_tau/kappa * log(z/z_not) * sin(theta) - sin(x)*cos(y)*sin(z)\n")
            self.write("#                              w = -sin(x)*sin(y)*( 1.0 + cos(z) )\n")
            self.write("#    ABCFlow                   u = cos(y) + sin(z)\n")
            self.write("#                              v = sin(x) + cos(z)\n")
            self.write("#                              w = cos(x) + sin(y)\n")
            self.write("# An option not listed or no option specified at all results in an\n")
            self.write("# initial field of zero velocity.\n")
            if self.active('ic_rb_sameAsInletProfile'):
                self.write("#<Same as Inlet Profile>\n")
            self.write("InitialCondition  ", cnfg_ic['InitialCondition'] + "\n")
            self.write("\n")

            # Flow Solver Parameters
            cnfg_fs = self.cnfg['FS']

            self.write("#-------------------------------------------------------------------------------\n")
            self.write("#                             Flow Solver Parameters\n")
            self.write("#-------------------------------------------------------------------------------\n")
            self.write("\n")
            self.write("#Choose the finite difference method for spatial derivatives\n")
            self.write("#0.00 means CDS, 1.00 means FOU, values between are allowed\n")
            self.write("AdvectionScheme ", cnfg_fs['AdvectionScheme'] + '\n')
            self.write("\n\n")

            self.write("# The time derivative method\n")
            self.write("#   Euler        (first order forward Euler)\n")
            self.write("#   AdamsBash    (second order Adams-Bashforth, aka 'AB2')\n")
            self.write("TimeMethod  ", cnfg_fs['TimeMethod'] + '\n')
            self.write("\n\n")

            self.write("# Specify a computational time step. Set a to negative number if dynamic time\n")
            self.write("# step adjustment is desired.\n")
            self.write("ConstTimeStep   ", cnfg_fs["ConstTimeStep"] + '\n')
            self.write("\n\n")

            self.write("#CFL*dz/velmax is the convective dt limit\n")
            self.write("CFL      ", cnfg_fs["CFL"] + '\n')
            self.write("\n\n")

            self.write("#DTStability*dz*dz/Nu is the viscous dt limit\n")
            self.write("DTStability   ", cnfg_fs['DTStability'] + '\n')
            self.write("\n\n")

            self.write("# Choose either Iterative or Multigrid\n")
            self.write("SolverMethod	          ", cnfg_fs['SolverMethod'] + '\n')
            self.write("\n")
            self.write("# Reasonable values:\n")
            if cnfg_fs['SolverMethod'] == 'Multigrid':
                self.write("#   Cycle         V, W, or W1, recommend V\n")
                self.write("#   Smoother      Jacobi or GS\n")
                self.write("#   Loops         1-10\n")
                self.write("#   Levels        2-N, recommend about 2-3 less than the maximum\n")
                self.write("#   SmoothIters   1-10 1-10. Good choices include '2 1' and '4 2'\n")
                self.write("#   JacobiWeight  0.55 - 0.95, recommend 0.86 (from Trottenberg says 0.857)\n")
                self.write("#   SORWeight     1.0 - 1.8\n")
                self.write("MG_Cycle                V\n")
                self.write("MG_Smoother             Jacobi\n")
                self.write("MG_Loops                ", cnfg_fs['MG_Loops'] + "\n")
                self.write("MG_Levels               ", cnfg_fs['MG_Levels'] + "\n")  # default should be 50
                self.write("MG_Smoother_Iterations  ", ' '.join(cnfg_fs['MG_Smoother_Iterations']) + "\n")
                self.write("MG_Jacobi_Weight        ", cnfg_fs['MG_Jacobi_Weight'] + "\n")
                self.write("MG_SOR_Weight           1.00\n")
            else:
                self.write("#            Solver   Jacobi     GS\n")
                self.write("#            Loops    10-40      5-30\n")
                self.write("#            Weight   1.0        1.0-1.8\n")
                self.write("Iterative_Solver        Jacobi \n")
                self.write("Iterative_Loops         ", cnfg_fs['Iterative_Loops'] + "\n")
                self.write("Iterative_Weight        ", cnfg_fs['Iterative_Weight'] + "\n")
                self.write("\n")

            # Solution Output Parameters
            cnfg_so = self.cnfg['SO']
            self.write("\n")
            self.write("#-------------------------------------------------------------------------------\n")
            self.write("#                           Solution Output Parameters\n")
            self.write("#-------------------------------------------------------------------------------\n")

            op = cnfg_so['OutputPath']
            if op is not None:
                self.write("# An optional path for the results\n")
                self.write("OutputPath     ", op + "\n")
            self.write("# Optional prefix (default is 'gin3d_soln')\n")
            self.write("OutputPrefix     ", cnfg_so['OutputPrefix'] + "\n")
            self.write("\n")

            self.write("# Which 2D plane to output (XY, XZ, YZ) or Mesh for the entire 3D domain.\n")
            self.write("# Can select multiple planes.\n")
            self.write("OutputPlane Mesh\n")
            
            if self.active('so_cb_screenDump'):
	      self.write("OutputPlane ", cnfg_so['OutputPlane'] + "\n")
	      self.write("\n")
	      self.write("# If a plane is selected, choose the location of the slice along the\n")
	      self.write("# direction perpendicular to the plane as a decimal (0.0-1.0)\n")
	      self.write("OutputSlice ", cnfg_so['OutputSlice'] + "\n")
	      self.write("\n")
	      self.write("# Data output for the Matrix output: u, v, w, p, phi\n")
	      self.write("MatrixData ", cnfg_so['MatrixData'] + "\n")
	      self.write("\n")

            try:
                rd = cnfg_so['RestartDirectory']
            except KeyError:
                rd = None
            if rd is not None:
                self.write("# To restart a simulation, simply give the directory to the ouptut from the previous\n")
                self.write("# simulation. The most recent timestep present in the directory will be chosen.\n")
                self.write("# Full path starting from / is recommended. Ensure the last / is also included. This\n")
                self.write("# does not automatically place in the / if missing.\n")
                self.write("RestartDirectory   ", rd + "\n")

            self.write("\n")
            self.write("# Output formats.  Output formats can be combined by giving one line per type.\n")
            self.write("# Supported format types are:\n")
            self.write("#    Matrix          display values on the screen\n")
            self.write("#    DataFile        writes data to an ASCII file with .dat suffix\n")
            self.write("#    Binary          writes a single binary file using MPI-IO with a .bin\n")
            self.write("#                       suffix\n")
            self.write("#    PBinary         writes parallel binary files using MPI-IO with a .bin\n")
            self.write("#                       suffix\n")
            self.write("#    VTK             Visualization Toolkit Image file (.vti)\n")
            self.write("#    PVTK            Parallel VTK Image file (.pvti and .vti). Number of vti\n")
            self.write("#                       files generated is the same as the number of GPUs.\n")
            self.write("#    HDF5_Multiple   Write out heavy data in HDF5 format with XDMF files for\n")
            self.write("#                       visualization software. Produces multiple HDF5 files,\n")
            self.write("#                       one for each time step.\n")
            self.write("#\n")
            self.write("# Multiple output formats can be requested (e.g. both Matrix and VTK). The\n")
            self.write("# PVTK option cannot be combined with its sequential counterpart.  Doing so\n")
            self.write("# defaults to the VTK option. Also, VTK cannot be chosen when using more than\n")
            self.write("# eight GPUs and will default to PVTK.\n")
            [self.write("OutputFormat  ", e + "\n") for e in cnfg_so['OutputFormat']]

            self.write("\n")
            self.write("# Two methods of controlling output frequency: iterations and physical time.\n")
            self.write("# Setting an interval negative will disable that method of control.\n")
            self.write("\n")
            self.write("# Output by iteration:\n")
            self.write("# Write files every N timesteps. OutputFrequencyVTK is for VTK and PVTK file\n")
            self.write("# types. OutputFrequencyText is for DataFile, Binary, and PBinary file types.\n")
            self.write("# The output frequency can be changed part way through a simulation by\n")
            self.write("# setting OutputFrequencyVTKSwitch or OutputFrequencyTextSwitch. WTK output\n")
            self.write("# frequency will be every OutputFrequencyVTK1 timesteps until\n")
            self.write("# OutputFrequencyVTKSwitch timesteps is reached, then the frequency will be\n")
            self.write("# every OutputFrequencyVTK2 timesteps.\n")
            self.write("\n")
            try:
                self.write("OutputFrequencyVTK1          ", cnfg_so['OutputFrequencyVTK1'] + "\n")
            except KeyError:
                pass
            try:
                self.write("OutputFrequencyVTK2          ", cnfg_so['OutputFrequencyVTK2'] + "\n")
                self.write("OutputFrequencyVTKSwitch     ", cnfg_so['OutputFrequencyVTKSwitch'] + "\n")
            except KeyError:
                pass
            try:
                self.write("OutputFrequencyText1         ", cnfg_so['OutputFrequencyText1'] + "\n")
            except KeyError:
                pass
            try:
                self.write("OutputFrequencyText2         ", cnfg_so['OutputFrequencyText2'] + "\n")
                self.write("OutputFrequencyTextSwitch    ", cnfg_so['OutputFrequencyTextSwitch'] + "\n")
            except KeyError:
                pass
            self.write("\n")
            self.write("# Output by physical time:\n")
            self.write("# Writes files every time interval (in seconds) rather than every N\n")
            self.write("# See avobe for description of interval switching and corresponding file\n")
            self.write("# types.\n")
            try:
                self.write("OutputTimeIntervalVTK1          ", cnfg_so['OutputTimeIntervalVTK1'] + "\n")
            except KeyError:
                pass
            try:
                self.write("OutputTimeIntervalVTK2          ", cnfg_so['OutputTimeIntervalVTK2'] + "\n")
                self.write("OutputTimeIntervalVTKSwitch     ", cnfg_so['OutputTimeIntervalVTKSwitch'] + "\n")
            except KeyError:
                pass
            try:
                self.write("OutputTimeIntervalText1         ", cnfg_so['OutputTimeIntervalText1'] + "\n")
            except KeyError:
                pass
            try:
                self.write("OutputTimeIntervalText2         ", cnfg_so['OutputTimeIntervalText2'] + "\n")
                self.write("OutputTimeIntervalTextSwitch    ", cnfg_so['OutputTimeIntervalTextSwitch'] + "\n")
            except KeyError:
                pass
            self.write("\n")

            self.write("# Output options for OutputVariablesVTK and OutputVariablesTXT tags:\n")
            self.write("#   *   All variables\n")
            self.write("#   @   Only variables needed for restarting simulation\n")
            self.write("#   u   Velocity\n")
            self.write("#   p   Pressure\n")
            self.write("#   t   Temperature\n")
            self.write("#   q   Q-criterion\n")
            self.write("#   d   Density (PGM Solid Representation)\n")
            self.write("#   g   IB Flags(Immersed Boundary Method)\n")
            self.write("#   z   Distance field\n")
            self.write("#   n   Eddy Viscosity\n")
            self.write("#   c   Smagorinsky model coefficient (Dynamic C_s)\n")
            self.write("#   i   Inflow Perturbation\n")
            self.write("# Choose any combination of these, each character separated by spaces.  The\n")
            self.write("# time-averaged values are output in the final files, if applicable.\n")
            self.write("# For example, OutputVariablesVTK u p t will output velocity, pressure,\n")
            self.write("# instantaneous temperature.\n")
            self.write("# If no or incorrect options are given, * is the default.\n")
            self.write("\n")
            self.write("# Choose output variables for VTK\n")
            vtk = ""
            try:
                vtk = " ".join(cnfg_so['OutputVariablesVTK'])
            except:
                pass
	    self.write("OutputVariablesVTK  ", vtk)
	    
            self.write("\n")
            self.write("# Choose output variables for DataFile and/or Binary\n")
            
            txt = ""
            try:
                txt = " ".join(cnfg_so['OutputVariablesTXT'])
            except:
                pass
	    self.write("OutputVariablesTXT  ", txt)
            self.write("\n")
            self.write("\n")

            self.write("# -------------------------------------------------------------------------------\n")
            self.write("#                               Time Series Options\n")
            self.write("# -------------------------------------------------------------------------------\n")
            self.write("# Set averaging to true or false.  Set the starting time for averaging to start.\n")
            self.write("# Starting time does not matter if PrecursorSim is set to True.\n")
            self.write("TimeAvg      " + cnfg_so['TimeAvg'] + "\n")
            self.write("StartingTime " + cnfg_so['StartingTime'] + "\n")
            self.write("\n")
            self.write("# Set statistical quantity to calculate with TurbStats tags:\n")
            self.write("#   m   time averaging of primitive quantities\n")
            self.write("#   s   symmetric components of Reynolds stress tensor (shear)\n")
            self.write("#   d   diagonal components of Reynolds stress tensor (normal)\n")
            self.write("#   i   two-point covariance of u component in x- and y-directions\n")
            self.write("TurbStats " + " ".join(cnfg_so['TurbStats']) + "\n")
            self.write("\n")
            self.write("# Pass in mean velocity profile to compute turbulent fluctuations in the\n")
            self.write("# z-direction.\n")
            try:
                mp_u = cnfg_so['MeanProfile_U']
                mp_v = cnfg_so['MeanProfile_V']
                mp_w = cnfg_so['MeanProfile_W']
                if mp_u is not None and mp_v is not None and mp_w is not None:
                    self.write("MeanProfile_U   ", mp_u + "\n")
                    self.write("MeanProfile_V   ", mp_v + "\n")
                    self.write("MeanProfile_W   ", mp_w + "\n")
            except KeyError:
                pass
            try:
                self.write("MeanProfile_P   ", cnfg_so['MeanProfile_P'] + "\n")
            except (KeyError, TypeError):
                pass
            try:
                self.write("MeanProfile_Phi ", cnfg_so['MeanProfile_Phi'] + "\n")
            except (KeyError, TypeError):
                pass
            self.write("\n")
	    
	    
	    tc_file = cnfg_so['TrackCoordFile']
	    
	    if tc_file is not None:
		self.write("# -------------------------------------------------------------------------------\n")
		self.write("#                          Variable Tracking\n")
		self.write("# -------------------------------------------------------------------------------\n")
		self.write("# Track all applicable primitive variables as the simulation progresses.\n")
		self.write("#\n")
		self.write("# Give the physical time interval to output the variables e.g., giving a value\n")
		self.write("# of 20.0 will output every 20 seconds of physical time.\n")
		self.write("TrackInterval", cnfg_so['TrackInterval'] + "\n")
		self.write("\n")
           
                self.write("# Give the integer number of coordinates to track after the TrackCoordinates tag.\n")
                self.write("# Important for allocation when loading in the tracked coordinates.\n")
                self.write("# List the coordinates where tracking will occur with the TrackCoord tag in\n")
                self.write("# front, one per line, in the format:\n")
                self.write("# TrackCoord x y z\n")
                self.write("# where x, y, and z are floating point numbers. Set number of coordinates to\n")
                self.write("# track to 0 if no tracking is desired.\n")
                self.write("#\n")
                self.write("#<TrackCoordFile> ", tc_file + "\n")
                self.write("TrackCoordinates ", str(cnfg_so['TrackCoordinates']) + "\n")

                try:
                    trackCoords = cnfg_so['TrackCoordinates']
                    if trackCoords > 0:
                        trackCoord = cnfg_so['TrackCoord']
                        for tc in trackCoord:
                            self.write("TrackCoord  ", " ".join(tc) + "\n")
                except Exception as err:
                    log(err, 'e')

            self.write("#\n")

            # Temperature Equation Parameters
            if cnfg_main['Temperature']:
                try:
                    cnfg_tp = self.cnfg['TP']

                    self.write("# -------------------------------------------------------------------------------\n")
                    self.write("#                              Temperature\n")
                    self.write("# -------------------------------------------------------------------------------\n")
                    self.write("\n")
                    self.write("# Notes:\n")
                    self.write("#    Re = (ReferenceVelocity * L / Nu)\n")
                    self.write("#    Ra = (g * Beta * (Tmax - Tmin) * L^3) / (Gamma * Nu)\n")
                    self.write("#    Pr = (Nu * RhoInf) / Gamma\n")
                    self.write("#    Gr = Ra / Pr\n")
                    self.write("#Choose the finite difference method for spatial derivatives\n")
		    self.write("#0.00 means CDS, 1.00 means FOU, values between are allowed\n")
                    self.write("TemperatureAdvectionScheme ", cnfg_tp['TemperatureAdvectionScheme'] + '\n')
                    self.write("\n")
                    self.write("# Gravity magnitude\n")
                    self.write("Gravity ", " ".join(cnfg_tp['Gravity']) + "\n")
                    self.write("\n")
                    self.write("# Thermal expansion coefficient (1/T for ideal gas)\n")
                    self.write("Beta ", cnfg_tp['Beta'] + "\n")
                    self.write("\n")
                    self.write("# Reference Temperature\n")
                    self.write("Temp_Infinity ", cnfg_tp['Temp_Infinity'] + "\n")
                    self.write("\n")
                    self.write("# The turbulent prandtl number\n")
                    self.write("TurbulentPrandtl ", cnfg_tp['TurbulentPrandtl'] + "\n")
                    self.write("\n")
                    self.write("# Thermal Diffusivity  (Nu / Prandtl number)\n")
                    self.write("Gamma ", cnfg_tp['Gamma'] + "\n")
                    self.write("\n")
                    self.write("# Scalar transport for temperature\n")
                    self.write("# Note that gravity is assumed to act in the negative Z direction.\n")
                    self.write("# The Boussinesq approximation is included in the w-momentum equation.\n")
                    self.write("Temp_West   ", " ".join(cnfg_tp['Temp_West']) + "\n")
                    self.write("Temp_East   ", " ".join(cnfg_tp['Temp_East']) + "\n")
                    self.write("Temp_South  ", " ".join(cnfg_tp['Temp_South']) + "\n")
                    self.write("Temp_North  ", " ".join(cnfg_tp['Temp_North']) + "\n")
                    self.write("Temp_Bottom ", " ".join(cnfg_tp['Temp_Bottom']) + "\n")
                    self.write("Temp_Top    ", " ".join(cnfg_tp['Temp_Top']) + "\n")
                    self.write("\n")
                    self.write("# Source term for temperature\n")
                    self.write("Source_PHI ", cnfg_tp['Source_PHI'] + "\n")
                    self.write("\n")

                except Exception as err:
                    log(err, 'e')

            # Turbulence Model Parameters
            cnfg_tm = self.cnfg['TM']
            
            self.write("# -------------------------------------------------------------------------------\n")
	    self.write("#                            Turbulence Parameters\n")
	    self.write("# -------------------------------------------------------------------------------\n")
	    self.write("\n")
	    self.write("# The filter width for LES models is the numerical grid.\n")
	    self.write("# Turbulence model:\n")
	    self.write("#   Laminar        Initialize simulation as laminar\n")
	    self.write("#   None           No turbulence model but initialized as turbulent flow\n")
	    self.write("#   Smagorinsky    The basic eddy viscosity model (J. Smagorinsky, 1963)\n")
	    self.write("#   LagDynamic     The Lagrangian dynamic model Ref:(Meneveau, Lund, Cabot; 1996)\n")
	    self.write("#   SmagRANS       A hybrid RANS/LES approach with the original Smagorinsky\n")
	    self.write("#                    (LES) and Prandtl mixing length (RANS) models\n")
	    self.write("#   LagDynRANS     A hybrid RANS/LES approach with the Lagrangian dynamic\n")
	    self.write("#                    (LES) and Prandtl mixing length (RANS)  models\n")
	    self.write("#   RSCSmag        Reynolds-stress-constrained Smagorinsky model.\n")
	    self.write("#                  See S. Chen et al., JFM, 2012 for RSC concept.\n")
	    self.write("#   RSCLagDyn      Reynolds-stress-constrained Lagrangian dynamic model.\n")
	    self.write("#                  See S. Chen et al., JFM, 2012 for RSC concept.\n")
	    self.write("#   RSCSmagIB      The RSCSmag option but with immersed boundary reconstruction.\n")
	    self.write("#   RSCLagDynIB    The RSCLagDyn option but with immersed boundary reconstruction.\n")
	    self.write("#   LinearForcing  Linear forcing for DNS of isotropic turbulence.\n")
	    self.write("#                  See T.S. Lundgren, CTR Briefs, 2003.\n")
	    self.write("#   MasonThomson   Mason Thomson model used with original Smagorinsky model.\n")
	    self.write("#                  See Mason, Thomson, JFM, 1992\n")
	    self.write("\n")
	    self.write("TurbulenceModel  ", cnfg_tm['TurbulenceModel'] + "\n")
	    self.write("\n")
	    
            if cnfg_main['Turbulence']:
                try:
                    self.write("# The default initialization is a quiescent flow field.  Certain combinations\n")
                    self.write("# of boundary conditions and flow regimes have preset initial conditions in\n")
                    self.write("# the code, see Domain Boundary Conditions section for more details.  Care\n")
                    self.write("# should be taken when using the LagDynamic model as a quiescent field will\n")
                    self.write("# result in Lij and Mij both being zero.  Cs in with this model is essentially\n")
                    self.write("# (LijMij) / (MijMij).\n")
                    self.write("# Cs for Smagorinsky model.  Typically 0.01 - 0.25.\n")
                    self.write("TurbCS ", cnfg_tm['TurbCS'] + "\n")
                    self.write("\n")
                    self.write("# For RANS/LES transition height, either set the tag\n")
                    self.write("# TransitionHeight to a value above 0.00001 or set the\n")
                    self.write("# dimensionless parameter Zeta,\n")
                    self.write("#   Zeta = h / (2*Delta),\n")
                    self.write("# where h is height normal to surface and Delta is filter width.\n")
                    self.write("# TransitionHeight overrides Zeta if both are set.\n")
                    self.write("TransitionHeight  ", cnfg_tm['TransitionHeight'] + "\n")
                    self.write("\n")

                    try:
                        df = cnfg_tm['DistanceField']
                        self.write("# File for loading distance field\n")
                        self.write("DistanceField  ", df + "\n")
                    except:
                        pass

                except Exception as err:
                    log(err, 'e')

            # Geometry Parameters
            if cnfg_main['SolidGeometry']:
                try:
                    cnfg_ge = self.cnfg['GE']

                    self.write("# -------------------------------------------------------------------------------\n")
                    self.write("#                            Geometry Parameters\n")
                    self.write("# -------------------------------------------------------------------------------\n")
                    self.write("\n")
                    
                    try:
		      ibr = cnfg_ge['IBReconstruction']
		      self.write("IBReconstruction  ", ibr + "\n")
		      if ibr == 'RoughLogLaw':
			roughness = cnfg_ge['Roughness']
			self.write("Roughness  ", roughness + "\n")
		    except:
		      pass
		    
                    ib_file = cnfg_ge['IBFiles']
                    if ib_file is not None:
                        self.write("# The following are only used for the Immersed Boundary Method.\n")
                        self.write("\n")
                        self.write("# The files to load the immersed boundary information from. See README for or\n")
                        self.write("# samples under ibdata directory for format. There will be a total of eight files\n")
                        self.write("# to load, four immersed boundary node files and four obstacles flag files.\n")
                        self.write("# The four files of each type correspond to the following locations:\n")
                        self.write("#     Scalar  (Cell center)\n")
                        self.write("#     U Face  (Location of u velocity component)\n")
                        self.write("#     V Face  (Location of v velocity component)\n")
                        self.write("#     W Face  (Location of w velocity component)\n")
                        self.write("# The number of nodes loaded for each location is set above with the Obstruction\n")
                        self.write("# tag.\n")
                        self.write("#\n")
                        self.write("# With the option of IBDF, a distance field file must be specified that is in\n")
                        self.write("# the same format as as the file used for restarting a simulation.\n")
                        self.write("\n")
                        self.write("#<IBFiles> ", ib_file + "\n")
                        self.write("# Files for loading immersed boundary nodes\n")
                        self.write("IBNodeFile_S ", cnfg_ge['IBNodeFile_S'] + "\n")
                        self.write("IBNodeFile_U ", cnfg_ge['IBNodeFile_U'] + "\n")
                        self.write("IBNodeFile_V ", cnfg_ge['IBNodeFile_V'] + "\n")
                        self.write("IBNodeFile_W ", cnfg_ge['IBNodeFile_W'] + "\n")
                        self.write("\n")
                        self.write("# Files for loading obstacle flags\n")
                        self.write("IBFlagFile_S ", cnfg_ge['IBFlagFile_S'] + "\n")
                        self.write("IBFlagFile_U ", cnfg_ge['IBFlagFile_U'] + "\n")
                        self.write("IBFlagFile_V ", cnfg_ge['IBFlagFile_V'] + "\n")
                        self.write("IBFlagFile_W ", cnfg_ge['IBFlagFile_W'] + "\n")
                        self.write("\n")

                        if not cnfg_main['Turbulence']:
                            try:
                                df = cnfg_ge['DistanceField']
                                self.write("# File for loading distance field\n")
                                self.write("DistanceField  ", df + "\n")
                            except:
                                pass


                except Exception as err:
                    log(err, 'e')

            self.fp.close()
        except IndexError as err:
            log(err, 'e')
            # display(err, 'e',)
