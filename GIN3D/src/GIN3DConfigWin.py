import sys, os
from utilities import *
from handler import Handler

try:
    import gi

    gi.require_version('Gtk', '3.0')
    from gi.repository import Gtk, Gdk
except Exception as e:
    log(e, 'e')
    sys.exit(1)


class GIN3DConfigWin:
    GLADEFILE = 'GIN3D_ConfigGUI.glade'

    def __init__(self, cfg=None):

        if cfg is not None:
            try:
                self.cnfg = self.parse_file(cfg)
            except FileNotFoundError as err:
                log(err, 'e')
                return

        self.builder = Gtk.Builder()
        self.builder.add_from_file(self.GLADEFILE)
        self.builder.connect_signals(Handler(self.builder))
        win = self.builder.get_object('MainWindow')
	
        if cfg is not None:
            self.apply_cfg()
            win.set_title(win.get_title() + ' - ' + cfg[cfg.rfind('/')+1:])

        win.show()

        Gtk.main()

    def parse_file(self, filename):

        if filename is None:
            return {}

        line_number = 0
        cnfg = {}
        with open(filename, "r") as fp:
            for line in fp:
                line_number += 1
                l = line.strip().replace('\t', '')
                if l.startswith('#<SolidGeometry True>'):
		  cnfg['SolidGeometry'] = 'True'
                elif l.startswith('#<Same as Inlet Profile>'):
                    cnfg['IC_InletProfile'] = 'True'
                elif l.startswith('#<TrackCoordFile>'):
                    cnfg['TrackCoordFile'] = l.replace('#<TrackCoordFile>','').strip()
                elif l.startswith('#<IBFiles>'):
                    cnfg['IBFiles'] = l.replace('#<IBFiles>','').strip()
                elif not (l.startswith('#') or '' == l):
                    split = l.split(' ', 1)
                    key = split[0]
                    val = ''
                    if len(split) == 2:
                        val = split[1]

                    if key == 'ForcingRegionBounds':
                        try:
                            if key in cnfg:
                                cnfg[key] = cnfg[key] + ' ' + val.strip().split(' ')[1]
                            else:
                                cnfg[key] = val.strip().split(' ')[1]
                        except IndexError:
                            pass
                    elif key in ['ForcingRegionFlowRate']:
                        if key in cnfg:
                            tmp = cnfg[key]
                            tmp.append(val.strip().split())
                            cnfg[key] = tmp
                        else:
                            cnfg[key] = [val.strip().split()]
                    elif key in ['OutputFormat']:
                        if key in cnfg:
                            tmp = cnfg[key]
                            tmp.append(val.strip())
                            cnfg[key] = tmp
                        else:
                            cnfg[key] = [val.strip()]
                    elif key == 'OutputPlane' and val == 'Mesh':
                        continue
                    else:
                        cnfg[key] = val
        return cnfg

    def cfgVal(self, key):
        try:
            val = self.cnfg[key]
            if isinstance(val, str):
                val_list = val.strip().split(' ')
                return val_list[0] if len(val_list) == 1 else val_list
            return val
        except KeyError as err:
            # log("'" + key + "' not found.", 'd')
            return None

    def apply_cfg(self):

        # -------- Main parameters --------
        # Domain Size
        try:
            domain = self.cfgVal('Dimensions')
            if domain is not None and len(domain) == 3:
                self.setVal('domain_lx', domain[0])
                self.setVal('domain_ly', domain[1])
                self.setVal('domain_lz', domain[2])
        except Exception as err:
            log(err, 'e')

        # Mesh Size
        try:
            dimensions = self.cfgVal('Grid')
            if dimensions is not None and len(dimensions) == 3:
                self.setVal('mesh_nx', dimensions[0])
                self.setVal('mesh_ny', dimensions[1])
                self.setVal('mesh_nz', dimensions[2])
        except Exception as err:
            log(err, 'e')

        # Reference Length
        try:
            reference_length = self.cfgVal('ReferenceLength')
            if reference_length is not None:
                self.setVal('refLength', reference_length)
        except Exception as err:
            log(err, 'e')

        # Laminar/Kinematic Viscosity
        try:
            nu = self.cfgVal('Nu')
            if nu is not None:
                self.setVal('laminar_viscosity', nu)
        except Exception as err:
            log(err, 'e')

        # Simulation End Time
        try:
            stop_time = self.cfgVal('StopTime')
            if stop_time is not None:
                self.setVal('sim_endTime', stop_time)
        except Exception as err:
            log(err, 'e')

        # Turbulence
        try:
            turbulence = self.cfgVal('TurbulenceModel')
            if turbulence is not None and turbulence != 'Laminar':
                turbulence = True
            else:
                turbulence = False
        except Exception as err:
            log(err, 'e')
            turbulence = False

        self.setVal('simParam_switch_turbulence', turbulence)

        # Temperature
        try:
            temperature = self.cfgVal('SolveTemperature')
            if temperature is not None:
                temperature = eval(temperature)
            else:
                temperature = False
        except Exception as err:
            log(err, 'e')
            temperature = False
        self.setVal('simParam_switch_temperature', temperature)

        # Solid Geometry
        try:
            solid_geometry = self.cfgVal('SolidGeometry')
            if solid_geometry is not None:
                solid_geometry = eval(solid_geometry)
            else:
                solid_geometry = False
        except Exception as err:
            log(err, 'e')
            solid_geometry = False
        self.setVal('simParam_switch_geometry', solid_geometry)

        # Forcing
        try:
            forcing = self.cfgVal('Forcing')
            if forcing is not None:
                self.setVal('simParam_rb_forcing_constPresGrad', None)
                self.setVal('simParam_txtbox_fx', forcing[0])
                self.setVal('simParam_txtbox_fy', forcing[1])
                self.setVal('simParam_txtbox_fz', forcing[2])
        except Exception as err:
            log(err, 'e')
            forcing = None

        if forcing is None:
            try:
                constMassFlowRate = self.cfgVal('ConstantMassFlowRate')
                if constMassFlowRate is not None and eval(constMassFlowRate):
                    self.setVal('simParam_rb_forcing_constMassFlowRate', None)
                    numforcings = self.cfgVal('NumberForcingRegions')
                    heights = self.cfgVal('ForcingRegionBounds')
                    mdot = self.cfgVal('ForcingRegionFlowRate')

                    if numforcings is not None and heights is not None and mdot is not None:
                        self.setVal('no_of_forcing_regions', str(numforcings))
                        for i in range(eval(numforcings)):
                            try:
                                h = first_child("simParam_forcing_cmfr_box", "height_", i, self.builder)
                                if h is not None:
                                    self.setVal(h, heights[i])
                                    x = first_child("simParam_forcing_cmfr_box", "frfr_x_", i, self.builder)
                                    if x is not None:
                                        self.setVal(x, mdot[i][0])

                                    y = first_child("simParam_forcing_cmfr_box", "frfr_y_", i, self.builder)
                                    if y is not None:
                                        self.setVal(y, mdot[i][1])

                                    z = first_child("simParam_forcing_cmfr_box", "frfr_z_", i, self.builder)
                                    if z is not None:
                                        self.setVal(z, mdot[i][2])
                            except Exception as err:
                                log(err, 'e')

            except Exception as err:
                log(err, 'e')

        # -------- Boundary Conditions ----
        # Face Values
        # Caution: the order of the values matter
        inlet_values = ['RoughLogLawInlet', 'SmoothLogLawInlet', 'OneSeventhsInlet', 'Parabolic', 'DNS2000']
        variable_profile = None

        face_west = self.cfgVal('Face_West')
        if face_west is None:
            face_west = self.builder.get_object('west_face').get_active_id()

        if face_west in inlet_values:
            variable_profile = face_west
            face_west = 'Inlet'
        self.setVal('west_face', face_west)

        face_east = self.cfgVal('Face_East')
        if face_east is None:
            face_east = self.builder.get_object('east_face').get_active_id()
        self.setVal('east_face', face_east)

        face_south = self.cfgVal('Face_South')
        if face_south is None:
            face_south = self.builder.get_object('south_face').get_active_id()

        if face_south in inlet_values:
            variable_profile = face_south
            face_south = 'Inlet'
        self.setVal('south_face', face_south)

        face_north = self.cfgVal('Face_North')
        if face_north is None:
            face_north = self.builder.get_object('north_face').get_active_id()
        self.setVal('north_face', face_north)

        face_bottom = self.cfgVal('Face_Bottom')
        if face_bottom is None:
            face_bottom = self.builder.get_object('bottom_face').get_active_id()
        self.setVal('bottom_face', face_bottom)

        face_top = self.cfgVal('Face_Top')
        if face_top is None:
            face_top = self.builder.get_object('top_face').get_active_id()

        self.setVal('top_face', face_top)

        # Inlet Profile Parameters
        # check for inlet
        if 'Inlet' in [face_west, face_south]:
            try:
                inlet_velocities = self.cfgVal('InletVelocities')
                if inlet_velocities is not None:
                    self.setVal('inletProfile_u', inlet_velocities[0])
                    self.setVal('inletProfile_v', inlet_velocities[1])
                    self.setVal('inletProfile_w', inlet_velocities[2])
            except Exception as err:
                log(err, 'e')

        # Variable Profile
        if variable_profile is not None:
            try:
                self.setVal('inletProfile_switch_variable', True)
                # Caution: the order of the values matter
                rb_inletProfile_id = ['roughWall', 'smoothWall', 'powerLaw', 'parabolic', 'dns']
                self.setVal('rb_inletProfile_' + rb_inletProfile_id[inlet_values.index(variable_profile)], None)
            except Exception as err:
                log(err, 'e')

        # Turbulent Inlet for LES
        try:
            perturb_inflow = self.cfgVal('PerturbTurbInflow')
            if perturb_inflow is not None and eval(perturb_inflow):
                self.setVal('inletProfile_switch_turbulent', True)
                self.setVal('inletProfile_offset', self.cfgVal('PerturbZoneOffset'))
                self.setVal('inletProfile_richardson', self.cfgVal('PerturbZoneRichardson'))
                self.setVal('inletProfile_depth', self.cfgVal('PerturbZoneDepth'))
                self.setVal('inletProfile_length', self.cfgVal('PerturbBoxLength'))
                self.setVal('inletProfile_width', self.cfgVal('PerturbBoxWidth'))
                self.setVal('inletProfile_height', self.cfgVal('PerturbBoxHeight'))
        except Exception as err:
            log(err, 'e')

        # --------- Initial Conditions ----
        try:
            # check if same as inlet profile was selected
            ic_inletProfile = self.cfgVal('IC_InletProfile')
            if ic_inletProfile is not None:
                ic_inletProfile = eval(ic_inletProfile)
            else:
                ic_inletProfile = False

            if ic_inletProfile:
                self.setVal('ic_rb_sameAsInletProfile', None)
            else:
                initial_condition = self.cfgVal('InitialCondition')
                if initial_condition is not None:
                    if isinstance(initial_condition, str):
                        if initial_condition == 'ZeroField':
                            self.setVal('ic_rb_zeroField', None)
                        # Special Cases
                        else:
                            self.setVal('ic_rb_specialCases', None)
                            if initial_condition == 'TurbulentChannelX':
                                self.setVal('ic_list_specialCases', '1')
                            elif initial_condition == 'ComplexTerrainSinPerturb':
                                self.setVal('ic_list_specialCases', '2')
                            elif initial_condition == 'ABCFlow':
                                self.setVal('ic_list_specialCases', '3')
                    elif isinstance(initial_condition, list):
                        # Uniform
                        try:
                            if initial_condition[0] == 'Uniform':
                                self.setVal('ic_rb_uniform', None)
                                self.setVal('ic_uniform_u', initial_condition[1])
                                self.setVal('ic_uniform_v', initial_condition[2])
                                self.setVal('ic_uniform_w', initial_condition[3])
                        except IndexError as err:
                            log(err, 'e')
        except Exception as err:
            log(err, 'e')

        # -------- Flow Solver Parameters ----
        try:
            # Advection Scheme
            # Momentum
            try:
                advection_scheme_momentum = self.cfgVal('AdvectionScheme')
                if advection_scheme_momentum is not None:
                    if advection_scheme_momentum == '0.00':
                        self.setVal('fsp_rb_momentum_cds', None)
                    elif advection_scheme_momentum == '1.00':
                        self.setVal('fsp_rb_momentum_fou', None)
                    elif advection_scheme_momentum == '-1.00':
                        self.setVal('fsp_rb_momentum_quick', None)
                    elif advection_scheme_momentum == '':
                        self.setVal('fsp_rb_momentum_kappa', None)
                    else:
                        self.setVal('fsp_rb_momentum_hybrid', None)
                        self.setVal('fsp_text_momentum_hybrid', advection_scheme_momentum)
            except Exception as err:
                log(err, 'e')

            # Time Marching
            try:
                self.setVal('fsp_timeMarching', self.cfgVal('TimeMethod'))
            except Exception as err:
                log(err, 'e')

            # Time Step Size
            try:
                # not implemented
                const_time_step = self.cfgVal('ConstTimeStep')
                if const_time_step is not None:
                    if const_time_step == '-0.1':
                        self.setVal('fsp_timeStepSize', 'Variable')
                        self.setVal('fsp_cfl', self.cfgVal('CFL'))
                        self.setVal('fsp_viscousLimit', self.cfgVal('DTStability'))
                    else:
                        self.setVal('fsp_timeStepSize', 'Constant')
                        self.setVal('fsp_dt', self.cfgVal('ConstTimeStep'))
            except Exception as err:
                log(err, 'e')

            # Poisson Solver
            try:
                poisson_solver = self.cfgVal('SolverMethod')
                self.setVal('fsp_poissonSolver', poisson_solver)
                if poisson_solver == 'Multigrid':
                    self.setVal('fsp_text_poissonSolver_outerLoops', self.cfgVal('MG_Loops'))
                    iterations = self.cfgVal('MG_Smoother_Iterations')
                    if len(iterations) == 2:
                        self.setVal('fsp_text_poissonSolver_innerLoops_down', iterations[0])
                        self.setVal('fsp_text_poissonSolver_innerLoops_up', iterations[1])
                    self.setVal('fsp_text_poissonSolver_mg_jacobiWeight', self.cfgVal('MG_Jacobi_Weight'))
                    self.setVal('fsp_text_poissonSolver_mgLevels', self.cfgVal('MG_Levels'))
                elif poisson_solver == 'Iterative':
                    self.setVal('fsp_text_poissonSolver_iterations', self.cfgVal('Iterative_Loops'))
                    self.setVal('fsp_text_poissonSolver_pj_jacobiWeight', self.cfgVal('Iterative_Weight'))
            except Exception as err:
                log(err, 'e')

        except Exception as err:
            log(err, 'e')

        # -------- Solution Output Parameters ----
        try:
            # Output path
            self.setVal('so_folder_outputPath', self.cfgVal('OutputPath'))
            # Output prefix
            self.setVal('so_text_outputPrefix', self.cfgVal('OutputPrefix'))

            # Restart Solution
            rd = self.cfgVal('RestartDirectory')
            if rd is not None:
                self.setVal('so_switch_restartSolution', True)
                self.setVal('so_folder_restartDirectory', rd)

            outputFormat = self.cfgVal('OutputFormat')
            # Formatted Output
            # File Formats
            if outputFormat is not None:
                if 'HDF5_Multiple' in outputFormat:
                    self.setVal('so_cb_fileFormat_hdf5', True)
                if 'VTK' in outputFormat:
                    self.setVal('so_cb_fileFormat_vti', True)

            # Time Step
            ts_vtk1 = self.cfgVal('OutputFrequencyVTK1')
            ts_vtk2 = self.cfgVal('OutputFrequencyVTK2')
            ts_vtks = self.cfgVal('OutputFrequencyVTKSwitch')

            if ts_vtk1 is not None:
                self.setVal('so_cb_formatOutput_timeStep', True)
                self.setVal('so_text_formatOutput_timeStep_interval1', ts_vtk1)
            if ts_vtk2 is not None:
                self.setVal('so_cb_formatOutput_timeStep_switch', True)
                self.setVal('so_text_formatOutput_timeStep_interval2', ts_vtk2)
            if ts_vtks is not None:
                self.setVal('so_cb_formatOutput_timeStep_switch', True)
                self.setVal('so_text_formatOutput_timeStep_switch', ts_vtks)

            # Physical Time
            pt_vtk1 = self.cfgVal('OutputTimeIntervalVTK1')
            pt_vtk2 = self.cfgVal('OutputTimeIntervalVTK2')
            pt_vtks = self.cfgVal('OutputTimeIntervalVTKSwitch')

            if pt_vtk1 is not None:
                self.setVal('so_cb_formatOutput_physicalTime', True)
                self.setVal('so_text_formatOutput_physicalTime_interval1', pt_vtk1)
            if pt_vtk2 is not None:
                self.setVal('so_cb_formatOutput_physicalTime_switch', True)
                self.setVal('so_text_formatOutput_physicalTime_interval2', pt_vtk2)
            if pt_vtks is not None:
                self.setVal('so_cb_formatOutput_physicalTime_switch', True)
                self.setVal('so_text_formatOutput_physicalTime_switch', pt_vtks)

            # Output Variables
            output_variables_vtk = self.cfgVal('OutputVariablesVTK')
            if output_variables_vtk is not None:
                self.setVal('so_cb_outputVTK_u', 'u' in output_variables_vtk)
                self.setVal('so_cb_outputVTK_p', 'p' in output_variables_vtk)
                self.setVal('so_cb_outputVTK_t', 't' in output_variables_vtk)
                self.setVal('so_cb_outputVTK_n', 'n' in output_variables_vtk)
                self.setVal('so_cb_outputVTK_z', 'z' in output_variables_vtk)
                self.setVal('so_cb_outputVTK_q', 'q' in output_variables_vtk)
                self.setVal('so_cb_outputVTK_g', 'g' in output_variables_vtk)
                self.setVal('so_cb_outputVTK_c', 'c' in output_variables_vtk)
                self.setVal('so_cb_outputVTK_i', 'i' in output_variables_vtk)

            # Raw Data
            # File Formats
            if outputFormat is not None:
                if 'Binary' in outputFormat:
                    self.setVal('so_cb_rawData_binary', True)
                elif 'PBinary' in outputFormat:
                    self.setVal('so_cb_rawData_parallelBinary', True)
                if 'DataFile' in outputFormat:
                    self.setVal('so_cb_rawData_ASCII', True)
                elif 'PDataFile' in outputFormat:
                    self.setVal('so_cb_rawData_parallelASCII', True)

            # Time Step
            ts_txt1 = self.cfgVal('OutputFrequencyText1')
            ts_txt2 = self.cfgVal('OutputFrequencyText2')
            ts_txts = self.cfgVal('OutputFrequencyTextSwitch')

            if ts_txt1 is not None:
                self.setVal('so_cb_rawData_timeStep', True)
                self.setVal('so_text_rawData_timeStep_interval1', ts_txt1)
            if ts_txt2 is not None:
                self.setVal('so_cb_rawData_timeStep_switch', True)
                self.setVal('so_text_rawData_timeStep_interval2', ts_txt2)
            if ts_txts is not None:
                self.setVal('so_cb_rawData_timeStep_switch', True)
                self.setVal('so_text_rawData_timeStep_switch', ts_txts)

            # Physical Time
            pt_txt1 = self.cfgVal('OutputTimeIntervalText1')
            pt_txt2 = self.cfgVal('OutputTimeIntervalText2')
            pt_txts = self.cfgVal('OutputTimeIntervalTextSwitch')

            if pt_txt1 is not None:
                self.setVal('so_cb_rawData_physicalTime', True)
                self.setVal('so_text_rawData_physicalTime_interval1', pt_txt1)
            if pt_txt2 is not None:
                self.setVal('so_cb_rawData_physicalTime_switch', True)
                self.setVal('so_text_rawData_physicalTime_interval2', pt_txt2)
            if pt_txts is not None:
                self.setVal('so_cb_rawData_physicalTime_switch', True)
                self.setVal('so_text_rawData_physicalTime_switch', pt_txts)

            # Output Variables
            output_variables_txt = self.cfgVal('OutputVariablesTXT')
            if output_variables_txt is not None:
                self.setVal('so_cb_outputTXT_u', 'u' in output_variables_txt)
                self.setVal('so_cb_outputTXT_p', 'p' in output_variables_txt)
                self.setVal('so_cb_outputTXT_t', 't' in output_variables_txt)
                self.setVal('so_cb_outputTXT_n', 'n' in output_variables_txt)
                self.setVal('so_cb_outputTXT_z', 'z' in output_variables_txt)
                self.setVal('so_cb_outputTXT_q', 'q' in output_variables_txt)
                self.setVal('so_cb_outputTXT_g', 'g' in output_variables_txt)
                self.setVal('so_cb_outputTXT_c', 'c' in output_variables_txt)
                self.setVal('so_cb_outputTXT_i', 'i' in output_variables_txt)

            # Screen Dump
            # File Formats
            if outputFormat is not None:
                if 'Matrix' in outputFormat:
                    self.setVal('so_cb_screenDump', True)

            output_plane = self.cfgVal('OutputPlane')
            if output_plane is not None:
                if output_plane.upper() == 'XY':
                    self.setVal('so_rb_screenDump_xy', None)
                    self.setVal('so_text_screenDump_xy', self.cfgVal('OutputSlice'))
                elif output_plane.upper() == 'YZ':
                    self.setVal('so_rb_screenDump_yz', None)
                    self.setVal('so_text_screenDump_yz', self.cfgVal('OutputSlice'))
                elif output_plane.upper() == 'XZ':
                    self.setVal('so_rb_screenDump_xz', None)
                    self.setVal('so_text_screenDump_xz', self.cfgVal('OutputSlice'))

            matrix_data = self.cfgVal('MatrixData')
            if matrix_data is not None:
                if matrix_data.lower() == 'u':
                    self.setVal('so_rb_screenDump_u', None)
                elif matrix_data.lower() == 'v':
                    self.setVal('so_rb_screenDump_v', None)
                elif matrix_data.lower() == 'w':
                    self.setVal('so_rb_screenDump_w', None)
                elif matrix_data.lower() == 'p':
                    self.setVal('so_rb_screenDump_p', None)
                elif matrix_data.lower() == 'phi':
                    self.setVal('so_rb_screenDump_temperature', None)

            # Sampling
            # Read the coordinate file (*.csv)
            self.setVal('so_file_coordFile', self.cfgVal('TrackCoordFile'))

            # Sampling Frequency
            # Physical Time - Interval
            self.setVal('so_text_physicalTime_interval', self.cfgVal('TrackInterval'))

            # Collect Statistics
            try:
                self.setVal('so_switch_collectStatistics', eval(self.cfgVal('TimeAvg')))
            except:
                pass

            # Sampling Time Begin
            try:
                self.setVal('so_text_samplingTimeBegin', self.cfgVal('StartingTime'))
            except:
                pass

            # Order Statistics
            turb_stats = self.cfgVal('TurbStats')
            if turb_stats is not None:
                if isinstance(turb_stats, str):
                    turb_stats = [turb_stats]

                lower_turb_stats = [o.lower() for o in turb_stats]
                if 'm' in lower_turb_stats:
                    self.setVal('so_rb_firstOrderStatistics', None)
                else:
                    self.setVal('so_rb_higherOrderStatistics', None)
                    self.setVal('so_cb_stats_d', 'd' in lower_turb_stats)
                    self.setVal('so_cb_stats_s', 's' in lower_turb_stats)
                    self.setVal('so_cb_stats_i', 'i' in lower_turb_stats)

                # Spatial Averaging
                self.setVal('so_switch_spatialAveraging', turb_stats[0] in ['M', 'S', 'D', 'I'])

            # Mean quantity files for fluctuation calculations
            mean_profile_u = self.cfgVal('MeanProfile_U')
            mean_profile_v = self.cfgVal('MeanProfile_V')
            mean_profile_w = self.cfgVal('MeanProfile_W')
            mean_profile_p = self.cfgVal('MeanProfile_P')
            mean_profile_phi = self.cfgVal('MeanProfile_Phi')

            try:
                if mean_profile_u is not None or mean_profile_v is not None or mean_profile_w is not None:
                    self.setVal('so_file_meanprofile_u', mean_profile_u)
                    self.setVal('so_file_meanprofile_v', mean_profile_v)
                    self.setVal('so_file_meanprofile_w', mean_profile_w)
                    self.setVal('so_cb_mean_velocity', True)
            except:
                pass

            if mean_profile_p is not None:
                self.setVal('so_cb_mean_pressure', True)
                self.setVal('so_file_meanprofile_p', mean_profile_p)
            if mean_profile_phi is not None:
                self.setVal('so_cb_mean_temperature', True)
                self.setVal('so_file_meanprofile_phi', mean_profile_p)

        except Exception as err:
            log(err, 'e')

        # -------- Temperature Equation Parameters ----
        try:
            if temperature:
                face_list = ['West', 'East', 'South', 'North', 'Bottom', 'Top']
                for face in face_list:
                    val = self.cfgVal('Temp_' + face)
                    try:
                        if isinstance(val, list):
                            self.setVal('temp_' + face.lower() + '_face', val[0])
                            self.setVal('temp_text_' + face.lower() + 'Face', val[1])
                        else:
                            self.setVal('temp_' + face.lower() + '_face', val)
                    except Exception as er:
                        log(er, 'e')

                # Advection Scheme Temperature
                try:
                    advection_scheme_temperature = self.cfgVal('TemperatureAdvectionScheme')
                    if advection_scheme_temperature is not None:
                        if advection_scheme_temperature == '0.00':
                            self.setVal('fsp_rb_temperature_cds', None)
                        elif advection_scheme_temperature == '1.00':
                            self.setVal('fsp_rb_temperature_fou', None)
                        elif advection_scheme_temperature == '-1.00':
                            self.setVal('fsp_rb_temperature_quick', None)
                        elif advection_scheme_temperature == '':
                            self.setVal('fsp_rb_temperature_kappa', None)
                        else:
                            self.setVal('fsp_rb_temperature_hybrid', None)
                            self.setVal('fsp_text_temperature_hybrid', advection_scheme_temperature)
                except Exception as err:
                    log(err, 'e')

                # Model Parameters
                # Thermal expansion coefficient
                self.setVal('temp_text_thermalExpansionCoefficient', self.cfgVal('Beta'))
                # Thermal diffusivity
                self.setVal('temp_text_thermalDiffusivity', self.cfgVal('Gamma'))

                # Prandtl number
                self.setVal('temp_text_prandtlNumber', self.cfgVal('TurbulentPrandtl'))

                # Reference temperature
                self.setVal('temp_text_refTemperature', self.cfgVal('Temp_Infinity'))

                # Gravity
                gravity = self.cfgVal('Gravity')
                if gravity is not None:
                    try:
                        self.setVal('temp_text_gravity_x', gravity[0])
                        self.setVal('temp_text_gravity_y', gravity[1])
                        self.setVal('temp_text_gravity_z', gravity[2])
                    except Exception as err:
                        log(err, 'e')

                # Source magnitude
                self.setVal('temp_text_sourceMagnitude', self.cfgVal('Source_PHI'))
        except Exception as err:
            log(err, 'e')

        # -------- Turbulence Model Parameters ----
        try:
            if turbulence:
                turbulenceModel = self.cfgVal('TurbulenceModel')
                if turbulenceModel is not None:
                    if turbulenceModel == 'Smagorinsky':
                        self.setVal('tm_rb_originalSmagorinsky', True)
                    elif turbulenceModel == 'LagDynamic':
                        self.setVal('tm_rb_lagrangianDynamicSmagorinsky', True)
                    elif turbulenceModel == 'SmagRANS':
                        self.setVal('tm_rb_prandtlSmagorinsky', True)
                    elif turbulenceModel == 'LagDynRANS':
                        self.setVal('tm_rb_prandtlLagrangianDynamicSmagorinsky', True)

                # Smagorinsly Coefficient
                self.setVal('tm_text_smagorinskyCoefficient', self.cfgVal('TurbCS'))

                # RANS-LES Blending Height
                self.setVal('tm_text_blendingHeight', self.cfgVal('TransitionHeight'))

                # Distance Field
                distanceField = self.cfgVal('DistanceField')
                if distanceField is not None:
                    if distanceField == 'TPC_X_WallNormal':
                        self.setVal('geom_list_distanceField', '1')
                        self.setVal('geom_rb_channel_x', True)
                    elif distanceField == 'TPC_Y_WallNormal':
                        self.setVal('geom_list_distanceField', '1')
                        self.setVal('geom_rb_channel_y', True)
                    elif distanceField == 'TPC_Z_WallNormal':
                        self.setVal('geom_list_distanceField', '1')
                        self.setVal('geom_rb_channel_z', True)
                    else:
                        self.setVal('geom_list_distanceField', '0')
                        self.setVal('geom_file_distanceField', distanceField)

        except Exception as err:
            log(err, 'e')

        # -------- Geometry Parameters ----
        try:
            if solid_geometry:

                # IB Reconstruction Scheme
                ibReconstruction = self.cfgVal('IBReconstruction')
                if ibReconstruction is not None:
                    if ibReconstruction == 'Linear':
                        self.setVal('geom_rb_ib_linear', True)
                    elif ibReconstruction == 'RoughLogLaw':
                        self.setVal('geom_rb_ib_roughLogLaw', True)
                        self.setVal('geom_text_roughness', self.cfgVal('Roughness'))
                    elif ibReconstruction == 'OneSevenths':
                        self.setVal('geom_rb_ib_powerLaw', True)

                if not turbulence:
                    # Distance Field
                    distanceField = self.cfgVal('DistanceField')
                    if distanceField is not None:
                        if distanceField == 'TPC_X_WallNormal':
                            self.setVal('tm_list_distanceField', '1')
                            self.setVal('tm_rb_channel_x', True)
                        elif distanceField == 'TPC_Y_WallNormal':
                            self.setVal('tm_list_distanceField', '1')
                            self.setVal('tm_rb_channel_y', True)
                        elif distanceField == 'TPC_Z_WallNormal':
                            self.setVal('tm_list_distanceField', '1')
                            self.setVal('tm_rb_channel_z', True)
                        else:
                            self.setVal('tm_list_distanceField', '0')
                            self.setVal('tm_file_distanceField', distanceField)
        except Exception as err:
            log(err, 'e')
        return

    def setVal(self, obj_id, val):

        if isinstance(obj_id, str):
            obj = self.builder.get_object(obj_id)
        else:
            obj = obj_id

        if obj is None:
            log('Object Id: ' + obj_id + ' not found.', 'e')
            return

        if isinstance(obj, Gtk.Entry):
            if val is None:
                return
            obj.set_text(val)
            #event = Gdk.Event(Gdk.EventType.FOCUS_CHANGE)
            event = Gdk.Event()
            event.send_event = True
            event.in_ = False
            obj.emit('focus-out-event', event)

        elif isinstance(obj, Gtk.Switch):
            if val is None:
                return
            obj.set_active(not val)
            #event = Gdk.Event(Gdk.EventType.BUTTON_RELEASE)
            event = Gdk.Event()
            event.send_event = True
            obj.emit('button-release-event', event)
            obj.set_active(val)

        elif isinstance(obj, Gtk.RadioButton):
            obj.set_active(True)
            obj.emit('toggled')

        elif isinstance(obj, Gtk.CheckButton):
            obj.set_active(val)
            obj.emit('toggled')

        elif isinstance(obj, Gtk.ComboBoxText):
            if val is None:
                return
            obj.set_active_id(val)
            obj.emit('changed')

        elif isinstance(obj, Gtk.FileChooserButton):
            if val is None:
                return
            obj.set_filename(val)
