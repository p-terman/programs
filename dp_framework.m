
function status = dp_framework(filename_evt)

addpath(genpath(['~/LUXcode/']));
pathbase= '/home/paul/LUXdata'

filename_prefix=filename_evt(1:19)

data_path_evt = [pathbase filesep filename_prefix filesep];
filename_rq = strrep(filename_evt,'evt','rq');
data_path_rq = [pathbase filesep];

 lug_iqs_xml_file = '~/MATLAB/lug_iqs_new4.xml';
iq_xml_path=lug_iqs_xml_file;

% xmlfid = fopen(lug_iqs_xml_file,'r');
% lug_iqs_xml = char(fread(xmlfid,inf,'uchar')');
% fclose(xmlfid);
%************check dp_settings_new for max num pulses ******* also check
%rubiks cube

%data_processing_settings_path = which('default_data_processing_settings.xml');
data_processing_settings_path ='~/MATLAB/dp_settings_new.xml';
data_processing_xml_path =data_processing_settings_path;
%might want to edit /home/paul/Matlab/dp_settings_new.xml for max num pulses


% Don't use % status = InitializeRQFile_Default(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_settings_path,lug_iqs_xml_file);

%% Here we start with some cpp modules
%in order to get the cpp modules to work, they must first be complied. Go
%to the dp framework folder for the relevant cpp files version and type
%'make' . That should have all the modules in the bin folder to run a
%module, just use its full path a give it the same 6 parameters. if you’re
%inside the CppModules folder, binocessing_settings_path/PhotonTiming path_to_evt_file . I have added the path to .bashrc. these
%can only be done in the linux terminal.

%use instead: 

system(['/home/paul/LUXcode/CppModules/bin/InitializeRQFile_Initialize ' filename_evt ' ' data_path_evt ' ' filename_rq ' ' pathbase ' 1 ' data_processing_settings_path ' ' lug_iqs_xml_file ])
%InitializeRQFile_Initialize luxsm_20170608T1512_f000000001.evt /home/paul/LUXdata/luxsm_20170608T1512 luxsm_20170608T1512_f000000001.rq /home/paul/LUXdata/ 1  /home/paul/Matlab/dp_settings_new.xml /home/paul/Matlab/lug_iqs_new4.xml

status = PulseCalibration_BaselineZen(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_settings_path,lug_iqs_xml_file);
status = PODSummer_LUXSumPOD(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_settings_path,lug_iqs_xml_file);

%{ 
I might need to compile perusePeeks.c before
PulseFinder_TransparentRubiksCube will work. This code should do it. 
cd ~/LUXcode/Stable_Releases/DataProcessingStableReleases/v2.0/DataProcessing
cd MatlabModules/PulseFinder_TransparentRubiksCube/
mex perusePeeks.c 
%}

status = PulseFinder_TransparentRubiksCube(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_settings_path,lug_iqs_xml_file);
system(['/home/paul/LUXcode/CppModules/bin/PulseTiming_HeightTiming ' filename_evt ' ' data_path_evt ' ' filename_rq ' ' pathbase ' 5 ' data_processing_settings_path ' ' lug_iqs_xml_file ])
%PulseTiming_HeightTiming  luxsm_20180410T0819_f000000001.evt /home/paul/LUXdata/luxsm_20180410T0819_f000000001/ luxsm_20180410T0819_f000000001.rq /home/paul/LUXdata/ 5  /home/paul/Matlab/dp_settings_new.xml /home/paul/Matlab/lug_iqs_new4.xml
%PulseTiming_HeightTiming lux10_20120210T1756_f000001_eb00067.evt /media/paul/TOSHIBA/lux10_20120210T1756 lux10_20120210T1756_f000001_eb00067.rq /media/paul/TOSHIBA/ 5  /home/paul/Matlab/dp_settings_new.xml /home/paul/Matlab/lug_iqs_new4.xml
 
%not used status = PulseTiming_PerusePeeksMatlab(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
%don't use: status = PulseTiming_BasicSet(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path); % make sure the the xml has the snipet added for this module

status = PulseQuantities_MinimumSet(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);

status = PulseQuantities_PhotonCounting(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);

status = PulseClassifier_MultiDimensional(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);

status = S1S2Pairing_Naive(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);

status = Event_Classification(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
%do we even need this? This looks like it is for golden events

%don't need% PositionReconstructizon_CorrCentroid lux10_20130506T2323_f000502_eb00020.evt /home/paul/LUXdata/lux10_20130506T2323/ lux10_20130506T2323_f000502_eb00020.rq /home/paul/LUXdata 13  /home/paul/Matlab/data_processing_settings.xml /home/paul/Matlab/iqs.xml
% hitmap module% don't need

PositionReconstruction_MercuryI (filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
Corrections_PositionCorrection(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
Corrections_ApplyCorrections(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)

% 170608 tomaz's talk recommended removal of EnergyReconstruction_Naive(
%EnergyReconstruction_Naive(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
PulseQuantities_TimeSince(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path) %deals with etrains, and looks at the time since the last 'big' event

%PulseQuantities_WaterPmtRQs lux10_20130506T2323_f000502_eb00020.evt /home/paul/LUXdata/lux10_20130506T2323/ lux10_20130506T2323_f000502_eb00020.rq /home/paul/LUXdata 20  /home/paul/Matlab/dp_settings_new.xml /home/paul/Matlab/lug_iqs_new4.xml
%TriggerRQModule(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%%
status =AdditionalFileFormat_SaveMatFile(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%this is to save the .rq file as a .mat file saved /home/paul/LUXdata/matfiles//lux10_20130506T2323_f000502_eb00020.rq.mat
