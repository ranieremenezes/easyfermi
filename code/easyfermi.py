import os
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QDialog, QApplication, QFileDialog
from PyQt5.uic import loadUi
import sys
import glob
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib
import numpy as np
import astropy.io.fits as pyfits  # We need astropy version 5.2.2 or earlier
from astropy.time import Time
from fermipy.gtanalysis import GTAnalysis
from fermipy.plotting import ROIPlotter
import platform
from astroquery import fermi  # It requires astroquery version 0.4.6
from astroquery.simbad import Simbad
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
import psutil  # Version: 5.9.8
from scipy import interpolate
from gammapy.modeling.models import EBLAbsorptionNormSpectralModel  # Version 0.20.1
import emcee  # Version: 3.1.4
import corner  # Version: 2.2.2
from pathlib import Path
import yaml
from yaml import Loader

plt.rcParams.update({'font.size': 12})
libpath = Path(__file__).parent.resolve() / Path("resources/images")
EBLpath = Path(__file__).parent.resolve() / Path("resources")
Working_directory = os.getcwd()
OS_name = platform.system()

matplotlib.interactive(True)
os.environ["LANG"] = 'C'

try:
    os.environ['LC_NAME'] = 'en_IN:en'
except:
    print("No 'LC_NAME' entry was found.")
    
try:
    os.environ["LANGUAGE"] = 'en_IN:en'
except:
    print("No 'LANGUAGE' entry was found.")

try:
    os.environ['LC_MEASUREMENT'] = 'en_IN:en'
except:
    pass

try:
    os.environ['LC_PAPER'] = 'en_IN:en'
except:
    pass
    
try: 
    os.environ['LC_MONETARY'] = 'en_IN:en'
except:
    pass
    
try:
    os.environ['LC_ADDRESS'] = 'en_IN:en'
except:
    pass
    
try:
    os.environ['LC_NUMERIC'] = 'en_IN:en'
except:
    pass
    
try:
    os.environ['LC_TELEPHONE'] = 'en_IN:en'
except:
    pass
    
try:
    os.environ['LC_IDENTIFICATION'] = 'en_IN:en'
except:
    pass

try:
    os.environ['LC_TIME'] = 'en_IN:en'
except:
    pass



class Worker(QtCore.QObject):

    """
    This class is used to run the functions from the main class (i.e. Ui_mainWindow) in a different thread, 
    so avoiding the freezing of the easyfermi window when the analysis is ongoing.
    """

    starting = QtCore.pyqtSignal()
    starting_download_photons = QtCore.pyqtSignal()
    finished = QtCore.pyqtSignal()
    finished_download_photons = QtCore.pyqtSignal()
    progress = QtCore.pyqtSignal(int)
    
    
    def run_gtsetup(self):
        """Long-running task."""
        self.starting.emit()
        ui.setFermipy()
        self.progress.emit(0)  # These numbers 0, 1, 2 etc defined by emit() enter as "n" in the function reportProgress(self,n)
        ui.gta.setup()
        self.progress.emit(1)
        N_iter_adaptive_LC = ui.analysisBasics()
        self.progress.emit(2)
        calculate_Sun = ui.fit_model()
        if calculate_Sun:
            self.progress.emit(3)
            ui.Sun_path()

        self.progress.emit(4)
        ui.relocalize_the_target()
        self.progress.emit(5)
        ui.compute_TSmap()
        self.progress.emit(6)
        ui.compute_Extension()
        self.progress.emit(7)
        ui.compute_SED()
        self.progress.emit(8)
        ui.EBL_and_MCMC()
        self.progress.emit(9)
        ui.compute_LC()
        ui.plot_LCs(adaptive=False)
        self.progress.emit(10)
        for i in range(N_iter_adaptive_LC):
            ui.compute_LC_adaptive()
        ui.plot_LCs(adaptive=True)
        self.progress.emit(11)
        self.finished.emit()


class Downloads(QtCore.QObject):

    """
    This class is used to download the Fermi-LAT data using in a different thread than that 
    of the easyfermi window, then avoiding that it freezes when the downloads are ongoing.
    """
    starting = QtCore.pyqtSignal()
    finished_downloads = QtCore.pyqtSignal()
    progress = QtCore.pyqtSignal(int)
    
    def run_download_SC(self):
        """Download spacecraft file."""
        self.starting.emit()
        self.progress.emit(-1)
        ui.download_SC()
        self.finished_downloads.emit()
    
    def run_download_Photons(self):
        """Download photon files."""
        self.starting.emit()
        self.progress.emit(-2)
        ui.download_Photons()
        self.finished_downloads.emit()

    def run_download_Diffuse(self):
        """Download diffuse models."""
        self.starting.emit()
        self.progress.emit(-3)
        ui.download_Diffuse()
        self.finished_downloads.emit()



class Ui_mainWindow(QDialog):

    """Class hosting the main window of easyfermi as well as all major functions."""

    def setupUi(self, mainWindow):

        """
        In this function we setup all the easyfermi buttons, checkboxes, toolboxes etc, i.e.
        this is the function that tells easyfermi how it must look like.

        Here we also set the tooltips and associate the buttons to their respective functions.  
        """
        mainWindow.setObjectName("mainWindow")
        mainWindow.resize(1110, 595)
        mainWindow.setWindowOpacity(1.0)
        self.centralwidget = QtWidgets.QWidget(mainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.progressBar = QtWidgets.QProgressBar(self.centralwidget)
        self.progressBar.setGeometry(QtCore.QRect(10, 530, 1090, 20))
        self.progressBar.setProperty("value", 0)
        self.progressBar.setObjectName("progressBar")
        self.call_download_window = None  # The donwload window is off
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.number_of_clicks_to_download = 1  # This variable will be called only if a download is requested

        self.label_configuration_file = QtWidgets.QLabel(self.centralwidget)
        self.label_configuration_file.setGeometry(QtCore.QRect(10, 0, 121, 21))
        self.label_configuration_file.setObjectName("label_configuration_file")

        self.radioButton_Standard = QtWidgets.QRadioButton(self.centralwidget)
        self.radioButton_Standard.setEnabled(True)
        self.radioButton_Standard.setGeometry(QtCore.QRect(20, 20, 91, 23))
        self.radioButton_Standard.setChecked(True)
        self.radioButton_Standard.setAutoExclusive(True)
        self.radioButton_Standard.setObjectName("radioButton_Standard")

        ##################################################################
        ######## Group box Science
        ##################################################################

        self.groupBox_Science = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_Science.setGeometry(QtCore.QRect(300, 190, 341, 331))
        self.groupBox_Science.setObjectName("groupBox_Science")

        #LC:
        self.label_N_time_bins_LC = QtWidgets.QLabel(self.groupBox_Science)
        self.label_N_time_bins_LC.setEnabled(False)
        self.label_N_time_bins_LC.setGeometry(QtCore.QRect(95, 40, 121, 21))
        self.label_N_time_bins_LC.setObjectName("label_N_time_bins_LC")
        self.label_N_cores_LC = QtWidgets.QLabel(self.groupBox_Science)
        self.label_N_cores_LC.setEnabled(False)
        self.label_N_cores_LC.setGeometry(QtCore.QRect(252, 40, 81, 21))
        self.label_N_cores_LC.setObjectName("label_N_cores_LC")
        self.checkBox_LC = QtWidgets.QCheckBox(self.groupBox_Science)
        self.checkBox_LC.setGeometry(QtCore.QRect(10, 20, 131, 23))
        self.checkBox_LC.setObjectName("checkBox_LC")
        self.spinBox_LC_N_time_bins = QtWidgets.QSpinBox(self.groupBox_Science)
        self.spinBox_LC_N_time_bins.setEnabled(False)
        self.spinBox_LC_N_time_bins.setGeometry(QtCore.QRect(40, 40, 48, 21))
        self.spinBox_LC_N_time_bins.setRange(3,999)
        self.spinBox_LC_N_time_bins.setProperty("value", 20)
        self.spinBox_LC_N_time_bins.setObjectName("spinBox_LC_N_time_bins")
        self.spinBox_N_cores_LC = QtWidgets.QSpinBox(self.groupBox_Science)
        self.spinBox_N_cores_LC.setEnabled(False)
        self.spinBox_N_cores_LC.setGeometry(QtCore.QRect(200, 40, 48, 21))
        self.spinBox_N_cores_LC.setMinimum(1)
        self.spinBox_N_cores_LC.setProperty("value", 1)
        self.spinBox_N_cores_LC.setObjectName("spinBox_N_cores_LC")
        self.checkBox_adaptive_binning = QtWidgets.QCheckBox(self.groupBox_Science)
        self.checkBox_adaptive_binning.setGeometry(QtCore.QRect(40, 65, 131, 23))
        self.checkBox_adaptive_binning.setObjectName("checkBox_adaptive_binning")
        self.checkBox_adaptive_binning.setChecked(False)
        self.checkBox_adaptive_binning.setEnabled(False)
        self.white_box_TS_threshold = QtWidgets.QLineEdit(self.groupBox_Science)
        self.white_box_TS_threshold.setGeometry(QtCore.QRect(60, 90, 48, 21))
        self.white_box_TS_threshold.setObjectName("white_box_TS_threshold")
        self.white_box_TS_threshold.setEnabled(False)
        self.label_TS_threshold = QtWidgets.QLabel(self.groupBox_Science)
        self.label_TS_threshold.setEnabled(False)
        self.label_TS_threshold.setGeometry(QtCore.QRect(112, 90, 121, 21))
        self.label_TS_threshold.setObjectName("label_TS_threshold")
        self.spinBox_N_iter = QtWidgets.QSpinBox(self.groupBox_Science)
        self.spinBox_N_iter.setEnabled(False)
        self.spinBox_N_iter.setGeometry(QtCore.QRect(200, 90, 48, 21))
        self.spinBox_N_iter.setRange(1,10)
        self.spinBox_N_iter.setProperty("value", 1)
        self.spinBox_N_iter.setObjectName("spinBox_N_iter")
        self.label_N_iter = QtWidgets.QLabel(self.groupBox_Science)
        self.label_N_iter.setEnabled(False)
        self.label_N_iter.setGeometry(QtCore.QRect(252, 90, 121, 21))
        self.label_N_iter.setObjectName("label_N_iter")


        #SED:
        self.checkBox_SED = QtWidgets.QCheckBox(self.groupBox_Science)
        self.checkBox_SED.setGeometry(QtCore.QRect(10, 110, 131, 23))
        self.checkBox_SED.setChecked(True)
        self.checkBox_SED.setObjectName("checkBox_SED")
        self.label_N_energy_bins = QtWidgets.QLabel(self.groupBox_Science)
        self.label_N_energy_bins.setGeometry(QtCore.QRect(95, 130, 121, 21))
        self.label_N_energy_bins.setObjectName("label_N_energy_bins")
        self.spinBox_SED_N_energy_bins = QtWidgets.QSpinBox(self.groupBox_Science)
        self.spinBox_SED_N_energy_bins.setGeometry(QtCore.QRect(40, 130, 48, 26))
        self.spinBox_SED_N_energy_bins.setMinimum(3)
        self.spinBox_SED_N_energy_bins.setProperty("value", 10)
        self.spinBox_SED_N_energy_bins.setObjectName("spinBox_SED_N_energy_bins")
        self.checkBox_use_local_index = QtWidgets.QCheckBox(self.groupBox_Science)
        self.checkBox_use_local_index.setGeometry(QtCore.QRect(210, 130, 131, 23))
        self.checkBox_use_local_index.setChecked(True)
        self.checkBox_use_local_index.setObjectName("checkBox_use_local_index")
        self.checkBox_use_local_index.setChecked(False)
        self.label_redshift = QtWidgets.QLabel(self.groupBox_Science)
        self.label_redshift.setGeometry(QtCore.QRect(95, 160, 50, 21))
        self.label_redshift.setObjectName("label_redshift")
        self.white_box_redshift = QtWidgets.QLineEdit(self.groupBox_Science)
        self.white_box_redshift.setGeometry(QtCore.QRect(40, 160, 48, 25))
        self.white_box_redshift.setObjectName("white_box_redshift")
        self.comboBox_redshift = QtWidgets.QComboBox(self.groupBox_Science)
        self.comboBox_redshift.setEnabled(True)
        self.comboBox_redshift.setGeometry(QtCore.QRect(150, 160, 185, 25))
        self.comboBox_redshift.setObjectName("comboBox_redshift")
        self.comboBox_redshift.addItem("")
        self.comboBox_redshift.addItem("")
        self.comboBox_redshift.addItem("")
        self.comboBox_redshift.addItem("")
        self.comboBox_redshift.addItem("")
        self.comboBox_redshift.setCurrentIndex(4)
        self.comboBox_MCMC = QtWidgets.QComboBox(self.groupBox_Science)
        self.comboBox_MCMC.setEnabled(True)
        self.comboBox_MCMC.setGeometry(QtCore.QRect(40, 190, 80, 25))
        self.comboBox_MCMC.setObjectName("comboBox_MCMC")
        self.comboBox_MCMC.addItem("")
        self.comboBox_MCMC.addItem("")
        self.comboBox_MCMC.addItem("")
        self.comboBox_MCMC.addItem("")
        self.comboBox_MCMC.addItem("")
        self.comboBox_MCMC.setCurrentIndex(1)
        self.label_MCMC = QtWidgets.QLabel(self.groupBox_Science)
        self.label_MCMC.setGeometry(QtCore.QRect(125, 190, 50, 21))
        self.label_MCMC.setObjectName("label_MCMC")
        self.white_box_VHE = QtWidgets.QLineEdit(self.groupBox_Science)
        self.white_box_VHE.setGeometry(QtCore.QRect(170, 190, 135, 25))
        self.white_box_VHE.setObjectName("white_box_VHE")
        self.white_box_VHE.setText("Add VHE data?")
        self.toolButton_VHE = QtWidgets.QToolButton(self.groupBox_Science)
        self.toolButton_VHE.setGeometry(QtCore.QRect(310, 190, 26, 24))
        self.toolButton_VHE.setWhatsThis("")
        self.toolButton_VHE.setObjectName("toolButton_VHE")


        #Extension:
        self.checkBox_extension = QtWidgets.QCheckBox(self.groupBox_Science)
        self.checkBox_extension.setGeometry(QtCore.QRect(10, 215, 131, 23))
        self.checkBox_extension.setObjectName("checkBox_extension")
        self.radioButton_disk = QtWidgets.QRadioButton(self.groupBox_Science)
        self.radioButton_disk.setEnabled(False)
        self.radioButton_disk.setGeometry(QtCore.QRect(40, 235, 61, 23))
        self.radioButton_disk.setChecked(True)
        self.radioButton_disk.setAutoExclusive(True)
        self.radioButton_disk.setObjectName("radioButton_disk")
        self.radioButton_2D_Gauss = QtWidgets.QRadioButton(self.groupBox_Science)
        self.radioButton_2D_Gauss.setEnabled(False)
        self.radioButton_2D_Gauss.setGeometry(QtCore.QRect(95, 235, 91, 23))
        self.radioButton_2D_Gauss.setChecked(False)
        self.radioButton_2D_Gauss.setAutoExclusive(True)
        self.radioButton_2D_Gauss.setObjectName("radioButton_2D_Gauss")
        self.label_Extension_max_size = QtWidgets.QLabel(self.groupBox_Science)
        self.label_Extension_max_size.setEnabled(False)
        self.label_Extension_max_size.setGeometry(QtCore.QRect(252, 230, 81, 31))
        self.label_Extension_max_size.setObjectName("label_Extension_max_size")
        self.doubleSpinBox_extension_max_size = QtWidgets.QDoubleSpinBox(self.groupBox_Science)
        self.doubleSpinBox_extension_max_size.setEnabled(False)
        self.doubleSpinBox_extension_max_size.setGeometry(QtCore.QRect(200, 235, 48, 21))
        self.doubleSpinBox_extension_max_size.setProperty("value", 1.0)
        self.doubleSpinBox_extension_max_size.setMinimum(0.1)
        self.doubleSpinBox_extension_max_size.setObjectName("doubleSpinBox_extension_max_size")
        
        
        #Relocalize and TS maps:
        self.checkBox_relocalize = QtWidgets.QCheckBox(self.groupBox_Science)
        self.checkBox_relocalize.setGeometry(QtCore.QRect(10, 255, 131, 23))
        self.checkBox_relocalize.setObjectName("checkBox_relocalize")

        self.checkBox_TSmap = QtWidgets.QCheckBox(self.groupBox_Science)
        self.checkBox_TSmap.setGeometry(QtCore.QRect(10, 280, 151, 23))
        self.checkBox_TSmap.setChecked(True)
        self.checkBox_TSmap.setObjectName("checkBox_TSmap")
        self.checkBox_residual_TSmap = QtWidgets.QCheckBox(self.groupBox_Science)
        self.checkBox_residual_TSmap.setGeometry(QtCore.QRect(210, 300, 131, 23))
        self.checkBox_residual_TSmap.setChecked(True)
        self.checkBox_residual_TSmap.setObjectName("checkBox_residual_TSmap")
        self.doubleSpinBox_Photon_index_TS = QtWidgets.QDoubleSpinBox(self.groupBox_Science)
        self.doubleSpinBox_Photon_index_TS.setGeometry(QtCore.QRect(40, 300, 59, 26))
        self.doubleSpinBox_Photon_index_TS.setMinimum(0.5)
        self.doubleSpinBox_Photon_index_TS.setProperty("value", 2.0)
        self.doubleSpinBox_Photon_index_TS.setObjectName("doubleSpinBox_Photon_index_TS")
        self.label_photon_index_TS = QtWidgets.QLabel(self.groupBox_Science)
        self.label_photon_index_TS.setGeometry(QtCore.QRect(105, 300, 90, 21))
        self.label_photon_index_TS.setObjectName("label_photon_index_TS")
        
        
        
        ##################################################################
        ######## Group box fit fine-tune
        ##################################################################
        
        self.groupBox_fit_finetune = QtWidgets.QGroupBox(self.centralwidget)
        self.groupBox_fit_finetune.setGeometry(QtCore.QRect(10, 190, 282, 331))
        self.groupBox_fit_finetune.setObjectName("groupBox_fit_finetune")
        
        self.checkBox_delete_sources = QtWidgets.QCheckBox(self.groupBox_fit_finetune)
        self.checkBox_delete_sources.setGeometry(QtCore.QRect(10, 150, 141, 23))
        self.checkBox_delete_sources.setObjectName("checkBox_delete_sources")
        self.comboBox_change_model = QtWidgets.QComboBox(self.groupBox_fit_finetune)
        self.comboBox_change_model.setEnabled(False)
        self.comboBox_change_model.setGeometry(QtCore.QRect(30, 110, 101, 25))
        self.comboBox_change_model.setObjectName("comboBox_change_model")
        self.comboBox_change_model.addItem("")
        self.comboBox_change_model.addItem("")
        self.comboBox_change_model.addItem("")
        self.comboBox_change_model.addItem("")
        self.comboBox_change_model.addItem("")
        self.comboBox_change_model.addItem("")
        self.comboBox_change_model.addItem("")
        self.comboBox_change_model.addItem("")
        self.comboBox_change_model.addItem("")
        self.comboBox_change_model.addItem("")
        self.white_box_list_of_sources_to_delete = QtWidgets.QLineEdit(self.groupBox_fit_finetune)
        self.white_box_list_of_sources_to_delete.setEnabled(False)
        self.white_box_list_of_sources_to_delete.setGeometry(QtCore.QRect(30, 180, 101, 21))
        self.white_box_list_of_sources_to_delete.setObjectName("white_box_list_of_sources_to_delete")
        self.checkBox_diagnostic_plots = QtWidgets.QCheckBox(self.groupBox_fit_finetune)
        self.checkBox_diagnostic_plots.setGeometry(QtCore.QRect(10, 300, 151, 23))
        self.checkBox_diagnostic_plots.setChecked(True)
        self.checkBox_diagnostic_plots.setObjectName("checkBox_diagnostic_plots")
        self.label_minimum_separation = QtWidgets.QLabel(self.groupBox_fit_finetune)
        self.label_minimum_separation.setGeometry(QtCore.QRect(110, 272, 141, 17))
        self.label_minimum_separation.setObjectName("label_minimum_separation")
        self.doubleSpinBox_min_significance = QtWidgets.QDoubleSpinBox(self.groupBox_fit_finetune)
        self.doubleSpinBox_min_significance.setEnabled(True)
        self.doubleSpinBox_min_significance.setGeometry(QtCore.QRect(30, 240, 69, 21))
        self.doubleSpinBox_min_significance.setMinimum(1.0)
        self.doubleSpinBox_min_significance.setProperty("value", 5.0)
        self.doubleSpinBox_min_significance.setObjectName("doubleSpinBox_min_significance")
        self.doubleSpinBox_min_separation = QtWidgets.QDoubleSpinBox(self.groupBox_fit_finetune)
        self.doubleSpinBox_min_separation.setEnabled(True)
        self.doubleSpinBox_min_separation.setGeometry(QtCore.QRect(30, 270, 69, 21))
        self.doubleSpinBox_min_separation.setMinimum(0.1)
        self.doubleSpinBox_min_separation.setProperty("value", 0.5)
        self.doubleSpinBox_min_separation.setObjectName("doubleSpinBox_min_separation")
        self.checkBox_find_extra_sources = QtWidgets.QCheckBox(self.groupBox_fit_finetune)
        self.checkBox_find_extra_sources.setGeometry(QtCore.QRect(10, 220, 221, 23))
        self.checkBox_find_extra_sources.setChecked(True)
        self.checkBox_find_extra_sources.setObjectName("checkBox_find_extra_sources")
        self.label_min_significance = QtWidgets.QLabel(self.groupBox_fit_finetune)
        self.label_min_significance.setGeometry(QtCore.QRect(110, 242, 141, 17))
        self.label_min_significance.setObjectName("label_min_significance")
        self.checkBox_change_model = QtWidgets.QCheckBox(self.groupBox_fit_finetune)
        self.checkBox_change_model.setGeometry(QtCore.QRect(10, 90, 131, 21))
        self.checkBox_change_model.setObjectName("checkBox_change_model")
        self.checkBox_minimizer = QtWidgets.QCheckBox(self.groupBox_fit_finetune)
        self.checkBox_minimizer.setGeometry(QtCore.QRect(10, 30, 131, 21))
        self.checkBox_minimizer.setObjectName("checkBox_minimizer")
        self.comboBox_minimizer = QtWidgets.QComboBox(self.groupBox_fit_finetune)
        self.comboBox_minimizer.setEnabled(False)
        self.comboBox_minimizer.setGeometry(QtCore.QRect(30, 50, 101, 25))
        self.comboBox_minimizer.setObjectName("comboBox_minimizer")
        self.comboBox_minimizer.addItem("")
        self.comboBox_minimizer.addItem("")
        self.comboBox_minimizer.addItem("")
        self.comboBox_minimizer.addItem("")
        self.comboBox_minimizer.addItem("")
        self.radioButton_free_source_radius = QtWidgets.QRadioButton(self.groupBox_fit_finetune)
        self.radioButton_free_source_radius.setGeometry(QtCore.QRect(150, 50, 81, 23))
        self.radioButton_free_source_radius.setChecked(True)
        self.radioButton_free_source_radius.setAutoExclusive(True)
        self.radioButton_free_source_radius.setObjectName("radioButton_free_source_radius")
        self.radioButton_free_source_radius_customized = QtWidgets.QRadioButton(self.groupBox_fit_finetune)
        self.radioButton_free_source_radius_customized.setGeometry(QtCore.QRect(150, 70, 112, 23))
        self.radioButton_free_source_radius_customized.setAutoExclusive(True)
        self.radioButton_free_source_radius_customized.setObjectName("radioButton_free_source_radius_customized")
        self.label_free_src_radius = QtWidgets.QLabel(self.groupBox_fit_finetune)
        self.label_free_src_radius.setGeometry(QtCore.QRect(150, 30, 131, 17))
        self.label_free_src_radius.setObjectName("label_free_src_radius")
        self.label_radius = QtWidgets.QLabel(self.groupBox_fit_finetune)
        self.label_radius.setEnabled(False)
        self.label_radius.setGeometry(QtCore.QRect(220, 100, 71, 21))
        self.label_radius.setObjectName("label_radius")
        self.checkBox_only_norm = QtWidgets.QCheckBox(self.groupBox_fit_finetune)
        self.checkBox_only_norm.setEnabled(False)
        self.checkBox_only_norm.setGeometry(QtCore.QRect(170, 130, 121, 21))
        self.checkBox_only_norm.setChecked(False)
        self.checkBox_only_norm.setObjectName("checkBox_only_norm")
        self.white_box_radius = QtWidgets.QLineEdit(self.groupBox_fit_finetune)
        self.white_box_radius.setEnabled(False)
        self.white_box_radius.setGeometry(QtCore.QRect(170, 100, 41, 21))
        self.white_box_radius.setObjectName("white_box_radius")
        self.line_4 = QtWidgets.QFrame(self.groupBox_fit_finetune)
        self.line_4.setGeometry(QtCore.QRect(0, 210, 281, 16))
        self.line_4.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_4.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_4.setObjectName("line_4")
        self.checkBox_freeze_gal = QtWidgets.QCheckBox(self.groupBox_fit_finetune)
        self.checkBox_freeze_gal.setEnabled(False)
        self.checkBox_freeze_gal.setGeometry(QtCore.QRect(170, 150, 121, 21))
        self.checkBox_freeze_gal.setChecked(False)
        self.checkBox_freeze_gal.setObjectName("checkBox_freeze_gal")
        self.checkBox_freeze_iso = QtWidgets.QCheckBox(self.groupBox_fit_finetune)
        self.checkBox_freeze_iso.setEnabled(False)
        self.checkBox_freeze_iso.setGeometry(QtCore.QRect(170, 170, 121, 21))
        self.checkBox_freeze_iso.setChecked(False)
        self.checkBox_freeze_iso.setObjectName("checkBox_freeze_iso")
        self.checkBox_freeze_spec_shape = QtWidgets.QCheckBox(self.groupBox_fit_finetune)
        self.checkBox_freeze_spec_shape.setEnabled(True)
        self.checkBox_freeze_spec_shape.setGeometry(QtCore.QRect(150, 190, 141, 21))
        self.checkBox_freeze_spec_shape.setChecked(False)
        self.checkBox_freeze_spec_shape.setObjectName("checkBox_freeze_spec_shape")
        self.label_output_format = QtWidgets.QLabel(self.groupBox_fit_finetune)
        self.label_output_format.setGeometry(QtCore.QRect(150, 290, 111, 17))
        self.label_output_format.setObjectName("label_output_format")
        self.comboBox_output_format = QtWidgets.QComboBox(self.groupBox_fit_finetune)
        self.comboBox_output_format.setEnabled(True)
        self.comboBox_output_format.setGeometry(QtCore.QRect(150, 308, 101, 20))
        self.comboBox_output_format.setObjectName("comboBox_output_format")
        self.comboBox_output_format.addItem("")
        self.comboBox_output_format.addItem("")
        
        

        ##########################################################################
        ###### Log-box region
        ##########################################################################

        logbox_position = 650
        self.large_white_box_Log = QtWidgets.QPlainTextEdit(self.centralwidget)
        self.large_white_box_Log.setGeometry(QtCore.QRect(logbox_position, 210, 450, 261))
        self.large_white_box_Log.setObjectName("large_white_box_Log")
        self.label_Log = QtWidgets.QLabel(self.centralwidget)
        self.label_Log.setGeometry(QtCore.QRect(logbox_position, 190, 71, 20))
        self.label_Log.setObjectName("label_Log")
        self.picture = QtWidgets.QLabel(self.centralwidget)
        self.picture.setEnabled(True)
        self.picture.setGeometry(QtCore.QRect(logbox_position+110, 230, 231, 231))
        self.picture.setFont(font)
        self.picture.setMouseTracking(False)
        self.picture.setAutoFillBackground(False)
        self.picture.setText("")
        self.picture.setPixmap(QtGui.QPixmap(str(libpath / "fermi.png")))
        self.picture.setScaledContents(True)
        self.picture.setWordWrap(False)
        self.picture.setObjectName("picture")
        self.pushButton_Go = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_Go.setGeometry(QtCore.QRect(logbox_position+260, 480, 191, 41))
        self.pushButton_Go.setObjectName("pushButton_Go")
        self.label_output_dir = QtWidgets.QLabel(self.centralwidget)
        self.label_output_dir.setGeometry(QtCore.QRect(logbox_position, 470, 131, 20))
        self.label_output_dir.setObjectName("label_output_dir")
        self.toolButton_output_dir = QtWidgets.QToolButton(self.centralwidget)
        self.toolButton_output_dir.setGeometry(QtCore.QRect(logbox_position+110, 490, 26, 24))
        self.toolButton_output_dir.setObjectName("toolButton_output_dir")
        self.white_box_output_dir = QtWidgets.QLineEdit(self.centralwidget)
        self.white_box_output_dir.setGeometry(QtCore.QRect(logbox_position, 490, 101, 25))
        self.white_box_output_dir.setObjectName("white_box_output_dir")
        

        ##########################################################################
        ###### Column 1 under "Standard"
        ##########################################################################

        self.label_Catalog = QtWidgets.QLabel(self.centralwidget)
        self.label_Catalog.setGeometry(QtCore.QRect(40, 40, 61, 21))
        self.label_Catalog.setObjectName("label_Catalog")        
        self.comboBox_Catalog = QtWidgets.QComboBox(self.centralwidget)
        self.comboBox_Catalog.setGeometry(QtCore.QRect(40, 60, 91, 25))
        self.comboBox_Catalog.setObjectName("comboBox_Catalog")
        self.comboBox_Catalog.addItem("")
        self.comboBox_Catalog.addItem("")
        self.comboBox_Catalog.addItem("")
        self.comboBox_Catalog.addItem("")
        self.label_is_it_cataloged = QtWidgets.QLabel(self.centralwidget)
        self.label_is_it_cataloged.setGeometry(QtCore.QRect(40, 90, 151, 21))
        self.label_is_it_cataloged.setObjectName("label_is_it_cataloged")
        self.comboBox_is_it_cataloged = QtWidgets.QComboBox(self.centralwidget)
        self.comboBox_is_it_cataloged.setGeometry(QtCore.QRect(40, 110, 47, 25))
        self.comboBox_is_it_cataloged.setObjectName("comboBox_is_it_cataloged")
        self.comboBox_is_it_cataloged.addItem("")
        self.comboBox_is_it_cataloged.addItem("")
        self.white_box_target_name = QtWidgets.QLineEdit(self.centralwidget)
        self.white_box_target_name.setEnabled(False)
        self.white_box_target_name.setGeometry(QtCore.QRect(92, 110, 90, 25))
        self.white_box_target_name.setObjectName("white_box_target_name")

        ##########################################################################
        ###### Column 2 under "Standard"
        ##########################################################################
        
        col2_position = 200
        self.label_energy = QtWidgets.QLabel(self.centralwidget)
        self.label_energy.setGeometry(QtCore.QRect(col2_position, 90, 111, 17))
        self.label_energy.setObjectName("label_energy")
        self.white_box_energy = QtWidgets.QLineEdit(self.centralwidget)
        self.white_box_energy.setGeometry(QtCore.QRect(col2_position, 110, 111, 25))
        self.white_box_energy.setObjectName("white_box_energy")
        self.label_RAandDec = QtWidgets.QLabel(self.centralwidget)
        self.label_RAandDec.setGeometry(QtCore.QRect(col2_position, 40, 81, 17))
        self.label_RAandDec.setObjectName("label_RAandDec")
        self.white_box_RAandDec = QtWidgets.QLineEdit(self.centralwidget)
        self.white_box_RAandDec.setGeometry(QtCore.QRect(col2_position, 60, 111, 25))
        self.white_box_RAandDec.setObjectName("white_box_RAandDec")
        self.vertical_line_col2 = QtWidgets.QFrame(self.centralwidget)
        self.vertical_line_col2.setGeometry(QtCore.QRect(col2_position-20, 40, 20, 101))
        self.vertical_line_col2.setFrameShape(QtWidgets.QFrame.VLine)
        self.vertical_line_col2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.vertical_line_col2.setObjectName("vertical_line_col2")

        ##########################################################################
        ###### Column 3 under "Standard"
        ##########################################################################

        col3_position = 327
        self.vertical_line_col3 = QtWidgets.QFrame(self.centralwidget)
        self.vertical_line_col3.setGeometry(QtCore.QRect(col3_position-20, 40, 20, 101))
        self.vertical_line_col3.setFrameShape(QtWidgets.QFrame.VLine)
        self.vertical_line_col3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.vertical_line_col3.setObjectName("vertical_line_col3")
        self.dateTimeEdit_2 = QtWidgets.QDateTimeEdit(self.centralwidget)
        self.dateTimeEdit_2.setGeometry(QtCore.QRect(col3_position, 110, 171, 21))
        self.dateTimeEdit_2.setDateTime(QtCore.QDateTime(QtCore.QDate(2008, 10, 14), QtCore.QTime(15, 43, 0)))
        self.dateTimeEdit_2.setMinimumDateTime(QtCore.QDateTime(QtCore.QDate(2008, 8, 5), QtCore.QTime(00, 00, 00)))
        self.dateTimeEdit_2.setCalendarPopup(False)
        self.dateTimeEdit_2.setObjectName("dateTimeEdit_2")
        self.dateTimeEdit = QtWidgets.QDateTimeEdit(self.centralwidget)
        self.dateTimeEdit.setGeometry(QtCore.QRect(col3_position, 60, 171, 21))
        self.dateTimeEdit.setDateTime(QtCore.QDateTime(QtCore.QDate(2008, 8, 4), QtCore.QTime(15, 43, 36)))
        self.dateTimeEdit.setDate(QtCore.QDate(2008, 8, 4))
        self.dateTimeEdit.setTime(QtCore.QTime(15, 43, 36))
        self.dateTimeEdit.setMinimumDateTime(QtCore.QDateTime(QtCore.QDate(2008, 8, 4), QtCore.QTime(15, 43, 36)))
        self.dateTimeEdit.setObjectName("dateTimeEdit")
        self.label_start_time = QtWidgets.QLabel(self.centralwidget)
        self.label_start_time.setGeometry(QtCore.QRect(col3_position, 40, 111, 17))
        self.label_start_time.setObjectName("label_start_time")
        self.label_end_time = QtWidgets.QLabel(self.centralwidget)
        self.label_end_time.setGeometry(QtCore.QRect(col3_position, 90, 111, 17))
        self.label_end_time.setObjectName("label_end_time")

        ##########################################################################
        ###### Column 4 under "Standard"
        ##########################################################################

        col4_position = 515
        self.line = QtWidgets.QFrame(self.centralwidget)
        self.line.setGeometry(QtCore.QRect(col4_position-20, 40, 20, 101))
        self.line.setFrameShape(QtWidgets.QFrame.VLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.white_box_photon_dir = QtWidgets.QLineEdit(self.centralwidget)
        self.white_box_photon_dir.setGeometry(QtCore.QRect(col4_position, 110, 131, 25))
        self.white_box_photon_dir.setObjectName("white_box_photon_dir")
        self.white_box_spacecraft_file = QtWidgets.QLineEdit(self.centralwidget)
        self.white_box_spacecraft_file.setGeometry(QtCore.QRect(col4_position, 60, 131, 25))
        self.white_box_spacecraft_file.setObjectName("white_box_spacecraft_file")
        self.label_dir_photons = QtWidgets.QLabel(self.centralwidget)
        self.label_dir_photons.setGeometry(QtCore.QRect(col4_position, 90, 131, 20))
        self.label_dir_photons.setObjectName("label_dir_photons")
        self.label_dir_spacecraft = QtWidgets.QLabel(self.centralwidget)
        self.label_dir_spacecraft.setGeometry(QtCore.QRect(col4_position, 40, 111, 17))
        self.label_dir_spacecraft.setObjectName("label")
        self.toolButton_spacecraft = QtWidgets.QToolButton(self.centralwidget)
        self.toolButton_spacecraft.setGeometry(QtCore.QRect(col4_position+140, 60, 26, 24))
        self.toolButton_spacecraft.setWhatsThis("")
        self.toolButton_spacecraft.setObjectName("toolButton_spacecraft")
        self.toolButton_photons = QtWidgets.QToolButton(self.centralwidget)
        self.toolButton_photons.setGeometry(QtCore.QRect(col4_position+140, 110, 26, 24))
        self.toolButton_photons.setObjectName("toolButton_photons")
        self.pushButton_Download_SC = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_Download_SC.setGeometry(QtCore.QRect(col4_position+170, 60, 26, 24))
        self.pushButton_Download_SC.setObjectName("pushButton_Download_SC")
        self.pushButton_Download_SC.setIcon(QtGui.QIcon(str(libpath/'download_logo.png')))
        self.pushButton_Download_Photons = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_Download_Photons.setGeometry(QtCore.QRect(col4_position+170, 110, 26, 24))
        self.pushButton_Download_Photons.setObjectName("pushButton_Download_Photons")
        self.pushButton_Download_Photons.setIcon(QtGui.QIcon(str(libpath/'download_logo.png')))

        ##########################################################################
        ###### Column 5 under "Standard"
        ##########################################################################

        col5_position = 730
        self.toolButton_External_ltcube = QtWidgets.QToolButton(self.centralwidget)
        self.toolButton_External_ltcube.setEnabled(False)
        self.toolButton_External_ltcube.setGeometry(QtCore.QRect(col5_position+160, 110, 26, 24))
        self.toolButton_External_ltcube.setObjectName("toolButton_External_ltcube")
        self.white_box_External_ltcube = QtWidgets.QLineEdit(self.centralwidget)
        self.white_box_External_ltcube.setEnabled(False)
        self.white_box_External_ltcube.setGeometry(QtCore.QRect(col5_position+20, 110, 131, 25))
        self.white_box_External_ltcube.setObjectName("white_box_External_ltcube")
        self.white_box_Diffuse_dir = QtWidgets.QLineEdit(self.centralwidget)
        self.white_box_Diffuse_dir.setGeometry(QtCore.QRect(col5_position, 60, 151, 25))
        self.white_box_Diffuse_dir.setObjectName("white_box_Diffuse_dir")
        self.toolButton_dir_diffuse = QtWidgets.QToolButton(self.centralwidget)
        self.toolButton_dir_diffuse.setGeometry(QtCore.QRect(col5_position+160, 60, 26, 24))
        self.toolButton_dir_diffuse.setObjectName("toolButton_dir_diffuse")
        self.label_dir_diffuse = QtWidgets.QLabel(self.centralwidget)
        self.label_dir_diffuse.setGeometry(QtCore.QRect(col5_position, 40, 161, 17))
        self.label_dir_diffuse.setObjectName("label_dir_diffuse")
        self.checkBox_External_ltcube = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_External_ltcube.setGeometry(QtCore.QRect(col5_position, 90, 161, 23))
        self.checkBox_External_ltcube.setObjectName("checkBox_External_ltcube")
        self.vertical_line_col5 = QtWidgets.QFrame(self.centralwidget)
        self.vertical_line_col5.setGeometry(QtCore.QRect(col5_position-20, 40, 20, 101))
        self.vertical_line_col5.setFrameShape(QtWidgets.QFrame.VLine)
        self.vertical_line_col5.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.vertical_line_col5.setObjectName("vertical_line_col5")
        self.pushButton_Download_diffuse = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_Download_diffuse.setGeometry(QtCore.QRect(col5_position+190, 60, 26, 24))
        self.pushButton_Download_diffuse.setObjectName("pushButton_Download_diffuse")
        self.pushButton_Download_diffuse.setIcon(QtGui.QIcon(str(libpath/'download_logo.png')))


        ##########################################################################
        ###### Column 6 under "Standard"
        ##########################################################################

        col6_position = 962
        self.vertical_line_col6 = QtWidgets.QFrame(self.centralwidget)
        self.vertical_line_col6.setGeometry(QtCore.QRect(col6_position-20, 40, 20, 101))
        self.vertical_line_col6.setFrameShape(QtWidgets.QFrame.VLine)
        self.vertical_line_col6.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.vertical_line_col6.setObjectName("vertical_line_col6")
        self.label_limits = QtWidgets.QLabel(self.centralwidget)
        self.label_limits.setGeometry(QtCore.QRect(col6_position, 40, 111, 17))
        self.label_limits.setObjectName("label_limits")
        self.checkBox_highest_resolution = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_highest_resolution.setGeometry(QtCore.QRect(col6_position, 60, 161, 23))
        self.checkBox_highest_resolution.setObjectName("checkBox_highest_resolution")
        self.checkBox_high_sensitivity = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_high_sensitivity.setGeometry(QtCore.QRect(col6_position, 80, 161, 23))
        self.checkBox_high_sensitivity.setObjectName("checkBox_high_sensitivity")
        self.pushButton_config = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_config.setGeometry(QtCore.QRect(col6_position, 110, 140, 24))
        self.pushButton_config.setObjectName("pushButton_config")



        
        ##########################################################################
        ###### Custom
        ##########################################################################

        self.radioButton_Custom = QtWidgets.QRadioButton(self.centralwidget)
        self.radioButton_Custom.setEnabled(True)
        self.radioButton_Custom.setGeometry(QtCore.QRect(20, 140, 91, 23))
        self.radioButton_Custom.setChecked(False)
        self.radioButton_Custom.setAutoExclusive(True)
        self.radioButton_Custom.setObjectName("radioButton_Custom")

        self.toolButton_Custom = QtWidgets.QToolButton(self.centralwidget)
        self.toolButton_Custom.setEnabled(False)
        self.toolButton_Custom.setGeometry(QtCore.QRect(260, 160, 26, 24))
        self.toolButton_Custom.setObjectName("toolButton_Custom")
        self.white_box_config_file = QtWidgets.QLineEdit(self.centralwidget)
        self.white_box_config_file.setEnabled(False)
        self.white_box_config_file.setGeometry(QtCore.QRect(40, 160, 211, 25))
        self.white_box_config_file.setObjectName("white_box_config_file")



        
        




        self.line_6 = QtWidgets.QFrame(self.centralwidget)
        self.line_6.setGeometry(QtCore.QRect(140, 220, 20, 181))
        self.line_6.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_6.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_6.setObjectName("line_6")
        self.line_7 = QtWidgets.QFrame(self.centralwidget)
        self.line_7.setGeometry(QtCore.QRect(140, 480, 20, 41))
        self.line_7.setFrameShape(QtWidgets.QFrame.VLine)
        self.line_7.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_7.setObjectName("line_7")



        self.pushButton_Go.raise_()
        self.pushButton_config.raise_()
        self.pushButton_Download_SC.raise_()
        self.pushButton_Download_Photons.raise_()
        self.pushButton_Download_diffuse.raise_()
        self.progressBar.raise_()
        self.label_Log.raise_()
        self.groupBox_Science.raise_()
        self.groupBox_fit_finetune.raise_()
        self.large_white_box_Log.raise_()
        self.toolButton_External_ltcube.raise_()
        self.white_box_External_ltcube.raise_()
        self.white_box_Diffuse_dir.raise_()
        self.toolButton_dir_diffuse.raise_()
        self.label_dir_diffuse.raise_()
        self.checkBox_External_ltcube.raise_()
        self.checkBox_highest_resolution.raise_()
        self.vertical_line_col3.raise_()
        self.vertical_line_col5.raise_()
        self.dateTimeEdit_2.raise_()
        self.dateTimeEdit.raise_()
        self.label_start_time.raise_()
        self.label_end_time.raise_()
        self.label_limits.raise_()
        self.line.raise_()
        self.toolButton_output_dir.raise_()
        self.white_box_output_dir.raise_()
        self.white_box_redshift.raise_()
        self.white_box_TS_threshold.raise_()
        self.white_box_photon_dir.raise_()
        self.white_box_spacecraft_file.raise_()
        self.white_box_VHE.raise_()
        self.label_dir_photons.raise_()
        self.label_dir_spacecraft.raise_()
        self.label_output_dir.raise_()
        self.toolButton_spacecraft.raise_()
        self.toolButton_photons.raise_()
        self.toolButton_Custom.raise_()
        self.toolButton_VHE.raise_()
        self.label_energy.raise_()
        self.white_box_energy.raise_()
        self.label_Catalog.raise_()
        self.label_RAandDec.raise_()
        self.comboBox_Catalog.raise_()
        self.white_box_config_file.raise_()
        self.white_box_RAandDec.raise_()
        self.label_configuration_file.raise_()
        self.picture.raise_()
        self.label_is_it_cataloged.raise_()
        self.comboBox_is_it_cataloged.raise_()
        self.vertical_line_col2.raise_()
        self.line_6.raise_()
        self.line_7.raise_()
        mainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(mainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 870, 22))
        self.menubar.setObjectName("menubar")
        self.menuTutorial = QtWidgets.QMenu(self.menubar)
        self.menuTutorial.setObjectName("menuTutorial")
        self.menuCredits = QtWidgets.QMenu(self.menubar)
        self.menuCredits.setObjectName("menuCredits")
        mainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(mainWindow)
        self.statusbar.setObjectName("statusbar")
        mainWindow.setStatusBar(self.statusbar)
        self.actionOpen = QtWidgets.QAction(mainWindow)
        self.actionOpen.setShortcut("")
        self.actionOpen.setObjectName("actionOpen")
        self.actionSave = QtWidgets.QAction(mainWindow)
        self.actionSave.setObjectName("actionSave")
        self.actionCopy = QtWidgets.QAction(mainWindow)
        self.actionCopy.setObjectName("actionCopy")
        self.actionPaste = QtWidgets.QAction(mainWindow)
        self.actionPaste.setObjectName("actionPaste")
        self.actionOpen_Tutorial = QtWidgets.QAction(mainWindow)
        self.actionOpen_Tutorial.setObjectName("actionOpen_Tutorial")
        self.actionSee_credits = QtWidgets.QAction(mainWindow)
        self.actionSee_credits.setObjectName("actionSee_credits")
        self.actionLoad_state = QtWidgets.QAction(mainWindow)
        self.actionLoad_state.setObjectName("actionLoad_state")
        self.menuTutorial.addAction(self.actionOpen_Tutorial)
        self.menuTutorial.addAction(self.actionLoad_state)
        self.menuCredits.addAction(self.actionSee_credits)
        self.menubar.addAction(self.menuTutorial.menuAction())
        self.menubar.addAction(self.menuCredits.menuAction())







        
        self.retranslateUi(mainWindow)
        QtCore.QMetaObject.connectSlotsByName(mainWindow)
        
        #Click GO!
        self.pushButton_Go.setToolTip('Click to start the analysis.')
        self.pushButton_Go.clicked.connect(self.runLongTask)
        self.pushButton_config.setToolTip('Click to generate the yaml configuration file (optional).')
        self.pushButton_config.clicked.connect(self.click_to_generateConfig)
        self.pushButton_Download_SC.setToolTip('Click to download the spacecraft data file.')
        self.pushButton_Download_SC.clicked.connect(self.popup_download_SC)
        self.pushButton_Download_Photons.setToolTip('Click to download the photon data files.')
        self.pushButton_Download_Photons.clicked.connect(self.popup_download_Photons)
        self.pushButton_Download_diffuse.setToolTip("Click to download the latest diffuse files.\nThese files can be used for all sources in the sky.\nYou don't need to download them multiple times.")
        self.pushButton_Download_diffuse.clicked.connect(self.popup_download_Diffuse)
        #self.actionOpen_Tutorial.connect(self.popup_tutorial)
        self.actionOpen_Tutorial.triggered.connect(self.popup_tutorial)
        self.actionSee_credits.triggered.connect(self.popup_credits)
        self.actionLoad_state.triggered.connect(self.load_GUIstate)
            
        ############ Loading main files:
        self.toolButton_spacecraft.setCheckable(True)
        self.toolButton_spacecraft.setToolTip('You can download the spacecraft file from https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi')
        self.toolButton_spacecraft.clicked.connect(self.browsefiles)
        self.toolButton_photons.setCheckable(True)
        self.toolButton_photons.setToolTip('You can download the photon files from https://fermi.gsfc.nasa.gov/cgi-bin/ssc/LAT/LATDataQuery.cgi')
        self.toolButton_photons.clicked.connect(self.browsefiles)
        self.toolButton_Custom.setCheckable(True)
        self.toolButton_Custom.setToolTip('Please upload your own config.yaml file.')
        self.toolButton_Custom.clicked.connect(self.browsefiles)
        self.toolButton_dir_diffuse.setCheckable(True)
        self.toolButton_dir_diffuse.setToolTip('You can download the background files from https://fermi.gsfc.nasa.gov/ssc/data/access/lat/BackgroundModels.html')
        self.toolButton_dir_diffuse.clicked.connect(self.browsefiles)
        self.toolButton_output_dir.setCheckable(True)
        self.toolButton_output_dir.clicked.connect(self.browsefiles)
        self.toolButton_External_ltcube.setCheckable(True)
        self.toolButton_External_ltcube.clicked.connect(self.browsefiles)
        self.toolButton_VHE.setCheckable(True)
        self.toolButton_VHE.clicked.connect(self.browsefiles)
        self.toolButton_VHE.setToolTip('Here you can select a fits table with the VHE data. This is not mandatory for the MCMC.\nCheck the GitHub of easyfermi for details on the format of this table.')
        self.toolButton_External_ltcube.setToolTip('Here you can select a txt file listing one or more ltcube files.\neasyfermi automatically saves this file as "ltcube_list.txt" once you run an analysis.')

        self.white_box_list_of_sources_to_delete.setToolTip("e.g.: 4FGL J1222.5+0414,4FGL J1219.7+0444,4FGL ...")
        
        
        ###### Activating/deactivating options
        self.checkBox_LC.clicked.connect(self.activate)
        self.checkBox_adaptive_binning.clicked.connect(self.activate)
        self.checkBox_SED.clicked.connect(self.activate)
        self.checkBox_extension.clicked.connect(self.activate)
        self.checkBox_TSmap.clicked.connect(self.activate)
        self.checkBox_find_extra_sources.clicked.connect(self.activate)
        self.checkBox_External_ltcube.clicked.connect(self.activate)
        self.checkBox_highest_resolution.clicked.connect(self.activate)
        self.checkBox_high_sensitivity.clicked.connect(self.activate)
        self.checkBox_delete_sources.clicked.connect(self.activate)
        self.checkBox_change_model.clicked.connect(self.activate)
        self.checkBox_minimizer.clicked.connect(self.activate)
        self.radioButton_disk.clicked.connect(self.activate)
        self.radioButton_2D_Gauss.clicked.connect(self.activate)
        self.radioButton_Standard.clicked.connect(self.activate)
        self.radioButton_Custom.clicked.connect(self.activate)
        self.radioButton_free_source_radius.clicked.connect(self.activate)
        self.radioButton_free_source_radius_customized.clicked.connect(self.activate)

        self.comboBox_is_it_cataloged.activated.connect(self.activate)
        self.comboBox_is_it_cataloged.setToolTip("Is your target listed in the catalog selected above?")
        self.comboBox_change_model.setToolTip("Select a spectral model for your target.")
        self.comboBox_minimizer.setToolTip("Select an optimizer for the fit.")
        self.comboBox_redshift.setToolTip("Select an EBL absorption model.")
        self.comboBox_MCMC.setToolTip("Choose a spectral model for the MCMC.\n\nWalkers = 300\nIterations = 500")
        self.comboBox_output_format.setToolTip("Select the output format for the main plots (i.e. SED, light curve etc).")
        self.checkBox_only_norm.setToolTip("Check if you wish that only the normalizations can vary.")
        self.checkBox_freeze_gal.setToolTip("Freeze the Galactic diffuse model.")
        self.checkBox_freeze_iso.setToolTip("Freeze the isotropic diffuse model.")
        self.checkBox_freeze_spec_shape.setToolTip("If checked, you will freeze the spectral shape of the target.")
        self.checkBox_find_extra_sources.setToolTip("Check to look for extra sources in the ROI, i.e. sources not listed in the adopted catalog.")
        self.checkBox_residual_TSmap.setToolTip("If checked, easyfermi will also compute the residuals TS map.")
        self.label_photon_index_TS.setToolTip("The photon index of the test source adopted in the TS and excess maps.")
        self.spinBox_N_cores_LC.setToolTip("Multiprocessing is available only for Linux OS.")
        self.checkBox_relocalize.setToolTip("Check this to find the optimal R.A. and Dec. of the target's gamma-ray emission.")
        self.white_box_radius.setToolTip("Only the sources within this radius will have free parameters during the fit.")
        self.label_dir_diffuse.setToolTip("Directory containing the diffuse emission files (e.g. gll_iem_v07.fits and iso_P8R3_SOURCE_V3_v1.txt).")
        self.label_RAandDec.setToolTip("J2000 coordinates in degrees.\nIf you are connected to the internet, you can type the target name, e.g.: M31 or NGC 1275.")
        self.checkBox_highest_resolution.setToolTip('Check this box to increase angular resolution at the cost of decreasing by nearly half the sensitivity.\nThis method will select only PSF 2 and 3 events.\nIf combined with the "Improve sensitivity" option below, it will select:\n- PSF 2 and 3 events below 500 MeV.\n- PSF 1, 2 and 3 events between 500 MeV and 1000 MeV.\n- All events above 1000 MeV.')
        self.checkBox_high_sensitivity.setToolTip("Check this box to slightly improve sensitivity (less than 5%) at high energies at the cost of a longer analysis.\nThis will be useful as long as E_min < 1000 MeV.")
        self.checkBox_diagnostic_plots.setToolTip("Check this box to compute a series of diagnostic plots (e.g.: distance between the target and the Sun).")
        self.checkBox_adaptive_binning.setToolTip("Check this box to compute an adaptive-binning light curve.\nThis process can take quite a long time.\nTry running it with more than 1 core.")
        self.checkBox_use_local_index.setToolTip("Check to use a power-law approximation to the shape of the global spectrum in each energy bin. If not checked, a constant index will be used.")
        self.label_redshift.setToolTip("If you want to take EBL absorption into account, you can put the cosmological redshift of the target here.\nIf z = 0.0, the EBL absorption correction is not applied.")
        self.white_box_redshift.setToolTip("If you want to take EBL absorption into account, you can put the cosmological redshift of the target here.\nIf z = 0.0, the EBL absorption correction is not applied.")
        self.white_box_VHE.setToolTip('Optional feature:\nHere you can select a fits table with the VHE data.\nCheck the GitHub of easyfermi for details on the format of this table.')
        self.white_box_TS_threshold.setToolTip("This is the average TS for each one of the adaptive time bins. We recommend keeping it > 50.")
        self.spinBox_N_iter.setToolTip("The greater the number of iterations, the greater the temporal resolution of the light curve.")
        self.white_box_target_name.setToolTip("Please insert a nickname for your target avoiding blank spaces.")


    def retranslateUi(self, mainWindow):

        """
        Here we tell the buttons/checkboxes/toolboxes/etc what they must display on the easyfermi window.
        
        """
        _translate = QtCore.QCoreApplication.translate
        mainWindow.setWindowIcon(QtGui.QIcon(str(libpath/'easyFermiIcon.png')))
        mainWindow.setWindowTitle(_translate("mainWindow", "easyfermi"))
        self.pushButton_Go.setText(_translate("mainWindow", "Go!"))
        self.pushButton_config.setText(_translate("mainWindow", "Generate config file"))
        self.pushButton_Download_SC.setText(_translate("mainWindow", ""))
        self.pushButton_Download_Photons.setText(_translate("mainWindow", ""))
        self.pushButton_Download_diffuse.setText(_translate("mainWindow", ""))
        self.label_Log.setText(_translate("mainWindow", "Log:"))
        self.groupBox_Science.setTitle(_translate("mainWindow", "Analyses selection:"))
        self.radioButton_disk.setText(_translate("mainWindow", "Disk"))
        self.label_N_time_bins_LC.setText(_translate("mainWindow", "N° of time bins"))
        self.checkBox_extension.setText(_translate("mainWindow", "Extension:"))
        self.checkBox_adaptive_binning.setText(_translate("mainWindow", "Adaptive-binning"))
        self.label_N_cores_LC.setText(_translate("mainWindow", "N° of cores"))
        self.label_TS_threshold.setText(_translate("mainWindow", "TS threshold"))
        self.label_N_iter.setText(_translate("mainWindow", "N° iterations"))
        self.label_N_energy_bins.setText(_translate("mainWindow", "N⁰ of energy bins"))
        self.label_redshift.setText(_translate("mainWindow", "Redshift"))
        self.label_MCMC.setText(_translate("mainWindow", "MCMC"))
        self.checkBox_LC.setText(_translate("mainWindow", "Light curve:"))
        self.label_Extension_max_size.setText(_translate("mainWindow", "Max. radius"))
        self.checkBox_SED.setText(_translate("mainWindow", "SED:"))
        self.checkBox_use_local_index.setText(_translate("mainWindow", "Use local index"))
        self.radioButton_2D_Gauss.setText(_translate("mainWindow", "2D-Gauss"))
        self.checkBox_residual_TSmap.setText(_translate("mainWindow", "Residuals TS map"))
        self.checkBox_TSmap.setText(_translate("mainWindow", "TS and excess maps:"))
        self.checkBox_relocalize.setText(_translate("mainWindow", "Relocalize"))
        self.label_photon_index_TS.setText(_translate("mainWindow", "Photon index"))
        self.groupBox_fit_finetune.setTitle(_translate("mainWindow", "Fine-tuning the fit:"))
        self.checkBox_delete_sources.setText(_translate("mainWindow", "Delete sources:"))
        self.comboBox_change_model.setAccessibleName(_translate("mainWindow", "4FGL"))
        self.comboBox_change_model.setAccessibleDescription(_translate("mainWindow", "4FGL"))
        self.comboBox_change_model.setItemText(0, _translate("mainWindow", "Select..."))
        self.comboBox_change_model.setItemText(1, _translate("mainWindow", "Power-law"))
        self.comboBox_change_model.setItemText(2, _translate("mainWindow", "Power-law2"))
        self.comboBox_change_model.setItemText(3, _translate("mainWindow", "LogPar"))
        self.comboBox_change_model.setItemText(4, _translate("mainWindow", "PLEC"))
        self.comboBox_change_model.setItemText(5, _translate("mainWindow", "PLEC2"))
        self.comboBox_change_model.setItemText(6, _translate("mainWindow", "PLEC3"))
        self.comboBox_change_model.setItemText(7, _translate("mainWindow", "PLEC4"))
        self.comboBox_change_model.setItemText(8, _translate("mainWindow", "BPL"))
        self.comboBox_change_model.setItemText(9, _translate("mainWindow", "ExpCutOff-EBL"))
        self.comboBox_minimizer.setAccessibleName(_translate("mainWindow", "4FGL"))
        self.comboBox_minimizer.setAccessibleDescription(_translate("mainWindow", "4FGL"))
        self.comboBox_minimizer.setItemText(0, _translate("mainWindow", "NEWMINUIT"))
        self.comboBox_minimizer.setItemText(1, _translate("mainWindow", "MINUIT"))
        self.comboBox_minimizer.setItemText(2, _translate("mainWindow", "DRMNGB"))
        self.comboBox_minimizer.setItemText(3, _translate("mainWindow", "DRMNFB"))
        self.comboBox_minimizer.setItemText(4, _translate("mainWindow", "LBFGS"))
        self.comboBox_redshift.setAccessibleName(_translate("mainWindow", "redshift_model"))
        self.comboBox_redshift.setAccessibleDescription(_translate("mainWindow", "redshift_model"))
        self.comboBox_redshift.setItemText(0, _translate("mainWindow", "Franceschini et al. (2008)"))
        self.comboBox_redshift.setItemText(1, _translate("mainWindow", "Finke et al. (2010)"))
        self.comboBox_redshift.setItemText(2, _translate("mainWindow", "Dominguez et al. (2011)"))
        self.comboBox_redshift.setItemText(3, _translate("mainWindow", "Franceschini & Rodighiero (2017)"))
        self.comboBox_redshift.setItemText(4, _translate("mainWindow", "Saldana-Lopez et al. (2021)"))
        self.comboBox_MCMC.setAccessibleName(_translate("mainWindow", "MCMC_model"))
        self.comboBox_MCMC.setAccessibleDescription(_translate("mainWindow", "MCMC_model"))
        self.comboBox_MCMC.setItemText(0, _translate("mainWindow", "PowerLaw"))
        self.comboBox_MCMC.setItemText(1, _translate("mainWindow", "LogPar"))
        self.comboBox_MCMC.setItemText(2, _translate("mainWindow", "LogPar2"))
        self.comboBox_MCMC.setItemText(3, _translate("mainWindow", "PLEC"))
        self.comboBox_MCMC.setItemText(4, _translate("mainWindow", "PLEC_bfix"))
        self.white_box_list_of_sources_to_delete.setText(_translate("mainWindow", ""))
        self.checkBox_diagnostic_plots.setText(_translate("mainWindow", "Diagnostic plots"))
        self.label_minimum_separation.setText(_translate("mainWindow", "Minimum separation (⁰)"))
        self.checkBox_find_extra_sources.setText(_translate("mainWindow", "Find extra sources in the ROI:"))
        self.label_min_significance.setText(_translate("mainWindow", "Minimum significance"))
        self.checkBox_change_model.setText(_translate("mainWindow", "Change model:"))
        self.checkBox_minimizer.setText(_translate("mainWindow", "Change optimizer:"))
        self.radioButton_free_source_radius.setText(_translate("mainWindow", "Defaut"))
        self.radioButton_free_source_radius_customized.setText(_translate("mainWindow", "Customized"))
        self.label_free_src_radius.setText(_translate("mainWindow", "Free source radius:"))
        self.label_radius.setText(_translate("mainWindow", "Radius (⁰)"))
        self.checkBox_only_norm.setText(_translate("mainWindow", "Only norm."))
        self.checkBox_freeze_gal.setText(_translate("mainWindow", "Freeze Gal."))
        self.checkBox_freeze_iso.setText(_translate("mainWindow", "Freeze Iso."))
        self.checkBox_freeze_spec_shape.setText(_translate("mainWindow", "Freeze shape targ."))
        self.label_output_format.setText(_translate("mainWindow", "Output format:"))
        self.comboBox_output_format.setAccessibleName(_translate("mainWindow", "4FGL"))
        self.comboBox_output_format.setAccessibleDescription(_translate("mainWindow", "4FGL"))
        self.comboBox_output_format.setItemText(0, _translate("mainWindow", "pdf"))
        self.comboBox_output_format.setItemText(1, _translate("mainWindow", "png"))
        self.large_white_box_Log.setPlainText(_translate("mainWindow", "\n"))

        self.toolButton_External_ltcube.setText(_translate("mainWindow", "..."))
        self.toolButton_dir_diffuse.setText(_translate("mainWindow", "..."))
        self.label_dir_diffuse.setText(_translate("mainWindow", "Diffuse emission directory:"))
        self.checkBox_External_ltcube.setText(_translate("mainWindow", "Use external ltcube:"))
        self.checkBox_highest_resolution.setText(_translate("mainWindow", "Improve resolution."))
        self.checkBox_high_sensitivity.setText(_translate("mainWindow", "Improve sensitivity."))
        self.dateTimeEdit_2.setDisplayFormat(_translate("mainWindow", "dd/MM/yyyy HH:mm:ss"))
        self.dateTimeEdit.setDisplayFormat(_translate("mainWindow", "dd/MM/yyyy HH:mm:ss"))
        self.label_start_time.setText(_translate("mainWindow", "Start time:"))
        self.label_end_time.setText(_translate("mainWindow", "Stop time:"))
        self.label_limits.setText(_translate("mainWindow", "Advanced options:"))
        self.toolButton_output_dir.setText(_translate("mainWindow", "..."))
        self.white_box_output_dir.setText(_translate("mainWindow", "./Output"))
        self.white_box_redshift.setText(_translate("mainWindow", "0.0"))
        self.white_box_TS_threshold.setText(_translate("mainWindow", "50.0"))
        self.label_dir_photons.setText(_translate("mainWindow", "Photon files directory:"))
        self.label_dir_spacecraft.setText(_translate("mainWindow", "Spacecraft file:"))
        self.label_output_dir.setText(_translate("mainWindow", "Output directory:"))
        self.toolButton_spacecraft.setAccessibleDescription(_translate("mainWindow", "spacecraft mission file"))
        self.toolButton_spacecraft.setText(_translate("mainWindow", "..."))
        self.toolButton_photons.setText(_translate("mainWindow", "..."))
        self.toolButton_Custom.setText(_translate("mainWindow", "..."))
        self.toolButton_VHE.setText((_translate("mainWindow", "...")))
        self.label_energy.setText(_translate("mainWindow", "<html><head/><body><p>E<span style=\" vertical-align:sub;\">min</span>, E<span style=\" vertical-align:sub;\">max</span> (MeV):</p></body></html>"))
        self.white_box_energy.setText(_translate("mainWindow", "100, 300000"))
        self.label_Catalog.setText(_translate("mainWindow", "Catalog:"))
        self.label_RAandDec.setText(_translate("mainWindow", "RA, Dec (⁰):"))
        self.comboBox_Catalog.setAccessibleName(_translate("mainWindow", "4FGL"))
        self.comboBox_Catalog.setAccessibleDescription(_translate("mainWindow", "4FGL"))
        self.comboBox_Catalog.setItemText(0, _translate("mainWindow", "4FGL-DR3"))
        self.comboBox_Catalog.setItemText(1, _translate("mainWindow", "4FGL-DR2"))
        self.comboBox_Catalog.setItemText(2, _translate("mainWindow", "4FGL"))
        self.comboBox_Catalog.setItemText(3, _translate("mainWindow", "3FGL"))
        self.white_box_config_file.setText(_translate("mainWindow", "Configuration file (yaml)"))
        self.white_box_target_name.setText(_translate("mainWindow", "Target name"))
        self.radioButton_Standard.setText(_translate("mainWindow", "Standard"))
        self.radioButton_Custom.setText(_translate("mainWindow", "Custom"))
        self.label_configuration_file.setText(_translate("mainWindow", "Configuration file:"))
        self.label_is_it_cataloged.setText(_translate("mainWindow", "Target cataloged/name?"))
        self.comboBox_is_it_cataloged.setAccessibleName(_translate("mainWindow", "4FGL"))
        self.comboBox_is_it_cataloged.setAccessibleDescription(_translate("mainWindow", "4FGL"))
        self.comboBox_is_it_cataloged.setItemText(0, _translate("mainWindow", "Yes"))
        self.comboBox_is_it_cataloged.setItemText(1, _translate("mainWindow", "No"))
        self.menuTutorial.setTitle(_translate("mainWindow", "Menu"))
        self.menuCredits.setTitle(_translate("mainWindow", "Credits"))
        self.actionOpen.setText(_translate("mainWindow", "Open..."))
        self.actionSave.setText(_translate("mainWindow", "Save"))
        self.actionSave.setIconText(_translate("mainWindow", "Save"))
        self.actionSave.setShortcut(_translate("mainWindow", "Ctrl+S"))
        self.actionCopy.setText(_translate("mainWindow", "Copy"))
        self.actionCopy.setShortcut(_translate("mainWindow", "Ctrl+C"))
        self.actionPaste.setText(_translate("mainWindow", "Paste"))
        self.actionPaste.setStatusTip(_translate("mainWindow", "Paste a file"))
        self.actionPaste.setShortcut(_translate("mainWindow", "Ctrl+V"))
        self.actionOpen_Tutorial.setText(_translate("mainWindow", "Open Tutorial"))
        self.actionSee_credits.setText(_translate("mainWindow", "See credits"))
        self.actionLoad_state.setText(_translate("mainWindow", "Load GUI state"))

    def reportProgress(self, n):

        """
        This function is called only from a secondary thread where the class Worker() is running.
        The goal of this function is to print updates about the Fermi-LAT analysis in the easyfermi Log box.

        Parameters
        ----------
        self: instance of the class Ui_mainWindow
            This parameter contains all the variables read from the Graphical interface.
        n: int
            An integer number emitted by pyqtsignal from the Worker() class.

        Returns
        -------
        """

        #Making Fermi logo transparent:
        new_pix = QtGui.QPixmap(str(libpath/"fermi.png"))
        new_pix.fill(QtCore.Qt.transparent)
        painter = QtGui.QPainter(new_pix)
        painter.setOpacity(0.2)
        painter.drawPixmap(QtCore.QPoint(), QtGui.QPixmap(str(libpath/"fermi.png")))
        painter.end()
        self.picture.setPixmap(new_pix)
        self.picture.setGeometry(QtCore.QRect(760, 230, 231, 231))

        if n == -3:
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Downloading Diffuse model files...\n")
        
        if n == -2:
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Downloading Photon files...\n")

        if n == -1:
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Downloading spacecraft file...\n")

        if n == 0:
            self.progressBar.setProperty("value", 5)

            if self.white_box_list_of_sources_to_delete.text() != '':
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Deleting source(s): "+self.white_box_list_of_sources_to_delete.text()+".\n")
                
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"------------------------\n- Running the setup.\n")

            try:
                # Here we check if the computer is connected to power
                if psutil.sensors_battery()[2] is False:
                    self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- WARNING: Your computer is running on battery. Plugging it to power will make the analysis significantly faster.\n")
            except:
                pass

            if (self.IsThereLtcube is None) & (self.IsThereLtcube3==0) & (self.checkBox_External_ltcube.isChecked() is False):
                if self.checkBox_high_sensitivity.isChecked():
                    if self.Emin < 500 and self.Emax > 1000:
                        multiplication_factor = 2
                    elif self.Emin < 500 and 500 <= self.Emax < 1000:
                        multiplication_factor = 2
                    elif 500 <= self.Emin < 1000 and self.Emax > 1000:
                        multiplication_factor = 2
                    else:
                        multiplication_factor = 1
                else:
                    multiplication_factor = 1
        
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- It will take about ~"+str(multiplication_factor*int(self.Time_intervMJD*206/(30.0*60)))+" min to run the ltcube.\n")
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- (Don't worry, this really takes some time...)\n")
                try:
                    if psutil.sensors_battery()[2] is False:
                        self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- WARNING: Since your computer is running on battery power, this process can actually take much longer.\n")
                except:
                    pass

            else:
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Using precomputed ltcube.\n")
            
            
        if n == 1:
            self.progressBar.setProperty("value", 25)
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Setup finished.\n")
            os.system('ls '+self.OutputDir+'ltcube_*.fits > '+self.OutputDir+'ltcube_list.txt')
            list_of_photon_files = glob.glob(self.white_box_output_dir.text()+"/ft1*.fits")
            max_photon_energy= 0
            for photon_file in list_of_photon_files:
                max_energy_in_the_file = pyfits.open(photon_file)[1].data["ENERGY"].max()
                if max_energy_in_the_file > max_photon_energy:
                    max_photon_energy = max_energy_in_the_file
            
            highest_energy_photon_RoI = round(max_photon_energy/1000,2)
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Highest energy photon in the RoI: "+str(highest_energy_photon_RoI)+" GeV.\n")
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Optimizing the RoI...\n")
        
        if n == 2:
            self.progressBar.setProperty("value", 30)
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Optimization is done.\n")
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+f"- Target spectral type: {self.gta.roi.sources[0]['SpectrumType']}.\n")

            if self.freeradiusalert != 'ok':
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+self.freeradiusalert+"\n")

            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Performing the fit...\n")
        
        if n == 3:
            self.progressBar.setProperty("value", 40)
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+self.fitquality+"\n")
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Main results (flux, spectral index, TS, etc) saved in Target_results.txt\n")
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Computing distance from the Sun...\n")

        if n == 4:
            self.progressBar.setProperty("value", 50)

            if self.checkBox_diagnostic_plots.isChecked():
                rounded_separation = round(self.Solar_separation.min(),2)
                if rounded_separation < 15:
                    self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+
                                                          "- Minimum separation between the target and the Sun: "+str(rounded_separation)+
                                                          "°.\nPlease be aware that the Sun can affect your observations.\nYou can check the time ranges when the Sun is nearby in the diagnostic plots.")
                else:
                    self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Minimum separation between the target and the Sun: "+
                                                          str(rounded_separation)+"°.\n")
                
                if len(self.Solar_separation) < 10:
                    self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+
                                                          '- Time window is too short for computing the target-Sun separation plot. Minimum window required is 4 days.\n')
            else:
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+self.fitquality+"\n")
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+
                                                      "- Main results (flux, spectral index, TS, etc) saved in Target_results.txt\n")
                    
            if self.checkBox_relocalize.isChecked():
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Relocalizing target...\n")

        if n == 5:
            self.progressBar.setProperty("value", 55)
            if self.checkBox_TSmap.isChecked():
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Computing TS maps...\n")

        if n == 6:
            self.progressBar.setProperty("value", 60)
            if self.checkBox_extension.isChecked():
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Looking for extended emission...\n")

        if n == 7:
            self.progressBar.setProperty("value", 70)
            if self.checkBox_SED.isChecked():
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Computing SED.\n")

        if n == 8:
            self.progressBar.setProperty("value", 75)
            if self.checkBox_SED.isChecked():
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Starting MCMC...\n")

        if n == 9:
            self.progressBar.setProperty("value", 80)
            try:
                if self.allow_MCMC is False:
                    self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+
                                                      "- MCMC not allowed. We require at least 3 LAT data points with TS > 9 in the SED to proceed.\n")
            except:
                pass

            if self.include_VHE is False:
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- VHE data not available.\n")

            if self.checkBox_LC.isChecked():
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Computing light curve.\n")
        
        if n == 10:
            self.progressBar.setProperty("value", 90)

            if self.Compute_LC is False:
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+f"- A light curve with {self.spinBox_LC_N_time_bins.value()} was already found. Skipping new light curve computation.\n")

            if self.adaptive:
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Computing adaptive-binning light curve.\n")
                
        if n == 11:
            self.progressBar.setProperty("value", 99)
            if self.checkBox_relocalize.isChecked():
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- New position: RA = "+str(round(self.locRA,3))+
                                                      ", Dec = "+str(round(self.locDec,3))+", r_95 = "+str(round(self.locr95,3))+"\n")
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Localization results saved in "+self.sourcename+"_loc.fits\n")        
            if self.checkBox_TSmap.isChecked():
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- TS maps saved as figures and fits files.\n")
            if self.checkBox_SED.isChecked():
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- SED data saved in "+self.sourcename+"_sed.fits\n")
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Preliminar SED shown in figure Quickplot_SED\n")
                try:
                    self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+self.redshift_error)
                    del self.redshift_error
                except:
                    pass
            if self.checkBox_extension.isChecked():
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Extension data saved in "+self.sourcename+"_ext.fits\n")
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Extension is shown in figure Quickplot_extension\n")
            if self.checkBox_LC.isChecked():
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- LC saved in file "+self.sourcename+"_lightcurve.fits\n")
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- LC is shown in figure Quickplot_LC\n")
                        
        
        
            
    def runLongTask(self):

        """
        This function connects the easyfermi window with the class Worker()
        and sends the heavy job to another thread, then leaving the easyfermi
        window active even when the analysis is running.
        """
        
        can_we_go = self.check_for_erros()
        
        if can_we_go:
            self.thread = QtCore.QThread()
            self.worker = Worker()
            self.worker.moveToThread(self.thread)
        
            # Standard or custom analysis?
            Standard = self.radioButton_Standard.isChecked()
            Custom = self.radioButton_Custom.isChecked()

            # Step 5: Connect signals and slots
            self.thread.started.connect(self.worker.run_gtsetup)
            self.worker.finished.connect(self.thread.quit)
            self.worker.finished.connect(self.worker.deleteLater)
            self.thread.finished.connect(self.thread.deleteLater)
            self.worker.starting.connect(self.readytogo)
            self.worker.progress.connect(self.reportProgress)
        
            # Step 6: Start the thread
            self.thread.start()
        
            # Final reset
            self.thread.finished.connect(  lambda: self.pushButton_Go.setEnabled(True)  )
            self.thread.finished.connect(  lambda: self.radioButton_Standard.setEnabled(True)  )
            self.thread.finished.connect(  lambda: self.radioButton_Custom.setEnabled(True)  )
            self.thread.finished.connect(  lambda: self.radioButton_Standard.setChecked(Standard)  )
            self.thread.finished.connect(  lambda: self.radioButton_Custom.setChecked(Custom)  )
            self.thread.finished.connect(  self.activate  )
            
            _translate = QtCore.QCoreApplication.translate
        
            self.thread.finished.connect(  lambda: self.pushButton_Go.setText(_translate("MainWindow", "Go!"))  )
        
            self.thread.finished.connect(  lambda: self.progressBar.setProperty("value", 100)  )
            
            self.thread.finished.connect(  lambda: self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Saving GUI status...\n")  )
            
            self.thread.finished.connect(  self.save_GUIstate   )
            
            self.thread.finished.connect(  lambda: self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Process finished!\n")  )
            
        else:
            self.popup_go()
    
    
    def save_GUIstate(self):

        """
        This function is called only at the end of the analysis.
        It saves the current state of the GUI in the file GUI_status.npy
        in the output directory. The user can latter use this file to recover the
        state of the GUI.
        """

        state = {}

        state["Standard"] = self.radioButton_Standard.isChecked()
        state["Coords"] = self.white_box_RAandDec.text()
        state["Energ"] = self.white_box_energy.text()
        state["date"] = self.dateTimeEdit.text()
        state["date2"] = self.dateTimeEdit_2.text()
        state["spacecraft"] = self.white_box_spacecraft_file.text()
        state["diffuse"] = self.white_box_Diffuse_dir.text()
        state["dir_photon"] = self.white_box_photon_dir.text()
        state["Use_external_ltcube"] = self.checkBox_External_ltcube.isChecked()
        state["external_ltcube"] = self.white_box_External_ltcube.text()
        state["catalog"] = self.comboBox_Catalog.currentText()
        state["cataloged"] = self.comboBox_is_it_cataloged.currentText()
        state["configfile"] = self.white_box_config_file.text()
        state["nickname"] = self.white_box_target_name.text()
        state["High_resolution"] = self.checkBox_highest_resolution.isChecked()
        state["High_sensitivity"] = self.checkBox_high_sensitivity.isChecked()

        state["change_minimizer"] = self.checkBox_minimizer.isChecked()
        state["which_minimizer"] = self.comboBox_minimizer.currentText()
        state["change_model"] = self.checkBox_change_model.isChecked()
        state["which_model"] = self.comboBox_change_model.currentText()
        state["delete_sources"] = self.checkBox_delete_sources.isChecked()
        state["which_sources_deleted"] = self.white_box_list_of_sources_to_delete.text()
        state["Free_radius_standard"] = self.radioButton_free_source_radius.isChecked()
        state["Free_radius_custom"] = self.radioButton_free_source_radius_customized.isChecked()
        state["free_radius"] = self.white_box_radius.text()
        state["Only_norm"] = self.checkBox_only_norm.isChecked()
        state["Freeze_Gal"] = self.checkBox_freeze_gal.isChecked()
        state["Freeze_Iso"] = self.checkBox_freeze_iso.isChecked()
        state["Freeze_targ_shape"] = self.checkBox_freeze_spec_shape.isChecked()
        state["find_sources"] = self.checkBox_find_extra_sources.isChecked()
        state["min_sig"] = self.doubleSpinBox_min_significance.value()
        state["min_sep"] = self.doubleSpinBox_min_separation.value()
        state["diagnostic"] = self.checkBox_diagnostic_plots.isChecked()
        state["output_format"] = self.comboBox_output_format.currentText()

        state["LC"] = self.checkBox_LC.isChecked()
        state["LC_Nbins"] = self.spinBox_LC_N_time_bins.value()
        state["LC_Ncores"] = self.spinBox_N_cores_LC.value()
        state["adaptive_binning"] = self.checkBox_adaptive_binning.isChecked()
        state["TS_threshold"] = self.white_box_TS_threshold.text()
        state["N_iter"] = self.spinBox_N_iter.value()
        state["SED"] = self.checkBox_SED.isChecked()
        state["SED_Nbins"] = self.spinBox_SED_N_energy_bins.value()
        state["VHE"] = self.white_box_VHE.text()
        state["use_local_index"] = self.checkBox_use_local_index.isChecked()
        state["which_MCMC_model"] = self.comboBox_MCMC.currentText()
        state["redshift_value"] = self.white_box_redshift.text()
        state["EBL_model"] = self.comboBox_redshift.currentText()
        state["extension"] = self.checkBox_extension.isChecked()
        state["Disk"] = self.radioButton_disk.isChecked()
        state["Gauss2D"] = self.radioButton_2D_Gauss.isChecked()
        state["max_size"] = self.doubleSpinBox_extension_max_size.value()
        state["reloc"] = self.checkBox_relocalize.isChecked()
        state["TS_map"] = self.checkBox_TSmap.isChecked()
        state["test_source_index"] = self.doubleSpinBox_Photon_index_TS.value()
        state["remove_targ_from_model"] = self.checkBox_residual_TSmap.isChecked()        
        state["output"] = self.white_box_output_dir.text()
        
        
        
        with open(self.OutputDir+"./GUI_status.yaml", 'w') as yaml_file:
            yaml_file.write( yaml.dump(state, default_flow_style=False))            
        
        
            
    def load_GUIstate(self):
    
        """
        This function allows the user to load all the entries previously used in an analysis.
        It can be called fro mthe main menu -> Load GUI state.
        """
        
        fname = QFileDialog.getOpenFileName(self, 'Open file', '', '(*.yaml)')
        stream = open(fname[0], 'r')
        keys = yaml.load(stream,Loader)
        
        if keys["Standard"] is True:
            self.radioButton_Standard.setChecked(True)
            self.radioButton_Custom.setChecked(False)
            self.white_box_RAandDec.setText(keys["Coords"])
            self.white_box_energy.setText(keys["Energ"])
            self.white_box_spacecraft_file.setText(keys["spacecraft"])
            self.white_box_Diffuse_dir.setText(keys["diffuse"])
            self.white_box_photon_dir.setText(keys["dir_photon"])
            self.white_box_target_name.setText(keys["nickname"])
            self.white_box_External_ltcube.setText(keys["external_ltcube"])
            self.checkBox_External_ltcube.setChecked(keys["Use_external_ltcube"])

            date1 = keys["date"].split(' ')
            date = np.asarray(date1[0].split('/')).astype(int)
            time = np.asarray(date1[1].split(':')).astype(int)
            self.dateTimeEdit.setDateTime(QtCore.QDateTime(QtCore.QDate(date[2], date[1], date[0]), QtCore.QTime(time[0], time[1], time[2])))
            date2 = keys["date2"].split(' ')
            date = np.asarray(date2[0].split('/')).astype(int)
            time = np.asarray(date2[1].split(':')).astype(int)
            self.dateTimeEdit_2.setDateTime(QtCore.QDateTime(QtCore.QDate(date[2], date[1], date[0]), QtCore.QTime(time[0], time[1], time[2])))

            catalog = keys["catalog"]
            cataloged = keys["cataloged"]
            aux = self.comboBox_Catalog.findText(catalog, QtCore.Qt.MatchFixedString)
            self.comboBox_Catalog.setCurrentIndex(aux)
            aux = self.comboBox_is_it_cataloged.findText(cataloged, QtCore.Qt.MatchFixedString)
            self.comboBox_is_it_cataloged.setCurrentIndex(aux)
            
        else:
            self.radioButton_Standard.setChecked(False)
            self.radioButton_Custom.setChecked(True)
            self.white_box_config_file.setText(keys["configfile"])
            
        
        # Advanced options:
        self.checkBox_highest_resolution.setChecked(keys["High_resolution"])
        self.checkBox_high_sensitivity.setChecked(keys["High_sensitivity"])
            
        # Change model:    
        self.checkBox_change_model.setChecked(keys["change_model"])
        aux = self.comboBox_change_model.findText(keys["which_model"], QtCore.Qt.MatchFixedString)
        self.comboBox_change_model.setCurrentIndex(aux)

        # Change minimizer:
        self.checkBox_minimizer.setChecked(keys["change_minimizer"])
        aux = self.comboBox_minimizer.findText(keys["which_minimizer"], QtCore.Qt.MatchFixedString)
        self.comboBox_minimizer.setCurrentIndex(aux)
        
        # Delete sources:    
        self.checkBox_delete_sources.setChecked(keys["delete_sources"])
        self.white_box_list_of_sources_to_delete.setText(keys["which_sources_deleted"])
        
        # Free source radius:
        self.radioButton_free_source_radius.setChecked(keys["Free_radius_standard"])
        self.radioButton_free_source_radius_customized.setChecked(keys["Free_radius_custom"])
        self.white_box_radius.setText(keys["free_radius"])
        self.checkBox_only_norm.setChecked(keys["Only_norm"])
        self.checkBox_freeze_gal.setChecked(keys["Freeze_Gal"])
        self.checkBox_freeze_iso.setChecked(keys["Freeze_Iso"])
        self.checkBox_freeze_spec_shape.setChecked(keys["Freeze_targ_shape"])
        
        # Find sources:
        self.checkBox_find_extra_sources.setChecked(keys["find_sources"])
        self.doubleSpinBox_min_significance.setProperty("value", float(keys["min_sig"]))
        self.doubleSpinBox_min_separation.setProperty("value", float(keys["min_sep"]))

        # Diagnostic plots:
        self.checkBox_diagnostic_plots.setChecked(keys["diagnostic"])
            
        # Output format:
        aux = self.comboBox_output_format.findText(keys["output_format"], QtCore.Qt.MatchFixedString)
        self.comboBox_output_format.setCurrentIndex(aux)
        
        # LC:
        self.checkBox_LC.setChecked(keys["LC"])
        self.checkBox_adaptive_binning.setChecked(keys["adaptive_binning"])
        self.spinBox_LC_N_time_bins.setProperty("value", int(keys["LC_Nbins"]))     
        self.spinBox_N_cores_LC.setProperty("value", int(keys["LC_Ncores"]))  
        self.spinBox_N_iter.setProperty("value", int(keys["N_iter"]))  
        self.white_box_TS_threshold.setText(keys["TS_threshold"])
        
        # SED:
        self.checkBox_SED.setChecked(keys["SED"])
        self.spinBox_SED_N_energy_bins.setProperty("value", int(keys["SED_Nbins"]))
        self.checkBox_use_local_index.setChecked(keys["use_local_index"])
        self.white_box_redshift.setText(keys["redshift_value"])
        aux = self.comboBox_redshift.findText(keys["EBL_model"], QtCore.Qt.MatchFixedString)
        self.comboBox_redshift.setCurrentIndex(aux)
        aux = self.comboBox_MCMC.findText(keys["which_MCMC_model"], QtCore.Qt.MatchFixedString)
        self.comboBox_MCMC.setCurrentIndex(aux)
        self.white_box_VHE.setText(keys["VHE"])

        # Extension:
        self.checkBox_extension.setChecked(keys["extension"])
        self.radioButton_disk.setChecked(keys["Disk"])
        self.radioButton_2D_Gauss.setChecked(keys["Gauss2D"])
        self.doubleSpinBox_extension_max_size.setProperty("value", float(keys["max_size"])) 
        
        # Relocalize:
        self.checkBox_relocalize.setChecked(keys["reloc"])
        
        # TSmap:
        self.checkBox_TSmap.setChecked(keys["TS_map"])
        self.doubleSpinBox_Photon_index_TS.setProperty("value", float(keys["test_source_index"]))     
        self.checkBox_residual_TSmap.setChecked(keys["remove_targ_from_model"])
        
        # Output dir:
        self.white_box_output_dir.setText(keys["output"])
        
        # Activate/deactivate boxes according to the loaded GUI_satatus file
        self.activate()

    
    def check_for_erros(self,jump_paths=False):

        """
        Here we check if all input parameters given by the user are ok.
        
        Parameters
        ----------
        jump_paths (optional): boolean
                               This is used only when downloading the data and, if True, will skip the checking of the data paths.

        Returns
        -------
            A boolean answer set as "True" if all inputs are ok, or set as "False" if any problem is detected.

        """

        if self.radioButton_Standard.isChecked():
            self.reportProgress(n=None)           
            check = 0
            Coords = self.white_box_RAandDec.text().split(',')
            Energ = self.white_box_energy.text().split(',')
            date = self.dateTimeEdit.text()
            date2 = self.dateTimeEdit_2.text()
            times = [str(date[6:10])+'-'+str(date[3:5])+'-'+str(date[:2])+'T'+str(date[11:19]),    str(date2[6:10])+'-'+str(date2[3:5])+'-'+str(date2[:2])+'T'+str(date2[11:19])              ]
            t = Time(times)
            t0 = t.mjd[0]
            t1 = t.mjd[1]
            
            
            if jump_paths:
                pass
            else:
                if (len(self.white_box_spacecraft_file.text()) > 0) & (len(self.white_box_Diffuse_dir.text()) > 0) & (len(self.white_box_photon_dir.text()) > 0):
                    pass
                else:
                    check = check + 1
                    self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Invalid path. Please double check the data paths above.\n")
            
                if self.checkBox_External_ltcube.isChecked():
                    if len(self.white_box_External_ltcube.text()) < 1:
                        check = check + 1
                        self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Invalid path for external ltcube.\n")
                
                if (len(self.white_box_spacecraft_file.text().split(' ')) > 1) or (len(self.white_box_Diffuse_dir.text().split(' ')) > 1) or (len(self.white_box_photon_dir.text().split(' ')) > 1):
                    check = check + 1
                    self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- The data paths cannot contain blank spaces. Please fix this before proceeding.\n")
            
            try:
                if len(Coords) == 1:
                    self.recover_coords_from_name = Simbad.query_object(Coords[0])
                    try:
                        if len(self.recover_coords_from_name) == 1:
                            pass
                    except:
                        print('Target name not found in Simbad.')
                        self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Target name not found in Simbad.\n")
                        check = check + 1

                else:
                    Coords = np.asarray(Coords).astype(float)
                    if (Coords[0] <= 360) & (Coords[0] >= 0) & (Coords[1] >= -90) & (Coords[1] <= 90):
                        pass
                    else:
                        print('Problems with the coordinates.')
                        check = check + 1
                        self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Invalid coordinates. Try 0 <= RA <= 360 and -90 <= Dec <= 90.\n")

                Energ = np.asarray(Energ).astype(float)
                if (Energ[0] < Energ[1]) & (Energ[0] > 20) & (Energ[1] <= 10000000):
                    pass
                else:
                    print("Problems in the energy range.")
                    check = check + 1
                    self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Invalid energy range. Try values from 20 to 10000000 MeV.\n")
                if t1 <= t0:
                    check = check + 1
                    self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Stop time must be set to a date after start time.\n")
            
                list_of_photon_files = glob.glob(self.white_box_photon_dir.text()+"/*.fits")
                max_photon_energy= 0
                for photon_file in list_of_photon_files:
                    max_energy_in_the_file = pyfits.open(photon_file)[1].data["ENERGY"].max()
                    if max_energy_in_the_file > max_photon_energy:
                        max_photon_energy = max_energy_in_the_file
                
                if max_photon_energy > 0 and Energ[1] > 1.1*max_photon_energy:
                    check = check + 1
                    self.max_energy_wanning = f"The maximum energy is set to {Energ[1]} MeV, however, the highest energy photon in the dataset has only {max_photon_energy} MeV.\n\nPlease set Emax <= {max_photon_energy} MeV before proceeding."

                if check == 0:
                    answer = True
                else:
                    answer = False
            except:
                answer = False
            
            
        
        else:
            check = 0
            Config = self.white_box_config_file.text().split('.')
            if Config[-1] == 'yaml':
                pass
            else:
                check = check + 1
                self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"- Please provide a yaml configuration file.\n")
            
            
            if check == 0:
                answer = True
            else:
                answer = False
        return answer
    

    def readytogo(self):

        """
        Here we update the "Go!" button when the analysis is started
        and set the progress bar to zero.
        """
        _translate = QtCore.QCoreApplication.translate
        self.pushButton_Go.setEnabled(False)
        self.pushButton_Go.setText(_translate("MainWindow", "Running..."))
        self.progressBar.setProperty("value", 0)

        self.radioButton_Standard.setEnabled(False)
        self.radioButton_Custom.setEnabled(False)
        self.radioButton_Standard.setChecked(False)
        self.radioButton_Custom.setChecked(False)
        self.activate()
        
        
        
        
        
    def setFermipy(self):
        
        """
        Here we setup the GTAnalysis class using the input information given by the user
        in the GUI, or in the customized configuration file.
        """
        
        if self.radioButton_Custom.isChecked(): 
            self.OutputDir = self.white_box_output_dir.text()+'/'
            self.gta = GTAnalysis(self.white_box_config_file.text(),logging={'verbosity': 3})
            self.roiwidth = self.gta.config.get('binning').get('roiwidth')
            self.Emin = self.gta.config.get('selection').get('emin')
            self.Emax = self.gta.config.get('selection').get('emax')
            self.RA = self.gta.config.get('selection').get('ra')
            self.Dec = self.gta.config.get('selection').get('dec')
            self.Time_intervMJD = (self.gta.config.get('selection').get('tmax') - self.gta.config.get('selection').get('tmin'))/86400
            self.spacecraft_file = self.gta.config.get('data').get('scfile')
            self.tmin = self.gta.config.get('selection').get('tmin')
            self.tmax = self.gta.config.get('selection').get('tmax')
            
        else:
            self.generateConfig()
            self.gta = GTAnalysis(self.OutputDir+'config.yaml',logging={'verbosity': 3})
        
        
        
        #Get target name:
        for n,s in enumerate(self.gta.roi.sources):
            if n == 0:
                self.sourcename = s.name
        
        
        #Checking for ltcube:
        self.IsThereLtcube = self.gta.config.get('data').get('ltcube')
        self.IsThereLtcube2 = glob.glob(self.OutputDir+'*.fits')  
        self.IsThereLtcube3 = 0
        if self.OutputDir+'ltcube_00.fits' in self.IsThereLtcube2:
            self.IsThereLtcube3 = 1


    def click_to_generateConfig(self):

        """
        Fucntion that calls the generateConfig() 
        """

        self.reportProgress(n=None)
        can_we_go = self.generateConfig()
        if can_we_go:
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+f"- The file config.yaml was saved in {self.OutputDir}\n")

        
    def generateConfig(self):

        """
        This function generates the yaml configuration file required for the analysis.
        
        Returns
        -------
            File named "config.yaml" saved in the chosen output directory and containing information on RA, Dec, energy range, time range etc.

        """
        
        output_dir = Path(self.white_box_output_dir.text())
        self.OutputDir = str(output_dir.parent.resolve())+'/'+output_dir.name+'/'

        if not os.path.exists(self.OutputDir):
            os.system(f"mkdir {self.OutputDir}")

        can_we_go = self.check_for_erros()

        if can_we_go:
            date = self.dateTimeEdit.text()
            date2 = self.dateTimeEdit_2.text()
            timeStart = ['2008-08-04T15:43:36']
            tStart = Time(timeStart)
            self.tStartMJD = tStart.mjd[0]
            self.METStart = 239557417.0
            times = [str(date[6:10])+'-'+str(date[3:5])+'-'+str(date[:2])+'T'+str(date[11:19]),    str(date2[6:10])+'-'+str(date2[3:5])+'-'+str(date2[:2])+'T'+str(date2[11:19])              ]
            t = Time(times)
            t0 = t.mjd[0]
            t1 = t.mjd[1]
            self.tmin = (t0-self.tStartMJD)*86400 + self.METStart
            self.tmax = (t1-self.tStartMJD)*86400 + self.METStart
            self.spacecraft_file = self.white_box_spacecraft_file.text()


            self.Time_intervMJD = t1-t0 
            
            catalog = self.comboBox_Catalog.currentText()
            
            f = open(self.OutputDir+'config.yaml','w')
            f.write('data:\n')
            os.system('ls '+self.white_box_photon_dir.text()+'/*PH*.fits > '+self.OutputDir+'list.txt')
            f.write('  evfile : '+self.OutputDir+'list.txt\n')   
            f.write('  scfile : '+self.spacecraft_file+'\n')
            if self.checkBox_External_ltcube.isChecked():
                ltcube_list = np.loadtxt(self.white_box_External_ltcube.text(),ndmin=2,dtype=str)
                ltcube_list = ltcube_list[:,0]
                if len(ltcube_list) == 1:
                    f.write('  ltcube : '+ltcube_list[0]+'\n')
            
            
            Coords = self.white_box_RAandDec.text().split(",")

            if len(Coords) == 1:
                recovered_RA = self.recover_coords_from_name["RA"][0]
                recovered_Dec = self.recover_coords_from_name["DEC"][0]
                c = SkyCoord(recovered_RA+' '+recovered_Dec, unit=(u.hourangle, u.deg))
                self.RA = str(c.ra.value)
                self.Dec = str(c.dec.value)
            else:
                self.RA = Coords[0]
                self.Dec = Coords[1]
            
            Energies = self.white_box_energy.text().split(',')
            
            self.Emin = float(Energies[0])
            self.Emax = float(Energies[1])
            if self.Emin < 100:
                zmax = 80
                self.roiwidth = 17
            elif 100 <= self.Emin < 500:
                zmax = 90
                self.roiwidth = 15
            elif 500 <= self.Emin < 1000:
                zmax = 100
                self.roiwidth = 12
            else:
                zmax = 105
                self.roiwidth = 10
            
            f.write('\nbinning:\n')
            f.write('  roiwidth   : '+str(self.roiwidth)+'\n')
            f.write('  binsz      : 0.1\n')
            f.write('  binsperdec : 8\n\n')
            f.write('selection :\n')
            f.write('  emin : '+str(self.Emin)+'\n')
            f.write('  emax : '+str(self.Emax)+'\n')
            f.write('  zmax    : '+str(zmax)+'\n')
            f.write('  evclass : 128\n')
            if self.checkBox_highest_resolution.isChecked():
                f.write('  evtype  : 48\n')
            else:
                f.write('  evtype  : 3\n')
            f.write('  ra: '+self.RA+'\n')
            f.write('  dec: '+self.Dec+'\n')
            f.write('  tmin: '+str(int(self.tmin))+'\n')
            f.write('  tmax: '+str(int(self.tmax))+'\n\n')
            f.write('gtlike:\n')
            f.write('  edisp : True\n')
            f.write("  irfs : 'P8R3_SOURCE_V3'\n")
            f.write("  edisp_disable : ['isodiff']\n")
            f.write('  edisp_bins : -2\n\n')
            f.write('model:\n')
            f.write('  src_roiwidth : '+str(self.roiwidth+10)+'\n')
            f.write("  galdiff  : '"+self.white_box_Diffuse_dir.text()+"/gll_iem_v07.fits'\n")
            f.write("  isodiff  : '"+self.white_box_Diffuse_dir.text()+"/iso_P8R3_SOURCE_V3_v1.txt'\n")
            f.write("  catalogs : ['"+catalog+"']\n")            
            
            if self.comboBox_is_it_cataloged.currentText() == "No": 
                f.write("  sources  :\n")
                f.write("    - { name: '"+self.white_box_target_name.text()+"', ra : "+self.RA+", dec : "+self.Dec+", SpectrumType : 'PowerLaw', SpatialModel: 'PointSource' }")


            if self.checkBox_high_sensitivity.isChecked():
                if self.Emax < 500:
                    n_components = 1
                elif self.Emin < 500 and 500 <= self.Emax < 1000:
                    n_components = 2
                    emin = [self.Emin,500]
                    emax = [500,self.Emax]
                    if self.checkBox_highest_resolution.isChecked():
                        evtype = [48,56]                        
                    else:
                        evtype = [3,3]
                    
                    if self.Emin < 100:
                        zmax = [80,100]
                    else:
                        zmax = [90,100]

                elif 500 <= self.Emin < 1000 and self.Emax > 1000:
                    n_components = 2
                    emin = [self.Emin,1000]
                    emax = [1000,self.Emax]
                    zmax = [100,105]
                    if self.checkBox_highest_resolution.isChecked():
                        evtype = [56,3]
                    else:
                        evtype = [3,3]
                elif self.Emin < 500 and self.Emax > 1000:
                    n_components = 3
                    emin = [self.Emin,500,1000]
                    emax = [500,1000,self.Emax]
                    if self.checkBox_highest_resolution.isChecked():
                        evtype = [48,56,3]                        
                    else:
                        evtype = [3,3,3]

                    if self.Emin < 100:
                        zmax = [80,100,105]
                    else:
                        zmax = [90,100,105]
                else:
                    n_components = 1

                if n_components > 1:
                    f.write('\ncomponents:\n')
                    for i in range(n_components):
                        f.write("  - model:\n")
                        f.write("      galdiff  : '"+self.white_box_Diffuse_dir.text()+"/gll_iem_v07.fits'\n")
                        f.write("      isodiff  : '"+self.white_box_Diffuse_dir.text()+"/iso_P8R3_SOURCE_V3_v1.txt'\n")
                        f.write('    selection:\n')
                        f.write(f'      emin : {emin[i]}\n')
                        f.write(f'      emax : {emax[i]}\n')
                        f.write(f'      zmax : {zmax[i]}\n')
                        f.write(f'      evtype : {evtype[i]}\n')
                        if self.checkBox_External_ltcube.isChecked():
                            f.write('    data:\n')
                            f.write(f'      ltcube : {ltcube_list[i]}\n')

                            

            f.close()
        else:
            self.popup_go()
        
        return can_we_go
            
        
    def popup_tutorial(self):

        """
        Popup with a link to the video tutorials of easyfermi.
        """
        msg = QtWidgets.QMessageBox()
        msg.setWindowTitle("Quick tutorial")
        msg.setText("Please check the tutorial on YouTube:\n<a href='https://www.youtube.com/channel/UCeLCfEoWasUKky6CPNN_opQ'>Video tutorial</a>")
        msg.setTextFormat(QtCore.Qt.RichText)
               
        msg.setIcon(QtWidgets.QMessageBox.Information) #Information, Critical, Warning
        
        msg.exec_()
        
        
    def popup_credits(self):

        """
        Popup showing the credits.
        """

        msg = QtWidgets.QMessageBox()
        msg.setWindowTitle("Credits")
        msg.setText("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n- I would like to thank Clodomir Vianna, Fabio Cafardo, Lucas Costa Campos and Raí Menezes for their help and strong support in this project.\n- A big thanks to Alessandra Azzollini, Douglas Carlos, Kaori Nakashima, Lucas Siconato, Matt Pui, and Romana Grossova, the first users/testers of easyfermi. \n ")
        msg.setInformativeText("To acknowledge easyfermi, please cite <a href='https://ui.adsabs.harvard.edu/abs/2022arXiv220611272D/abstract'>de Menezes, R. (2022)</a>. Since easyfermi relies on fermipy, gammapy, and astropy, please also cite <a href='https://ui.adsabs.harvard.edu/abs/2017ICRC...35..824W/abstract'>Wood et al. (2017)</a>, <a href='https://ui.adsabs.harvard.edu/abs/2023A%26A...678A.157D/abstract'>Donath et al. (2023)</a>, and <a href='https://ui.adsabs.harvard.edu/abs/2018AJ....156..123A/abstract'>Astropy Collaboration (2018)</a>.\n\n")
        msg.setIcon(QtWidgets.QMessageBox.Information) #Information, Critical, Warning
        
        msg.exec_()
    
    def popup_go(self):

        """
        Warning popup called by check_for_errors().
        """

        msg = QtWidgets.QMessageBox()
        msg.setWindowTitle("Something is wrong")
        try:
            msg.setText(self.max_energy_wanning)
            del self.max_energy_wanning
        except:
            if self.radioButton_Standard.isChecked():
                msg.setText("One of the mandatory fields under 'Standard' is not properly filled. Check the Log for more details.")
            else:
                msg.setText("The mandatory field under 'Custom' is not properly filled.")
               
        msg.setIcon(QtWidgets.QMessageBox.Warning) #Information, Critical, Warning
        
        
        msg.exec_()
            
    def popup_block_download(self):

        """
        Warning popup if the download was already performed. 
        """

        msg = QtWidgets.QMessageBox()
        msg.setWindowTitle(f"{self.which_download} file(s) already found")
        msg.setText(f"We found one or more {self.which_download} files in the download directory. The download is canceled.\n")
        msg.setIcon(QtWidgets.QMessageBox.Warning) #Information, Critical, Warning
        msg.exec_()

    def read_yes_or_no_button(self, i):
        self.yes_or_no = i.text()


    def download_SC(self):

        """
        Function to donwload the spacecraft data.
        """

        self.download_spacecraft_is_over = False
        Coords = self.white_box_RAandDec.text()
        date = self.dateTimeEdit.text()
        date2 = self.dateTimeEdit_2.text()
        times = [str(date[6:10])+'-'+str(date[3:5])+'-'+str(date[:2])+' '+str(date[11:19]),    str(date2[6:10])+'-'+str(date2[3:5])+'-'+str(date2[:2])+' '+str(date2[11:19])   ]
            
        print("\nQuerying spacecraft data...")
        query = fermi.FermiLAT.query_object(Coords,searchradius="10",
                obsdates=f'{times[0]}, {times[1]}',LATdatatype="None",spacecraftdata=True)
        print("\nDownloading spacecraft data...")
        if not os.path.exists("./spacecraft"):
            os.mkdir("./spacecraft")
            
        dir_path = os.path.dirname(os.path.realpath("./spacecraft/*"))
        for i in range(len(query)):
            if query[i].split("/")[-1][-9:] == 'SC00.fits':
                query_spacecraft = query[i]

        self.white_box_spacecraft_file.setText(dir_path+'/'+query_spacecraft.split("/")[-1])
        if OS_name != "Darwin":
            os.system(f"wget -P ./spacecraft/ {query_spacecraft}")
        else:
            os.system(f"curl -o ./spacecraft/{query_spacecraft.split('/')[-1]} {query_spacecraft}")

        self.download_spacecraft_is_over = True

    def onAndOff_spacecraft_dowload_button(self):

        """
        Function that enables/disables the spacecraft download button.
        """

        if self.pushButton_Download_SC.isEnabled():
            self.pushButton_Download_SC.setEnabled(False)
        elif self.radioButton_Standard.isChecked():
            self.pushButton_Download_SC.setEnabled(True)
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+'- Spacecraft file was downloaded to ./spacecraft/\n')
        else:
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+'- Spacecraft file was downloaded to ./spacecraft/\n')

    def popup_download_SC(self):

        """
        Poput with details about the download of the spacecraft file.
        """
        can_we_go = self.check_for_erros(jump_paths=True)
        self.which_download = "Spacecraft"

        if can_we_go:

            try:
                data_path = glob.glob("./spacecraft/*.fits*")
                if len(data_path[0]) > 0:
                    self.popup_block_download()
            except:
                self.popup_SC = QtWidgets.QMessageBox()
                self.popup_SC.setWindowTitle("Preparing donwload")
                self.popup_SC.setText("Please make sure you have internet conection.\n\nThe following proceedure will be done:\n1) Query spacecraft data. This may take a few minutes.\n2) Download the spacecraft file to ./spacecraft. This may take a while depending on your internet conection.\n\nWould you like to proceed?")
                self.popup_SC.setIcon(QtWidgets.QMessageBox.Information) #Information, Critical, Warning
                self.popup_SC.setStandardButtons(QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes) # seperate buttons with "|"
                self.popup_SC.buttonClicked.connect(self.read_yes_or_no_button)
                self.popup_SC.exec_()
                
                if self.yes_or_no == "&Yes":
                    
                    exec(f"self.thread_Download{self.number_of_clicks_to_download} = QtCore.QThread()")
                    exec(f"self.Download{self.number_of_clicks_to_download} = Downloads()")
                    exec(f"self.Download{self.number_of_clicks_to_download}.moveToThread(self.thread_Download{self.number_of_clicks_to_download})")
                    exec(f"self.thread_Download{self.number_of_clicks_to_download}.started.connect( self.onAndOff_spacecraft_dowload_button )")
                    exec(f"self.thread_Download{self.number_of_clicks_to_download}.started.connect(self.Download{self.number_of_clicks_to_download}.run_download_SC)")
                    exec(f"self.Download{self.number_of_clicks_to_download}.finished_downloads.connect( self.onAndOff_spacecraft_dowload_button )")
                    exec(f"self.thread_Download{self.number_of_clicks_to_download}.finished.connect(self.thread_Download{self.number_of_clicks_to_download}.deleteLater)")
                    exec(f"self.Download{self.number_of_clicks_to_download}.progress.connect(self.reportProgress)")
                    exec(f"self.thread_Download{self.number_of_clicks_to_download}.start()")
                    
                    self.number_of_clicks_to_download = self.number_of_clicks_to_download + 1

        else:
            self.popup_go()

    def download_Photons(self):

        """
        Function to donwload the photon data.
        """

        self.download_photons_is_over = False
        Coords = self.white_box_RAandDec.text()
        Energies = self.white_box_energy.text().split(',')
        self.Emin = float(Energies[0])
        self.Emax = float(Energies[1])
        if self.Emin < 100:
            radius = 14
        elif 100 <= self.Emin < 500:
            radius = 12
        elif 500 <= self.Emin < 1000:
            radius = 11
        else:
            radius = 9
        date = self.dateTimeEdit.text()
        date2 = self.dateTimeEdit_2.text()
        times = [str(date[6:10])+'-'+str(date[3:5])+'-'+str(date[:2])+' '+str(date[11:19]),    str(date2[6:10])+'-'+str(date2[3:5])+'-'+str(date2[:2])+' '+str(date2[11:19])   ]
            
        print("\nQuerying Photon data...")
        query_photons = fermi.FermiLAT.query_object(Coords,searchradius=f"{radius}", energyrange_MeV=f'{self.Emin}, {self.Emax}',
                obsdates=f'{times[0]}, {times[1]}',LATdatatype="Photon",spacecraftdata=False)
        print("\nDownloading Photon data...")
        if not os.path.exists("./Photons"):
            os.mkdir("./Photons")

        dir_path = os.path.dirname(os.path.realpath("./Photons/*"))
        self.white_box_photon_dir.setText(dir_path)
        if OS_name != "Darwin":
            for i in range(len(query_photons)):
                if query_photons[i].split("/")[-1][-9:] != 'SC00.fits':
                    os.system(f"wget -P ./Photons/ {query_photons[i]}")
        else:
            for i in range(len(query_photons)):
                if query_photons[i].split("/")[-1][-9:] != 'SC00.fits':
                    os.system(f"curl -o ./Photons/{query_photons[i].split('/')[-1]} {query_photons[i]}")
        
        self.download_photons_is_over = True

    def onAndOff_photon_dowload_button(self):

        """
        Function that enables/disables the photon download button.
        """

        if self.pushButton_Download_Photons.isEnabled():
            self.pushButton_Download_Photons.setEnabled(False)
        elif self.radioButton_Standard.isChecked():
            self.pushButton_Download_Photons.setEnabled(True)
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+'- Photon files were downloaded to ./Photons/\n')
        else:
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+'- Photon files were downloaded to ./Photons/\n')

    def popup_download_Photons(self):

        """
        Poput with details about the download of the photon files.
        """

        can_we_go = self.check_for_erros(jump_paths=True)
        self.which_download = "Photon"

        if can_we_go:

            try:
                data_path = glob.glob("./Photons/*.fits*")
                if len(data_path[0]) > 0:
                    self.popup_block_download()
            except:
                self.popup_Photons = QtWidgets.QMessageBox()
                self.popup_Photons.setWindowTitle("Preparing donwload")
                self.popup_Photons.setText("Please make sure you have internet conection.\n\nThe following proceedure will be done:\n1) Query Photon data. This may take a few minutes.\n2) Download the Photon data files to ./Photons/. This may take a while depending on your internet conection.\n\nWould you like to proceed?")
                self.popup_Photons.setIcon(QtWidgets.QMessageBox.Information) #Information, Critical, Warning
                self.popup_Photons.setStandardButtons(QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes) # seperate buttons with "|"
                self.popup_Photons.buttonClicked.connect(self.read_yes_or_no_button)
                self.popup_Photons.exec_()
                
                if self.yes_or_no == "&Yes":
                    
                    exec(f"self.thread_Download{self.number_of_clicks_to_download} = QtCore.QThread()")
                    exec(f"self.Download{self.number_of_clicks_to_download} = Downloads()")
                    exec(f"self.Download{self.number_of_clicks_to_download}.moveToThread(self.thread_Download{self.number_of_clicks_to_download})")
                    exec(f"self.thread_Download{self.number_of_clicks_to_download}.started.connect( self.onAndOff_photon_dowload_button )")
                    exec(f"self.thread_Download{self.number_of_clicks_to_download}.started.connect(self.Download{self.number_of_clicks_to_download}.run_download_Photons)")
                    exec(f"self.Download{self.number_of_clicks_to_download}.finished_downloads.connect( self.onAndOff_photon_dowload_button )")
                    exec(f"self.thread_Download{self.number_of_clicks_to_download}.finished.connect(self.thread_Download{self.number_of_clicks_to_download}.deleteLater)")
                    exec(f"self.Download{self.number_of_clicks_to_download}.progress.connect(self.reportProgress)")
                    exec(f"self.thread_Download{self.number_of_clicks_to_download}.start()")
                    
                    self.number_of_clicks_to_download = self.number_of_clicks_to_download + 1
            
        else:
            self.popup_go()

    def download_Diffuse(self):

        """
        Function to donwload the diffuse data.
        """

        self.download_diffuse_is_over = False
        print("\nDownloading latest diffuse models...")
        if not os.path.exists("./Diffuse"):
            os.mkdir("./Diffuse")

        dir_path = os.path.dirname(os.path.realpath("./Diffuse/*"))
        self.white_box_Diffuse_dir.setText(dir_path)
        if OS_name != "Darwin":
            os.system("wget -P ./Diffuse/ https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/iso_P8R3_SOURCE_V3_v1.txt")
            os.system("wget -P ./Diffuse/ https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/4fgl/gll_iem_v07.fits")
        else:
            os.system("curl -o ./Diffuse/iso_P8R3_SOURCE_V3_v1.txt https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/iso_P8R3_SOURCE_V3_v1.txt")
            os.system("curl -o ./Diffuse/gll_iem_v07.fits https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/4fgl/gll_iem_v07.fits")

        self.download_diffuse_is_over = True

    def onAndOff_diffuse_dowload_button(self):

        """
        Function that enables/disables the diffuse download button.
        """

        if self.pushButton_Download_diffuse.isEnabled():
            self.pushButton_Download_diffuse.setEnabled(False)
        elif self.radioButton_Standard.isChecked():
            self.pushButton_Download_diffuse.setEnabled(True)
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+'- Diffise models were downloaded to ./Diffuse/\n')
        else:
            self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+'- Diffise models were downloaded to ./Diffuse/\n')

    def popup_download_Diffuse(self):

        """
        Poput with details about the download of the diffuse files.
        """

        self.which_download = "Diffuse"

        try:
            data_path = glob.glob("./Diffuse/*.fits*")
            if len(data_path[0]) > 0:
                self.popup_block_download()
        except:
            self.popup_Diffuse = QtWidgets.QMessageBox()
            self.popup_Diffuse.setWindowTitle("Preparing donwload")
            self.popup_Diffuse.setText("Please make sure you have internet conection.\n\nThe following files will be downloaded:\n1) The Galactic gll_iem_v07.fits model.\n2) The isotropic template iso_P8R3_SOURCE_V3_v1.txt.\n\nBoth files will be saved in ./Diffuse/\nThis may take a while depending on your internet conection.\n\nWould you like to proceed?")
            self.popup_Diffuse.setIcon(QtWidgets.QMessageBox.Information) #Information, Critical, Warning
            self.popup_Diffuse.setStandardButtons(QtWidgets.QMessageBox.No | QtWidgets.QMessageBox.Yes) # seperate buttons with "|"
            self.popup_Diffuse.buttonClicked.connect(self.read_yes_or_no_button)
            self.popup_Diffuse.exec_()
            
            if self.yes_or_no == "&Yes":

                exec(f"self.thread_Download{self.number_of_clicks_to_download} = QtCore.QThread()")
                exec(f"self.Download{self.number_of_clicks_to_download} = Downloads()")
                exec(f"self.Download{self.number_of_clicks_to_download}.moveToThread(self.thread_Download{self.number_of_clicks_to_download})")
                exec(f"self.thread_Download{self.number_of_clicks_to_download}.started.connect( self.onAndOff_diffuse_dowload_button )")
                exec(f"self.thread_Download{self.number_of_clicks_to_download}.started.connect(self.Download{self.number_of_clicks_to_download}.run_download_Diffuse)")
                exec(f"self.Download{self.number_of_clicks_to_download}.finished_downloads.connect( self.onAndOff_diffuse_dowload_button )")
                exec(f"self.thread_Download{self.number_of_clicks_to_download}.finished.connect(self.thread_Download{self.number_of_clicks_to_download}.deleteLater)")
                exec(f"self.Download{self.number_of_clicks_to_download}.progress.connect(self.reportProgress)")
                exec(f"self.thread_Download{self.number_of_clicks_to_download}.start()")
                
                self.number_of_clicks_to_download = self.number_of_clicks_to_download + 1




    def browsefiles(self):

        """
        This function is called by the toolbuttons in easyfermi and
        allows the browse of files such as the spacecraft file, the VHE table, the diffuse model files etc.
        """
        if self.toolButton_spacecraft.isChecked():
            fname = QFileDialog.getOpenFileName(self, 'Open file', '', '(*.fits)')
            self.white_box_spacecraft_file.setText(fname[0])
            self.toolButton_spacecraft.setChecked(False)
        if self.toolButton_photons.isChecked():
            #fname = QFileDialog.getOpenFileName(self, 'Open file', './', '(*.txt *.list *.dat)')
            dir_ = QFileDialog.getExistingDirectory(self, 'Open directory:', './', QtWidgets.QFileDialog.ShowDirsOnly)
            self.white_box_photon_dir.setText(dir_)
            self.toolButton_photons.setChecked(False)
        if self.toolButton_dir_diffuse.isChecked():
            dir_ = QFileDialog.getExistingDirectory(self, 'Open directory:', './', QtWidgets.QFileDialog.ShowDirsOnly)
            self.white_box_Diffuse_dir.setText(dir_)
            self.toolButton_dir_diffuse.setChecked(False)
        if self.toolButton_output_dir.isChecked():
            dir_ = QFileDialog.getExistingDirectory(self, 'Open directory:', './', QtWidgets.QFileDialog.ShowDirsOnly)
            self.white_box_output_dir.setText(dir_)
            self.toolButton_output_dir.setChecked(False)
        if self.toolButton_External_ltcube.isChecked():
            fname = QFileDialog.getOpenFileName(self, 'Open file', './', '(*.txt)')
            self.white_box_External_ltcube.setText(fname[0])
            self.toolButton_External_ltcube.setChecked(False)
        if self.toolButton_Custom.isChecked():
            fname = QFileDialog.getOpenFileName(self, 'Open file', './', '(*.yaml)')
            self.white_box_config_file.setText(fname[0])
            self.toolButton_Custom.setChecked(False)
        if self.toolButton_VHE.isChecked():
            fname = QFileDialog.getOpenFileName(self, 'Open file', './', '(*.fits)')
            self.white_box_VHE.setText(fname[0])
            self.toolButton_VHE.setChecked(False)
        
       
    def activate(self):

        """ This function activates/deactivates the entries in the main window"""

        if self.checkBox_LC.isChecked():
            self.spinBox_LC_N_time_bins.setEnabled(True)
            if OS_name != "Darwin":
                self.spinBox_N_cores_LC.setEnabled(True)
                self.label_N_cores_LC.setEnabled(True)
            self.label_N_time_bins_LC.setEnabled(True)
            self.checkBox_adaptive_binning.setEnabled(True)
            self.white_box_TS_threshold.setEnabled(True)
            if self.checkBox_adaptive_binning.isChecked():
                self.white_box_TS_threshold.setEnabled(True)
                self.label_TS_threshold.setEnabled(True)
                self.spinBox_N_iter.setEnabled(True)
                self.label_N_iter.setEnabled(True)
            else:
                self.white_box_TS_threshold.setEnabled(False)
                self.label_TS_threshold.setEnabled(False)
                self.spinBox_N_iter.setEnabled(False)
                self.label_N_iter.setEnabled(False)
            
        else:
            self.spinBox_LC_N_time_bins.setEnabled(False)
            self.spinBox_N_cores_LC.setEnabled(False)
            self.label_N_time_bins_LC.setEnabled(False)
            self.label_N_cores_LC.setEnabled(False)
            self.checkBox_adaptive_binning.setEnabled(False)
            self.white_box_TS_threshold.setEnabled(False)
            self.label_TS_threshold.setEnabled(False)
            self.spinBox_N_iter.setEnabled(False)
            self.label_N_iter.setEnabled(False)
            
        if self.checkBox_SED.isChecked():
            self.spinBox_SED_N_energy_bins.setEnabled(True)
            self.label_N_energy_bins.setEnabled(True)
            self.checkBox_use_local_index.setEnabled(True)
            self.white_box_redshift.setEnabled(True)
            self.label_redshift.setEnabled(True)
            self.comboBox_redshift.setEnabled(True)
            self.comboBox_MCMC.setEnabled(True)
            self.label_MCMC.setEnabled(True)
            self.white_box_VHE.setEnabled(True)
            self.toolButton_VHE.setEnabled(True)
        else:
            self.spinBox_SED_N_energy_bins.setEnabled(False)
            self.label_N_energy_bins.setEnabled(False)
            self.checkBox_use_local_index.setEnabled(False)
            self.label_redshift.setEnabled(False)
            self.white_box_redshift.setEnabled(False)
            self.comboBox_redshift.setEnabled(False)
            self.comboBox_MCMC.setEnabled(False)
            self.label_MCMC.setEnabled(False)
            self.white_box_VHE.setEnabled(False)
            self.toolButton_VHE.setEnabled(False)

            
        if self.checkBox_extension.isChecked():
            self.doubleSpinBox_extension_max_size.setEnabled(True)
            self.radioButton_disk.setEnabled(True)
            self.radioButton_2D_Gauss.setEnabled(True)
            self.label_Extension_max_size.setEnabled(True)
        else:
            self.doubleSpinBox_extension_max_size.setEnabled(False)
            self.radioButton_disk.setEnabled(False)
            self.radioButton_2D_Gauss.setEnabled(False)
            self.label_Extension_max_size.setEnabled(False)
            
        if self.checkBox_TSmap.isChecked():
            self.doubleSpinBox_Photon_index_TS.setEnabled(True)
            self.checkBox_residual_TSmap.setEnabled(True)
            self.label_photon_index_TS.setEnabled(True)
        else:
            self.doubleSpinBox_Photon_index_TS.setEnabled(False)
            self.checkBox_residual_TSmap.setEnabled(False)
            self.label_photon_index_TS.setEnabled(False)
        
                
        if self.checkBox_find_extra_sources.isChecked():
            self.doubleSpinBox_min_significance.setEnabled(True)
            self.doubleSpinBox_min_separation.setEnabled(True)
            self.label_min_significance.setEnabled(True)
            self.label_minimum_separation.setEnabled(True)
        else:
            self.doubleSpinBox_min_significance.setEnabled(False)
            self.doubleSpinBox_min_separation.setEnabled(False)
            self.label_min_significance.setEnabled(False)
            self.label_minimum_separation.setEnabled(False)
            
        if self.checkBox_External_ltcube.isChecked():
            self.white_box_External_ltcube.setEnabled(True)
            self.toolButton_External_ltcube.setEnabled(True)
        else:
            self.white_box_External_ltcube.setEnabled(False)
            self.toolButton_External_ltcube.setEnabled(False)
        
        if self.checkBox_delete_sources.isChecked():
            self.white_box_list_of_sources_to_delete.setEnabled(True)
        else:
            self.white_box_list_of_sources_to_delete.setEnabled(False)
            
        if self.checkBox_change_model.isChecked():
            self.comboBox_change_model.setEnabled(True)
        else:
            self.comboBox_change_model.setEnabled(False)
        
        if self.checkBox_minimizer.isChecked():
            self.comboBox_minimizer.setEnabled(True)
        else:
            self.comboBox_minimizer.setEnabled(False)
            
        if self.radioButton_Custom.isChecked():
            self.toolButton_Custom.setEnabled(True)
            self.white_box_config_file.setEnabled(True)
            
        else:
            self.toolButton_Custom.setEnabled(False)
            self.white_box_config_file.setEnabled(False)
        
        if self.radioButton_free_source_radius_customized.isChecked():
            self.white_box_radius.setEnabled(True)
            self.label_radius.setEnabled(True)
            self.checkBox_only_norm.setEnabled(True)
            self.checkBox_freeze_gal.setEnabled(True)
            self.checkBox_freeze_iso.setEnabled(True)
        else:
            self.white_box_radius.setEnabled(False)
            self.label_radius.setEnabled(False)
            self.checkBox_only_norm.setEnabled(False)
            self.checkBox_freeze_gal.setEnabled(False)
            self.checkBox_freeze_iso.setEnabled(False)
        
            
        if self.radioButton_Standard.isChecked():
            self.toolButton_spacecraft.setEnabled(True)
            self.toolButton_photons.setEnabled(True)
            self.toolButton_dir_diffuse.setEnabled(True)
            self.white_box_RAandDec.setEnabled(True)
            self.white_box_energy.setEnabled(True)
            self.white_box_spacecraft_file.setEnabled(True)
            self.white_box_Diffuse_dir.setEnabled(True)
            self.white_box_photon_dir.setEnabled(True)
            self.label_dir_spacecraft.setEnabled(True)
            self.label_dir_photons.setEnabled(True)
            self.label_start_time.setEnabled(True)
            self.label_end_time.setEnabled(True)
            self.label_limits.setEnabled(True)
            self.label_dir_diffuse.setEnabled(True)
            self.label_RAandDec.setEnabled(True)
            self.label_energy.setEnabled(True)
            self.label_Catalog.setEnabled(True)
            self.label_is_it_cataloged.setEnabled(True)
            self.comboBox_Catalog.setEnabled(True)
            self.comboBox_is_it_cataloged.setEnabled(True)
            self.dateTimeEdit.setEnabled(True)
            self.dateTimeEdit_2.setEnabled(True)
            self.checkBox_External_ltcube.setEnabled(True)
            self.checkBox_highest_resolution.setEnabled(True)
            self.checkBox_high_sensitivity.setEnabled(True)
            self.pushButton_config.setEnabled(True)
            try:
                if self.download_spacecraft_is_over:
                    self.pushButton_Download_SC.setEnabled(True)
            except:
                self.pushButton_Download_SC.setEnabled(True)
            try:
                if self.download_photons_is_over:
                    self.pushButton_Download_Photons.setEnabled(True)
            except:
                self.pushButton_Download_Photons.setEnabled(True)
            try:
                if self.download_diffuse_is_over:
                    self.pushButton_Download_diffuse.setEnabled(True)
            except:
                self.pushButton_Download_diffuse.setEnabled(True)
        else:
            self.toolButton_spacecraft.setEnabled(False)
            self.toolButton_photons.setEnabled(False)
            self.toolButton_dir_diffuse.setEnabled(False)
            self.white_box_RAandDec.setEnabled(False)
            self.white_box_energy.setEnabled(False)
            self.white_box_spacecraft_file.setEnabled(False)
            self.white_box_Diffuse_dir.setEnabled(False)
            self.white_box_photon_dir.setEnabled(False)
            self.label_dir_spacecraft.setEnabled(False)
            self.label_dir_photons.setEnabled(False)
            self.label_start_time.setEnabled(False)
            self.label_end_time.setEnabled(False)
            self.label_limits.setEnabled(False)
            self.label_dir_diffuse.setEnabled(False)
            self.label_RAandDec.setEnabled(False)
            self.label_energy.setEnabled(False)
            self.label_Catalog.setEnabled(False)
            self.label_is_it_cataloged.setEnabled(False)
            self.comboBox_Catalog.setEnabled(False)
            self.comboBox_is_it_cataloged.setEnabled(False)
            self.dateTimeEdit.setEnabled(False)
            self.dateTimeEdit_2.setEnabled(False)
            self.checkBox_External_ltcube.setEnabled(False)
            self.checkBox_highest_resolution.setEnabled(False)
            self.checkBox_high_sensitivity.setEnabled(False)
            self.pushButton_config.setEnabled(False)
            self.pushButton_Download_SC.setEnabled(False)
            self.pushButton_Download_Photons.setEnabled(False)
            self.pushButton_Download_diffuse.setEnabled(False)
            if self.checkBox_External_ltcube.isChecked():
                self.white_box_External_ltcube.setEnabled(False)
                self.toolButton_External_ltcube.setEnabled(False)
            

        if self.comboBox_is_it_cataloged.currentText() == 'Yes':
            self.white_box_target_name.setEnabled(False)
        elif self.radioButton_Standard.isChecked():
            self.white_box_target_name.setEnabled(True)
        else:
            self.white_box_target_name.setEnabled(False)

    def find_nearest(self,array, value):

        """
        This function finds the nearest value in a numpy array.
        
        Parameters
        ----------
        array: numpy array
            Integer or float array 
        value: float
            The number that you are looking for.

        Returns
        -------
        idx: int
            the index corresponding to the element in the 'array' which is closest to 'value'
        """

        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    def Sun_path(self):

        """
        Function that calculates the target-Sun distance over
        the time window set in the main window of easyfermi.

        Returns
        -------
        Quickplot_Sun: pdf
            pdf plot in the output directory.

        """

        timeStart = ['2008-08-04T15:43:36']
        tStart = Time(timeStart)
        self.tStartMJD = tStart.mjd[0]
        self.METStart = 239557417.0
        
        output_format = self.comboBox_output_format.currentText()

        spacecraft_data_file = pyfits.open(self.spacecraft_file)
        spacecraft_data_file = spacecraft_data_file[1].data

        START = spacecraft_data_file["START"]
        STOP = spacecraft_data_file["STOP"]
        time_window_min_index = self.find_nearest(START, self.tmin)
        time_window_max_index = self.find_nearest(STOP, self.tmax)


        # Selecting the Sun RA and Dec within the chosen time window:
        Sun_RA = spacecraft_data_file["RA_SUN"][time_window_min_index:time_window_max_index]
        Sun_DEC = spacecraft_data_file["DEC_SUN"][time_window_min_index:time_window_max_index]
        Sun_RA = Sun_RA[::1000]
        Sun_DEC = Sun_DEC[::1000]
        target_RA = float(self.RA)
        target_Dec = float(self.Dec)


        Coords_Sun = SkyCoord(Sun_RA, Sun_DEC, frame='icrs', unit='deg')
        Coords_target = SkyCoord(target_RA, target_Dec, frame='icrs', unit='deg')
        self.Solar_separation = Coords_Sun.separation(Coords_target).value

        t0_MJD = self.tStartMJD + (START[time_window_min_index] - self.METStart)/86400  # Computing the MJD in the beginning of observations
        t1_MJD = self.tStartMJD + (STOP[time_window_max_index] - self.METStart)/86400

        if len(self.Solar_separation) >= 10:
            time_range = np.linspace(t0_MJD,t1_MJD,len(self.Solar_separation))

            f = plt.figure(figsize=(9,4),dpi=250) 
            ax = f.add_subplot(1,1,1)
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
            ax.tick_params(which='major', length=5, direction='in')
            ax.tick_params(which='minor', length=2.5, direction='in',bottom=True, top=True, left=True, right=True)
            ax.tick_params(bottom=True, top=True, left=True, right=True)
            ax.grid(linestyle=':',which='both')
            
            plt.plot(time_range,15*np.ones(len(time_range)))
            plt.plot(time_range,self.Solar_separation)
            plt.fill_between(time_range,15,alpha=0.4,color="gray")
            shift = t1_MJD - t0_MJD
            shift = 0.05*shift
            plt.text(t0_MJD + shift,15,"Solar contamination is possible below this line")
            plt.xlabel('Time [MJD]')
            plt.ylabel('Angular separation [$^{\circ}$]')
            plt.title('Angular separation between '+self.sourcename+' and the Sun',fontsize=11)
            
            plt.xlim(time_range[0],time_range[-1])
            plt.ylim(0,None)
            plt.tight_layout()
            plt.savefig(self.OutputDir+'Quickplot_Sun.'+output_format,bbox_inches='tight')

        else:
            print("Time window is too short for computing the target-Sun separation plot. Minimum window is 4 days.")




    
    def analysisBasics(self):
    
        """
        This function calls fermipy to optimize the RoI model in the Fermi-LAT data after the setup is done (i.e. after the computation
        of ltcube, exposure map, srcmaps etc). This function is also in charge of generating a counts map, changing the target model
        (if requested by the user), and deleting sources from the model (if requested).
        
        Parameters
        ----------
        self: instance of the class Ui_mainWindow
            This parameter contains all the variables read from the Graphical interface.

        Returns
        -------
        cmap.fits: data file saved in the output directory.
            This is the counts map centered on the target.
        """

        #Cmap:
        h = pyfits.open(self.OutputDir+'ccube.fits')
        counts = h[0].data
        counts.shape
        h[0].data = np.sum(counts,axis=0)
        h.writeto(self.OutputDir+'cmap.fits',overwrite=True)

        #Change model if requested:
        if self.checkBox_change_model.isChecked():
            aux = self.comboBox_change_model.currentText()
            if aux == 'Select...':
                pass
            elif aux == 'Power-law':
                aux = 'PowerLaw'
            elif aux == 'Power-law2':
                aux = 'PowerLaw2'
            elif aux == 'LogPar':
                aux = 'LogParabola'
            elif aux == 'PLEC':
                aux = 'PLSuperExpCutoff'
            elif aux == 'PLEC2':
                aux = 'PLSuperExpCutoff2'
            elif aux == 'PLEC3':
                aux = 'PLSuperExpCutoff3'
            elif aux == 'PLEC4':
                aux = 'PLSuperExpCutoff4'
            elif aux == 'BPL':
                aux = 'BrokenPowerLaw'
            elif aux == 'ExpCutOff-EBL':
                aux = 'ExpCutoff'
            
            if aux != 'Select...':    
                self.gta.delete_source(self.sourcename)
                skip_list = list(self.gta.roi.create_source_table()["source_name"])
                self.gta.add_source(self.sourcename, src_dict={'ra' : float(self.RA), 'dec' : float(self.Dec) , 'SpatialModel' : 'PointSource', 'SpectrumType' : aux})
                self.gta.optimize(skip=skip_list, npred_frac=0)  # We optimize only the new source added above. All other sources are skipped.
        
        # Optimizing RoI:
        self.gta.optimize(npred_threshold=50,shape_ts_threshold=30)

        if self.checkBox_find_extra_sources.isChecked():
            self.gta.find_sources(sqrt_ts_threshold=self.doubleSpinBox_min_significance.value(), min_separation=self.doubleSpinBox_min_separation.value(), multithread=True)
        
        #Delete sources:
        if self.checkBox_delete_sources.isChecked():
            a = self.white_box_list_of_sources_to_delete.text().split(',')
            if a[0] == '':
                pass
            else:    
                for i in a:
                    self.gta.delete_source(i)
        
        
        self.freeradius = self.roiwidth/2.
        self.freeradiusalert = 'ok'
        if self.radioButton_free_source_radius_customized.isChecked():
            if self.white_box_radius.text() != '':
                self.freeradius = float(self.white_box_radius.text())
            else:
                self.freeradiusalert = '- No free source radius available. Using default free source radius: '+str(self.freeradius)+"°."
            
            if self.checkBox_only_norm.isChecked():
                """Free only the normalizations:"""
                self.gta.free_sources(distance=self.freeradius,pars='norm')
            else:
                self.gta.free_sources(distance=self.freeradius)
                
            if self.checkBox_freeze_gal.isChecked():
                """Freeze Galactic diffuse model:"""
                self.gta.free_source('galdiff',free=False)
            else:
                pass
                
            if self.checkBox_freeze_iso.isChecked():
                """Freeze Isotropic diffuse model:"""
                self.gta.free_source('isodiff',free=False)
            else:
                pass
        else:    
            self.gta.free_sources(distance=self.freeradius)
            #self.gta.free_source('galdiff')
            #self.gta.free_source('isodiff')
            #self.gta.free_source(self.sourcename)
        
        if self.checkBox_freeze_spec_shape.isChecked():
            self.gta.free_source(self.sourcename,free=False)
            self.gta.free_source(self.sourcename,pars='norm')

        N_iter_adaptive_LC = self.spinBox_LC_N_time_bins.value()
        
        return N_iter_adaptive_LC
        
    def fit_model(self):
        
        """
        This function calls fermipy to fit the RoI model in the Fermi-LAT data after the optimization is done.
        Furthermore, it also generates the diagnostic plots (if requested by the user).
        
        Parameters
        ----------
        self: instance of the class Ui_mainWindow
            This parameter contains all the variables read from the Graphical interface.

        Returns
        -------
        do_diagnostic_plots: boolean
            This is the status of the diagnostic plot check box. If True, the Solar separation will be computed in the Class Worker.
        
        Target_results.txt: data file saved in the output directory.
            This file contains the results of the RoI fit only for the target, such that the user can have a quick look at it.
        Results.fits and Results.npy: data files saved in the output directory.
            These files contain the same information. They give the user the full set of results regarding the fit of the RoI.
        """

        if self.checkBox_minimizer.isChecked():
            optimizer = self.comboBox_minimizer.currentText()
        else:
            optimizer = 'NEWMINUIT'

        fit_results = self.gta.fit(optimizer=optimizer)
        original_fit_quality = fit_results['fit_quality']

        if original_fit_quality == 3:
            self.fitquality = '- Fit quality: 3. Excellent fit! Full accurate covariance matrix.'
        else:
            if self.gta.get_sources()[0]["ts"] < 25:
                self.gta.delete_sources(minmax_ts=[None,0.99*self.gta.get_sources()[0]["ts"]])  # Here we delete the sources that have TS < TS_target and refit the model
                fit_results = self.gta.fit(optimizer='NEWMINUIT')
                if fit_results['fit_quality'] == 3:
                    self.fitquality = f'- Original fit quality: {original_fit_quality}. As the TS of the target was only {self.gta.get_sources()[0]["ts"]}, we deleted all sources with TS < {self.gta.get_sources()[0]["ts"]} and performed the fit again.\n- Fit quality of the second try: 3. Excellent fit! Full accurate covariance matrix.'
                elif fit_results['fit_quality'] == 2:
                    self.fitquality = f'- Original fit quality: {original_fit_quality}. As the TS of the target was only {self.gta.get_sources()[0]["ts"]}, we deleted all sources with TS < {self.gta.get_sources()[0]["ts"]} and performed the fit again.\n- Fit quality of the second try: 2. Reasonable fit. Full matrix, but forced positive-definite (i.e. not accurate).'
                elif fit_results['fit_quality'] == 1:
                    self.fitquality = f'- Original fit quality: {original_fit_quality}. As the TS of the target was only {self.gta.get_sources()[0]["ts"]}, we deleted all sources with TS < {self.gta.get_sources()[0]["ts"]} and performed the fit again.\n- Fit quality of the second try: 1. Poor fit. Diagonal approximation only, not accurate.'
                else:
                    self.fitquality = f'- Original fit quality: {original_fit_quality}. As the TS of the target was only {self.gta.get_sources()[0]["ts"]}, we deleted all sources with TS < {self.gta.get_sources()[0]["ts"]} and performed the fit again.\n- Fit quality of the second try: 0. Bad fit. Error matrix not calculated.'
                
            else:
                self.gta.delete_sources(minmax_ts=[None,25])  # Here we delete the sources that have TS < 25 and refit the model
                fit_results = self.gta.fit(optimizer='NEWMINUIT')
                if fit_results['fit_quality'] == 3:
                    self.fitquality = f'- Original fit quality: {original_fit_quality}. We deleted all sources with TS < 25 and performed the fit again.\n- Fit quality of the second try: 3. Excellent fit! Full accurate covariance matrix.'
                elif fit_results['fit_quality'] == 2:
                    self.fitquality = f'- Original fit quality: {original_fit_quality}. We deleted all sources with TS < 25 and performed the fit again.\n- Fit quality of the second try: 2. Reasonable fit. Full matrix, but forced positive-definite (i.e. not accurate).'
                elif fit_results['fit_quality'] == 1:
                    self.fitquality = f'- Original fit quality: {original_fit_quality}. We deleted all sources with TS < 25 and performed the fit again.\n- Fit quality of the second try: 1. Poor fit. Diagonal approximation only, not accurate.'
                else:
                    self.fitquality = f'- Original fit quality: {original_fit_quality}. We deleted all sources with TS < 25 and performed the fit again.\n- Fit quality of the second try: 0. Bad fit. Error matrix not calculated.'
                

        #Do plots:
        do_diagnostic_plots = self.checkBox_diagnostic_plots.isChecked()
        if do_diagnostic_plots:
            self.gta.write_roi('Results',make_plots=True)
        else:
            self.gta.write_roi('Results',make_plots=False)
        
        #Saving results
        f = open(self.OutputDir+'Target_results.txt','w')
        f.write(str(self.gta.roi[self.sourcename]))
        f.write('\nEnergy flux upper limit (MeV cm-2 s-1): '+str(self.gta.roi.sources[0]['eflux_ul95']))
        f.write('\nPhoton flux upper limit (cm-2 s-1): '+str(self.gta.roi.sources[0]['flux_ul95']))
        f.close()
    
        return do_diagnostic_plots
    

    def relocalize_the_target(self):

        """This function finds the best-fit position for the gamma-ray target.

        Parameters
        ----------
        self: instance of the class Ui_mainWindow
            This parameter contains all the variables read from the Graphical interface.

        Returns
        -------
        TARGET_NAME_loc.fits:
            Table containing the old and new coordinates of the gamma-ray target, as well as the r68, r95, and r99 uncertainty radii. File is saved in the output directory read from the graphical interface.
        TARGET_NAME_loc.npy:
            Same as above but in the npy format.
        TARGET_NAME_localize.png and TARGET_NAME_localize_peak.png:
            A couple of plots showing the shift in the target coordinates.

        """

        if self.checkBox_relocalize.isChecked():
            loc = self.gta.localize(self.sourcename, make_plots=True,update=True)
            self.locRA = loc['ra']
            self.locDec = loc['dec']
            self.locr68 = loc['pos_r68']
            self.locr95 = loc['pos_r95']
            #print("New position. RA = ",self.locRA,", Dec = ",self.locDec,", r_68 = ",self.locr68,", r_95 = ",self.locr95)
            #self.large_white_box_Log.setPlainText(self.large_white_box_Log.toPlainText()+"Localization results saved in the file "+self.sourcename+"_loc.fits\n")

            #Saving results
            f = open(self.OutputDir+'Target_results.txt','w')
            f.write(str(self.gta.roi[self.sourcename]))
            f.write('\nEnergy flux upper limit (MeV cm-2 s-1): '+str(self.gta.roi.sources[0]['eflux_ul95']))
            f.write('\nPhoton flux upper limit (cm-2 s-1): '+str(self.gta.roi.sources[0]['flux_ul95']))
            f.close()

    def compute_TSmap(self):
        
        """This function computes the TS maps.

        Parameters
        ----------
        self: instance of the class Ui_mainWindow
            This parameter contains all the variables read from the Graphical interface.

        Returns
        -------
        Target_TS_map_pointsource_powerlaw_2.00_tsmap.fits:
            TS map highlighting the target significance computed for a point-source with power-law index defined by the user (default = 2). File is saved in the output directory read from the graphical interface.
        Residuals_TS_map__pointsource_powerlaw_2.00_tsmap.fits:
            Residuals TS map assuming the power-law index defined by the user (default = 2). File is saved in the output directory read from the graphical interface.
        
        """

        output_format = self.comboBox_output_format.currentText()

        if self.checkBox_TSmap.isChecked():
            model = {'Index' : self.doubleSpinBox_Photon_index_TS.value(), 'SpatialModel' : 'PointSource'}
            if self.checkBox_residual_TSmap.isChecked():
                TSmap_res = self.gta.tsmap('Residuals_TS_map_', model=model)
                plt.figure(figsize=(8,8))
                ROIPlotter(TSmap_res['sqrt_ts'],roi=self.gta.roi).plot(vmin=0,vmax=5,levels=[3,5,7,9],subplot=111,cmap='magma')
                plt.gca().set_title('sqrt(TS)')
                plt.savefig(self.OutputDir+'TSmap_residuals.'+output_format,bbox_inches='tight')
                          

            TSmap = self.gta.tsmap('Target_TS_map_', model=model, exclude=self.sourcename)
            plt.figure(figsize=(8,8))
            ROIPlotter(TSmap['sqrt_ts'],roi=self.gta.roi).plot(vmin=0,vmax=5,levels=[3,5,7,9],subplot=111,cmap='magma')
            plt.gca().set_title('sqrt(TS)')
            plt.savefig(self.OutputDir+'TSmap_target_highlighted.'+output_format,bbox_inches='tight')
            
            #Below we compute the excess, significance, model, and data maps:
            self.gta.residmap('Excess_'+self.sourcename,model=model,make_plots=True, write_fits=True, write_npy=False)
        
        
    def compute_SED(self):

        """This function computes the SED for the target.

        Parameters
        ----------
        self: instance of the class Ui_mainWindow
            This parameter contains all the variables read from the Graphical interface.

        Returns
        -------
        TARGET_NAME_sed.fits
            Table containing the SED data. File is saved in the output directory read from the graphical interface.
        Quickplot_SED.png or Quickplot_SED.pdf:
            A plot of the SED for quick visualization. File is saved in the output directory read from the graphical interface.

        """

        output_format = self.comboBox_output_format.currentText()
        if self.checkBox_SED.isChecked():

            try:
                self.redshift = float(self.white_box_redshift.text())
            except:
                print("- WARNING: easyfermi could not read the redshift. Please insert a valid float number. Setting redshift to 0.0.\n")
                self.redshift = 0.0
                self.redshift_error = "- WARNING: easyfermi could not read the redshift. Please insert a valid float number. Setting redshift to 0.0.\n"

            c = np.load(self.OutputDir+'Results.npy',allow_pickle=True).flat[0]
            self.E = np.array(c['sources'][self.sourcename]['model_flux']['energies'])
            self.dnde = np.array(c['sources'][self.sourcename]['model_flux']['dnde'])
            dnde_hi = np.array(c['sources'][self.sourcename]['model_flux']['dnde_hi'])
            dnde_lo = np.array(c['sources'][self.sourcename]['model_flux']['dnde_lo'])
            TSmin = 9                

            Nbins = self.spinBox_SED_N_energy_bins.value() + 1
            ebins = np.linspace(np.log10(self.Emin),np.log10(self.Emax),num=Nbins)
            ebins = ebins.tolist()

            use_local_index = False
            if self.checkBox_use_local_index.isChecked():
                use_local_index = True

            self.sed = self.gta.sed(self.sourcename,loge_bins=ebins,make_plots=False,use_local_index=use_local_index,write_npy=False)
            
            self.Energy_data_points = self.sed['e_ctr'][self.sed['ts']>TSmin]
            self.e2dnde_data_points = self.sed['e2dnde'][self.sed['ts']>TSmin]
            self.xerr_data_points = [self.sed['e_ctr'][self.sed['ts']>TSmin] - self.sed['e_min'][self.sed['ts']>TSmin], self.sed['e_max'][self.sed['ts']>TSmin] - self.sed['e_ctr'][self.sed['ts']>TSmin]]
            self.yerr_data_points = self.sed['e2dnde_err'][self.sed['ts']>TSmin]
            self.Energy_uplims = self.sed['e_ctr'][self.sed['ts']<=TSmin]
            self.e2dnde_uplims = self.sed['e2dnde_ul95'][self.sed['ts']<=TSmin]
            self.xerr_uplims = [self.sed['e_ctr'][self.sed['ts']<=TSmin] - self.sed['e_min'][self.sed['ts']<=TSmin], self.sed['e_max'][self.sed['ts']<=TSmin] - self.sed['e_ctr'][self.sed['ts']<=TSmin]]
            self.yerr_uplims = 0.3*self.sed['e2dnde_ul95'][self.sed['ts']<=TSmin]


            f = plt.figure(figsize=(6,5),dpi=250)
            ax = f.add_subplot(1,1,1)
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
            ax.tick_params(which='major', length=5, direction='in')
            ax.tick_params(which='minor', length=2.5, direction='in',bottom=True, top=True, left=True, right=True)
            ax.tick_params(bottom=True, top=True, left=True, right=True)
            ax.grid(linestyle=':',which='both')

            plt.loglog(self.E, (self.E**2)*self.dnde, 'k--')
            plt.loglog(self.E, (self.E**2)*dnde_hi, 'k')
            plt.loglog(self.E, (self.E**2)*dnde_lo, 'k')
            plt.errorbar(self.Energy_data_points, self.e2dnde_data_points, xerr = self.xerr_data_points, yerr = self.yerr_data_points, color = "C0", markeredgecolor='black', fmt ='o', capsize=4)
            plt.errorbar(self.Energy_uplims, self.e2dnde_uplims, xerr = self.xerr_uplims, yerr = self.yerr_uplims, markeredgecolor='black', fmt='o', uplims=True, color='C0', capsize=4)
            plt.xlabel('Energy [MeV]')
            plt.ylabel(r'E$^{2}$ dN/dE [MeV cm$^{-2}$ s$^{-1}$]')
            plt.title(self.sourcename+' - SED')
            
            if len(self.sed['e2dnde'][self.sed['ts']>TSmin]) > 0:
                ymax = 2*(self.sed['e2dnde'][self.sed['ts']>TSmin] + self.sed['e2dnde_err'][self.sed['ts']>TSmin]).max()
                ymin = 0.5*(self.sed['e2dnde'][self.sed['ts']>TSmin] - self.sed['e2dnde_err'][self.sed['ts']>TSmin]).min()
                if ymax > 4*self.sed['e2dnde'][self.sed['ts']>TSmin].max():
                    ymax = 4*self.sed['e2dnde'][self.sed['ts']>TSmin].max()
                
                if ymin < self.sed['e2dnde'][self.sed['ts']>TSmin].min()/5.0:
                    ymin = self.sed['e2dnde'][self.sed['ts']>TSmin].min()/5.0
                    
                if len(self.sed['e2dnde_ul95'][self.sed['ts']<=TSmin]) > 0:
                    yaux = (self.sed['e2dnde_ul95'][self.sed['ts']<=TSmin]).max()
                    if yaux > ymax:
                        ymax = 2*yaux
                        
                    yaux = (self.sed['e2dnde_ul95'][self.sed['ts']<=TSmin]).min()
                    if yaux < ymin:
                        ymin = 0.5*yaux
            else:
                ymax = 2*self.sed['e2dnde_ul95'][self.sed['ts']<=TSmin].max()
                ymin = 0.5*self.sed['e2dnde_ul95'][self.sed['ts']<=TSmin].min()
                
            
            plt.ylim(ymin,ymax)
            plt.xlim(self.Emin*0.8,self.Emax*1.2)
            plt.tight_layout()
            plt.savefig(self.OutputDir+'Quickplot_SED_fermipy.'+output_format,bbox_inches='tight')


    def EBL_and_MCMC(self):

        """This function corrects the data for EBL absorption and applies a MCMC to the corrected data.
           If the redshift is 0.0, the MCMC is applied to the non-corrected LAT data.

        Parameters
        ----------
        self: instance of the class Ui_mainWindow
            This parameter contains all the variables read from the Graphical interface.

        Returns
        -------
        TARGET_NAME_sed.fits
            Table containing the SED data + MCMC results. File is saved in the output directory read from the graphical interface.
        Quickplot_SED_MCMC.png or Quickplot_SED_MCMC.pdf:
            A plot of the SED for quick visualization. File is saved in the output directory read from the graphical interface.

        """

        TSmin = 9 
        self.include_VHE = False
        if len(self.Energy_data_points) > 2:
            self.allow_MCMC = True
        else:
            self.allow_MCMC = False

        if self.checkBox_SED.isChecked() and self.allow_MCMC:
            if self.white_box_VHE.text().split('.')[-1] == 'fits':
                try:
                    self.include_VHE = True
                    VHE = pyfits.open(self.white_box_VHE.text())[1].data

                    energy_VHE = VHE["e_ref"] * 1000000  # TeV to MeV
                    energy_VHE_min = VHE["e_min"] * 1000000
                    energy_VHE_max = VHE["e_max"] * 1000000
                    e2dnde_VHE = VHE["e2dnde"] * 1e6  # TeV cm-2 s-1 to MeV cm-2 s-1
                    e2dnde_err_VHE = VHE["e2dnde_err"] * 1e6
                    e2dnde_UL_VHE = VHE["e2dnde_ul"] * 1e6
                    TS_VHE = VHE["ts"]

                    delete_nan_elements = np.isnan(e2dnde_UL_VHE)

                    TS_VHE = np.delete(TS_VHE, delete_nan_elements)
                    e2dnde_UL_VHE = np.delete(e2dnde_UL_VHE, delete_nan_elements)[TS_VHE <= TSmin]
                    energy_UL_VHE = np.delete(energy_VHE, delete_nan_elements)[TS_VHE <= TSmin]
                    energy_UL_VHE_min = np.delete(energy_VHE_min, delete_nan_elements)[TS_VHE <= TSmin]
                    energy_UL_VHE_max = np.delete(energy_VHE_max, delete_nan_elements)[TS_VHE <= TSmin]
                    energy_VHE = np.delete(energy_VHE, delete_nan_elements)[TS_VHE > TSmin]
                    energy_VHE_min = np.delete(energy_VHE_min, delete_nan_elements)[TS_VHE > TSmin]
                    energy_VHE_max = np.delete(energy_VHE_max, delete_nan_elements)[TS_VHE > TSmin]
                    e2dnde_VHE = np.delete(e2dnde_VHE, delete_nan_elements)[TS_VHE > TSmin]
                    e2dnde_err_VHE = np.delete(e2dnde_err_VHE, delete_nan_elements)[TS_VHE > TSmin]

                    x_err_VHE = [energy_VHE - energy_VHE_min, energy_VHE_max - energy_VHE]
                    x_err_UL_VHE =  [energy_UL_VHE - energy_UL_VHE_min, energy_UL_VHE_max - energy_UL_VHE]

                except:
                    self.include_VHE = False
            else:
                self.include_VHE = False

            if self.include_VHE:
                Energy_SED = np.concatenate([self.Energy_data_points, energy_VHE])
                Energy_err_SED = [np.concatenate([self.xerr_data_points[0],x_err_VHE[0]]), np.concatenate([self.xerr_data_points[1],x_err_VHE[1]])]
                dnde_SED = np.concatenate([self.e2dnde_data_points, e2dnde_VHE]) / (Energy_SED**2)  # y in dnde
                dnde_err_SED = np.concatenate([self.yerr_data_points, e2dnde_err_VHE]) / (Energy_SED**2)  # yerr in dnde
                self.Energy_uplims = np.concatenate([self.Energy_uplims, energy_UL_VHE])
                self.xerr_uplims = [np.concatenate([self.xerr_uplims[0],x_err_UL_VHE[0]]), np.concatenate([self.xerr_uplims[1],x_err_UL_VHE[1]])]
                self.e2dnde_uplims = np.concatenate([self.e2dnde_uplims, e2dnde_UL_VHE])
                self.yerr_uplims = np.concatenate([self.yerr_uplims, 0.3*e2dnde_UL_VHE])
                N_bins_VHE = len(energy_VHE)
                N_bins_VHE_UL = len(energy_UL_VHE)
            else:
                Energy_SED = self.Energy_data_points
                Energy_err_SED = self.xerr_data_points
                dnde_SED = self.e2dnde_data_points / (Energy_SED**2)
                dnde_err_SED = self.yerr_data_points / (Energy_SED**2)

            if self.redshift > 0.0:
                # Here we compute the EBL absorption model for a given redshift:
                EBL_model = self.comboBox_redshift.currentText().split(" ")[0]
                if EBL_model == "Dominguez":
                    EBL_model = "dominguez"
                elif EBL_model == "Franceschini":
                    if self.comboBox_redshift.currentText().split(" ")[1] == "et":
                        EBL_model = "franceschini"
                    else:
                        EBL_model = "franceschini17"
                elif EBL_model == "Saldana-Lopez":
                    EBL_model = "saldana-lopez21"
                else:
                    EBL_model = "finke"

                os.environ["GAMMAPY_DATA"] = EBLpath.name  # gammapy will look for EBL models in this directory
                absorption = EBLAbsorptionNormSpectralModel.read_builtin(EBL_model, redshift=self.redshift)
                abs_data = absorption.evaluate(Energy_SED*u.MeV,self.redshift, alpha_norm=1)
                dnde_data_points_deabsorbed = dnde_SED/abs_data
                dnde_error_deabsorbed = dnde_err_SED/abs_data


                abs_uplims = absorption.evaluate(self.Energy_uplims*u.MeV,self.redshift, alpha_norm=1)
                e2dnde_uplims_deabsorbed = self.e2dnde_uplims/abs_uplims
                yerr_uplims_deabsorbed = self.yerr_uplims/abs_uplims

                Energy_SED_deabsorbed_log = np.log10(Energy_SED)
                dnde_data_points_deabsorbed_log = np.log10(dnde_data_points_deabsorbed)
                dnde_error_deabsorbed_log = (1 / np.log(10)) * (1 / dnde_data_points_deabsorbed) * dnde_error_deabsorbed  # converting the errors to log scale

                # data to be saved later:
                absorption_all_data = absorption.evaluate(self.sed['e_ctr']*u.MeV,self.redshift, alpha_norm=1)
                e2dnde_all_data_unabsorbed = self.sed['e2dnde']/absorption_all_data
                e2dnde_all_data_errors_unabsorbed = self.sed['e2dnde_err']/absorption_all_data
                e2dnde_all_ul95_unabsorbed = self.sed['e2dnde_ul95']/absorption_all_data
            else:
                Energy_SED_log = np.log10(Energy_SED)
                dnde_data_points_log = np.log10(dnde_SED)
                dnde_error_log = (1 / np.log(10)) * (1 / dnde_SED) * dnde_err_SED  # converting the errors to log scale


            def plotter(sampler, x):

                """
                Function to plot the distribution of possible models around the maximum likelihood model.
                
                Parameters
                ----------
                sampler: array
                    table where each column contains the distribution of possible values for a given parameter
                x: array
                    the x-axis of the SED to be plotted.

                Returns
                -------
                """

                samples = sampler.flatchain
                for theta in samples[np.random.randint(len(samples), size=100)]:
                    if self.comboBox_MCMC.currentText() == "LogPar":
                        model = 10 ** (2 * x + LogPar(theta, x))
                    elif self.comboBox_MCMC.currentText() == "LogPar2":
                        model = 10 ** LogPar2(theta, x)
                    elif self.comboBox_MCMC.currentText() == "PLEC":
                        model = 10 ** (2 * x + PLEC(theta, x))
                    elif self.comboBox_MCMC.currentText() == "PLEC_bfix":
                        model = 10 ** (2 * x + PLEC_bfix(theta, x))
                    elif self.comboBox_MCMC.currentText() == "PowerLaw":
                        model = 10 ** (2 * x + PowerLaw(theta, x))
                    plt.plot(10**x, model, color="r", zorder=0, alpha=0.1)  # plotting with parameters in the posterior distribution
                plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
                plt.xlabel("Energy [MeV]")
                plt.ylabel("E$^2dN/dE$ [MeV cm$^{-2}$ s$^{-1}$]")
                plt.grid(linestyle=":")
                plt.legend()

            def PowerLaw(theta, x):
                N0, alpha = theta
                Ep = np.log10(2*self.Emin)
                return N0 - alpha*(x - Ep)
            
            def LogPar(theta, x):
                N0, alpha, beta = theta
                Ep = np.log10(2*self.Emin)
                return N0 + (-alpha - beta * np.log((10**x) / (10**Ep))) * (x - Ep)
            
            def LogPar2(theta, x):
                Splog, alpha, Ep = theta
                return Splog + (-alpha*(np.log10((10**x)/(10**Ep))**2))

            def PLEC(theta, x):
                Ep = np.log10(2*self.Emin)
                N0, alpha, Ec, b = theta
                # N0, alpha, Ec = theta
                return N0 - alpha * (x - Ep) + np.log10(np.exp(-((10**x / (10**Ec)) ** b)))

            def PLEC_bfix(theta, x):
                Ep = np.log10(2*self.Emin)
                b = 1
                N0, alpha, Ec = theta
                # N0, alpha, Ec = theta
                return N0 - alpha * (x - Ep) + np.log10(np.exp(-((10**x / (10**Ec)) ** b)))

            def lnlike(theta, x, y, yerr):
                if self.comboBox_MCMC.currentText() == "LogPar":
                    likelihood = -0.5 * np.sum(((y - LogPar(theta, x)) / yerr) ** 2)
                elif self.comboBox_MCMC.currentText() == "LogPar2":
                    likelihood = -0.5 * np.sum(((y - LogPar2(theta, x)) / yerr) ** 2)
                elif self.comboBox_MCMC.currentText() == "PLEC":
                    likelihood = -0.5 * np.sum(((y - PLEC(theta, x)) / yerr) ** 2)
                elif self.comboBox_MCMC.currentText() == "PLEC_bfix":
                    likelihood = -0.5 * np.sum(((y - PLEC_bfix(theta, x)) / yerr) ** 2)
                elif self.comboBox_MCMC.currentText() == "PowerLaw":
                    likelihood = -0.5 * np.sum(((y - PowerLaw(theta, x)) / yerr) ** 2)

                return likelihood

            def lnprior_PowerLaw(theta):
                N0, alpha = theta
                if -15 < N0 < -8 and 0.5 < alpha < 5.0:
                    return 0.0
                return -np.inf
            
            def lnprior_LogPar(theta):
                N0, alpha, beta = theta
                if -15 < N0 < -8 and 1.0 < alpha < 4.0 and -1 < beta < 1.0:
                    return 0.0
                return -np.inf
            
            def lnprior_LogPar2(theta):
                Splog, alpha, Ep = theta
                if -7 < Splog < -2 and -1.0 < alpha < 1.0 and 2 < Ep < 7:
                    return 0.0
                return -np.inf
            
            def lnprior_PLEC(theta):
                N0, alpha, Ec, b = theta
                # N0, alpha, Ec = theta
                if -15 < N0 < -8 and 1.0 < alpha < 4.0 and 3.0 < Ec < 7.0 and 0.2 < b < 3.0:
                    return 0.0
                return -np.inf
            
            def lnprior_PLEC_bfix(theta):
                N0, alpha, Ec = theta
                # N0, alpha, Ec = theta
                if -15 < N0 < -8 and 1.0 < alpha < 4.0 and 3.0 < Ec < 7.0:
                    return 0.0
                return -np.inf

            def lnprob(theta, x, y, yerr):
                if self.comboBox_MCMC.currentText() == "LogPar":
                    lp = lnprior_LogPar(theta)
                elif self.comboBox_MCMC.currentText() == "LogPar2":
                    lp = lnprior_LogPar2(theta)
                elif self.comboBox_MCMC.currentText() == "PLEC":
                    lp = lnprior_PLEC(theta)
                elif self.comboBox_MCMC.currentText() == "PLEC_bfix":
                    lp = lnprior_PLEC_bfix(theta)
                elif self.comboBox_MCMC.currentText() == "PowerLaw":
                    lp = lnprior_PowerLaw(theta)

                if not np.isfinite(lp):
                    return -np.inf
                return lp + lnlike(theta, x, y, yerr)


            if self.redshift > 0.0:
                if self.comboBox_MCMC.currentText() == "LogPar2":
                    e2dnde_err = (1 / np.log(10)) * (1 / (dnde_data_points_deabsorbed*(Energy_SED**2))) * (dnde_error_deabsorbed*(Energy_SED**2))  # converting the errors to log scale
                    data = (Energy_SED_deabsorbed_log, dnde_data_points_deabsorbed_log + 2*Energy_SED_deabsorbed_log, e2dnde_err)
                else:
                    data = (Energy_SED_deabsorbed_log, dnde_data_points_deabsorbed_log, dnde_error_deabsorbed_log)
            else:
                if self.comboBox_MCMC.currentText() == "LogPar2":
                    e2dnde_err = (1 / np.log(10)) * (1 / (dnde_SED*(Energy_SED**2))) * (dnde_err_SED*(Energy_SED**2))  # converting the errors to log scale
                    data = (Energy_SED_log, dnde_data_points_log + 2*Energy_SED_log, e2dnde_err)
                else:
                    data = (Energy_SED_log, dnde_data_points_log, dnde_error_log)

            nwalkers = 300
            niter = 500
            if self.comboBox_MCMC.currentText() == "LogPar":
                initial = np.array([-13, 1.7, 0.2])
            elif self.comboBox_MCMC.currentText() == "LogPar2":
                initial = np.array([-4.5, 0.2, 3.5])
            elif self.comboBox_MCMC.currentText() == "PLEC":
                initial = np.array([-13, 1.7, 5, 1])
            elif self.comboBox_MCMC.currentText() == "PLEC_bfix":
                initial = np.array([-13, 1.7, 5])
            elif self.comboBox_MCMC.currentText() == "PowerLaw":
                initial = np.array([-13, 2.0])

            ndim = len(initial)
            p0 = [
                np.array(initial) + 0.3 * np.random.randn(ndim) for i in range(nwalkers)
            ]  # p0 is the methodology of stepping from one place on a grid to the next.


            def main(p0, nwalkers, niter, ndim, lnprob, data):
                sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=data)

                print("Running burn-in...")
                p0, _, _ = sampler.run_mcmc(p0, 100)
                sampler.reset()

                print("Running production...")
                pos, prob, state = sampler.run_mcmc(p0, niter)

                return sampler, pos, prob, state

            sampler, pos, prob, state = main(p0, nwalkers, niter, ndim, lnprob, data)

            # Setting the x limits:
            xmin = np.log10(Energy_SED[0]-Energy_err_SED[0][0])
            xmax = np.log10(Energy_SED[-1]+Energy_err_SED[1][-1])
            
            try:
                xmin2 = np.log10(self.Energy_uplims[0]-self.xerr_uplims[0][0])
                if xmin2 < xmin:
                    xmin = xmin2
            except:
                pass
            try:
                xmax2 = np.log10(self.Energy_uplims[-1]+self.xerr_uplims[1][-1])
                if xmax2 > xmax:
                    xmax = xmax2
            except:
                pass

            x = np.linspace(xmin, xmax, 1000)  # defining a high resolution x for the plots

            # Saving MCMC parameters and data points corrected for EBL:
            SED_file = glob.glob(self.OutputDir+"*_sed.fits")[0]
            hdul = pyfits.open(SED_file)
            if self.redshift > 0.0:
                original_data_cols = hdul[1].data.columns
                new_col_e2dnde_unabsorbed_data = pyfits.Column(name="e2dnde_EBL_corrected",array=e2dnde_all_data_unabsorbed,format="D",unit="cm-2 MeV s-1")
                new_col_e2dnde_unabsorbed_errors = pyfits.Column(name="e2dnde_err_EBL_corrected",array=e2dnde_all_data_errors_unabsorbed,format="D",unit="cm-2 MeV s-1")
                new_col_e2dnde_ul95_unabsorbed = pyfits.Column(name="e2dnde_ul95_EBL_corrected",array=e2dnde_all_ul95_unabsorbed,format="D",unit="cm-2 MeV s-1")                
                all_cols = pyfits.BinTableHDU.from_columns(original_data_cols + new_col_e2dnde_unabsorbed_data + new_col_e2dnde_unabsorbed_errors + new_col_e2dnde_ul95_unabsorbed)
                hdul[1].data = all_cols.data
                hdul[1].name = "SED"

                if self.white_box_VHE.text().split('.')[-1] == 'fits':
                    
                    VHE_col_energy = VHE["e_ref"] * 1000000  # in MeV
                    VHE_col_energy = np.delete(VHE_col_energy, delete_nan_elements)
                    VHE_absorption = absorption.evaluate(VHE_col_energy*u.MeV,self.redshift, alpha_norm=1)
                    VHE_col_energy = pyfits.Column(name="energy",array= VHE_col_energy,format="D",unit="MeV")

                    VHE_col_energy_min = VHE["e_min"] * 1000000
                    VHE_col_energy_min = np.delete(VHE_col_energy_min, delete_nan_elements)
                    VHE_col_energy_min = pyfits.Column(name="energy_min",array= VHE_col_energy_min,format="D",unit="MeV")

                    VHE_col_energy_max = VHE["e_max"] * 1000000
                    VHE_col_energy_max = np.delete(VHE_col_energy_max, delete_nan_elements)
                    VHE_col_energy_max = pyfits.Column(name="energy_max",array= VHE_col_energy_max,format="D",unit="MeV")

                    VHE_col_e2dnde = VHE["e2dnde"] * 1e6  # TeV cm-2 s-1 to MeV cm-2 s-1
                    VHE_col_e2dnde = np.delete(VHE_col_e2dnde, delete_nan_elements)/VHE_absorption  # Correcting for EBL absorption
                    VHE_col_e2dnde = pyfits.Column(name="e2dnde_VHE",array= VHE_col_e2dnde,format="D",unit="MeV cm-2 s-1")

                    VHE_col_e2dnde_err = VHE["e2dnde_err"] * 1e6
                    VHE_col_e2dnde_err = np.delete(VHE_col_e2dnde_err, delete_nan_elements)/VHE_absorption  # Correcting for EBL absorption
                    VHE_col_e2dnde_err = pyfits.Column(name="e2dnde_VHE_err",array= VHE_col_e2dnde_err,format="D",unit="MeV cm-2 s-1")

                    VHE_col_e2dnde_UL = VHE["e2dnde_ul"] * 1e6
                    VHE_col_e2dnde_UL = np.delete(VHE_col_e2dnde_UL, delete_nan_elements)/VHE_absorption  # Correcting for EBL absorption
                    VHE_col_e2dnde_UL = pyfits.Column(name="e2dnde_VHE_UL95",array= VHE_col_e2dnde_UL,format="D",unit="MeV cm-2 s-1")

                    VHE_col_TS = VHE["ts"]
                    VHE_col_TS = np.delete(VHE_col_TS, delete_nan_elements)
                    VHE_col_TS = pyfits.Column(name="TS",array= VHE_col_TS,format="D")

                    VHE_table_EBL = pyfits.BinTableHDU.from_columns([VHE_col_energy, VHE_col_energy_min, VHE_col_energy_max, VHE_col_e2dnde, VHE_col_e2dnde_err, VHE_col_e2dnde_UL, VHE_col_TS])

            add_results = open(self.OutputDir+"Target_results.txt","a")
            add_results.write("\n\nMCMC results:\n")
            add_results.write("Model: "+self.comboBox_MCMC.currentText()+"\n")
            samples = sampler.flatchain
            theta_max = samples[np.argmax(sampler.flatlnprobability)]
            if self.comboBox_MCMC.currentText() == "LogPar":
                best_fit_model = LogPar(theta_max, x)
                N0 = np.quantile(samples[:,0],q=[0.16,0.5,0.84])
                Alpha = np.quantile(samples[:,1],q=[0.16,0.5,0.84])
                Beta = np.quantile(samples[:,2],q=[0.16,0.5,0.84])
                add_results.write(f"N0 (log scale): {N0[1]} - {N0[1]-N0[0]} + {N0[2]-N0[1]}\n")
                add_results.write(f"Alpha: {Alpha[1]} - {Alpha[1]-Alpha[0]} + {Alpha[2]-Alpha[1]}\n")
                add_results.write(f"Beta: {Beta[1]} - {Beta[1]-Beta[0]} + {Beta[2]-Beta[1]}\n")
                c1 = np.array(["N0 (log scale)","Alpha","Beta","Ep=2*Emin (log scale)"])
                c2 = np.array([N0[1],Alpha[1],Beta[1],np.log10(2*self.Emin)])
                c3 = np.array([N0[1]-N0[0],Alpha[1]-Alpha[0],Beta[1]-Beta[0],0])
                c4 = np.array([N0[2]-N0[1],Alpha[2]-Alpha[1],Beta[2]-Beta[1],0])
                c1 = pyfits.Column(name='Parameter', array=c1, format='22A')
                c2 = pyfits.Column(name='Value', array=c2, format='D')
                c3 = pyfits.Column(name='error_minus', array=c3, format='D')
                c4 = pyfits.Column(name='error_plus', array=c4, format='D')
                table_hdu = pyfits.BinTableHDU.from_columns([c1, c2, c3, c4])
                c1_posterior = pyfits.Column(name='N0 distribution (log scale)', array=samples[:,0], format='D')
                c2_posterior = pyfits.Column(name='Alpha distribution', array=samples[:,1], format='D')
                c3_posterior = pyfits.Column(name='Beta distribution', array=samples[:,2], format='D')
                table_hdu_posterior = pyfits.BinTableHDU.from_columns([c1_posterior, c2_posterior, c3_posterior])

            elif self.comboBox_MCMC.currentText() == "LogPar2":
                best_fit_model = LogPar2(theta_max, x)
                N0 = np.quantile(samples[:,0],q=[0.16,0.5,0.84])
                Alpha = np.quantile(samples[:,1],q=[0.16,0.5,0.84])
                Ep = np.quantile(samples[:,2],q=[0.16,0.5,0.84])
                add_results.write(f"N0 (log scale): {N0[1]} - {N0[1]-N0[0]} + {N0[2]-N0[1]}\n")
                add_results.write(f"Alpha: {Alpha[1]} - {Alpha[1]-Alpha[0]} + {Alpha[2]-Alpha[1]}\n")
                add_results.write(f"Ep (log scale): {Ep[1]} - {Ep[1]-Ep[0]} + {Ep[2]-Ep[1]}\n")
                c1 = np.array(["N0 (log scale)","Alpha","Ep (log scale)"])
                c2 = np.array([N0[1],Alpha[1],Ep[1]])
                c3 = np.array([N0[1]-N0[0],Alpha[1]-Alpha[0],Ep[1]-Ep[0]])
                c4 = np.array([N0[2]-N0[1],Alpha[2]-Alpha[1],Ep[2]-Ep[1]])
                c1 = pyfits.Column(name='Parameter', array=c1, format='22A')
                c2 = pyfits.Column(name='Value', array=c2, format='D')
                c3 = pyfits.Column(name='error_minus', array=c3, format='D')
                c4 = pyfits.Column(name='error_plus', array=c4, format='D')
                table_hdu = pyfits.BinTableHDU.from_columns([c1, c2, c3, c4])
                c1_posterior = pyfits.Column(name='N0 distribution (log scale)', array=samples[:,0], format='D')
                c2_posterior = pyfits.Column(name='Alpha distribution', array=samples[:,1], format='D')
                c3_posterior = pyfits.Column(name='Ep distribution', array=samples[:,2], format='D')
                table_hdu_posterior = pyfits.BinTableHDU.from_columns([c1_posterior, c2_posterior, c3_posterior])

            elif self.comboBox_MCMC.currentText() == "PLEC":
                best_fit_model = PLEC(theta_max, x)
                N0 = np.quantile(samples[:,0],q=[0.16,0.5,0.84])
                Alpha = np.quantile(samples[:,1],q=[0.16,0.5,0.84])
                Ec = np.quantile(samples[:,2],q=[0.16,0.5,0.84])
                b = np.quantile(samples[:,3],q=[0.16,0.5,0.84])
                add_results.write(f"N0 (log scale): {N0[1]} - {N0[1]-N0[0]} + {N0[2]-N0[1]}\n")
                add_results.write(f"Alpha: {Alpha[1]} - {Alpha[1]-Alpha[0]} + {Alpha[2]-Alpha[1]}\n")
                add_results.write(f"Ec: {Ec[1]} - {Ec[1]-Ec[0]} + {Ec[2]-Ec[1]}\n")
                add_results.write(f"b: {b[1]} - {b[1]-b[0]} + {b[2]-b[1]}\n")
                c1 = np.array(["N0 (log scale)","Alpha","Ec","b","Ep=2*Emin (log scale)"])
                c2 = np.array([N0[1],Alpha[1],Ec[1],b[1],np.log10(2*self.Emin)])
                c3 = np.array([N0[1]-N0[0],Alpha[1]-Alpha[0],Ec[1]-Ec[0],b[1]-b[0],0])
                c4 = np.array([N0[2]-N0[1],Alpha[2]-Alpha[1],Ec[2]-Ec[1],b[2]-b[1],0])
                c1 = pyfits.Column(name='Parameter', array=c1, format='22A')
                c2 = pyfits.Column(name='Value', array=c2, format='D')
                c3 = pyfits.Column(name='error_minus', array=c3, format='D')
                c4 = pyfits.Column(name='error_plus', array=c4, format='D')
                table_hdu = pyfits.BinTableHDU.from_columns([c1, c2, c3, c4])
                c1_posterior = pyfits.Column(name='N0 distribution (log scale)', array=samples[:,0], format='D')
                c2_posterior = pyfits.Column(name='Alpha distribution', array=samples[:,1], format='D')
                c3_posterior = pyfits.Column(name='Ec distribution', array=samples[:,2], format='D')
                c4_posterior = pyfits.Column(name='b distribution', array=samples[:,3], format='D')
                table_hdu_posterior = pyfits.BinTableHDU.from_columns([c1_posterior, c2_posterior, c3_posterior, c4_posterior])
            
            elif self.comboBox_MCMC.currentText() == "PLEC_bfix":
                best_fit_model = PLEC_bfix(theta_max, x)
                N0 = np.quantile(samples[:,0],q=[0.16,0.5,0.84])
                Alpha = np.quantile(samples[:,1],q=[0.16,0.5,0.84])
                Ec = np.quantile(samples[:,2],q=[0.16,0.5,0.84])
                add_results.write(f"N0 (log scale): {N0[1]} - {N0[1]-N0[0]} + {N0[2]-N0[1]}\n")
                add_results.write(f"Alpha: {Alpha[1]} - {Alpha[1]-Alpha[0]} + {Alpha[2]-Alpha[1]}\n")
                add_results.write(f"Ec: {Ec[1]} - {Ec[1]-Ec[0]} + {Ec[2]-Ec[1]}\n")
                add_results.write("b: 1\n")
                c1 = np.array(["N0 (log scale)","Alpha","Ec","b","Ep=2*Emin (log scale)"])
                c2 = np.array([N0[1],Alpha[1],Ec[1],1,np.log10(2*self.Emin)])
                c3 = np.array([N0[1]-N0[0],Alpha[1]-Alpha[0],Ec[1]-Ec[0],0,0])
                c4 = np.array([N0[2]-N0[1],Alpha[2]-Alpha[1],Ec[2]-Ec[1],0,0])
                c1 = pyfits.Column(name='Parameter', array=c1, format='22A')
                c2 = pyfits.Column(name='Value', array=c2, format='D')
                c3 = pyfits.Column(name='error_minus', array=c3, format='D')
                c4 = pyfits.Column(name='error_plus', array=c4, format='D')
                table_hdu = pyfits.BinTableHDU.from_columns([c1, c2, c3, c4])
                c1_posterior = pyfits.Column(name='N0 distribution (log scale)', array=samples[:,0], format='D')
                c2_posterior = pyfits.Column(name='Alpha distribution', array=samples[:,1], format='D')
                c3_posterior = pyfits.Column(name='Ec distribution', array=samples[:,2], format='D')
                table_hdu_posterior = pyfits.BinTableHDU.from_columns([c1_posterior, c2_posterior, c3_posterior])

            elif self.comboBox_MCMC.currentText() == "PowerLaw":
                best_fit_model = PowerLaw(theta_max, x)
                N0 = np.quantile(samples[:,0],q=[0.16,0.5,0.84])
                Alpha = np.quantile(samples[:,1],q=[0.16,0.5,0.84])
                add_results.write(f"N0 (log scale): {N0[1]} - {N0[1]-N0[0]} + {N0[2]-N0[1]}\n")
                add_results.write(f"Alpha: {Alpha[1]} - {Alpha[1]-Alpha[0]} + {Alpha[2]-Alpha[1]}\n")
                c1 = np.array(["N0 (log scale)","Alpha","Ep=2*Emin (log scale)"])
                c2 = np.array([N0[1],Alpha[1],np.log10(2*self.Emin)])
                c3 = np.array([N0[1]-N0[0],Alpha[1]-Alpha[0],0])
                c4 = np.array([N0[2]-N0[1],Alpha[2]-Alpha[1],0])
                c1 = pyfits.Column(name='Parameter', array=c1, format='22A')
                c2 = pyfits.Column(name='Value', array=c2, format='D')
                c3 = pyfits.Column(name='error_minus', array=c3, format='D')
                c4 = pyfits.Column(name='error_plus', array=c4, format='D')
                table_hdu = pyfits.BinTableHDU.from_columns([c1, c2, c3, c4])
                c1_posterior = pyfits.Column(name='N0 distribution (log scale)', array=samples[:,0], format='D')
                c2_posterior = pyfits.Column(name='Alpha distribution', array=samples[:,1], format='D')
                table_hdu_posterior = pyfits.BinTableHDU.from_columns([c1_posterior, c2_posterior])

            add_results.close()
            hdul.append(table_hdu)
            hdul[4].name = "MCMC Parameters"
            hdul.append(table_hdu_posterior)
            hdul[5].name = "MCMC Posterior dist."
            if self.redshift > 0.0 and self.white_box_VHE.text().split('.')[-1] == 'fits':
                hdul.append(VHE_table_EBL)
                hdul[6].name = "VHE data corrected for EBL"

            hdul.writeto(SED_file,overwrite=True)
            hdul.close()


            # Corner plot:
            if self.comboBox_MCMC.currentText() == "LogPar":
                labels = ["N0", "alpha", "beta"]
            elif self.comboBox_MCMC.currentText() == "LogPar2":
                labels = ["Sp_log", "alpha", "Ep_log"]
            elif self.comboBox_MCMC.currentText() == "PLEC":
                labels = ["N0", "alpha", "Ec", "b"]
            elif self.comboBox_MCMC.currentText() == "PLEC_bfix":
                labels = ["N0", "alpha", "Ec"]
            elif self.comboBox_MCMC.currentText() == "PowerLaw":
                labels = ["N0", "alpha"]

            corner.corner(
                samples,
                show_titles=True,
                labels=labels,
                plot_datapoints=True,
                quantiles=[0.16, 0.5, 0.84],
            )

            plt.savefig(self.OutputDir+'Quickplot_MCMC_SED_pars.png',bbox_inches='tight')

            # SED plot:
            f = plt.figure(figsize=(6,5),dpi=250)
            ax = f.add_subplot(1,1,1)
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
            ax.tick_params(which='major', length=5, direction='in')
            ax.tick_params(which='minor', length=2.5, direction='in',bottom=True, top=True, left=True, right=True)
            ax.tick_params(bottom=True, top=True, left=True, right=True)
            ax.grid(linestyle=':',which='both')

            alpha_original_data = 1
            color_original_data = "C0"
            color_original_VHE_data = "C1"
            if self.redshift > 0.0:
                plt.errorbar(Energy_SED, dnde_data_points_deabsorbed*(Energy_SED**2), xerr=Energy_err_SED, yerr=dnde_error_deabsorbed*(Energy_SED**2), color="C0", markeredgecolor="black", ecolor="black", zorder = 130, fmt="o", label="Fermi-LAT EBL corrected")
                plt.errorbar(self.Energy_uplims, e2dnde_uplims_deabsorbed, xerr=self.xerr_uplims, yerr=yerr_uplims_deabsorbed, uplims=True, color="C0", markeredgecolor="black", ecolor="black", zorder = 130, fmt="o")
                if self.include_VHE:
                    plt.errorbar(Energy_SED[-N_bins_VHE:], dnde_data_points_deabsorbed[-N_bins_VHE:]*(Energy_SED[-N_bins_VHE:]**2), xerr=[Energy_err_SED[0][-N_bins_VHE:], Energy_err_SED[1][-N_bins_VHE:]], yerr=dnde_error_deabsorbed[-N_bins_VHE:]*(Energy_SED[-N_bins_VHE:]**2), color="C1", markeredgecolor="black", ecolor="black", zorder = 130, fmt="o", label="VHE EBL corrected")
                    plt.errorbar(self.Energy_uplims[-N_bins_VHE_UL:], e2dnde_uplims_deabsorbed[-N_bins_VHE_UL:], xerr=[self.xerr_uplims[0][-N_bins_VHE_UL:], self.xerr_uplims[1][-N_bins_VHE_UL:]], yerr=yerr_uplims_deabsorbed[-N_bins_VHE_UL:], uplims=True, color="C1", markeredgecolor="black", ecolor="black", zorder = 130, fmt="o")
                
                alpha_original_data = 0.4
                #color_original_data = "gray"
                #color_original_VHE_data = "gray"

            plt.errorbar(Energy_SED, dnde_SED*(Energy_SED**2), xerr=Energy_err_SED, yerr=dnde_err_SED*(Energy_SED**2), color=color_original_data, alpha = alpha_original_data, zorder = 120, fmt="o", label="Fermi-LAT")
            if self.include_VHE:
                plt.errorbar(Energy_SED[-N_bins_VHE:], dnde_SED[-N_bins_VHE:]*(Energy_SED[-N_bins_VHE:]**2), xerr=[Energy_err_SED[0][-N_bins_VHE:], Energy_err_SED[1][-N_bins_VHE:]], yerr=dnde_err_SED[-N_bins_VHE:]*(Energy_SED[-N_bins_VHE:]**2), color=color_original_VHE_data, alpha = alpha_original_data, zorder = 120, fmt="o", label="VHE instrument")
            
            plt.errorbar(self.Energy_uplims, self.e2dnde_uplims, xerr=self.xerr_uplims, yerr=self.yerr_uplims, uplims=True, color=color_original_data, alpha = alpha_original_data, zorder = 120, fmt="o")
            if self.include_VHE:
                plt.errorbar(self.Energy_uplims[-N_bins_VHE_UL:], self.e2dnde_uplims[-N_bins_VHE_UL:], xerr=[self.xerr_uplims[0][-N_bins_VHE_UL:], self.xerr_uplims[1][-N_bins_VHE_UL:]], yerr=self.yerr_uplims[-N_bins_VHE_UL:], uplims=True, color=color_original_VHE_data, alpha = alpha_original_data, zorder = 120, fmt="o")

            plotter(sampler,x)

            plt.plot(self.E, self.dnde*(self.E**2), color="gray",alpha=0.9, zorder = 119, label="Fermipy fit")
            if self.comboBox_MCMC.currentText() == "LogPar2":
                plt.plot(10**x, 10 ** best_fit_model, color="black", zorder = 119, label="Highest Likelihood MCMC")
            else:
                plt.plot(10**x, 10 ** (2 * x + best_fit_model), color="black", zorder = 119, label="Highest Likelihood MCMC")

            plt.xscale("log")
            plt.yscale("log")
            plt.xlabel("Energy [MeV]")
            plt.ylabel("E$^2dN/dE$ [MeV cm$^{-2}$ s$^{-1}$]")
            plt.title(self.sourcename+' - SED - MCMC - '+self.comboBox_MCMC.currentText())
            plt.grid(which="both", linestyle=":")
                
            if len(dnde_SED) > 0:
                ymax = 2*((dnde_SED + dnde_err_SED)*(Energy_SED**2)).max()
                ymin = 0.5*((dnde_SED - dnde_err_SED)*(Energy_SED**2)).min()
                if ymax > 4*(dnde_SED*(Energy_SED**2)).max():
                    ymax = 4*(dnde_SED*(Energy_SED**2)).max()
                
                if self.redshift > 0.0:
                    EBL_corrected_max = 2*np.max(dnde_data_points_deabsorbed*(Energy_SED**2) + dnde_error_deabsorbed*(Energy_SED**2))
                    if EBL_corrected_max > ymax:
                        ymax = EBL_corrected_max

                if ymin < (dnde_SED*(Energy_SED**2)).min()/5.0:
                    ymin = (dnde_SED*(Energy_SED**2)).min()/5.0
                    
                if len(self.e2dnde_uplims) > 0:
                    yaux = self.e2dnde_uplims.max()
                    if yaux > ymax:
                        ymax = 2*yaux
                        
                    yaux = self.e2dnde_uplims.min()
                    if yaux < ymin:
                        ymin = 0.5*yaux
            else:
                ymax = 2*self.e2dnde_uplims.max()
                ymin = 0.5*self.e2dnde_uplims.min()
                
            
            plt.ylim(ymin,ymax)
            plt.xlim(0.8 * 10**x.min(),1.2 * 10**x.max())
            plt.legend(fontsize=11)
            plt.tight_layout()
            plt.savefig(self.OutputDir+'Quickplot_SED_MCMC.'+self.comboBox_output_format.currentText(),bbox_inches='tight')


    def compute_Extension(self):

        """This function searches for extended gamma-ray emission from the target. It can use a disk or a 2D Gaussian model.

        Parameters
        ----------
        self: instance of the class Ui_mainWindow
            This parameter contains all the variables read from the Graphical interface.

        Returns
        -------
        TARGET_NAME_ext.fits
            Table containing the extendend emission fit results. File is saved in the output directory read from the graphical interface.
        TARGET_NAME_ext.npy:
            Same as above but in the npy format.
        Quickplot_extension.png or Quickplot_extension.pdf:
            A plot for quick visualization of the extension radius. File is saved in the output directory read from the graphical interface.

        """

        output_format = self.comboBox_output_format.currentText()

        #Extension:    
        if self.checkBox_extension.isChecked():
            self.gta.config['extension']['width_min'] = 0.01
            if self.radioButton_disk.isChecked(): 
                exten = self.gta.extension(self.sourcename,width=np.linspace(0.01,self.doubleSpinBox_extension_max_size.value(),20).tolist(), spatial_model='RadialDisk')
            else:
                exten = self.gta.extension(self.sourcename,width=np.linspace(0.01,self.doubleSpinBox_extension_max_size.value(),20).tolist(), spatial_model='RadialGaussian')
                
            self.gta.write_roi(self.sourcename+'_extension')            
            f = plt.figure(figsize=(6,5),dpi=250)
            ax = f.add_subplot(1,1,1)
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
            ax.tick_params(which='major', length=5, direction='in')
            ax.tick_params(which='minor', length=2.5, direction='in',bottom=True, top=True, left=True, right=True)
            ax.tick_params(bottom=True, top=True, left=True, right=True)
            ax.grid(linestyle=':',which='both')
            
            plt.plot(exten['width'],exten['dloglike'],marker='o', markeredgecolor='black')
            plt.grid(linestyle=':',which='both')
            plt.gca().set_xlabel('Radius [$^{\circ}$]')
            plt.gca().set_ylabel('Delta Log-Likelihood')
            plt.gca().axvline(exten['ext'], color='black', label='Extension radius')
            plt.gca().axvspan(exten['ext']-exten['ext_err_lo'],exten['ext']+exten['ext_err_hi'], alpha=0.2,color='b')

            plt.annotate(' TS$_{\mathrm{ext}}$ = %.2f\n R$_{68}$ = %.3f $\pm$ %.3f'%
                (exten['ts_ext'],exten['ext'],exten['ext_err']),xy=(0.05,0.05),xycoords='axes fraction')
            plt.gca().legend(frameon=False)
            plt.title(self.sourcename+' - Extension')
            plt.tight_layout()
            plt.savefig(self.OutputDir+'Quickplot_extension.'+output_format,bbox_inches='tight')


    def compute_LC(self):
        
        """
        This function calls fermipy to compute the target light curve (under request by the user).
        
        Parameters
        ----------
        self: instance of the class Ui_mainWindow
            This parameter contains all the variables read from the Graphical interface.

        Returns
        -------
        TARGET_NAME_lightcurve.fits and TARGET_NAME_lightcurve.npy:
            Data files with all the information regarding the energy-flux and photon-flux light curves.
            
        """
        
        self.Compute_LC = True   

        if os.path.exists(self.OutputDir+"Light_curve_001"):
            last_LC_directory = np.sort(glob.glob(self.OutputDir+"Light_curve_*"))[-1]
            last_number_of_bins = np.sort(glob.glob(last_LC_directory+"/lightcurve_*"))
            if len(last_number_of_bins) == self.spinBox_LC_N_time_bins.value():
                self.Compute_LC = False


            
        if self.checkBox_LC.isChecked() and self.Compute_LC:
            #Running the LC in parallel cores is possible only in Linux systems:
            if OS_name != "Darwin":
                self.gta.lightcurve(self.sourcename, nbins=self.spinBox_LC_N_time_bins.value(), free_radius=self.roiwidth/2,use_local_ltcube=True, use_scaled_srcmap=True, free_params=['norm','shape'], shape_ts_threshold=9, multithread=True, nthread=self.spinBox_N_cores_LC.value())
            else:
                self.gta.lightcurve(self.sourcename, nbins=self.spinBox_LC_N_time_bins.value(), free_radius=self.roiwidth/2,use_local_ltcube=True, use_scaled_srcmap=True, free_params=['norm','shape'], shape_ts_threshold=9, multithread=False)
            
            if not os.path.exists(self.OutputDir+"Light_curve_001"):
                os.mkdir(self.OutputDir+"Light_curve_001")
                os.system(f"mv {self.OutputDir}*lightcurve* {self.OutputDir}Light_curve_001")
            else:
                last_LC = int(np.sort(glob.glob(self.OutputDir+"Light_curve_*"))[-1].split("_")[-1]) + 1  # Taking the number corresponding to the last light curve computed and adding 1
                last_LC = "{:03d}".format(last_LC)
                os.mkdir(self.OutputDir+f"Light_curve_{last_LC}")
                os.system(f"mv {self.OutputDir}*lightcurve* {self.OutputDir}Light_curve_{last_LC}")


    def compute_LC_adaptive(self):
        
        """
        This function calls fermipy to compute the target light curve with adaptive binning.
        
        Parameters
        ----------
        self: instance of the class Ui_mainWindow
            This parameter contains all the variables read from the Graphical interface.

        Returns
        -------
        Adaptive-binning_lightcurve.fits:
            Data file with all the information regarding the energy-flux and photon-flux adaptive-binning light curves.
            
        """
        
        
        if self.checkBox_adaptive_binning.isChecked() and self.checkBox_LC.isChecked():
            try:
                TS_Threshold = float(self.white_box_TS_threshold.text())
                if TS_Threshold > 0.0:
                    self.adaptive = True
            except:
                print("Invalid TS threshold value.")
        else:
            self.adaptive = False


        if self.adaptive:

            Adaptive_binning_LC_tables = {}

            if os.path.exists(self.OutputDir+"Adaptive-binning_light_curve_001"):
                last_LC_directory = np.sort(glob.glob(self.OutputDir+"Adaptive-binning_light_curve_*"))[-1]
            else:
                last_LC_directory = np.sort(glob.glob(self.OutputDir+"Light_curve_*"))[-1]

            last_LC_bins_directories = np.sort(glob.glob(last_LC_directory+"/lightcurve_*"))

            LC_last_data_file = glob.glob(last_LC_directory+"/*_lightcurve.fits")[0]
            hdul = pyfits.open(LC_last_data_file)
            TS_at_each_time_bin = hdul[1].data["ts"]
            High_TS_bins = np.where(TS_at_each_time_bin > 2*TS_Threshold)[0]
            if len(High_TS_bins) > 0:
                for high_TS_index in High_TS_bins:
                    number_of_subBins = int(TS_at_each_time_bin[high_TS_index]/TS_Threshold)
                    if number_of_subBins > 1:
                        os.chdir(last_LC_bins_directories[high_TS_index])  # Here we enter in the directory of each LC bin with TS > 2*TS_Threshold
                        stream = open("./config.yaml", 'r')
                        data = yaml.load(stream,Loader)
                        data["fileio"]["workdir"] = "./"
                        data["fileio"]["outdir"] = "./Output"
                        data["fileio"]["logfile"] = "./Output"
                        data["components"][0]["gtlike"]["srcmap_base"] = "./srcmap_00.fits"
                        data["components"][0]["gtlike"]["bexpmap_roi_base"] = "./bexpmap_roi_00.fits"
                        data["components"][0]["gtlike"]["bexpmap_base"] = "./bexpmap_00.fits"
                        data["components"][0]["data"]["evfile"] = "./ft1_00.fits"
                        with open("./config.yaml", 'w') as yaml_file:
                            yaml_file.write( yaml.dump(data, default_flow_style=False))

                        gta = GTAnalysis("./config.yaml",logging={'verbosity': 3})
                        gta.setup()
                        if OS_name != "Darwin":
                            gta.lightcurve(self.sourcename, nbins=number_of_subBins, free_radius=self.roiwidth/2,use_local_ltcube=True, 
                                                     use_scaled_srcmap=True, free_params=['norm','shape'], shape_ts_threshold=9, multithread=True, nthread=self.spinBox_N_cores_LC.value())
                        else:
                            gta.lightcurve(self.sourcename, nbins=number_of_subBins, free_radius=self.roiwidth/2,use_local_ltcube=True, use_scaled_srcmap=True, free_params=['norm','shape'], shape_ts_threshold=9, multithread=False)
                        
                        os.chdir(Working_directory)

                # Here we copy/move the adaptive bins to the Adaptive-binning_light_curve directory, leaving a copy of the standard LC files in the Light_curve_XXX directory.
                if not os.path.exists(self.OutputDir+"Adaptive-binning_light_curve_001"):
                    os.mkdir(self.OutputDir+"Adaptive-binning_light_curve_001")
                    last_adaptive_LC = "{:03d}".format(1)
                    for n,new_lc_bins in enumerate(last_LC_bins_directories):
                        if n in High_TS_bins:
                            os.system(f"mv {new_lc_bins}/Output/lightcurve_* {self.OutputDir}Adaptive-binning_light_curve_001")
                            local_lc_table = glob.glob(f"{new_lc_bins}/Output/*_lightcurve.fits")[0]
                            Adaptive_binning_LC_tables[n] = Table.read(local_lc_table,format="fits", hdu=1)
                        else:
                            os.system(f"cp -r {new_lc_bins} {self.OutputDir}Adaptive-binning_light_curve_001")
                            local_lc_table = glob.glob(f"{last_LC_directory}/*_lightcurve.fits")[0]
                            Adaptive_binning_LC_tables[n] = Table.read(local_lc_table,format="fits", hdu=1)[n]
                else:
                    last_adaptive_LC = int(np.sort(glob.glob(self.OutputDir+"Adaptive-binning_light_curve_*"))[-1].split("_")[-1]) + 1  # Taking the number corresponding to the last adaptive binning light curve computed and adding 1
                    last_adaptive_LC = "{:03d}".format(last_adaptive_LC)
                    os.mkdir(self.OutputDir+f"Adaptive-binning_light_curve_{last_adaptive_LC}")
                    for n,new_lc_bins in enumerate(last_LC_bins_directories):
                        if n in High_TS_bins:
                            os.system(f"mv {new_lc_bins}/Output/lightcurve_* {self.OutputDir}Adaptive-binning_light_curve_{last_adaptive_LC}")
                            local_lc_table = glob.glob(f"{new_lc_bins}/Output/*_lightcurve.fits")[0]
                            Adaptive_binning_LC_tables[n] = Table.read(local_lc_table,format="fits", hdu=1)
                        else:
                            os.system(f"cp -r {new_lc_bins} {self.OutputDir}Adaptive-binning_light_curve_{last_adaptive_LC}")
                            local_lc_table = glob.glob(f"{last_LC_directory}/*_lightcurve.fits")[0]
                            Adaptive_binning_LC_tables[n] = Table.read(local_lc_table,format="fits", hdu=1)[n]
                
                
                # Saving the final table with the adaptive-binning LC:
                final_table = Adaptive_binning_LC_tables[0]
                for n in range((len(Adaptive_binning_LC_tables) - 1)):
                    final_table = vstack([final_table,Adaptive_binning_LC_tables[n+1]])

                final_table
                final_table.write(f"{self.OutputDir}Adaptive-binning_light_curve_{last_adaptive_LC}/Adaptive-binning_lightcurve.fits", format="fits",overwrite=True)
                


            else:
                print(f"No time bins with TS > 2x{TS_Threshold}. Current iteration of adaptive binning LC was not performed.")


    def plot_LCs(self,adaptive):

        """
        This function plots the latest light curves computed.
        
        Parameters
        ----------
        self: instance of the class Ui_mainWindow
            This parameter contains all the variables read from the Graphical interface.

        Returns
        -------
        Quickplot_adaptive_LC.pdf and Quickplot_adaptive_eLC.pdf:
            Plots showing the flux light curve and the energy flux light curve. File is saved in the output directory read from the graphical interface.

            
        """

        # Checking if there is something to plot:
        if adaptive:
            if os.path.exists(self.OutputDir+"Adaptive-binning_light_curve_001"):
                last_LC_directory = np.sort(glob.glob(self.OutputDir+"Adaptive-binning_light_curve_*"))[-1]
                we_can_plot = True
            else:
                we_can_plot = False
        else:
            if os.path.exists(self.OutputDir+"Light_curve_001"):
                last_LC_directory = np.sort(glob.glob(self.OutputDir+"Light_curve_*"))[-1]
                we_can_plot = True
            else:
                we_can_plot = False


        # Reading the data and plotting:
        if we_can_plot and self.checkBox_LC.isChecked():

            TSmin = 9
            output_format = self.comboBox_output_format.currentText() 

            LC_last_data_file = glob.glob(last_LC_directory+"/*_lightcurve.fits")[0]
            hdul = pyfits.open(LC_last_data_file)
            lc = hdul[1].data

            ################################################
            ########## LC energy flux
            ################################################
            f = plt.figure(figsize=(9,4),dpi=250)
            ax = f.add_subplot(1,1,1)
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
            ax.tick_params(which='major', length=5, direction='in')
            ax.tick_params(which='minor', length=2.5, direction='in',bottom=True, top=True, left=True, right=True)
            ax.tick_params(bottom=True, top=True, left=True, right=True)
            ax.grid(linestyle=':',which='both')

            
            if len(lc['eflux'][lc['ts']>TSmin]) > 0:
                scale = int(np.log10(lc['eflux'][lc['ts']>TSmin].max())) -2
            else:
                scale = int(np.log10(lc['eflux_ul95'].max())) -2
                
            tmean = (lc['tmin_mjd'] + lc['tmax_mjd'])/2


            if len(tmean[lc['ts']>TSmin]) > 9:
                time_continuum = np.linspace(np.min(lc['tmin_mjd']),np.max(lc['tmax_mjd']),10*len(lc['tmin_mjd']))
                tck_eflux = interpolate.splrep(tmean[lc['ts']>TSmin], (10**-scale)*lc['eflux'][lc['ts']>TSmin], k=3)
                tck_eflux_error = interpolate.splrep(tmean[lc['ts']>TSmin], (10**-scale)*lc['eflux_err'][lc['ts']>TSmin],k=3)

                eflux_continuum = interpolate.splev(time_continuum, tck_eflux)
                eflux_continuum_err = interpolate.splev(time_continuum, tck_eflux_error)
                ax.plot(time_continuum,eflux_continuum,color="C1", label="Spline")
                
                # Saving Spline
                col_time = pyfits.Column(name="time [MJD]",array=time_continuum,format="D",unit="MJD")
                col_eflux = pyfits.Column(name="eflux_continuum",array=eflux_continuum,format="D",unit="10^"+str(scale)+" MeV cm-2 s-1")
                col_eflux_err = pyfits.Column(name="eflux_err_continuum",array=eflux_continuum_err,format="D",unit="10^"+str(scale)+" MeV cm-2 s-1")

            plt.errorbar(tmean[lc['ts']>TSmin], (10**-scale)*lc['eflux'][lc['ts']>TSmin], xerr = [ tmean[lc['ts']>TSmin]- lc['tmin_mjd'][lc['ts']>TSmin], lc['tmax_mjd'][lc['ts']>TSmin] - tmean[lc['ts']>TSmin] ], yerr=(10**-scale)*lc['eflux_err'][lc['ts']>TSmin], markeredgecolor='black', fmt='o', capsize=4)
            plt.errorbar(tmean[lc['ts']<=TSmin], (10**-scale)*lc['eflux_ul95'][lc['ts']<=TSmin], xerr = [ tmean[lc['ts']<=TSmin]- lc['tmin_mjd'][lc['ts']<=TSmin], lc['tmax_mjd'][lc['ts']<=TSmin] - tmean[lc['ts']<=TSmin] ], yerr=5*np.ones(len(lc['eflux_err'][lc['ts']<=TSmin])), markeredgecolor='black', fmt='o', uplims=True, color='orange', capsize=4)
            plt.ylabel(r'Energy flux [$10^{'+str(scale)+'}$ MeV cm$^{-2}$ s$^{-1}$]')
            plt.xlabel('Time [MJD]')
            plt.title(self.sourcename+' - Energy light curve (free shape)')
            
            
            if len(lc['eflux'][lc['ts']>TSmin]) > 0:
                y0 = (lc['eflux'][lc['ts']>TSmin]).max()
                y1 = (lc['eflux'][lc['ts']>TSmin] + lc['eflux_err'][lc['ts']>TSmin]).max()
                if y1 > 4*y0:
                    y1 = 4*y0
                
                if len(lc['eflux_ul95'][lc['ts']<=TSmin]) > 0:
                    y2 = (lc['eflux_ul95'][lc['ts']<=TSmin]).max()
                    if y2 > y1:
                        y1 = y2
                
            else:
                y1 = (lc['eflux_ul95'][lc['ts']<=TSmin]).max()
                
            ymin = -(10**-scale)*0.1*y1
            plt.ylim(ymin,(10**-scale)*1.1*y1)
            plt.legend()
            if adaptive:
                plt.savefig(self.OutputDir+f'Quickplot_adaptive-binning_eLC_{len(lc["ts"])}_bins.'+output_format,bbox_inches='tight')
            else:
                plt.savefig(self.OutputDir+f'Quickplot_eLC_{len(lc["ts"])}_bins.'+output_format,bbox_inches='tight')
            

            ################################################
            ########## LC photon flux
            ################################################
            f = plt.figure(figsize=(9,4),dpi=250)
            ax = f.add_subplot(1,1,1)
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
            ax.tick_params(which='major', length=5, direction='in')
            ax.tick_params(which='minor', length=2.5, direction='in',bottom=True, top=True, left=True, right=True)
            ax.tick_params(bottom=True, top=True, left=True, right=True)
            ax.grid(linestyle=':',which='both')
            
            if len(lc['flux'][lc['ts']>TSmin]) > 0:
                scale = int(np.log10(lc['flux'][lc['ts']>TSmin].max())) -2
            else:
                scale = int(np.log10(lc['flux_ul95'].max())) -2
            
            if len(tmean[lc['ts']>TSmin]) > 9:
                tck_flux = interpolate.splrep(tmean[lc['ts']>TSmin], (10**-scale)*lc['flux'][lc['ts']>TSmin], k=3)
                tck_flux_error = interpolate.splrep(tmean[lc['ts']>TSmin], (10**-scale)*lc['flux_err'][lc['ts']>TSmin],k=3)

                flux_continuum = interpolate.splev(time_continuum, tck_flux)
                flux_continuum_err = interpolate.splev(time_continuum, tck_flux_error)
                ax.plot(time_continuum,flux_continuum,color="C1", label="Spline")
                
                col_flux = pyfits.Column(name="flux_continuum",array=flux_continuum,format="D",unit="10^"+str(scale)+" ph cm-2 s-1")
                col_flux_err = pyfits.Column(name="flux_err_continuum",array=flux_continuum_err,format="D",unit="10^"+str(scale)+" ph cm-2 s-1")
                all_cols = pyfits.BinTableHDU.from_columns([col_time, col_eflux, col_eflux_err, col_flux, col_flux_err])
                if len(hdul) < 3: 
                    hdul.append(all_cols)
                    hdul[2].name = "LC spline data"
                    hdul.writeto(LC_last_data_file,overwrite=True)
                    hdul.close()                

                
            plt.errorbar(tmean[lc['ts']>TSmin], (10**-scale)*lc['flux'][lc['ts']>TSmin], xerr = [ tmean[lc['ts']>TSmin]- lc['tmin_mjd'][lc['ts']>TSmin], lc['tmax_mjd'][lc['ts']>TSmin] - tmean[lc['ts']>TSmin] ], yerr=(10**-scale)*lc['flux_err'][lc['ts']>TSmin], markeredgecolor='black', fmt='o', capsize=4)
            plt.errorbar(tmean[lc['ts']<=TSmin], (10**-scale)*lc['flux_ul95'][lc['ts']<=TSmin], xerr = [ tmean[lc['ts']<=TSmin]- lc['tmin_mjd'][lc['ts']<=TSmin], lc['tmax_mjd'][lc['ts']<=TSmin] - tmean[lc['ts']<=TSmin] ], yerr=5*np.ones(len(lc['flux_err'][lc['ts']<=TSmin])), markeredgecolor='black', fmt='o', uplims=True, color='orange', capsize=4)
            plt.ylabel(r'Flux [$10^{'+str(scale)+'}$ cm$^{-2}$ s$^{-1}$]')
            plt.xlabel('Time [MJD]')
            plt.title(self.sourcename+' - Light curve (free shape)')
            
            
            if len(lc['flux'][lc['ts']>TSmin]) > 0:
                y0 = (lc['flux'][lc['ts']>TSmin]).max()
                y1 = (lc['flux'][lc['ts']>TSmin] + lc['flux_err'][lc['ts']>TSmin]).max()
                if y1 > 4*y0:
                    y1 = 4*y0
                
                if len(lc['flux_ul95'][lc['ts']<=TSmin]) > 0:
                    y2 = (lc['flux_ul95'][lc['ts']<=TSmin]).max()
                    if y2 > y1:
                        y1 = y2
                
            else:
                y1 = (lc['flux_ul95'][lc['ts']<=TSmin]).max()
            
            ymin = -(10**-scale)*0.1*y1               
            plt.ylim(ymin,(10**-scale)*1.1*y1)
            plt.legend()
            
            if adaptive:
                plt.savefig(self.OutputDir+f'Quickplot_adaptive-binning_LC_{len(lc["ts"])}_bins.'+output_format,bbox_inches='tight')
            else:
                plt.savefig(self.OutputDir+f'Quickplot_LC_{len(lc["ts"])}_bins.'+output_format,bbox_inches='tight')

    





#if __name__ == "__main__":
app = QtWidgets.QApplication(sys.argv)
mainWindow = QtWidgets.QMainWindow()
ui = Ui_mainWindow()
ui.setupUi(mainWindow)
mainWindow.show()
sys.exit(app.exec_())





