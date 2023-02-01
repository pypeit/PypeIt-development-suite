import pytest

import glob
from collections import namedtuple
import shutil
import re
import os

from qtpy.QtCore import Qt, QSettings

from pypeit.setup_gui import view, model, controller
from pypeit.tests.tstutils import data_path
from pypeit.spectrographs.util import load_spectrograph
from pypeit.metadata import PypeItMetaData
from pypeit.inputfiles import PypeItFile

# Override pytest-qt's QApplication arguments
@pytest.fixture(scope="session")
def qapp_args():
    return ["pytest", "-platform", "offscreen"]

def get_test_metadata():
    # Create a PyepItMetadata to test with
    file_pattern = data_path("b*.fits.gz")
    files = glob.glob(file_pattern)
    spec = load_spectrograph("shane_kast_blue")
    par = spec.default_pypeit_par()
    return PypeItMetaData(spec, par, files)

def get_mock_setup(metadata=None):
    MockSetup = namedtuple('MockSetup', ['par', 'user_cfg', 'fitstbl'])


    if metadata is not None:
        spec = metadata.spectrograph
        par = metadata.par
    else:
        spec = load_spectrograph("shane_kast_blue")
        par = spec.default_pypeit_par()

    user_cfg = ["[rdx]", # Basic entry in config
                "spectrograph = shane_kast_blue", 
                "[reduce]",  # Nested entries in config, with no values at top level
                "[[findobj]]", 
                "snr_thresh = 5.0",
                "[calibrations]", # Double nested entries in config
                "bpm_usebias = true",
                "[[biasframe]]",
                "frametype = bias",
                "[[[process]]]",
                "combine = median", ]

    return MockSetup(par, user_cfg, metadata)



def test_metadata_proxy(qtbot, qtmodeltester):

    metadata = get_test_metadata()
    proxy = model.PypeItMetadataProxy()

    # Test a model with data in it
    with qtbot.waitSignals([proxy.modelAboutToBeReset, proxy.modelReset], raising=True, order='strict', timeout=5000):
        proxy.setMetadata(metadata)

    # Run Qt tester to cover the basics
    qtmodeltester.check(proxy)

    # Test some things not covered by the above check
    
    # Floating point rounding    
    # Sorting and floating point rounding
    dec_column = metadata.set_pypeit_cols(write_bkg_pairs=True).index('dec')
    proxy.sort(dec_column)
    indx = proxy.createIndex(0, dec_column)
    assert proxy.data(indx, Qt.DisplayRole) == '37.324'

    # Column names and sort order
    assert proxy.headerData(dec_column, Qt.Orientation.Horizontal, Qt.DisplayRole) == 'dec'
    assert proxy.headerData(dec_column, Qt.Orientation.Vertical, Qt.DecorationRole) is None
    assert proxy.headerData(dec_column, Qt.Orientation.Vertical, Qt.DisplayRole) == str(dec_column+1)
    assert proxy.headerData(dec_column, Qt.Orientation.Horizontal, Qt.InitialSortOrderRole) == Qt.AscendingOrder
    

def test_pypeit_params_proxy(qtmodeltester):
    # Mock setup object
    mock_setup = get_mock_setup()
    proxy = model.PypeItParamsProxy(mock_setup)

    # Run Qt tester to cover the basics
    qtmodeltester.check(proxy)
    idx = proxy.index(0,1)
    assert proxy.data(idx) is None
    parent = proxy.index(0,0)
    idx = proxy.index(0,1, parent)
    assert proxy.data(idx) == "shane_kast_blue"
    parent = proxy.index(1,0) # [reduce]
    parent = proxy.index(0,0,parent) #[[findobj]]
    idx = proxy.index(0,1, parent) # snr_thresh value
    assert proxy.data(idx) == '5.0'
    parent = proxy.index(1,1)
    idx = proxy.index(0,0, parent)
    assert not idx.isValid()

    # Header data
    assert proxy.headerData(0, Qt.Orientation.Horizontal, Qt.DisplayRole) == "Setting"
    assert proxy.headerData(1, Qt.Orientation.Horizontal, Qt.DisplayRole) == "Value"
    assert proxy.headerData(1, Qt.Orientation.Vertical, Qt.DisplayRole) == 2 # Row #s on side
    assert proxy.headerData(1, Qt.Orientation.Vertical, Qt.DecorationRole) is None
    
def verify_config_name(name):
    return name == 'B'

def test_setup_config_model(qtbot, tmp_path):
    # Get a mock pypeit setup
    metadata = get_test_metadata()
    # This will set "setup" to B for every file except b1.fits.gz
    metadata.get_frame_types()
    metadata.table['setup'] = 'A'
    metadata.table['setup'][metadata.table['decker'] == '2.0 arcsec'] = 'B'
    mock_setup = get_mock_setup(metadata)
    config_dict = {'B': {'dispname': '600/4310', 'dichroic': 'd55'}}
    config_model = model.SetupConfigurationModel(mock_setup, "B", config_dict, model.ModelState.CHANGED)

    with qtbot.waitSignal(config_model.stateChanged, raising=True, check_params_cb=verify_config_name, timeout=5000):
        config_model.save(str(tmp_path))
    
    # Make sure state changed
    assert config_model.state == model.ModelState.UNCHANGED

    # Make sure we can open the created file
    pf = PypeItFile.from_file(str(tmp_path / "shane_kast_blue_B" / "shane_kast_blue_B.pypeit"))
    filenames = pf.filenames

    # Make sure one known file is there, and that b1.fits.gz is NOT
    assert data_path("b11.fits.gz") in filenames
    assert data_path("b1.fits.gz") not in filenames

def verify_failed_gen_obslog(message):
    return message.startswith("Could not read metadata")

def verify_canceled_gen_obslog(message):
    return message == "CANCEL"

def raise_cancel_exception(*args):
    raise model.OpCanceledError()

def test_pypeit_setup_proxy(qtbot, tmp_path):

    proxy = model.PypeItSetupProxy()
    logname = str(tmp_path / "setup_gui_model_test.log")
    proxy.setup_logging(logname, 2)

    assert len(proxy.available_spectrographs) != 0
    assert proxy.state == model.ModelState.NEW

    # Copy only those fits files we want to run on
    shutil.copy2(data_path("b1.fits.gz"), tmp_path)
    shutil.copy2(data_path("b27.fits.gz"), tmp_path)

    # Test a failed setup run with the wrong spectrograph
    with qtbot.waitSignal(proxy.operation_complete, raising=True, check_params_cb=verify_failed_gen_obslog, timeout=10000):
        proxy.set_spectrograph("keck_deimos")
        proxy.set_raw_data_directory(str(tmp_path))
        assert proxy.get_num_raw_files() == 2
        proxy.generate_obslog()

    # The controller normally does this reset after an error
    proxy.reset()

    # Test a canceled setup run 
    with qtbot.waitSignal(proxy.operation_complete, raising=True, check_params_cb=verify_canceled_gen_obslog, timeout=10000):
        proxy.set_spectrograph("shane_kast_blue")
        proxy.set_raw_data_directory(str(tmp_path))
        assert proxy.get_num_raw_files() == 2
        # Use the proxy's internal log watcher to trigger the canceled exception
        proxy._log_watcher.watch("test_cancel", re.compile("Adding metadata for .*$"), raise_cancel_exception)
        proxy.generate_obslog()
        proxy._log_watcher.unwatch("test_cancel")

    # The controller normally does this reset after an error
    proxy.reset()

    # Run setup on those files
    signal_list = [(proxy.spectrograph_changed, "set_spec"), 
                   (proxy.raw_data_dir_changed, "reset_raw_data"), 
                   (proxy.spectrograph_changed, "reset_spec"), 
                   (proxy.raw_data_dir_changed, "set_raw_data"), 
                   (proxy.spectrograph_changed, "set_raw_data_spec"), 
                   (proxy.operation_progress, "op_progress_file1"),
                   (proxy.operation_progress, "op_progress_file2"),
                   (proxy.configs_added, "configs_added"),
                   (proxy.operation_complete, "op_complete"),]

    with qtbot.waitSignals(signal_list, raising=True, order='strict', timeout=10000):
        proxy.set_spectrograph("shane_kast_blue")
        proxy.set_raw_data_directory(str(tmp_path))
        assert proxy.get_num_raw_files() == 2
        proxy.generate_obslog()
        assert proxy.state == model.ModelState.CHANGED
        assert proxy._pypeit_setup is not None

    # Save the setup
    proxy.save_all(str(tmp_path))
    dest_file = tmp_path / "shane_kast_blue_A" / "shane_kast_blue_A.pypeit"
    assert dest_file.exists()
    assert proxy.state == model.ModelState.UNCHANGED
    
    # Test a failed save
    dest_file.unlink()
    dest_file.mkdir(parents=True)
    assert dest_file.is_dir()

    with pytest.raises(RuntimeError, match="Failed saving setup"):
        proxy.save_config("A", str(tmp_path))

    # Save a specific setup
    dest_file.rmdir()
    dest_file.parent.rmdir()
    assert not dest_file.exists()
    proxy.save_config("A", str(tmp_path))
    assert dest_file.exists()


    # Reset the proxy
    signal_list = [(proxy.raw_data_dir_changed, "reset_raw_data"), 
                   (proxy.spectrograph_changed, "reset_spec"), 
                   (proxy.configs_deleted, "configs_deleted"),]
    with qtbot.waitSignals(signal_list, raising=True, order='strict', timeout=1000):
        proxy.reset()

    assert proxy.state == model.ModelState.NEW

    # Now open the previously created file
    signal_list = [(proxy.raw_data_dir_changed, "set_raw_data"), 
                   (proxy.spectrograph_changed, "set_spec"),
                   (proxy.configs_added,        "configs_added")]

    with qtbot.waitSignals(signal_list, raising=True, order='strict', timeout=1000):
        proxy.open_pypeit_file(str(dest_file))

    assert proxy.state == model.ModelState.UNCHANGED

    assert proxy.spectrograph == "shane_kast_blue"
    assert proxy.raw_data_directory == str(tmp_path)
    assert "A" in proxy.configs
    assert "b1.fits.gz" in proxy.configs["A"].metadata_model._metadata.table['filename']

def test_gui_basic_usage(qtbot, tmp_path, monkeypatch):
    
    # Change to tmp_path so that log goes there
    monkeypatch.chdir(tmp_path)

    # Set config location and format so that history information is saved to a local ini file
    # instead of the user's home directory (linux and mac) or the registry (windows)
    QSettings.setDefaultFormat(QSettings.Format.IniFormat)
    QSettings.setPath(QSettings.Format.IniFormat, QSettings.Scope.UserScope, str(tmp_path / "config_user.ini"))
    QSettings.setPath(QSettings.Format.IniFormat, QSettings.Scope.SystemScope, str(tmp_path / "config_system.ini"))
    
    # Start the gui with with the "offscreen" platform for headless testing.
    Args = namedtuple("Args", ['verbosity'])
    args = Args(verbosity=2)
    c = controller.SetupGUIController(args)

    main_window = c.view
    pypeit_setup_model = c.model
    qtbot.addWidget(main_window)
    main_window.show()

    # Validate initial state
    assert main_window.spectrograph.currentIndex() == -1
    assert main_window.spectrograph.isEnabled()

    assert main_window.raw_data.current_location is None
    assert not main_window.raw_data.isEnabled()

    assert main_window.outdir.current_location is None
    assert main_window.outdir.isEnabled()

    assert main_window.openButton.isEnabled()
    assert not main_window.clearButton.isEnabled()
    assert not main_window.setupButton.isEnabled()
    assert not main_window.saveTabButton.isEnabled()
    assert not main_window.save_all_button.isEnabled()

    assert main_window.setup_view.state == model.ModelState.NEW
    assert main_window.setup_view.get_tab(0).state == model.ModelState.UNCHANGED
    assert main_window.setup_view.get_tab(0).tab_name == "ObsLog"

    # Select the "keck_mosfire" spectrograph 
    keck_mosfire_offset = main_window.spectrograph.findText("keck_mosfire", Qt.MatchExactly)
    assert keck_mosfire_offset != -1
    with qtbot.waitSignal(pypeit_setup_model.spectrograph_changed, timeout=1000):
        main_window.spectrograph.setCurrentIndex(keck_mosfire_offset)

    assert pypeit_setup_model.spectrograph == "keck_mosfire"

    assert main_window.raw_data.current_location is None
    assert main_window.raw_data.isEnabled()

    # set raw data
    raw_dir = os.path.join(os.getenv('PYPEIT_DEV'), 
                           'RAW_DATA', 'keck_mosfire')
    j_multi = os.path.join(raw_dir, "J_multi")
    with qtbot.waitSignal(pypeit_setup_model.operation_complete, raising=True, timeout=10000):
        main_window.raw_data._location.setCurrentText(j_multi)
        qtbot.keyClick(main_window.raw_data._location, Qt.Key_Enter)

    
    assert main_window.raw_data.current_location == j_multi
    assert main_window.raw_data.isEnabled()

    assert main_window.openButton.isEnabled()
    assert main_window.clearButton.isEnabled()
    assert main_window.setupButton.isEnabled()
    assert not main_window.saveTabButton.isEnabled()
    assert main_window.save_all_button.isEnabled()

    assert main_window.setup_view.state == model.ModelState.CHANGED
    assert main_window.setup_view.get_tab(0).state == model.ModelState.UNCHANGED
    assert main_window.setup_view.get_tab(0).tab_name == "ObsLog"
    
    # Select Tab A
    with qtbot.waitSignal(main_window.setup_view.currentChanged, raising=True, timeout=10000):
        main_window.setup_view.setCurrentIndex(1)

    tab_a = main_window.setup_view.get_tab(1) 
    tab_b = main_window.setup_view.get_tab(2) 
    assert tab_a.name == "A"
    assert tab_a.tab_name == "*A"
    assert tab_a.state == model.ModelState.CHANGED
    assert tab_b.name == "B"
    assert tab_b.tab_name == "*B"
    assert tab_b.state == model.ModelState.CHANGED
    assert main_window.setup_view.state == model.ModelState.CHANGED

    assert main_window.saveTabButton.isEnabled()

    # Set output directory
    with qtbot.waitSignal(main_window.outdir.location_changed, raising=True, timeout=10000):
        main_window.outdir.set_current_location(str(tmp_path))
        qtbot.keyClick(main_window.outdir._location, Qt.Key_Enter)

    # Save the tab
    with qtbot.waitSignal(tab_a._model.stateChanged, raising=True, timeout=10000):
        qtbot.mouseClick(main_window.saveTabButton, Qt.MouseButton.LeftButton)

    assert tab_a.state == model.ModelState.UNCHANGED
    assert tab_b.state == model.ModelState.CHANGED
    assert main_window.setup_view.state == model.ModelState.CHANGED
    assert not main_window.saveTabButton.isEnabled()
    assert main_window.save_all_button.isEnabled()

    tab_a_file = tmp_path / "keck_mosfire_A/keck_mosfire_A.pypeit"
    assert tab_a_file.is_file()
    tab_a_file.unlink()

    # Do a "save all"
    with qtbot.waitSignal(tab_b._model.stateChanged, raising=True, timeout=10000):
        qtbot.mouseClick(main_window.save_all_button, Qt.MouseButton.LeftButton)
    
    assert tab_a.state == model.ModelState.UNCHANGED
    assert tab_b.state == model.ModelState.UNCHANGED
    assert main_window.setup_view.state == model.ModelState.UNCHANGED
    assert not main_window.saveTabButton.isEnabled()
    assert not main_window.save_all_button.isEnabled()

    tab_a_file = tmp_path / "keck_mosfire_A/keck_mosfire_A.pypeit"
    assert tab_a_file.exists()

    tab_b_file = tmp_path / "keck_mosfire_B/keck_mosfire_B.pypeit"
    assert tab_b_file.is_file()

    # clear
    with qtbot.waitSignal(pypeit_setup_model.configs_deleted, raising=True, timeout=10000):
        qtbot.mouseClick(main_window.clearButton, Qt.MouseButton.LeftButton)

    # Verify the GUI has been reset to its initial state
    assert main_window.openButton.isEnabled()
    assert not main_window.clearButton.isEnabled()
    assert not main_window.setupButton.isEnabled()
    assert not main_window.saveTabButton.isEnabled()
    assert not main_window.save_all_button.isEnabled()

    assert main_window.setup_view.state == model.ModelState.NEW
    assert main_window.setup_view.get_tab(0).state == model.ModelState.UNCHANGED
    assert main_window.setup_view.get_tab(0).tab_name == "ObsLog"
    assert main_window.setup_view.count() == 1

    # open saved
    def mock_prompt_for_file(*args, **kwargs):
        return tab_b_file

    monkeypatch.setattr(main_window, "prompt_for_file", mock_prompt_for_file)
    with qtbot.waitSignal(pypeit_setup_model.configs_added, raising=True, timeout=10000):
        qtbot.mouseClick(main_window.openButton, Qt.MouseButton.LeftButton)

    # Verify a single tab was opened and everthing is the correct state
    assert main_window.openButton.isEnabled()
    assert main_window.clearButton.isEnabled()
    assert main_window.setupButton.isEnabled()
    assert not main_window.saveTabButton.isEnabled()
    assert not main_window.save_all_button.isEnabled()

    assert main_window.setup_view.state == model.ModelState.UNCHANGED
    assert main_window.setup_view.count() == 2
    assert main_window.setup_view.currentIndex() == 0
    assert main_window.setup_view.get_tab(0).state == model.ModelState.UNCHANGED
    assert main_window.setup_view.get_tab(0).tab_name == "ObsLog"
    assert main_window.setup_view.get_tab(1).state == model.ModelState.UNCHANGED
    assert main_window.setup_view.get_tab(1).tab_name == "B"

    # re-run setup
    with qtbot.waitSignal(pypeit_setup_model.operation_complete, raising=True, timeout=10000):
        qtbot.mouseClick(main_window.setupButton, Qt.MouseButton.LeftButton)

    # Verify a both and "A" and "B" tab was created and everthing is the correct state
    assert main_window.openButton.isEnabled()
    assert main_window.clearButton.isEnabled()
    assert main_window.setupButton.isEnabled()
    assert not main_window.saveTabButton.isEnabled()
    assert main_window.save_all_button.isEnabled()

    assert main_window.setup_view.state == model.ModelState.CHANGED
    assert main_window.setup_view.count() == 3
    assert main_window.setup_view.currentIndex() == 0
    assert main_window.setup_view.get_tab(0).state == model.ModelState.UNCHANGED
    assert main_window.setup_view.get_tab(0).tab_name == "ObsLog"
    assert main_window.setup_view.get_tab(1).state == model.ModelState.CHANGED
    assert main_window.setup_view.get_tab(1).tab_name == "*A"
    assert main_window.setup_view.get_tab(2).state == model.ModelState.CHANGED
    assert main_window.setup_view.get_tab(2).tab_name == "*B"

    # verify history
    settings = QSettings()
    settings.beginGroup("RawDataDirectory")
    history = settings.value('History')
    assert history == [str(j_multi)]
    settings = QSettings()
    settings.beginGroup("OutputDirectory")
    history = settings.value('History')
    assert history == [str(tmp_path)]
    settings = QSettings()
    settings.beginGroup("OpenPypeItFile")
    history = settings.value('History')
    assert history == [str(tab_b_file.parent)]
