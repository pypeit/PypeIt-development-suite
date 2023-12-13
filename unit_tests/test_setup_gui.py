import pytest

import glob
from collections import namedtuple
import shutil
import re
import os

from qtpy.QtCore import Qt, QSettings,QObject
from qtpy.QtWidgets import QFileDialog

from pypeit.setup_gui import view, model, controller
from pypeit.setup_gui.view import FileDialog
from pypeit.tests.tstutils import data_path
from pypeit.spectrographs.util import load_spectrograph
from pypeit.metadata import PypeItMetaData
from pypeit.inputfiles import PypeItFile
from pypeit import msgs

# Override pytest-qt's QApplication arguments
@pytest.fixture(scope="session")
def qapp_args():
    return ["pytest", "-platform", "offscreen"]

def get_test_metadata():
    # Create a PypeItMetadata to test with
    file_pattern = data_path("b*.fits.gz")
    files = glob.glob(file_pattern)
    spec = load_spectrograph("shane_kast_blue")
    par = spec.default_pypeit_par()
    return PypeItMetaData(spec, par, files)

def get_mock_setup(metadata=None):
    MockSetup = namedtuple('MockSetup', ['par', 'user_cfg', 'fitstbl', 'spectrograph'])


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

    return MockSetup(par, user_cfg, metadata, spec)


class MockFileDialog():
    mock_response = view.DialogResponses.CANCEL
    mock_location = None
    call_count = 0
    def __init__(parent, *args, **kwargs):
        pass

    def show(self):
        self.selected_path = self.mock_location
        self.__class__.call_count += 1
        print(f"MockFileDialog Returning mock response: {self.mock_response}")
        return self.mock_response


def test_metadata_model(qtbot, qtmodeltester):

    metadata = get_test_metadata()
    metadata_model = model.PypeItMetadataModel(metadata=None)
    spec_name = metadata.spectrograph.name
    spec = metadata.spectrograph

    # Sort by dec so we get a consistent order
    dec_column = metadata.set_pypeit_cols(write_bkg_pairs=True).index('dec')
    metadata.sort('dec')
    
    # Test setting the metadata
    with qtbot.waitSignals([metadata_model.modelAboutToBeReset, metadata_model.modelReset], raising=True, order='strict', timeout=5000):
        metadata_model.setMetadata(metadata)

    # Run Qt tester to cover the basics of what Qt expects a table model to do
    qtmodeltester.check(metadata_model)

    # Test some things not covered by the above check

    # Floating point rounding    
    indx = metadata_model.createIndex(0, dec_column)
    assert metadata_model.data(indx, Qt.DisplayRole) == '24.996'

    # Column names and sort order
    assert metadata_model.headerData(dec_column, Qt.Orientation.Horizontal, Qt.DisplayRole) == 'dec'
    assert metadata_model.headerData(dec_column, Qt.Orientation.Vertical, Qt.DecorationRole) is None
    
    # TODO: More tests?
    # Test removing metadata rows
    # Test adding metadata rows
    # Test paste
    # Test copyForRows
    # Test CopyFromConfig?
    # Test comment/uncomment/iscommented

def test_pypeit_params_proxy(qtmodeltester):

    # Mock setup object
    mock_setup = get_mock_setup()
    
    proxy = model.PypeItParamsProxy(mock_setup.par, mock_setup.user_cfg)

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
    
def verify_state_change(name, state):
    return name == 'B' and state == model.ModelState.UNCHANGED

def test_pypeit_file_model(qtbot, tmp_path):

    # Get a mock pypeit setup
    metadata = get_test_metadata()
    # This will set "setup" to B for every file except b1.fits.gz
    metadata.get_frame_types()
    metadata.table['setup'] = 'A'
    metadata.table['setup'][metadata.table['decker'] == '2.0 arcsec'] = 'B'
    metadata_model = model.PypeItMetadataModel(metadata=metadata)
    mock_setup = get_mock_setup(metadata)
    params_model = model.PypeItParamsProxy(mock_setup.par, mock_setup.user_cfg)
    config_dict = {'B': {'dispname': '600/4310', 'dichroic': 'd55'}}
    file_model = model.PypeItFileModel(metadata_model = metadata_model, name_stem="B", config_dict=config_dict, state=model.ModelState.CHANGED, params_model=params_model)
    file_model.save_location = str(tmp_path)
    with qtbot.waitSignal(file_model.stateChanged, raising=True, check_params_cb=verify_state_change, timeout=5000):
        file_model.save()
    
    # Make sure state changed
    assert file_model.state == model.ModelState.UNCHANGED

    # Make sure we can open the created file
    pf = PypeItFile.from_file(str(tmp_path / "shane_kast_blue_B.pypeit"))
    filenames = pf.filenames

    # Make sure one setup A and one setup B is there. This would probably be invalid
    # but the GUI trusts the user to know what they're doing
    assert data_path("b11.fits.gz") in filenames
    assert data_path("b1.fits.gz") in filenames

    # TODO
    # test commented out files being saved,reloaded correctly
    # test removed files not being saved
    # test comments in configuration section being saved

def verify_failed_setup(message):
    return message.startswith("Could not read metadata")

def verify_canceled_setup(message):
    return message == "CANCEL"

def raise_cancel_exception(*args):
    raise model.OpCanceledError()

def test_pypeit_setup_model(qtbot, tmp_path):

    setup_model = model.PypeItSetupModel()
    logname = str(tmp_path / "setup_gui_model_test.log")
    setup_model.setup_logging(logname, 2)

    assert setup_model.state == model.ModelState.NEW

    # Copy only those fits files we want to run on
    shutil.copy2(data_path("b1.fits.gz"), tmp_path)
    shutil.copy2(data_path("b27.fits.gz"), tmp_path)

    # Test a failed setup run with the wrong spectrograph
    with qtbot.waitSignal(setup_model.operation_complete, raising=True, check_params_cb=verify_failed_setup, timeout=10000):
        setup_model.set_spectrograph("keck_deimos")
        setup_model.add_raw_data_directory(str(tmp_path))
        assert setup_model.scan_raw_data_directories() == 2
        setup_model.run_setup()

    # The controller normally does this reset after an error
    setup_model.reset()

    # Test a canceled setup run 
    with qtbot.waitSignal(setup_model.operation_complete, raising=True, check_params_cb=verify_canceled_setup, timeout=10000):
        setup_model.set_spectrograph("shane_kast_blue")
        setup_model.add_raw_data_directory(str(tmp_path))
        assert setup_model.scan_raw_data_directories() == 2
        # Use the proxy's internal log watcher to trigger the canceled exception
        setup_model.log_buffer.watch("test_cancel", re.compile("Adding metadata for .*$"), raise_cancel_exception)
        setup_model.run_setup()
        setup_model.log_buffer.unwatch("test_cancel")

    # The controller normally does this reset after an error
    setup_model.reset()

    # Run setup on those files
    signal_list = [(setup_model.spectrograph_changed, "set_spec"), 
                   (setup_model.raw_data_dirs_changed, "reset_raw_data"), 
                   (setup_model.operation_progress, "op_progress_file1"),
                   (setup_model.operation_progress, "op_progress_file2"),
                   (setup_model.configs_added, "configs_added"),
                   (setup_model.operation_complete, "op_complete"),]

    with qtbot.waitSignals(signal_list, raising=True, order='strict', timeout=10000):
        setup_model.set_spectrograph("shane_kast_blue")
        setup_model.add_raw_data_directory(str(tmp_path))
        assert setup_model.scan_raw_data_directories() == 2
        setup_model.run_setup()
        assert setup_model.state == model.ModelState.CHANGED
        assert setup_model._pypeit_setup is not None

    # Save the setup
    pypeit_file_model = setup_model.pypeit_files["A"]
    pypeit_file_model.save_location = str(tmp_path)
    pypeit_file_model.save()
    dest_file = tmp_path / "shane_kast_blue_A" / "shane_kast_blue_A.pypeit"
    assert dest_file.exists()
    assert setup_model.state == model.ModelState.UNCHANGED

    # Reset the proxy
    signal_list = [(setup_model.raw_data_dirs_changed, "reset_raw_data"), 
                   (setup_model.spectrograph_changed, "reset_spec"), 
                   (setup_model.configs_deleted, "configs_deleted"),]
    with qtbot.waitSignals(signal_list, raising=True, order='strict', timeout=1000):
        setup_model.reset()
    assert setup_model.state == model.ModelState.NEW
    assert setup_model.spectrograph == None
    assert len(setup_model.raw_data_directories) == 0
    assert len(setup_model.pypeit_files.values()) == 0

    # Now open the previously created file
    signal_list = [(setup_model.raw_data_dirs_changed, "set_raw_data"), 
                   (setup_model.spectrograph_changed, "set_spec"),
                   (setup_model.configs_added,        "configs_added")]

    with qtbot.waitSignals(signal_list, raising=True, order='strict', timeout=1000):
        setup_model.open_pypeit_file(str(dest_file))

    assert setup_model.state == model.ModelState.UNCHANGED

    assert setup_model.spectrograph == "shane_kast_blue"
    assert str(tmp_path) in setup_model.raw_data_directories
    assert "A" in setup_model.pypeit_files
    assert "b1.fits.gz" in setup_model.pypeit_files["A"].metadata_model._metadata.table['filename']


    # Test a failed save
    dest_file.unlink()
    dest_file.mkdir(parents=True)
    assert dest_file.is_dir()

    with pytest.raises(RuntimeError, match="Failed saving setup"):
        pypeit_file_model.save()


    # Reset the proxy
    signal_list = [(setup_model.raw_data_dirs_changed, "reset_raw_data"), 
                   (setup_model.spectrograph_changed, "reset_spec"), 
                   (setup_model.configs_deleted, "configs_deleted"),]
    with qtbot.waitSignals(signal_list, raising=True, order='strict', timeout=1000):
        setup_model.reset()

    assert setup_model.state == model.ModelState.NEW


def setup_offscreen_gui(tmp_path, monkeypatch, qtbot, verbosity=2, logfile="test_log.log"):
    """Helper function to setup the gui to run in offscreen mode"""
    # Change to tmp_path so that log goes there
    monkeypatch.chdir(tmp_path)

    # Set config location and format so that history information is saved to a local ini file
    # instead of the user's home directory (linux and mac) or the registry (windows)
    QSettings.setDefaultFormat(QSettings.Format.IniFormat)
    QSettings.setPath(QSettings.Format.IniFormat, QSettings.Scope.UserScope, str(tmp_path / "config_user.ini"))
    QSettings.setPath(QSettings.Format.IniFormat, QSettings.Scope.SystemScope, str(tmp_path / "config_system.ini"))
    
    # Start the gui with with the "offscreen" platform for headless testing.
    Args = namedtuple("Args", ['verbosity','logfile'])
    args = Args(verbosity=verbosity,logfile=logfile)
    c = controller.SetupGUIController(args)
    main_window = c.view

    # Patch display error so we don't get hung on it's modal dialog box
    def mock_display_error(message):
        raise ValueError(message)
    monkeypatch.setattr(main_window, "display_error", mock_display_error)
    qtbot.addWidget(main_window)
    main_window.show()
    return c, main_window

def run_setup(spectrograph, test_setup, main_window, qtbot):
    """Run setup in an already started GUI. This is a helper method intended to set things up
    for tests that focus on functionality that require a populated GUI."""
    
    spectrograph_widget = main_window.setup_view._obs_log_tab.spectrograph
    raw_data_paths_widget = main_window.setup_view._obs_log_tab.raw_data_paths._location

    with qtbot.waitSignal(main_window.model.spectrograph_changed, timeout=1000):
        qtbot.keyClicks(spectrograph_widget, "keck_mosfire")
        qtbot.keyClick(spectrograph_widget, Qt.Key.Key_Enter)

    # set raw data
    raw_dir = os.path.join(os.getenv('PYPEIT_DEV'), 
                           'RAW_DATA', spectrograph, test_setup)

    with qtbot.waitSignal(main_window.model.raw_data_dirs_changed, raising=True, timeout=1000):
        raw_data_paths_widget.setCurrentText(raw_dir)
        qtbot.keyClick(raw_data_paths_widget, Qt.Key_Enter)

    # Click the Run Setup button
    with qtbot.waitSignal(main_window.model.operation_complete, raising=True, timeout=10000):
        qtbot.mouseClick(main_window.setupButton, Qt.MouseButton.LeftButton)

def test_run_setup(qtbot, tmp_path, monkeypatch):
    """Test the basic process of setting a spectrograph, a raw data directory, and
    running setup."""
    c, main_window = setup_offscreen_gui(tmp_path, monkeypatch, qtbot)    
    pypeit_setup_model = c.model

    # Validate initial state
    obs_log_tab = main_window.setup_view._obs_log_tab
    assert obs_log_tab.spectrograph.currentIndex() == -1
    assert obs_log_tab.spectrograph.isEnabled()

    assert len(obs_log_tab.raw_data_paths.locations) == 0
    assert not obs_log_tab.raw_data_paths.isEnabled()


    assert main_window.openButton.isEnabled()
    assert not main_window.clearButton.isEnabled()
    assert not main_window.setupButton.isEnabled()
    assert not main_window.saveTabButton.isEnabled()
    assert not main_window.saveAllButton.isEnabled()

    assert main_window.setup_view.state == model.ModelState.NEW
    assert main_window.setup_view.currentWidget().state == model.ModelState.UNCHANGED
    assert main_window.setup_view.currentWidget().tab_name == "ObsLog"

    # Select the "keck_mosfire" spectrograph 
    with qtbot.waitSignal(pypeit_setup_model.spectrograph_changed, timeout=1000):
        qtbot.keyClicks(obs_log_tab.spectrograph, "keck_mosfire")
        qtbot.keyClick(obs_log_tab.spectrograph, Qt.Key.Key_Enter)

    assert pypeit_setup_model.spectrograph == "keck_mosfire"

    assert len(obs_log_tab.raw_data_paths.locations) == 0
    assert obs_log_tab.raw_data_paths.isEnabled()

    # set raw data
    raw_dir = os.path.join(os.getenv('PYPEIT_DEV'), 
                           'RAW_DATA', 'keck_mosfire')
    j_multi = os.path.join(raw_dir, "J_multi")
    with qtbot.waitSignal(pypeit_setup_model.raw_data_dirs_changed, raising=True, timeout=1000):
        obs_log_tab.raw_data_paths._location.setCurrentText(j_multi)
        qtbot.keyClick(obs_log_tab.raw_data_paths._location, Qt.Key_Enter)

    
    assert len(obs_log_tab.raw_data_paths.locations) == 1
    assert j_multi in obs_log_tab.raw_data_paths.locations
    assert obs_log_tab.raw_data_paths.isEnabled()
    assert obs_log_tab.spectrograph.isEnabled()
    assert main_window.openButton.isEnabled()
    assert not main_window.clearButton.isEnabled()
    assert main_window.setupButton.isEnabled()
    assert not main_window.saveTabButton.isEnabled()
    assert not main_window.saveAllButton.isEnabled()

    assert main_window.setup_view.state == model.ModelState.NEW
    assert main_window.setup_view.currentWidget().state == model.ModelState.UNCHANGED
    assert main_window.setup_view.currentWidget().tab_name == "ObsLog"

    # Click the Run Setup button
    with qtbot.waitSignal(pypeit_setup_model.operation_complete, raising=True, timeout=10000):
        qtbot.mouseClick(main_window.setupButton, Qt.MouseButton.LeftButton)

    assert main_window.saveAllButton.isEnabled()
    assert not main_window.saveTabButton.isEnabled()
    assert main_window.clearButton.isEnabled()

    assert main_window.setup_view.state == model.ModelState.CHANGED
    assert main_window.setup_view.currentWidget().state == model.ModelState.UNCHANGED
    assert main_window.setup_view.currentWidget().tab_name == "ObsLog"
    assert not obs_log_tab.spectrograph.isEnabled()

    tab_a = main_window.setup_view.widget(1) 
    tab_b = main_window.setup_view.widget(2) 
    assert tab_a.name == "A"
    assert tab_a.tab_name == "*A"
    assert tab_a.state == model.ModelState.CHANGED
    assert tab_a.file_metadata_table.model().rowCount() == 69
    assert tab_b.name == "B"
    assert tab_b.tab_name == "*B"
    assert tab_b.state == model.ModelState.CHANGED
    assert tab_b.file_metadata_table.model().rowCount() == 4

    # Verify history
    settings = QSettings()
    settings.beginGroup("RawDataDirectory")
    history = settings.value('History')
    assert history == [str(j_multi)]


def test_multi_paths(qtbot, tmp_path, monkeypatch):
    """Test re-running setup on setting multiple paths"""
    c, main_window = setup_offscreen_gui(tmp_path, monkeypatch, qtbot)    
    pypeit_setup_model = c.model
    obs_log_tab = main_window.setup_view._obs_log_tab

    # run_setup helper to run setup on one keck_mosfire setup
    run_setup("keck_mosfire", "Y_long", main_window, qtbot)

    # add a second raw data path
    raw_dir = os.path.join(os.getenv('PYPEIT_DEV'), 
                           'RAW_DATA', 'keck_mosfire')
    j_multi = os.path.join(raw_dir, "J_multi")
    y_long = os.path.join(raw_dir, "Y_long")
    with qtbot.waitSignal(pypeit_setup_model.raw_data_dirs_changed, raising=True, timeout=1000):
        obs_log_tab.raw_data_paths._location.setCurrentText(j_multi)
        qtbot.keyClick(obs_log_tab.raw_data_paths._location, Qt.Key_Enter)

    assert len(obs_log_tab.raw_data_paths.locations) == 2
    assert j_multi in obs_log_tab.raw_data_paths.locations
    assert y_long in obs_log_tab.raw_data_paths.locations
    assert obs_log_tab.raw_data_paths.isEnabled()
    assert not obs_log_tab.spectrograph.isEnabled()
    assert main_window.openButton.isEnabled()
    assert main_window.clearButton.isEnabled()
    assert main_window.setupButton.isEnabled()
    assert not main_window.saveTabButton.isEnabled()
    assert main_window.saveAllButton.isEnabled()

    assert main_window.setup_view.state == model.ModelState.CHANGED

    # It will ask us to discard if we run setup now, so we monkeypatch that dialog to say yes
    with monkeypatch.context() as m:
        m.setattr(main_window, "prompt_for_save", lambda: view.DialogResponses.ACCEPT)

        # Click the Run Setup button to re-run setup
        with qtbot.waitSignal(pypeit_setup_model.operation_complete, raising=True, timeout=10000):
            qtbot.mouseClick(main_window.setupButton, Qt.MouseButton.LeftButton)

    assert main_window.saveAllButton.isEnabled()
    assert not main_window.saveTabButton.isEnabled()
    assert main_window.clearButton.isEnabled()

    assert main_window.setup_view.state == model.ModelState.CHANGED
    assert main_window.setup_view.currentWidget().state == model.ModelState.UNCHANGED
    assert main_window.setup_view.currentWidget().tab_name == "ObsLog"
    assert not obs_log_tab.spectrograph.isEnabled()

    tab_a = main_window.setup_view.widget(1) 
    tab_b = main_window.setup_view.widget(2) 
    tab_c = main_window.setup_view.widget(3) 
    assert tab_a.name == "A"
    assert tab_a.tab_name == "*A"
    assert tab_a.state == model.ModelState.CHANGED
    assert tab_a.file_metadata_table.model().rowCount() == 69
    assert tab_b.name == "B"
    assert tab_b.tab_name == "*B"
    assert tab_b.state == model.ModelState.CHANGED
    assert tab_b.file_metadata_table.model().rowCount() == 4
    assert tab_c.name == "C"
    assert tab_c.tab_name == "*C"
    assert tab_c.state == model.ModelState.CHANGED
    assert tab_c.file_metadata_table.model().rowCount() == 16

    # Verify history
    settings = QSettings()
    settings.beginGroup("RawDataDirectory")
    history = settings.value('History')
    assert len(history) == 2
    assert str(j_multi) in history
    assert str(y_long) in history


def test_save_and_open(qtbot, tmp_path, monkeypatch):
    c, main_window = setup_offscreen_gui(tmp_path, monkeypatch, qtbot)    

    # Use the run_setup helper to run setup on J_muilti, which will create two tabs
    run_setup("keck_mosfire", "J_multi", main_window, qtbot)

    # Select Tab B, test that the save tab button is enabled
    with qtbot.waitSignal(main_window.setup_view.currentChanged, raising=True, timeout=1000):
        main_window.setup_view.setCurrentIndex(2)

    tab_b = main_window.setup_view.currentWidget()
    assert tab_b.tab_name == "*B"
    # The filename displayed should be just the base name until a location is set
    assert tab_b._model.filename == "keck_mosfire_B.pypeit"
    assert main_window.saveTabButton.isEnabled()


    # Set output directory. using a mock file dialog
    with monkeypatch.context() as m:
        m.setattr("pypeit.setup_gui.view.FileDialog", MockFileDialog)
        MockFileDialog.call_count = 0
        MockFileDialog.mock_response = view.DialogResponses.ACCEPT
        MockFileDialog.mock_location = str(tmp_path)
        # Save the tab
        with qtbot.waitSignal(tab_b._model.state_changed, raising=True, timeout=10000):
            qtbot.mouseClick(main_window.saveTabButton, Qt.MouseButton.LeftButton)

    # Assert the save file dialog was called once and only once
    assert MockFileDialog.call_count == 1

    # Assert that tab B's state changed, and the save tab button is disabled,
    # while tab A's state is unchanged and Save All is still enabled
    assert tab_b.state == model.ModelState.UNCHANGED
    assert tab_b.tab_name == "B"

    tab_a = main_window.setup_view.widget(1)
    assert tab_a.state == model.ModelState.CHANGED
    assert main_window.setup_view.state == model.ModelState.CHANGED
    assert not main_window.saveTabButton.isEnabled()
    assert main_window.saveAllButton.isEnabled()

    # Assert the file was created
    tab_b_file = tmp_path / "keck_mosfire_B/keck_mosfire_B.pypeit"
    assert tab_b_file.is_file()

    # Make sure tab_b properly report it's name
    assert tab_b._model.filename == str(tab_b_file)

    # Now test opening that file

    # It will ask us to discard the unchanged tab, and for the file to open, 
    # so we monkeypatch the dialogs

    with monkeypatch.context() as m:
        m.setattr(main_window, "prompt_for_save", lambda: view.DialogResponses.ACCEPT)
        m.setattr("pypeit.setup_gui.view.FileDialog", MockFileDialog)
        MockFileDialog.call_count = 0
        MockFileDialog.mock_response = view.DialogResponses.ACCEPT
        MockFileDialog.mock_location = str(tab_b_file)

        with qtbot.waitSignal(c.model.state_changed, raising=True, timeout=10000):
            qtbot.mouseClick(main_window.openButton, Qt.MouseButton.LeftButton)

    assert MockFileDialog.call_count == 1

    # Verify a single tab was opened and everthing is the correct state
    assert main_window.openButton.isEnabled()
    assert main_window.clearButton.isEnabled()
    assert main_window.setupButton.isEnabled()
    assert not main_window.saveTabButton.isEnabled()
    assert not main_window.saveAllButton.isEnabled()

    assert main_window.setup_view.state == model.ModelState.UNCHANGED
    assert main_window.setup_view.count() == 2
    assert main_window.setup_view.currentIndex() == 0
    assert main_window.setup_view.widget(0).state == model.ModelState.UNCHANGED
    assert main_window.setup_view.widget(0).tab_name == "ObsLog"
    assert main_window.setup_view.widget(1).state == model.ModelState.UNCHANGED
    assert main_window.setup_view.widget(1).tab_name == "B"
    assert main_window.setup_view.widget(1).file_metadata_table.model().rowCount() == 4

    # Verify history
    settings = QSettings()
    settings.beginGroup("OutputDirectory")
    history = settings.value('History')
    assert history == [str(tmp_path)]
    settings.endGroup()
    settings.beginGroup("OpenFile")
    history = settings.value('History')
    assert history == [str(tab_b_file.parent)]


def test_save_all(qtbot, tmp_path, monkeypatch):
    c, main_window = setup_offscreen_gui(tmp_path, monkeypatch, qtbot)    

    # Use the run_setup helper to run setup on J_muilti, which will create two tabs
    run_setup("keck_mosfire", "J_multi", main_window, qtbot)

    assert main_window.saveAllButton.isEnabled()

    # Select Tab B, test that the save tab button is enabled
    with qtbot.waitSignal(main_window.setup_view.currentChanged, raising=True, timeout=1000):
        main_window.setup_view.setCurrentIndex(2)

    tab_b = main_window.setup_view.currentWidget()
    assert tab_b.tab_name == "*B"
    assert main_window.saveTabButton.isEnabled()
    assert main_window.saveAllButton.isEnabled()

    # Click the save all button
    with monkeypatch.context() as m:
        m.setattr("pypeit.setup_gui.view.FileDialog", MockFileDialog)
        MockFileDialog.call_count = 0
        MockFileDialog.mock_response = view.DialogResponses.ACCEPT
        MockFileDialog.mock_location = str(tmp_path)

        with qtbot.waitSignal(tab_b._model.state_changed, raising=True, timeout=10000):
            qtbot.mouseClick(main_window.saveAllButton, Qt.MouseButton.LeftButton)

    # Assert the save file dialog was called twice, once for each tab
    assert MockFileDialog.call_count == 2

    # Assert both tabs are saved now, and that both the Save Tab and Save All buttons
    # are disabled
    tab_a = main_window.setup_view.widget(1)
    assert tab_a.state == model.ModelState.UNCHANGED
    assert tab_b.state == model.ModelState.UNCHANGED
    assert main_window.setup_view.state == model.ModelState.UNCHANGED
    assert not main_window.saveTabButton.isEnabled()
    assert not main_window.saveAllButton.isEnabled()

    tab_a_file = tmp_path / "keck_mosfire_A/keck_mosfire_A.pypeit"
    tab_b_file = tmp_path / "keck_mosfire_B/keck_mosfire_B.pypeit"

    assert tab_a_file.is_file()
    assert tab_b_file.is_file()

    # Retest with the ACCEPT_FOR_ALL response. First re-run setup. It should
    # *not* prompt us to save
    def raise_on_prompt():
        raise RuntimeError("The GUI should *not* be prompting to save.")

    with monkeypatch.context() as m:
        m.setattr(main_window, "prompt_for_save",  raise_on_prompt)

        # Click the Run Setup button to re-run setup
        with qtbot.waitSignal(c.model.operation_complete, raising=True, timeout=10000):
            qtbot.mouseClick(main_window.setupButton, Qt.MouseButton.LeftButton)

    tab_a_file.unlink()
    tab_b_file.unlink()

    # Re-test with ACCEPT_FOR_ALL

    tab_a = main_window.setup_view.widget(1)
    tab_b = main_window.setup_view.widget(2)

    # Click the save all button
    with monkeypatch.context() as m:
        m.setattr("pypeit.setup_gui.view.FileDialog", MockFileDialog)
        MockFileDialog.call_count = 0
        MockFileDialog.mock_response = view.DialogResponses.ACCEPT_FOR_ALL
        MockFileDialog.mock_location = str(tmp_path)
        with qtbot.waitSignal(tab_b._model.state_changed, raising=True, timeout=10000):
            qtbot.mouseClick(main_window.saveAllButton, Qt.MouseButton.LeftButton)

    # Assert the save file dialog was called once, with ACCEPT_FOR_ALL preventing the second tab
    # from prompting for a location
    assert MockFileDialog.call_count == 1

    # Verify the save worked
    assert tab_a.state == model.ModelState.UNCHANGED
    assert tab_b.state == model.ModelState.UNCHANGED
    assert main_window.setup_view.state == model.ModelState.UNCHANGED
    assert not main_window.saveTabButton.isEnabled()
    assert not main_window.saveAllButton.isEnabled()

    assert tab_a_file.is_file()
    assert tab_b_file.is_file()

    # Retest with the CANCEL response. 
    def raise_on_prompt():
        raise RuntimeError("The GUI should *not* be prompting to save.")

    with monkeypatch.context() as m:
        m.setattr(main_window, "prompt_for_save",  raise_on_prompt)

        # Click the Run Setup button to re-run setup
        with qtbot.waitSignal(c.model.operation_complete, raising=True, timeout=10000):
            qtbot.mouseClick(main_window.setupButton, Qt.MouseButton.LeftButton)

    tab_a_file.unlink()
    tab_b_file.unlink()

    tab_a = main_window.setup_view.widget(1)
    tab_b = main_window.setup_view.widget(2)

    # Re-test with CANCEL, with a location for tab 1. Testing that we only prompt for save
    # if there's no location, and that cancel works
    tab_a._model.save_location = str(tmp_path)
    # Click the save all button
    with monkeypatch.context() as m:
        m.setattr("pypeit.setup_gui.view.FileDialog", MockFileDialog)
        MockFileDialog.call_count = 0
        MockFileDialog.mock_response = view.DialogResponses.CANCEL
        MockFileDialog.mock_location = str(tmp_path)

        qtbot.mouseClick(main_window.saveAllButton, Qt.MouseButton.LeftButton)


    # Assert the save file dialog was called once, with CANCEL preventing the second tab
    # from saving
    assert MockFileDialog.call_count == 1

    # Verify tab_a was saved but tab_b was not.
    tab_b = main_window.setup_view.widget(2)
    assert tab_a.state == model.ModelState.UNCHANGED
    assert tab_b.state == model.ModelState.CHANGED
    assert main_window.setup_view.state == model.ModelState.CHANGED
    assert not main_window.saveTabButton.isEnabled()
    assert main_window.saveAllButton.isEnabled()

    assert tab_a_file.is_file()
    assert not tab_b_file.is_file()

def test_clear(qtbot, tmp_path, monkeypatch):
    """
    Test the "clear" button
    """

    c, main_window = setup_offscreen_gui(tmp_path, monkeypatch, qtbot)    

    # Use the run_setup helper to run setup on J_muilti, which will create two tabs
    run_setup("keck_mosfire", "J_multi", main_window, qtbot)

    pypeit_setup_model = c.model
    obs_log_tab = main_window.setup_view._obs_log_tab

    # Suppress the prompt for save dialog
    def mock_prompt_for_save():
        mock_prompt_for_save.call_count += 1
        return view.DialogResponses.ACCEPT
    mock_prompt_for_save.call_count = 0

    with monkeypatch.context() as m:
        m.setattr(main_window, "prompt_for_save", mock_prompt_for_save)
        
        # clear
        with qtbot.waitSignal(pypeit_setup_model.configs_deleted, raising=True, timeout=10000):
            qtbot.mouseClick(main_window.clearButton, Qt.MouseButton.LeftButton)

    # Make sure the prompt for save dialog was actually called
    assert mock_prompt_for_save.call_count == 1

    # Verify the GUI has been reset to its initial state
    assert obs_log_tab.spectrograph.currentIndex() == -1
    assert obs_log_tab.spectrograph.isEnabled()

    assert len(obs_log_tab.raw_data_paths.locations) == 0
    assert not obs_log_tab.raw_data_paths.isEnabled()

    assert main_window.openButton.isEnabled()
    assert not main_window.clearButton.isEnabled()
    assert not main_window.setupButton.isEnabled()
    assert not main_window.saveTabButton.isEnabled()
    assert not main_window.saveAllButton.isEnabled()

    assert main_window.setup_view.state == model.ModelState.NEW
    assert main_window.setup_view.widget(0).state == model.ModelState.UNCHANGED
    assert main_window.setup_view.widget(0).tab_name == "ObsLog"
    assert main_window.setup_view.count() == 1

def test_log_window(qtbot, tmp_path, monkeypatch):
    c, main_window = setup_offscreen_gui(tmp_path, monkeypatch, qtbot, verbosity=1, logfile=None)

    # Open the log window
    qtbot.mouseClick(main_window.logButton, Qt.MouseButton.LeftButton)

    # Use the run_setup helper to run setup on J_muilti, which will create two tabs
    run_setup("keck_mosfire", "J_multi", main_window, qtbot)

    logWindow = main_window._logWindow
    assert logWindow is not None

    # Log our own message for testing    
    info_msg ="Info message should appear in verbosity=1" 
    wip_msg = "WIP message should not appear in verbosity=1"

    with qtbot.waitSignal(logWindow.logViewer.textChanged, raising=True, timeout=1000):
        msgs.work(wip_msg)
        msgs.info(info_msg)

    log_text = logWindow.logViewer.toPlainText()
    assert info_msg in log_text
    assert wip_msg not in log_text

    # save the log
    log_file = tmp_path / "test_log_window.log"

    with monkeypatch.context() as m:
        m.setattr("pypeit.setup_gui.view.FileDialog", MockFileDialog)
        MockFileDialog.call_count = 0
        MockFileDialog.mock_response = view.DialogResponses.ACCEPT
        MockFileDialog.mock_location = str(log_file)

        # Wait for the "log saved" message
        with qtbot.waitSignal(logWindow.logViewer.textChanged, raising=True, timeout=10000):
            qtbot.mouseClick(main_window._logWindow.saveButton, Qt.MouseButton.LeftButton)

    # Assert the save file dialog was called once and only once
    assert MockFileDialog.call_count == 1

    assert log_file.is_file()
    # Look for our log message in the file
    found_info_msg = False
    found_wip_msg = False
    with open(log_file, "r") as f:
        for line in f:
            if info_msg in line:
                found_info_msg = True
            if wip_msg in line:
                found_wip_msg = True
    assert found_info_msg is True
    assert found_wip_msg is False

    # Close the log window
    with qtbot.waitSignal(logWindow.closed, raising=True, timeout=1000):
        qtbot.mouseClick(main_window._logWindow.closeButton, Qt.MouseButton.LeftButton)

    assert main_window._logWindow is None

