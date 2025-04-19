import pytest

from collections import namedtuple
import shutil
import os
from pathlib import Path
import re
import gc
from qtpy.QtCore import Qt, QSettings, QDeadlineTimer

import numpy as np

from pypeit.setup_gui import dialog_helpers
from pypeit.setup_gui import view, model, controller
from pypeit.spectrographs.util import load_spectrograph
from pypeit.metadata import PypeItMetaData
from pypeit.inputfiles import PypeItFile
from pypeit import msgs

# Override pytest-qt's QApplication arguments
@pytest.fixture(scope="session")
def qapp_args():
    return ["pytest", "-platform", "offscreen"]

@pytest.fixture
def raw_data_path():
    return Path(os.environ['PYPEIT_DEV'], 'RAW_DATA')

def get_test_metadata(raw_data_path):
    # Create a PypeItMetadata to test with
    files = list((raw_data_path / "shane_kast_blue" / "600_4310_d55").absolute().glob("*.fits.gz"))
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
                "bpm_usebias = True",
                "[[biasframe]]",
                "frametype = bias",
                "[[[process]]]",
                "combine = median", ]

    return MockSetup(par, user_cfg, metadata, spec)


class MockFileDialog():
    mock_response = view.DialogResponses.CANCEL
    selected_path = None
    call_count = 0
    def __init__(parent, *args, **kwargs):
        pass

    def show(self):
        self.__class__.call_count += 1
        print(f"MockFileDialog Returning mock response: {self.mock_response}")
        return self.mock_response

    @classmethod
    def create_save_location_dialog(*args, **kwargs):
        return MockFileDialog()

    @classmethod
    def create_open_file_dialog(*args, **kwargs):
        return MockFileDialog()

def test_metadata_model(raw_data_path, qtbot, qtmodeltester):

    metadata = get_test_metadata(raw_data_path)
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
    

    # Test removing metadata rows
    rows_to_remove = np.where(np.isin(metadata_model.metadata['filename'], ['b22.fits.gz', 'b23.fits.gz', 'b24.fits.gz']))[0]
    orig_count = metadata_model.rowCount()
    with qtbot.waitSignals([metadata_model.rowsAboutToBeRemoved, metadata_model.rowsRemoved]*len(rows_to_remove), raising=True, order='strict', timeout=5000):
        metadata_model.removeMetadataRows(rows_to_remove)
    assert orig_count - len(rows_to_remove) == metadata_model.rowCount()
    assert not np.any(np.isin(metadata_model.metadata['filename'], ['b22.fits.gz', 'b23.fits.gz', 'b24.fits.gz']))

    # TODO: More tests?
    # Test paste
    # Test copyForRows
    # Test CopyFromConfig?

def test_commenting_out_files(raw_data_path, tmp_path, qtbot):

    metadata = get_test_metadata(raw_data_path)
    metadata.get_frame_types() # Needed so the metadata can be saved later
    metadata_model = model.PypeItMetadataModel(metadata=metadata)

    # Test commenting out 3 known files, and make sure the appropriate signal is sent
    files_to_comment = ['b22.fits.gz', 'b23.fits.gz', 'b24.fits.gz']
    rows_to_comment = [np.where(metadata['filename'] == file)[0][0] for file in files_to_comment]
    with qtbot.waitSignal(metadata_model.dataChanged, raising=True, timeout=1000):
        metadata_model.commentMetadataRows(rows_to_comment)

    filename_col = metadata_model.colnames.index('filename')
    for row, file in zip(rows_to_comment, files_to_comment):
        # Make sure file is commented out
        assert metadata_model.metadata[row]['filename'] == "# " + file

        # Make sure the comment character does not appear in the GUI
        indx =  metadata_model.createIndex(row, filename_col)
        assert metadata_model.data(indx, Qt.DisplayRole) == file

        # Make sure the model reports it as commented out
        assert metadata_model.isCommentedOut(indx) is True

    # Try commenting out the files again, to make sure double comments aren't added
    metadata_model.commentMetadataRows(rows_to_comment)
    for row in rows_to_comment:
        assert metadata_model.metadata[row]['filename'].count("#") == 1

    # Now save and re-open the file to test that these are saved correctly
    mock_setup = get_mock_setup(metadata)
    params_model = model.PypeItParamsProxy(mock_setup.par, mock_setup.user_cfg)
    config_dict = {'A': {'dispname': '600/4310', 'dichroic': 'd55'}}
    file_model = model.PypeItFileModel(metadata_model = metadata_model, name_stem="A", config_dict=config_dict, state=model.ModelState.CHANGED, params_model=params_model)
    file_model.save_location = str(tmp_path)
    file_model.save()
    assert file_model.state == model.ModelState.UNCHANGED
    pypeit_file_path = tmp_path / "shane_kast_blue_A.pypeit"
    assert pypeit_file_path.exists()

    gui_model = model.PypeItSetupGUIModel()
    with qtbot.waitSignals([gui_model.filesAdded, gui_model.stateChanged], order='strict', raising=True, timeout=5000):
        gui_model.open_pypeit_file(pypeit_file_path)

    opened_file_model = gui_model.pypeit_files['A']
    opened_metadata_model = opened_file_model.metadata_model

    # Make sure files are still commented out, without assuming the same order
    # while gathering the row numbers for later use
    commented_files = ["# " + file for file in files_to_comment]

    # We concatenate because the np.where will return an empty result if the file isn't there
    commented_rows = np.concatenate([np.where(opened_metadata_model.metadata['filename'] == file)[0] for file in commented_files])
    assert len(commented_rows) == len(commented_files)

    # Test uncommenting out the files
    with qtbot.waitSignal(opened_metadata_model.dataChanged, raising=True, timeout=1000):
        opened_metadata_model.uncommentMetadataRows(commented_rows)
    
    for row,file in zip(commented_rows,files_to_comment):
        assert opened_metadata_model.metadata[row]['filename'] == file

        # Make sure model reports it as uncommented
        indx =  metadata_model.createIndex(row, filename_col)
        assert opened_metadata_model.isCommentedOut(indx) is False


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

def test_pypeit_file_model(qtbot, raw_data_path, tmp_path):

    # Get a mock pypeit setup
    metadata = get_test_metadata(raw_data_path)
    # This will set "setup" to B for every file except b1.fits.gz
    metadata.get_frame_types()
    metadata.table['setup'] = 'A'
    metadata.table['setup'][metadata.table['decker'] == '2.0 arcsec'] = 'B'
    metadata_model = model.PypeItMetadataModel(metadata=metadata)
    mock_setup = get_mock_setup(metadata)
    params_model = model.PypeItParamsProxy(mock_setup.par, mock_setup.user_cfg)
    config_dict = {'B': {'dispname': '600/4310', 'dichroic': 'd55'}}
    file_model = model.PypeItFileModel(metadata_model = metadata_model, name_stem="B",
                                       config_dict=config_dict, state=model.ModelState.CHANGED,
                                       params_model=params_model)
    file_model.save_location = str(tmp_path)
    with qtbot.waitSignal(file_model.stateChanged, raising=True,
                          check_params_cb=verify_state_change, timeout=5000):
        file_model.save()
    
    # Make sure state changed
    assert file_model.state == model.ModelState.UNCHANGED

    # Make sure we can open the created file
    pf = PypeItFile.from_file(str(tmp_path / "shane_kast_blue_B.pypeit"))
    filenames = pf.filenames

    # Make sure one setup A and one setup B is there. This would probably be invalid
    # but the GUI trusts the user to know what they're doing
    shane_kast_blue_path = (raw_data_path / "shane_kast_blue" / "600_4310_d55").absolute()
    assert str((shane_kast_blue_path / "b11.fits.gz").absolute()) in filenames
    assert str((shane_kast_blue_path / "b1.fits.gz").absolute()) in filenames

    # TODO
    # test removed files not being saved
    # test comments in configuration section being saved

 
def test_pypeit_obslog_model(qtbot, raw_data_path, tmp_path):
    obslog_model = model.PypeItObsLogModel()
    assert obslog_model.state == model.ModelState.NEW

    # Copy only those fits files we want to run on
    shane_kast_blue_path = (raw_data_path / "shane_kast_blue" / "600_4310_d55").absolute()
    shutil.copy2(shane_kast_blue_path / "b1.fits.gz", tmp_path)
    shutil.copy2(shane_kast_blue_path / "b27.fits.gz", tmp_path)

    # Test setting the spectrograph, scanning the raw data directories, and then
    # setting metadata. THe same that would happen when running setup
    with qtbot.waitSignal(obslog_model.spectrograph_changed, raising=True, timeout=5000):
        obslog_model.set_spectrograph("shane_kast_blue")
    assert obslog_model.spec_name == "shane_kast_blue"
    assert obslog_model.spectrograph.name == "shane_kast_blue"
    with qtbot.waitSignals([(obslog_model.paths_model.rowsInserted, "add_row"),
                            (obslog_model.paths_model.dataChanged, "data_changed")],
                            raising=True, order='strict', timeout=5000): 
        obslog_model.add_raw_data_directory(str(tmp_path))

    raw_data_dirs = obslog_model.raw_data_directories
    assert len(raw_data_dirs) == 1
    assert raw_data_dirs[0] == str(tmp_path)
    
    raw_data_files = obslog_model.scan_raw_data_directories()
    assert len(raw_data_files) == 2
    assert str(tmp_path/"b1.fits.gz") in raw_data_files
    assert str(tmp_path/"b27.fits.gz") in raw_data_files

    # Reset and test with a path that has a pattern
    with qtbot.waitSignals([(obslog_model.paths_model.modelReset, "paths_reset"),
                            (obslog_model.spectrograph_changed, "spec_changed")],
                            raising=True, timeout=5000):
        obslog_model.reset()

    # Make sure the reset worked
    assert obslog_model.state == model.ModelState.NEW
    assert len(obslog_model.raw_data_directories) == 0
    assert obslog_model.spec_name is None
    assert obslog_model.spectrograph is None

    obslog_model.set_spectrograph("shane_kast_blue")
    obslog_model.add_raw_data_directory(str(tmp_path /"b1*"))

    raw_data_dirs = obslog_model.raw_data_directories
    assert len(raw_data_dirs) == 1
    assert raw_data_dirs[0] == str(tmp_path / "b1*")

    # The scan should now find only one file
    raw_data_files = obslog_model.scan_raw_data_directories()
    assert len(raw_data_files) == 1
    assert str(tmp_path/"b1.fits.gz") in raw_data_files

    # Test setting the metadata
    metadata = get_test_metadata(raw_data_path)

    with qtbot.waitSignals([(obslog_model.paths_model.modelReset, "paths_reset"),
                            (obslog_model.spectrograph_changed, "spec_changed")],
                            raising=True, timeout=5000):
        obslog_model.setMetadata(metadata)

    raw_data_dirs = obslog_model.raw_data_directories
    assert len(raw_data_dirs) == 1
    assert raw_data_dirs[0] == str(shane_kast_blue_path.absolute())
    
    assert obslog_model.spec_name == "shane_kast_blue"
    assert obslog_model.spectrograph.name == "shane_kast_blue"
    assert obslog_model.state == model.ModelState.UNCHANGED


def test_pypeit_setup_gui_model(qtbot, tmp_path):

    """
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
    """

def setup_offscreen_gui(tmp_path, monkeypatch, qapp, qtbot, verbosity=2):
    """Helper function to setup the gui to run in offscreen mode"""
    # Change to tmp_path so that log goes there
    monkeypatch.chdir(tmp_path)

    # Set config location and format so that history information is saved to a local ini file
    # instead of the user's home directory (linux and mac) or the registry (windows)
    QSettings.setDefaultFormat(QSettings.Format.IniFormat)
    QSettings.setPath(QSettings.Format.IniFormat, QSettings.Scope.UserScope, str(tmp_path / "config_user.ini"))
    QSettings.setPath(QSettings.Format.IniFormat, QSettings.Scope.SystemScope, str(tmp_path / "config_system.ini"))
    
    # Start the gui with with the "offscreen" platform for headless testing.
    c = controller.SetupGUIController(qapp, verbosity=verbosity)
    main_window = c.main_window

    # Patch display error so we don't get hung on it's modal dialog box
    def mock_display_error(parent, message):
        pass
    # We have to patch it on controller, even though it lives in dialog_helpers, because otherwise
    # controller will only see its own imported version and no patched one
    monkeypatch.setattr(controller, "display_error", mock_display_error)

    # qtbot will close the window for us, but we have to take care of the
    # class variables in SetupGUIController
    qtbot.addWidget(main_window)
    main_window.show()
    return c, main_window

def run_setup(spectrograph_name, setup_raw_data_path, main_controller, qtbot):
    """Run setup in an already started GUI. This is a helper method intended to set things up
    for tests that focus on functionality that require a populated GUI."""
    main_window = main_controller.main_window
    spectrograph_widget = main_window._obs_log_tab.spectrograph
    raw_data_paths_widget = main_window._obs_log_tab.paths_editor._path
    obslog_model = main_window.model.obslog_model
    with qtbot.waitSignal(obslog_model.spectrograph_changed, timeout=1000):
        qtbot.keyClicks(spectrograph_widget, spectrograph_name)
        qtbot.keyClick(spectrograph_widget, Qt.Key.Key_Enter)

    # set raw data
    with qtbot.waitSignals([(obslog_model.paths_model.rowsInserted, "path inserted"),
                            (obslog_model.paths_model.dataChanged, "path data set")], 
                            order = 'strict', raising=True, timeout=1000):
        raw_data_paths_widget.setCurrentText(str(setup_raw_data_path))
        qtbot.keyClick(raw_data_paths_widget, Qt.Key_Enter)

    # Click the Run Setup button
    with qtbot.waitSignal(main_window.model.stateChanged, raising=True, timeout=10000):
        main_window.setupAction.trigger()

    # Make sure the thread finishes before continuing
    assert main_controller.operation_thread.wait(QDeadlineTimer(1000)) is True

def test_run_setup(qapp, qtbot, raw_data_path, tmp_path, monkeypatch):
    """Test the basic process of setting a spectrograph, a raw data directory, and
    running setup."""
    c, main_window = setup_offscreen_gui(tmp_path, monkeypatch, qapp, qtbot)
    setup_gui_model = c.model

    # Validate initial state
    obs_log_tab = main_window._obs_log_tab
    assert obs_log_tab.spectrograph.currentIndex() == -1
    assert obs_log_tab.spectrograph.isEnabled()

    assert obs_log_tab._paths_viewer.model().rowCount() == 0
    assert obs_log_tab.paths_editor.history().rowCount() == 0
    assert not obs_log_tab.paths_editor.isEnabled()


    assert main_window.openAction.isEnabled()
    assert main_window.viewLogAction.isEnabled()
    assert not main_window.clearAction.isEnabled()
    assert not main_window.setupAction.isEnabled()
    assert not main_window.saveTabAction.isEnabled()
    assert not main_window.saveAllAction.isEnabled()

    assert obs_log_tab.state == model.ModelState.UNCHANGED
    assert main_window._current_tab is obs_log_tab
    assert main_window.tab_widget.currentWidget().name == "ObsLog"
    assert obs_log_tab.closeable is False
    assert setup_gui_model.state == model.ModelState.NEW

    # Select the "keck_mosfire" spectrograph 
    with qtbot.waitSignal(setup_gui_model.obslog_model.spectrograph_changed, timeout=1000):
        qtbot.keyClicks(obs_log_tab.spectrograph, "keck_mosfire")
        qtbot.keyClick(obs_log_tab.spectrograph, Qt.Key.Key_Enter)

    assert setup_gui_model.obslog_model.spectrograph.name == "keck_mosfire"

    # The paths editor should now be enabled
    assert obs_log_tab._paths_viewer.model().rowCount() == 0
    assert obs_log_tab.paths_editor.history().rowCount() == 0
    assert obs_log_tab.paths_editor.isEnabled()

    # set raw data
    raw_dir = raw_data_path / 'keck_mosfire'
    j_multi = (raw_dir / "J_multi").absolute()
    with qtbot.waitSignals([(setup_gui_model.obslog_model.paths_model.rowsInserted, "path inserted"),
                            (setup_gui_model.obslog_model.paths_model.dataChanged, "path data set")], 
                            order = 'strict', raising=True, timeout=1000):
        obs_log_tab.paths_editor._path.setCurrentText(str(j_multi))
        qtbot.keyClick(obs_log_tab.paths_editor._path, Qt.Key_Enter)

    
    assert obs_log_tab._paths_viewer.model().rowCount() == 1
    assert obs_log_tab._paths_viewer.model().stringList()[0] == str(j_multi)
    assert obs_log_tab.paths_editor.isEnabled()
    assert obs_log_tab.spectrograph.isEnabled()
    assert main_window.openAction.isEnabled()
    assert not main_window.clearAction.isEnabled()
    assert main_window.setupAction.isEnabled()
    assert not main_window.saveTabAction.isEnabled()
    assert not main_window.saveAllAction.isEnabled()

    assert setup_gui_model.state == model.ModelState.NEW
    assert obs_log_tab.state == model.ModelState.UNCHANGED
    assert main_window._current_tab is obs_log_tab
    assert main_window.tab_widget.currentWidget().name == "ObsLog"

    # Click the Run Setup button and wait for the setup operation to complete

    with qtbot.waitSignal(setup_gui_model.stateChanged, raising=True, timeout=10000):
        main_window.setupAction.trigger()
    assert c.operation_thread.wait(QDeadlineTimer(1000)) is True

    assert main_window.saveAllAction.isEnabled()
    assert not main_window.saveTabAction.isEnabled()
    assert main_window.clearAction.isEnabled()

    assert setup_gui_model.state == model.ModelState.CHANGED
    assert obs_log_tab.state == model.ModelState.UNCHANGED
    assert main_window._current_tab is obs_log_tab
    assert main_window.tab_widget.currentWidget().name == "ObsLog"
    assert not obs_log_tab.spectrograph.isEnabled()

    tab_a = main_window.tab_widget.widget(1) 
    tab_b = main_window.tab_widget.widget(2) 
    assert tab_a.name == "A"
    assert main_window.tab_widget.tabText(1) == "*A"
    assert tab_a.state == model.ModelState.CHANGED
    assert tab_a.closeable is True
    assert tab_a.file_metadata_table.model().rowCount() == 69
    assert tab_b.name == "B"
    assert main_window.tab_widget.tabText(2) == "*B"
    assert tab_b.state == model.ModelState.CHANGED
    assert tab_b.closeable is True
    assert tab_b.file_metadata_table.model().rowCount() == 4

    # Verify history
    assert obs_log_tab.paths_editor.history().rowCount() == 1
    assert obs_log_tab.paths_editor.history().stringList()[0] == str(j_multi)

    settings = QSettings()
    settings.beginGroup("RawDataDirectory")
    history = settings.value('History')
    assert history == [str(j_multi)]

def test_run_setup_failure(qapp, qtbot, raw_data_path, tmp_path, monkeypatch):
    c, main_window = setup_offscreen_gui(tmp_path, monkeypatch, qapp, qtbot)    
    setup_gui_model = c.model
    obs_log_tab = main_window._obs_log_tab

    # Select the spectrograph
    with qtbot.waitSignal(setup_gui_model.obslog_model.spectrograph_changed, timeout=1000):
        qtbot.keyClicks(obs_log_tab.spectrograph, "shane_kast_blue")
        qtbot.keyClick(obs_log_tab.spectrograph, Qt.Key.Key_Enter)

    # set raw data
    raw_dir = raw_data_path / 'keck_mosfire'
    j_multi = (raw_dir / "J_multi").absolute()
    with qtbot.waitSignals([(setup_gui_model.obslog_model.paths_model.rowsInserted, "path inserted"),
                            (setup_gui_model.obslog_model.paths_model.dataChanged, "path data set")], 
                            order = 'strict', raising=True, timeout=1000):
        obs_log_tab.paths_editor._path.setCurrentText(str(j_multi))
        qtbot.keyClick(obs_log_tab.paths_editor._path, Qt.Key_Enter)
    # Click the Run Setup button and wait for the background operation to complete, verifying it was canceled
    with qtbot.waitSignal(setup_gui_model.stateChanged, raising=True, timeout=10000):
        main_window.setupAction.trigger()
    assert c.operation_thread.wait(QDeadlineTimer(1000)) is True

    # Verify all of the files are commented out
    assert np.all([filename.startswith("#") for filename in setup_gui_model.obslog_model.metadata_model.metadata['filename']])

def test_run_setup_cancel(qapp, qtbot, raw_data_path, tmp_path, monkeypatch):
    c, main_window = setup_offscreen_gui(tmp_path, monkeypatch, qapp, qtbot)    
    setup_gui_model = c.model
    obs_log_tab = main_window._obs_log_tab

    # Select the wrong spectrograph for the test spectrograph 
    with qtbot.waitSignal(setup_gui_model.obslog_model.spectrograph_changed, timeout=1000):
        qtbot.keyClicks(obs_log_tab.spectrograph, "keck_mosfire")
        qtbot.keyClick(obs_log_tab.spectrograph, Qt.Key.Key_Enter)

    # set raw data
    raw_dir = raw_data_path / 'keck_mosfire'
    j_multi = (raw_dir / "J_multi").absolute()
    with qtbot.waitSignals([(setup_gui_model.obslog_model.paths_model.rowsInserted, "path inserted"),
                            (setup_gui_model.obslog_model.paths_model.dataChanged, "path data set")], 
                            order = 'strict', raising=True, timeout=1000):
        obs_log_tab.paths_editor._path.setCurrentText(str(j_multi))
        qtbot.keyClick(obs_log_tab.paths_editor._path, Qt.Key_Enter)


    # Helpers to verify and cause a cancel
    def verify_canceled_setup(canceled, exc_info):
        return canceled is True and (exc_info is None or all([value is None for value in exc_info]))

    def raise_cancel_exception(*args):
        raise controller.OpCanceledError()

    # Click the Run Setup button and trigger a cancel exception

    with qtbot.waitSignal(c.operation_thread.completed, check_params_cb=verify_canceled_setup, raising=True, timeout=10000):
        # Use the GUI's internal log watcher to trigger the canceled exception
        setup_gui_model.log_buffer.watch("test_cancel", re.compile("Adding metadata for .*$"), raise_cancel_exception)
        main_window.setupAction.trigger()

    # Make sure the thread finishes before continuing
    assert c.operation_thread.wait(QDeadlineTimer(1000)) is True

    assert setup_gui_model.state == model.ModelState.NEW
 
    setup_gui_model.log_buffer.unwatch("test_cancel")

def test_multi_paths(qapp, qtbot, raw_data_path, tmp_path, monkeypatch):
    """Test re-running setup on setting multiple paths"""
    c, main_window = setup_offscreen_gui(tmp_path, monkeypatch, qapp, qtbot)    
    setup_gui_model = c.model
    obs_log_tab = main_window._obs_log_tab

    raw_dir = raw_data_path / 'keck_mosfire'

    # run_setup helper to run setup on one keck_mosfire setup
    y_long = str((raw_dir / "Y_long").absolute())
    run_setup("keck_mosfire", y_long, c, qtbot)

    # add a second raw data path
    j_multi = str((raw_dir / "J_multi").absolute())

    with qtbot.waitSignals([(setup_gui_model.obslog_model.paths_model.rowsInserted, "path inserted"),
                            (setup_gui_model.obslog_model.paths_model.dataChanged, "path data set")], 
                            order = 'strict', raising=True, timeout=1000):
        obs_log_tab.paths_editor._path.setCurrentText(j_multi)
        qtbot.keyClick(obs_log_tab.paths_editor._path, Qt.Key_Enter)


    assert obs_log_tab._paths_viewer.model().rowCount() == 2
    assert j_multi in obs_log_tab._paths_viewer.model().stringList()
    assert y_long in obs_log_tab._paths_viewer.model().stringList()
    assert obs_log_tab.paths_editor.isEnabled()
    assert not obs_log_tab.spectrograph.isEnabled()
    assert main_window.openAction.isEnabled()
    assert main_window.clearAction.isEnabled()
    assert main_window.setupAction.isEnabled()
    assert not main_window.saveTabAction.isEnabled()
    assert main_window.saveAllAction.isEnabled()

    assert setup_gui_model.state == model.ModelState.CHANGED

    # It will ask us to discard if we run setup now, so we monkeypatch that dialog to say yes
    with monkeypatch.context() as m:
        m.setattr(controller, "prompt_to_save", lambda x: dialog_helpers.DialogResponses.ACCEPT)

        # Click the Run Setup button to re-run setup
        with qtbot.waitSignal(setup_gui_model.stateChanged, raising=True, timeout=10000):
            main_window.setupAction.trigger()

        # Make sure the thread finishes before continuing
        assert c.operation_thread.wait(QDeadlineTimer(1000)) is True

    assert main_window.saveAllAction.isEnabled()
    assert not main_window.saveTabAction.isEnabled()
    assert main_window.clearAction.isEnabled()

    assert setup_gui_model.state == model.ModelState.CHANGED
    assert main_window.tab_widget.currentWidget().state == model.ModelState.UNCHANGED
    assert main_window.tab_widget.currentWidget() is obs_log_tab
    assert not obs_log_tab.spectrograph.isEnabled()

    tab_widget = main_window.tab_widget
    tab_a = tab_widget.widget(1) 
    tab_b = tab_widget.widget(2) 
    tab_c = tab_widget.widget(3) 
    assert tab_a.name == "A"
    assert tab_widget.tabText(1) == "*A"
    assert tab_a.state == model.ModelState.CHANGED
    assert tab_a.file_metadata_table.model().rowCount() == 69
    assert tab_b.name == "B"
    assert tab_widget.tabText(2) == "*B"
    assert tab_b.state == model.ModelState.CHANGED
    assert tab_b.file_metadata_table.model().rowCount() == 4
    assert tab_c.name == "C"
    assert tab_widget.tabText(3) == "*C"
    assert tab_c.state == model.ModelState.CHANGED
    assert tab_c.file_metadata_table.model().rowCount() == 16

    # Verify history
    settings = QSettings()
    settings.beginGroup("RawDataDirectory")
    history = settings.value('History')
    assert len(history) == 2
    assert str(j_multi) in history
    assert str(y_long) in history


def test_save_and_open(qapp, qtbot, raw_data_path, tmp_path, monkeypatch):
    c, main_window = setup_offscreen_gui(tmp_path, monkeypatch, qapp, qtbot)    

    # Use the run_setup helper to run setup on J_muilti, which will create two tabs
    run_setup("keck_mosfire", (raw_data_path / "keck_mosfire" / "J_multi").absolute(), c, qtbot)

    tab_widget = main_window.tab_widget

    # Select Tab B, test that the save tab button is enabled
    with qtbot.waitSignal(tab_widget.currentChanged, raising=True, timeout=1000):
        tab_widget.setCurrentIndex(2)

    tab_b = tab_widget.currentWidget()
    assert tab_widget.tabText(2) == "*B"
    # The filename displayed should be just the base name until a location is set
    assert tab_b.model.filename == "keck_mosfire_B.pypeit"
    assert main_window.saveTabAction.isEnabled()


    # Set output directory. using a mock file dialog
    with monkeypatch.context() as m:
        m.setattr(controller, "FileDialog", MockFileDialog)
        MockFileDialog.call_count = 0
        MockFileDialog.mock_response = dialog_helpers.DialogResponses.ACCEPT
        MockFileDialog.selected_path = str(tmp_path)
        # Save the tab
        with qtbot.waitSignal(tab_b.model.stateChanged, raising=True, timeout=10000):
            main_window.saveTabAction.trigger()

    # Assert the save file dialog was called once and only once
    assert MockFileDialog.call_count == 1

    # Assert that tab B's state changed, and the save tab button is disabled,
    # while tab A's state is unchanged and Save All is still enabled
    assert tab_b.state == model.ModelState.UNCHANGED
    assert tab_widget.tabText(2) == "B"

    tab_a = tab_widget.widget(1)
    assert tab_a.state == model.ModelState.CHANGED
    assert tab_widget.tabText(1) == "*A"
    assert main_window.model.state == model.ModelState.CHANGED
    assert not main_window.saveTabAction.isEnabled()
    assert main_window.saveAllAction.isEnabled()

    # Assert the file was created
    tab_b_file = tmp_path / "keck_mosfire_B.pypeit"
    assert tab_b_file.is_file()

    # Make sure tab_b properly report it's name
    assert tab_b.model.filename == str(tab_b_file)

    # Now test opening that file

    # It will ask us to discard the unchanged tab, and for the file to open, 
    # so we monkeypatch the dialogs

    with monkeypatch.context() as m:
        m.setattr(controller, "prompt_to_save", lambda x: view.DialogResponses.ACCEPT)
        m.setattr(controller, "FileDialog", MockFileDialog)
        MockFileDialog.call_count = 0
        MockFileDialog.mock_response = view.DialogResponses.ACCEPT
        MockFileDialog.selected_path = str(tab_b_file)

        with qtbot.waitSignal(c.model.stateChanged, raising=True, timeout=10000):
            main_window.openAction.trigger()

        # Make sure the thread finishes before continuing
        assert c.operation_thread.wait(QDeadlineTimer(1000)) is True

    assert MockFileDialog.call_count == 1

    # Verify a single tab was opened and everthing is the correct state
    assert main_window.openAction.isEnabled()
    assert main_window.clearAction.isEnabled()
    assert main_window.setupAction.isEnabled()
    assert not main_window.saveTabAction.isEnabled()
    assert not main_window.saveAllAction.isEnabled()

    assert main_window.model.state == model.ModelState.UNCHANGED
    assert tab_widget.count() == 3 # The tab, the obslog tab, and the "new tab" tab
    assert tab_widget.currentIndex() == 0
    assert tab_widget.widget(0).state == model.ModelState.UNCHANGED
    assert tab_widget.widget(0).name == "ObsLog"
    assert tab_widget.widget(1).state == model.ModelState.UNCHANGED
    assert tab_widget.widget(1).name == "B"
    assert tab_widget.tabText(1) == "B"
    assert tab_widget.widget(1).file_metadata_table.model().rowCount() == 4

# TODO test dialog_helpers directly and verify it's history abilities

def test_save_all(qapp, qtbot, raw_data_path, tmp_path, monkeypatch):
    c, main_window = setup_offscreen_gui(tmp_path, monkeypatch, qapp, qtbot)    

    # Use the run_setup helper to run setup on J_muilti, which will create two tabs
    run_setup("keck_mosfire", (raw_data_path / "keck_mosfire" / "J_multi").absolute(), c, qtbot)

    assert main_window.saveAllAction.isEnabled()

    # Select Tab B, test that the save tab button is enabled
    with qtbot.waitSignal(main_window.tab_widget.currentChanged, raising=True, timeout=1000):
        main_window.tab_widget.setCurrentIndex(2)

    tab_b = main_window.tab_widget.currentWidget()
    assert tab_b.name == "B"
    assert main_window.tab_widget.tabText(2) == "*B"
    assert main_window.saveTabAction.isEnabled()
    assert main_window.saveAllAction.isEnabled()

    # Click the save all button
    with monkeypatch.context() as m:
        m.setattr(controller, "FileDialog", MockFileDialog)
        MockFileDialog.call_count = 0
        MockFileDialog.mock_response = dialog_helpers.DialogResponses.ACCEPT
        MockFileDialog.selected_path = str(tmp_path)

        with qtbot.waitSignal(tab_b.model.stateChanged, raising=True, timeout=10000):
            main_window.saveAllAction.trigger()

    # Assert the save file dialog was called twice, once for each tab
    assert MockFileDialog.call_count == 2

    # Assert both tabs are saved now, and that both the Save Tab and Save All buttons
    # are disabled
    tab_a = main_window.tab_widget.widget(1)
    assert tab_a.state == model.ModelState.UNCHANGED
    assert tab_b.state == model.ModelState.UNCHANGED
    assert c.model.state == model.ModelState.UNCHANGED
    assert not main_window.saveTabAction.isEnabled()
    assert not main_window.saveAllAction.isEnabled()

    tab_a_file = tmp_path / "keck_mosfire_A.pypeit"
    tab_b_file = tmp_path / "keck_mosfire_B.pypeit"

    assert tab_a_file.is_file()
    assert tab_b_file.is_file()

    # Retest with the ACCEPT_FOR_ALL response. First re-run setup. It should
    # *not* prompt us to save
    def raise_on_prompt():
        raise RuntimeError("The GUI should *not* be prompting to save.")

    with monkeypatch.context() as m:
        m.setattr(controller, "prompt_to_save",  raise_on_prompt)

        # Click the Run Setup button to re-run setup
        with qtbot.waitSignal(c.model.stateChanged, raising=True, timeout=10000):
            main_window.setupAction.trigger()

    tab_a_file.unlink()
    tab_b_file.unlink()

    # Re-test with ACCEPT_FOR_ALL

    tab_a = main_window.tab_widget.widget(1)
    tab_b = main_window.tab_widget.widget(2)

    # Click the save all button
    with monkeypatch.context() as m:
        m.setattr(controller, "FileDialog", MockFileDialog)
        MockFileDialog.call_count = 0
        MockFileDialog.mock_response = view.DialogResponses.ACCEPT_FOR_ALL
        MockFileDialog.selected_path = str(tmp_path)
        with qtbot.waitSignal(tab_b.model.stateChanged, raising=True, timeout=10000):
            main_window.saveAllAction.trigger()

    # Assert the save file dialog was called once, with ACCEPT_FOR_ALL preventing the second tab
    # from prompting for a location
    assert MockFileDialog.call_count == 1

    # Verify the save worked
    assert tab_a.state == model.ModelState.UNCHANGED
    assert tab_b.state == model.ModelState.UNCHANGED
    assert c.model.state == model.ModelState.UNCHANGED
    assert not main_window.saveTabAction.isEnabled()
    assert not main_window.saveAllAction.isEnabled()

    assert tab_a_file.is_file()
    assert tab_b_file.is_file()

    # Retest with the CANCEL response. 
    def raise_on_prompt():
        raise RuntimeError("The GUI should *not* be prompting to save.")

    with monkeypatch.context() as m:
        m.setattr(controller, "prompt_to_save",  raise_on_prompt)

        # Click the Run Setup button to re-run setup
        with qtbot.waitSignal(c.model.stateChanged, raising=True, timeout=10000):
            main_window.setupAction.trigger()

        # Make sure the thread finishes before continuing
        assert c.operation_thread.wait(QDeadlineTimer(1000)) is True

    tab_a_file.unlink()
    tab_b_file.unlink()

    tab_a = main_window.tab_widget.widget(1)
    tab_b = main_window.tab_widget.widget(2)

    # Re-test with CANCEL, with a location for tab 1. Testing that we only prompt for save
    # if there's no location, and that cancel works
    tab_a.model.save_location = str(tmp_path)

    # Click the save all button
    with monkeypatch.context() as m:
        m.setattr(controller, "FileDialog", MockFileDialog)
        MockFileDialog.call_count = 0
        MockFileDialog.mock_response = view.DialogResponses.CANCEL
        MockFileDialog.selected_path = str(tmp_path)

        main_window.saveAllAction.trigger()
        qtbot.waitUntil(lambda: MockFileDialog.call_count > 0 ,timeout=10000)

    assert MockFileDialog.call_count == 1

    # Verify tab_a was saved but tab_b was not.
    tab_b = main_window.tab_widget.widget(2)
    assert tab_a.state == model.ModelState.UNCHANGED
    assert tab_b.state == model.ModelState.CHANGED
    assert c.model.state == model.ModelState.CHANGED
    assert not main_window.saveTabAction.isEnabled()
    assert main_window.saveAllAction.isEnabled()

    assert tab_a_file.is_file()
    assert not tab_b_file.is_file()

def test_clear(qapp, qtbot, raw_data_path, tmp_path, monkeypatch):
    """
    Test the "clear" button
    """

    c, main_window = setup_offscreen_gui(tmp_path, monkeypatch, qapp, qtbot)    

    # Use the run_setup helper to run setup on J_muilti, which will create two tabs
    run_setup("keck_mosfire", (raw_data_path / "keck_mosfire" / "J_multi").absolute(), c, qtbot)

    setup_gui_model = c.model
    obs_log_tab = main_window._obs_log_tab

    # Suppress the prompt for save dialog
    def mock_prompt_to_save(*args, **kwargs):
        mock_prompt_to_save.call_count += 1
        return dialog_helpers.DialogResponses.ACCEPT
    mock_prompt_to_save.call_count = 0

    with monkeypatch.context() as m:
        m.setattr(controller, "prompt_to_save", mock_prompt_to_save)
        
        # clear
        with qtbot.waitSignal(setup_gui_model.filesDeleted, raising=True, timeout=10000):
            main_window.clearAction.trigger()

    # Make sure the prompt for save dialog was actually called
    assert mock_prompt_to_save.call_count == 1

    # Verify the GUI has been reset to its initial state
    assert obs_log_tab.spectrograph.currentIndex() == -1
    assert obs_log_tab.spectrograph.isEnabled()

    assert obs_log_tab._paths_viewer.model().rowCount() == 0
    assert not obs_log_tab.paths_editor.isEnabled()

    assert main_window.openAction.isEnabled()
    assert not main_window.clearAction.isEnabled()
    assert not main_window.setupAction.isEnabled()
    assert not main_window.saveTabAction.isEnabled()
    assert not main_window.saveAllAction.isEnabled()

    assert c.model.state == model.ModelState.NEW
    assert main_window.tab_widget.widget(0).state == model.ModelState.UNCHANGED
    assert main_window.tab_widget.widget(0).name == "ObsLog"
    assert main_window.tab_widget.count() == 2

def test_log_window(qapp, qtbot, tmp_path, monkeypatch):
    c, main_window = setup_offscreen_gui(tmp_path, monkeypatch, qapp, qtbot, verbosity=1)

    # Open the log window
    main_window.viewLogAction.trigger()
    qtbot.waitUntil(lambda: main_window._logWindow is not None and main_window._logWindow.isVisible(),timeout=5000)

    logWindow = main_window._logWindow

    # Log our own message for testing    
    info_msg ="Info message should appear in verbosity=1" 
    wip_msg = "WIP message should not appear in verbosity=1"

    with qtbot.waitSignal(logWindow.textViewer.textChanged, raising=True, timeout=1000):
        msgs.work(wip_msg)
        msgs.info(info_msg)

    log_text = logWindow.textViewer.toPlainText()
    assert info_msg in log_text
    assert wip_msg not in log_text

    # save the log
    log_file = tmp_path / "test_log_window.log"

    with monkeypatch.context() as m:
        m.setattr("pypeit.setup_gui.text_viewer.FileDialog", MockFileDialog)
        MockFileDialog.call_count = 0
        MockFileDialog.mock_response = view.DialogResponses.ACCEPT
        MockFileDialog.selected_path = str(log_file)

        # Wait for the "log saved" message
        with qtbot.waitSignal(logWindow.textViewer.textChanged, raising=True, timeout=10000):
            qtbot.mouseClick(logWindow.saveButton, Qt.MouseButton.LeftButton)

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
        qtbot.mouseClick(logWindow.closeButton, Qt.MouseButton.LeftButton)

    assert main_window._logWindow is None

