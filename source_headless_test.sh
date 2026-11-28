# Source this file to enable running the dev suite in a headless linux environment
# It requires a Qt backend; by default PypeIt installs pyqt6, whose wheels bundle
# the Qt libraries.  See https://pypeit.readthedocs.io/en/latest/installing.html
# e.g.
#
# $ source source_headless_test.sh
# $ ./pypeit_test all
#
export QT_QPA_PLATFORM=offscreen