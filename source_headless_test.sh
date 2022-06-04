# Source this file to enable running the dev suite in a headless linux environment
# It requires the QT libraries, which on Ubuntu 20.04 can be done with apt install -y qt5-default
# e.g. 
#
# $ source source_headless_test.sh
# $ ./pypeit_test all
#
export QT_QPA_PLATFORM=offscreen