# Source this file to enable running the dev suite in a headless linux environment
# It requires installing Xvfb
# e.g. 
#
# $ source source_headless_test.sh
# $ ./pypeit_test all
#
Xvfb :1 -screen 0 800x600x24 &
export DISPLAY=1:0
