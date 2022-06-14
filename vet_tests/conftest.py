# Local pytest plugin to get the REDUX_OUT location from the pytest command line
# using a "redux_out" fixture

import pytest
import os 

def pytest_addoption(parser):
    parser.addoption("--redux_out", action="store", 
                     default=os.path.join(os.getenv('PYPEIT_DEV'), "REDUX_OUT"),
                     help="Location of dev-suite REDUX_OUT directory")

@pytest.fixture
def redux_out(request):
    return request.config.getoption("--redux_out")