import pickle
import numpy as np
from pypeit.core.gui import object_find as gui_object_find


def load_pkl(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


det = 1
frame = np.load("frame.npy")
trace_dict = load_pkl("trace_dict")
sobjs = load_pkl("sobjs")
# Make some updates
trace_dict["slit_left"] = trace_dict["edges_l"]
trace_dict["slit_righ"] = trace_dict["edges_r"]
gui_object_find.initialise(det, frame, trace_dict, sobjs)
