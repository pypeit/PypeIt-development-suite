import pickle
import numpy as np
from pypeit.core.gui import object_find as gui_object_find


def load_pkl(name):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)


frame = np.load("frame.npy")
trace_dict = load_pkl("trace_dict")
sobjs = load_pkl("sobjs")
gui_object_find.initialise(frame, trace_dict, sobjs)
