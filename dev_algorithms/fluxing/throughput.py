
import numpy as np
from astropy import units as u
from astropy import constants as const
from pypeit import sensfunc



def sensfunc_to_thru(sensfile):


    wave, sensfunction, meta_table, out_table, header_sens = sensfunc.SensFunc.load(sensfile)
