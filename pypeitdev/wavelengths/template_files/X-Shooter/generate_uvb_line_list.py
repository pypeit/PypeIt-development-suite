""" Quick script ot generate the linelist for X-Shooter
  Assumes a 1'' slit (approximate) """

from pypeit.wavemodel import create_ThArlinelist

# File generate is in vacuum
create_ThArlinelist(6000.)