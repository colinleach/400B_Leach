import numpy as np
from numpy.linalg import norm
import pandas as pd
# from galaxy.db import DB
from galaxy.galaxies import Galaxies
from galaxy.timecourse import TimeCourse

tc = TimeCourse(usesql=True)

tc.write_total_angmom()