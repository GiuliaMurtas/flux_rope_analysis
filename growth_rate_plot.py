#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pipreadmods
import h5py

filename = "growth_rate_PIP_1.h5"
with h5py.File(filename, "r") as f:
    print("Keys: %s" % f.keys())