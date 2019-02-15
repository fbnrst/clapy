import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.transforms import Bbox, TransformedBbox, blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import BboxConnector
import matplotlib.lines as pltlines
import matplotlib.transforms as plttransformsimport 

import pandas as pd
idx = pd.IndexSlice


import warnings


def idx_name_levels(data,name):
    pos = data.index.names.index(name)
    return data.index.levels[pos]

def idx_name_labels(data,name):
    pos = data.index.names.index(name)
    return data.index.labels[pos]

def idx_name(df,name):
    return idx_name_levels(df,name)[idx_name_labels(df,name)]


def plot2d_pandas(df,column,levelx,levely,return_data=False,sym_limits=None,sym_scale=1,axis=None,**kwargs):
    xindx = idx_name_levels(df,levelx)
    yindx = idx_name_levels(df,levely)
    datag = df[column].groupby(level=[levely,levelx])
    if np.sum(datag.count() == 1) != len(yindx)*len(xindx):
        warnings.warn("DataFrame has more than 2 dimension. Average is plotted")
    data = datag.mean().sort_index().values.reshape((len(yindx),len(xindx)))
    dx = (xindx[1] - xindx[0])
    dy = (yindx[1] - yindx[0])
    ext = [xindx[0] - dx/2,xindx[-1] + dx/2,yindx[0] - dy/2,yindx[-1] + dy/2]
    args = dict()
    if sym_limits is not None:
        limit = np.absolute((data - sym_limits)).max()*sym_scale
        args.update({'vmin':sym_limits-limit,'vmax':sym_limits+limit})
    if axis is not None:
        im = axis.imshow(data,origin='lower',extent=ext,**args,**kwargs)
    else:
        im = plt.imshow(data,origin='lower',extent=ext,**args,**kwargs)

    if return_data:
        return im,data
    return im

def zoom_effect02(ax1, ax2,loc1a=3, loc2a=2, loc1b=4, loc2b=1, **kwargs):
    """
    ax1 : the main axes
    ax1 : the zoomed axes

    Similar to zoom_effect01.  The xmin & xmax will be taken from the
    ax1.viewLim.
    """

    tt = ax1.transScale + (ax1.transLimits + ax2.transAxes)
    trans = blended_transform_factory(ax2.transData, tt)

    bbox1 = ax1.bbox
    bbox2 = TransformedBbox(ax1.viewLim, trans)
    c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **kwargs)
    c1.set_clip_on(False)
    c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **kwargs)
    c2.set_clip_on(False)

    ax2.add_patch(c1)
    ax2.add_patch(c2)

    return c1, c2
