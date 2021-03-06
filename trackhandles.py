
# -*- coding: utf-8 -*-
"""
Track Handles
--------------

This separate pieces of code, though repeated at first glance, each determine a different Python-Object, which
would require very formal Python to condense. For clarity, we have chosen the *long* version of explictly stating all objects. 
"""
import matplotlib.patches as mpatches
class MObject(object):
    pass
class Hurricane1Object(object):
    pass
class Hurricane3Object(object):
    pass
class MajorHurricaneObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        patch = mpatches.Circle([13, 4],6.5,facecolor='red',
                                   transform=handlebox.get_transform())
        handlebox.add_artist(patch)
        return patch
class Hurricane1ObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        patch = mpatches.Circle([13, 4],4.5,facecolor='yellow',
                                   transform=handlebox.get_transform())
        handlebox.add_artist(patch)
        return patch
class Hurricane3ObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        patch = mpatches.Circle([13, 4],5.25,facecolor='orange',
                                   transform=handlebox.get_transform())
        handlebox.add_artist(patch)
        return patch
class TSObject(object):
    pass
class TDObject(object):
    pass
class EXObject(object):
    pass
class TWObject(object):
    pass
class TSObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        patch = mpatches.Circle([13, 4],4,facecolor='aqua',
                                   transform=handlebox.get_transform())
        handlebox.add_artist(patch)
        return patch
class TDObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        center = 0.5 * width - 0.5 * handlebox.xdescent, 0.5 * height - 0.5 * handlebox.ydescent
        patch = mpatches.Circle([13,4],3,facecolor='blue',
                                   transform=handlebox.get_transform())
        handlebox.add_artist(patch)
        return patch
class EXObjectHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        center = 0.5 * width - 0.5 * handlebox.xdescent, 0.5 * height - 0.5 * handlebox.ydescent
        patch = mpatches.Circle([13,4],4,facecolor='magenta',
                                   transform=handlebox.get_transform())
        handlebox.add_artist(patch)
        return patch
