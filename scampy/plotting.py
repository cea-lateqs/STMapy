import numpy as np
from .processing import levelTopo


def waterfallPlot(ax, x_indices, yz_data, offset):
    zf = len(x_indices)
    for z, z_data in enumerate(yz_data):
        z_offset = z * offset
        ax.plot(x_indices, z_data + z_offset, "k", zorder=(z + 1) * 2)
        # White filling under the curves
#        ax.fill_between(
#            x_indices,
#            z_data + z_offset,
#            offset,
#            facecolor="w",
#            lw=0,
#            zorder=(z + 1) * 2 - 1,
#        )
    ax.set_xlim([x_indices[0], x_indices[-1]])


def lineTopoPlot(fig, max_x, topo):
    """ Draws the topography read while opening a Line Spectro """
    yPx, xPx = topo
    # yPx is 1 (line spectro)
    assert yPx == 1

    ax = fig.add_subplot(111)
    x_data = np.linspace(0, max_x, xPx)

    ax.plot(x_data, topo[0], label="Without line leveling")
    ax.plot(x_data, levelTopo(topo)[0], label="With line leveling")
    ax.set_xlim([x_data[0], x_data[-1]])
    ax.legend(loc=0)
    return ax


def imageTopoPlot(fig, max_x, max_y, topo, colormap, line_fit):
    """ Draws the topography read while opening a 2D CITS."""
    # Get parameters
    yPx, xPx = topo.shape

    ax = fig.add_subplot(111)

    # pcolormesh takes *vertices* in arguments
    # so the X (Y) array need to be from 0 to W (H) INCLUDED
    XYmap = ax.pcolormesh(
        np.linspace(0, max_x, xPx + 1),
        np.linspace(0, max_y, yPx + 1),
        levelTopo(topo) if line_fit else topo,
        cmap=colormap,
    )

    # Colorbar stuff
    cbar = fig.colorbar(XYmap, shrink=0.9, pad=0.05, aspect=15)
    cbar.ax.yaxis.set_ticks_position("both")
    cbar.ax.tick_params(axis="y", direction="in")
    return ax
