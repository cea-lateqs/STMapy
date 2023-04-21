import numpy as np
from .processing import levelTopoLine, levelTopoPlane


def waterfallPlot(ax, x_indices, yz_data, offset):
    for z, z_data in enumerate(yz_data):
        z_offset = z * offset
        ax.plot(x_indices, z_data + z_offset, "k", zorder=(z + 1) * 2)
    ax.set_xlim([x_indices[0], x_indices[-1]])


def lineTopoPlot(fig, max_x, topo):
    """ Draws the topography read while opening a Line Spectro """
    yPx, xPx = topo
    # yPx is 1 (line spectro)
    assert yPx == 1

    ax = fig.add_subplot(111)
    x_data = np.linspace(0, max_x, xPx)

    ax.plot(x_data, topo[0], label="Without line leveling")
    ax.plot(x_data, levelTopoLine(topo)[0], label="With line leveling")
    ax.set_xlim([x_data[0], x_data[-1]])
    ax.legend(loc=0)
    return ax


def imageTopoPlot(fig, max_x, max_y, topo, colormap, level_algo):
    """ Draws the topography read while opening a 2D CITS."""
    if level_algo == "line":
        topoToPlot = levelTopoLine(topo)
        fig.suptitle("Leveled topo (line fit)")
    elif level_algo == "plane":
        topoToPlot = levelTopoPlane(topo)
        fig.suptitle("Leveled topo (plane fit)")
    else:
        topoToPlot = topo
        fig.suptitle("Raw topo data")

    # Get parameters
    yPx, xPx = topo.shape

    ax = fig.add_subplot(1,1,1,aspect=float(max_y)/max_x)

    # pcolormesh takes *vertices* in arguments
    # so the X (Y) array need to be from 0 to W (H) INCLUDED
    XYmap = ax.pcolormesh(
        np.linspace(0, max_x, xPx + 1),
        np.linspace(0, max_y, yPx + 1),
        topoToPlot,
        cmap=colormap,
    )

    # Colorbar stuff
    cbar = fig.colorbar(XYmap, shrink=0.9, pad=0.05, aspect=15)
    cbar.ax.yaxis.set_ticks_position("both")
    cbar.ax.tick_params(axis="y", direction="in")
    return ax

def visualize_statistics(fig, peak_pos):
        # Plot the in a new figure
        ax = fig.add_subplot(1, 1, 1)
        Egap = peak_pos[:,1]-peak_pos[:,0]
        center_gap = Egap/2+peak_pos[:,0]
        # ax.plot(Egap*1000, center_gap*1000, 'x')
        # ax.set_ylabel("Gap position (meV)")
        ax.plot(Egap*1000, peak_pos[:,0]*1000, 'x')
        ax.set_xlabel("Gap Value (meV)")
        ax.set_ylabel("Lower peak position (meV)")

def visualize_statistics_tests(fig, peak_pos):
        # Plot the in a new figure
        ax = fig.add_subplot(1, 1, 1)
        RefEgap = 26 #meV
        Refcenter_gap = -100 #meV
        d=0.422 #nm
        delta = 0.33
        eVg = 33 #meV
        Egap = peak_pos[:,1]*1000-peak_pos[:,0]*1000
        center_gap = Egap/2+peak_pos[:,0]*1000
        deltad = ((Refcenter_gap-center_gap)+d*(RefEgap-Egap)/delta)/((RefEgap-Egap)/delta-center_gap/d)
        deltad1 = (RefEgap-Egap)*d**2/(eVg-d*(RefEgap-Egap))
        deltad2 = ((eVg-center_gap)/(eVg-Refcenter_gap)-1)*d
        ax.plot(deltad1,deltad2, 'x')
        ax.set_xlabel("deltad1 (nm)")
        ax.set_ylabel("deltad1 (nm)")
        # ax.set_ylabel("Gap value (meV)")

def visualize_statistics_topo(fig, peak_pos, Zarray):
        # Plot the in a new figure
        ax = fig.add_subplot(1, 1, 1)
        ax2 = ax.twinx()
        Egap = peak_pos[:,1]-peak_pos[:,0]
        center_gap = Egap/2+peak_pos[:,0]
        ax.plot(Zarray, peak_pos[:,0]*1000, 'x')
        # ax.plot(Zarray, center_gap, 'x')
        ax2.plot(Zarray, Egap*1000, 'x', color='k')
        ax2.set_ylabel("Gap Value (meV)", color='k')
        ax.set_ylabel("Lower peak position (meV)")
        ax.set_xlabel("Z (nm)")