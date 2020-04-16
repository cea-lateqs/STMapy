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
