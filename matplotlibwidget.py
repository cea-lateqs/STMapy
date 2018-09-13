from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as Canvas
from matplotlib.figure import Figure


class MatplotlibWidget(Canvas):
    def __init__(self, parent=None, dpi=100):
        super(MatplotlibWidget, self).__init__(Figure())

        self.setParent(parent)
        self.figure = Figure(dpi=dpi)
        self.canvas = Canvas(self.figure)
