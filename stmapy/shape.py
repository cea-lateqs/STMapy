# -*- coding: utf-8 -*-
from matplotlib import patches


def generateShape(event, fig_map, fig_topo, color, ratio_x, ratio_y):
    if event.button == 0:
        return Dot(
            int(event.xdata),
            int(event.ydata),
            fig_map,
            fig_topo,
            color,
            ratio_x,
            ratio_y,
        )
    elif event.button == 1:
        return Line(
            int(event.xdata),
            int(event.ydata),
            fig_map,
            fig_topo,
            color,
            ratio_x,
            ratio_y,
        )
    elif event.button == 3:
        return Rectangle(
            int(event.xdata),
            int(event.ydata),
            fig_map,
            fig_topo,
            color,
            ratio_x,
            ratio_y,
        )
    else:
        raise ValueError("Event button not recognized in generateShape !")


def changeToDot(shape):
    return Dot(
        shape.xi,
        shape.yi,
        shape.map1,
        shape.map2,
        shape.color,
        shape.ratio_x,
        shape.ratio_y,
    )


class Shape:
    """Shape class that is used to draw the shapes when clicking on the map"""

    def __init__(
        self, xi, yi, fig_map, fig_topo, color, ratio_x, ratio_y, xf=None, yf=None
    ):
        self.xi = xi
        self.yi = yi
        self.xf = xi if xf is None else xf
        self.yf = yi if yf is None else yf
        self.shape1 = None
        self.shape2 = None
        self.color = color
        self.map1 = fig_map
        self.map2 = fig_topo
        self.ratio_x = ratio_x
        self.ratio_y = ratio_y
        self.create()

    def update(self, event):
        if event.xdata is not None and event.ydata is not None:
            self.xf = int(event.xdata)
            self.yf = int(event.ydata)
            self.draw()

    def forceUpdate(self, x, y, xp, yp):
        self.xi = x
        self.yi = y
        self.xf = xp
        self.yf = yp
        self.draw()

    def create(self):
        raise NotImplementedError("Create not implemented !")

    def recreate(self, fig_map, fig_topo):
        self.updateMaps(fig_map, fig_topo)
        self.create()

    def draw(self):
        raise NotImplementedError("Draw not implemented !")

    def drawMaps(self):
        if self.map1 is not None:
            self.map1.canvas.draw()
        if self.map2 is not None:
            self.map2.canvas.draw()

    def updateMaps(self, fig_map, fig_topo):
        self.map1 = fig_map
        self.map2 = fig_topo

    def remove(self):
        self.shape1.remove()
        del self.shape1
        if self.map2 is not None:
            self.shape2.remove()
            del self.shape2
        self.drawMaps()


class Dot(Shape):
    def create(self):
        (self.shape1,) = self.map1.axes[0].plot(
            [self.xi + 0.5, self.xf + 0.5],
            [self.yi + 0.5, self.yf + 0.5],
            c=self.color,
            ls="None",
            marker="o",
        )
        if self.map2 is not None:
            (self.shape2,) = self.map2.axes[0].plot(
                [(self.xi + 0.5) * self.ratio_x, (self.xf + 0.5) * self.ratio_x],
                [(self.yi + 0.5) * self.ratio_y, (self.yf + 0.5) * self.ratio_y],
                c=self.color,
                ls="None",
                marker="o",
            )
        self.drawMaps()

    def draw(self):
        self.shape1.set_data(
            [self.xi + 0.5, self.xf + 0.5], [self.yi + 0.5, self.yf + 0.5]
        )
        if self.map2 is not None:
            self.shape2.set_data(
                [(self.xi + 0.5) * self.ratio_x, (self.xf + 0.5) * self.ratio_x],
                [(self.yi + 0.5) * self.ratio_y, (self.yf + 0.5) * self.ratio_y],
            )
        self.drawMaps()


class Line(Shape):
    def create(self):
        (self.shape1,) = self.map1.axes[0].plot(
            [self.xi + 0.5, self.xf + 0.5], [self.yi + 0.5, self.yf + 0.5], "k--"
        )
        if self.map2 is not None:
            (self.shape2,) = self.map2.axes[0].plot(
                [(self.xi + 0.5) * self.ratio_x, (self.xf + 0.5) * self.ratio_x],
                [(self.yi + 0.5) * self.ratio_y, (self.yf + 0.5) * self.ratio_y],
                "k--",
            )
        self.drawMaps()

    def draw(self):
        self.shape1.set_data(
            [self.xi + 0.5, self.xf + 0.5], [self.yi + 0.5, self.yf + 0.5]
        )
        if self.map2 is not None:
            self.shape2.set_data(
                [(self.xi + 0.5) * self.ratio_x, (self.xf + 0.5) * self.ratio_x],
                [(self.yi + 0.5) * self.ratio_y, (self.yf + 0.5) * self.ratio_y],
            )
        self.drawMaps()


class Rectangle(Shape):
    def create(self):
        self.shape1 = self.map1.axes[0].add_patch(
            patches.Rectangle(
                (self.xi, self.yi),
                self.xf - self.xi,
                self.yf - self.yi,
                color=self.color,
                alpha=0.5,
            )
        )
        if self.map2 is not None:
            self.shape2 = self.map2.axes[0].add_patch(
                patches.Rectangle(
                    (self.xi * self.ratio_x, self.yi * self.ratio_y),
                    (self.xf - self.xi) * self.ratio_x,
                    (self.yf - self.yi) * self.ratio_y,
                    color=self.color,
                    alpha=0.5,
                )
            )
        self.drawMaps()

    def draw(self):
        self.shape1.set_xy([self.xi, self.yi])
        self.shape1.set_width(self.xf - self.xi)
        self.shape1.set_height(self.yf - self.yi)
        if self.map2 is not None:
            self.shape2.set_xy([self.xi * self.ratio_x, self.yi * self.ratio_y])
            self.shape2.set_width(self.xf * self.ratio_x - self.xi * self.ratio_x)
            self.shape2.set_height(self.yf * self.ratio_y - self.yi * self.ratio_y)

        self.drawMaps()
