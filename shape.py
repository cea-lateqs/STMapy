#-*- coding: utf-8 -*-
from matplotlib import patches


class Shape:
    """ Shape class that is used to draw the shapes when clicking on the map """
    def __init__(self, event, fig_map, fig_topo, color, ratioX, ratioY):
        self.xi = int(event.xdata)
        self.yi = int(event.ydata)
        self.xf = int(event.xdata)
        self.yf = int(event.ydata)
        self.d = int(event.button)
        self.color = color
        self.map1 = fig_map
        self.map2 = fig_topo
        self.ratioX = ratioX
        self.ratioY = ratioY
        self.create()

    def changeToPt(self):
        self.d = 0
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
        if self.d == 1:
            self.shape1, = self.map1.axes[0].plot([self.xi + 0.5, self.xf + 0.5],
                                                  [self.yi + 0.5, self.yf + 0.5], 'k--')
            if self.map2 != 0:
                self.shape2, = self.map2.axes[0].plot([(self.xi + 0.5) * self.ratioX, (self.xf + 0.5) * self.ratioX],
                                                      [(self.yi + 0.5) * self.ratioY, (self.yf + 0.5) * self.ratioY], 'k--')
        elif self.d == 3:
            self.shape1 = self.map1.axes[0].add_patch(patches.Rectangle((self.xi, self.yi),
                                                                        self.xf - self.xi,
                                                                        self.yf - self.yi, color=self.color, alpha=0.5))
            if self.map2 != 0:
                self.shape2 = self.map2.axes[0].add_patch(patches.Rectangle((self.xi * self.ratioX, self.yi * self.ratioY),
                                                                            (self.xf - self.xi) * self.ratioX,
                                                                            (self.yf - self.yi) * self.ratioY, color=self.color, alpha=0.5))
        elif self.d == 0:
            self.shape1, = self.map1.axes[0].plot([self.xi + 0.5, self.xf + 0.5], [self.yi + 0.5, self.yf + 0.5],
                                                  c=self.color, ls='None', marker='o')
            if self.map2 != 0:
                self.shape2, = self.map2.axes[0].plot([(self.xi + 0.5) * self.ratioX, (self.xf + 0.5) * self.ratioX],
                                                      [(self.yi + 0.5) * self.ratioY, (self.yf + 0.5) * self.ratioY],
                                                      c=self.color, ls='None', marker='o')
        self.drawMaps()

    def recreate(self, fig_map, fig_topo):
        self.updateMaps(fig_map, fig_topo)
        self.create()

    def draw(self):
        if self.d == 1 or self.d == 0:
            self.shape1.set_data([self.xi + 0.5, self.xf + 0.5], [self.yi + 0.5, self.yf + 0.5])
            if self.map2 != 0:
                self.shape2.set_data([(self.xi + 0.5) * self.ratioX, (self.xf + 0.5) * self.ratioX],
                                     [(self.yi + 0.5) * self.ratioY, (self.yf + 0.5) * self.ratioY])
        elif self.d == 3:
            self.shape1.set_xy([self.xi, self.yi])
            self.shape1.set_width(self.xf - self.xi)
            self.shape1.set_height(self.yf - self.yi)
            if self.map2 != 0:
                self.shape2.set_xy([self.xi * self.ratioX, self.yi * self.ratioY])
                self.shape2.set_width(self.xf * self.ratioX - self.xi * self.ratioX)
                self.shape2.set_height(self.yf * self.ratioY - self.yi * self.ratioY)
        self.drawMaps()

    def drawMaps(self):
        if self.map1 != 0: self.map1.canvas.draw()
        if self.map2 != 0: self.map2.canvas.draw()

    def updateMaps(self, fig_map, fig_topo):
        self.map1 = fig_map
        self.map2 = fig_topo
        self.shape1.set_axes(self.map1.axes[0])
        if self.map2 != 0:
            self.shape2.set_axes(self.map2.axes[0])

    def remove(self):
        self.shape1.remove()
        del self.shape1
        if self.map2 != 0:
            self.shape2.remove()
            del self.shape2
        self.drawMaps()
