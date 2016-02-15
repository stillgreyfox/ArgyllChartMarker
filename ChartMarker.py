#!venv/bin/python
# #!/usr/bin/env python
#
# Copyright 2005-2010, Karljohan Lundin Palmerius
#
# ChartMarker is a Graphical User Interface (GUI) for use together
# with the Argyll Color Management System for finding the coordinates
# of the colour chart (target) in a photo. These can then be used
# together with the scanin tool to extract colour information that can
# then be used with the colprof tool to create a colour profile (ICC
# profile).
#
# 2010-03-15, Version 1.1
# Added correct projective mapping that seems to be identical to the
# one used in the scanin tool.
#
# 2010-01-26, Version 1.0
# First version with basic functionality. The code looks like an ugly
# hack, but works. Requires Python and wxPython.
#
#
# ChartMarker is distributed under GNU GPL version 2.
#
# ChartMarker is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# ChartMarker is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ChartMarker; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
# USA

import math
import sys
import wx
import numpy


class support:
    have = {}

    @staticmethod
    def fetch(libname, warning=None):
        try:
            lib = __import__(libname, globals(), locals(), [], -1)
            globals()[libname] = lib
            support.have[libname] = True
            return True
        except:
            support.have[libname] = False
            if warning is not None:
                sys.stderr.write(warning)
            return False
'''
if not support.fetch("wx", "Error: wxPython is required for this tool!"):
    exit(1)
support.fetch("numpy",
              "Warning: NumPy is needed for drawing target - function disabled")
'''

from optparse import OptionParser

parser = OptionParser(usage="usage: %prog [options] <tiff image>")
parser.set_defaults(handle_size=20,
                    border_size=20,
                    show_area=False,
                    area_shrink=2,
                    scale=1.0,
                    fmt="x1,y1,x2,y2,x3,y3,x4,y4",
                    returnValue=0,
                    target=None)
parser.add_option("-H", "--handle", type="int", dest="handle_size",
                  help="Set size of polygon handles",
                  metavar="PIX")
parser.add_option("-B", "--border", type="float", dest="border_size",
                  help="Set size of border that activates " +
                       "scrolling during dragging",
                  metavar="PIX")
parser.add_option("-F", "--format", dest="fmt",
                  help="The output format of the chart polygon. " +
                       "Default is x1,y1,x2,y2,x3,y3,x4,y4",
                  metavar="STRING")
# if support.have["numpy"]:
if True:
    parser.add_option("-t", "--target", dest="target",
                      help="Set target to show chart structure",
                      metavar="FILE")
    parser.add_option("-a", "--show-area",
                      action="store_true", dest="show_area",
                      help="Show the area that will be sampled")
parser.add_option("-s", "--scale",
                  type="float", dest="scale",
                  help="Scale the image to this ratio")

(opts, cmd_args) = parser.parse_args()

if len(cmd_args) < 1:
    parser.print_help()
    exit(1)


class Vec3:

    def __init__(self):
        self.x = 0
        self.y = 0
        self.z = 0

    def __mult__(self, o2):
        print(self, o2)


class Mark:

    def __init__(self, isSample):
        self.Nx = 0
        self.Ny = 0
        self.origin = wx.Point(0, 0)
        self.size = wx.Point(0, 0)
        self.offset = wx.Point(0, 0)
        self.isSample = isSample


def dist(a, b):
    return math.sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y))


def countSteps(str1, str2):
    if len(str1) != len(str2):
        n1 = int(str1)
        n2 = int(str2)
        return n2 - n1 + 1
    count = 0
    for (c1, c2) in zip(str1, str2):
        count = 10 * count + (ord(c2) - ord(c1))
    return count + 1


class ChartMarker(wx.Frame):

    def __init__(self):
        wx.Frame.__init__(self, None, -1, "Chart Marker",
                          wx.DefaultPosition, wx.Size(700, 500))

        wx.EvtHandler.Bind(self, wx.EVT_CLOSE, self.OnClose)
        wx.EvtHandler.Bind(self, wx.EVT_PAINT, self.OnPaint)
        wx.EvtHandler.Bind(self, wx.EVT_MOTION, self.OnMouseMotion)
        wx.EvtHandler.Bind(self, wx.EVT_LEFT_DOWN, self.OnMousePress)
        wx.EvtHandler.Bind(self, wx.EVT_LEFT_UP, self.OnMousePress)

        self.image = wx.Image(cmd_args[0], wx.BITMAP_TYPE_ANY)
        if math.fabs(1 - opts.scale) > 0.001:
            self.image = \
                self.image.Scale(int(self.image.GetWidth() * opts.scale),
                                 int(self.image.GetHeight() * opts.scale))
        self.bitmap = wx.Bitmap(self.image)
        self.image_position = wx.Point(0, 0)
        self.image_size = [self.image.GetWidth(), self.image.GetHeight()]
        self.target_area = [wx.Point(0, 0),
                            wx.Point(self.image_size[0], 0),
                            wx.Point(self.image_size[0], self.image_size[1]),
                            wx.Point(0, self.image_size[1])]
        '''
        size = [self.GetSize().GetWidth(), self.GetSize().GetHeight()]
        self.target_area = [wx.Point(1 * size[0] / 4, 1 * size[1] / 4),
                            wx.Point(3 * size[0] / 4, 1 * size[1] / 4),
                            wx.Point(3 * size[0] / 4, 3 * size[1] / 4),
                            wx.Point(1 * size[0] / 4, 3 * size[1] / 4)]
        '''

        '''
        # this is now done after instantiating class. wouldn't work this way
        self.SetSize(wx.Size(min(1920, max(200, self.image_size[0])),
                             min(1080, max(200, self.image_size[1]))))
        '''

        self.mouseState_drag = False
        self.selected_handle = None
        self.button_status = 0

        self.fiducials = None
        self.marks = []
        if opts.target is not None:
            with open(opts.target, "r") as tin:
                for line in tin:
                    if line.find("BOXES") > -1:
                        found_boxes = True
                        break
                    else:
                        found_boxes = False

                if found_boxes:
                    print(line)
                    line = tin.__next__().strip()
                    while len(line) > 0:
                        print(line)
                        print('here')
                        T = line.split()
                        if T[0] == "F":
                            self.fiducials = [
                                wx.RealPoint(float(T[3]), float(T[4])),
                                wx.RealPoint(float(T[5]), float(T[6])),
                                wx.RealPoint(float(T[7]), float(T[8])),
                                wx.RealPoint(float(T[9]), float(T[10]))]
                            print(self.fiducials)

                        elif (T[0] == "X" or
                              T[0] == "Y" or
                              T[0] == "D"):
                            mark = Mark(T[0] != "D")
                            if T[1] == "_" or \
                               T[1] == "ALL" or \
                               T[1] == "MARK":
                                mark.Nx = 1
                            else:
                                mark.Nx = countSteps(T[1], T[2])

                            if T[3] == "_" or T[3] == "ALL" or T[3] == "MARK":
                                mark.Ny = 1
                            else:
                                mark.Ny = countSteps(T[3], T[4])

                            mark.size = wx.RealPoint(float(T[5]), float(T[6]))
                            mark.origin = \
                                wx.RealPoint(float(T[7]), float(T[8]))
                            mark.offset = \
                                wx.RealPoint(float(T[9]), float(T[10]))

                            self.marks.append(mark)

                        else:
                            sys.stderr.write("Parse error in target file.")

                        line = tin.__next__().strip()

                    maxX = 0.0
                    maxY = 0.0
                    for mark in self.fiducials:
                        maxX = max(maxX, mark.x)
                        maxY = max(maxY, mark.y)
                    for mark in self.marks:
                        maxX = max(maxX, mark.origin.x)
                        maxY = max(maxY, mark.origin.y)
                        maxX = max(maxX, mark.origin.x +
                                   mark.Nx * mark.offset.x)
                        maxY = max(maxY, mark.origin.y +
                                   mark.Ny * mark.offset.y)

                    maxA = max(maxX, maxY)
                    maxX = maxY = maxA

                    for mark in self.fiducials:
                        mark.x = mark.x / maxX
                        mark.y = mark.y / maxY
                    for mark in self.marks:
                        mark.origin.x = mark.origin.x / maxX
                        mark.origin.y = mark.origin.y / maxY
                        mark.offset.x = mark.offset.x / maxX
                        mark.offset.y = mark.offset.y / maxY
                        mark.size.x = mark.size.x / maxX
                        mark.size.y = mark.size.y / maxY

                    self.UpdateMarkerData()

                else:
                    sys.stderr.write("Could not find boxes in target file.")

                for line in tin:
                    if line.find("BOX_SHRINK") > -1:
                        opts.area_shrink = float(line.split()[1])
                        opts.area_shrink = float(line.split()[1])
                        break
                opts.area_shrink_x = opts.area_shrink / maxX
                opts.area_shrink_y = opts.area_shrink / maxY

        self.Fit()

    def UpdateMarkerData(self):
        if self.fiducials is not None:
            U = self.fiducials
            X = self.target_area

            M = numpy.array(
                [[U[0].x, U[0].y, 1, 0, 0, 0, -U[0].x*X[0].x, -U[0].y*X[0].x],
                 [U[1].x, U[1].y, 1, 0, 0, 0, -U[1].x*X[1].x, -U[1].y*X[1].x],
                 [U[2].x, U[2].y, 1, 0, 0, 0, -U[2].x*X[2].x, -U[2].y*X[2].x],
                 [U[3].x, U[3].y, 1, 0, 0, 0, -U[3].x*X[3].x, -U[3].y*X[3].x],
                 [0, 0, 0, U[0].x, U[0].y, 1, -U[0].x*X[0].y, -U[0].y*X[0].y],
                 [0, 0, 0, U[1].x, U[1].y, 1, -U[1].x*X[1].y, -U[1].y*X[1].y],
                 [0, 0, 0, U[2].x, U[2].y, 1, -U[2].x*X[2].y, -U[2].y*X[2].y],
                 [0, 0, 0, U[3].x, U[3].y, 1, -U[3].x*X[3].y, -U[3].y*X[3].y]])

            D = numpy.array(
                [X[0].x,
                 X[1].x,
                 X[2].x,
                 X[3].x,
                 X[0].y,
                 X[1].y,
                 X[2].y,
                 X[3].y])

            Y = numpy.linalg.solve(M, D)

            self.Msd = numpy.array([[Y[0], Y[1], Y[2]],
                                    [Y[3], Y[4], Y[5]],
                                    [Y[6], Y[7], 1]])

    def GetABCDPoint(self, pos):
        U = numpy.array([pos[0], pos[1], 1])
        X = numpy.dot(self.Msd, U)
        return wx.Point(int(X[0] / X[2]), int(X[1] / X[2]))

    def OnClose(self, event):
        self.OnExit(False)

    def OnExit(self, success):
        if success:
            fmt = opts.fmt
            fmt = fmt.replace("x1", str(int(self.target_area[0].x/opts.scale)))
            fmt = fmt.replace("y1", str(int(self.target_area[0].y/opts.scale)))
            fmt = fmt.replace("x2", str(int(self.target_area[1].x/opts.scale)))
            fmt = fmt.replace("y2", str(int(self.target_area[1].y/opts.scale)))
            fmt = fmt.replace("x3", str(int(self.target_area[2].x/opts.scale)))
            fmt = fmt.replace("y3", str(int(self.target_area[2].y/opts.scale)))
            fmt = fmt.replace("x4", str(int(self.target_area[3].x/opts.scale)))
            fmt = fmt.replace("y4", str(int(self.target_area[3].y/opts.scale)))
            print(fmt)
        else:
            opts.returnValue = -1
        self.Destroy()

    def OnMouseMotion(self, event):
        mousepos = event.GetPosition()
        needsRefresh = False
        if self.mouseState_drag:
            if self.selected_handle is None:
                if self.button_status == 0:
                    self.image_position = self.image_position + mousepos - \
                        self.mouseState_zpos
            else:
                if self.selected_handle < 4:  # Corner handle
                    dx = mousepos - self.mouseState_zpos
                    self.target_area[self.selected_handle] += dx
                elif self.selected_handle < 8:  # Edge handle
                    dx = mousepos - self.mouseState_zpos
                    self.target_area[self.selected_handle - 4] += dx
                    self.target_area[(self.selected_handle + 1) % 4] += dx
                elif self.selected_handle < 9:  # Middle handle
                    dx = mousepos - self.mouseState_zpos
                    for idx in range(4):
                        self.target_area[idx] += dx
                else:  # Rotation handle
                    v = mousepos - self.rotation_center - self.image_position
                    angle = math.atan2(v.y, v.x) - self.rotation_start
                    for idx in range(4):
                        X = [self.target_R[idx] *
                             math.cos(self.target_A[idx] + angle),
                             self.target_R[idx] *
                             math.sin(self.target_A[idx] + angle)]
                        self.target_area[idx] = \
                            wx.Point(int(X[0]), int(X[1])) + \
                            self.rotation_center
                self.UpdateMarkerData()

                size = self.GetSize()
                if mousepos.x < opts.border_size:
                    self.image_position = self.image_position + \
                        wx.Point(opts.border_size - mousepos.x, 0)
                    mousepos = mousepos + \
                        wx.Point(opts.border_size - mousepos.x, 0)
                if mousepos.y < opts.border_size:
                    self.image_position = self.image_position + \
                        wx.Point(0, opts.border_size - mousepos.y)
                    mousepos = mousepos + \
                        wx.Point(0, opts.border_size - mousepos.y)
                if mousepos.x > size.GetWidth() - opts.border_size:
                    self.image_position = self.image_position + \
                        wx.Point(size.GetWidth() - opts.border_size -
                                 mousepos.x, 0)
                    mousepos = mousepos + \
                        wx.Point(size.GetWidth() - opts.border_size -
                                 mousepos.x, 0)
                if mousepos.y > size.GetHeight() - opts.border_size:
                    self.image_position = self.image_position + \
                        wx.Point(0, size.GetHeight() - opts.border_size -
                                 mousepos.y)
                    mousepos = mousepos + \
                        wx.Point(0, size.GetHeight() - opts.border_size -
                                 mousepos.y)
            self.mouseState_zpos = mousepos
            needsRefresh = True

        else:

            if opts.handle_size < mousepos.x and \
                    mousepos.x < 2 * opts.handle_size and \
                    opts.handle_size < mousepos.y and \
                    mousepos.y < 2 * opts.handle_size:
                if self.button_status != 1:
                    needsRefresh = True
                self.button_status = 1
            else:
                if self.button_status != 0:
                    needsRefresh = True
                self.button_status = 0

            old_selected_handle = self.selected_handle
            self.selected_handle = None
            for idx in range(len(self.target_area)):
                corner = self.target_area[idx]
                d = mousepos - corner - self.image_position
                if -opts.handle_size / 2 < d.x and d.x < opts.handle_size and \
                   -opts.handle_size / 2 < d.y and d.y < opts.handle_size:
                    self.selected_handle = idx
                    break

                corner2 = self.target_area[(idx + 1) % 4]
                d = mousepos - \
                    wx.Point((corner.x + corner2.x) / 2,
                             (corner.y + corner2.y) / 2) - self.image_position
                if -opts.handle_size / 2 < d.x and d.x < opts.handle_size and \
                   -opts.handle_size / 2 < d.y and d.y < opts.handle_size:
                    self.selected_handle = idx + 4
                    break

            center = wx.Point(0, 0)
            for corner in self.target_area:
                center = center + corner
            center.x /= 4
            center.y /= 4

            d = mousepos - center - self.image_position
            if -opts.handle_size / 2 < d.x and d.x < opts.handle_size and \
                    -opts.handle_size / 2 < d.y and d.y < opts.handle_size:
                self.selected_handle = 8

            center2 = center + center + \
                self.target_area[0] + self.target_area[1]
            center2.x /= 4
            center2.y /= 4

            d = mousepos - center2 - self.image_position
            if -opts.handle_size / 2 < d.x and d.x < opts.handle_size and \
                    -opts.handle_size / 2 < d.y and d.y < opts.handle_size:

                self.selected_handle = 9

                v = mousepos - center - self.image_position
                self.rotation_center = center
                self.rotation_start = math.atan2(v.y, v.x)
                self.target_R = []
                self.target_A = []
                for idx in range(4):
                    corner = self.target_area[idx] - center
                    R = math.sqrt(corner.x * corner.x + corner.y * corner.y)
                    self.target_R.append(R)
                    A = math.atan2(corner.y, corner.x)
                    self.target_A.append(A)

            if old_selected_handle != self.selected_handle:
                needsRefresh = True

        if needsRefresh:
            self.Refresh()

    def OnMousePress(self, event):
        self.mouseState_drag = event.ButtonDown()
        self.mouseState_zpos = event.GetPosition()
        if event.ButtonDown():
            if self.button_status == 1:
                self.button_status = 2
                self.Refresh()
        else:
            if self.button_status == 2:
                mousepos = event.GetPosition()
                if opts.handle_size < mousepos.x and \
                        mousepos.x < 2 * opts.handle_size and \
                        opts.handle_size < mousepos.y and \
                        mousepos.y < 2 * opts.handle_size:
                    self.OnExit(True)
                else:
                    self.button_status = 0
                    self.Refresh()

    def OnPaint(self, event):
        g = wx.PaintDC(self)
        g.Clear()

        if self.bitmap is None:
            return

        g.SetLogicalFunction(wx.COPY)
        g.DrawBitmap(self.bitmap,
                     self.image_position.x, self.image_position.y,
                     False)

        g.SetPen(wx.GREEN_PEN)
        g.SetBrush(wx.TRANSPARENT_BRUSH)
        if self.fiducials is not None:
            for mark in self.marks:
                if opts.show_area:
                    shrunksizex = mark.size.x - opts.area_shrink_x
                    shrunksizey = mark.size.y - opts.area_shrink_y

                posx = mark.origin.x
                posy = mark.origin.y
                for py in range(mark.Ny):
                    for px in range(mark.Nx):
                        corners = [[posx,
                                    posy],
                                   [posx + mark.size.x,
                                    posy],
                                   [posx + mark.size.x,
                                    posy + mark.size.y],
                                   [posx,
                                    posy + mark.size.y]]
                        points = []
                        for corner in corners:
                            points.append(self.GetABCDPoint(corner) +
                                          self.image_position)
                        g.DrawPolygon(points)

                        if opts.show_area and mark.isSample:
                            corners = [[posx + opts.area_shrink_x,
                                        posy + opts.area_shrink_y],
                                       [posx + shrunksizex,
                                        posy + opts.area_shrink_y],
                                       [posx + shrunksizex,
                                        posy + shrunksizey],
                                       [posx + opts.area_shrink_x,
                                        posy + shrunksizey]]
                            points = []
                            for corner in corners:
                                points.append(self.GetABCDPoint(corner) +
                                              self.image_position)
                            g.DrawPolygon(points)
                            g.DrawLine(points[0].x, points[0].y,
                                       points[2].x, points[2].y)
                            g.DrawLine(points[1].x, points[1].y,
                                       points[3].x, points[3].y)

                        posx = posx + mark.offset.x

                    posx = mark.origin.x
                    posy = posy + mark.offset.y

        g.SetPen(wx.RED_PEN)
        g.DrawPolygon(self.target_area,
                      self.image_position.x, self.image_position.y)

        g.SetPen(wx.BLACK_PEN)
        g.SetBrush(wx.WHITE_BRUSH)

        for idx in range(len(self.target_area)):
            corner = self.target_area[idx]
            g.DrawRectangle(corner.x +
                            self.image_position.x - opts.handle_size / 2,
                            corner.y +
                            self.image_position.y - opts.handle_size / 2,
                            opts.handle_size, opts.handle_size)
            corner2 = self.target_area[(idx + 1) % 4]
            g.DrawRectangle((corner.x + corner2.x) / 2 +
                            self.image_position.x - opts.handle_size / 2,
                            (corner.y + corner2.y) / 2 +
                            self.image_position.y - opts.handle_size / 2,
                            opts.handle_size, opts.handle_size)

        center = wx.Point(0, 0)
        for corner in self.target_area:
            center = center + corner
        center.x /= 4
        center.y /= 4
        g.DrawRectangle(center.x +
                        self.image_position.x - opts.handle_size / 2,
                        center.y +
                        self.image_position.y - opts.handle_size / 2,
                        opts.handle_size, opts.handle_size)

        center2 = center + center + self.target_area[0] + self.target_area[1]
        center2.x /= 4
        center2.y /= 4
        g.DrawRectangle(center2.x +
                        self.image_position.x - opts.handle_size / 2,
                        center2.y +
                        self.image_position.y - opts.handle_size / 2,
                        opts.handle_size, opts.handle_size)

        if self.selected_handle is not None:
            g.SetBrush(wx.BLACK_BRUSH)
            if self.selected_handle < 4:
                corner = self.target_area[self.selected_handle]
                g.DrawRectangle(corner.x +
                                self.image_position.x - opts.handle_size / 2,
                                corner.y +
                                self.image_position.y - opts.handle_size / 2,
                                opts.handle_size, opts.handle_size)
            elif self.selected_handle < 8:
                corner1 = self.target_area[self.selected_handle - 4]
                corner2 = self.target_area[(self.selected_handle - 3) % 4]
                corner = corner1 + corner2
                g.DrawRectangle(corner.x / 2 +
                                self.image_position.x - opts.handle_size / 2,
                                corner.y / 2 +
                                self.image_position.y - opts.handle_size / 2,
                                opts.handle_size, opts.handle_size)
            elif self.selected_handle < 9:
                g.DrawRectangle(center.x +
                                self.image_position.x - opts.handle_size / 2,
                                center.y +
                                self.image_position.y - opts.handle_size / 2,
                                opts.handle_size, opts.handle_size)
            else:
                g.DrawRectangle(center2.x +
                                self.image_position.x - opts.handle_size / 2,
                                center2.y +
                                self.image_position.y - opts.handle_size / 2,
                                opts.handle_size, opts.handle_size)

            g.SetBrush(wx.TRANSPARENT_BRUSH)

        g.SetPen(wx.BLACK_PEN)
        if self.button_status == 0:
            g.SetBrush(wx.WHITE_BRUSH)
        elif self.button_status == 1:
            g.SetBrush(wx.GREEN_BRUSH)
        else:
            g.SetBrush(wx.RED_BRUSH)
        g.DrawRectangle(opts.handle_size, opts.handle_size,
                        opts.handle_size, opts.handle_size)


class MarkerApp(wx.App):

    def __init__(self, id):
        wx.App.__init__(self, id)

    def OnInit(self):
        frame = ChartMarker()
        frame.SetSize(wx.Size(min(1920, max(200, frame.image_size[0])),
                              min(1080, max(200, frame.image_size[1] + 27))))
        frame.Show(True)
        self.SetTopWindow(frame)
        print(frame.image_size)
        return True


def main():
    app = MarkerApp(0)
    app.MainLoop()
    return opts.returnValue

if __name__ == "__main__":
    sys.exit(main())
