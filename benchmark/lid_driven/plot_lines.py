#!/usr/bin/env pvpython
"""Sample u/v/pressure along x = 0.5 and y = 0.5 for the lid-driven cavity and save as CSV.

vel_x/vel_y/pressure are stored as elemental (cell) variables, so a plain
PlotOverLine on the raw reader gives a piecewise-constant staircase (one flat
step per triangle). CellDatatoPointData averages them onto the mesh nodes
first, so the line plot is a smooth, linearly-interpolated profile instead.

This must run under ParaView's own Python, not a system python3 -- paraview.simple
is not a pip package here. Point this at whichever ParaView install you use
(they are not on PATH by default on macOS):

Usage:
    /Applications/ParaView-6.1.1.app/Contents/bin/pvpython benchmark/lid_driven/plot_lines.py \\
        benchmark/lid_driven/cavity_Re100_out.e benchmark/lid_driven
"""
import sys

from paraview.simple import CellDatatoPointData, OpenDataFile, PlotOverLine, SaveData

RESOLUTION = 100  # sample points along each line

LINES = {
    # name -> (Point1, Point2)
    "u_mid": ((0.5, 0.0, 0.0), (0.5, 1.0, 0.0)),  # vertical line at x = 0.5
    "v_mid": ((0.0, 0.5, 0.0), (1.0, 0.5, 0.0)),  # horizontal line at y = 0.5
}


def main(exodus_path, out_dir):
    reader = OpenDataFile(exodus_path)
    reader.UpdatePipeline()
    # Without this the pipeline evaluates at the first timestep, not the converged
    # one -- SaveData/PlotOverLine otherwise silently sample the transient start-up
    # state instead of steady state.
    last_time = reader.TimestepValues[-1]

    c2p = CellDatatoPointData(Input=reader)

    for name, (point1, point2) in LINES.items():
        line = PlotOverLine(Input=c2p, Point1=point1, Point2=point2, Resolution=RESOLUTION)
        line.UpdatePipeline(last_time)
        out_path = "%s/%s.csv" % (out_dir, name)
        SaveData(out_path, proxy=line)
        print("wrote %s" % out_path)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
