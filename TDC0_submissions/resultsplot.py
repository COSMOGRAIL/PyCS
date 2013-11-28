
import pycs
import glob

dtfiles = glob.glob("*.dt")

estimateslists = [e for dtfile in dtfiles for e in pycs.tdc.est.readsubmission(dtfile)]

pycs.tdc.est.bigplot(estimateslists, plotpath = "compa.pdf")


