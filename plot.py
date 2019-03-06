from ROOT import *


gInterpreter.LoadFile("GaussSeidel_Manuel_plot.cpp")



n = 4;

length = 101;

binsWidth = 5./(length - 1)

stencil = Stencil(n, length, length, 10, 10);


matrix = MMatrix(stencil.get_values());



timeLoop = timeloop(stencil, matrix);


c = TCanvas()

gStyle.SetOptStat(0);

h = TH2D("pot", "potential;x;y;#phi", length, -binsWidth, 10 + binsWidth, length, -binsWidth, 10 + binsWidth)


for x in range(0, length):
    
    for y in range(0, length):
        
        xx = 10.*float(x) / float(length - 1)
        
        yy = 10.*float(y) / float(length - 1)
        
        binnr = h.FindBin(xx, yy)
        
        h.SetBinContent(binnr, timeLoop[x][y])



h.Draw("SURF2")



c.Update()

c.SaveAs("potentialT.pdf")

