from ROOT import *


gInterpreter.LoadFile("./GaussSeidel_plot.cpp")



n = 2;

length = 51;

stencil = Stencil(n, length, length, 10, 10);



matrix = MMatrix(stencil.get_values());



timeLoop = timeloop(stencil, matrix);



c = TCanvas()






h = TH2D("pot", "potential;x;y;#phi", length, -0.1, 10.1, length, -0.1, 10.1)



for x in range(0, length):

    for y in range(0, length):

        xx = 10.*float(x) / float(length - 1)

        yy = 10.*float(y) / float(length - 1)

        binnr = h.FindBin(xx, yy)

        h.SetBinContent(binnr, timeLoop[x][y])



h.Draw("SURF2")



c.Update()

c.SaveAs("potential2.pdf")

