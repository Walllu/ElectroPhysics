from ROOT import *



#c = TCanvas()

#h = TH2D("pot", "potential;x;y;#phi", 4, -1., 1., 4, -1., 1.)



l = [[1,2,3,2],[2,3,3,2],[2,3,3,3],[1,2,4,2]]







#h.Draw("CONT")



#c.Update()

#c.SaveAs("potential.pdf")



gInterpreter.LoadFile("./GaussSeidel_plot.cpp")



n = 2;

length = 101;

stencil = Stencil(n, length, length, 10, 10);



matrix = MMatrix(stencil.get_values());



timeLoop = timeloop(stencil, matrix);



c = TCanvas()






h = TH2D("pot", "potential;x;y;#phi", length, -0.5, 5.5, length, -0.5, 5.5)



for x in range(0, length):

    for y in range(0, length):

        xx = float(x) / 10.

        yy = float(y) / 10.

        binnr = h.FindBin(xx, yy)

        h.SetBinContent(binnr, timeLoop[x][y])



h.Draw("SURF2")



c.Update()

c.SaveAs("potential2.pdf")

