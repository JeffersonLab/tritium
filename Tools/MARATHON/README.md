# MARATHON DIS Model and Rate Calculator
##  --- created on 11/10/2016 by Zhihong Ye
##  --- also can be found on www.jlab.org/~yez/Work/DIS/

## Initial Commit:
The main code is "DIS.h", which includes three different models, "F2ALLM97", "F1F2IN09" and "CTEQ-JLAB", to calculate F1 and F2 for protons, neutrons and deutrons. The DIS cross section can be calculated with F2.

In addition to these, there is another model, Peter Bosted's DIS model, which more focuses on the residenance region. When compared with other models, it gives a factor between 2 and 3 larger cross section values. The reason is still unknown.

"F2ALLM97" only provides F2p calculations. A "F2n/F2p" ratio was used to get the F2n. Hence, the F2d calculation is not precise as no Fermi smearing is embedded. When compring with world data, this model works very well at high Q2.

"F1F2IN09" gives F2p, F2n and F2d. It directly calls Eric Christy and Peter Bosted's 2009-Fit code (in Fortran). It describes CLAS's data at high-x very well and gives consistance results as "F2ALLM97" at certain Q2 range.
F2d goes to negative at Q2 roughly larger than 10 GeV. Overall, it may work well at moderate Q2 range like Jlab-6eV data.

"CTEQ-JLab" provides PDFs. I simply combined the quark flavor and charge to get F2. Hence, there is no high-twist, target-mass
or other effects. It gives the shape of the F2 distribution roughly similar to other two model but with smaller values which is expected.
Please treat it as a reference as it won't give precise F2 without these effects included.

