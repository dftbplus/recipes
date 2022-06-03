import os

#set here broadening parameter
gamma = 0.2
#set here spin contamination threshold
spin_cont = 0.5

#constants
pi = 3.141592654
C = 0.000035


exc_dat = os.path.join('EXC.DAT')
output = open(exc_dat, "r")
data = output.readlines()
output.close()

closed_shell = False
if 'Sym.' in data[1]:
    closed_shell = True

ener = []
osc = []
for line in data[5:]:
    line_split = line.split()
    if (closed_shell):
        ener.append(float(line_split[0]))
        osc.append(float(line_split[1]))
    else:
        if (abs(float(line_split[7])) < spin_cont):
            ener.append(float(line_split[0]))
            osc.append(float(line_split[1]))

wstart = int(ener[0])
wend = int(ener[-1]) + 1
deltaw = 0.01
npoint = int((wend - wstart)/deltaw)

w = wstart
for j in range(1,npoint):
    abs = 0.0
    for i in range(len(ener)):
        abs = abs + 0.5*osc[i]*gamma/(C*pi*((w-ener[i])**2 + 0.25*gamma**2))

    print("{:.2f}".format(w), "{:.2f}".format(abs))
    w = w + deltaw
