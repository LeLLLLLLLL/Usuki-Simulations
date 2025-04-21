import numpy as np
import math
import cmath
import time
import os

if os.path.exists("waves.txt"): os.remove("waves.txt")
if os.path.exists("tr_b.txt"): os.remove("tr_b.txt")

# Constants
rows = 119
cols = 68
islmax = 300
a = 2.50e-9  # Grid size
b = 0.000    # Magnetic field in Tesla
EF0 = 0.014338  # Fermi energy in eV
alpha = 68.214  # Spreading factor in nm/V

# Derived Constants
m0 = 9.10938356e-31  # Electron mass (kg)
hbar = 1.0545718e-34  # Reduced Planck constant (J·s)
q = 1.602176634e-19   # Elementary charge (C)
mass = 0.45 * m0
rmass = .067
angfac=0.2626*rmass
hb2o2m = (hbar / (2 * mass)) * (hbar / q)
ci = 1j
c1 = complex(1.0, 0.0)
c0 = complex(0.0, 0.0)

#thop = hb2o2m / (a ** 2)
xmax=4000.0
ymin=0.0
ymax=7000.0
dely=(ymax-ymin)/(rows+1)
delx=dely
thop =1.0/(delx*delx*angfac)
nsl=int(xmax/delx)
# nsl should be used instead of cols

emax=0.00305
emin=0.00
ne=70
de=(emax-emin)/(ne-1)

bet = q * b * (a ** 2) / hbar

# Matrices and arrays
#pot = np.zeros((rows, nsl), dtype=np.float64)
vel = np.zeros(rows)
psil = np.zeros(rows)
pot = np.full((nsl+1,rows), 3)
evec = np.zeros((rows, rows), dtype=np.complex128)
up = np.zeros((rows, rows), dtype=np.complex128)
upl = np.zeros((rows, rows), dtype=np.complex128)
um = np.zeros((rows, rows), dtype=np.complex128)
uml = np.zeros((rows, rows), dtype=np.complex128)
c1l1 = np.zeros((rows, rows), dtype=np.complex128)
c2l1 = np.zeros((rows, rows), dtype=np.complex128)
d1l1 = np.zeros((rows, rows), dtype=np.complex128)
d2l1 = np.zeros((rows, rows), dtype=np.complex128)
p2i = np.zeros((rows, rows), dtype=np.complex128)
evs = np.zeros(rows, dtype=np.complex128)
T21 = np.zeros((rows, rows), dtype=np.complex128)
T22 = np.zeros((rows, rows), dtype=np.complex128)
Tl = np.zeros((2*rows, 2*rows), dtype=np.complex128)  # Initialize Tl with c0
pl1 = np.zeros((rows, rows, islmax), dtype=np.complex128)
pl2 = np.zeros((rows, rows, islmax), dtype=np.complex128)
pl2i = np.zeros((rows, rows, islmax), dtype=np.complex128)
psimode = np.zeros((islmax, rows, rows), dtype=np.complex128)
psipm = np.zeros((islmax, rows), dtype=np.double)
phi1new = np.zeros((rows, rows), dtype=np.complex128)
phi2new = np.zeros((rows, rows), dtype=np.complex128)
phi1old = np.zeros((rows, rows), dtype=np.complex128)
phi2old = np.zeros((rows, rows), dtype=np.complex128)

for i in range((rows//2)-15, (rows//2)+15):
    for j in range(nsl+1):
        pot[j, i] = 0.0
        pot[j, rows - i - 1] = 0.0

# energy loop
trans = 0.0
#for ef in range(30,31):
for ef in range(1,31):

    en=emin+ef*de
    ehop=en/thop

    # Propagation through columns
    for i in range(rows):
        # Note: Fortran is 1-based; Python is 0-based
        Pmi = -np.exp(ci * bet)
        Pmi1 = -np.exp(ci * bet)
        # Zero the i-th row
        T21[i, :] = c0
        T22[i, :] = c0
        # Diagonal entry of T21
        T21[i, i] = -Pmi * Pmi1
        # Diagonal entry of T22 with complex potential
        T22[i, i] = (ehop - 4.0 - (pot[1, i] )/thop) * Pmi
        # Off-diagonal entries of T22
        if i < rows - 1:
            T22[i, i + 1] = -Pmi
        if i > 0:
            T22[i, i - 1] = -Pmi

    for i in range(rows):
        Tl[i, i + rows] = c1
        Tl[i + rows, i] = T21[i, i]
        for j in range(rows):
            Tl[i + rows, j + rows] = T22[i, j]

    evals, evec = np.linalg.eig(Tl)
# Preallocate arrays
    xp = np.zeros(rows)
    xpe = np.zeros(rows)
    xm = np.zeros(rows)
    xme = np.zeros(rows)

    ipv = np.zeros(rows, dtype=int)
    ipev = np.zeros(rows, dtype=int)
    imv = np.zeros(rows, dtype=int)
    imev = np.zeros(rows, dtype=int)

    rnorm = np.zeros(2*rows)
    cur = np.zeros(2*rows)

    im = ip = ime = ipe = 0
    nprop2 = 0

# Compute rnorm
    for j in range(2*rows):
        rnorm[j] = sum(abs(evec[i, j])**2 for i in range(rows))

# Loop through all eigenvalues
    for j in range(2*rows):
        x = abs(evals[j])

        if 0.9999 < x < 1.0001:
            nprop2 += 1
            rk = (cmath.log(evals[j] ) / complex(0.0, 1.0)).real
            rnorm[j] = 0.0
            cur[j] = 0.0

            for i in range(rows):
                add = abs(evec[i, j])
                cur[j] += math.sin(rk + 2.0 * math.pi * bet * (i + 1)) * add ** 2
                rnorm[j] += add ** 2

            vf = cur[j] / rnorm[j] if rnorm[j] != 0 else 0.0

            if vf > 0.0:
                ipv[ip] = j
                xp[ip] = rk
                ip += 1
            else:
                imv[im] = j
                xm[im] = rk
                im += 1

        elif x < 0.9999:
            ipev[ipe] = j
            xpe[ipe] = x
            ipe += 1
            rnorm[j] = sum(abs(evec[i, j])**2 for i in range(rows))

        elif x > 1.0001:
            imev[ime] = j
            xme[ime] = x
            ime += 1
            rnorm[j] = sum(abs(evec[i, j])**2 for i in range(rows))

# Final value
    nprop = nprop2 // 2
    print("ip,ipe,im,ime",ip,ipe,im,ime)

# Sort xp/ipv ascending
    if ip > 1:
        sorted_indices = np.argsort(xp[:ip])
        xp[:ip] = xp[sorted_indices]
        ipv[:ip] = ipv[sorted_indices]

# Sort xm/imv descending
    if im > 1:
        sorted_indices = np.argsort(-xm[:im])
        xm[:im] = xm[sorted_indices]
        imv[:im] = imv[sorted_indices]

    for j in range(ip):
        idx = ipv[j]
        vel[j] = cur[idx] / rnorm[idx]

        rkval = xp[j] / dely
        evalp = ehop * thop - (rkval) ** 2.0 / angfac
        abs_eval_diff = abs(evals[idx]) - 1.0

        for i in range(rows):
            up[i, j] = evec[i, idx] / np.sqrt(delx * delx * rnorm[idx])
            upl[i, j] = evec[i + rows, idx] / np.sqrt(delx * delx * rnorm[idx])

    # Evanescent positive modes
    for j in range(ipe):
        idx = ipev[j]
        for i in range(rows):
            up[i, j + ip] = evec[i, idx] / np.sqrt(delx * delx * rnorm[idx])
            upl[i, j + ip] = evec[i + rows, idx] / np.sqrt(delx * delx * rnorm[idx])

    # Negative propagating modes
    for j in range(im):
        idx = imv[j]
        for i in range(rows):
            um[i, j] = evec[i, idx] / np.sqrt(delx * delx * rnorm[idx])
            uml[i, j] = evec[i + rows, idx] / np.sqrt(delx * delx * rnorm[idx])

    # Evanescent negative modes
    for j in range(ime):
        idx = imev[j]
        for i in range(rows):
            um[i, j + im] = evec[i, idx] / np.sqrt(delx * delx * rnorm[idx])
            uml[i, j + im] = evec[i + rows, idx] / np.sqrt(delx * delx * rnorm[idx])

    print(en,nprop)
    #print(um(:,range(im)))
    if nprop >= 1:
        d2l1 = np.linalg.inv(uml)
        c2l1 = um @ d2l1
        p2i = np.linalg.inv(c2l1)
        d1l1 = -d2l1 @ upl
        c1l1 = up - c2l1 @ upl
        iii = 0
        for ii in range(1, nsl + 1):  
            istart = ii
            if istart >= 0:
                iii += 1

    # Copy matrices
            c1l = c1l1.copy()
            c2l = c2l1.copy()
            d1l = d1l1.copy()
            d2l = d2l1.copy()
    
            if istart >= 0:
                pl1[:, :, iii - 1] = c1l1
                pl2[:, :, iii - 1] = c2l1
                pl2i[:, :, iii - 1] = p2i
    
            for i in range(rows):
                Pmi = -np.exp(ci * bet)
                Pmi1 = -np.exp(ci * bet)
                T21[i, :] = c0
                T22[i, :] = c0
                T21[i, i] = -Pmi * Pmi1
                T22[i, i] = (ehop - 4.0 - (pot[ii, i] )) * Pmi
                if i < rows - 1:
                    T22[i, i + 1] = -Pmi
                if i > 0:
                    T22[i, i - 1] = -Pmi
            save = np.diag(T21)[:, np.newaxis] * c2l  
            p2i = save + T22
            c2l1 = np.linalg.inv(p2i)
            save = np.diag(T21)[:, np.newaxis] * c1l 
            save2 = c2l1 @ save                     
            c1l1 = -save2                          
            d2l1 = d2l @ c2l1
            save = d2l @ c1l1
            d1l1 = d1l + save

        c1l = c1l1.copy()
        c2l = c2l1.copy()
        d1l = d1l1.copy()
        d2l = d2l1.copy()
        pl1[:, :, nsl + 1] = c1l1
        pl2[:, :, nsl + 1] = c2l1
        pl2i[:, :, nsl + 1] = p2i
        # final slice
        upli = np.linalg.inv(upl)
        save = up @ upli
        save2 = c2l - save
        p2 = np.linalg.inv(save2)
        save2 = p2 @ c1l
        p1 = -save2
        save = upli @ save2
        c1l1 = -save
        save = d2l @ save2
        d1l1 = d1l - save
        pl1[:, :, nsl + 2] = p1
        pl2[:, :, nsl + 2] = p2
        pl2i[:, :, nsl + 2] = p2i
    
        trans = 0.0
        ref = 0.0
    
        for j in range(nprop):
            for i in range(nprop):
                term = (vel[i] / vel[j]) * abs(c1l1[i, j])**2
                trans += term
                ref += (vel[i] / vel[j]) * abs(d1l1[i, j])**2
    
        corrinel = (float(nprop) - (trans + ref)) / 2.0
    
        print('e,bmag,trans,ref,error')
        print(en, bet, trans, ref,float(nprop) - trans - ref)
    
        phi1new[:, :] = c0
        phi2new[:, :] = c0
        np.fill_diagonal(phi2new, c1)
        ii = nsl + 1
        phi1old[:, :] = phi1new
        phi2old[:, :] = phi2new
        p1[:, :] = pl1[:, :, ii + 1]
        p2[:, :] = pl2[:, :, ii + 1]
        phi2new = phi2old @ p2
        save = phi2old @ p1
        phi1new = phi1old + save
    
        nslice = 0
    
        # Backward propagation
        for lplot in range(nsl + 1, 0, -1):
            nslice += 1
            p2[:, :] = pl2[:, :, lplot]
            p1[:, :] = pl1[:, :, lplot]
            phi1old[:, :] = phi1new
    
            save = p2 @ phi1old
            phi1new = p1 + save
            phi1old[:, :] = phi1new
    
            for m in range(rows):
                for i in range(nprop):
                    psimode[lplot, m, i] = phi1new[m, i]
    
        # Accumulate probabilities
        for lplot in range(1, nsl + 2):  # Fortran 1 to nsl+1
            for m in range(rows):
                psil[m] = 0.0
                psipm[lplot, m] = 0.0
                for i in range(nprop):
                    amp2 = abs(psimode[lplot, m, i])**2
                    psil[m] += amp2
                    psipm[lplot, m] += amp2

    ttot = trans

    with open("tr_b.txt", "a") as trb:
        trb.write(f"{en:.8e} {bet:.8e} {ttot:.8e}\n")

    with open("waves.txt", "a") as wvs:
        for lplot in range(1, nsl + 2):  
            for m in range(rows):
                if psipm[lplot, m] < 1e-10:
                    psipm[lplot, m] = 1e-10

                wvs.write(f"{psipm[lplot, m]:.8e} \n")
