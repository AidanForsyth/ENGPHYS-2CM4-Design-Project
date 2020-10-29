TITLE 'PPMT strap optimization'    
COORDINATES cartesian2  { coordinate system, 1D,2D,3D, etc }
VARIABLES        { system variables }
    Voltage              
	u
	v
    rho_free
SELECT         { method controls }
errlim = 1e-6

DEFINITIONS    { parameter definitions }

Lx = 10e-6 !length of copper PZT sandwhich (left or right side)
Lxapp = 10e-6 !diameter of the apperature hole
LyCU = 2e-6 !height of copper layer
LyPZT = 3e-6 !height of PZT layer
LyStrap = 5.999999999999999e-06 

!omega = 20*2*pi {not the actual resonant frequency of the PZT material  just was a test case for a time dependant simulation}

rhoe 
rho

Efield = -grad(Voltage)
Efieldx = xcomp(Efield)
Efieldy = ycomp(Efield)
Efieldz =0 !no z-component being modelled 

J = 1/rhoe*Efield

Vmax = 1 !assuming that the PPMT is run at a steady state/low frequency


E = 110e9 !parameter for CU
nu = 0.34 !parameter for CU

G = E/(2*(1+nu))

C11 =E*(1-nu)/((1+nu)*(1-2*nu)) C12 = E*nu/((1+nu)*(1-2*nu))	C13 = C12		C14 = 0			C15 = 0			C16 = 0
C21 = C12	C22 = C11				C23 = C13			C24 = 0			C25 = 0 C26 = 0
C31 = C13	C32 = C23		C33 =C11	C34 = 0		C35 = 0	C36 = 0
C41 = C14	C42 = C24		C43 = C34	C44 = G		C45 = 0	C46 =0
C51 = C15	C52 = C25		C53 = C35	C54 = C45	C55 = G		C56 = 0 
C61 = C16	C62 = C26		C63 = C36	C64 = C46	C65 = C56		C66 = G


DeltaTemp = 0 !Thermal Expansion is turned off

alpha = 1e-6

!thermal expansion coefficients
alphax = alpha
alphay = alpha
alphaz = 0
alphayz = 0 !0 unless monoclinic or triclinic
alphaxz = 0
alphaxy = 0

!piezoelectric coupling coefficients
d11 = 0				d12 = 0		d13 = 0			d14 = 0		d15 = 584e-12	d16 = 0
d21 = 0				d22 = 0		d23 = 0			d24 = d15	d25 = 0			d26 = 0
d31 = -171e-12	d32 = d31	d33 = 374e-12	d34 = 0		d35 = 0			d36 = 0

!Strain definitions from displacements
ex = dx(u)
ey = dy(v)
ez = 0!dz(w)
gyz = 0!dy(w) + dz(v)
gxz = 0!dx(w) + dz(u)
gxy = dx(v) + dy(u)

!Mechanical strain
exm  = ex  - alphax *DeltaTemp - (d11*Efieldx+d21*Efieldy+d31*Efieldz)
eym  = ey  - alphay *DeltaTemp - (d12*Efieldx+d22*Efieldy+d32*Efieldz)
ezm  = ez  - alphaz *DeltaTemp - (d13*Efieldx+d23*Efieldy+d33*Efieldz)
gyzm = gyz - alphayz*DeltaTemp - (d14*Efieldx+d24*Efieldy+d34*Efieldz)
gxzm = gxz - alphaxz*DeltaTemp - (d15*Efieldx+d25*Efieldy+d35*Efieldz)
gxym = gxy - alphaxy*DeltaTemp - (d16*Efieldx+d26*Efieldy+d36*Efieldz)

!Hookes Law
sx  = C11*exm + C12*eym + C13*ezm + C14*gyzm + C15*gxzm + C16*gxym
sy  = C21*exm + C22*eym + C23*ezm + C24*gyzm + C25*gxzm + C26*gxym
sz  = C31*exm + C32*eym + C33*ezm + C34*gyzm + C35*gxzm + C36*gxym
syz = C41*exm + C42*eym + C43*ezm + C44*gyzm + C45*gxzm + C46*gxym
sxz = C51*exm + C52*eym + C53*ezm + C54*gyzm + C55*gxzm + C56*gxym
sxy = C61*exm + C62*eym + C63*ezm + C64*gyzm + C65*gxzm + C66*gxym

epsilonrel = 0

epsilon0 = 8.854e-12
epsilonT = epsilonrel*epsilon0

DfieldPiezox = d11*sx+d12*sy+d13*sz+d14*syz+d15*sxz+d16*sxy
DfieldPiezoy = d21*sx+d22*sy+d23*sz+d24*syz+d25*sxz+d26*sxy
DfieldPiezoz = d31*sx+d32*sy+d33*sz+d34*syz+d35*sxz+d36*sxy

DfieldPiezo = vector(DfieldPiezox, DfieldPiezoy, DfieldPiezoz)
Dfield = epsilonT*Efield + DfieldPiezo

{ scaling factor for displacement plots } 
   Mt =0.1*globalmax(magnitude(x,y))/globalmax(magnitude(U,V)) 

EQUATIONS        { PDE's, one for each variable }
Voltage:  	div(J)=0 
u: dx(sx) + dy(sxy) = 0
v: dx(sxy) + dy(sy) = 0
rho_free: div(Dfield) = rho_free

BOUNDARIES       

REGION 'Aluminum Strap'

E = 68.3e9
nu = 0.33
rho = 2700
rhoe = 36e-9

 {Not a piezoelectric, therefore, no d-matrix}
 d11 =0		d12 = 0	d13 = 0	d14 = 0	d15 = 0	d16 = 0
d21 = 0	d22 = 0	d23 = 0	d24 = 0	d25 = 0	d26 =0
d31 = 0	d32 = 0	d33 = 0	d34 = 0	d35 = 0	d36 = 0

C11 =E*(1-nu)/((1+nu)*(1-2*nu)) C12 = E*nu/((1+nu)*(1-2*nu))	C13 = C12		C14 = 0			C15 = 0			C16 = 0
C21 = C12	C22 = C11				C23 = C13			C24 = 0			C25 = 0 C26 = 0
C31 = C13	C32 = C23		C33 =C11	C34 = 0		C35 = 0	C36 = 0
C41 = C14	C42 = C24		C43 = C34	C44 = G		C45 = 0	C46 =0
C51 = C15	C52 = C25		C53 = C35	C54 = C45	C55 = G		C56 = 0 
C61 = C16	C62 = C26		C63 = C36	C64 = C46	C65 = C56		C66 = G

START (0, LyCU*2 + LyPZT)
LINE TO (Lx*2 + Lxapp, LyCU*2 + LyPZT)
LINE TO (Lx*2 + Lxapp, LyCU*2 + LyPZT + LyStrap)
LINE TO (0, LyCU*2 + LyPZT + LyStrap)
value(u) = 0
value(v) = 0
LINE TO CLOSE
 
 REGION 'Copper Sandwich'
 
 rho = 8960
 rhoe = 16.78e-9
 
 {Not a piezoelectric, therefore, no d-matrix}
 d11 =0		d12 = 0	d13 = 0	d14 = 0	d15 = 0	d16 = 0
d21 = 0	d22 = 0	d23 = 0	d24 = 0	d25 = 0	d26 =0
d31 = 0	d32 = 0	d33 = 0	d34 = 0	d35 = 0	d36 = 0

 !Bottom
 START (0,0)
  load(u) = 0
 load(v) = 0
 load(Voltage) = 0
 LINE TO (Lx*2 + Lxapp, 0)
 LINE TO (Lx*2 + Lxapp, LyCU)
 LINE TO (0, LyCU)
 value(u) = 0
 value(v) = 0
value(Voltage) = 0
 LINE TO CLOSE
 !Upper
 START (0, LyCU + LyPZT)
 load(u) = 0
 load(v) = 0
 load(Voltage) = 0
 LINE TO (Lx*2 + Lxapp, LyCU + LyPZT)
 LINE TO (Lx*2 + Lxapp, LyCU*2 + LyPZT)
 LINE TO (0, LyCU*2 + LyPZT)
 value(Voltage) = Vmax !*sin(omega*t)
 value(u) = 0
 value(v) = 0
 LINE TO CLOSE
 
 REGION 'PZT Layer' 
 
 rho = 0
 rhoe = 1 !this represents a very high rhoe value, therefore, inhibiting the current flow throught the material 
 epsilonrel = 1730
 
C11 =82e9	C12 = 35e9		C13 = 34e9	C14 = 0		C15 = -2e9	C16 = 0
C21 = C12	C22 = 63e9		C23 = 35e9	C24 = 0		C25 = -8e9	C26 = 0
C31 = C13	C32 = C23		C33 = 58e9	C34 = 0		C35 = -3e9	C36 = 0
C41 = C14	C42 = C24		C43 = C34	C44 = 21e9	C45 = 0		C46 = -5e9
C51 = C15	C52 = C25		C53 = C35	C54 = C45	C55 = 28e9	C56 = 0
C61 = C16	C62 = C26		C63 = C36	C64 = C46	C65 = C56	C66 = 29e9
 
 START (0, LyCU)
 load(u) = 0
 load(v) = 0
 load(Voltage) = 0
 LINE TO (Lx*2 + Lxapp, LyCU)
 LINE TO (Lx*2 + Lxapp, LyCU + LyPZT) 
 LINE TO (0, LyCU + LyPZT)
 value(u) = 0
 value(v) = 0
 LINE TO CLOSE 

!TIME 0 TO 3    { if time dependent }

MONITORS         { show progress }

PLOTS      
	!for t = endtime     
 grid(x+Mt*u, y+Mt*v) 
 	report val(v, Lx*2 + Lxapp, (LyCU*2 + LyPZT + LyStrap)/2)
  !CONTOUR(Voltage) painted
	!vector(J) norm
	!CONTOUR(rho_free) painted
	!vector(Dfield) norm
	!vector(Efield) norm

elevation(v) from (10e-6,2e-6) to (10e-6,3e-6)  !export format '#x#b#1' file='2CM4ProjectForsya5.txt' 

SUMMARY
SUMMARY EXPORT FILE =  "2CM4ProjectForsya5.txt" !a way to export the tip displacement, as exporting data from the elevation plot didnt work
    report val(v, Lx*2 + Lxapp, (LyCU*2 + LyPZT + LyStrap)/2) as "Tip Displacement "
	

	!report val(DFieldPiezox, Lx/2, hb+hp/2) as "DFieldX"
	!report val(DFieldPiezoy, Lx/2, hb+hp/2) as "DFieldY"

END

