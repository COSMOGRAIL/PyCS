
PyCS_D3CScombi_v1.dt
====================

- pour les courbes avec bcp de points et du bon signal (soit rung 0 et 3) :
Les mesures à l'oeil sont des mesures bruitées quasi-indépendantes car on regarde tous des bouts différents.
Donc erreur sur la médiane = deviation standard / sqrt(N)

- pour les autres courbes (saisons, faible signal, microlensing) :
Pas de réduction de l'erreur en sqrt(N).
On prend l'erreur maximale entre
       - dev std des délais
       - médiane des erreurs individuelles de D3CS


bottomline :
"scatter" for rungs 0 and 3
"d3cscombi1" for the other rungs




PyCS_regdiff2_v1.dt
===================

pycs.tdc.run.multirun(iniests, 
	sploptfct = pycs.tdc.splopt.spl2, optfct = pycs.tdc.optfct.regdiff2,
	addmlfct = None,
	nobs = 10, nsim = 20,
	ncpu=0,
	diagnostics = True,
	saveplots = True,
	tdcpath = "./tdc0", outdir="./test2_regdiff2", method="regdiff2", methodpar="First test")


PyCS_spl2_v1.dt
===============

pycs.tdc.run.multirun(iniests, 
	sploptfct = pycs.tdc.splopt.spl2, optfct = pycs.tdc.splopt.spl2,
	addmlfct = None but using pycs.tdc.splopt.splml1 for curve 6 6 !
	nobs = 10, nsim = 20,
	ncpu=0,
	diagnostics = True,
	saveplots = True,
	tdcpath = "./tdc0", outdir="./test2_spl2", method="spl2", methodpar="First test")



