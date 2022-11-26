This is the code of the diploma thesis of my undergraduate degree in Mechanical Engineering. 

You can find the full text here (incl. ENG abstract): https://dspace.lib.ntua.gr/xmlui/bitstream/handle/123456789/51815/Vavoula.pdf?sequence=4&isAllowed=y

The objective is to create a population–based stochastic optimization method that is based on the concept of the optimization method Harmony Search. 
New operators are applied for the improvement of the classical Harmony Search method, based on the theory of dance composition, and as a result, the method is named Choreography Composition-Like
Optimization (CCLO). 

The proposed method uses machine learning and can solve multi-objective problems with or without constraints by computing Pareto fronts.

For its assessment, the method is applied to:
> Several computational minimization problems 
> An aerodynamics problem, where the shape of an isolated airfoil is optimized, either for maximum lift coefficient or for both maximum lift and minimum drag coefficients. The shape of the airfoil is also examined under constraints on the moment coefficient.
> A finance problem, where the economic order quantity of a product in a company is calculated. The company uses the shortage model and the supplier applies tiered pricing to the
product.


Notes:
> After compiling the code, it should be called as: CLO.exe .\eas.cfg, where eas.cfg is the file with the initial values of the variables.
> The files eas.cfg and ann.conf, which are required, are in the TestFunctions folder, for each Test Function. 
> In the same folder you can also find the objective.exe for the calculation of the objective functions.
> You can select to run the CCLO code without machine learning, through the initialization of the eas.cfg file.