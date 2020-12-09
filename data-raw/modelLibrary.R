### The nlmixrUI object contains several values that
### are random at every generation
### Saving an nlmixrUI object is therefore not reproducible,
### and results in changes in git at every regeneration.
###
### We explored many solutions to this problem. The optimal solution depends on
### your requirements:
### A) Do you want the models to be listed in the `data()` section?
### B) Do you want the models to be available through lazy-loading ?
### C) Are you willing to adapt tdmore code?
###
### Our aim is to allow a user to simply use the demo models. They should be attached lazily.
### It may also be interesting to have them in the `data()` section, but this is not a requirement IMHO.
### Therefore, the best solution is to store the R source code in inst/models/ and attach this
### through a lazy-load mechanism. The `data()` section is too limited to include lots of meta-data anyway.
###
### Solutions explored:
### 1) Store a serializable version of the object (e.g. nlmixr or RxODE),
### and change it back to the full version upon load.
### Unfortunately, this circumvents the `data()` system from R.
### 2) Adapt the nlmixrUI object to remove all random values.
### This is deemed impossible, because the random values are also in function environments.
### This is a wild goose chase. See e.g. x$nmodel$pred
### 3) Store the tdmore object, which will include an RxODE object.
### We can adapt the RxODE object to remove all random values (which are generated upon compiling).
### Unfortunately, adapting the RxODE object is also very difficult...
### 4) R supports loading data as .R source files.
### However, it sources these files without keep.source=TRUE. This is a requirement for nlmixr,
### so no dice yet again!!
###
