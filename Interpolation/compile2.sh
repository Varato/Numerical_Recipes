gcc -fPIC -shared polint_mod.c ../nrlib/nrutil.c ../nrlib/polint.c -o polint_mod.so -I/usr/include/python2.7 -lpython2.7
