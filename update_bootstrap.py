from sys import argv
from write_isobars import update_bootstrap

path = argv[1]
if len(argv)>2:
	datadir = argv[2]
	update_bootstrap(path,datadir)
else:
	update_bootstrap(path,"./data/")




