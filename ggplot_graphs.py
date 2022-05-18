##import subprocess
##
##res = subprocess.call("Rscript /Users/dradecic/Desktop/script.R", shell=True)
##res

import rpy2

##from rpy2 import robjects

pi = rpy2.robjects.r['pi']
pi


##from rpy2.robjects.packages import importr, data
##import rpy2.robjects.lib.ggplot2 as ggplot2
##
##grdevices = importr('grDevices')
##grdevices.png(file="/Users/dradecic/Desktop/mtcars.png", width=1024, height=512)
##datasets = importr('datasets')
##mtcars = data(datasets).fetch('mtcars')['mtcars']
##
##pp = (ggplot2.ggplot(mtcars) +
##      ggplot2.aes_string(x='wt', y='mpg', col='factor(cyl)') +
##      ggplot2.geom_point())
##pp.plot()
##
##grdevices.dev_off()
