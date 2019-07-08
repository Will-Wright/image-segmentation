import numpy as np

def PrintIteration(iter, pr_obj, du_obj, pr_res, du_res, \
                   time_eig, time_SDP, time_main, pr_dim, tol_print, \
                   y, verbosity=1):
    if verbosity>=1:
        time_total = time_eig + time_SDP + time_main
        if verbosity==2:
            print("%3i" % iter, \
                  "|% 7.5e"%pr_obj, \
                  "% 5.2e" % du_res, \
                  "%4.2e" % pr_res, \
                  "%4.2e |" % np.absolute(pr_obj-du_obj), \
                  "%4.2f"%(time_eig/time_total), \
                  "%4.2f"%(time_SDP/time_total), \
                  "%4.2f"%(time_main/time_total), \
                  "| %3.1e" %tol_print, \
                  "%5i"%pr_dim, \
                  "%6.2f"%time_total, \
                  "| % 3.5e"% y[0], \
                  "% 3.5e" % y[1] ) 
        else:
            print("%3i" % iter, \
              "|% 7.5e"%pr_obj, \
              "% 5.2e" % du_res, \
              "%4.2e" % pr_res, \
              "%4.2e |" % np.absolute(pr_obj-du_obj), \
              "%4.2f"%(time_eig/time_total), \
              "%4.2f"%(time_SDP/time_total), \
              "%4.2f"%(time_main/time_total), \
              "| %3.1e" %tol_print, \
              "%5i"%pr_dim, \
              "%6.2f"%time_total)
    return