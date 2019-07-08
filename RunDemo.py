import os.path
import numpy as np
from PIL import Image, ImageTk
import tkinter
import copy
from matplotlib import interactive
import matplotlib.pyplot as plt

from ImSeg import ImSeg


def main():

    demo_cache_file = 'demo_cache.npz'
    
    if os.path.exists(demo_cache_file):
        demo_cached = True
    else:
        demo_cached = False
 

    im_path = 'test/person_walking.jpg'
    must_link_list1 = [[273, 33], [279, 88], [235, 116], [265, 154], [318, 114], 
                       [281, 192], [264, 266], [275, 300]]
    must_link_list2 = [[222, 286], [164, 73], [127, 81], [36, 32], 
                       [40, 106], [34, 163], [227, 36], [378, 25], 
                       [394, 144], [391, 281], [217, 141], [335, 119]]
    
    resize_ratio = 0.90
    resize_min_num_pix = 200
    resize_filter = Image.ANTIALIAS
    
    
    im_RGB = Image.open(im_path)
    im_RGB_mod = copy.deepcopy(im_RGB)
    im_width = im_RGB.size[0]
    im_height = im_RGB.size[1]
    im_num_pix = im_width*im_height
    
    
    sig_dist=4 # default: 3-5
    sig_feat=7 # default: 5-7
    max_dist=3 # default: 3-5

    use_RGB = True
    is_demo_mode=False

    verbosity=0

    attempt_resize=False

    return_data = True

    eigsh_tol_init = 1e-4
    eigsh_tol_update = 5e-1
    eigsh_tol_min = 1e-4

    lobpcg_tol_init = 1e-4
    lobpcg_tol_update = 5e-1
    lobpcg_tol_min = 1e-4

    pr_tol=1e-5
    du_tol=1e-5
    gap_tol=1e-5

    #eig_solver = "eigsh"
    eig_solver = "lobpcg"
    SDPSS_max_iters = 50
    SDPSS_pr_update = 1

    use_NewtonEigSolver = True
    #use_NewtonEigSolver = False
    #use_Newton_step = True
    use_Newton_step = False
    Newton_max_iters = 50
    Newton_eig_solver="eigsh"
    #Newton_eig_solver="lobpcg"

    cvx_verbose = False
    
    
   
    window = tkinter.Tk(className='Initial image')
    frame = tkinter.Canvas(window, width=im_width, height=im_height)
    frame.focus_set()
    frame.pack()

    num_return_pressed = list()
    num_return_pressed.append(0)
    im_tk = ImageTk.PhotoImage(im_RGB)
    im_tk_frame = frame.create_image(im_width//2, im_height//2, image=im_tk)
    print("ImSeg demo called.  Displaying full-size demo image.")
    print("To view must-link constraints, hover mouse over image and press <Return>.")

    def key(event):
        global im_tk
        if (event.keysym == 'Return') and (num_return_pressed[0] == 0):
            dot_len = 5
            for i in range(0,len(must_link_list1)):
                (w_val, h_val) = must_link_list1[i]
                w_min = max(0, w_val-dot_len); w_max = min(im_width-1, w_val+dot_len);
                h_min = max(0, h_val-dot_len); h_max = min(im_width-1, h_val+dot_len);        
                for w in range(w_min, w_max+1):
                    for h in range(h_min, h_max+1):
                        if (w_val - w)**2 + (h_val - h)**2 <= dot_len**2 - 5:
                            dot_pix_val = 255
                            im_RGB_mod.putpixel((w,h), (dot_pix_val, dot_pix_val, dot_pix_val) )
                        elif (w_val - w)**2 + (h_val - h)**2 <= dot_len**2 + 4:
                            dot_pix_val = 0
                            im_RGB_mod.putpixel((w,h), (dot_pix_val, dot_pix_val, dot_pix_val) )
            for i in range(0,len(must_link_list2)):
                (w_val, h_val) = must_link_list2[i]
                w_min = max(0, w_val-dot_len); w_max = min(im_width-1, w_val+dot_len);
                h_min = max(0, h_val-dot_len); h_max = min(im_width-1, h_val+dot_len);        
                for w in range(w_min, w_max+1):
                    for h in range(h_min, h_max+1):
                        if (w_val - w)**2 + (h_val - h)**2 <= dot_len**2 - 5:
                            dot_pix_val = 0
                            im_RGB_mod.putpixel((w,h), (dot_pix_val, dot_pix_val, dot_pix_val) )
                        elif (w_val - w)**2 + (h_val - h)**2 <= dot_len**2 + 4:
                            dot_pix_val = 255
                            im_RGB_mod.putpixel((w,h), (dot_pix_val, dot_pix_val, dot_pix_val) ) 
            im_tk = ImageTk.PhotoImage(im_RGB_mod)
            frame.itemconfigure(im_tk_frame, image=im_tk)
            num_return_pressed[0] = 1
            print("To continue demo, hover over image and press <Return> again.")
        elif (event.keysym == 'Return') and (num_return_pressed[0] == 1):
            window.destroy()
    frame.bind("<Key>", key)
    tkinter.mainloop()
    

    
    resize_dims_list = list()
    resize_ratio_list = list()
        
    resize_dims_list.append((im_width, im_height))
    resize_ratio_list.append(1)

    im_num_pix_temp = resize_dims_list[-1][0]*resize_dims_list[-1][1]

    while (im_num_pix_temp >= resize_min_num_pix):
        num_resized_images = len(resize_ratio_list)
        im_size_resize = (int(round(  im_width*(resize_ratio**num_resized_images)  )), \
                          int(round( im_height*(resize_ratio**num_resized_images)  )))
        resize_dims_list.append(im_size_resize)
        resize_ratio_list.append(resize_ratio**num_resized_images)
        im_num_pix_temp = im_num_pix * (resize_ratio**(2*(num_resized_images+1)))        
    
    
    #print(resize_dims_list)
    #print(resize_ratio_list)

    
    num_exp_total = len(resize_ratio_list)
    SDPSS_time = np.zeros((num_exp_total,1))
    NES_time = np.zeros((num_exp_total,1))
    SDPSS_eig_calls = np.zeros((num_exp_total,1))
    NES_eig_calls = np.zeros((num_exp_total,1))
    num_pixels = np.zeros((num_exp_total,1))

    if not demo_cached:
        for num_exp in range(0, num_exp_total):
            im_size_resize = resize_dims_list[num_exp_total - (num_exp + 1)]
            resize_ratio = resize_ratio_list[num_exp_total - (num_exp + 1)]
            
            im_RGB_resize = Image.open(im_path)
            im_RGB_resize = im_RGB_resize.resize(im_size_resize, resize_filter)
            
            print("Running experiment %i" %(num_exp+1), "of %i" %num_exp_total\
                  , "with image size", im_size_resize)

            if (len(must_link_list1) + len(must_link_list2) >= 1):
                mll_resize1 = [[int(resize_ratio*i) for i in tups] \
                             for tups in must_link_list1]
                mll_resize2 = [[int(resize_ratio*i) for i in tups] \
                             for tups in must_link_list2]
            else:
                mll_resize1 = []
                mll_resize2 = []

            #print(mll_resize1)
            #print(mll_resize2)

            [im_RGB_arr, adj_mat, mll_resize1, mll_resize2, \
             X_sub, V_sub, im_seg_RGB_arr, im_RGB, N, data] \
                = ImSeg(im_path = [], Tk_im_file = im_RGB_resize, \
                        attempt_resize=attempt_resize, \
                        return_data=return_data, \
                        cvx_verbose=cvx_verbose, eig_solver = eig_solver, \
                        pr_tol=pr_tol, du_tol=du_tol, gap_tol=gap_tol, \
                        sig_dist=sig_dist, sig_feat=sig_feat, max_dist=max_dist, \
                        eigsh_tol_init=eigsh_tol_init, \
                        eigsh_tol_update=eigsh_tol_update, eigsh_tol_min=eigsh_tol_min, \
                        lobpcg_tol_init=lobpcg_tol_init, \
                        lobpcg_tol_update=lobpcg_tol_update, \
                        lobpcg_tol_min=lobpcg_tol_min, \
                        SDPSS_max_iters=SDPSS_max_iters, \
                        SDPSS_pr_update=SDPSS_pr_update, \
                        must_link_list1 = mll_resize1, \
                        must_link_list2 = mll_resize2, \
                        use_RGB=use_RGB, \
                        use_NewtonEigSolver=use_NewtonEigSolver, \
                        Newton_eig_solver=Newton_eig_solver, \
                        use_Newton_step=use_Newton_step, \
                        Newton_max_iters=Newton_max_iters, \
                        is_demo_mode=is_demo_mode, verbosity=verbosity);
            SDPSS_time[num_exp] = data['SDPSS_data']['time_total']
            NES_time[num_exp] = data['NES_data']['time_total']
            SDPSS_eig_calls[num_exp] = data['SDPSS_data']['num_eig_calls']
            NES_eig_calls[num_exp] = data['NES_data']['num_eig_calls']
            num_pixels[num_exp] = im_size_resize[0]*im_size_resize[1]

        
        # Saves figure of performance comparison 
        fig = plt.figure()
        fig.canvas.draw()
        ax = fig.add_subplot(111)
        plt.plot(num_pixels, NES_time, 'r--', label='Newton')
        plt.plot(num_pixels, SDPSS_time, 'b', label='SDP-Subspace')

        fig.suptitle('SDPSS vs Newton for various sized images')
        plt.xlabel('Number pixels (x1000)')
        plt.ylabel('Runtime (secs)')
        plt.legend()
        
        labels=ax.get_xticks().tolist()
        labels = [np.int16(l/1000) for l in labels]
        ax.set_xticklabels(labels)
        plt.savefig('demo_fig1.png')
        
        
        np.savez_compressed(demo_cache_file, SDPSS_time=SDPSS_time, \
                            NES_time=NES_time, SDPSS_eig_calls=SDPSS_eig_calls,\
                            NES_eig_calls=NES_eig_calls, num_pixels=num_pixels,\
                            im_seg_RGB_arr=im_seg_RGB_arr)
        cache_data = np.load(demo_cache_file)

    else:
        cache_data = np.load(demo_cache_file)
        SDPSS_time = cache_data['SDPSS_time']
        NES_time = cache_data['NES_time']
        num_pixels = cache_data['num_pixels']
        SDPSS_eig_calls = cache_data['SDPSS_eig_calls']
        NES_eig_calls = cache_data['NES_eig_calls']
        im_seg_RGB_arr = cache_data['im_seg_RGB_arr']
        
    print("All experiments complete.\nPlotting segmented image.")
        
    im_seg_RGB = Image.fromarray(im_seg_RGB_arr)
    im_seg_RGB.show(title="Newton Image")
    
    return cache_data

