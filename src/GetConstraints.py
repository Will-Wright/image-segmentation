import tkinter
import copy
from PIL import Image, ImageTk

def GetConstraints(im_width, im_height, im_RGB, \
                   must_link_list1=[], must_link_list2=[], \
                   is_demo_mode=False):

    im_RGB_mod = copy.deepcopy(im_RGB)

           
    # Prints image to get must_link and cannot_link inputs
    if (not must_link_list1) or (not must_link_list2):    
        window = tkinter.Tk(className='Initial image')
        frame = tkinter.Canvas(window, width=im_width, height=im_height)
        frame.focus_set()
        frame.pack()
        
        im_tk = ImageTk.PhotoImage(im_RGB)
        im_tk_frame = frame.create_image(im_width//2, im_height//2, image=im_tk)
        must_link_list1 = list()
        must_link_list2 = list()
        num_return_pressed = list()
        num_return_pressed.append(0)
        print("Click on pixels in region 1.  Press <Return> to change to region 2.")
        def key(event):
            if (event.keysym == 'Return') and (num_return_pressed[0] == 0):
                num_return_pressed[0] = 1
                print("Click on pixels in region 2.  Press <Return> when done.")
            elif (event.keysym == 'Return') and (num_return_pressed[0] == 1):
                window.destroy()            
        def callback(event):
            global im_tk
            dot_len = 5
            if (0 <= event.x <= im_width-1) and (0 <= event.y <= im_height-1):
                w_min = max(0, event.x-dot_len); w_max = min(im_width-1, event.x+dot_len);
                h_min = max(0, event.y-dot_len); h_max = min(im_height-1, event.y+dot_len);
                for w in range(w_min, w_max+1):
                    for h in range(h_min, h_max+1):
                        if (event.x - w)**2 + (event.y - h)**2 <= dot_len**2 - 5:
                            dot_pix_val = 255*num_return_pressed[0]
                            im_RGB_mod.putpixel((w,h), (dot_pix_val, dot_pix_val, dot_pix_val) )
                        elif (event.x - w)**2 + (event.y - h)**2 <= dot_len**2 + 4:
                            dot_pix_val = 255*(1-num_return_pressed[0])
                            im_RGB_mod.putpixel((w,h), (dot_pix_val, dot_pix_val, dot_pix_val) )
                im_tk = ImageTk.PhotoImage(im_RGB_mod)
                frame.itemconfigure(im_tk_frame, image=im_tk)
                if (num_return_pressed[0] == 0):
                    must_link_list1.append([event.x, event.y])
                    print("Clicked at (%i," %event.x, "%i)" %event.y)
                elif (num_return_pressed[0] == 1):
                    must_link_list2.append([event.x, event.y])
                    print("Clicked at (%i," %event.x, "%i)" %event.y)
        frame.bind("<Button-1>", callback)
        frame.bind("<Key>", key)
        tkinter.mainloop()
        
    elif is_demo_mode: 
        window = tkinter.Tk(className='Initial image')
        frame = tkinter.Canvas(window, width=im_width, height=im_height)
        frame.focus_set()
        frame.pack()
    
        num_clicked = list()
        num_clicked.append(0)
        im_tk = ImageTk.PhotoImage(im_RGB)
        im_tk_frame = frame.create_image(im_width//2, im_height//2, image=im_tk)
        print("Must-link pixels passed by user.\nClick on image to see pixels.")
        def callback(event):
            global im_tk
            if (num_clicked[0] == 0):
                dot_len = 5
                #global im_tk
                for i in range(0,len(must_link_list1)):
                    (w_val, h_val) = must_link_list1[i]
                    w_min = max(0, w_val-dot_len); w_max = min(im_width-1, w_val+dot_len);
                    h_min = max(0, h_val-dot_len); h_max = min(im_height-1, h_val+dot_len);        
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
                    h_min = max(0, h_val-dot_len); h_max = min(im_height-1, h_val+dot_len);        
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
                num_clicked[0] = 1
                print("Click on image again to continue.")
            elif (num_clicked[0] == 1):
                window.destroy()
        frame.bind("<Button-1>", callback)
        tkinter.mainloop()
      
    return (must_link_list1, must_link_list2)