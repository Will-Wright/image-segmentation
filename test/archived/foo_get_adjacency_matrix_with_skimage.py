import Tkinter
from PIL import Image, ImageTk
#from skimage import data, segmentation
#from skimage.future import graph
import time
from sys import argv
import numpy as np

def main(*argv):
    window = Tkinter.Tk(className="Original image")

    image = Image.open(argv[0] if len(argv) >=1 else 'mellon_collie_and_the_infinite_sadness_test.jpg').convert('LA')
    image_width = image.size[0]
    image_height = image.size[1]
    num_pixels = image_height*image_width

    canvas = Tkinter.Canvas(window, width=image_width, height=image_height)
    canvas.pack()

    image_array = np.array(image).astype(np.float16)[0:image_height, 0:image_width, 0]
#    image_array_scaled = StandardScaler().fit_transform(image_array)

    image_tk = ImageTk.PhotoImage(image)
    canvas.create_image(image_width//2, image_height//2, image=image_tk)

#    labels = segmentation.slic(image_array)

#    timer_start = time.time()
#    rag = graph.rag_mean_color(image_array, labels)
#    timer_end = time.time()
#    runtime_skimage = timer_end - timer_start
#    print(runtime_skimage)

    def callback(event):
	print "clicked at: ", event.x, event.y

    canvas.bind("<Button-1>", callback)
    Tkinter.mainloop()

    return (image_array)






