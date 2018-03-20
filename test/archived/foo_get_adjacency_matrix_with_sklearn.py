import Tkinter
from PIL import Image, ImageTk
from sklearn.metrics import pairwise
from sklearn.preprocessing import StandardScaler
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
    image_array_scaled = StandardScaler().fit_transform(image_array)

    image_tk = ImageTk.PhotoImage(image)
    canvas.create_image(image_width//2, image_height//2, image=image_tk)

    image_vector = np.reshape(image_array_scaled, [num_pixels, 1])
    affinity_type = 'rbf' #'rbf'
    params = {}
    params['gamma'] = 1.0
    params['degree'] = 3
    params['coef0'] = 1

    timer_start = time.time()
    affinity_matrix = pairwise.pairwise_kernels(image_vector, metric=affinity_type, filter_params=True, **params)
    timer_end = time.time()
    runtime_sklearn = timer_end - timer_start
    print(runtime_sklearn)

    def callback(event):
	print "clicked at: ", event.x, event.y

    canvas.bind("<Button-1>", callback)
    Tkinter.mainloop()

    return (image_array, image_array_scaled, affinity_matrix)






