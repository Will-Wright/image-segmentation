import tkinter
from PIL import Image, ImageTk
from sys import argv
import numpy as np

def main():
    window = tkinter.Tk(className="bla")

    image = Image.open(argv[1] if len(argv) >=2 else 'mellon_collie_and_the_infinite_sadness_test.jpg')
    canvas = tkinter.Canvas(window, width=image.size[0], height=image.size[1])
    canvas.pack()
    image_tk = ImageTk.PhotoImage(image)
    canvas.create_image(image.size[0]//2, image.size[1]//2, image=image_tk)

    def callback(event):
	print("clicked at: ", event.x, event.y)

    canvas.bind("<Button-1>", callback)
    tkinter.mainloop()

