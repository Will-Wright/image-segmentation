from cv2 import *
import numpy as np

def main():
    img_array = cv2.imread('mellon_collie_and_the_infinite_sadness_test.jpg')
    cv2.namedWindow('Test Image', cv2.WINDOW_NORMAL)


    cv2.startWindowThread()


    cv2.imshow('Test Image', img_array)
    res = cv2.waitKey(1000)
#    print 'You pressed %d (0x%x), LSB: %d (%s)' % (res, res, res % 256,
#        repr(chr(res%256)) if res%256 < 128 else '?')

    cv2.waitKey(1)
#    cv2.destroyWindow('Test Image')
    cv2.destroyAllWindows()
    cv2.waitKey(0)
    cv2.waitKey(1)
    cv2.waitKey(1)
    cv2.waitKey(1)
    cv2.destroyAllWindows()





#    cv2.namedWindow('image',cv2.WINDOW_NORMAL)
#    cv2.resizeWindow('image', self.w, self.h)
#    cv2.setMouseCallback('image',self.mouse_callback)
#    while(not self.finished):
#        cv2.imshow('image',self.img)
#        k = cv2.waitKey(4) & 0xFF
#        if k == 27:
#             breakim
#    cv2.destroyAllWindows()

#    image_gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
#    cv2.imwrite('graytest.jpg',image_gray)


# mouse callback function
#    def mouse_callback(self,event,x,y,flags,param):
#        if event == cv2.EVENT_LBUTTONDOWN:
#             print x, y
