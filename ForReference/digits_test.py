import matplotlib.pyplot as plt
from sklearn import datasets, svm, metrics
import sys


def main(argv=None):
    digits = datasets.load_digits()
    print(type(digits))
    for key, value in digits.items():
        try:
            print(key, ": ", value.shape)
        except:
            print(key)
    images_and_labels = list(zip(digits.images, digits.target))
    print(images_and_labels)
    for index, (image, label) in enumerate(images_and_labels[:4]):
        plt.subplot(2, 4, index + 1)
        plt.axis('off')
        plt.imshow(image, cmap=plt.cm.gray_r, interpolation='nearest')
        plt.title('Training: %i' % label)
    plt.show()
if __name__ == "__main__":
    sys.exit(main())
