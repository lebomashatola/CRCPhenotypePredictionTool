

import tensorflow as tf
from tensorflow import keras
import matplotlib.pyplot as plt
from tensorflow.keras import layers
from keras.layers import GaussianNoise

from keras.layers import Activation, Dropout, Flatten, Dense
from keras.preprocessing.image import ImageDataGenerator
from keras.layers import Conv2D, MaxPooling2D
from keras.models import Sequential


input_shape = (224, 224)
batch_size = 64

train_ds = tf.keras.preprocessing.image_dataset_from_directory(
    "/Users/lebmash",
    label_mode="binary",
    shuffle=True,
    subset=None,
    image_size=input_shape,
    batch_size=batch_size)


val_ds = tf.keras.preprocessing.image_dataset_from_directory(
    "",
    label_mode="binary",
    shuffle=True,
    subset=None,
    image_size=input_shape,
    batch_size=batch_size)


model = Sequential()
model.add(layers.Conv2D(16, (2, 2), input_shape=(224,224,3)))
model.add(layers.MaxPooling2D())
model.add(GaussianNoise(0.1))
model.add(Activation('relu'))

model.add(layers.Conv2D(8,(2, 2),activation='relu'))
model.add(layers.MaxPooling2D())

model.add(Dropout(0.5))
model.add(layers.Dense(1, activation='sigmoid'))

model.compile(optimizer=keras.optimizers.Adam(1e-3),
                loss="binary_crossentropy",
                metrics=["accuracy"])

history = model.fit(train_ds, epochs=150, validation_data=val_ds, verbose=1)
plot_model(model, to_file='/Users/lebomash/documents/phd_2022/results/cnn.pdf', show_shapes=True)

def plot(label1, label2, label3, history):

    plt.plot(history.history[label1], color='black', linewidth=3)
    plt.plot(history.history[label2], color='red', linewidth=3)

    plt.title(label3, fontsize=16)

    plt.legend(['training', 'validation'], loc='upper left')
    plt.xlabel('epoch', fontsize=14)
    plt.ylabel(label1, fontsize=14)

    plt.rc('xtick',labelsize=13)
    plt.rc('ytick',labelsize=13)

    plt.savefig('/Users/lebomash/documents/phd_2022/results/' + label1 + '.pdf')


plot('accuracy', 'val_accuracy', 'model accuracy', history)
plot('loss', 'val_loss', 'model loss', history)
