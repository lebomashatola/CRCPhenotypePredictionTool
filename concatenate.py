

from keras.layers import Conv2D, MaxPooling2D, Activation, Dropout, Flatten, Dense
from tensorflow.keras.utils import image_dataset_from_directory
from sklearn.model_selection import train_test_split
from tensorflow.keras.layers import concatenate
from sklearn.preprocessing import LabelEncoder
from keras.utils.vis_utils import plot_model
from sklearn.preprocessing import scale
from keras.models import Sequential
from sklearn import preprocessing
from tensorflow.keras import Model
from sklearn.utils import shuffle
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras
from keras import models
import pandas as pd
import numpy as np



################################################################################


from keras.layers import Input, Concatenate, Conv2D, Flatten, Dense
from keras.models import Model

train_ds = tf.keras.preprocessing.image_dataset_from_directory(
    "/Users/lebomash/downloads/histology_images/training",
    label_mode="binary",
    shuffle=True,
    subset=None,
    image_size=(224, 224),
    batch_size=16)

val_ds = tf.keras.preprocessing.image_dataset_from_directory(
    "/Users/lebomash/downloads/histology_images/validation",
    label_mode="binary",
    shuffle=True,
    subset=None,
    image_size=(224, 224),
    batch_size=16)

def process(path, label):

    df = pd.read_csv(path, low_memory=False, index_col=0)
    norm = np.linalg.norm(df)
    normal_array = df/norm
    df = normal_array.transpose()
    df['Y'] = label

    return(df)

def processed():

    R = process('/Users/lebomash/documents/phd_2022/counts/processed_files/Resistant.csv', 'R')
    S = process('/Users/lebomash/documents/phd_2022/counts/processed_files/Sensitive.csv', 'S')

    df_all = R.append(pd.DataFrame(data = S), ignore_index=True)

    return(df_all)


df_counts = processed()
df_counts = shuffle(df_counts)

encoder = LabelEncoder()
df_counts['Y'] = encoder.fit_transform(df_counts['Y'])

x = df_counts.drop(['Y'], axis=1)
y = df_counts['Y']

x = scale(x)
X_train, X_test, y_train, y_test = train_test_split(x, y, test_size = 0.2)



vector_shape = (X_train.shape[1], )
image_shape = (224, 224, 3)

image_input = Input(image_shape)
vector_input = Input(vector_shape)

x = Conv2D(32, (3,3))(image_input)
x = Flatten()(x)

x = Dense(8, activation="relu")(image_input)
x = Dense(4, activation="relu")(x)
x = Model(inputs=image_input, outputs=x)


y = Dense(64, activation="relu")(vector_input)
y = Dense(32, activation="relu")(y)
y = Dense(4, activation="relu")(y)
y = Model(inputs=vector_input, outputs=y)


combined = concatenate([x.output, y.output])


z = Dense(2, activation="relu")(combined)
z = Dense(1, activation="sigmoid")(z)


model = Model(inputs=[x.input, y.input], outputs=z)

plot_model(model, to_file='demo.pdf', show_shapes=True)


opt = Adam(lr=1e-3, decay=1e-3 / 200)
model.compile(loss="mean_absolute_percentage_error", optimizer=opt)



model.fit(
	x=[X_train, train_ds], y=y_train,
	validation_data=([X_test, val_ds], y_test),
	epochs=10, batch_size=16)

################################################################################


def CNN():
    model = Sequential()

    model.add(Conv2D(128, (2, 2), activation='relu', input_shape=(224,224,1)))
    model.add(MaxPooling2D((2, 2)))

    model.add(Flatten())
    model.add(Dense(64, activation='relu'))

    model.add(Dropout(0.5))
    model.add(Dense(16, activation='sigmoid'))

    return(model)


def MLP(X_train):

    model = Sequential()

    model.add(Dense(128, activation='relu', input_shape=(X_train.shape[1], )))
    model.add(Dense(64, activation='relu'))
    model.add(Dense(16, activation='relu'))

    return(model)




input_shape = (224, 224)
batch_size = 8

train_ds = image_dataset_from_directory(
    "/Users/lebomash/downloads/histology_images/training",
    label_mode="binary",
    color_mode="grayscale",
    shuffle=True,
    subset=None,
    image_size=input_shape,
    batch_size=batch_size)

val_ds = image_dataset_from_directory(
    "/Users/lebomash/downloads/histology_images/validation",
    label_mode="binary",
    color_mode="grayscale",
    shuffle=True,
    subset=None,
    image_size=input_shape,
    batch_size=batch_size)


def process(path, label):

    df = pd.read_csv(path, low_memory=False, index_col=0)
    norm = np.linalg.norm(df)
    normal_array = df/norm
    df = normal_array.transpose()
    df['Y'] = label

    return(df)

def processed():

    R = process('/Users/lebomash/documents/phd_2022/counts/processed_files/Resistant.csv', 'R')
    S = process('/Users/lebomash/documents/phd_2022/counts/processed_files/Sensitive.csv', 'S')
    df_all = R.append(pd.DataFrame(data = S), ignore_index=True)

    return(df_all)




df_counts = processed()
df_counts = shuffle(df_counts)

encoder = LabelEncoder()
df_counts['Y'] = encoder.fit_transform(df_counts['Y'])

x = df_counts.drop(['Y'], axis=1)
y = df_counts['Y']
x = scale(x)

X_train, X_test, y_train, y_test = train_test_split(x, y, test_size = 0.2)


mlp_model = MLP(X_train)
cnn_model = CNN()

combinedInput = concatenate([mlp_model.output, cnn_model.output])

x_layer = Dense(9, activation="relu")(combinedInput)
x_layer = Dense(1, activation="sigmoid")(x_layer)

model = Model(inputs=[mlp_model.input, cnn_model.input], outputs=x_layer)
plot_model(model, to_file='demo.pdf', show_shapes=True)


model.compile(optimizer=keras.optimizers.Adam(1e-3),
                loss="binary_crossentropy",
                metrics=["accuracy"])


model.fit([X_train, train_np], y_train.to_numpy(), epochs=5, batch_size=8)



unseen = model.evaluate([X_test, val_ds], y_test)
seen = model.evaluate([X_train, train_ds], y_train)

print(' Accuracy test data: ', round(unseen[1] * 100, 2),'%', '\n',
        'Accuracy train data: ',
        round(seen[1] * 100, 2),'%')
