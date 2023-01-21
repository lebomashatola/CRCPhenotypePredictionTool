

from sklearn.model_selection import train_test_split
from tensorflow.keras.layers import Activation, Dense
from tensorflow.keras.models import Sequential
from tensorflow.keras.optimizers import Adam
import tensorflow as tf


#Transcription

df_blue = '/Users/lebo/Documents/PhD_2022/ML/GeneExpression/Blue.csv'

#Histology

train_path = '<dir>'
valid_path = '<dir>'
test_path = '<dir>'

train_batch = ImageDataGenerator(preprocessing_function=self.image_gen) \
.flow_from_directory(directory=train_path, target_size=(224,224),
classes=['resistant', 'sensitive'], batch_size=30)

test_batch = ImageDataGenerator(preprocessing_function=self.image_gen) \
.flow_from_directory(directory=test_path, target_size=(224,224),
classes=['resistant', 'sensitive'], batch_size=30)

validate_batch = ImageDataGenerator(preprocessing_function=self.image_gen) \
.flow_from_directory(directory=valid_path, target_size=(224,224),
classes=['resistant', 'sensitive'], batch_size=30)
