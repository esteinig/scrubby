import os
import typer
import logging
import random
import numpy as np
import tensorflow as tf

from Bio import SeqIO
from typing import List

from tensorflow.keras import layers, models, optimizers
from tensorflow.keras.callbacks import EarlyStopping


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'


# Optionally, set the logging level for the Python logging module
logging.getLogger('tensorflow').setLevel(logging.ERROR)

# Set TensorFlow logging level
tf.get_logger().setLevel('ERROR')


INPUT_SIZE = 150  # Length of the DNA sequence
NUM_CLASSES = 5  # Number of classes
NUM_CHROMOSOMES = 25  # Assuming chromosomes 1-22, X, Y, and MT
DROPOUT_PROB = 0.5  # Dropout probability

class HybridModel(tf.keras.Model):
    def __init__(self, input_size, hidden_size, num_classes, aux_input_size=None, use_lstm=True, train=True):
        super(HybridModel, self).__init__()

        self.cnn = models.Sequential([
            layers.Conv1D(32, 3, activation='relu', input_shape=(input_size, 1)),
            layers.MaxPooling1D(2),
            layers.Conv1D(64, 3, activation='relu'),
            layers.MaxPooling1D(2)
        ])

        self.use_lstm = use_lstm

        if use_lstm:
            self.lstm = layers.Bidirectional(layers.LSTM(hidden_size, dropout=DROPOUT_PROB, return_sequences=True))
        else:
            self.lstm = None

        if aux_input_size is not None:
            self.aux_fc = layers.Dense(hidden_size)
        else:
            self.aux_fc = None

        self.fc = layers.Dense(num_classes)

    def call(self, inputs, aux_input=None):
        x = self.cnn(inputs)
        if self.use_lstm:
            x = self.lstm(x)

        if self.aux_fc is not None and aux_input is not None:
            aux_out = self.aux_fc(aux_input)
            x = tf.concat([tf.reduce_mean(x, axis=1), aux_out], axis=1)
        else:
            x = tf.reduce_mean(x, axis=1)

        logits = self.fc(x)
        return logits
    

# Helper functions
def one_hot_encode(labels, num_classes):
    return np.eye(num_classes)[labels]

def normalize_sequence(seq):
    """Normalize the DNA sequence to a numerical representation."""
    mapping = { 'A': 1.0, 'C': 2.0, 'G': 3.0, 'T': 4.0 }
    return np.array([mapping.get(base, 0.0) for base in seq], dtype=np.float32)

def get_label_from_filename(file_path):
    """Extract label from the filename."""
    filename = os.path.splitext(os.path.splitext(os.path.basename(file_path))[0])[0]
    label_pos = filename.rfind("__")
    if label_pos != -1:
        return int(filename[label_pos + 2:])
    raise ValueError("Invalid filename format for extracting label")

def load_alignment_info(file_path):
    """Load alignment information from a file."""
    alignment_info = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            read_id = parts[0]
            chromosome = int(parts[1])
            start = int(parts[2])
            end = int(parts[3])
            alignment_info[read_id] = (chromosome, start, end)
    return alignment_info

def load_sequences(file_path, alignment_info=None, num_chromosomes=NUM_CHROMOSOMES):
    seqs = []
    labels = []
    aux_inputs = []

    seq_label = get_label_from_filename(file_path)
    print(f"Label from filename is: {seq_label}")

    excluded = 0
    total = 0

    for record in SeqIO.parse(file_path, "fastq"):
        total += 1
        seq = normalize_sequence(str(record.seq))

        if len(seq) < INPUT_SIZE:
            excluded += 1
            continue


        # Correcting the shape to match the model input shape (None, 150, 1)
        seq_tensor = tf.convert_to_tensor(seq[:INPUT_SIZE], dtype=tf.float32)
        seq_tensor = tf.expand_dims(seq_tensor, axis=-1)  # Shape: (150, 1)
        

        if alignment_info:
            read_id = record.id
            if read_id in alignment_info:
                chromosome, start, end = alignment_info[read_id]
                chrom_tensor = tf.zeros([num_chromosomes], dtype=tf.float32)
                chrom_tensor = tf.tensor_scatter_nd_update(chrom_tensor, [[chromosome]], [1.0])

                start_tensor = tf.convert_to_tensor([start], dtype=tf.float32)
                end_tensor = tf.convert_to_tensor([end], dtype=tf.float32)
                aux_input = tf.concat([chrom_tensor, start_tensor, end_tensor], axis=0)
                aux_input = tf.expand_dims(aux_input, axis=0)

                aux_inputs.append(aux_input)

        seqs.append(seq_tensor)
        labels.append(tf.convert_to_tensor([seq_label], dtype=tf.int64))

    print(f"Excluded {excluded}/{total} sequences for not matching input size {INPUT_SIZE}")

    if aux_inputs:
        return np.array(seqs), np.array(labels), np.array(aux_inputs)
    else:
        return np.array(seqs), np.array(labels), None

def train_test_val_split(data_len, train_ratio, test_ratio):
    indices = list(range(data_len))
    random.shuffle(indices)
    train_end = int(data_len * train_ratio)
    test_end = train_end + int(data_len * test_ratio)
    return indices[:train_end], indices[train_end:test_end], indices[test_end:]

def evaluate(model, sequences, labels, aux_inputs=None):
    predictions = model(sequences, aux_inputs)
    predicted_classes = tf.argmax(predictions, axis=1)
    accuracy = tf.reduce_mean(tf.cast(predicted_classes == labels, tf.float32))
    return accuracy.numpy()

def train_nn(fastq_files, model_weights, alignment_data=None, epochs=10, batch_size=32):
    
    # Create a MirroredStrategy
    strategy = tf.distribute.MirroredStrategy()

    # Load sequences and labels
    all_sequences = []
    all_labels = []
    all_aux_inputs = []
    has_aux_inputs = False

    for file_path in fastq_files:
        sequences, labels, aux_inputs = load_sequences(file_path, alignment_data)
        all_sequences.append(sequences)
        all_labels.append(labels)
        if aux_inputs is not None:
            all_aux_inputs.append(aux_inputs)
            has_aux_inputs = True

    all_sequences = np.concatenate(all_sequences)
    all_labels = np.concatenate(all_labels)
    if has_aux_inputs:
        all_aux_inputs = np.concatenate(all_aux_inputs)
    else:
        all_aux_inputs = None

    train_indices, test_indices, val_indices = train_test_val_split(len(all_sequences), 0.7, 0.15)

    train_sequences = all_sequences[train_indices]
    test_sequences = all_sequences[test_indices]
    val_sequences = all_sequences[val_indices]

    train_labels = all_labels[train_indices]
    test_labels = all_labels[test_indices]
    val_labels = all_labels[val_indices]

    if has_aux_inputs:
        train_aux_inputs = all_aux_inputs[train_indices]
        test_aux_inputs = all_aux_inputs[test_indices]
        val_aux_inputs = all_aux_inputs[val_indices]
    else:
        train_aux_inputs = test_aux_inputs = val_aux_inputs = None

    aux_input_size = NUM_CHROMOSOMES + 2 if has_aux_inputs else None

    with strategy.scope():
        model = HybridModel(INPUT_SIZE, 128, NUM_CLASSES, aux_input_size)
        model.compile(optimizer=optimizers.Adam(1e-4), loss='sparse_categorical_crossentropy', metrics=['accuracy'])

    if train_aux_inputs is not None:
        train_data = (train_sequences, train_aux_inputs)
    else:
        train_data = train_sequences

    model.fit(train_data, train_labels, epochs=epochs, batch_size=batch_size, validation_data=(val_sequences, val_labels),
              callbacks=[tf.keras.callbacks.EarlyStopping(patience=3)])

    val_accuracy = evaluate(model, val_sequences, val_labels, val_aux_inputs)
    print(f'Final Validation Accuracy: {val_accuracy * 100:.2f}%')

    model.save_weights(model_weights)


# Prediction function
# def predict_nn(model_weights: str, fastq_files: List[str], alignment_data: str = None):
#     aux_input_size = NUM_CHROMOSOMES + 2 if alignment_data is not None else None
#     model = HybridModel(INPUT_SIZE, 128, NUM_CLASSES, aux_input_size)
#     model.load_weights(model_weights)

#     for fastq_path in fastq_files:
#         sequences, _, aux_inputs = load_sequences(fastq_path, alignment_data)
#         predictions = model(sequences, aux_inputs)
#         predicted_classes = tf.argmax(predictions, axis=1)
#         print(f'Predicted classes for {fastq_path}: {predicted_classes.numpy()}')


def predict_nn(model_weights: str, fastq_files: List[str], alignment_data: str = None):
    aux_input_size = NUM_CHROMOSOMES + 2 if alignment_data is not None else None
    model = HybridModel(INPUT_SIZE, 128, NUM_CLASSES, aux_input_size, True, False)
    model.load_weights(model_weights)

    for fastq_path in fastq_files:
        sequences, _, aux_inputs = load_sequences(fastq_path, alignment_info=None if alignment_data is None else load_alignment_info(alignment_data))
    

        # Make predictions on the aggregated data
        predictions = model(sequences, aux_inputs)
        predicted_classes = tf.argmax(predictions, axis=1).numpy()

        # Count the occurrences of each class
        class_counts = np.bincount(predicted_classes, minlength=NUM_CLASSES)

        total_predictions = len(predicted_classes)
        class_distribution = class_counts / total_predictions * 100  # Convert to percentage

        print(f"\nPredictions for {fastq_path}")
        for cls, count in enumerate(class_counts):
            print(f"Class {cls}: {count} predictions ({class_distribution[cls]:.2f}%)")
