use tch::{nn, nn::Module, nn::OptimizerConfig, Device, Tensor, Kind, no_grad};
use needletail::{parse_fastx_file, Sequence};
use rand::seq::SliceRandom;
use rand::thread_rng;
use std::cmp;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use tch::nn::RNN;

use crate::error::ScrubbyError;

const INPUT_SIZE: i64 = 150; // Length of the DNA sequence
const NUM_CLASSES: i64 = 5; // Number of classes
const NUM_CHROMOSOMES: usize = 25; // Assuming chromosomes 1-22, X, Y, and MT
const DROPOUT_PROB: f64 = 0.5; // Dropout probability

enum AuxDataOption {
    Include,
    Exclude,
}


#[derive(Debug)]
struct HybridModel {
    cnn: nn::Sequential,
    _cnn_layers: Vec<String>,
    lstm: Option<nn::LSTM>,
    fc: nn::Linear,
    aux_fc: Option<nn::Linear>,
}

impl HybridModel {
    fn new(
        vs: &nn::Path,
        input_size: i64,
        hidden_size: i64,
        num_classes: i64,
        aux_input_size: Option<i64>,
        use_lstm: bool,
        train: bool
    ) -> Self {

        let _cnn_layers = Vec::from([
            "Conv1D(1, 32, 3)".to_string(),
            "ReLU".to_string(),
            "MaxPool1D(2)".to_string(),
            "Conv1D(32, 64, 3)".to_string(),
            "ReLU".to_string(),
            "MaxPool1D(2)".to_string()
        ]);

        let cnn = nn::seq()
            .add(nn::conv1d(vs / "cnn1", 1, 32, 3, Default::default()))
            .add_fn(|x| x.relu())
            .add_fn(|x| x.max_pool1d(2, 2, 0, 1, false))
            .add(nn::conv1d(vs / "cnn2", 32, 64, 3, Default::default()))
            .add_fn(|x| x.relu())
            .add_fn(|x| x.max_pool1d(2, 2, 0, 1, false));

        // Get the output length after the CNN layers
        let cnn_output_size = {
            let input = Tensor::zeros(&[1, 1, input_size], (tch::Kind::Float, Device::cuda_if_available()));
            let output = cnn.forward(&input);
            output.size()[1]  
        };

        log::info!("CNN output size for input dimension for LSTM: {}", cnn_output_size);

        let lstm = if use_lstm {
            Some(nn::lstm(
                vs / "lstm",
                cnn_output_size,
                hidden_size,
                nn::RNNConfig {
                    batch_first: true,
                    dropout: DROPOUT_PROB,
                    train,
                    bidirectional: true,
                    ..Default::default()
                },
            ))
        } else {
            None
        };

        let fc_input_size = if use_lstm {
            hidden_size * 2 + aux_input_size.unwrap_or(0) // bidirectional LSTM
        } else {
            cnn_output_size + aux_input_size.unwrap_or(0)
        };

        let fc = nn::linear(vs / "fc", fc_input_size, num_classes, Default::default());

        let aux_fc = aux_input_size.map(|size| nn::linear(vs / "aux_fc", size, hidden_size, Default::default()));

        HybridModel {
            cnn,
            _cnn_layers,
            lstm,
            fc,
            aux_fc,
        }
    }
    fn forward(&self, xs: &Tensor, aux_input: Option<&Tensor>) -> Tensor {
        let cnn_out = self.cnn.forward(&xs.view([xs.size()[0], 1, INPUT_SIZE]));

        // Print the shape of cnn_out
        // log::info!("CNN output shape before view: {:?}", cnn_out.size());

        // let batch_size = cnn_out.size()[0];
        // let features = cnn_out.size()[1];
        // let seq_len = cnn_out.size()[2];

        let cnn_out = cnn_out.permute(&[0, 2, 1]); //  reshape to [batch_size, seq_len, features] for LSTM

        // log::info!("CNN output shape after view: {:?}", cnn_out.size());
        
        let lstm_out = if let Some(lstm) = &self.lstm {
            let (output, _) = lstm.seq(&cnn_out);
            output
        } else {
            cnn_out
        };

        if let Some(aux_fc) = &self.aux_fc {
            if let Some(aux_input) = aux_input {
                let aux_out = aux_fc.forward(aux_input);
                let combined_input = Tensor::cat(&[lstm_out.squeeze(), aux_out], 1);

                let logits =  self.fc.forward(&combined_input);
                return logits
            }
        }
        
        let logits = self.fc.forward(&lstm_out.squeeze());
        return logits
    }

    fn forward_with_softmax(&self, xs: &Tensor, aux_input: Option<&Tensor>) -> Tensor {
        let logits = self.forward(xs, aux_input);
        logits.softmax(-1, Kind::Float)
    }
}


fn load_alignment_info(file_path: &Path) -> HashMap<String, (i64, i64, i64)> {
    let mut alignment_info = HashMap::new();
    let file = File::open(file_path).expect("Failed to open alignment file");
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line.expect("Failed to read line");
        let parts: Vec<&str> = line.split(',').collect();
        let read_id = parts[0].to_string();
        let chromosome = parts[1].parse::<i64>().expect("Invalid chromosome");
        let start = parts[2].parse::<i64>().expect("Invalid start position");
        let end = parts[3].parse::<i64>().expect("Invalid end position");
        alignment_info.insert(read_id, (chromosome, start, end));
    }

    alignment_info
}


fn get_label_from_filename(file_path: &PathBuf) -> Result<i64, ScrubbyError> {
    let filename = file_path.with_extension("").with_extension("");
    let filename = filename.file_name()
        .ok_or(ScrubbyError::ReadNeuralNetworkModelLabel)?
        .to_str()
        .ok_or(ScrubbyError::ReadNeuralNetworkModelLabel)?;
    
    if let Some(label_pos) = filename.rfind("__") {
        Ok(filename[label_pos + 2..].parse::<i64>().expect("Invalid numeric label <i64>"))
    } else {
        Err(ScrubbyError::ReadNeuralNetworkModelLabel)
    }
}

fn load_sequences(device: Device, file_path: &PathBuf, alignment_info: Option<&HashMap<String, (i64, i64, i64)>>, num_chromosomes: usize) -> Result<(Vec<Tensor>, Vec<Tensor>, Option<Vec<Tensor>>), ScrubbyError> {
    
    let mut seqs = Vec::new();
    let mut labels = Vec::new();
    let mut aux_inputs = Vec::new();

    let mut reader = parse_fastx_file(file_path).map_err(
        |_| ScrubbyError::ReadNeuralNetworkFastq(file_path.to_path_buf())
    )?;

    let seq_label = get_label_from_filename(file_path)?;

    while let Some(record) = reader.next() {
        let record = record?;
        let seq = record.normalize(false);

        if record.num_bases() < INPUT_SIZE as usize {
            log::warn!("Read is smaller with {} bp than expected input size of {} bp", record.num_bases(), INPUT_SIZE);
            continue;
        }
            
        let seq_tensor = Tensor::from_slice(&seq)
            .to_device(device)
            .to_kind(tch::Kind::Float)
            .unsqueeze(0)
            .unsqueeze(0);
        
        if let Some(alignment_info) = alignment_info {
            let read_id = std::str::from_utf8(record.id())?.to_string();

            if let Some(&(chromosome, start, end)) = alignment_info.get(&read_id) {
                let chrom_tensor = Tensor::zeros(
                    &[num_chromosomes as i64], 
                    (Kind::Float, device)
                )
                .narrow(0, chromosome, 1)
                .fill_(1.0);

                let start_tensor = Tensor::from_slice(&[start as f32])
                    .to_device(device);
                let end_tensor = Tensor::from_slice(&[end as f32])
                    .to_device(device);
                let aux_input = Tensor::cat(&[chrom_tensor, start_tensor, end_tensor], 0)
                    .unsqueeze(0);

                aux_inputs.push(aux_input);
            }
        }

        seqs.push(seq_tensor);

        labels.push(
            Tensor::from_slice(&[seq_label])
                .to_device(device)
                .to_kind(tch::Kind::Int64)
        );
    }

    if !aux_inputs.is_empty() {
        Ok((seqs, labels, Some(aux_inputs)))
    } else {
        Ok((seqs, labels, None))
    }
}

fn predict(model: &HybridModel, seqs: Vec<Tensor>, aux_inputs: Option<Vec<Tensor>>) -> i64 {
    let mut all_predictions = Vec::new();

    no_grad(|| {
        for (i, seq) in seqs.into_iter().enumerate() {
            let aux_input = aux_inputs.as_ref().map(|aux| &aux[i]);
            let logits = model.forward_with_softmax(&seq, aux_input);
            all_predictions.push(logits);
        }
    });

    let all_predictions = Tensor::cat(&all_predictions, 0);
    let average_predictions = all_predictions.mean_dim(Some([0].as_ref()), true, Kind::Float);
    let probabilities = average_predictions.softmax(1, Kind::Float);
    let final_prediction = probabilities.argmax(-1, false);

    log::info!("Average predictions: {}", average_predictions);
    log::info!("Final prediction: {}", final_prediction);

    final_prediction.int64_value(&[])
}


fn one_hot_encode(device: Device, labels: &Tensor, num_classes: i64, kind: Kind) -> Tensor {
    let batch_size = labels.size()[0];
    let mut one_hot = Tensor::zeros(&[batch_size, num_classes], (kind, device));
    one_hot = one_hot.scatter_(
        1, 
        &labels.unsqueeze(1), 
        &Tensor::ones(
            &[batch_size, 1], 
            (kind, device)
        )
    );
    one_hot
}

fn train(
    model: &HybridModel,
    device: Device,
    vs: &nn::VarStore,
    sequences: &[Tensor],
    labels: &[Tensor],
    aux_inputs: Option<&[Tensor]>,
    test_sequences: &[Tensor],
    test_labels: &[Tensor],
    test_aux_inputs: Option<&[Tensor]>,
    epochs: i64,
    batch_size: usize,
) {
    let mut optimizer = nn::Adam::default().build(&vs, 1e-4).unwrap();

    for epoch in 0..epochs {
        let mut batch_indices: Vec<usize> = (0..sequences.len()).collect();
        batch_indices.shuffle(&mut thread_rng());

        for batch_start in (0..sequences.len()).step_by(batch_size) {
            let batch_end = cmp::min(batch_start + batch_size, sequences.len());

            let batch_seqs: Vec<_> = batch_indices[batch_start..batch_end]
                .iter()
                .map(|&i| sequences[i].unsqueeze(0))
                .collect();
            let batch_labels: Vec<_> = batch_indices[batch_start..batch_end]
                .iter()
                .map(|&i| labels[i].unsqueeze(0))
                .collect();


            let batch_seqs = Tensor::cat(&batch_seqs, 0);
            let batch_labels = Tensor::cat(&batch_labels, 0).squeeze_dim(1);

            let output = if let Some(aux_inputs) = aux_inputs {
                let batch_aux: Vec<_> = batch_indices[batch_start..batch_end]
                    .iter()
                    .map(|&i| aux_inputs[i].unsqueeze(0))
                    .collect();
                let batch_aux = Tensor::cat(&batch_aux, 0);
                model.forward(&batch_seqs, Some(&batch_aux))
            } else {
                model.forward(&batch_seqs, None)
            };

            let loss = output.cross_entropy_loss(&one_hot_encode(device, &batch_labels, NUM_CLASSES, Kind::Int64), None::<&Tensor>, tch::Reduction::Mean, -100, 0.0);

            optimizer.zero_grad();
            loss.backward();
            optimizer.step();

            log::info!("Epoch: {}, Loss: {}", epoch, loss.double_value(&[]));
        }

        // Evaluate on the test set after each epoch
        let test_accuracy = evaluate(model, test_sequences, test_labels, test_aux_inputs);
        log::info!("Epoch: {}, Test Accuracy: {:.2}%", epoch, test_accuracy * 100.0);
    }
}

pub fn train_nn(
    device: usize,
    fastq_files: Vec<PathBuf>,
    model_weights: PathBuf,
    alignment_data: Option<PathBuf>,
    epochs: i64,
    batch_size: usize,
) -> Result<(), ScrubbyError> {

    let device = Device::Cuda(device);
    
    let vs = nn::VarStore::new(device);

    log::info!("Device is CUDA: {:?}", device.is_cuda());

    let aux_data_option = match alignment_data {
        Some(_) => AuxDataOption::Include,
        None => AuxDataOption::Exclude,
    };
    let aux_input_size = match aux_data_option {
        AuxDataOption::Include => Some((NUM_CHROMOSOMES + 2) as i64),
        AuxDataOption::Exclude => None,
    };

    let model = HybridModel::new(&vs.root(), INPUT_SIZE, 128, NUM_CLASSES, aux_input_size, true, true);

    let alignment_info = if matches!(aux_data_option, AuxDataOption::Include) {
        Some(load_alignment_info(&alignment_data.ok_or(ScrubbyError::ReadNeuralNetworkModel)?))
    } else {
        None
    };

    let mut all_sequences = Vec::new();
    let mut all_labels = Vec::new();
    let mut all_aux_inputs = Vec::new();
    let mut has_aux_inputs = false;

    for file_path in fastq_files {
        let (sequences, labels, aux_inputs) = load_sequences(device, &file_path, alignment_info.as_ref(), NUM_CHROMOSOMES)?;
        all_sequences.extend(sequences);
        all_labels.extend(labels);
        if let Some(aux) = aux_inputs {
            all_aux_inputs.extend(aux);
            has_aux_inputs = true;
        }
    }

    let aux_inputs = if has_aux_inputs { Some(all_aux_inputs) } else { None };

    // Get the indices for train, test, and validation splits
    let (train_indices, test_indices, val_indices) = train_test_val_split(all_sequences.len(), 0.7, 0.15);

    // Function to gather tensors based on indices
    let gather_tensors = |indices: &[usize], data: &[Tensor]| -> Vec<Tensor> {
        indices.iter().map(|&i| data[i].shallow_clone()).collect()
    };

    let train_sequences = gather_tensors(&train_indices, &all_sequences);
    let test_sequences = gather_tensors(&test_indices, &all_sequences);
    let val_sequences = gather_tensors(&val_indices, &all_sequences);

    let train_labels = gather_tensors(&train_indices, &all_labels);
    let test_labels = gather_tensors(&test_indices, &all_labels);
    let val_labels = gather_tensors(&val_indices, &all_labels);

    let train_aux_inputs = aux_inputs.as_ref().map(|aux| gather_tensors(&train_indices, aux));
    let test_aux_inputs = aux_inputs.as_ref().map(|aux| gather_tensors(&test_indices, aux));
    let val_aux_inputs = aux_inputs.as_ref().map(|aux| gather_tensors(&val_indices, aux));

    train(
        &model,
        device,
        &vs,
        &train_sequences,
        &train_labels,
        train_aux_inputs.as_deref(),
        &test_sequences,
        &test_labels,
        test_aux_inputs.as_deref(),
        epochs,
        batch_size,
    );

    // Evaluate on the validation set
    let val_accuracy = evaluate(&model, &val_sequences, &val_labels, val_aux_inputs.as_deref());
    log::info!("Final Validation Accuracy: {:.2}%", val_accuracy * 100.0);

    vs.save(model_weights).map_err(|_| ScrubbyError::SaveNeuralNetworkModel)?;

    Ok(())
}



pub fn predict_nn(device: usize, model_weights: PathBuf, fastq: Vec<PathBuf>, alignment_data: Option<PathBuf>) -> Result<(), ScrubbyError>{

    let device = Device::Cuda(device);
    let mut vs = nn::VarStore::new(device);

    let aux_data_option = AuxDataOption::Exclude; // Set to AuxDataOption::Exclude to exclude auxiliary alignment data
    let aux_input_size = match aux_data_option {
        AuxDataOption::Include => Some((NUM_CHROMOSOMES + 2) as i64),
        AuxDataOption::Exclude => None,
    };

    let model = HybridModel::new(&vs.root(), INPUT_SIZE, 128, NUM_CLASSES, aux_input_size, true, false);
    
    vs.load(model_weights).expect("Failed to load model weights");
    
    let alignment_info = if matches!(aux_data_option, AuxDataOption::Include) {
        Some(load_alignment_info(&alignment_data.ok_or(ScrubbyError::ReadNeuralNetworkModel)?))
    } else {
        None
    };

    for fastq_path in fastq {

        log::info!("Loading read tensors: {}", fastq_path.display());
        let (seqs, _, aux_inputs) = load_sequences(device, &fastq_path, alignment_info.as_ref(), NUM_CHROMOSOMES)?;
        let final_class = predict(&model, seqs, aux_inputs);
        log::info!("Predicted class: {}", final_class);

    }

    Ok(())
}


pub fn check_gpu_connectivity() -> bool {
    // Create a simple tensor
    let tensor = Tensor::ones(&[1], (tch::Kind::Float, Device::Cpu));
    
    // Try to move the tensor to the GPU
    let result = tensor.to_device(Device::Cuda(0));
    
    // Check if the device of the result tensor is GPU
    result.device() == Device::Cuda(0)
}

fn train_test_val_split(data_len: usize, train_ratio: f64, test_ratio: f64) -> (Vec<usize>, Vec<usize>, Vec<usize>) {
    let mut rng = thread_rng();
    let mut indices: Vec<usize> = (0..data_len).collect();
    indices.shuffle(&mut rng);

    let train_end = (data_len as f64 * train_ratio).round() as usize;
    let test_end = train_end + (data_len as f64 * test_ratio).round() as usize;

    let train_indices = indices[..train_end].to_vec();
    let test_indices = indices[train_end..test_end].to_vec();
    let val_indices = indices[test_end..].to_vec();

    (train_indices, test_indices, val_indices)
}

fn evaluate(model: &HybridModel, sequences: &[Tensor], labels: &[Tensor], aux_inputs: Option<&[Tensor]>) -> f64 {
    let mut correct = 0;
    let mut total = 0;

    no_grad(|| {
        for (i, seq) in sequences.iter().enumerate() {
            let aux_input = aux_inputs.map(|aux| &aux[i]);

            let logits = model.forward(seq, aux_input);
            log::debug!("Output: {}", logits);
            
            // Average the logits across the batch for each sequence
            let average_logits = logits.mean_dim(Some(&[0_i64][..]), false, tch::Kind::Float);
            log::debug!("Average Logits: {}", average_logits);

            // Apply softmax to get probabilities
            let probabilities = average_logits.softmax(0, Kind::Float);
            log::debug!("Probabilities: {}", probabilities);

            // Get the predicted class by taking argmax
            let predicted = probabilities.argmax(0, true);
            log::debug!("Predicted: {}", predicted);

            // Ensure target is a single value tensor
            let target = &labels[i];
            log::debug!("Target: {}", target);

            if predicted == *target {
                correct += 1;
            }
            total += 1;
        }
    });

    correct as f64 / total as f64
}